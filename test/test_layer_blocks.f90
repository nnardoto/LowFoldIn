! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Portão 1 (fase 1): blocos de camada principal.
!>  - H00 + H01 e^{2 pi i k_s} + H10 e^{-2 pi i k_s} = compute_Hk_at_kpoint a 1e-12,
!>    em modelos sintéticos (três eixos de empilhamento) e nos exemplos .hsio;
!>  - erro explícito quando |n_s| > pl_cells; reagrupamento de P células numa
!>    cadeia com segundo vizinho preserva as bandas;
!>  - fases do gauge atômico com a_3 inclinado;
!>  - zlinalg: LU com vários lados direitos contra multiplicação direta.
program test_layer_blocks
    use types, only: dp
    use HS_io, only: HSDataFormat
    use structure_io, only: structure_t
    use core, only: compute_Hk_at_kpoint
    use layer_blocks, only: pl_blocks_t, build_pl_blocks, pl_blocks_at_k, attach_pl_geometry, &
                            layer_gauge_phases, kpar_cartesian
    use zlinalg, only: lu_factor, lu_solve, zmm, herm_check
    use green_models, only: make_model, make_structure, rnd, set_seed, slab_eigen
    use testing, only: check, check_close, test_summary
    implicit none

    real(dp), parameter :: TWO_PI = 6.283185307179586476925286766559_dp
    type(HSDataFormat) :: hh, hs
    type(pl_blocks_t) :: pl
    integer :: ierr, axis
    character(len=:), allocatable :: errmsg

    call test_zlinalg()

    ! --- modelo sintético não ortogonal complexo, três eixos ---
    call make_model([2, 1, 3], [1, 1, 1], 0.3_dp, 0.03_dp, .true., 7, hh, hs)
    do axis = 1, 3
        call build_pl_blocks(hh, axis, 1, pl, ierr, errmsg, hs=hs)
        call check(ierr == 0, 'build_pl_blocks synthetic axis')
        call check(pl%norb == 6 .and. pl%ncell == 1, 'PL size for pl_cells = 1')
        call gate1(hh, hs, pl, 'synthetic')
    end do

    ! --- exemplos existentes ---
    call example('examples/UnfoldChain1D/Chain_H.hsio', 'examples/UnfoldChain1D/Chain_S.hsio')
    call example('examples/UnfoldDiatomic/Diatomic_H.hsio', 'examples/UnfoldDiatomic/Diatomic_S.hsio')
    call example('examples/UnfoldSpinorNonCollinear/Spinor_H.hsio', 'examples/UnfoldSpinorNonCollinear/Spinor_S.hsio')

    ! --- segundo vizinho: erro explícito e reagrupamento ---
    call make_model([1, 1], [0, 0, 2], 0.4_dp, 0.02_dp, .true., 11, hh, hs)
    call build_pl_blocks(hh, 3, 1, pl, ierr, errmsg, hs=hs)
    call check(ierr == 2, 'pl_cells = 1 with |n_s| = 2 is an error')
    call check(index(errmsg, 'pl_cells >= 2') > 0, 'error message suggests pl_cells >= 2')
    call build_pl_blocks(hh, 3, 0, pl, ierr, errmsg, hs=hs)
    call check(ierr == 0 .and. pl%ncell == 2 .and. pl%norb == 4, 'automatic pl_cells = 2')
    call regroup_bands(hh, hs, pl)
    call build_pl_blocks(hh, 3, 3, pl, ierr, errmsg, hs=hs)
    call check(ierr == 0 .and. pl%ncell == 3, 'explicit pl_cells = 3')
    call regroup_bands(hh, hs, pl)

    ! --- R com blocos nulos (cópias de célula do OpenMX além do raio de corte) ---
    call make_model([2, 1], [1, 1, 3], 0.3_dp, 0.03_dp, .true., 13, hh, hs)
    block
        integer :: ir
        do ir = 1, hh%Neighbors
            if (abs(hh%IndexList(3, ir)) > 1) then
                hh%Data(:, :, ir) = (0.0_dp, 0.0_dp)
                hs%Data(:, :, ir) = (0.0_dp, 0.0_dp)
            end if
        end do
    end block
    call build_pl_blocks(hh, 3, 0, pl, ierr, errmsg, hs=hs)
    call check(ierr == 0 .and. pl%ncell == 1 .and. pl%max_ns == 1, 'zero blocks at |n_s| = 3 ignored')
    call check(pl%nr == 9, 'zero blocks do not enter the R_par list')
    call gate1(hh, hs, pl, 'zero-padded')

    ! --- base ortonormal (sem S) ---
    call make_model([2], [1, 0, 1], 0.3_dp, 0.0_dp, .true., 5, hh, hs)
    call build_pl_blocks(hh, 3, 1, pl, ierr, errmsg)
    call check(ierr == 0 .and. .not. pl%has_overlap, 'orthonormal PL blocks')
    block
        complex(dp) :: h(2, 2, -1:1), s(2, 2, -1:1)
        call pl_blocks_at_k(pl, [0.1_dp, 0.0_dp], h, s)
        call check(abs(s(1, 1, 0) - 1.0_dp) < 1e-15_dp .and. maxval(abs(s(:, :, 1))) == 0.0_dp, &
                   'S = I without overlap file')
    end block

    call test_gauge()

    call test_summary('test_layer_blocks')

contains

    subroutine test_zlinalg()
        complex(dp) :: a(5, 5), lu(5, 5), x(5, 3), b(5, 3), ax(5, 3)
        integer, allocatable :: ipiv(:)
        integer :: i, j
        call set_seed(3)
        do j = 1, 5
            do i = 1, 5
                a(i, j) = cmplx(rnd(), rnd(), dp)
            end do
            a(j, j) = a(j, j) + 3.0_dp
        end do
        do j = 1, 3
            do i = 1, 5
                b(i, j) = cmplx(rnd(), rnd(), dp)
            end do
        end do
        lu = a
        call lu_factor(lu, ipiv, ierr, errmsg)
        call check(ierr == 0, 'lu_factor')
        x = b
        call lu_solve(lu, ipiv, x, ierr, errmsg)
        call zmm('N', 'N', a, x, ax)
        call check_close(maxval(abs(ax - b)), 0.0_dp, 1e-13_dp, 'LU solve, 3 right-hand sides')
        x = b
        call lu_solve(lu, ipiv, x, ierr, errmsg, trans='C')
        call zmm('C', 'N', a, x, ax)
        call check_close(maxval(abs(ax - b)), 0.0_dp, 1e-13_dp, 'LU solve with A^H')
        lu = (0.0_dp, 0.0_dp)
        call lu_factor(lu, ipiv, ierr, errmsg)
        call check(ierr == 2, 'singular matrix detected')
        call check(herm_check(a + conjg(transpose(a))) < 1e-15_dp, 'herm_check')
    end subroutine test_zlinalg

    subroutine gate1(hh, hs, pl, tag)
        type(HSDataFormat), intent(in) :: hh, hs
        type(pl_blocks_t),  intent(in) :: pl
        character(len=*),   intent(in) :: tag
        complex(dp), allocatable :: hk(:, :), sk(:, :), h(:, :, :), s(:, :, :), hb(:, :), sb(:, :)
        real(dp) :: k(3), dev, devs
        integer :: it, n
        complex(dp) :: ph
        n = pl%norb
        allocate(hk(n, n), sk(n, n), h(n, n, -1:1), s(n, n, -1:1))
        dev = 0.0_dp; devs = 0.0_dp
        call set_seed(99)
        do it = 1, 6
            k = 0.5_dp * [rnd(), rnd(), rnd()]
            call compute_Hk_at_kpoint(hh%Data, hh%IndexList, hk, k)
            call compute_Hk_at_kpoint(hs%Data, hs%IndexList, sk, k)
            call pl_blocks_at_k(pl, [k(pl%inplane(1)), k(pl%inplane(2))], h, s)
            ph = exp(cmplx(0.0_dp, TWO_PI * k(pl%stack_axis), dp))
            hb = h(:, :, 0) + h(:, :, 1) * ph + h(:, :, -1) / ph
            sb = s(:, :, 0) + s(:, :, 1) * ph + s(:, :, -1) / ph
            dev = max(dev, maxval(abs(hb - hk)))
            devs = max(devs, maxval(abs(sb - sk)))
            call check(maxval(abs(h(:, :, -1) - conjg(transpose(h(:, :, 1))))) < 1e-13_dp, &
                       tag // ': H10 = H01^dagger')
        end do
        call check_close(dev, 0.0_dp, 1e-12_dp, tag // ': H00 + H01 e^{ik} + H10 e^{-ik} = H(k)')
        call check_close(devs, 0.0_dp, 1e-12_dp, tag // ': same for S')
    end subroutine gate1

    subroutine example(hfile, sfile)
        character(len=*), intent(in) :: hfile, sfile
        type(HSDataFormat) :: eh, es
        type(pl_blocks_t) :: epl
        integer :: ax, nused
        call eh%ReadFile(hfile)
        call es%ReadFile(sfile)
        nused = 0
        do ax = 1, 3
            if (maxval(abs(eh%IndexList(ax, :))) > 1) cycle
            call build_pl_blocks(eh, ax, 1, epl, ierr, errmsg, hs=es)
            call check(ierr == 0, 'example build ' // hfile)
            call gate1(eh, es, epl, hfile)
            nused = nused + 1
        end do
        call check(nused > 0, 'example has an axis with |n| <= 1: ' // hfile)
    end subroutine example

    !> Bandas de PL com P células em K = união das bandas em k_s = (K + j)/P.
    subroutine regroup_bands(hh, hs, pl)
        type(HSDataFormat), intent(in) :: hh, hs
        type(pl_blocks_t),  intent(in) :: pl
        complex(dp), allocatable :: h(:, :, :), s(:, :, :), hb(:, :), sb(:, :), hk(:, :), sk(:, :), c(:, :)
        real(dp), allocatable :: epl(:), eref(:), e1(:)
        real(dp) :: kk, k3
        integer :: j, n, nb, P
        complex(dp) :: ph
        n = pl%norb; nb = pl%nbasis; P = pl%ncell
        allocate(h(n, n, -1:1), s(n, n, -1:1), epl(n), eref(n), e1(nb), c(n, n), hk(nb, nb), sk(nb, nb))
        kk = 0.137_dp
        call pl_blocks_at_k(pl, [0.21_dp, -0.05_dp], h, s)
        ph = exp(cmplx(0.0_dp, TWO_PI * kk, dp))
        hb = h(:, :, 0) + h(:, :, 1) * ph + h(:, :, -1) / ph
        sb = s(:, :, 0) + s(:, :, 1) * ph + s(:, :, -1) / ph
        call slab_eigen(hb, sb, epl, c)
        do j = 0, P - 1
            k3 = (kk + real(j, dp)) / real(P, dp)
            call compute_Hk_at_kpoint(hh%Data, hh%IndexList, hk, [0.21_dp, -0.05_dp, k3])
            call compute_Hk_at_kpoint(hs%Data, hs%IndexList, sk, [0.21_dp, -0.05_dp, k3])
            deallocate(c); allocate(c(nb, nb))
            call slab_eigen(hk, sk, e1, c)
            eref(j * nb + 1:(j + 1) * nb) = e1
            deallocate(c); allocate(c(n, n))
        end do
        call sort(eref)
        call check_close(maxval(abs(epl - eref)), 0.0_dp, 1e-12_dp, 'regrouped PL bands = folded bands')
    end subroutine regroup_bands

    subroutine sort(a)
        real(dp), intent(inout) :: a(:)
        integer :: i, j
        real(dp) :: t
        do i = 2, size(a)
            t = a(i)
            j = i - 1
            do while (j >= 1)
                if (a(j) <= t) exit
                a(j + 1) = a(j)
                j = j - 1
            end do
            a(j + 1) = t
        end do
    end subroutine sort

    !> Gauge atômico com a_3 inclinado: e^{-ik.r_i} H_{0,p} e^{ik.r_j} (fases de
    !> layer_gauge_phases) = sum_R H(R) e^{i k_par.(R_par + tau_j - tau_i)}.
    subroutine test_gauge()
        type(structure_t) :: st
        type(pl_blocks_t) :: gp
        real(dp) :: lat(3, 3), kp(2), kc(3), rv(3), dev
        complex(dp), allocatable :: h(:, :, :), s(:, :, :), ph0(:), php(:), ref(:, :)
        integer :: p, ir, i, j, n
        call make_model([2, 1], [1, 1, 1], 0.3_dp, 0.03_dp, .true., 21, hh, hs)
        lat(1, :) = [3.0_dp, 0.0_dp, 0.0_dp]
        lat(2, :) = [-1.5_dp, 2.6_dp, 0.0_dp]
        lat(3, :) = [1.1_dp, 0.7_dp, 6.0_dp]           ! a_3 inclinado
        call make_structure(lat, 2, 'X', 4, st)
        call build_pl_blocks(hh, 3, 1, gp, ierr, errmsg, hs=hs)
        call attach_pl_geometry(gp, st, ierr, errmsg)
        call check(ierr == 0 .and. gp%has_geometry, 'attach_pl_geometry')
        call check(abs(dot_product(gp%normal, [0.0_dp, 0.0_dp, 1.0_dp]) - 1.0_dp) < 1e-14_dp, 'normal')
        n = gp%norb
        allocate(h(n, n, -1:1), s(n, n, -1:1), ph0(n), php(n), ref(n, n))
        kp = [0.17_dp, -0.31_dp]
        kc = kpar_cartesian(gp, kp)
        call check(abs(dot_product(kc, gp%normal)) < 1e-14_dp, 'k_par lies in the plane')
        call check_close(dot_product(kc, lat(1, :)), TWO_PI * kp(1), 1e-13_dp, 'k_par . a1 = 2 pi k1')
        call pl_blocks_at_k(gp, kp, h, s)
        call layer_gauge_phases(gp, kp, 0, ph0)
        dev = 0.0_dp
        do p = -1, 1
            call layer_gauge_phases(gp, kp, p, php)
            ref = (0.0_dp, 0.0_dp)
            do ir = 1, hh%Neighbors
                if (hh%IndexList(3, ir) /= p) cycle
                rv = matmul(real(hh%IndexList(:, ir), dp), lat)
                do j = 1, n
                    do i = 1, n
                        ref(i, j) = ref(i, j) + hh%Data(i, j, ir) * &
                            exp(cmplx(0.0_dp, dot_product(kc, rv + gp%tau(:, j) - gp%tau(:, i)), dp))
                    end do
                end do
            end do
            do j = 1, n
                do i = 1, n
                    dev = max(dev, abs(conjg(ph0(i)) * h(i, j, p) * php(j) - ref(i, j)))
                end do
            end do
        end do
        call check_close(dev, 0.0_dp, 1e-12_dp, 'atomic-gauge blocks with tilted a3')
    end subroutine test_gauge

end program test_layer_blocks
