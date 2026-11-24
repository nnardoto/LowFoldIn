! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Portão 3 (fase 3), nível de módulo: região central, RGF e projeções.
!>  - RGF (G_mm, G_{m,m±1}, blocos de borda) e K x = w contra a inversa densa;
!>  - teste A/A: camada central entre dois meios do mesmo bulk recupera a DOS de
!>    bulk (referência Python, integração em k_s) a 1e-8; só Tr(G_mm S_mm) não;
!>  - slab espesso: G_11 de 400 PLs = g_s semi-infinita;
!>  - filme desacoplado (H_AB = S_AB = 0) reproduz os autovalores do slab finito;
!>  - deslocamento rígido H + Delta S move o espectro exatamente por Delta;
!>  - mapa coerente: I_coh >= 0 e polos nos autovalores do slab para qualquer w;
!>  - Löwdin e Mulliken têm o mesmo traço.
program test_central_green
    use types, only: dp
    use zlinalg, only: lu_factor, lu_solve_identity, zmm
    use surface_green, only: surface_gf, SIDE_TOP, SIDE_BOTTOM
    use central_green, only: central_t, rgf_t, central_init, assemble_central, attach_lead, rgf_solve, &
                             solve_projection
    use green_projection, only: orbital_dos_mulliken, orbital_dos_lowdin, coherent_map, dense_central, &
                                min_neg_im_gdiag, region_dos
    use spectral_projection, only: orbital_groups_t
    use green_models, only: read_cmat, rnd, set_seed, slab_eigen
    use testing, only: check, check_close, test_summary
    implicit none

    real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp
    integer :: n, ne
    real(dp) :: eta
    complex(dp), allocatable :: h00(:, :), h01(:, :), s00(:, :), s01(:, :)
    real(dp), allocatable :: e_ref(:), rho_ref(:), rho_diag(:)
    integer :: ierr
    character(len=:), allocatable :: errmsg

    call load_chain()
    call rgf_vs_dense()
    call aa_test()
    call thick_slab()
    call decoupled_film_and_coherent()
    call rigid_shift()
    call lowdin_trace()
    call test_summary('test_central_green')

contains

    subroutine load_chain()
        integer :: u, ie
        character(len=256) :: line
        complex(dp), allocatable :: tmp(:, :)
        real(dp) :: e, remb
        open(newunit=u, file='test/data/green/nonorth_chain.txt', status='old', action='read')
        do
            read(u, '(A)') line
            if (line(1:1) /= '#') exit
        end do
        read(line, *) n, ne, eta
        allocate(h00(n, n), h01(n, n), s00(n, n), s01(n, n), tmp(n, n))
        call read_cmat(u, h00)
        call read_cmat(u, h01)
        call read_cmat(u, s00)
        call read_cmat(u, s01)
        close(u)
        open(newunit=u, file='test/data/green/aa_reference.txt', status='old', action='read')
        do
            read(u, '(A)') line
            if (line(1:1) /= '#') exit
        end do
        read(line, *) ne
        allocate(e_ref(ne), rho_ref(ne), rho_diag(ne))
        do ie = 1, ne
            read(u, *) e, rho_ref(ie), remb, rho_diag(ie)
            e_ref(ie) = e
            rho_ref(ie) = 0.5_dp * (rho_ref(ie) + remb)   ! média k-int / embedding (iguais a 1e-15)
        end do
        close(u)
    end subroutine load_chain

    !> Região central de nl PLs da cadeia; leads: 'none', 'bot', 'top', 'both';
    !> v: potencial uniforme; coupled = .false. desacopla os meios (contato 'none').
    subroutine chain_central(nl, z, v, leads, sys, coupled)
        integer,          intent(in)  :: nl
        complex(dp),      intent(in)  :: z
        real(dp),         intent(in)  :: v
        character(len=*), intent(in)  :: leads
        type(central_t),  intent(out) :: sys
        logical, optional, intent(in) :: coupled
        complex(dp), allocatable :: g(:, :), k00(:, :), k01(:, :), k10(:, :), hc(:, :), sc(:, :)
        integer :: m, it
        real(dp) :: res
        logical :: cpl
        cpl = .true.
        if (present(coupled)) cpl = coupled
        call central_init(sys, [(n, m = 1, nl)])
        do m = 1, nl
            if (m < nl) then
                call assemble_central(sys, m, z, h00, s00, v, hu=h01, su=s01, vnext=v)
            else
                call assemble_central(sys, m, z, h00, s00, v)
            end if
        end do
        k00 = (z - v) * s00 - h00
        k01 = (z - v) * s01 - h01
        k10 = (z - v) * conjg(transpose(s01)) - conjg(transpose(h01))
        allocate(g(n, n))
        hc = h01; sc = s01
        if (.not. cpl) then
            hc = (0.0_dp, 0.0_dp); sc = (0.0_dp, 0.0_dp)
        end if
        if (leads == 'bot' .or. leads == 'both') then
            call surface_gf(k00, k01, k10, SIDE_TOP, 1.0e-14_dp, 300, g, it, res, ierr, errmsg)
            call attach_lead(sys, 'bot', z, hc, sc, v, v, g)
        end if
        if (leads == 'top' .or. leads == 'both') then
            call surface_gf(k00, k01, k10, SIDE_BOTTOM, 1.0e-14_dp, 300, g, it, res, ierr, errmsg)
            call attach_lead(sys, 'top', z, hc, sc, v, v, g)
        end if
    end subroutine chain_central

    subroutine rgf_vs_dense()
        type(central_t) :: sys
        type(rgf_t) :: res
        integer :: bs(4) = [3, 2, 4, 3], m, i, j, nt, i0, i1, j0, j1
        complex(dp), allocatable :: keff(:, :), sc(:, :), g(:, :), w(:, :), x(:, :), gx(:, :)
        complex(dp), allocatable :: hd(:, :), sd(:, :), hu(:, :), su(:, :)
        integer, allocatable :: ipiv(:)
        complex(dp) :: z
        real(dp) :: dev
        call set_seed(77)
        z = cmplx(0.3_dp, 0.05_dp, dp)
        call central_init(sys, bs)
        do m = 1, 4
            hd = rmat(bs(m), bs(m)); hd = 0.5_dp * (hd + conjg(transpose(hd)))
            sd = 0.05_dp * rmat(bs(m), bs(m)); sd = 0.5_dp * (sd + conjg(transpose(sd)))
            do i = 1, bs(m)
                sd(i, i) = sd(i, i) + 1.0_dp
            end do
            if (m < 4) then
                hu = 0.4_dp * rmat(bs(m), bs(m + 1)); su = 0.03_dp * rmat(bs(m), bs(m + 1))
                call assemble_central(sys, m, z, hd, sd, 0.1_dp * m, hu=hu, su=su, vnext=0.1_dp * (m + 1))
            else
                call assemble_central(sys, m, z, hd, sd, 0.4_dp)
            end if
        end do
        ! meios com g_X arbitrária (só a álgebra é testada aqui)
        gx = 0.3_dp * rmat(2, 2)
        call attach_lead(sys, 'bot', z, 0.3_dp * rmat(2, 3), 0.02_dp * rmat(2, 3), 0.1_dp, 0.0_dp, gx)
        gx = 0.3_dp * rmat(3, 3)
        call attach_lead(sys, 'top', z, 0.3_dp * rmat(3, 3), 0.02_dp * rmat(3, 3), 0.4_dp, 0.2_dp, gx)
        call rgf_solve(sys, res, ierr, errmsg)
        call check(ierr == 0, 'rgf_solve')
        call dense_central(sys, keff, sc)
        nt = sys%ntot
        allocate(g(nt, nt))
        call lu_factor(keff, ipiv, ierr, errmsg)
        call lu_solve_identity(keff, ipiv, g, ierr, errmsg)
        dev = 0.0_dp
        do m = 1, 4
            i0 = sys%boff(m) + 1; i1 = sys%boff(m + 1)
            dev = max(dev, maxval(abs(res%gd(m)%a - g(i0:i1, i0:i1))))
            if (m < 4) then
                j0 = sys%boff(m + 1) + 1; j1 = sys%boff(m + 2)
                dev = max(dev, maxval(abs(res%gu(m)%a - g(i0:i1, j0:j1))))
                dev = max(dev, maxval(abs(res%gl(m)%a - g(j0:j1, i0:i1))))
            end if
        end do
        call check_close(dev, 0.0_dp, 1.0e-12_dp, 'RGF blocks = dense inverse')
        ! G_{1,A} = -G_{1,:} K_{:,A} g_A envolve só o bloco 1 (K_{C,A} só acopla ao bloco 1)
        block
            complex(dp), allocatable :: ref(:, :), tmp(:, :)
            allocate(tmp(3, 2), ref(3, 2))
            call zmm('N', 'N', sys%kcx_bot%a, sys%gx_bot%a, tmp)
            call zmm('N', 'N', g(1:3, 1:3), tmp, ref, alpha=(-1.0_dp, 0.0_dp))
            call check_close(maxval(abs(res%glead_bot%a - ref)), 0.0_dp, 1.0e-12_dp, 'G_{1,A} boundary block')
        end block
        ! K x = w
        allocate(w(nt, 3), x(nt, 3))
        do j = 1, 3
            do i = 1, nt
                w(i, j) = cmplx(rnd(), rnd(), dp)
            end do
        end do
        x = w
        call solve_projection(sys, res, x, ierr, errmsg)
        call check(ierr == 0, 'solve_projection')
        call zmm('N', 'N', g, w, keff(:, 1:3))
        call check_close(maxval(abs(x - keff(:, 1:3))), 0.0_dp, 1.0e-12_dp, 'K x = w by block elimination')
    end subroutine rgf_vs_dense

    function rmat(a, b) result(m)
        integer, intent(in) :: a, b
        complex(dp) :: m(a, b)
        integer :: i, j
        do j = 1, b
            do i = 1, a
                m(i, j) = cmplx(rnd(), rnd(), dp)
            end do
        end do
    end function rmat

    subroutine aa_test()
        type(central_t) :: sys
        type(rgf_t) :: res
        real(dp), allocatable :: d(:)
        real(dp) :: rho, rhod
        integer :: ie, m
        do ie = 1, ne
            call chain_central(3, cmplx(e_ref(ie), eta, dp), 0.0_dp, 'both', sys)
            call rgf_solve(sys, res, ierr, errmsg)
            allocate(d(sys%ntot))
            call orbital_dos_mulliken(sys, res, d)
            do m = 1, 3
                rho = sum(d(sys%boff(m) + 1:sys%boff(m + 1)))
                call check_close(rho, rho_ref(ie), 1.0e-8_dp, 'A/A: layer DOS = bulk DOS')
            end do
            ! só o termo diagonal da camada do meio: erro de ~1% (o plano v4 mede 2%)
            rhod = 0.0_dp
            do m = 1, n
                rhod = rhod - aimag(sum(res%gd(2)%a(m, :) * sys%sd(2)%a(:, m))) / PI
            end do
            call check_close(rhod, rho_diag(ie), 1.0e-8_dp, 'A/A: diagonal-only term = Python value')
            call check(abs(rhod - rho_ref(ie)) > 1.0e-3_dp * rho_ref(ie), 'A/A: diagonal-only term is not the DOS')
            call check(min_neg_im_gdiag(res) >= 0.0_dp, 'A/A: -Im G_ii >= 0')
            deallocate(d)
        end do
    end subroutine aa_test

    subroutine thick_slab()
        type(central_t) :: sys
        type(rgf_t) :: res
        complex(dp), allocatable :: g(:, :), k00(:, :), k01(:, :), k10(:, :)
        complex(dp) :: z
        integer :: it
        real(dp) :: r
        z = cmplx(0.2_dp, eta, dp)
        call chain_central(400, z, 0.0_dp, 'none', sys)
        call rgf_solve(sys, res, ierr, errmsg)
        call check(ierr == 0, 'thick slab RGF')
        k00 = z * s00 - h00
        k01 = z * s01 - h01
        k10 = z * conjg(transpose(s01)) - conjg(transpose(h01))
        allocate(g(n, n))
        call surface_gf(k00, k01, k10, SIDE_BOTTOM, 1.0e-14_dp, 300, g, it, r, ierr, errmsg)
        call check_close(maxval(abs(res%gd(1)%a - g)), 0.0_dp, 1.0e-10_dp, 'thick slab G_11 = semi-infinite (bottom)')
        call surface_gf(k00, k01, k10, SIDE_TOP, 1.0e-14_dp, 300, g, it, r, ierr, errmsg)
        call check_close(maxval(abs(res%gd(400)%a - g)), 0.0_dp, 1.0e-10_dp, 'thick slab G_NN = semi-infinite (top)')
    end subroutine thick_slab

    subroutine decoupled_film_and_coherent()
        type(central_t) :: sys, sys0
        type(rgf_t) :: res
        complex(dp), allocatable :: keff(:, :), sc(:, :), hs(:, :), c(:, :), w(:, :)
        real(dp), allocatable :: d(:), ev(:), icoh(:), iinc(:), rd(:)
        type(orbital_groups_t) :: groups
        real(dp) :: energies(3) = [-1.2_dp, 0.05_dp, 1.7_dp], e, ref, dev, devc, minc
        complex(dp) :: z, acc, amp
        integer :: ie, nt, j, i, iq
        integer, parameter :: nl = 3, nq = 4

        nt = nl * n
        ! H e S densos do slab finito de 3 PLs: K(z=0) = -H, S
        call chain_central(nl, (0.0_dp, 0.0_dp), 0.0_dp, 'none', sys0)
        call dense_central(sys0, keff, sc)
        hs = -keff
        allocate(ev(nt), c(nt, nt), d(nt), icoh(nq), iinc(nq), w(nt, nq), rd(nl))
        call slab_eigen(hs, sc, ev, c)
        groups%ngroup = nl
        groups%norb = nt
        groups%start = [(1 + j * n, j = 0, nl)]
        groups%orb = [(j, j = 1, nt)]
        call set_seed(5)
        do iq = 1, nq
            do i = 1, nt
                w(i, iq) = cmplx(rnd(), rnd(), dp)
            end do
        end do
        dev = 0.0_dp; devc = 0.0_dp; minc = huge(1.0_dp)
        do ie = 1, size(energies)
            e = energies(ie)
            z = cmplx(e, eta, dp)
            ! filme sobre um substrato semi-infinito com contato nulo
            call chain_central(nl, z, 0.0_dp, 'bot', sys, coupled=.false.)
            call rgf_solve(sys, res, ierr, errmsg)
            call orbital_dos_mulliken(sys, res, d)
            ref = 0.0_dp
            do j = 1, nt
                ref = ref + (eta / PI) / ((e - ev(j))**2 + eta**2)
            end do
            dev = max(dev, abs(sum(d) - ref) / ref)
            call region_dos(d, groups, rd)
            call check_close(sum(rd), sum(d), 1.0e-12_dp, 'sum of layer regions = total')
            ! mapa coerente: polos nos autovalores do slab para qualquer w
            call coherent_map(sys, res, w, icoh, iinc, ierr, errmsg)
            do iq = 1, nq
                acc = (0.0_dp, 0.0_dp)
                do j = 1, nt
                    amp = dot_product(c(:, j), w(:, iq))      ! c_j^dagger w
                    acc = acc + abs(amp)**2 / (z - ev(j))
                end do
                devc = max(devc, abs(icoh(iq) + aimag(acc) / PI) / max(abs(aimag(acc)) / PI, 1.0e-300_dp))
                minc = min(minc, icoh(iq))
            end do
        end do
        call check_close(dev, 0.0_dp, 1.0e-10_dp, 'decoupled film: DOS = finite-slab Lorentzians')
        call check_close(devc, 0.0_dp, 1.0e-10_dp, 'coherent map: poles at slab eigenvalues for any w')
        call check(minc >= 0.0_dp, 'coherent map I_coh >= 0')
    end subroutine decoupled_film_and_coherent

    subroutine rigid_shift()
        type(central_t) :: sys
        type(rgf_t) :: res
        real(dp), allocatable :: d0(:), d1(:)
        real(dp), parameter :: delta = 0.37_dp, e = 0.2_dp
        call chain_central(2, cmplx(e - delta, eta, dp), 0.0_dp, 'both', sys)
        call rgf_solve(sys, res, ierr, errmsg)
        allocate(d0(sys%ntot), d1(sys%ntot))
        call orbital_dos_mulliken(sys, res, d0)
        call chain_central(2, cmplx(e, eta, dp), delta, 'both', sys)
        call rgf_solve(sys, res, ierr, errmsg)
        call orbital_dos_mulliken(sys, res, d1)
        call check_close(maxval(abs(d1 - d0)), 0.0_dp, 1.0e-11_dp, 'rigid shift H + Delta S moves spectrum by Delta')
    end subroutine rigid_shift

    subroutine lowdin_trace()
        type(central_t) :: sys
        type(rgf_t) :: res
        real(dp), allocatable :: dm(:), dl(:)
        real(dp) :: smin
        call chain_central(3, cmplx(0.4_dp, eta, dp), 0.0_dp, 'none', sys)
        call rgf_solve(sys, res, ierr, errmsg)
        allocate(dm(sys%ntot), dl(sys%ntot))
        call orbital_dos_mulliken(sys, res, dm)
        call orbital_dos_lowdin(sys, dl, smin, ierr, errmsg)
        call check(ierr == 0 .and. smin > 0.0_dp, 'Lowdin projection')
        call check_close(sum(dl), sum(dm), 1.0e-11_dp, 'Lowdin and Mulliken totals agree (isolated slab)')
        call check(maxval(abs(dl - dm)) > 1.0e-6_dp, 'Lowdin and Mulliken orbital weights differ')
        call check(minval(dl) >= 0.0_dp, 'Lowdin orbital DOS is non-negative')
    end subroutine lowdin_trace

end program test_central_green
