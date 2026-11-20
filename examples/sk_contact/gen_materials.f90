! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Dois materiais SINTÉTICOS (não DFT) para o exemplo do contato SK:
!>   A: 2 átomos (Ax, Ay, base s1p1), a_s inclinado, 8 orbitais por PL;
!>   B: 1 átomo (Bz, base s1p1d1), a_s ao longo de z, 9 orbitais por PL.
!> Mesma rede hexagonal no plano (a = 3.2 Å). O acoplamento entre camadas de
!> cada bulk é SK com decaimento exponencial; S não ortogonal com S01 pequeno.
!> As PLs têm tamanhos diferentes: o contato 'scaled' não se aplica, o 'sk' sim.
program gen_materials
    use types, only: dp
    use HS_io, only: HSDataFormat
    use sk_table, only: sk_block, species_basis_t, parse_basis
    implicit none

    real(dp), parameter :: A1(3) = [3.2_dp, 0.0_dp, 0.0_dp], A2(3) = [-1.6_dp, 2.7712812921102_dp, 0.0_dp]
    real(dp), parameter :: RC = 5.5_dp
    integer :: st = 7

    call material('A', [A1, A2, [0.35_dp, 0.20_dp, 4.4_dp]], 2, &
                  reshape([0.0_dp, 0.0_dp, 0.3_dp, 1.6_dp, 0.92_dp, 1.9_dp], [3, 2]), &
                  ['Ax s1p1  ', 'Ay s1p1  '], [-0.6_dp, 0.4_dp], 1.6_dp, -0.25_dp)
    call material('B', [A1, A2, [0.0_dp, 0.0_dp, 3.7_dp]], 1, &
                  reshape([0.0_dp, 0.0_dp, 0.5_dp], [3, 1]), &
                  ['Bz s1p1d1'], [0.1_dp], 1.3_dp, 0.15_dp)
    print '(A)', 'Wrote A_H.hsio, A_S.hsio, A.nml, B_H.hsio, B_S.hsio, B.nml (synthetic, NOT DFT)'

contains

    real(dp) function rnd()
        st = modulo(1103515245 * st + 12345, 2147483647)
        rnd = 2.0_dp * real(st, dp) / 2147483647.0_dp - 1.0_dp
    end function rnd

    !> lat = [a1, a2, a3]; nat átomos em tau; bases; on-site; beta do
    !> acoplamento; ef (E_F gravado nos metadados, escala absoluta).
    subroutine material(name, latv, nat, tau, btxt, eon, beta, ef)
        character(len=*), intent(in) :: name
        integer, intent(in) :: nat
        real(dp), intent(in) :: latv(9), tau(3, nat), eon(nat), beta, ef
        character(len=*), intent(in) :: btxt(nat)
        type(species_basis_t) :: b(nat)
        type(HSDataFormat) :: hh, hs
        real(dp) :: lat(3, 3), d(3), r(3), vs(3, 0:2, 0:2), vh(3, 0:2, 0:2), vo(3, 0:2, 0:2)
        real(dp), allocatable :: bh(:, :), bs(:, :)
        integer, allocatable :: off(:)
        integer :: i, j, n, nr, ir, n1, n2, n3, ns, ierr, ir0, u, la, lb, k
        character(len=:), allocatable :: msg

        lat = transpose(reshape(latv, [3, 3]))
        allocate(off(nat))
        n = 0
        do i = 1, nat
            call parse_basis(btxt(i), b(i), ierr, msg)
            off(i) = n
            n = n + b(i)%nspatial
        end do
        ! parâmetros por (la, lb): entre PLs (vh), dentro da PL (vs), overlap (vo)
        do la = 0, 2
            do lb = la, 2
                do k = 1, 3
                    vh(k, la, lb) = merge(0.9_dp * rnd(), 0.0_dp, k <= la + 1)
                    vs(k, la, lb) = merge(0.7_dp * rnd(), 0.0_dp, k <= la + 1)
                    vo(k, la, lb) = merge(0.03_dp * rnd(), 0.0_dp, k <= la + 1)
                end do
                vh(:, lb, la) = vh(:, la, lb); vs(:, lb, la) = vs(:, la, lb); vo(:, lb, la) = vo(:, la, lb)
            end do
        end do
        ns = 4
        nr = (2 * ns + 1)**2 * 3
        hh%nAtoms = nat; hh%NOrbs = n; hh%Neighbors = nr; hh%Title = 'synthetic SK material ' // name
        allocate(hh%Orbs(nat), hh%IndexList(3, nr), hh%Degenerecency(nr), hh%Data(n, n, nr))
        hh%Orbs = [(b(i)%nspatial, i = 1, nat)]
        hh%Degenerecency = 1
        hh%Data = (0.0_dp, 0.0_dp)
        ir = 0
        do n3 = -1, 1
            do n2 = -ns, ns
                do n1 = -ns, ns
                    ir = ir + 1
                    hh%IndexList(:, ir) = [n1, n2, n3]
                    if (all([n1, n2, n3] == 0)) ir0 = ir
                end do
            end do
        end do
        hs = hh
        do ir = 1, nr
            r = matmul(real(hh%IndexList(:, ir), dp), lat)
            n3 = hh%IndexList(3, ir)
            do i = 1, nat
                do j = 1, nat
                    d = tau(:, j) + r - tau(:, i)
                    if (norm2(d) > RC .or. norm2(d) < 1.0e-8_dp) cycle
                    if (n3 == 0) then
                        call blk(b(i), b(j), d, vs, 0.9_dp, bh)
                        call blk(b(i), b(j), d, vo, 1.0_dp, bs)
                    else
                        call blk(b(i), b(j), d, vh, beta, bh)
                        call blk(b(i), b(j), d, vo, 1.0_dp, bs)
                        bs = 0.3_dp * bs
                    end if
                    hh%Data(off(i) + 1:off(i) + b(i)%nspatial, off(j) + 1:off(j) + b(j)%nspatial, ir) = &
                        cmplx(bh, 0.0_dp, dp)
                    hs%Data(off(i) + 1:off(i) + b(i)%nspatial, off(j) + 1:off(j) + b(j)%nspatial, ir) = &
                        cmplx(bs, 0.0_dp, dp)
                end do
            end do
        end do
        do i = 1, nat
            do k = 1, b(i)%nspatial
                hh%Data(off(i) + k, off(i) + k, ir0) = hh%Data(off(i) + k, off(i) + k, ir0) + &
                    cmplx(ef + eon(i) + 0.5_dp * real(k - 1, dp) / real(b(i)%nspatial, dp), 0.0_dp, dp)
                hs%Data(off(i) + k, off(i) + k, ir0) = hs%Data(off(i) + k, off(i) + k, ir0) + (1.0_dp, 0.0_dp)
            end do
        end do
        ! H(R) acima está na escala absoluta: E_F = ef, sem deslocamento
        hh%has_meta = .true.; hh%fermi_ev = ef; hh%fermi_shift_applied = .false.
        hs%has_meta = .true.; hs%fermi_ev = ef; hs%fermi_shift_applied = .false.
        call hh%WriteFile(name // '_H.hsio')
        call hs%WriteFile(name // '_S.hsio')
        open(newunit=u, file=name // '.nml', status='replace', action='write')
        write(u, '(A)') '! synthetic material ' // name // ' (NOT DFT)'
        write(u, '(A)') '&STRUCTURE'
        do i = 1, 3
            write(u, '(A,I0,A,3(F20.14,:,","))') '  LatticeVectors(', i, ',:) = ', lat(i, :)
        end do
        write(u, '(A)') '/'
        write(u, '(A)') '&ATOMS'
        write(u, '(I0)') nat
        do i = 1, nat
            write(u, '(A,3F20.14)') trim(b(i)%name), tau(:, i)
        end do
        write(u, '(A)') '/'
        close(u)
    end subroutine material

    !> Bloco SK espacial entre os átomos de bases bi, bj à distância d.
    subroutine blk(bi, bj, d, v, bt, o)
        type(species_basis_t), intent(in) :: bi, bj
        real(dp), intent(in) :: d(3), v(3, 0:2, 0:2), bt
        real(dp), allocatable, intent(out) :: o(:, :)
        integer :: si, sj, li, lj
        allocate(o(bi%nspatial, bj%nspatial))
        o = 0.0_dp
        do si = 1, bi%nshell
            do sj = 1, bj%nshell
                li = bi%shell(si)%l; lj = bj%shell(sj)%l
                o(bi%shell(si)%offset + 1:bi%shell(si)%offset + 2 * li + 1, &
                  bj%shell(sj)%offset + 1:bj%shell(sj)%offset + 2 * lj + 1) = &
                    sk_block(li, lj, d, v(:, li, lj) * exp(-bt * (norm2(d) - 3.0_dp)))
            end do
        end do
    end subroutine blk

end program gen_materials
