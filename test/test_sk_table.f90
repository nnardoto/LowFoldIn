! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Tabela de Slater–Koster (sk_table):
!>  - blocos iguais à referência independente (rotação de harmônicos reais,
!>    test/data/green/gen_sk_reference.py) a 1e-12;
!>  - paridade E_ab(n) = (-1)^(la+lb) E_ba(n)^T e E_ab(-n) = (-1)^(la+lb) E_ab(n);
!>  - ligação ao longo de z: só sigma/pi/delta na diagonal de |m|;
!>  - leitura de bases do OpenMX ('Te s3p2d2f1').
program test_sk_table
    use types, only: dp
    use testing, only: check, check_close, test_summary
    use sk_table, only: sk_block, parse_basis, species_basis_t, basis_shell_index
    implicit none

    call reference_blocks()
    call parity()
    call bond_along_z()
    call basis_strings()
    call test_summary('test_sk_table')

contains

    subroutine reference_blocks()
        integer :: u, ios, nc, c, la, lb, a
        real(dp) :: n(3), v(3), ref(5, 5), dmax
        real(dp), allocatable :: e(:, :)
        open(newunit=u, file='test/data/green/sk_reference.txt', status='old', action='read', iostat=ios)
        call check(ios == 0, 'open sk_reference.txt')
        if (ios /= 0) return
        read(u, *) nc
        dmax = 0.0_dp
        do c = 1, nc
            read(u, *) la, lb, n, v
            do a = 1, 2 * la + 1
                read(u, *) ref(a, 1:2 * lb + 1)
            end do
            e = sk_block(la, lb, n, v)
            dmax = max(dmax, maxval(abs(e - ref(1:2 * la + 1, 1:2 * lb + 1))))
        end do
        close(u)
        call check(nc == 72, 'reference has 72 blocks')
        call check_close(dmax, 0.0_dp, 1.0e-12_dp, 'SK table = rotation reference (s, p, d)')
    end subroutine reference_blocks

    subroutine parity()
        integer :: la, lb
        real(dp) :: n(3), v(3), dmax, dinv
        n = [0.31_dp, -0.72_dp, 0.45_dp]
        v = [1.1_dp, -0.37_dp, 0.21_dp]
        dmax = 0.0_dp
        dinv = 0.0_dp
        do la = 0, 2
            do lb = 0, 2
                dmax = max(dmax, maxval(abs(sk_block(la, lb, n, v) - &
                           real((-1)**(la + lb), dp) * transpose(sk_block(lb, la, n, v)))))
                dinv = max(dinv, maxval(abs(sk_block(la, lb, -n, v) - &
                           real((-1)**(la + lb), dp) * sk_block(la, lb, n, v))))
            end do
        end do
        call check_close(dmax, 0.0_dp, 1.0e-14_dp, 'E_ab(n) = (-1)^(la+lb) E_ba(n)^T')
        call check_close(dinv, 0.0_dp, 1.0e-14_dp, 'E_ab(-n) = (-1)^(la+lb) E_ab(n)')
    end subroutine parity

    subroutine bond_along_z()
        real(dp) :: v(3), e(5, 5), p(3, 3)
        real(dp), allocatable :: x(:, :)
        v = [0.9_dp, -0.4_dp, 0.15_dp]
        e = sk_block(2, 2, [0.0_dp, 0.0_dp, 2.0_dp], v)
        ! 3z2-r2 sigma; x2-y2, xy delta; xz, yz pi
        call check_close(e(1, 1), v(1), 1.0e-14_dp, 'dd along z: 3z2-r2 = sigma')
        call check_close(e(2, 2), v(3), 1.0e-14_dp, 'dd along z: x2-y2 = delta')
        call check_close(e(3, 3), v(3), 1.0e-14_dp, 'dd along z: xy = delta')
        call check_close(e(4, 4), v(2), 1.0e-14_dp, 'dd along z: xz = pi')
        call check_close(e(5, 5), v(2), 1.0e-14_dp, 'dd along z: yz = pi')
        call check_close(sum(abs(e)) - (abs(e(1, 1)) + abs(e(2, 2)) + abs(e(3, 3)) + abs(e(4, 4)) + &
                         abs(e(5, 5))), 0.0_dp, 1.0e-14_dp, 'dd along z: diagonal')
        p = sk_block(1, 1, [0.0_dp, 0.0_dp, 1.0_dp], v)
        call check_close(p(3, 3), v(1), 1.0e-14_dp, 'pp along z: pz = sigma')
        call check_close(p(1, 1), v(2), 1.0e-14_dp, 'pp along z: px = pi')
        x = sk_block(0, 1, [0.0_dp, 0.0_dp, 1.0_dp], v)
        call check_close(x(1, 3), v(1), 1.0e-14_dp, 'E(s, pz) along +z = V_sp_sigma')
        x = sk_block(1, 0, [0.0_dp, 0.0_dp, 1.0_dp], v)
        call check_close(x(3, 1), -v(1), 1.0e-14_dp, 'E(pz, s) along +z = -V_sp_sigma')
    end subroutine bond_along_z

    subroutine basis_strings()
        type(species_basis_t) :: b
        integer :: ierr
        character(len=:), allocatable :: msg
        call parse_basis('Te s3p2d2f1', b, ierr, msg)
        call check(ierr == 0, 'parse Te s3p2d2f1')
        call check(b%nspatial == 3 + 6 + 10 + 7, 's3p2d2f1 has 26 spatial orbitals')
        call check(b%nshell == 8, 's3p2d2f1 has 8 shells')
        call check(b%shell(basis_shell_index(b, 1, 2))%offset == 6, 'p zeta 2 starts at offset 6')
        call check(b%shell(basis_shell_index(b, 2, 1))%offset == 9, 'd zeta 1 starts at offset 9')
        call check(basis_shell_index(b, 2, 3) == 0, 'd zeta 3 absent')
        call parse_basis('Sb s2p2d1', b, ierr, msg)
        call check(ierr == 0 .and. b%nspatial == 2 + 6 + 5, 'parse Sb s2p2d1')
        call parse_basis('Sb s2x1', b, ierr, msg)
        call check(ierr /= 0, 'invalid letter is an error')
        call parse_basis('Sb', b, ierr, msg)
        call check(ierr /= 0, 'missing basis string is an error')
    end subroutine basis_strings

end program test_sk_table
