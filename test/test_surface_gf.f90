! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Portão 2 (fase 2): Green de superfície semi-infinita.
!>  - cadeia analítica (1 orbital) a 1e-10 nos dois lados;
!>  - cadeia não ortogonal complexa igual à referência Python a 1e-12;
!>  - resíduo da equação de Dyson pequeno e -Im g_ii >= 0;
!>  - L(z) /= U(z)^dagger para eta > 0 (o código não supõe hermiticidade de K).
program test_surface_gf
    use types, only: dp
    use surface_green, only: surface_gf, surface_gf_residual, SIDE_TOP, SIDE_BOTTOM
    use green_models, only: read_cmat
    use testing, only: check, check_close, test_summary
    implicit none

    call analytic_chain()
    call python_reference()
    call test_summary('test_surface_gf')

contains

    subroutine analytic_chain()
        real(dp), parameter :: eps = 0.3_dp, t = -1.0_dp, eta = 1.0e-3_dp
        real(dp) :: energies(4) = [-1.5_dp, 0.2_dp, 2.5_dp, -2.9_dp]
        complex(dp) :: k00(1, 1), k01(1, 1), k10(1, 1), gs(1, 1), z, g1, g2, gex
        integer :: ie, niter, ierr, side
        real(dp) :: res
        character(len=:), allocatable :: errmsg
        do ie = 1, size(energies)
            z = cmplx(energies(ie), eta, dp)
            k00 = z - eps
            k01 = -t
            k10 = -t
            g1 = ((z - eps) + sqrt((z - eps)**2 - 4.0_dp * t**2)) / (2.0_dp * t**2)
            g2 = ((z - eps) - sqrt((z - eps)**2 - 4.0_dp * t**2)) / (2.0_dp * t**2)
            gex = merge(g1, g2, abs(g1 * t) < 1.0_dp)
            do side = -1, 1, 2
                call surface_gf(k00, k01, k10, side, 1.0e-13_dp, 400, gs, niter, res, ierr, errmsg)
                call check(ierr == 0, 'analytic chain: decimation converged')
                call check_close(abs(gs(1, 1) - gex), 0.0_dp, 1.0e-10_dp, 'analytic chain g_s')
                call check(res < 1.0e-10_dp, 'analytic chain: Dyson residual')
                call check(-aimag(gs(1, 1)) >= 0.0_dp, 'analytic chain: -Im g >= 0')
            end do
        end do
    end subroutine analytic_chain

    subroutine python_reference()
        integer :: u, n, ne, ie, niter, ierr, i
        real(dp) :: eta, e, res
        complex(dp), allocatable :: h00(:, :), h01(:, :), s00(:, :), s01(:, :), gtop(:, :), gbot(:, :)
        complex(dp), allocatable :: k00(:, :), k01(:, :), k10(:, :), gs(:, :)
        complex(dp) :: z
        character(len=256) :: line
        character(len=:), allocatable :: errmsg

        open(newunit=u, file='test/data/green/nonorth_chain.txt', status='old', action='read')
        do
            read(u, '(A)') line
            if (line(1:1) /= '#') exit
        end do
        read(line, *) n, ne, eta
        allocate(h00(n, n), h01(n, n), s00(n, n), s01(n, n), gtop(n, n), gbot(n, n), gs(n, n))
        call read_cmat(u, h00)
        call read_cmat(u, h01)
        call read_cmat(u, s00)
        call read_cmat(u, s01)
        do ie = 1, ne
            read(u, *) e
            call read_cmat(u, gtop)
            call read_cmat(u, gbot)
            z = cmplx(e, eta, dp)
            k00 = z * s00 - h00
            k01 = z * s01 - h01
            k10 = z * conjg(transpose(s01)) - conjg(transpose(h01))
            call check(maxval(abs(k10 - conjg(transpose(k01)))) > 1.0e-3_dp, 'L(z) /= U(z)^dagger for eta > 0')

            call surface_gf(k00, k01, k10, SIDE_TOP, 1.0e-14_dp, 200, gs, niter, res, ierr, errmsg)
            call check(ierr == 0, 'nonorth chain top: converged')
            call check_close(maxval(abs(gs - gtop)), 0.0_dp, 1.0e-12_dp, 'nonorth chain top = Python reference')
            call check(res < 1.0e-12_dp, 'nonorth chain top: residual')
            call check(niter < 40, 'nonorth chain top: iteration count')
            do i = 1, n
                call check(-aimag(gs(i, i)) >= 0.0_dp, 'nonorth chain top: -Im g_ii >= 0')
            end do
            call check_close(surface_gf_residual(k00, k01, k10, SIDE_TOP, gs), res, 1.0e-15_dp, &
                             'residual routine consistent')

            call surface_gf(k00, k01, k10, SIDE_BOTTOM, 1.0e-14_dp, 200, gs, niter, res, ierr, errmsg)
            call check(ierr == 0, 'nonorth chain bottom: converged')
            call check_close(maxval(abs(gs - gbot)), 0.0_dp, 1.0e-12_dp, 'nonorth chain bottom = Python reference')
            call check(res < 1.0e-12_dp, 'nonorth chain bottom: residual')
            do i = 1, n
                call check(-aimag(gs(i, i)) >= 0.0_dp, 'nonorth chain bottom: -Im g_ii >= 0')
            end do
            ! A superfície oposta é de fato diferente (cadeia sem simetria de inversão)
            call check(maxval(abs(gtop - gbot)) > 1.0e-3_dp, 'top and bottom surfaces differ')
        end do
        close(u)

        ! Sem convergência: ierr = 3, mas g_s ainda devolvida
        z = cmplx(0.2_dp, 1.0e-6_dp, dp)
        k00 = z * s00 - h00
        k01 = z * s01 - h01
        k10 = z * conjg(transpose(s01)) - conjg(transpose(h01))
        call surface_gf(k00, k01, k10, SIDE_TOP, 1.0e-14_dp, 3, gs, niter, res, ierr, errmsg)
        call check(ierr == 3 .and. niter == 3, 'non-convergence reported with ierr = 3')
        call surface_gf(k00, k01, k10, 0, 1.0e-14_dp, 3, gs, niter, res, ierr, errmsg)
        call check(ierr == 1, 'invalid side rejected')
    end subroutine python_reference

end program test_surface_gf
