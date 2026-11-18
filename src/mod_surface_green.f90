! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Green retartada de uma superficie semi-infinita em base ortonormal.
!>
!> `build_principal_layer` agrupa celulas consecutivas na direcao normal ate
!> que a Hamiltoniana tenha apenas acoplamento entre camadas principais
!> vizinhas. `surface_greens_function` usa a decimacao iterativa de Lopez
!> Sancho. Os lados `left` e `right` trocam H01 por H01^dagger e permitem
!> representar as duas terminacoes de um mesmo bulk.
module surface_green
    use types, only: dp
    implicit none
    private

    public :: build_principal_layer, surface_greens_function
    public :: spectral_matrix, trace_surface_dos

    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            use types, only: dp
            integer     :: n, nrhs, lda, ldb, info
            integer     :: ipiv(*)
            complex(dp) :: a(lda, *), b(ldb, *)
        end subroutine zgesv
    end interface

contains

    subroutine build_principal_layer(hr, rpt, kfrac, normal_axis, h00, h01, thickness, ierr, errmsg)
        complex(dp), intent(in)  :: hr(:, :, :)
        integer,     intent(in)  :: rpt(:, :)
        real(dp),    intent(in)  :: kfrac(3)
        integer,     intent(in)  :: normal_axis
        complex(dp), allocatable, intent(out) :: h00(:, :), h01(:, :)
        integer,     intent(out) :: thickness, ierr
        character(len=:), allocatable, intent(out) :: errmsg

        real(dp), parameter :: TWO_PI = 6.283185307179586476925286766559_dp
        integer :: n, nr, ir, a, b, rz, row0, col0, p
        real(dp) :: phase_arg
        complex(dp) :: phase

        ierr = 0
        errmsg = ''
        thickness = 0
        if (normal_axis < 1 .or. normal_axis > 3) then
            ierr = 1; errmsg = 'build_principal_layer: normal_axis must be 1, 2 or 3'; return
        end if
        n = size(hr, 1)
        nr = size(hr, 3)
        if (size(hr, 2) /= n .or. size(rpt, 1) /= 3 .or. size(rpt, 2) /= nr) then
            ierr = 1; errmsg = 'build_principal_layer: inconsistent H(R)/R dimensions'; return
        end if

        thickness = max(1, maxval(abs(rpt(normal_axis, :))))
        p = n * thickness
        allocate(h00(p, p), h01(p, p))
        h00 = cmplx(0.0_dp, 0.0_dp, dp)
        h01 = cmplx(0.0_dp, 0.0_dp, dp)

        do ir = 1, nr
            rz = rpt(normal_axis, ir)
            phase_arg = TWO_PI * (dot_product(kfrac, real(rpt(:, ir), dp)) - &
                                   kfrac(normal_axis) * real(rz, dp))
            phase = exp(cmplx(0.0_dp, phase_arg, dp))
            do a = 0, thickness - 1
                b = a + rz
                if (b >= 0 .and. b < thickness) then
                    row0 = a * n
                    col0 = b * n
                    h00(row0 + 1:row0 + n, col0 + 1:col0 + n) = &
                        h00(row0 + 1:row0 + n, col0 + 1:col0 + n) + phase * hr(:, :, ir)
                else if (b >= thickness .and. b < 2 * thickness) then
                    row0 = a * n
                    col0 = (b - thickness) * n
                    h01(row0 + 1:row0 + n, col0 + 1:col0 + n) = &
                        h01(row0 + 1:row0 + n, col0 + 1:col0 + n) + phase * hr(:, :, ir)
                end if
            end do
        end do
    end subroutine build_principal_layer

    subroutine surface_greens_function(h00, h01, energy, eta, side, tol, max_iter, gs, niter, ierr, errmsg)
        complex(dp), intent(in)  :: h00(:, :), h01(:, :)
        real(dp),    intent(in)  :: energy, eta, tol
        character(len=*), intent(in) :: side
        integer,     intent(in)  :: max_iter
        complex(dp), intent(out) :: gs(:, :)
        integer,     intent(out) :: niter, ierr
        character(len=:), allocatable, intent(out) :: errmsg

        integer :: n, info
        logical :: converged
        complex(dp) :: z
        complex(dp), allocatable :: eps(:, :), eps_s(:, :), alpha(:, :), beta(:, :)
        complex(dp), allocatable :: g(:, :), agb(:, :), bga(:, :), next_a(:, :), next_b(:, :)

        ierr = 0
        errmsg = ''
        niter = 0
        n = size(h00, 1)
        if (size(h00, 2) /= n .or. any(shape(h01) /= [n, n]) .or. any(shape(gs) /= [n, n])) then
            ierr = 1; errmsg = 'surface_greens_function: inconsistent matrix dimensions'; return
        end if
        if (eta <= 0.0_dp .or. tol <= 0.0_dp .or. max_iter < 1) then
            ierr = 1; errmsg = 'surface_greens_function: eta, tol and max_iter must be positive'; return
        end if

        allocate(eps(n,n), eps_s(n,n), alpha(n,n), beta(n,n), g(n,n))
        allocate(agb(n,n), bga(n,n), next_a(n,n), next_b(n,n))
        eps = h00
        eps_s = h00
        select case (trim(side))
        case ('left')
            alpha = h01
            beta = conjg(transpose(h01))
        case ('right')
            alpha = conjg(transpose(h01))
            beta = h01
        case default
            ierr = 1; errmsg = "surface_greens_function: side must be 'left' or 'right'"; return
        end select
        z = cmplx(energy, eta, dp)
        converged = .false.

        do niter = 1, max_iter
            call resolvent(z, eps, g, info)
            if (info /= 0) then
                ierr = 2; errmsg = 'surface_greens_function: singular decimation resolvent'; return
            end if
            agb = matmul(alpha, matmul(g, beta))
            bga = matmul(beta, matmul(g, alpha))
            next_a = matmul(alpha, matmul(g, alpha))
            next_b = matmul(beta, matmul(g, beta))
            eps_s = eps_s + agb
            eps = eps + agb + bga
            alpha = next_a
            beta = next_b
            if (max(maxval(abs(alpha)), maxval(abs(beta))) < tol) then
                converged = .true.
                exit
            end if
        end do
        if (.not. converged) then
            ierr = 3; errmsg = 'surface_greens_function: decimation did not converge'; return
        end if

        call resolvent(z, eps_s, gs, info)
        if (info /= 0) then
            ierr = 2; errmsg = 'surface_greens_function: singular surface resolvent'; return
        end if

    contains

        subroutine resolvent(zin, h, inverse, lapack_info)
            complex(dp), intent(in)  :: zin, h(:, :)
            complex(dp), intent(out) :: inverse(:, :)
            integer,     intent(out) :: lapack_info
            complex(dp), allocatable :: a(:, :)
            integer, allocatable :: ipiv(:)
            integer :: k, nn
            nn = size(h, 1)
            allocate(a(nn, nn), ipiv(nn))
            a = -h
            do k = 1, nn
                a(k, k) = a(k, k) + zin
            end do
            inverse = cmplx(0.0_dp, 0.0_dp, dp)
            do k = 1, nn
                inverse(k, k) = cmplx(1.0_dp, 0.0_dp, dp)
            end do
            call zgesv(nn, nn, a, nn, ipiv, inverse, nn, lapack_info)
        end subroutine resolvent
    end subroutine surface_greens_function

    pure subroutine spectral_matrix(gs, a)
        complex(dp), intent(in)  :: gs(:, :)
        complex(dp), intent(out) :: a(:, :)
        real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp
        a = (cmplx(0.0_dp, 1.0_dp, dp) / (2.0_dp * PI)) * &
            (gs - conjg(transpose(gs)))
    end subroutine spectral_matrix

    pure real(dp) function trace_surface_dos(gs) result(dos)
        complex(dp), intent(in) :: gs(:, :)
        real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp
        integer :: i
        dos = 0.0_dp
        do i = 1, min(size(gs, 1), size(gs, 2))
            dos = dos - aimag(gs(i, i)) / PI
        end do
    end function trace_surface_dos

end module surface_green
