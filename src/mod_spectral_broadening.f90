! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Broadening of discrete (energy, weight) spectral data into continuous maps
!> A(k,E) on a uniform energy grid, plus Fermi-Dirac occupation
!> (docs/design/contracts.md §7).
!>
!> Conventions:
!>   - Lorentzian: width = HWHM eta,              L(x) = (1/pi) eta/(x^2+eta^2)
!>   - Gaussian:   width = standard deviation sigma (NOT HWHM),
!>                 G(x) = exp(-x^2/(2 sigma^2)) / (sigma sqrt(2 pi))
!>   - Both kernels integrate to 1 over the real line; weights are used with
!>     their sign (no abs, no clipping).
!> Library routines never stop. Callers must validate width > 0, ne >= 2 and
!> emax > emin beforehand (specmap does it at input validation).
module spectral_broadening
    use types, only: dp
    implicit none
    private

    integer, parameter, public :: KERNEL_LORENTZIAN = 1, KERNEL_GAUSSIAN = 2

    !> Boltzmann constant in eV/K.
    real(dp), parameter, public :: KB_EV_PER_K = 8.617333262e-5_dp

    !> Gaussian contributions beyond this many sigma are skipped
    !> (exp(-72) ~ 5e-32, far below double-precision relevance).
    real(dp), parameter :: GAUSS_CUTOFF_NSIGMA = 12.0_dp

    real(dp), parameter :: PI = acos(-1.0_dp)

    public :: lorentzian, gaussian, fermi_dirac, energy_grid, broaden

contains

    !> Normalised Lorentzian with half width at half maximum eta.
    elemental function lorentzian(x, eta) result(y)
        real(dp), intent(in) :: x, eta
        real(dp) :: y
        y = (eta / PI) / (x*x + eta*eta)
    end function lorentzian

    !> Normalised Gaussian with standard deviation sigma.
    elemental function gaussian(x, sigma) result(y)
        real(dp), intent(in) :: x, sigma
        real(dp) :: y
        y = exp(-x*x / (2.0_dp*sigma*sigma)) / (sigma * sqrt(2.0_dp*PI))
    end function gaussian

    !> Fermi-Dirac occupation f = 1/(1+exp((e-mu)/(kB T))), overflow-safe.
    !> T == 0 gives the step function with f(mu) = 1/2. T < 0 is not
    !> meaningful and must be rejected by the caller; it is treated as T = 0.
    elemental function fermi_dirac(e, mu, temperature_k) result(f)
        real(dp), intent(in) :: e, mu, temperature_k
        real(dp) :: f
        real(dp) :: z, ez

        if (temperature_k <= 0.0_dp) then
            if (e < mu) then
                f = 1.0_dp
            else if (e > mu) then
                f = 0.0_dp
            else
                f = 0.5_dp
            end if
            return
        end if

        z = (e - mu) / (KB_EV_PER_K * temperature_k)
        if (z > 700.0_dp) then
            f = 0.0_dp
        else if (z < -700.0_dp) then
            f = 1.0_dp
        else if (z > 0.0_dp) then
            ez = exp(-z)
            f = ez / (1.0_dp + ez)
        else
            f = 1.0_dp / (1.0_dp + exp(z))
        end if
    end function fermi_dirac

    !> Uniform grid E_i = emin + (i-1)(emax-emin)/(ne-1), i = 1..ne (ne >= 2).
    subroutine energy_grid(emin, emax, ne, grid)
        real(dp), intent(in)  :: emin, emax
        integer,  intent(in)  :: ne
        real(dp), intent(out) :: grid(:)
        integer :: i
        real(dp) :: de

        de = (emax - emin) / real(max(ne - 1, 1), dp)
        do i = 1, min(ne, size(grid))
            grid(i) = emin + real(i - 1, dp) * de
        end do
        if (ne >= 2 .and. size(grid) >= ne) grid(ne) = emax
    end subroutine energy_grid

    !> A(:,k) = sum_n weight(n,k) * K(grid - energy(n,k)).
    !> energy(nband,nk), weight(nband,nk), grid(ne) -> A(ne,nk).
    !> Lorentzian tails are always computed in full (no cutoff); Gaussian
    !> contributions with |x| > 12 sigma are skipped. Unknown kernel -> A = 0.
    subroutine broaden(energy, weight, grid, kernel, width, A)
        real(dp), intent(in)  :: energy(:,:), weight(:,:), grid(:)
        integer,  intent(in)  :: kernel
        real(dp), intent(in)  :: width
        real(dp), intent(out) :: A(:,:)

        integer  :: ik, ib, ie, nband, nk, ne
        real(dp) :: x, w, cut

        nband = size(energy, 1)
        nk    = size(energy, 2)
        ne    = size(grid)
        A = 0.0_dp
        cut = GAUSS_CUTOFF_NSIGMA * width

        do ik = 1, nk
            do ib = 1, nband
                w = weight(ib, ik)
                if (w == 0.0_dp) cycle
                select case (kernel)
                case (KERNEL_LORENTZIAN)
                    do ie = 1, ne
                        x = grid(ie) - energy(ib, ik)
                        A(ie, ik) = A(ie, ik) + w * lorentzian(x, width)
                    end do
                case (KERNEL_GAUSSIAN)
                    do ie = 1, ne
                        x = grid(ie) - energy(ib, ik)
                        if (abs(x) > cut) cycle
                        A(ie, ik) = A(ie, ik) + w * gaussian(x, width)
                    end do
                end select
            end do
        end do
    end subroutine broaden

end module spectral_broadening
