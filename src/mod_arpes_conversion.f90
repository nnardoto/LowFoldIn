! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Conversão de emissão normal entre energia de fóton e k_perp no modelo de
!> estado final de elétron quase livre. V0 é positivo quando o fundo da banda
!> interna está abaixo do vácuo; a energia inicial é fornecida relativa a E_F.
module arpes_conversion
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use types, only: dp
    implicit none
    private

    public :: HBAR2_OVER_2M_EV_ANG2, photon_to_kperp, kperp_to_photon

    real(dp), parameter :: HBAR2_OVER_2M_EV_ANG2 = 3.809982116154859_dp

contains

    subroutine photon_to_kperp(photon_ev, energy_relative_fermi_ev, work_function_ev, &
                               inner_potential_ev, kperp_inv_ang, kinetic_energy_ev, ierr, errmsg)
        real(dp), intent(in) :: photon_ev, energy_relative_fermi_ev
        real(dp), intent(in) :: work_function_ev, inner_potential_ev
        real(dp), intent(out) :: kperp_inv_ang, kinetic_energy_ev
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        real(dp) :: final_energy

        kperp_inv_ang = 0.0_dp
        kinetic_energy_ev = 0.0_dp
        ierr = 0
        errmsg = ''
        if (.not. all(ieee_is_finite([photon_ev, energy_relative_fermi_ev, &
                                      work_function_ev, inner_potential_ev]))) then
            ierr = 1; errmsg = 'ARPES conversion parameters must be finite'; return
        end if
        if (photon_ev <= 0.0_dp .or. work_function_ev < 0.0_dp .or. inner_potential_ev < 0.0_dp) then
            ierr = 1
            errmsg = 'require photon_energy > 0, work_function >= 0 and inner_potential >= 0'
            return
        end if
        kinetic_energy_ev = photon_ev - work_function_ev + energy_relative_fermi_ev
        if (.not. ieee_is_finite(kinetic_energy_ev)) then
            ierr = 1; errmsg = 'non-finite photoelectron kinetic energy'; return
        end if
        if (kinetic_energy_ev < 0.0_dp) then
            ierr = 2; errmsg = 'initial state is below the photoemission threshold'; return
        end if
        final_energy = kinetic_energy_ev + inner_potential_ev
        if (final_energy < 0.0_dp) then
            ierr = 2; errmsg = 'no real normal-emission final state'; return
        end if
        kperp_inv_ang = sqrt(final_energy / HBAR2_OVER_2M_EV_ANG2)
    end subroutine photon_to_kperp

    subroutine kperp_to_photon(kperp_inv_ang, energy_relative_fermi_ev, work_function_ev, &
                               inner_potential_ev, photon_ev, kinetic_energy_ev, ierr, errmsg)
        real(dp), intent(in) :: kperp_inv_ang, energy_relative_fermi_ev
        real(dp), intent(in) :: work_function_ev, inner_potential_ev
        real(dp), intent(out) :: photon_ev, kinetic_energy_ev
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg

        photon_ev = 0.0_dp
        kinetic_energy_ev = 0.0_dp
        ierr = 0
        errmsg = ''
        if (.not. all(ieee_is_finite([kperp_inv_ang, energy_relative_fermi_ev, &
                                      work_function_ev, inner_potential_ev]))) then
            ierr = 1; errmsg = 'ARPES conversion parameters must be finite'; return
        end if
        if (kperp_inv_ang < 0.0_dp .or. work_function_ev < 0.0_dp .or. inner_potential_ev < 0.0_dp) then
            ierr = 1
            errmsg = 'require k_perp >= 0, work_function >= 0 and inner_potential >= 0'
            return
        end if
        kinetic_energy_ev = HBAR2_OVER_2M_EV_ANG2 * kperp_inv_ang**2 - inner_potential_ev
        if (.not. ieee_is_finite(kinetic_energy_ev)) then
            ierr = 1; errmsg = 'non-finite photoelectron kinetic energy'; return
        end if
        if (kinetic_energy_ev < 0.0_dp) then
            ierr = 2; errmsg = 'k_perp corresponds to a negative vacuum kinetic energy'; return
        end if
        photon_ev = kinetic_energy_ev + work_function_ev - energy_relative_fermi_ev
        if (.not. ieee_is_finite(photon_ev) .or. photon_ev <= 0.0_dp) then
            ierr = 2; errmsg = 'k_perp and initial energy do not give a physical photon energy'; return
        end if
    end subroutine kperp_to_photon

end module arpes_conversion
