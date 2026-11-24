! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

program test_arpes_conversion
    use types, only: dp
    use arpes_conversion, only: HBAR2_OVER_2M_EV_ANG2, photon_to_kperp, kperp_to_photon
    use testing, only: check, check_close, test_summary
    implicit none

    real(dp) :: kperp, photon, kinetic
    integer :: ierr
    character(len=:), allocatable :: errmsg

    call photon_to_kperp(21.2_dp, 0.0_dp, 4.5_dp, 10.0_dp, kperp, kinetic, ierr, errmsg)
    call check(ierr == 0, 'He-I conversion at E_F: ' // errmsg)
    call check_close(kperp, 2.6472451549338456_dp, 2.0e-15_dp, 'He-I k_perp at E_F')
    call check_close(kinetic, 16.7_dp, 1.0e-14_dp, 'He-I kinetic energy at E_F')

    call photon_to_kperp(21.2_dp, -0.5_dp, 4.5_dp, 10.0_dp, kperp, kinetic, ierr, errmsg)
    call check(ierr == 0, 'He-I conversion below E_F: ' // errmsg)
    call check_close(kperp, 2.6223410721179965_dp, 2.0e-15_dp, 'bound initial state increases binding, not kinetic energy')
    call kperp_to_photon(kperp, -0.5_dp, 4.5_dp, 10.0_dp, photon, kinetic, ierr, errmsg)
    call check(ierr == 0, 'inverse conversion: ' // errmsg)
    call check_close(photon, 21.2_dp, 2.0e-14_dp, 'photon/k_perp round trip')
    call check_close(kinetic, 16.2_dp, 1.0e-14_dp, 'inverse kinetic energy')

    call kperp_to_photon(1.0_dp, 0.0_dp, 4.5_dp, 0.0_dp, photon, kinetic, ierr, errmsg)
    call check(ierr == 0, 'zero-inner-potential limit: ' // errmsg)
    call check_close(kinetic, HBAR2_OVER_2M_EV_ANG2, 0.0_dp, 'k=1 analytic kinetic energy')
    call check_close(photon, 8.309982116154859_dp, 2.0e-15_dp, 'k=1 analytic photon energy')

    call photon_to_kperp(4.0_dp, 0.0_dp, 4.5_dp, 10.0_dp, kperp, kinetic, ierr, errmsg)
    call check(ierr /= 0 .and. index(errmsg, 'threshold') > 0, 'reject energy below photoemission threshold')
    call kperp_to_photon(0.0_dp, 0.0_dp, 4.5_dp, 10.0_dp, photon, kinetic, ierr, errmsg)
    call check(ierr /= 0 .and. index(errmsg, 'negative vacuum') > 0, 'reject negative vacuum kinetic energy')
    call photon_to_kperp(21.2_dp, 0.0_dp, 4.5_dp, -1.0_dp, kperp, kinetic, ierr, errmsg)
    call check(ierr /= 0 .and. index(errmsg, 'inner_potential') > 0, 'reject ambiguous negative V0 convention')
    call photon_to_kperp(-21.2_dp, 0.0_dp, 4.5_dp, 10.0_dp, kperp, kinetic, ierr, errmsg)
    call check(ierr /= 0 .and. index(errmsg, 'photon_energy') > 0, 'reject negative photon energy')
    call kperp_to_photon(-1.0_dp, 0.0_dp, 4.5_dp, 10.0_dp, photon, kinetic, ierr, errmsg)
    call check(ierr /= 0 .and. index(errmsg, 'k_perp') > 0, 'reject negative outward k_perp')

    call test_summary('test_arpes_conversion')
end program test_arpes_conversion
