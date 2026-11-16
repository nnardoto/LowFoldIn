! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! tools/specmap.f90
!
! Broadens the discrete (energy, weight) data of a spectral data file into
! continuous maps A_P(k,E) for regions or +/- combinations of regions
! (docs/design/contracts.md §7). Reads ONLY the spectral data file: no HSIO,
! no structure, no diagonalisation. Changing kernel/width/window/occupation
! therefore only requires re-running specmap.
!
! Usage:
!   specmap <specmap.nml>
!
! Input : namelist &SPECMAP (options in docs/specmap.md) pointing to the
!         <prefix>.spectral.bin file written by `spectral`.
! Output: <output_prefix>.<tag>.map.dat for every map expression and
!         <output_prefix>.kticks.gp (gnuplot include with the k-path ticks).
! Exit  : 0 on success; 1 on wrong usage or any error (message prefixed
!         "[specmap error]" on stderr). All inputs and map expressions are
!         validated before the first file is written.
!
! User documentation: docs/specmap.md
! ----------------------------------------------------------------------------------------------------------
program specmap_main
    use, intrinsic :: iso_fortran_env, only: int64, error_unit
    use types,               only: dp
    use spectral_data,       only: spectral_data_t, read_spectral_file, &
                                   METHOD_MULLIKEN, METHOD_LOWDIN, &
                                   EREF_UNKNOWN, EREF_FERMI, EREF_ABSOLUTE
    use spectral_broadening, only: KERNEL_LORENTZIAN, KERNEL_GAUSSIAN, &
                                   energy_grid, broaden, fermi_dirac
    use spectral_map,        only: parse_map_expression, write_map_file, &
                                    write_kticks_file, map_file_tag
    use terminal_ui,         only: print_box, print_box_close
    implicit none

    integer, parameter :: MAX_MAPS = 32

    ! Namelist variables (defaults here are the documented defaults)
    character(len=512) :: spectral_file = ''
    integer            :: spin_channel  = 1
    character(len=32)  :: kernel        = 'lorentzian'
    real(dp)           :: width_ev      = 0.05_dp
    real(dp)           :: emin = -3.0_dp, emax = 3.0_dp
    integer            :: ne = 601
    character(len=256) :: maps(MAX_MAPS) = ''
    logical            :: occupation    = .false.
    real(dp)           :: temperature_k = 0.0_dp
    real(dp)           :: mu_ev         = 0.0_dp
    logical            :: mu_set        = .false.
    character(len=512) :: output_prefix = 'specmap'

    namelist /SPECMAP/ spectral_file, spin_channel, kernel, width_ev, emin, emax, ne, &
                       maps, occupation, temperature_k, mu_ev, mu_set, output_prefix

    type(spectral_data_t) :: sd
    character(len=1024)   :: input_path, errmsg, iomsg
    character(len=:), allocatable :: tag, fname, kname, wmeaning, occ_line, eref_line
    character(len=640)    :: hdr(12)     ! header lines of the .map.dat file (11 used)
    integer  :: u, ios, ierr, im, jm, it, ie, nmaps, kern, map_number
    integer(int64) :: total_start, map_start, map_end, clock_rate
    real(dp) :: mu, map_seconds, total_seconds
    real(dp), allocatable :: grid(:), wcomb(:,:), A(:,:), focc(:)
    integer,  allocatable :: ridx(:)
    real(dp), allocatable :: coef(:)

    ! ================================================================
    ! 1. Read and validate inputs
    ! ================================================================
    if (command_argument_count() /= 1) then
        write(error_unit, '(A)') "Usage: specmap <specmap.nml>"
        stop 1
    end if
    call get_command_argument(1, input_path)
    call system_clock(total_start, clock_rate)

    call print_box('specmap | broadened spectral maps')
    print '(A,A)', '  Input file       : ', trim(input_path)
    call print_box_close()
    print '(A)', ''
    print '(A)', '[1/3] Reading and validating inputs'

    open(newunit=u, file=trim(input_path), status='old', action='read', iostat=ios)
    if (ios /= 0) call die("cannot open input file: " // trim(input_path))
    read(u, nml=SPECMAP, iostat=ios, iomsg=iomsg)
    close(u)
    if (ios /= 0) call die("error reading &SPECMAP namelist from " // trim(input_path) // &
                           ": " // trim(iomsg))

    ! Option validation (cheap checks first, before reading the data file)
    if (len_trim(spectral_file) == 0) call die("spectral_file is required")
    select case (lower(trim(adjustl(kernel))))
    case ('lorentzian')
        kern = KERNEL_LORENTZIAN
        kname = 'lorentzian'
        wmeaning = 'HWHM eta'
    case ('gaussian')
        kern = KERNEL_GAUSSIAN
        kname = 'gaussian'
        wmeaning = 'standard deviation sigma'
    case default
        call die("kernel must be 'lorentzian' or 'gaussian', got '" // trim(kernel) // "'")
    end select
    if (.not. (width_ev > 0.0_dp)) call die("width_ev must be > 0")
    if (ne < 2) call die("ne must be >= 2")
    if (.not. (emax > emin)) call die("emax must be > emin")
    if (temperature_k < 0.0_dp) call die("temperature_k must be >= 0")
    if (len_trim(output_prefix) == 0) call die("output_prefix must not be empty")
    nmaps = count(len_trim(maps) > 0)
    if (nmaps == 0) call die("no map expression given in 'maps'")

    call read_spectral_file(trim(spectral_file), sd, ierr, errmsg)
    if (ierr /= 0) call die(trim(errmsg))

    if (spin_channel < 1 .or. spin_channel > sd%nspin) &
        call die("spin_channel must be in 1.." // itoa(sd%nspin))

    ! Chemical potential for the Fermi-Dirac occupation. Only a 'fermi'
    ! reference has a known zero (E_F = 0); otherwise the user must give mu.
    mu = 0.0_dp
    if (occupation) then
        if (mu_set) then
            mu = mu_ev
        else if (sd%energy_ref == EREF_FERMI) then
            mu = 0.0_dp
        else if (sd%energy_ref == EREF_UNKNOWN) then
            call die("occupation with energy reference 'unknown' requires mu_set = .true. and mu_ev")
        else
            call die("occupation with an absolute energy reference requires mu_set = .true. and mu_ev")
        end if
    end if

    ! Validate every map expression before writing anything, and reject
    ! expressions that map to the same output file name (would overwrite).
    do im = 1, MAX_MAPS
        if (len_trim(maps(im)) == 0) cycle
        call parse_map_expression(trim(maps(im)), sd, ridx, coef, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))
        do jm = 1, im - 1
            if (len_trim(maps(jm)) == 0) cycle
            if (map_file_tag(trim(maps(jm))) == map_file_tag(trim(maps(im)))) &
                call die("maps '" // trim(maps(jm)) // "' and '" // trim(maps(im)) // &
                         "' would write the same output file")
        end do
    end do

    ! ================================================================
    ! 2. Set up energy grid, occupation and header text
    ! ================================================================
    allocate(grid(ne), A(ne, sd%nk), wcomb(sd%nband, sd%nk), focc(ne))
    call energy_grid(emin, emax, ne, grid)
    if (occupation) then
        focc = fermi_dirac(grid, mu, temperature_k)
        occ_line = 'occupation: Fermi-Dirac, T = ' // rtoa(temperature_k) // ' K, mu = ' // &
                   rtoa(mu) // ' eV (intensity = f(E) * A)'
    else
        focc = 1.0_dp
        occ_line = 'occupation: none'
    end if

    select case (sd%energy_ref)
    case (EREF_FERMI)
        eref_line = 'energy reference: fermi (E - E_F)'
    case (EREF_ABSOLUTE)
        eref_line = 'energy reference: absolute'
    case default
        eref_line = 'energy reference: unknown'
    end select
    if (sd%fermi_known) then
        eref_line = eref_line // ', E_F = ' // rtoa(sd%fermi_energy_ev) // ' eV'
    else
        eref_line = eref_line // ', E_F not known'
    end if

    print '(A,I0,A,I0,A,I0,A)', '      source: ', sd%nk, ' k-points, ', sd%nband, &
        ' bands, ', sd%nregion, ' regions'
    print '(A,A,A,A,A,I0)', '      grid: ', trim(kname), ', width ', rtoa(width_ev), &
        ' eV, energy points ', ne
    print '(A)', '[2/3] Broadening and writing maps'

    ! ================================================================
    ! 3. One map per expression: combine, broaden, weight, write
    !    (expressions are parsed again here; they were validated above)
    ! ================================================================
    map_number = 0
    do im = 1, MAX_MAPS
        if (len_trim(maps(im)) == 0) cycle
        map_number = map_number + 1
        call system_clock(map_start)
        print '(A,I0,A,I0,A,A,A)', '      [', map_number, '/', nmaps, '] ', trim(maps(im)), ' ...'
        call parse_map_expression(trim(maps(im)), sd, ridx, coef, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))

        ! Combine region weights first (linear: identical to combining the
        ! broadened maps), then broaden once.
        wcomb = 0.0_dp
        do it = 1, size(ridx)
            wcomb = wcomb + coef(it) * sd%weight(:, :, ridx(it), spin_channel)
        end do
        call broaden(sd%energy(:, :, spin_channel), wcomb, grid, kern, width_ev, A)
        if (occupation) then
            do ie = 1, ne
                A(ie, :) = A(ie, :) * focc(ie)
            end do
        end if

        ! Metadata header written at the top of the .map.dat file
        hdr = ''
        hdr(1)  = 'specmap spectral map (broadened projection weights; NOT a photoemission simulation)'
        hdr(2)  = 'source file: ' // trim(spectral_file)
        hdr(3)  = 'source (from file): ' // trim(sd%source)
        hdr(4)  = 'map expression: ' // trim(maps(im)) // '   spin channel: ' // itoa(spin_channel)
        hdr(5)  = 'method: ' // method_label(sd%method)
        hdr(6)  = 'kernel: ' // kname // ', width = ' // rtoa(width_ev) // ' eV (' // wmeaning // ')'
        hdr(7)  = 'energy window: [' // rtoa(emin) // ', ' // rtoa(emax) // '] eV, ne = ' // itoa(ne)
        hdr(8)  = occ_line
        hdr(9)  = 'units: k_distance 1/Angstrom, energy eV, intensity states/eV ' // &
                  '(Mulliken/Lowdin-weighted)'
        hdr(10) = eref_line
        hdr(11) = 'regions combined with coefficients +/-1, no normalisation; ' // &
                  'overlapping regions in a sum may double-count atoms'

        tag = map_file_tag(trim(maps(im)))
        fname = trim(output_prefix) // '.' // tag // '.map.dat'
        call write_map_file(fname, sd, A, grid, hdr(1:11), ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))
        call system_clock(map_end)
        map_seconds = real(map_end - map_start, dp) / real(max(clock_rate, 1_int64), dp)
        print '(A,A,A,F0.3,A)', '          wrote ', fname, ' in ', map_seconds, ' s'
    end do

    ! ================================================================
    ! 4. k-path ticks for gnuplot, then summary
    ! ================================================================
    print '(A)', '[3/3] Writing k-path metadata'
    fname = trim(output_prefix) // '.kticks.gp'
    call write_kticks_file(fname, sd, ierr, errmsg)
    if (ierr /= 0) call die(trim(errmsg))
    print '(A,A)', '      wrote ', fname

    call system_clock(map_end)
    total_seconds = real(map_end - total_start, dp) / real(max(clock_rate, 1_int64), dp)
    print '(A)', ''
    call print_box('specmap | summary')
    print '(A,I0)',   '  Maps written     : ', nmaps
    print '(A,A)',    '  Output prefix    : ', trim(output_prefix)
    print '(A,F10.3,A)', '  Total time       : ', total_seconds, ' s'
    call print_box_close()

contains

    !> Print the message to stderr and exit with code 1.
    subroutine die(msg)
        character(len=*), intent(in) :: msg
        write(error_unit, '(A, A)') "[specmap error] ", msg
        stop 1
    end subroutine die

    !> ASCII lower-casing (kernel name is case-insensitive).
    pure function lower(s) result(r)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: r
        integer :: i, c
        r = s
        do i = 1, len(s)
            c = iachar(s(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) r(i:i) = achar(c + 32)
        end do
    end function lower

    function itoa(i) result(r)
        integer, intent(in) :: i
        character(len=:), allocatable :: r
        character(len=32) :: buf
        write(buf, '(I0)') i
        r = trim(buf)
    end function itoa

    !> Real to compact string for headers/log: fixed notation with trailing
    !> zeros trimmed (one decimal kept, e.g. '3.0') for 1e-3 <= |x| < 1e6,
    !> scientific notation otherwise.
    function rtoa(x) result(r)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: r
        character(len=32) :: buf
        integer :: n
        if (x == 0.0_dp .or. (abs(x) >= 1.0e-3_dp .and. abs(x) < 1.0e6_dp)) then
            write(buf, '(F0.6)') x
            r = trim(adjustl(buf))
            n = len(r)
            do while (n > 1)
                if (r(n:n) /= '0') exit
                n = n - 1
            end do
            if (r(n:n) == '.') n = n + 1    ! keep one decimal digit: '3.0'
            r = r(1:n)
            if (r(1:1) == '.') r = '0' // r
            if (len(r) > 1) then
                if (r(1:2) == '-.') r = '-0' // r(2:)
            end if
        else
            write(buf, '(ES14.6E3)') x
            r = trim(adjustl(buf))
        end if
    end function rtoa

    !> Name of the projection method stored in the spectral file.
    function method_label(m) result(r)
        integer, intent(in) :: m
        character(len=:), allocatable :: r
        select case (m)
        case (METHOD_MULLIKEN); r = 'mulliken'
        case (METHOD_LOWDIN);   r = 'lowdin'
        case default;           r = 'unknown'
        end select
    end function method_label

end program specmap_main
