! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> skfit — parâmetros de Slater–Koster para o contato do surfgreen.
!>
!>   skfit FIT.nml      (&SKFIT)      ajusta SK ao acoplamento entre camadas
!>                                    principais de um bulk e testa as bandas
!>   skfit COMBINE.nml  (&SKCOMBINE)  combina dois ajustes (X-X de A, Y-Y de B)
!>                                    num arquivo de contato X-Y
!>
!> Ver docs/sk_contact.md.
program skfit
    use, intrinsic :: iso_fortran_env, only: error_unit, iostat_end
    use types, only: dp
    use HS_io, only: HSDataFormat
    use structure_io, only: structure_t, read_structure
    use layer_blocks, only: pl_blocks_t, build_pl_blocks, attach_pl_geometry
    use sk_contact, only: sk_model_t, sk_entry_t, read_sk_file, write_sk_file, sk_set_bases
    use sk_fit, only: fit_shell_set_t, sk_fit_report_t, parse_fit_shells, sk_fit_material, sk_band_check, &
                      sk_combine
    implicit none

    integer, parameter :: PLEN = 512, NB = 32
    character(len=1024) :: arg
    integer :: u, ios
    character(len=512) :: iomsg
    logical :: ex

    if (command_argument_count() /= 1) call usage()
    call get_command_argument(1, arg)
    if (trim(arg) == '-h' .or. trim(arg) == '--help') call usage()
    inquire(file=trim(arg), exist=ex)
    if (.not. ex) call die("input file not found: '" // trim(arg) // "'")
    open(newunit=u, file=trim(arg), status='old', action='read')
    if (has_namelist('SKCOMBINE')) then
        call do_combine()
    else
        call do_fit()
    end if
    close(u)

contains

    subroutine usage()
        write(error_unit, '(A)') 'Usage: skfit INPUT.nml   (namelist &SKFIT or &SKCOMBINE; see docs/sk_contact.md)'
        stop 1
    end subroutine usage

    subroutine die(msg)
        character(len=*), intent(in) :: msg
        write(error_unit, '(A,A)') '[skfit error] ', trim(msg)
        stop 1
    end subroutine die

    logical function has_namelist(name)
        character(len=*), intent(in) :: name
        character(len=1024) :: line
        has_namelist = .false.
        rewind(u)
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (line(1:1) == '&') then
                if (index(line, '&' // name) == 1) has_namelist = .true.
            end if
        end do
        rewind(u)
    end function has_namelist

    subroutine do_fit()
        character(len=PLEN) :: h_file, s_file, structure_file, output
        integer :: stacking_axis, pl_cells, check_nk
        character(len=64) :: basis(NB), fit_shells(NB)
        real(dp) :: cutoff_ang, anchor_ev, beta, beta_min, beta_max, check_kpar(2), check_window_ev
        namelist /SKFIT/ h_file, s_file, structure_file, stacking_axis, pl_cells, basis, fit_shells, &
                         cutoff_ang, anchor_ev, beta, beta_min, beta_max, output, check_kpar, check_nk, &
                         check_window_ev
        type(HSDataFormat) :: hh, hs
        type(structure_t) :: st
        type(pl_blocks_t) :: pl
        type(sk_model_t) :: model
        type(fit_shell_set_t), allocatable :: sets(:)
        type(sk_fit_report_t) :: rep
        character(len=:), allocatable :: errmsg
        character(len=160) :: hdr(8)
        character(len=512) :: msg
        integer :: ierr, nval, k
        real(dp) :: anchor, ma, ra, ms, rs

        h_file = ''; s_file = ''; structure_file = ''; output = 'fit.sk'
        stacking_axis = 3; pl_cells = 0
        basis = ''; fit_shells = ''
        cutoff_ang = 6.0_dp; anchor_ev = huge(1.0_dp)
        beta = -1.0_dp; beta_min = 0.0_dp; beta_max = 4.0_dp
        check_kpar = 0.0_dp; check_nk = 41; check_window_ev = 1.0_dp
        read(u, nml=SKFIT, iostat=ios, iomsg=iomsg)
        if (ios == iostat_end) call die('namelist &SKFIT not found')
        if (ios /= 0) call die('error reading &SKFIT: ' // trim(iomsg))
        if (len_trim(h_file) == 0 .or. len_trim(structure_file) == 0) call die('&SKFIT h_file and structure_file are required')
        if (count(len_trim(basis) > 0) == 0) call die("&SKFIT basis is required (e.g. basis = 'Te s3p2d2', 'Bi s3p2d2')")

        inquire(file=trim(h_file), exist=ex)
        if (.not. ex) call die("H file not found: '" // trim(h_file) // "'")
        call hh%ReadFile(trim(h_file))
        if (len_trim(s_file) > 0) then
            inquire(file=trim(s_file), exist=ex)
            if (.not. ex) call die("S file not found: '" // trim(s_file) // "'")
            call hs%ReadFile(trim(s_file))
            call build_pl_blocks(hh, stacking_axis, pl_cells, pl, ierr, errmsg, hs=hs)
        else
            call build_pl_blocks(hh, stacking_axis, pl_cells, pl, ierr, errmsg)
        end if
        if (ierr /= 0) call die(errmsg)
        call read_structure(trim(structure_file), st, ierr, msg)
        if (ierr /= 0) call die('structure: ' // trim(msg))
        call attach_pl_geometry(pl, st, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)

        ! âncora: E_F absoluto se o .hsio não foi deslocado; 0 se foi (E_F = 0)
        if (anchor_ev /= huge(1.0_dp)) then
            anchor = anchor_ev
        else if (hh%has_meta) then
            anchor = merge(0.0_dp, hh%fermi_ev, hh%fermi_shift_applied)
        else
            anchor = 0.0_dp
            print '(A)', '  WARNING: no energy metadata in the .hsio; anchor_ev = 0 (set anchor_ev explicitly)'
        end if

        model%cutoff = cutoff_ang
        call sk_set_bases(model, basis, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        call parse_fit_shells(fit_shells, sets, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        call sk_fit_material(pl, st%label, model, sets, anchor, beta, beta_min, beta_max, rep, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        call sk_band_check(pl, st%label, model, check_kpar, check_nk, anchor, check_window_ev, ma, ra, ms, rs, &
                           nval, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)

        hdr = ''
        hdr(1) = 'skfit: Slater-Koster fit of the coupling between principal layers'
        hdr(2) = 'H file: ' // trim(h_file)
        write(hdr(3), '(A,F10.5,A,F8.3,A,I0,A)') 'anchor E0 = ', anchor, ' eV, cutoff = ', cutoff_ang, &
            ' A, PL = ', pl%ncell, ' cell(s)'
        write(hdr(4), '(A,F8.4,A,F8.4,A,F8.4)') 'beta = ', rep%beta, ' 1/A; rel. residual (fitted shells) = ', &
            rep%rel_res_fitted, '; rel. error full H01 = ', rep%rel_err_full
        write(hdr(5), '(A,F10.2,A,F10.2,A,F5.2,A)') 'bulk bands, max |dE| (meV): anchored ', 1000.0_dp * ma, &
            ', SK ', 1000.0_dp * ms, ' (window +-', check_window_ev, ' eV)'
        hdr(6) = 'S_AB = 0 in the contact; parameters refer to the SK table with the lower l first'
        call write_sk_file(trim(output), model%e(1:model%nentry), hdr, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)

        print '(A)', 'skfit | Slater-Koster fit of the inter-layer coupling'
        print '(A,A)',        '  H file                  : ', trim(h_file)
        print '(A,I0,A,I0,A)', '  Principal layer         : ', pl%ncell, ' cell(s), ', pl%norb, ' orbitals'
        print '(A,F10.5,A)',  '  Anchor E0               : ', anchor, ' eV  (target H01 - E0 S01, S01 -> 0)'
        print '(A,F8.3,A)',   '  Cutoff                  : ', cutoff_ang, ' A'
        print '(A,I0,A,I0,A,I0)', '  Pairs / data / params   : ', rep%npairs, ' / ', rep%nrows, ' / ', rep%nparams
        print '(A,F8.4,A)',   '  beta                    : ', rep%beta, ' 1/A'
        print '(A,F8.4)',     '  Rel. residual (fitted)  : ', rep%rel_res_fitted
        print '(A,F8.4)',     '  Rel. error, full H01    : ', rep%rel_err_full
        print '(A,F8.4)',     '  Dropped (Im, spin-off)  : ', rep%dropped_spin
        print '(A,I0,A,F5.2,A)', '  Bulk bands check        : ', nval, ' eigenvalues within +-', check_window_ev, &
            ' eV of E0 along k_s'
        print '(A,2F10.2)',   '    anchored  max/rms meV : ', 1000.0_dp * ma, 1000.0_dp * ra
        print '(A,2F10.2)',   '    SK        max/rms meV : ', 1000.0_dp * ms, 1000.0_dp * rs
        print '(A,I0,A,A)',   '  Wrote ', model%nentry, ' lines to ', trim(output)
        do k = 1, min(model%nentry, 12)
            print '(4X,A,1X,A,I0,2X,A,1X,A,I0,F9.3,3F11.5)', trim(model%e(k)%spec_a), &
                'spdf'(model%e(k)%la + 1:model%e(k)%la + 1), model%e(k)%za, trim(model%e(k)%spec_b), &
                'spdf'(model%e(k)%lb + 1:model%e(k)%lb + 1), model%e(k)%zb, model%e(k)%d0, model%e(k)%v
        end do
        if (model%nentry > 12) print '(4X,A)', '...'
    end subroutine do_fit

    subroutine do_combine()
        character(len=PLEN) :: file_a, file_b, output
        character(len=16) :: species_a(NB), species_b(NB)
        namelist /SKCOMBINE/ file_a, file_b, species_a, species_b, output
        type(sk_model_t) :: ma, mb
        type(sk_entry_t), allocatable :: out(:)
        character(len=:), allocatable :: errmsg
        character(len=160) :: hdr(4)
        integer :: ierr, nwarn
        file_a = ''; file_b = ''; output = 'contact.sk'; species_a = ''; species_b = ''
        read(u, nml=SKCOMBINE, iostat=ios, iomsg=iomsg)
        if (ios /= 0) call die('error reading &SKCOMBINE: ' // trim(iomsg))
        call read_sk_file(trim(file_a), ma, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        call read_sk_file(trim(file_b), mb, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        call sk_combine(ma%e, mb%e, species_a, species_b, out, nwarn)
        if (size(out) == 0) call die('no homonuclear X-X lines in file_a matching Y-Y lines in file_b ' // &
                                     '(check species_a/species_b and the shells of the two fits)')
        hdr = ''
        hdr(1) = 'skfit --combine: geometric mean of homonuclear lines, at the mean d0'
        hdr(2) = 'from ' // trim(file_a) // ' and ' // trim(file_b)
        write(hdr(3), '(I0,A)') nwarn, ' channel(s) set to 0 because the two signs differ'
        hdr(4) = 'S_AB = 0; a starting point for the scan, not a prediction'
        call write_sk_file(trim(output), out, hdr, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        print '(A,I0,A,A)', 'skfit | combined ', size(out), ' lines -> ', trim(output)
        if (nwarn > 0) print '(A,I0,A)', '  WARNING: ', nwarn, ' channel(s) with opposite signs set to 0'
    end subroutine do_combine

end program skfit
