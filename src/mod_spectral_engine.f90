! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Spectral engine (contracts §3, §3.2, §4, §5): reads spectral.nml, the HSIO
!> H/S files, the structure file and regions.nml, diagonalizes H(k) along a
!> k-path and computes per-region projection weights (Mulliken).
!> The per-k work lives in `process_kpoint`, which only touches the arrays and
!> workspaces passed to it (ready for a later OpenMP loop).
!> Library routines never stop: errors are returned through ierr/errmsg.
!> (Exception: HSDataFormat%ReadFile stops on a corrupt HSIO file; missing
!> files are detected here before calling it.)
module spectral_engine
    use, intrinsic :: iso_fortran_env, only: int64, iostat_end
    use types, only: dp
    use hs_io, only: HSDataFormat
    use core, only: compute_Hk_at_kpoint
    use structure_io, only: structure_t, read_structure, structure_fingerprint, &
                            reciprocal_lattice, species_norb_of
    use atomic_regions, only: region_set_t, read_regions_file, validate_regions, &
                              regions_are_partition
    use spectral_projection, only: orbital_groups_t, build_orbital_groups, &
                                   overlap_times_vectors, mulliken_weights, lowdin_weights, &
                                   generalized_eigensolve
    use kpath, only: kpath_t, build_kpath
    use spectral_data, only: spectral_data_t, SPIN_NONE, SPIN_COLLINEAR, SPIN_NONCOLLINEAR, &
                             METHOD_MULLIKEN, METHOD_LOWDIN, EREF_UNKNOWN, EREF_FERMI, EREF_ABSOLUTE
    implicit none
    private

    integer, parameter, public :: MAX_KVERTEX = 100
    integer, parameter, public :: SPEC_PATH_LEN = 512
    integer, parameter :: MAX_WARN_LINES = 20

    type, public :: spectral_input_t
        character(len=SPEC_PATH_LEN) :: h_file = ''
        character(len=SPEC_PATH_LEN) :: h_up_file = ''
        character(len=SPEC_PATH_LEN) :: h_dn_file = ''
        character(len=SPEC_PATH_LEN) :: s_file = ''
        character(len=SPEC_PATH_LEN) :: structure_file = ''
        character(len=SPEC_PATH_LEN) :: regions_file = ''
        logical :: allow_no_fingerprint = .false.
        character(len=16) :: spin_mode = 'none'
        character(len=16) :: method = 'mulliken'
        character(len=16) :: energy_reference = 'unknown'
        real(dp) :: fermi_energy_ev = 0.0_dp
        logical :: fermi_known = .false.
        character(len=SPEC_PATH_LEN) :: output_prefix = 'slab'
        logical :: write_text = .false.
        real(dp) :: tol_norm = 1.0e-8_dp
        real(dp) :: tol_negative_warn = 1.0e-3_dp
        real(dp) :: tol_s_min = 1.0e-10_dp
        ! &KPATH
        integer :: nvertex = 0
        real(dp) :: vertex(3, MAX_KVERTEX) = 0.0_dp
        character(len=16) :: label(MAX_KVERTEX) = ''
        integer :: npoints(MAX_KVERTEX - 1) = 0
    end type spectral_input_t

    type, public :: spectral_timing_t
        real(dp) :: read_s  = 0.0_dp   ! structure, regions, HSIO
        real(dp) :: build_s = 0.0_dp   ! H(k), S(k) construction
        real(dp) :: diag_s  = 0.0_dp   ! generalized eigensolver
        real(dp) :: proj_s  = 0.0_dp   ! S*C + weights
        real(dp) :: write_s = 0.0_dp   ! filled by the driver
        real(dp) :: total_s = 0.0_dp
    end type spectral_timing_t

    public :: read_spectral_input, run_spectral, process_kpoint

contains

    ! ==================================================================
    ! Input
    ! ==================================================================

    subroutine read_spectral_input(filename, inp, ierr, errmsg)
        character(len=*),       intent(in)  :: filename
        type(spectral_input_t), intent(out) :: inp
        integer,                intent(out) :: ierr
        character(len=*),       intent(out) :: errmsg

        character(len=SPEC_PATH_LEN) :: h_file, h_up_file, h_dn_file, s_file, structure_file, &
                                        regions_file, output_prefix
        logical :: allow_no_fingerprint, fermi_known, write_text
        character(len=16) :: spin_mode, method, energy_reference
        real(dp) :: fermi_energy_ev, tol_norm, tol_negative_warn, tol_s_min
        integer :: nvertex
        real(dp) :: vertex(3, MAX_KVERTEX)
        character(len=16) :: label(MAX_KVERTEX)
        integer :: npoints(MAX_KVERTEX - 1)

        namelist /SPECTRAL/ h_file, h_up_file, h_dn_file, s_file, structure_file, regions_file, &
                            allow_no_fingerprint, spin_mode, method, energy_reference, &
                            fermi_energy_ev, fermi_known, output_prefix, write_text, &
                            tol_norm, tol_negative_warn, tol_s_min
        namelist /KPATH/ nvertex, vertex, label, npoints

        integer :: iu, ios, is
        character(len=512) :: iomsg
        logical :: exists

        ierr = 0
        errmsg = ''

        h_file = inp%h_file;               h_up_file = inp%h_up_file
        h_dn_file = inp%h_dn_file;         s_file = inp%s_file
        structure_file = inp%structure_file; regions_file = inp%regions_file
        allow_no_fingerprint = inp%allow_no_fingerprint
        spin_mode = inp%spin_mode;         method = inp%method
        energy_reference = inp%energy_reference
        fermi_energy_ev = inp%fermi_energy_ev; fermi_known = inp%fermi_known
        output_prefix = inp%output_prefix; write_text = inp%write_text
        tol_norm = inp%tol_norm;           tol_negative_warn = inp%tol_negative_warn
        tol_s_min = inp%tol_s_min
        nvertex = inp%nvertex;             vertex = inp%vertex
        label = inp%label;                 npoints = inp%npoints

        inquire(file=filename, exist=exists)
        if (.not. exists) then
            call fail("input file not found: '" // trim(filename) // "'")
            return
        end if
        open(newunit=iu, file=filename, status='old', action='read', iostat=ios, iomsg=iomsg)
        if (ios /= 0) then
            call fail("cannot open '" // trim(filename) // "': " // trim(iomsg))
            return
        end if

        read(iu, nml=SPECTRAL, iostat=ios, iomsg=iomsg)
        if (ios == iostat_end) then
            close(iu)
            call fail("namelist &SPECTRAL not found in '" // trim(filename) // "'")
            return
        else if (ios /= 0) then
            close(iu)
            call fail("error reading &SPECTRAL: " // trim(iomsg))
            return
        end if
        rewind(iu)
        read(iu, nml=KPATH, iostat=ios, iomsg=iomsg)
        if (ios == iostat_end) then
            close(iu)
            call fail("namelist &KPATH not found in '" // trim(filename) // "'")
            return
        else if (ios /= 0) then
            close(iu)
            call fail("error reading &KPATH: " // trim(iomsg))
            return
        end if
        close(iu)

        inp%h_file = h_file;               inp%h_up_file = h_up_file
        inp%h_dn_file = h_dn_file;         inp%s_file = s_file
        inp%structure_file = structure_file; inp%regions_file = regions_file
        inp%allow_no_fingerprint = allow_no_fingerprint
        inp%spin_mode = lower(adjustl(spin_mode))
        inp%method = lower(adjustl(method))
        inp%energy_reference = lower(adjustl(energy_reference))
        inp%fermi_energy_ev = fermi_energy_ev; inp%fermi_known = fermi_known
        inp%output_prefix = output_prefix; inp%write_text = write_text
        inp%tol_norm = tol_norm;           inp%tol_negative_warn = tol_negative_warn
        inp%tol_s_min = tol_s_min
        inp%nvertex = nvertex;             inp%vertex = vertex
        inp%label = label;                 inp%npoints = npoints

        ! --- enums -------------------------------------------------------
        select case (trim(inp%spin_mode))
        case ('none', 'noncollinear')
            if (len_trim(inp%h_file) == 0) then
                call fail("h_file is required for spin_mode='" // trim(inp%spin_mode) // "'")
                return
            end if
        case ('collinear')
            if (len_trim(inp%h_up_file) == 0 .or. len_trim(inp%h_dn_file) == 0) then
                call fail("h_up_file and h_dn_file are required for spin_mode='collinear'")
                return
            end if
        case default
            call fail("invalid spin_mode '" // trim(inp%spin_mode) // &
                      "' (expected 'none', 'collinear' or 'noncollinear')")
            return
        end select

        select case (trim(inp%method))
        case ('mulliken', 'lowdin')
        case default
            call fail("invalid method '" // trim(inp%method) // "' (expected 'mulliken' or 'lowdin')")
            return
        end select

        select case (trim(inp%energy_reference))
        case ('fermi')
            inp%fermi_known = .true.
            inp%fermi_energy_ev = 0.0_dp
        case ('absolute', 'unknown')
        case default
            call fail("invalid energy_reference '" // trim(inp%energy_reference) // &
                      "' (expected 'fermi', 'absolute' or 'unknown')")
            return
        end select

        ! --- required files ----------------------------------------------
        if (len_trim(inp%s_file) == 0) then
            call fail('s_file is required')
            return
        end if
        if (len_trim(inp%structure_file) == 0) then
            call fail('structure_file is required')
            return
        end if
        if (len_trim(inp%output_prefix) == 0) then
            call fail('output_prefix must not be empty')
            return
        end if

        ! --- tolerances --------------------------------------------------
        if (.not. (inp%tol_norm > 0.0_dp)) then
            call fail('tol_norm must be > 0')
            return
        end if
        if (.not. (inp%tol_negative_warn > 0.0_dp)) then
            call fail('tol_negative_warn must be > 0')
            return
        end if
        if (.not. (inp%tol_s_min > 0.0_dp)) then
            call fail('tol_s_min must be > 0')
            return
        end if

        ! --- k-path ------------------------------------------------------
        if (inp%nvertex < 2 .or. inp%nvertex > MAX_KVERTEX) then
            call fail('&KPATH nvertex must be in 2..' // itoa(MAX_KVERTEX) // &
                      ' (got ' // itoa(inp%nvertex) // ')')
            return
        end if
        do is = 1, inp%nvertex - 1
            if (inp%npoints(is) < 1) then
                call fail('&KPATH npoints(' // itoa(is) // ') must be >= 1 (got ' // &
                          itoa(inp%npoints(is)) // ')')
                return
            end if
        end do

    contains

        subroutine fail(m)
            character(len=*), intent(in) :: m
            ierr = 1
            errmsg = m
        end subroutine fail

    end subroutine read_spectral_input

    ! ==================================================================
    ! Driver of the calculation
    ! ==================================================================

    subroutine run_spectral(inp, sd, timing, ierr, errmsg)
        type(spectral_input_t),  intent(in)  :: inp
        type(spectral_data_t),   intent(out) :: sd
        type(spectral_timing_t), intent(out) :: timing
        integer,                 intent(out) :: ierr
        character(len=*),        intent(out) :: errmsg

        type(structure_t)  :: st
        type(region_set_t) :: rs
        type(HSDataFormat), allocatable :: hdat(:)
        type(HSDataFormat) :: sdat
        type(orbital_groups_t) :: groups
        type(kpath_t) :: path
        integer, allocatable :: list_start(:), list_atoms(:)
        complex(dp), allocatable :: Hk(:,:), Sk(:,:), C(:,:), SC(:,:)
        real(dp), allocatable :: E(:), W(:,:)
        character(len=SPEC_PATH_LEN) :: hname(2)
        integer(int64) :: c_start, c0, c1, crate
        integer :: nspin, nreg, ngroup, nat, n, ispin, ik, ia, ir, g, pos, smode, kerr, perr
        logical :: have_regions, partition, use_lowdin
        real(dp) :: smin
        character(len=512) :: msg

        ierr = 0
        errmsg = ''
        call system_clock(c_start, crate)
        use_lowdin = trim(inp%method) == 'lowdin'

        if (trim(inp%method) /= 'mulliken' .and. trim(inp%method) /= 'lowdin') then
            call fail("invalid method '" // trim(inp%method) // "'")
            return
        end if

        select case (trim(inp%spin_mode))
        case ('none')
            smode = SPIN_NONE;          nspin = 1; hname(1) = inp%h_file
        case ('noncollinear')
            smode = SPIN_NONCOLLINEAR;  nspin = 1; hname(1) = inp%h_file
        case ('collinear')
            smode = SPIN_COLLINEAR;     nspin = 2
            hname(1) = inp%h_up_file;   hname(2) = inp%h_dn_file
        case default
            call fail("invalid spin_mode '" // trim(inp%spin_mode) // "'")
            return
        end select

        ! ---------------- 1. read inputs --------------------------------
        call system_clock(c0)
        call read_structure(inp%structure_file, st, ierr, msg)
        if (ierr /= 0) then
            call fail('structure: ' // trim(msg))
            return
        end if

        have_regions = len_trim(inp%regions_file) > 0
        if (have_regions) then
            call read_regions_file(inp%regions_file, rs, ierr, msg)
            if (ierr /= 0) then
                call fail('regions: ' // trim(msg))
                return
            end if
            call validate_regions(rs, st, inp%allow_no_fingerprint, ierr, msg)
            if (ierr /= 0) then
                call fail('regions: ' // trim(msg))
                return
            end if
        end if

        if (.not. file_exists(inp%s_file)) then
            call fail("S file not found: '" // trim(inp%s_file) // "'")
            return
        end if
        do ispin = 1, nspin
            if (.not. file_exists(hname(ispin))) then
                call fail("H file not found: '" // trim(hname(ispin)) // "'")
                return
            end if
        end do
        allocate(hdat(nspin))
        call sdat%ReadFile(trim(inp%s_file))
        do ispin = 1, nspin
            call hdat(ispin)%ReadFile(trim(hname(ispin)))
        end do
        call system_clock(c1)
        timing%read_s = elapsed(c0, c1, crate)

        ! ---------------- 2. consistency --------------------------------
        do ispin = 1, nspin
            call compare_hs(hdat(ispin), sdat, trim(hname(ispin)), trim(inp%s_file), ierr, msg)
            if (ierr /= 0) then
                call fail(trim(msg))
                return
            end if
        end do
        if (.not. allocated(sdat%Orbs) .or. .not. allocated(sdat%Data) .or. &
            .not. allocated(sdat%IndexList)) then
            call fail("S file '" // trim(inp%s_file) // "' has no data")
            return
        end if
        nat = sdat%nAtoms
        n = sdat%NOrbs
        if (size(sdat%Orbs) /= nat .or. sum(sdat%Orbs) /= n) then
            call fail("HSIO Orbs inconsistent with nAtoms/NOrbs in '" // trim(inp%s_file) // "'")
            return
        end if
        if (nat /= st%nat) then
            call fail('HSIO nAtoms (' // itoa(nat) // ') differs from structure nat (' // &
                      itoa(st%nat) // ") in '" // trim(inp%structure_file) // "'")
            return
        end if
        if (st%has_species) then
            do ia = 1, nat
                if (species_norb_of(st, ia) /= sdat%Orbs(ia)) then
                    call fail('atom ' // itoa(ia) // " (species '" // trim(st%label(ia)) // &
                              "'): structure nOrb = " // itoa(species_norb_of(st, ia)) // &
                              ' but HSIO Orbs = ' // itoa(sdat%Orbs(ia)))
                    return
                end if
            end do
        end if

        ! ---------------- 3. orbital groups -----------------------------
        nreg = 0
        if (have_regions) nreg = rs%nregion
        ngroup = 1 + nreg
        allocate(list_start(ngroup + 1))
        list_start(1) = 1
        list_start(2) = 1 + nat
        do ir = 1, nreg
            list_start(ir + 2) = list_start(ir + 1) + size(rs%regions(ir)%atoms)
        end do
        allocate(list_atoms(list_start(ngroup + 1) - 1))
        list_atoms(1:nat) = [(ia, ia = 1, nat)]
        pos = nat
        do ir = 1, nreg
            list_atoms(pos + 1:pos + size(rs%regions(ir)%atoms)) = rs%regions(ir)%atoms
            pos = pos + size(rs%regions(ir)%atoms)
        end do
        call build_orbital_groups(sdat%Orbs, ngroup, list_start, list_atoms, groups)
        partition = .false.
        if (have_regions) partition = regions_are_partition(rs)

        ! ---------------- 4. k-path -------------------------------------
        call build_kpath(inp%vertex(:, 1:inp%nvertex), inp%npoints(1:inp%nvertex - 1), &
                         inp%label(1:inp%nvertex), reciprocal_lattice(st), path, ierr, msg)
        if (ierr /= 0) then
            call fail('k-path: ' // trim(msg))
            return
        end if

        ! ---------------- 5. spectral data metadata ---------------------
        sd%nk = path%nk
        sd%nband = n
        sd%nspin = nspin
        sd%nregion = ngroup
        sd%nvertex = path%nvertex
        sd%spin_mode = smode
        sd%method = METHOD_MULLIKEN
        if (use_lowdin) sd%method = METHOD_LOWDIN
        select case (trim(inp%energy_reference))
        case ('fermi');    sd%energy_ref = EREF_FERMI
        case ('absolute'); sd%energy_ref = EREF_ABSOLUTE
        case default;      sd%energy_ref = EREF_UNKNOWN
        end select
        sd%fermi_known = inp%fermi_known
        sd%fermi_energy_ev = 0.0_dp
        if (inp%fermi_known) sd%fermi_energy_ev = inp%fermi_energy_ev
        if (sd%energy_ref == EREF_FERMI) then
            sd%fermi_known = .true.
            sd%fermi_energy_ev = 0.0_dp
        end if
        sd%fingerprint = structure_fingerprint(st)
        if (smode == SPIN_COLLINEAR) then
            sd%source = 'H_up=' // trim(inp%h_up_file) // ' H_dn=' // trim(inp%h_dn_file) // &
                        ' S=' // trim(inp%s_file)
        else
            sd%source = 'H=' // trim(inp%h_file) // ' S=' // trim(inp%s_file)
        end if
        sd%kfrac = path%kfrac
        sd%kdist = path%kdist
        sd%vertex_k = path%vertex_k
        allocate(sd%vertex_label(path%nvertex))
        sd%vertex_label = path%label
        allocate(sd%region_name(ngroup), sd%region_natoms(ngroup))
        sd%region_name(1) = 'total'
        sd%region_natoms(1) = nat
        do ir = 1, nreg
            sd%region_name(ir + 1) = rs%regions(ir)%name
            sd%region_natoms(ir + 1) = size(rs%regions(ir)%atoms)
        end do
        allocate(sd%energy(n, path%nk, nspin), source=0.0_dp)
        allocate(sd%weight(n, path%nk, ngroup, nspin), source=0.0_dp)

        ! ---------------- 6. k loop -------------------------------------
        allocate(Hk(n, n), Sk(n, n), C(n, n), SC(n, n), E(n), W(n, ngroup))
        do ispin = 1, nspin
            do ik = 1, path%nk
                call process_kpoint(n, ngroup, hdat(ispin)%Data, sdat%Data, sdat%IndexList, &
                                    path%kfrac(:, ik), groups, Hk, Sk, C, SC, E, W, &
                                    timing%build_s, timing%diag_s, timing%proj_s, kerr, &
                                    lowdin=use_lowdin, tol_s_min=inp%tol_s_min, smin=smin, &
                                    proj_ierr=perr)
                if (perr == 1) then
                    write(msg, '(A,I0,A,I0,A,3F12.7,A,ES12.4,A,ES10.3,A)') &
                        'overlap matrix not positive definite above tol_s_min at ispin=', ispin, &
                        ' ik=', ik, ' kfrac=(', path%kfrac(:, ik), '): min eigenvalue smin=', smin, &
                        ' (tol_s_min=', inp%tol_s_min, ')'
                    call fail(trim(msg))
                    return
                else if (perr /= 0) then
                    write(msg, '(A,I0,A,I0,A,I0,A,3F12.7,A)') 'Lowdin projection failed (ierr=', &
                        perr, ') at ispin=', ispin, ' ik=', ik, ' kfrac=(', path%kfrac(:, ik), ')'
                    call fail(trim(msg))
                    return
                end if
                if (kerr /= 0) then
                    write(msg, '(A,I0,A,I0,A,I0,A,3F12.7,A)') 'generalized eigensolver failed (info=', &
                        kerr, ') at ispin=', ispin, ' ik=', ik, ' kfrac=(', path%kfrac(:, ik), ')'
                    call fail(trim(msg))
                    return
                end if
                sd%energy(:, ik, ispin) = E
                do g = 1, ngroup
                    sd%weight(:, ik, g, ispin) = W(:, g)
                end do
            end do
        end do

        ! ---------------- 7. checks -------------------------------------
        call check_weights(sd, partition, inp%tol_norm, inp%tol_negative_warn, use_lowdin)

        call system_clock(c1)
        timing%total_s = elapsed(c_start, c1, crate)

    contains

        subroutine fail(m)
            character(len=*), intent(in) :: m
            ierr = 1
            errmsg = m
        end subroutine fail

    end subroutine run_spectral

    ! ==================================================================
    ! Per-k body: touches only its arguments (thread-safe by construction).
    ! ==================================================================

    !> Optional Löwdin mode (lowdin=.true.): weights from lowdin_weights on the
    !> ORIGINAL Sk (SC is not used); tol_s_min, smin and proj_ierr are then
    !> required. proj_ierr = 1 means min eig(Sk) = smin <= tol_s_min.
    !> Without the optional arguments the Mulliken path is used (unchanged).
    subroutine process_kpoint(n, ngroup, Hr, Sr, Rpt, kfrac, groups, Hk, Sk, C, SC, E, W, &
                              t_build, t_diag, t_proj, ierr, lowdin, tol_s_min, smin, proj_ierr)
        integer,                intent(in)    :: n, ngroup
        complex(dp), contiguous, intent(in)   :: Hr(:,:,:), Sr(:,:,:)
        integer,     contiguous, intent(in)   :: Rpt(:,:)
        real(dp),               intent(in)    :: kfrac(3)
        type(orbital_groups_t), intent(in)    :: groups
        complex(dp), contiguous, intent(inout) :: Hk(:,:), Sk(:,:), C(:,:), SC(:,:)
        real(dp),               intent(inout) :: E(n), W(n, ngroup)
        real(dp),               intent(inout) :: t_build, t_diag, t_proj
        integer,                intent(out)   :: ierr
        logical,  optional,     intent(in)    :: lowdin
        real(dp), optional,     intent(in)    :: tol_s_min
        real(dp), optional,     intent(out)   :: smin
        integer,  optional,     intent(out)   :: proj_ierr
        integer(int64) :: c0, c1, c2, c3, crate
        logical :: do_lowdin
        real(dp) :: smin_loc
        integer :: perr

        do_lowdin = .false.
        if (present(lowdin)) do_lowdin = lowdin
        if (do_lowdin .and. .not. (present(tol_s_min) .and. present(proj_ierr))) &
            error stop 'process_kpoint: lowdin requires tol_s_min and proj_ierr'
        if (present(proj_ierr)) proj_ierr = 0
        if (present(smin)) smin = 0.0_dp

        call system_clock(c0, crate)
        call compute_Hk_at_kpoint(Hr, Rpt, Hk, kfrac)
        call compute_Hk_at_kpoint(Sr, Rpt, Sk, kfrac)
        call system_clock(c1)
        ! generalized_eigensolve works on internal copies: Sk stays the original S(k).
        call generalized_eigensolve(n, Hk, Sk, E, C, ierr)
        call system_clock(c2)
        if (ierr /= 0) then
            t_build = t_build + elapsed(c0, c1, crate)
            t_diag  = t_diag  + elapsed(c1, c2, crate)
            return
        end if
        if (do_lowdin) then
            call lowdin_weights(n, n, Sk, C, groups, W, tol_s_min, smin_loc, perr)
            proj_ierr = perr
            if (present(smin)) smin = smin_loc
        else
            call overlap_times_vectors(n, n, Sk, C, SC)
            call mulliken_weights(n, n, C, SC, groups, W)
        end if
        call system_clock(c3)
        t_build = t_build + elapsed(c0, c1, crate)
        t_diag  = t_diag  + elapsed(c1, c2, crate)
        t_proj  = t_proj  + elapsed(c2, c3, crate)
    end subroutine process_kpoint

    ! ==================================================================
    ! Helpers
    ! ==================================================================

    !> For Löwdin (lowdin=.true.) weights are non-negative by construction; the
    !> negative check then uses -tol_norm instead of -tol_negative_warn.
    subroutine check_weights(sd, partition, tol_norm, tol_neg_in, lowdin)
        type(spectral_data_t), intent(in) :: sd
        logical,               intent(in) :: partition
        real(dp),              intent(in) :: tol_norm, tol_neg_in
        logical,               intent(in) :: lowdin
        real(dp) :: tol_neg
        character(len=8) :: mname
        integer :: ispin, ik, ib, g, n_tot, n_part, n_neg
        real(dp) :: dev, max_tot, max_part, wmin, w

        if (lowdin) then
            tol_neg = tol_norm;   mname = 'Lowdin'
        else
            tol_neg = tol_neg_in; mname = 'Mulliken'
        end if
        n_tot = 0; n_part = 0; n_neg = 0
        max_tot = 0.0_dp; max_part = 0.0_dp; wmin = huge(1.0_dp)
        do ispin = 1, sd%nspin
            do ik = 1, sd%nk
                do ib = 1, sd%nband
                    dev = abs(sd%weight(ib, ik, 1, ispin) - 1.0_dp)
                    max_tot = max(max_tot, dev)
                    if (dev > tol_norm) then
                        n_tot = n_tot + 1
                        if (n_tot <= MAX_WARN_LINES) print '(A,I0,A,I0,A,I0,A,ES10.3)', &
                            '[spectral warning] total weight /= 1: spin=', ispin, ' k=', ik, &
                            ' band=', ib, ' |w-1|=', dev
                    end if
                    if (partition .and. sd%nregion > 1) then
                        dev = abs(sum(sd%weight(ib, ik, 2:sd%nregion, ispin)) - 1.0_dp)
                        max_part = max(max_part, dev)
                        if (dev > tol_norm) then
                            n_part = n_part + 1
                            if (n_part <= MAX_WARN_LINES) print '(A,I0,A,I0,A,I0,A,ES10.3)', &
                                '[spectral warning] sum over regions /= 1: spin=', ispin, ' k=', ik, &
                                ' band=', ib, ' |sum-1|=', dev
                        end if
                    end if
                    do g = 1, sd%nregion
                        w = sd%weight(ib, ik, g, ispin)
                        wmin = min(wmin, w)
                        if (w < -tol_neg) then
                            n_neg = n_neg + 1
                            if (n_neg <= MAX_WARN_LINES) print '(A,I0,A,I0,A,I0,A,A,A,ES10.3)', &
                                '[spectral warning] negative ' // trim(mname) // ' weight: spin=', ispin, ' k=', ik, &
                                ' band=', ib, ' region=', trim(sd%region_name(g)), ' w=', w
                        end if
                    end do
                end do
            end do
        end do
        if (n_tot > MAX_WARN_LINES) print '(A,I0,A)', '[spectral warning] ... ', &
            n_tot - MAX_WARN_LINES, ' more total-normalisation warnings'
        if (n_part > MAX_WARN_LINES) print '(A,I0,A)', '[spectral warning] ... ', &
            n_part - MAX_WARN_LINES, ' more region-sum warnings'
        if (n_neg > MAX_WARN_LINES) print '(A,I0,A)', '[spectral warning] ... ', &
            n_neg - MAX_WARN_LINES, ' more negative-weight warnings'

        print '(A)', 'Normalisation checks:'
        print '(A,ES10.3,A,I0,A,ES9.2,A)', '  max |total - 1|        = ', max_tot, &
            '  (', n_tot, ' above tol_norm=', tol_norm, ')'
        if (partition .and. sd%nregion > 1) then
            print '(A,ES10.3,A,I0,A)', '  max |sum regions - 1|  = ', max_part, &
                '  (', n_part, ' above tol_norm; regions form a partition)'
        else if (sd%nregion > 1) then
            print '(A)', '  sum over regions       : not checked (regions are not a partition)'
        end if
        if (lowdin) then
            print '(A,ES10.3,A,I0,A,ES9.2,A)', '  min Lowdin weight      = ', wmin, &
                '  (', n_neg, ' below -tol_norm=', -tol_neg, ')'
        else
            print '(A,ES10.3,A,I0,A,ES9.2,A)', '  min Mulliken weight    = ', wmin, &
                '  (', n_neg, ' below -tol_negative_warn=', -tol_neg, ')'
        end if
    end subroutine check_weights

    subroutine compare_hs(A, B, nameA, nameB, ierr, errmsg)
        type(HSDataFormat), intent(in)  :: A, B
        character(len=*),   intent(in)  :: nameA, nameB
        integer,            intent(out) :: ierr
        character(len=*),   intent(out) :: errmsg
        character(len=:), allocatable :: what

        ierr = 0
        errmsg = ''
        what = ''
        if (A%nAtoms /= B%nAtoms) then
            what = 'nAtoms'
        else if (A%NOrbs /= B%NOrbs) then
            what = 'NOrbs'
        else if (A%Neighbors /= B%Neighbors) then
            what = 'Neighbors'
        else if (allocated(A%Orbs) .neqv. allocated(B%Orbs)) then
            what = 'Orbs'
        else if (allocated(A%IndexList) .neqv. allocated(B%IndexList)) then
            what = 'IndexList'
        end if
        if (len(what) == 0 .and. allocated(A%Orbs)) then
            if (size(A%Orbs) /= size(B%Orbs)) then
                what = 'Orbs'
            else if (any(A%Orbs /= B%Orbs)) then
                what = 'Orbs'
            end if
        end if
        if (len(what) == 0 .and. allocated(A%IndexList)) then
            if (any(shape(A%IndexList) /= shape(B%IndexList))) then
                what = 'IndexList'
            else if (any(A%IndexList /= B%IndexList)) then
                what = 'IndexList'
            end if
        end if
        if (len(what) > 0) then
            ierr = 1
            errmsg = "HSIO files disagree in " // what // ": '" // nameA // "' vs '" // nameB // "'"
        end if
    end subroutine compare_hs

    logical function file_exists(fname)
        character(len=*), intent(in) :: fname
        inquire(file=trim(fname), exist=file_exists)
    end function file_exists

    pure real(dp) function elapsed(c0, c1, rate)
        integer(int64), intent(in) :: c0, c1, rate
        elapsed = real(c1 - c0, dp) / real(max(rate, 1_int64), dp)
    end function elapsed

    pure function lower(s) result(r)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: r
        integer :: i, ic
        r = s
        do i = 1, len(s)
            ic = iachar(s(i:i))
            if (ic >= iachar('A') .and. ic <= iachar('Z')) r(i:i) = achar(ic + 32)
        end do
    end function lower

    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=24) :: buf
        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

end module spectral_engine
