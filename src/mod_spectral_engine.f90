! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Spectral engine (docs/design/contracts.md §3, §3.2; uses §4 and §5): reads
!> spectral.nml, the HSIO H/S files, the structure file and regions.nml,
!> diagonalizes H(k) along a k-path and computes per-region projection weights
!> (Mulliken or Löwdin, selected by `method`).
!>
!> Pipeline position: backend of the `spectral` executable. Output is a
!> spectral_data_t (module spectral_data), written to disk by the driver and
!> later broadened/plotted by specmap.
!>
!> Conventions:
!>   - Energies in eV, k-points in fractional reciprocal coordinates, all
!>     indices 1-based; atom order = structure file = HSIO order.
!>   - Region group 1 is always 'total' (all atoms); user regions are groups
!>     2..nregion+1 (so sd%nregion = 1 + number of regions in regions.nml).
!>   - Array layout: sd%energy(nband, nk, nspin),
!>     sd%weight(nband, nk, group, nspin). nspin = 2 only for collinear.
!>
!> The per-k work lives in `process_kpoint`, which only touches the arrays and
!> workspaces passed to it. The (spin, k) loop of `run_spectral` is
!> parallelized with OpenMP (when compiled with -fopenmp): each thread owns
!> one set of workspaces, allocated once, and each iteration writes only its
!> own sd%energy / sd%weight slices, so results do not depend on the thread
!> count. Without -fopenmp the `!$` lines vanish and the loop runs serially.
!> Library routines never stop: errors are returned through ierr/errmsg.
!> (Exceptions: HSDataFormat%ReadFile stops on a corrupt HSIO file; missing
!> files are detected here before calling it. process_kpoint does `error stop`
!> on a programming error: Löwdin requested without its required optional
!> arguments.)
module spectral_engine
    use, intrinsic :: iso_fortran_env, only: int64, iostat_end, output_unit
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
    !$ use omp_lib, only: omp_get_num_threads
    use spectral_data, only: spectral_data_t, SPIN_NONE, SPIN_COLLINEAR, SPIN_NONCOLLINEAR, &
                             METHOD_MULLIKEN, METHOD_LOWDIN, EREF_UNKNOWN, EREF_FERMI, EREF_ABSOLUTE
    implicit none
    private

    ! ==================================================================
    ! Parameters and types
    ! ==================================================================

    integer, parameter, public :: MAX_KVERTEX = 100
    integer, parameter, public :: SPEC_PATH_LEN = 512
    !> Cap on printed warnings per category in check_weights.
    integer, parameter :: MAX_WARN_LINES = 20

    !> Mirror of the &SPECTRAL and &KPATH namelists (contracts §3). The
    !> defaults here are the defaults of the input file. Fields are plain
    !> strings/numbers; validation happens in read_spectral_input.
    !> vertex(:, i) = i-th k-path vertex (fractional); npoints(i) = points
    !> of the segment from vertex i to i+1 (exact meaning: see module kpath).
    type, public :: spectral_input_t
        character(len=SPEC_PATH_LEN) :: h_file = ''
        character(len=SPEC_PATH_LEN) :: h_up_file = ''
        character(len=SPEC_PATH_LEN) :: h_dn_file = ''
        character(len=SPEC_PATH_LEN) :: s_file = ''
        character(len=SPEC_PATH_LEN) :: structure_file = ''
        character(len=SPEC_PATH_LEN) :: regions_file = ''
        logical :: allow_no_fingerprint = .false.   ! accept regions.nml without structure fingerprint
        character(len=16) :: spin_mode = 'none'
        character(len=16) :: method = 'mulliken'
        character(len=16) :: energy_reference = 'unknown'
        real(dp) :: fermi_energy_ev = 0.0_dp
        logical :: fermi_known = .false.
        character(len=SPEC_PATH_LEN) :: output_prefix = 'slab'
        logical :: write_text = .false.
        real(dp) :: tol_norm = 1.0e-8_dp             ! normalisation-check tolerance (check_weights)
        real(dp) :: tol_negative_warn = 1.0e-3_dp    ! Mulliken negative-weight warning threshold
        real(dp) :: tol_s_min = 1.0e-10_dp           ! Löwdin: min eig(S) must exceed this
        ! &KPATH
        integer :: nvertex = 0
        real(dp) :: vertex(3, MAX_KVERTEX) = 0.0_dp
        character(len=16) :: label(MAX_KVERTEX) = ''
        integer :: npoints(MAX_KVERTEX - 1) = 0
    end type spectral_input_t

    !> All values are measured with system_clock (wall clock), in seconds.
    !> build_s, diag_s and proj_s are ACCUMULATED over all threads of the
    !> k-loop (thread-seconds: sum of each thread's wall time in that phase);
    !> loop_wall_s is the wall-clock time of the whole k-loop. With one thread
    !> build_s + diag_s + proj_s <= loop_wall_s.
    type, public :: spectral_timing_t
        real(dp) :: read_s  = 0.0_dp   ! structure, regions, HSIO (wall)
        real(dp) :: build_s = 0.0_dp   ! H(k), S(k) construction (accumulated thread-s)
        real(dp) :: diag_s  = 0.0_dp   ! generalized eigensolver (accumulated thread-s)
        real(dp) :: proj_s  = 0.0_dp   ! S*C + weights (accumulated thread-s)
        real(dp) :: loop_wall_s = 0.0_dp ! k-loop wall-clock time
        real(dp) :: write_s = 0.0_dp   ! filled by the driver (wall)
        real(dp) :: total_s = 0.0_dp   ! wall
        integer  :: nthreads = 1       ! threads that ran the k-loop (1 without OpenMP)
    end type spectral_timing_t

    public :: read_spectral_input, run_spectral, process_kpoint

contains

    ! ==================================================================
    ! Input
    ! ==================================================================

    !> Read and validate the &SPECTRAL and &KPATH namelists from `filename`
    !> (the two groups may appear in any order in the file).
    !> inp is initialised to its defaults first, so keys absent from the file
    !> keep them. spin_mode / method / energy_reference are lower-cased.
    !> energy_reference='fermi' forces fermi_known=.true. and E_F = 0.
    !> Never stops: on failure ierr = 1 and errmsg describes the first problem
    !> found (inp is then only partially filled).
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
        character(len=64) :: label(MAX_KVERTEX)   ! wider than stored (16) to detect truncation
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

        ! 1. Namelist reads only overwrite variables present in the file, so
        ! seed the local namelist variables with the defaults of spectral_input_t
        ! (inp is intent(out): its default initialisation applies on entry).
        h_file = inp%h_file;             h_up_file = inp%h_up_file
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

        ! 2. Open the file and read both namelists (rewind: order-independent).
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

        ! 3. Copy back into inp (normalising enum strings) and validate.
        inp%h_file = h_file;             inp%h_up_file = h_up_file
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
        if (any(len_trim(label) > len(inp%label))) then
            call fail('&KPATH label longer than 16 characters')
            return
        end if
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
            ! Energies are already relative to E_F: E_F = 0 by definition.
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

    !> Full calculation: read inputs, check consistency, build orbital groups
    !> and the k-path, diagonalize H(k) at every (spin, k) and fill `sd`.
    !>   inp           validated input (see read_spectral_input)
    !>   sd            output spectral data; sd%energy(nband,nk,nspin) in eV,
    !>                 sd%weight(nband,nk,group,nspin); group 1 = 'total'
    !>   timing        wall/thread timings; write_s is left for the driver
    !>   show_progress optional (default .false.): print [1/3]..[3/3] stage
    !>                 lines and k-loop progress/ETA to stdout
    !> Never stops: on failure ierr = 1 and errmsg holds the message (for k-loop
    !> failures it names ispin, ik and kfrac of the lowest failing iteration).
    !> Warnings from check_weights are printed, not returned. Exception:
    !> HSDataFormat%ReadFile may stop on a corrupt HSIO file.
    subroutine run_spectral(inp, sd, timing, ierr, errmsg, show_progress)
        type(spectral_input_t),  intent(in)  :: inp
        type(spectral_data_t),   intent(out) :: sd
        type(spectral_timing_t), intent(out) :: timing
        integer,                 intent(out) :: ierr
        character(len=*),        intent(out) :: errmsg
        logical, optional,       intent(in)  :: show_progress

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
        ! k-loop (OpenMP) state
        integer :: nk, nit, it, first_fail, ff, nthr
        integer, allocatable :: it_kerr(:), it_perr(:)
        real(dp), allocatable :: it_smin(:)
        real(dp) :: t_b, t_d, t_p
        integer(int64) :: cl0, cl1, cp
        integer :: completed, next_report, report_step
        logical :: report_progress
        real(dp) :: progress_elapsed, progress_eta

        ierr = 0
        errmsg = ''
        call system_clock(c_start, crate)
        use_lowdin = trim(inp%method) == 'lowdin'
        report_progress = .false.
        if (present(show_progress)) report_progress = show_progress

        if (report_progress) then
            print '(A)', '[1/3] Loading and validating input data'
            flush(output_unit)
        end if

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
        ! Group 1 = all atoms ('total'); groups 2.. = user regions, in file
        ! order. Atom lists are concatenated in CSR form (list_start/list_atoms).
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
        ! partition: user regions cover every atom exactly once, so their
        ! weights must sum to 1 (verified in check_weights).
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
        ! Fill everything in sd except the k-loop results (energy/weight are
        ! allocated zeroed here and written in step 6).
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

        ! ---------------- 6. k loop (OpenMP over the flat (ispin, ik) index) --
        ! it = (ispin-1)*nk + ik (inverted below as ispin = (it-1)/nk + 1).
        ! t_b/t_d/t_p are OpenMP reductions (thread-seconds). Each iteration writes only sd%energy(:,ik,ispin),
        ! sd%weight(:,ik,:,ispin) and its own entry of the status arrays, so the
        ! result is independent of the thread count and schedule. Workspaces
        ! are private and allocated once per thread. Progress I/O, when enabled,
        ! is serialized; there is no stop inside the region. The first failing
        ! iteration (lowest it) is reported after the parallel region.
        ! Iterations above the lowest failure seen so far are skipped; those
        ! below it always run, so the reported failure is deterministic.
        nk = path%nk
        nit = nspin * nk
        allocate(it_kerr(nit), it_perr(nit), it_smin(nit))
        it_kerr = 0
        it_perr = 0
        it_smin = 0.0_dp
        first_fail = nit + 1
        nthr = 1
        completed = 0
        report_step = max(1, (nit + 9) / 10)
        next_report = report_step
        t_b = 0.0_dp; t_d = 0.0_dp; t_p = 0.0_dp
        if (report_progress) then
            print '(A,I0,A)', '[2/3] Solving ', nit, ' (spin, k) points'
        end if
        call system_clock(cl0)

        !$omp parallel default(none) &
        !$omp shared(n, ngroup, nk, nit, hdat, sdat, path, groups, use_lowdin, inp, sd, &
        !$omp        it_kerr, it_perr, it_smin, first_fail, nthr, report_progress, completed, &
        !$omp        next_report, report_step, cl0, crate, cp, progress_elapsed, progress_eta) &
        !$omp private(it, ispin, ik, g, ff, kerr, perr, smin, Hk, Sk, C, SC, E, W) &
        !$omp reduction(+: t_b, t_d, t_p)
        allocate(Hk(n, n), Sk(n, n), C(n, n), SC(n, n), E(n), W(n, ngroup))
        !$omp single
        !$ nthr = omp_get_num_threads()
        !$omp end single
        !$omp do schedule(dynamic, 1)
        do it = 1, nit
            !$omp atomic read
            ff = first_fail
            if (it > ff) cycle
            ispin = (it - 1) / nk + 1
            ik = it - (ispin - 1) * nk
            call process_kpoint(n, ngroup, hdat(ispin)%Data, sdat%Data, sdat%IndexList, &
                                path%kfrac(:, ik), groups, Hk, Sk, C, SC, E, W, &
                                t_b, t_d, t_p, kerr, &
                                lowdin=use_lowdin, tol_s_min=inp%tol_s_min, smin=smin, &
                                proj_ierr=perr)
            if (kerr /= 0 .or. perr /= 0) then
                it_kerr(it) = kerr
                it_perr(it) = perr
                it_smin(it) = smin
                !$omp atomic update
                first_fail = min(first_fail, it)
                cycle
            end if
            sd%energy(:, ik, ispin) = E
            do g = 1, ngroup
                sd%weight(:, ik, g, ispin) = W(:, g)
            end do
            if (report_progress) then
                !$omp critical (spectral_progress)
                completed = completed + 1
                if (completed >= next_report .or. completed == nit) then
                    call system_clock(cp)
                    progress_elapsed = elapsed(cl0, cp, crate)
                    progress_eta = progress_elapsed * real(nit - completed, dp) / real(completed, dp)
                    write(*, '(A,I0,A,I0,A,F6.1,A,F9.2,A,F9.2,A)') '      ', completed, '/', nit, &
                        '  ', 100.0_dp * real(completed, dp) / real(nit, dp), &
                        '%  elapsed ', progress_elapsed, ' s  ETA ', progress_eta, ' s'
                    flush(output_unit)
                    next_report = min(nit, next_report + report_step)
                end if
                !$omp end critical (spectral_progress)
            end if
        end do
        !$omp end do
        deallocate(Hk, Sk, C, SC, E, W)
        !$omp end parallel

        call system_clock(cl1)
        timing%loop_wall_s = elapsed(cl0, cl1, crate)
        timing%build_s = t_b
        timing%diag_s  = t_d
        timing%proj_s  = t_p
        timing%nthreads = nthr

        ! Report the lowest failing iteration (deterministic across thread counts).
        ! perr == 1 is the Löwdin S-not-positive-definite code of lowdin_weights.
        if (first_fail <= nit) then
            it = first_fail
            ispin = (it - 1) / nk + 1
            ik = it - (ispin - 1) * nk
            perr = it_perr(it)
            kerr = it_kerr(it)
            smin = it_smin(it)
            if (perr == 1) then
                write(msg, '(A,I0,A,I0,A,3F12.7,A,ES12.4,A,ES10.3,A)') &
                    'overlap matrix not positive definite above tol_s_min at ispin=', ispin, &
                    ' ik=', ik, ' kfrac=(', path%kfrac(:, ik), '): min eigenvalue smin=', smin, &
                    ' (tol_s_min=', inp%tol_s_min, ')'
            else if (perr /= 0) then
                write(msg, '(A,I0,A,I0,A,I0,A,3F12.7,A)') 'Lowdin projection failed (ierr=', &
                    perr, ') at ispin=', ispin, ' ik=', ik, ' kfrac=(', path%kfrac(:, ik), ')'
            else if (kerr > n) then
                ! LAPACK zhegv convention: info > n means the leading minor of S(k) of
                ! order info-n is not positive definite.
                write(msg, '(A,I0,A,I0,A,I0,A,3F12.7,A)') &
                    'overlap S(k) not positive definite (zhegv info=', &
                    kerr, ') at ispin=', ispin, ' ik=', ik, ' kfrac=(', path%kfrac(:, ik), ')'
            else
                write(msg, '(A,I0,A,I0,A,I0,A,3F12.7,A)') 'generalized eigensolver failed (info=', &
                    kerr, ') at ispin=', ispin, ' ik=', ik, ' kfrac=(', path%kfrac(:, ik), ')'
            end if
            call fail(trim(msg))
            return
        end if

        ! ---------------- 7. checks -------------------------------------
        if (report_progress) print '(A)', '[3/3] Checking projection weights'
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

    !> One (spin, k) point: build H(k), S(k), solve H C = S C E, and compute the
    !> projection weights W(band, group).
    !>   Hr, Sr, Rpt  real-space H(R), S(R) and lattice vectors n_R (integer
    !>                triples), as stored in HSDataFormat (Data / IndexList);
    !>                Fourier sum H(k) = sum_R H(R) exp(+i 2 pi k.n_R)
    !>                (core::compute_Hk_at_kpoint, lattice gauge)
    !>   kfrac        k-point in fractional reciprocal coordinates
    !>   Hk,Sk,C,SC   caller-owned (n,n) workspaces (one set per thread);
    !>                on return Sk = original S(k), C = eigenvectors (columns)
    !>   E(n), W(n,ngroup)  eigenvalues in eV (ascending) and weights
    !>   t_build/t_diag/t_proj  wall seconds ADDED to the passed values
    !>   ierr         info of the generalized eigensolver (0 = ok); on ierr /= 0
    !>                E, C, W are undefined and t_proj is not updated
    !> Optional Löwdin mode (lowdin=.true.): weights from lowdin_weights on the
    !> ORIGINAL Sk (SC is not used); tol_s_min, smin and proj_ierr are then
    !> required. proj_ierr = 1 means min eig(Sk) = smin <= tol_s_min.
    !> Without the optional arguments the Mulliken path is used (unchanged).
    !> Violating the Löwdin argument requirement is a programming error and
    !> triggers `error stop` (the only stop in this module's library API).
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

        ! 1. H(k), S(k) by Fourier sum over R.
        call system_clock(c0, crate)
        call compute_Hk_at_kpoint(Hr, Rpt, Hk, kfrac)
        call compute_Hk_at_kpoint(Sr, Rpt, Sk, kfrac)
        call system_clock(c1)
        ! 2. Generalized eigenproblem. generalized_eigensolve works on internal
        ! copies: Sk stays the original S(k), as the projections require.
        call generalized_eigensolve(n, Hk, Sk, E, C, ierr)
        call system_clock(c2)
        if (ierr /= 0) then
            t_build = t_build + elapsed(c0, c1, crate)
            t_diag  = t_diag  + elapsed(c1, c2, crate)
            return
        end if
        ! 3. Projection weights (Löwdin from S^{1/2}C, Mulliken from C and S·C).
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

    !> Diagnostics on the finished spectral data; prints warnings (capped at
    !> MAX_WARN_LINES per category) and a summary, never fails. Checks:
    !>   - group 1 ('total') weight = 1 per band within tol_norm
    !>   - if `partition`, sum over user regions (groups 2..) = 1 within tol_norm
    !>   - no weight below -tol_neg (Mulliken populations may be slightly negative)
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

    !> Check that two HSIO files (e.g. H and S) describe the same system:
    !> nAtoms, NOrbs, Neighbors, Orbs and IndexList (R-vector list) must agree.
    !> ierr = 1 and errmsg names the first mismatching item.
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

    !> Seconds between two system_clock counts (rate guarded against 0).
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
