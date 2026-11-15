! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Tests for module spectral_engine (contracts §3, §3.2, §4, §5) on the
!> example HSIO fixtures (Diatomic and Spinor non-collinear chains).
program test_spectral_engine
    use types, only: dp
    use hs_io, only: HSDataFormat
    use core, only: compute_Hk_at_kpoint
    use structure_io, only: structure_t, read_structure, structure_fingerprint, species_norb_of
    use atomic_regions, only: region_set_t, add_region, write_regions_file
    use spectral_projection, only: orbital_groups_t, build_orbital_groups
    use spectral_data, only: spectral_data_t, write_spectral_file, read_spectral_file, &
                             write_spectral_text, find_region, &
                             SPIN_NONE, SPIN_COLLINEAR, SPIN_NONCOLLINEAR, METHOD_MULLIKEN, METHOD_LOWDIN, &
                             EREF_FERMI, EREF_UNKNOWN
    use spectral_engine, only: spectral_input_t, spectral_timing_t, read_spectral_input, run_spectral
    use testing, only: check, check_close, test_summary
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    implicit none

    interface
        subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
            use types, only: dp
            integer,          intent(in)    :: itype, n, lda, ldb, lwork
            character(len=1), intent(in)    :: jobz, uplo
            complex(dp),      intent(inout) :: a(lda, *), b(ldb, *)
            real(dp),         intent(out)   :: w(*), rwork(*)
            complex(dp),      intent(out)   :: work(*)
            integer,          intent(out)   :: info
        end subroutine zhegv
    end interface

    character(len=*), parameter :: DIA_H = 'examples/UnfoldDiatomic/Diatomic_H.hsio'
    character(len=*), parameter :: DIA_S = 'examples/UnfoldDiatomic/Diatomic_S.hsio'
    character(len=*), parameter :: DIA_ST = 'examples/UnfoldDiatomic/sc.nml'
    character(len=*), parameter :: SPI_H = 'examples/UnfoldSpinorNonCollinear/Spinor_H.hsio'
    character(len=*), parameter :: SPI_S = 'examples/UnfoldSpinorNonCollinear/Spinor_S.hsio'
    character(len=*), parameter :: SPI_ST = 'examples/UnfoldSpinorNonCollinear/sc.nml'
    character(len=*), parameter :: TMP = 'build/tse_'

    call test_read_input()
    call test_consistency_diatomic()
    call test_lowdin_diatomic()
    call test_regions_partition()
    call test_spin_noncollinear()
    call test_spin_collinear()
    call test_errors()
    call test_roundtrip_file()
    call test_determinism()

    call test_summary("test_spectral_engine")

contains

    ! ------------------------------------------------------------------
    ! Helpers
    ! ------------------------------------------------------------------

    subroutine base_input(inp, h, s, st)
        type(spectral_input_t), intent(out) :: inp
        character(len=*),       intent(in)  :: h, s, st
        inp%h_file = h
        inp%s_file = s
        inp%structure_file = st
        inp%nvertex = 2
        inp%vertex(:, 1) = [-0.5_dp, 0.0_dp, 0.0_dp]
        inp%vertex(:, 2) = [ 0.5_dp, 0.0_dp, 0.0_dp]
        inp%label(1:2) = ['X-', 'X+']
        inp%npoints(1) = 20
    end subroutine base_input

    subroutine write_lines(fname, lines)
        character(len=*), intent(in) :: fname
        character(len=*), intent(in) :: lines(:)
        integer :: iu, i
        open(newunit=iu, file=fname, status='replace', action='write')
        do i = 1, size(lines)
            write(iu, '(A)') trim(lines(i))
        end do
        close(iu)
    end subroutine write_lines

    subroutine remove_file(fname)
        character(len=*), intent(in) :: fname
        integer :: iu, ios
        logical :: ex
        inquire(file=fname, exist=ex)
        if (.not. ex) return
        open(newunit=iu, file=fname, status='old', iostat=ios)
        if (ios == 0) close(iu, status='delete')
    end subroutine remove_file

    !> Writes a regions file for the structure in st_file with the given
    !> CSR lists; fingerprint: 'auto' = correct one, else used verbatim.
    subroutine make_regions(fname, st_file, names, start, atoms, fp)
        character(len=*), intent(in) :: fname, st_file
        character(len=*), intent(in) :: names(:)
        integer,          intent(in) :: start(:), atoms(:)
        character(len=*), intent(in) :: fp
        type(structure_t) :: st
        type(region_set_t) :: rs
        integer :: ierr, i
        character(len=256) :: msg
        call read_structure(st_file, st, ierr, msg)
        rs%natoms = st%nat
        if (fp == 'auto') then
            rs%fingerprint = structure_fingerprint(st)
        else
            rs%fingerprint = fp
        end if
        do i = 1, size(names)
            call add_region(rs, trim(names(i)), atoms(start(i):start(i + 1) - 1), ierr, msg)
            call check(ierr == 0, 'add_region ' // trim(names(i)) // ': ' // trim(msg))
        end do
        call write_regions_file(fname, rs, ierr, msg)
        call check(ierr == 0, 'write_regions_file: ' // trim(msg))
    end subroutine make_regions

    real(dp) function max_total_dev(sd)
        type(spectral_data_t), intent(in) :: sd
        max_total_dev = maxval(abs(sd%weight(:, :, 1, :) - 1.0_dp))
    end function max_total_dev

    real(dp) function max_regsum_dev(sd)
        type(spectral_data_t), intent(in) :: sd
        max_regsum_dev = maxval(abs(sum(sd%weight(:, :, 2:sd%nregion, :), dim=3) - 1.0_dp))
    end function max_regsum_dev

    ! ------------------------------------------------------------------
    ! read_spectral_input
    ! ------------------------------------------------------------------

    subroutine test_read_input()
        type(spectral_input_t) :: inp
        integer :: ierr
        character(len=512) :: msg
        character(len=*), parameter :: f = TMP // 'input.nml'

        call write_lines(f, [character(len=120) :: &
            '&KPATH', &
            '  nvertex = 3', &
            '  vertex(:,1) = 0.3333333333, 0.3333333333, 0.0', &
            '  vertex(:,2) = 0.0, 0.0, 0.0', &
            '  vertex(:,3) = 0.5, 0.0, 0.0', &
            "  label = 'K', 'G', 'M'", &
            '  npoints = 40, 30', &
            '/', &
            '&SPECTRAL', &
            "  h_file = 'h.hsio'", &
            "  s_file = 's.hsio'", &
            "  structure_file = 'st.nml'", &
            "  spin_mode = 'NoneCollinear'", &
            '/'])
        call read_spectral_input(f, inp, ierr, msg)
        call check(ierr /= 0, 'invalid spin_mode rejected')
        call check(index(msg, 'spin_mode') > 0, 'invalid spin_mode message: ' // trim(msg))

        call write_lines(f, [character(len=120) :: &
            '&KPATH', &
            '  nvertex = 3', &
            '  vertex(:,1) = 0.3333333333, 0.3333333333, 0.0', &
            '  vertex(:,2) = 0.0, 0.0, 0.0', &
            '  vertex(:,3) = 0.5, 0.0, 0.0', &
            "  label = 'K', 'G', 'M'", &
            '  npoints = 40, 30', &
            '/', &
            '&SPECTRAL', &
            "  h_file = 'h.hsio'", &
            "  s_file = 's.hsio'", &
            "  structure_file = 'st.nml'", &
            "  spin_mode = 'noncollinear'", &
            "  energy_reference = 'fermi'", &
            '  fermi_energy_ev = 3.5', &
            '  tol_norm = 1e-6', &
            '/'])
        call read_spectral_input(f, inp, ierr, msg)
        call check(ierr == 0, 'read_spectral_input (KPATH first): ' // trim(msg))
        call check(inp%nvertex == 3, 'nvertex = 3')
        call check(all(inp%npoints(1:2) == [40, 30]), 'npoints')
        call check_close(inp%vertex(1, 1), 0.3333333333_dp, 1e-12_dp, 'vertex(1,1)')
        call check(trim(inp%label(2)) == 'G', 'label(2)')
        call check(trim(inp%h_file) == 'h.hsio', 'h_file')
        call check(trim(inp%spin_mode) == 'noncollinear', 'spin_mode')
        call check(trim(inp%method) == 'mulliken', 'default method')
        call check(inp%fermi_known, 'fermi forces fermi_known')
        call check_close(inp%fermi_energy_ev, 0.0_dp, 0.0_dp, 'fermi forces E_F = 0')
        call check_close(inp%tol_norm, 1e-6_dp, 1e-20_dp, 'tol_norm read')
        call check_close(inp%tol_negative_warn, 1e-3_dp, 1e-20_dp, 'tol_negative_warn default')
        call check(.not. inp%allow_no_fingerprint, 'allow_no_fingerprint default')
        call check(trim(inp%output_prefix) == 'slab', 'output_prefix default')

        ! collinear without h_dn_file
        call write_lines(f, [character(len=120) :: &
            '&SPECTRAL', "  h_up_file = 'u.hsio'", "  s_file = 's.hsio'", &
            "  structure_file = 'st.nml'", "  spin_mode = 'collinear'", '/', &
            '&KPATH', '  nvertex = 2', '  npoints = 4', '/'])
        call read_spectral_input(f, inp, ierr, msg)
        call check(ierr /= 0, 'collinear without h_dn_file rejected')

        ! nvertex out of range
        call write_lines(f, [character(len=120) :: &
            '&SPECTRAL', "  h_file = 'h.hsio'", "  s_file = 's.hsio'", &
            "  structure_file = 'st.nml'", '/', &
            '&KPATH', '  nvertex = 1', '/'])
        call read_spectral_input(f, inp, ierr, msg)
        call check(ierr /= 0, 'nvertex = 1 rejected')

        ! lowdin accepted by the parser; bad tolerance rejected
        call write_lines(f, [character(len=120) :: &
            '&SPECTRAL', "  h_file = 'h.hsio'", "  s_file = 's.hsio'", &
            "  structure_file = 'st.nml'", "  method = 'lowdin'", '/', &
            '&KPATH', '  nvertex = 2', '  npoints = 4', '/'])
        call read_spectral_input(f, inp, ierr, msg)
        call check(ierr == 0, 'lowdin accepted by parser: ' // trim(msg))
        call write_lines(f, [character(len=120) :: &
            '&SPECTRAL', "  h_file = 'h.hsio'", "  s_file = 's.hsio'", &
            "  structure_file = 'st.nml'", "  tol_norm = 0.0", '/', &
            '&KPATH', '  nvertex = 2', '  npoints = 4', '/'])
        call read_spectral_input(f, inp, ierr, msg)
        call check(ierr /= 0, 'tol_norm = 0 rejected')

        call remove_file(f)
    end subroutine test_read_input

    ! ------------------------------------------------------------------
    ! Test 7: energies vs independent zhegv('N'), total = 1
    ! ------------------------------------------------------------------

    subroutine test_consistency_diatomic()
        type(spectral_input_t) :: inp
        type(spectral_data_t) :: sd
        type(spectral_timing_t) :: tm
        type(HSDataFormat) :: hd, sdat
        complex(dp), allocatable :: Hk(:,:), Sk(:,:), work(:)
        real(dp), allocatable :: w(:), rwork(:)
        integer :: ierr, ik, n, info, lwork
        real(dp) :: maxdiff
        character(len=512) :: msg

        call base_input(inp, DIA_H, DIA_S, DIA_ST)
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr == 0, 'run_spectral Diatomic: ' // trim(msg))
        if (ierr /= 0) return
        call check(sd%nk == 21, 'nk = 21')
        call check(sd%nband == 9, 'nband = NOrbs = 9')
        call check(sd%nspin == 1, 'none: nspin = 1')
        call check(sd%spin_mode == SPIN_NONE, 'spin_mode none')
        call check(sd%method == METHOD_MULLIKEN, 'method mulliken')
        call check(sd%energy_ref == EREF_UNKNOWN, 'energy_ref unknown')
        call check(sd%nregion == 1, "only 'total' without regions file")
        call check(trim(sd%region_name(1)) == 'total', "region 1 = 'total'")
        call check(sd%region_natoms(1) == 6, 'total natoms = 6')
        call check(sd%nvertex == 2 .and. all(sd%vertex_k == [1, 21]), 'vertex_k')
        call check(len_trim(sd%fingerprint) == 16, 'fingerprint set')
        call check(index(sd%source, 'Diatomic_H.hsio') > 0, 'source names H file')
        call check(tm%total_s >= 0.0_dp .and. tm%diag_s >= 0.0_dp, 'timings non-negative')

        call hd%ReadFile(DIA_H)
        call sdat%ReadFile(DIA_S)
        n = hd%NOrbs
        lwork = max(1, 2 * n - 1)
        allocate(Hk(n, n), Sk(n, n), w(n), work(lwork), rwork(max(1, 3 * n - 2)))
        maxdiff = 0.0_dp
        do ik = 1, sd%nk
            call compute_Hk_at_kpoint(hd%Data, hd%IndexList, Hk, sd%kfrac(:, ik))
            call compute_Hk_at_kpoint(sdat%Data, sdat%IndexList, Sk, sd%kfrac(:, ik))
            call zhegv(1, 'N', 'U', n, Hk, n, Sk, n, w, work, lwork, rwork, info)
            call check(info == 0, 'reference zhegv info = 0')
            maxdiff = max(maxdiff, maxval(abs(w - sd%energy(:, ik, 1))))
        end do
        call check_close(maxdiff, 0.0_dp, 1e-10_dp, 'energies vs reference zhegv (max diff)')
        call check_close(max_total_dev(sd), 0.0_dp, 1e-10_dp, 'total weight = 1')
    end subroutine test_consistency_diatomic

    ! ------------------------------------------------------------------
    ! Lowdin on Diatomic: total = 1, partition = 1, energies = Mulliken run
    ! ------------------------------------------------------------------

    subroutine test_lowdin_diatomic()
        type(spectral_input_t) :: inp
        type(spectral_data_t) :: sdm, sdl
        type(spectral_timing_t) :: tm
        integer :: ierr
        character(len=512) :: msg
        character(len=*), parameter :: rf = TMP // 'regions_lowdin.nml'

        call make_regions(rf, DIA_ST, ['A_atoms', 'B_atoms'], [1, 4, 7], [1, 3, 5, 2, 4, 6], 'auto')
        call base_input(inp, DIA_H, DIA_S, DIA_ST)
        inp%regions_file = rf
        call run_spectral(inp, sdm, tm, ierr, msg)
        call check(ierr == 0, 'lowdin ref: Mulliken run: ' // trim(msg))
        inp%method = 'lowdin'
        call run_spectral(inp, sdl, tm, ierr, msg)
        call check(ierr == 0, 'run_spectral Diatomic lowdin: ' // trim(msg))
        if (ierr == 0) then
            call check(sdl%method == METHOD_LOWDIN, 'lowdin: sd%method = METHOD_LOWDIN')
            call check_close(max_total_dev(sdl), 0.0_dp, 1e-10_dp, 'lowdin: total weight = 1')
            call check_close(max_regsum_dev(sdl), 0.0_dp, 1e-10_dp, 'lowdin: partition sums to 1')
            call check(minval(sdl%weight) >= -1e-12_dp, 'lowdin: weights >= -1e-12')
            if (allocated(sdm%energy)) &
                call check(all(sdl%energy == sdm%energy), 'lowdin: energies identical to Mulliken run')
        end if
        call remove_file(rf)
    end subroutine test_lowdin_diatomic

    ! ------------------------------------------------------------------
    ! Test 2b: regions partition / overlap
    ! ------------------------------------------------------------------

    subroutine test_regions_partition()
        type(spectral_input_t) :: inp
        type(spectral_data_t) :: sd
        type(spectral_timing_t) :: tm
        integer :: ierr
        character(len=512) :: msg
        character(len=*), parameter :: rf = TMP // 'regions_ab.nml'

        call make_regions(rf, DIA_ST, ['A_atoms', 'B_atoms'], [1, 4, 7], [1, 3, 5, 2, 4, 6], 'auto')
        call base_input(inp, DIA_H, DIA_S, DIA_ST)
        inp%regions_file = rf
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr == 0, 'run_spectral with partition regions: ' // trim(msg))
        if (ierr == 0) then
            call check(sd%nregion == 3, 'nregion = total + 2')
            call check(find_region(sd, 'A_atoms') == 2 .and. find_region(sd, 'B_atoms') == 3, &
                       'region order = file order')
            call check(all(sd%region_natoms == [6, 3, 3]), 'region_natoms')
            call check_close(max_regsum_dev(sd), 0.0_dp, 1e-10_dp, 'sum over partition regions = 1')
            call check_close(max_total_dev(sd), 0.0_dp, 1e-10_dp, 'total = 1 (with regions)')
        end if

        ! overlapping regions are accepted
        call make_regions(rf, DIA_ST, ['left ', 'right'], [1, 4, 7], [1, 2, 3, 3, 4, 5], 'auto')
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr == 0, 'run_spectral with overlapping regions: ' // trim(msg))
        if (ierr == 0) call check(sd%nregion == 3, 'overlap: nregion = 3')

        call remove_file(rf)
    end subroutine test_regions_partition

    ! ------------------------------------------------------------------
    ! Test 8a: non-collinear spinor
    ! ------------------------------------------------------------------

    subroutine test_spin_noncollinear()
        type(spectral_input_t) :: inp
        type(spectral_data_t) :: sd
        type(spectral_timing_t) :: tm
        type(HSDataFormat) :: sdat
        type(structure_t) :: st
        type(orbital_groups_t) :: groups
        integer :: ierr, ia
        character(len=512) :: msg
        character(len=*), parameter :: rf = TMP // 'regions_atoms.nml'

        call make_regions(rf, SPI_ST, ['at1', 'at2', 'at3', 'at4', 'at5', 'at6'], &
                          [1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6], 'auto')
        call base_input(inp, SPI_H, SPI_S, SPI_ST)
        inp%spin_mode = 'noncollinear'
        inp%regions_file = rf
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr == 0, 'run_spectral noncollinear: ' // trim(msg))
        if (ierr == 0) then
            call check(sd%nspin == 1, 'noncollinear: nspin = 1')
            call check(sd%spin_mode == SPIN_NONCOLLINEAR, 'spin_mode noncollinear')
            call check(sd%nband == 18, 'spinor nband = 18')
            call check_close(max_total_dev(sd), 0.0_dp, 1e-10_dp, 'spinor total = 1')
            call check_close(max_regsum_dev(sd), 0.0_dp, 1e-10_dp, 'spinor per-atom regions sum = 1')
        end if

        ! orbital counts per atom = species table (4/2 alternating)
        call sdat%ReadFile(SPI_S)
        call read_structure(SPI_ST, st, ierr, msg)
        call build_orbital_groups(sdat%Orbs, 6, [1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6], groups)
        do ia = 1, 6
            call check(groups%start(ia + 1) - groups%start(ia) == sdat%Orbs(ia), 'group size = Orbs(a)')
            call check(sdat%Orbs(ia) == species_norb_of(st, ia), 'Orbs(a) = species nOrb')
        end do
        call check(all(groups%start(2:7) - groups%start(1:6) == [4, 2, 4, 2, 4, 2]), &
                   'spinor group sizes 4/2 alternating')

        call remove_file(rf)
    end subroutine test_spin_noncollinear

    ! ------------------------------------------------------------------
    ! Test 8b: collinear with identical up/dn
    ! ------------------------------------------------------------------

    subroutine test_spin_collinear()
        type(spectral_input_t) :: inp
        type(spectral_data_t) :: sd
        type(spectral_timing_t) :: tm
        integer :: ierr
        character(len=512) :: msg

        call base_input(inp, '', DIA_S, DIA_ST)
        inp%spin_mode = 'collinear'
        inp%h_up_file = DIA_H
        inp%h_dn_file = DIA_H
        inp%energy_reference = 'fermi'
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr == 0, 'run_spectral collinear: ' // trim(msg))
        if (ierr /= 0) return
        call check(sd%nspin == 2, 'collinear: nspin = 2')
        call check(sd%spin_mode == SPIN_COLLINEAR, 'spin_mode collinear')
        call check(sd%energy_ref == EREF_FERMI .and. sd%fermi_known, 'fermi reference metadata')
        call check(all(sd%energy(:, :, 1) == sd%energy(:, :, 2)), 'up/dn energies identical')
        call check(all(sd%weight(:, :, :, 1) == sd%weight(:, :, :, 2)), 'up/dn weights identical')
        call check(index(sd%source, 'H_up=') > 0, 'collinear source')
    end subroutine test_spin_collinear

    ! ------------------------------------------------------------------
    ! Error paths (ierr /= 0, no stop)
    ! ------------------------------------------------------------------

    subroutine test_errors()
        type(spectral_input_t) :: inp
        type(spectral_data_t) :: sd
        type(spectral_timing_t) :: tm
        integer :: ierr
        character(len=512) :: msg
        character(len=*), parameter :: sf = TMP // 'struct5.nml', rf = TMP // 'regions_bad.nml'

        ! 5-atom structure vs 6-atom HSIO
        call write_lines(sf, [character(len=80) :: &
            '&STRUCTURE', &
            '    LatticeVectors(1,:) = 3.0,0.0,0.0', &
            '    LatticeVectors(2,:) = 0.0,10.0,0.0', &
            '    LatticeVectors(3,:) = 0.0,0.0,10.0', &
            '/', '&ATOMS', '5', &
            'A 0.0 0.0 0.0', 'B 0.5 0.0 0.0', 'A 1.0 0.0 0.0', 'B 1.5 0.0 0.0', 'A 2.0 0.0 0.0', &
            '/'])
        call base_input(inp, DIA_H, DIA_S, sf)
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr /= 0, 'structure atom count mismatch rejected')
        call check(index(msg, 'nAtoms') > 0, 'atom count message: ' // trim(msg))

        ! wrong species nOrb (Spinor structure with Diatomic HSIO)
        call base_input(inp, DIA_H, DIA_S, SPI_ST)
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr /= 0, 'species nOrb mismatch rejected')
        call check(index(msg, 'atom 1') > 0, 'species mismatch names atom: ' // trim(msg))

        ! H/S mismatch
        call base_input(inp, SPI_H, DIA_S, DIA_ST)
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr /= 0, 'H vs S mismatch rejected')

        ! fingerprint mismatch
        call make_regions(rf, DIA_ST, ['r1'], [1, 3], [1, 2], '0123456789abcdef')
        call base_input(inp, DIA_H, DIA_S, DIA_ST)
        inp%regions_file = rf
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr /= 0, 'regions fingerprint mismatch rejected')

        ! missing fingerprint: rejected unless allowed
        call make_regions(rf, DIA_ST, ['r1'], [1, 3], [1, 2], '')
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr /= 0, 'missing fingerprint rejected without allow flag')
        inp%allow_no_fingerprint = .true.
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr == 0, 'missing fingerprint accepted with allow flag: ' // trim(msg))

        ! missing H file
        call base_input(inp, TMP // 'does_not_exist.hsio', DIA_S, DIA_ST)
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr /= 0, 'missing H file rejected')

        call remove_file(sf)
        call remove_file(rf)
    end subroutine test_errors

    ! ------------------------------------------------------------------
    ! run_spectral -> write -> read round trip
    ! ------------------------------------------------------------------

    subroutine test_roundtrip_file()
        type(spectral_input_t) :: inp
        type(spectral_data_t) :: sd, sd2
        type(spectral_timing_t) :: tm
        integer :: ierr
        character(len=512) :: msg
        character(len=*), parameter :: rf = TMP // 'regions_rt.nml', bf = TMP // 'rt.spectral.bin', &
                                       tf = TMP // 'rt.spectral.txt'

        call make_regions(rf, DIA_ST, ['A_atoms', 'B_atoms'], [1, 4, 7], [1, 3, 5, 2, 4, 6], 'auto')
        call base_input(inp, DIA_H, DIA_S, DIA_ST)
        inp%regions_file = rf
        call run_spectral(inp, sd, tm, ierr, msg)
        call check(ierr == 0, 'roundtrip run: ' // trim(msg))
        if (ierr /= 0) return
        call write_spectral_file(bf, sd, ierr, msg)
        call check(ierr == 0, 'write_spectral_file: ' // trim(msg))
        call write_spectral_text(tf, sd, ierr, msg)
        call check(ierr == 0, 'write_spectral_text: ' // trim(msg))
        call read_spectral_file(bf, sd2, ierr, msg)
        call check(ierr == 0, 'read_spectral_file: ' // trim(msg))
        if (ierr == 0) then
            call check(sd2%nk == sd%nk .and. sd2%nband == sd%nband .and. sd2%nregion == 3, 'dims')
            call check(all(sd2%energy == sd%energy), 'energy round trip')
            call check(all(sd2%weight == sd%weight), 'weight round trip')
            call check(sd2%fingerprint == sd%fingerprint, 'fingerprint round trip')
            call check(all(sd2%region_name == sd%region_name), 'region names round trip')
        end if
        call remove_file(rf)
        call remove_file(bf)
        call remove_file(tf)
    end subroutine test_roundtrip_file

    ! ------------------------------------------------------------------
    ! Determinism: repeated runs are bitwise identical; under OpenMP, 1 vs 4
    ! threads are bitwise identical as well (each (spin,k) is computed by one
    ! thread with the same serial code path; no reductions over k).
    ! ------------------------------------------------------------------

    subroutine test_determinism()
        type(spectral_input_t) :: inp
        type(spectral_data_t) :: sd1, sd2
        type(spectral_timing_t) :: tm
        integer :: ierr, max_thr
        character(len=512) :: msg

        max_thr = 1
        !$ max_thr = omp_get_max_threads()
        call base_input(inp, SPI_H, SPI_S, SPI_ST)
        inp%spin_mode = 'noncollinear'
        inp%npoints(1) = 60
        call run_spectral(inp, sd1, tm, ierr, msg)
        call check(ierr == 0, 'determinism: run 1: ' // trim(msg))
        call check(tm%nthreads >= 1, 'determinism: timing%nthreads >= 1')
        call check(tm%loop_wall_s >= 0.0_dp, 'determinism: loop_wall_s >= 0')
        call run_spectral(inp, sd2, tm, ierr, msg)
        call check(ierr == 0, 'determinism: run 2: ' // trim(msg))
        if (allocated(sd1%energy) .and. allocated(sd2%energy)) then
            call check(all(sd1%energy == sd2%energy), 'determinism: energies bitwise equal (repeat)')
            call check(all(sd1%weight == sd2%weight), 'determinism: weights bitwise equal (repeat)')
        end if

        ! OpenMP builds only: 1 thread vs 4 threads.
        !$ call omp_set_num_threads(1)
        !$ call run_spectral(inp, sd1, tm, ierr, msg)
        !$ call check(ierr == 0 .and. tm%nthreads == 1, 'determinism: 1-thread run: ' // trim(msg))
        !$ call omp_set_num_threads(4)
        !$ call run_spectral(inp, sd2, tm, ierr, msg)
        !$ call check(ierr == 0 .and. tm%nthreads == 4, 'determinism: 4-thread run: ' // trim(msg))
        !$ if (allocated(sd1%energy) .and. allocated(sd2%energy)) then
        !$     call check(all(sd1%energy == sd2%energy), 'determinism: energies bitwise equal (1 vs 4 threads)')
        !$     call check(all(sd1%weight == sd2%weight), 'determinism: weights bitwise equal (1 vs 4 threads)')
        !$ end if
        !$ call omp_set_num_threads(max_thr)
    end subroutine test_determinism

end program test_spectral_engine
