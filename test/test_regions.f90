! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Tests for atomic_regions (regions.nml) and atomic_selection (index lists).
program test_regions
    use types, only: dp
    use testing, only: check, test_summary
    use structure_io, only: structure_t, structure_fingerprint
    use atomic_regions, only: region_set_t, read_regions_file, write_regions_file, &
                              validate_regions, add_region, regions_are_partition, &
                              is_valid_region_name
    use atomic_selection, only: parse_index_list
    implicit none

    character(len=*), parameter :: TMPDIR = 'build/'
    type(structure_t) :: st
    character(len=16) :: fp
    integer :: ntmp = 0
    character(len=256) :: tmpfiles(64)

    call make_structure(st)
    fp = structure_fingerprint(st)

    call test_parse_index_list()
    call test_roundtrip()
    call test_overlap()
    call test_add_region_errors()
    call test_validate_errors()
    call test_file_errors()

    call cleanup()
    call test_summary('test_regions')

contains

    ! ------------------------------------------------------------------
    subroutine make_structure(s)
        type(structure_t), intent(out) :: s
        integer :: i
        s%nat = 8
        s%lattice = 0.0_dp
        s%lattice(1,1) = 3.2_dp; s%lattice(2,2) = 3.2_dp; s%lattice(3,3) = 20.0_dp
        allocate(s%cart(3, 8), s%label(8))
        do i = 1, 8
            s%cart(:, i) = [0.1_dp * i, 0.2_dp * i, 1.5_dp * i]
            if (mod(i, 2) == 1) then
                s%label(i) = 'Ga'
            else
                s%label(i) = 'N'
            end if
        end do
    end subroutine make_structure

    subroutine base_set(rs, with_species)
        type(region_set_t), intent(out) :: rs
        logical,            intent(in)  :: with_species
        rs%format_version = 1
        rs%natoms = 8
        rs%fingerprint = fp
        if (with_species) then
            rs%has_species = .true.
            rs%species = st%label
        end if
    end subroutine base_set

    function tmpname(tag) result(f)
        character(len=*), intent(in) :: tag
        character(len=256) :: f
        f = TMPDIR // 'test_regions_' // tag // '.nml'
        ntmp = ntmp + 1
        tmpfiles(ntmp) = f
    end function tmpname

    subroutine write_text(fname, lines)
        character(len=*), intent(in) :: fname
        character(len=*), intent(in) :: lines(:)
        integer :: iu, i
        open(newunit=iu, file=fname, status='replace', action='write')
        do i = 1, size(lines)
            write(iu, '(A)') trim(lines(i))
        end do
        close(iu)
    end subroutine write_text

    subroutine cleanup()
        integer :: i, iu, ios
        logical :: ex
        do i = 1, ntmp
            inquire(file=tmpfiles(i), exist=ex)
            if (ex) then
                open(newunit=iu, file=tmpfiles(i), status='old', iostat=ios)
                if (ios == 0) close(iu, status='delete')
            end if
        end do
    end subroutine cleanup

    ! ------------------------------------------------------------------
    subroutine test_parse_index_list()
        integer, allocatable :: l(:)
        integer :: ierr, i
        character(len=256) :: msg

        call parse_index_list('1-12', 12, l, ierr, msg)
        call check(ierr == 0 .and. size(l) == 12, 'parse "1-12"')
        if (ierr == 0) call check(all(l == [(i, i = 1, 12)]), 'parse "1-12" values')

        call parse_index_list('15, 17, 22', 30, l, ierr, msg)
        call check(ierr == 0, 'parse "15, 17, 22"')
        if (ierr == 0) call check(size(l) == 3 .and. all(l == [15, 17, 22]), 'parse "15, 17, 22" values')

        call parse_index_list('1-10, 15-20', 20, l, ierr, msg)
        call check(ierr == 0 .and. size(l) == 16, 'parse "1-10, 15-20"')
        if (ierr == 0) call check(all(l == [(i, i = 1, 10), (i, i = 15, 20)]), 'parse "1-10, 15-20" values')

        call parse_index_list(' 3 ,1-2', 12, l, ierr, msg)
        call check(ierr == 0, 'parse " 3 ,1-2"')
        if (ierr == 0) call check(size(l) == 3 .and. all(l == [1, 2, 3]), 'parse " 3 ,1-2" sorted')

        call parse_index_list(' 4 - 6 ', 12, l, ierr, msg)
        call check(ierr == 0 .and. size(l) == 3, 'parse " 4 - 6 " with spaces')

        call parse_index_list('', 12, l, ierr, msg)
        call check(ierr /= 0, 'parse "" -> error')
        call parse_index_list('1,,2', 12, l, ierr, msg)
        call check(ierr /= 0, 'parse "1,,2" -> error')
        call parse_index_list('1,', 12, l, ierr, msg)
        call check(ierr /= 0, 'parse "1," -> error')
        call parse_index_list('a', 12, l, ierr, msg)
        call check(ierr /= 0, 'parse "a" -> error')
        call parse_index_list('5-3', 12, l, ierr, msg)
        call check(ierr /= 0, 'parse "5-3" -> error')
        call parse_index_list('0', 12, l, ierr, msg)
        call check(ierr /= 0, 'parse "0" -> error')
        call parse_index_list('13', 12, l, ierr, msg)
        call check(ierr /= 0, 'parse "13" (nmax 12) -> error')
        call parse_index_list('1-5,3', 12, l, ierr, msg)
        call check(ierr /= 0, 'parse "1-5,3" duplicates -> error')
        call check(index(msg, '3') > 0 .and. index(msg, 'duplicat') > 0, 'duplicate message names index 3')
    end subroutine test_parse_index_list

    ! ------------------------------------------------------------------
    subroutine test_roundtrip()
        type(region_set_t) :: rs, rs2
        integer :: ierr, ir
        character(len=256) :: msg, f
        logical :: same

        call base_set(rs, .true.)
        call add_region(rs, 'top', [1, 2, 3, 4], ierr, msg)
        call check(ierr == 0, 'add_region top')
        call add_region(rs, 'bottom', [5, 6, 7, 8], ierr, msg)
        call check(ierr == 0, 'add_region bottom')
        call check(rs%nregion == 2, 'nregion == 2')

        f = tmpname('roundtrip')
        call write_regions_file(trim(f), rs, ierr, msg)
        call check(ierr == 0, 'write_regions_file')
        call read_regions_file(trim(f), rs2, ierr, msg)
        call check(ierr == 0, 'read_regions_file roundtrip: ' // trim(msg))
        if (ierr /= 0) return

        same = rs2%format_version == 1 .and. rs2%natoms == 8 .and. rs2%fingerprint == fp &
               .and. rs2%has_species .and. rs2%nregion == 2
        call check(same, 'roundtrip header fields')
        if (rs2%has_species) call check(all(rs2%species == rs%species), 'roundtrip species')
        if (rs2%nregion == 2) then
            do ir = 1, 2
                call check(rs2%regions(ir)%name == rs%regions(ir)%name, 'roundtrip region name')
                call check(size(rs2%regions(ir)%atoms) == size(rs%regions(ir)%atoms), 'roundtrip atoms size')
                if (size(rs2%regions(ir)%atoms) == size(rs%regions(ir)%atoms)) &
                    call check(all(rs2%regions(ir)%atoms == rs%regions(ir)%atoms), 'roundtrip atoms')
            end do
        end if

        call validate_regions(rs2, st, .false., ierr, msg)
        call check(ierr == 0, 'validate roundtrip set: ' // trim(msg))
        call check(regions_are_partition(rs2), 'top+bottom is a partition')

        call test_long_roundtrip()
    end subroutine test_roundtrip

    !> Larger set to exercise line wrapping of species and atoms.
    subroutine test_long_roundtrip()
        type(region_set_t) :: rs, rs2
        integer :: ierr, i
        character(len=256) :: msg, f
        rs%natoms = 45
        rs%has_species = .true.
        allocate(rs%species(45))
        do i = 1, 45
            write(rs%species(i), '(A,I0)') 'X', mod(i, 7)
        end do
        call add_region(rs, 'odd.part_1', [(i, i = 45, 1, -2)], ierr, msg)
        call check(ierr == 0, 'long add_region')
        f = tmpname('long')
        call write_regions_file(trim(f), rs, ierr, msg)
        call read_regions_file(trim(f), rs2, ierr, msg)
        call check(ierr == 0, 'long roundtrip read: ' // trim(msg))
        if (ierr /= 0) return
        call check(rs2%fingerprint == '' .and. all(rs2%species == rs%species), 'long roundtrip species')
        call check(all(rs2%regions(1)%atoms == rs%regions(1)%atoms), 'long roundtrip atoms (user order kept)')
        call check(.not. regions_are_partition(rs2), 'odd only is not a partition')
    end subroutine test_long_roundtrip

    ! ------------------------------------------------------------------
    subroutine test_overlap()
        type(region_set_t) :: rs, rs2
        integer :: ierr
        character(len=256) :: msg, f
        call base_set(rs, .false.)
        call add_region(rs, 'top', [1, 2, 3, 4], ierr, msg)
        call add_region(rs, 'bottom', [5, 6, 7, 8], ierr, msg)
        call add_region(rs, 'mid', [3, 4, 5, 6], ierr, msg)
        call check(ierr == 0 .and. rs%nregion == 3, 'add overlapping region mid')
        f = tmpname('overlap')
        call write_regions_file(trim(f), rs, ierr, msg)
        call read_regions_file(trim(f), rs2, ierr, msg)
        call check(ierr == 0 .and. rs2%nregion == 3, 'read overlapping set')
        call validate_regions(rs2, st, .false., ierr, msg)
        call check(ierr == 0, 'validate overlapping set')
        call check(.not. regions_are_partition(rs2), 'overlapping set is not a partition')
    end subroutine test_overlap

    ! ------------------------------------------------------------------
    subroutine test_add_region_errors()
        type(region_set_t) :: rs
        integer :: ierr
        character(len=256) :: msg
        integer :: empty(0)

        call base_set(rs, .false.)
        call add_region(rs, 'top', [1, 2, 3, 4], ierr, msg)

        call add_region(rs, 'a', [0, 1], ierr, msg)
        call check(ierr /= 0, 'add_region index 0 -> error')
        call add_region(rs, 'b', [1, 9], ierr, msg)
        call check(ierr /= 0, 'add_region index 9 -> error')
        call add_region(rs, 'c', [1, 2, 1], ierr, msg)
        call check(ierr /= 0 .and. index(msg, '1') > 0, 'add_region duplicate index -> error')
        call add_region(rs, 'd', empty, ierr, msg)
        call check(ierr /= 0, 'add_region empty -> error')
        call add_region(rs, 'top', [5], ierr, msg)
        call check(ierr /= 0, 'add_region duplicated name -> error')
        call add_region(rs, 'total', [5], ierr, msg)
        call check(ierr /= 0, "add_region name 'total' -> error")
        call add_region(rs, 'a-b', [5], ierr, msg)
        call check(ierr /= 0, "add_region name with '-' -> error")
        call add_region(rs, 'a b', [5], ierr, msg)
        call check(ierr /= 0, 'add_region name with space -> error')
        call add_region(rs, repeat('x', 65), [5], ierr, msg)
        call check(ierr /= 0, 'add_region name too long -> error')
        call check(rs%nregion == 1, 'failed add_region leaves set unchanged')

        call check(is_valid_region_name('Top_1.a'), 'valid name Top_1.a')
        call check(.not. is_valid_region_name(''), 'empty name invalid')
    end subroutine test_add_region_errors

    ! ------------------------------------------------------------------
    subroutine test_validate_errors()
        type(region_set_t) :: rs
        type(structure_t) :: st2
        integer :: ierr
        character(len=256) :: msg

        ! natoms mismatch
        call base_set(rs, .false.)
        rs%natoms = 7
        call add_region(rs, 'top', [1, 2], ierr, msg)
        call validate_regions(rs, st, .true., ierr, msg)
        call check(ierr /= 0, 'validate natoms mismatch -> error')

        ! wrong fingerprint
        call base_set(rs, .false.)
        rs%fingerprint = '0123456789abcdef'
        call add_region(rs, 'top', [1, 2], ierr, msg)
        call validate_regions(rs, st, .true., ierr, msg)
        call check(ierr /= 0, 'validate wrong fingerprint -> error')

        ! modified structure -> fingerprint differs
        st2 = st
        st2%cart(3, 8) = st2%cart(3, 8) + 0.01_dp
        call base_set(rs, .false.)
        call add_region(rs, 'top', [1, 2], ierr, msg)
        call validate_regions(rs, st2, .true., ierr, msg)
        call check(ierr /= 0, 'validate vs moved structure -> error')

        ! missing fingerprint
        call base_set(rs, .false.)
        rs%fingerprint = ''
        call add_region(rs, 'top', [1, 2], ierr, msg)
        call validate_regions(rs, st, .false., ierr, msg)
        call check(ierr /= 0, 'validate missing fingerprint (no allow) -> error')
        call validate_regions(rs, st, .true., ierr, msg)
        call check(ierr == 0, 'validate missing fingerprint (allow) -> ok')

        ! species mismatch
        call base_set(rs, .true.)
        rs%species(3) = 'N'
        call add_region(rs, 'top', [1, 2], ierr, msg)
        call validate_regions(rs, st, .false., ierr, msg)
        call check(ierr /= 0, 'validate species mismatch -> error')
    end subroutine test_validate_errors

    ! ------------------------------------------------------------------
    subroutine test_file_errors()
        type(region_set_t) :: rs
        integer :: ierr
        character(len=256) :: msg, f
        character(len=80) :: hdr(4)

        hdr(1) = '&REGION_FILE'
        hdr(2) = '    format_version = 1'
        hdr(3) = '    natoms = 8'
        hdr(4) = '/'

        f = tmpname('gap')
        call write_text(trim(f), [character(len=80) :: hdr, &
            '&REGION', "    name = 'top'", '    atoms(1) = 1, atoms(3) = 3', '/'])
        call read_regions_file(trim(f), rs, ierr, msg)
        call check(ierr /= 0 .and. index(msg, 'gap') > 0, 'gap in atoms list -> error')

        f = tmpname('syntax')
        call write_text(trim(f), [character(len=80) :: hdr, &
            '&REGION', "    name = 'top'", '    atoms = 1, 2', '    bogus = 3', '/'])
        call read_regions_file(trim(f), rs, ierr, msg)
        call check(ierr /= 0 .and. index(msg, 'block 1') > 0, 'unknown variable in &REGION -> error')

        f = tmpname('noname')
        call write_text(trim(f), [character(len=80) :: hdr, &
            '&REGION', "    name = 'top'", '    atoms = 1, 2', '/', &
            '&REGION', '    atoms = 3, 4', '/'])
        call read_regions_file(trim(f), rs, ierr, msg)
        call check(ierr /= 0 .and. index(msg, 'block 2') > 0, 'block without name -> error')

        f = tmpname('idx9')
        call write_text(trim(f), [character(len=80) :: hdr, &
            '&REGION', "    name = 'top'", '    atoms = 1, 9', '/'])
        call read_regions_file(trim(f), rs, ierr, msg)
        call check(ierr /= 0, 'file index 9 -> error')

        f = tmpname('dupname')
        call write_text(trim(f), [character(len=80) :: hdr, &
            '&REGION', "    name = 'top'", '    atoms = 1', '/', &
            '&REGION', "    name = 'top'", '    atoms = 2', '/'])
        call read_regions_file(trim(f), rs, ierr, msg)
        call check(ierr /= 0, 'file duplicated name -> error')

        f = tmpname('version')
        call write_text(trim(f), [character(len=80) :: '&REGION_FILE', '    format_version = 2', &
            '    natoms = 8', '/'])
        call read_regions_file(trim(f), rs, ierr, msg)
        call check(ierr /= 0, 'format_version 2 -> error')

        f = tmpname('nspecies')
        call write_text(trim(f), [character(len=80) :: '&REGION_FILE', '    format_version = 1', &
            '    natoms = 8', "    species = 'Ga', 'N'", '/'])
        call read_regions_file(trim(f), rs, ierr, msg)
        call check(ierr /= 0, 'species count /= natoms -> error')

        f = tmpname('order')
        call write_text(trim(f), [character(len=80) :: &
            '&REGION', "    name = 'top'", '    atoms = 1', '/', hdr])
        call read_regions_file(trim(f), rs, ierr, msg)
        call check(ierr /= 0, '&REGION_FILE not first -> error')

        f = tmpname('noregions')
        call write_text(trim(f), hdr)
        call read_regions_file(trim(f), rs, ierr, msg)
        call check(ierr == 0 .and. rs%nregion == 0, 'header only -> ok, zero regions')

        call read_regions_file(TMPDIR // 'test_regions_does_not_exist.nml', rs, ierr, msg)
        call check(ierr /= 0, 'missing file -> error')
    end subroutine test_file_errors

end program test_regions
