! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! tools/atomselect.f90
!
! Command-line tool for named atomic regions (docs/design/contracts.md §6):
! inspects a structure file, prints its fingerprint, creates/extends a
! regions.nml file with numeric atom selections and validates a regions file
! against a structure.
!
! Usage:
!   atomselect inspect     <structure.nml> [--xyz out.xyz] [--table out.txt]
!   atomselect fingerprint <structure.nml>
!   atomselect select      <structure.nml> --regions <regions.nml> --name NAME --atoms "LIST"
!   atomselect validate    <structure.nml> <regions.nml> [--allow-no-fingerprint]
!   atomselect help
!
! Exit code 0 on success, 1 on any error (message prefixed "[atomselect error]").
! ----------------------------------------------------------------------------------------------------------
program atomselect
    use structure_io,     only: structure_t, read_structure, structure_fingerprint
    use atomic_regions,   only: region_set_t, read_regions_file, write_regions_file, &
                                validate_regions, add_region, regions_are_partition
    use atomic_selection, only: parse_index_list
    implicit none

    type :: str_t
        character(len=:), allocatable :: s
    end type str_t

    integer, parameter :: MSG_LEN = 1024
    integer, parameter :: OPT_NAME_LEN = 32

    type(str_t), allocatable :: args(:)       ! all command-line arguments
    type(str_t), allocatable :: opt_key(:)    ! parsed options (after positionals)
    type(str_t), allocatable :: opt_val(:)
    integer :: nargs, nopt
    character(len=:), allocatable :: cmd

    call collect_arguments()

    if (nargs < 1) then
        call print_usage()
        stop
    end if

    cmd = args(1)%s
    select case (cmd)
    case ('help', '-h', '--help')
        call print_usage()
    case ('inspect')
        call cmd_inspect()
    case ('fingerprint')
        call cmd_fingerprint()
    case ('select')
        call cmd_select()
    case ('validate')
        call cmd_validate()
    case default
        call die("unknown subcommand '" // cmd // "' (run 'atomselect help')")
    end select

contains

    ! ==================================================================
    ! Subcommands
    ! ==================================================================

    subroutine cmd_inspect()
        type(structure_t) :: st
        character(len=16) :: fp
        integer :: iu, ios

        call require_positionals(1, 'inspect <structure.nml>')
        call parse_options(3, [character(len=OPT_NAME_LEN) :: '--xyz', '--table'], &
                           [character(len=OPT_NAME_LEN) :: ''])
        call load_structure(args(2)%s, st)
        fp = structure_fingerprint(st)

        call write_header(6, args(2)%s, st, fp)
        call write_table(6, st)

        if (has_opt('--table')) then
            open(newunit=iu, file=opt_value('--table'), status='replace', action='write', iostat=ios)
            if (ios /= 0) call die("cannot open table output file '" // opt_value('--table') // "'")
            call write_header(iu, args(2)%s, st, fp)
            call write_table(iu, st)
            close(iu)
            print '(A)', 'table written to ' // opt_value('--table')
        end if

        if (has_opt('--xyz')) then
            open(newunit=iu, file=opt_value('--xyz'), status='replace', action='write', iostat=ios)
            if (ios /= 0) call die("cannot open xyz output file '" // opt_value('--xyz') // "'")
            call write_xyz(iu, st, fp)
            close(iu)
            print '(A)', 'xyz written to ' // opt_value('--xyz')
        end if
    end subroutine cmd_inspect

    ! ------------------------------------------------------------------
    subroutine cmd_fingerprint()
        type(structure_t) :: st

        call require_positionals(1, 'fingerprint <structure.nml>')
        call parse_options(3, [character(len=OPT_NAME_LEN) :: ''], [character(len=OPT_NAME_LEN) :: ''])
        call load_structure(args(2)%s, st)
        print '(A)', structure_fingerprint(st)
    end subroutine cmd_fingerprint

    ! ------------------------------------------------------------------
    subroutine cmd_select()
        type(structure_t)  :: st
        type(region_set_t) :: rs
        integer, allocatable :: list(:)
        character(len=:), allocatable :: rfile, rname
        character(len=MSG_LEN) :: errmsg
        integer :: ierr
        logical :: exists

        call require_positionals(1, 'select <structure.nml> --regions <regions.nml> --name NAME --atoms "LIST"')
        call parse_options(3, [character(len=OPT_NAME_LEN) :: '--regions', '--name', '--atoms'], &
                           [character(len=OPT_NAME_LEN) :: ''])
        if (.not. has_opt('--regions')) call die('select: missing required option --regions')
        if (.not. has_opt('--name'))    call die('select: missing required option --name')
        if (.not. has_opt('--atoms'))   call die('select: missing required option --atoms')
        rfile = opt_value('--regions')
        rname = opt_value('--name')

        call load_structure(args(2)%s, st)

        call parse_index_list(opt_value('--atoms'), st%nat, list, ierr, errmsg)
        if (ierr /= 0) call die('select: --atoms: ' // trim(errmsg))

        inquire(file=rfile, exist=exists)
        if (exists) then
            call read_regions_file(rfile, rs, ierr, errmsg)
            if (ierr /= 0) call die(trim(errmsg))
            call validate_regions(rs, st, .false., ierr, errmsg)
            if (ierr /= 0) call die(trim(errmsg))
        else
            rs%format_version = 1
            rs%natoms = st%nat
            rs%fingerprint = structure_fingerprint(st)
            rs%has_species = .true.
            rs%species = st%label
            rs%nregion = 0
            allocate(rs%regions(0))
        end if

        call add_region(rs, rname, list, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))

        call write_regions_file(rfile, rs, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))

        if (exists) then
            print '(A)', "region '" // rname // "': " // itoa(size(list)) // ' atoms [' // &
                         compress_list(list) // '] added to ' // rfile // ' (updated)'
        else
            print '(A)', "region '" // rname // "': " // itoa(size(list)) // ' atoms [' // &
                         compress_list(list) // '] added to ' // rfile // ' (created)'
        end if
    end subroutine cmd_select

    ! ------------------------------------------------------------------
    subroutine cmd_validate()
        type(structure_t)  :: st
        type(region_set_t) :: rs
        character(len=MSG_LEN) :: errmsg
        integer, allocatable :: cnt(:)
        logical, allocatable :: in_i(:)
        integer :: ierr, ir, jr, a, nshared, noverlap, nuncovered

        call require_positionals(2, 'validate <structure.nml> <regions.nml>')
        call parse_options(4, [character(len=OPT_NAME_LEN) :: ''], &
                           [character(len=OPT_NAME_LEN) :: '--allow-no-fingerprint'])

        call load_structure(args(2)%s, st)
        call read_regions_file(args(3)%s, rs, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))
        call validate_regions(rs, st, has_opt('--allow-no-fingerprint'), ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))

        print '(A)', 'structure   : ' // args(2)%s
        print '(A)', 'regions     : ' // args(3)%s
        print '(A)', 'nat         : ' // itoa(st%nat)
        if (len_trim(rs%fingerprint) > 0) then
            print '(A)', 'fingerprint : ' // trim(rs%fingerprint) // ' (matches)'
        else
            print '(A)', 'fingerprint : none'
        end if
        print '(A)', 'nregion     : ' // itoa(rs%nregion)

        allocate(cnt(st%nat), source=0)
        allocate(in_i(st%nat))
        do ir = 1, rs%nregion
            print '(A)', "region '" // trim(rs%regions(ir)%name) // "': natoms = " // &
                         itoa(size(rs%regions(ir)%atoms)) // ' atoms = [' // &
                         compress_list(sorted(rs%regions(ir)%atoms)) // ']'
            do a = 1, size(rs%regions(ir)%atoms)
                cnt(rs%regions(ir)%atoms(a)) = cnt(rs%regions(ir)%atoms(a)) + 1
            end do
        end do

        noverlap = 0
        do ir = 1, rs%nregion
            in_i = .false.
            in_i(rs%regions(ir)%atoms) = .true.
            do jr = ir + 1, rs%nregion
                nshared = count(in_i(rs%regions(jr)%atoms))
                if (nshared > 0) then
                    noverlap = noverlap + 1
                    print '(A)', "overlap: '" // trim(rs%regions(ir)%name) // "' and '" // &
                                 trim(rs%regions(jr)%name) // "' share " // itoa(nshared) // ' atoms'
                end if
            end do
        end do
        if (noverlap == 0) print '(A)', 'overlaps    : none'

        nuncovered = count(cnt == 0)
        print '(A)', 'uncovered   : ' // itoa(nuncovered) // ' atoms'
        if (regions_are_partition(rs)) then
            print '(A)', 'partition   : yes'
        else
            print '(A)', 'partition   : no'
        end if
        print '(A)', 'valid       : yes'
    end subroutine cmd_validate

    ! ==================================================================
    ! Output helpers
    ! ==================================================================

    subroutine write_header(iu, filename, st, fp)
        integer,           intent(in) :: iu
        character(len=*),  intent(in) :: filename
        type(structure_t), intent(in) :: st
        character(len=16), intent(in) :: fp
        integer :: i

        write(iu, '(A)') '# file        : ' // filename
        write(iu, '(A)') '# nat         : ' // itoa(st%nat)
        do i = 1, 3
            write(iu, '(A,I1,A,3F12.6)') '# a', i, ' (A)      : ', st%lattice(i, :)
        end do
        write(iu, '(A)') '# fingerprint : ' // fp
    end subroutine write_header

    subroutine write_table(iu, st)
        integer,           intent(in) :: iu
        type(structure_t), intent(in) :: st
        integer :: i

        write(iu, '(A6,2X,A8,3A12)') 'Index', 'Element', 'X', 'Y', 'Z'
        do i = 1, st%nat
            write(iu, '(I6,2X,A8,3F12.6)') i, adjustl(st%label(i)), st%cart(:, i)
        end do
    end subroutine write_table

    subroutine write_xyz(iu, st, fp)
        integer,           intent(in) :: iu
        type(structure_t), intent(in) :: st
        character(len=16), intent(in) :: fp
        character(len=:), allocatable :: lat
        character(len=24) :: buf
        integer :: i, j

        lat = ''
        do i = 1, 3
            do j = 1, 3
                write(buf, '(F12.6)') st%lattice(i, j)
                if (len(lat) > 0) lat = lat // ' '
                lat = lat // trim(adjustl(buf))
            end do
        end do
        write(iu, '(I0)') st%nat
        write(iu, '(A)') 'Lattice="' // lat // '" fingerprint=' // fp
        do i = 1, st%nat
            write(iu, '(A,3(1X,F12.6))') trim(adjustl(st%label(i))), st%cart(:, i)
        end do
    end subroutine write_xyz

    subroutine print_usage()
        print '(A)', 'atomselect - inspect structures and manage named atomic regions (regions.nml)'
        print '(A)', ''
        print '(A)', 'Usage:'
        print '(A)', '  atomselect inspect     <structure.nml> [--xyz out.xyz] [--table out.txt]'
        print '(A)', '  atomselect fingerprint <structure.nml>'
        print '(A)', '  atomselect select      <structure.nml> --regions <regions.nml> --name NAME --atoms "LIST"'
        print '(A)', '  atomselect validate    <structure.nml> <regions.nml> [--allow-no-fingerprint]'
        print '(A)', '  atomselect help'
        print '(A)', ''
        print '(A)', 'LIST: comma-separated 1-based indices and inclusive ranges, e.g. "1-4, 7, 10-12".'
    end subroutine print_usage

    ! ==================================================================
    ! Argument handling
    ! ==================================================================

    subroutine collect_arguments()
        integer :: i, n

        nargs = command_argument_count()
        allocate(args(max(nargs, 0)))
        do i = 1, nargs
            call get_command_argument(i, length=n)
            allocate(character(len=n) :: args(i)%s)
            if (n > 0) call get_command_argument(i, args(i)%s)
        end do
        nopt = 0
        allocate(opt_key(0), opt_val(0))
    end subroutine collect_arguments

    !> Requires `npos` positional arguments after the subcommand (args(2:npos+1)),
    !> none of which may look like an option.
    subroutine require_positionals(npos, synopsis)
        integer,          intent(in) :: npos
        character(len=*), intent(in) :: synopsis
        integer :: i

        if (nargs < npos + 1) call die('missing arguments; usage: atomselect ' // synopsis)
        do i = 2, npos + 1
            if (is_option(args(i)%s)) call die('missing arguments; usage: atomselect ' // synopsis)
        end do
    end subroutine require_positionals

    !> Parses args(istart:) as options. `valued` lists options taking one value,
    !> `flags` lists value-less switches (blank entries ignored). Any other
    !> argument, or a repeated option, is an error.
    subroutine parse_options(istart, valued, flags)
        integer,          intent(in) :: istart
        character(len=*), intent(in) :: valued(:), flags(:)
        integer :: i
        character(len=:), allocatable :: a

        i = istart
        do while (i <= nargs)
            a = args(i)%s
            if (.not. is_option(a)) call die("unexpected argument '" // a // "'")
            if (has_opt(a)) call die("option '" // a // "' given more than once")
            if (in_list(a, valued)) then
                if (i + 1 > nargs) call die("option '" // a // "' requires a value")
                call push_opt(a, args(i + 1)%s)
                i = i + 2
            else if (in_list(a, flags)) then
                call push_opt(a, '')
                i = i + 1
            else
                call die("unknown option '" // a // "' for subcommand '" // cmd // "'")
            end if
        end do
    end subroutine parse_options

    subroutine push_opt(key, val)
        character(len=*), intent(in) :: key, val
        type(str_t), allocatable :: tk(:), tv(:)

        allocate(tk(nopt + 1), tv(nopt + 1))
        tk(1:nopt) = opt_key
        tv(1:nopt) = opt_val
        tk(nopt + 1)%s = key
        tv(nopt + 1)%s = val
        call move_alloc(tk, opt_key)
        call move_alloc(tv, opt_val)
        nopt = nopt + 1
    end subroutine push_opt

    logical function has_opt(key)
        character(len=*), intent(in) :: key
        integer :: i
        has_opt = .false.
        do i = 1, nopt
            if (opt_key(i)%s == key) then
                has_opt = .true.
                return
            end if
        end do
    end function has_opt

    function opt_value(key) result(val)
        character(len=*), intent(in) :: key
        character(len=:), allocatable :: val
        integer :: i
        val = ''
        do i = 1, nopt
            if (opt_key(i)%s == key) then
                val = opt_val(i)%s
                return
            end if
        end do
    end function opt_value

    logical function in_list(a, list)
        character(len=*), intent(in) :: a, list(:)
        integer :: i
        in_list = .false.
        do i = 1, size(list)
            if (len_trim(list(i)) == 0) cycle
            if (a == trim(list(i))) then
                in_list = .true.
                return
            end if
        end do
    end function in_list

    logical function is_option(a)
        character(len=*), intent(in) :: a
        is_option = .false.
        if (len(a) >= 2) is_option = (a(1:2) == '--')
    end function is_option

    ! ==================================================================
    ! Misc helpers
    ! ==================================================================

    subroutine load_structure(filename, st)
        character(len=*),  intent(in)  :: filename
        type(structure_t), intent(out) :: st
        character(len=MSG_LEN) :: errmsg
        integer :: ierr

        call read_structure(filename, st, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))
    end subroutine load_structure

    !> Compresses an ascending index list into ranges: [1,2,3,5] -> "1-3,5".
    function compress_list(list) result(s)
        integer, intent(in) :: list(:)
        character(len=:), allocatable :: s
        integer :: i, j

        s = ''
        i = 1
        do while (i <= size(list))
            j = i
            do while (j < size(list))
                if (list(j + 1) /= list(j) + 1) exit
                j = j + 1
            end do
            if (len(s) > 0) s = s // ','
            if (j > i) then
                s = s // itoa(list(i)) // '-' // itoa(list(j))
            else
                s = s // itoa(list(i))
            end if
            i = j + 1
        end do
    end function compress_list

    !> Returns an ascending copy (insertion sort; lists are small and mostly sorted).
    function sorted(list) result(r)
        integer, intent(in) :: list(:)
        integer, allocatable :: r(:)
        integer :: i, j, v

        r = list
        do i = 2, size(r)
            v = r(i)
            j = i - 1
            do while (j >= 1)
                if (r(j) <= v) exit
                r(j + 1) = r(j)
                j = j - 1
            end do
            r(j + 1) = v
        end do
    end function sorted

    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=16) :: buf
        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

    subroutine die(msg)
        character(len=*), intent(in) :: msg
        print '(A)', '[atomselect error] ' // msg
        stop 1
    end subroutine die

end program atomselect
