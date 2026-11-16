! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Map expressions over named regions and text output of spectral maps
!> (docs/design/contracts.md §7).
!>
!> Pipeline position: used by the specmap tool after the broadening step
!> (spectral_broadening) has turned discrete (energy, weight) pairs into an
!> intensity map A(ne,nk); this module selects/combines regions and writes
!> the gnuplot-ready outputs. It does no physics itself.
!>
!>   - parse_map_expression: 'top-bottom', 'a+b-c', '-top' -> region indices
!>     and +/-1 coefficients (a repeated region has its coefficients added)
!>   - map_file_tag:         'top-bottom' -> 'top_minus_bottom'
!>   - write_map_file:       pm3d-ready text map (one block per k)
!>   - write_kticks_file:    gnuplot include with xtics/xrange/EF/ELABEL
!>
!> Conventions: energies in eV, k distance in 1/Angstrom (as in
!> spectral_data); region indices are 1-based into sd%region_name.
!> Library routines never stop: errors are returned through ierr/errmsg
!> (ierr = 0 means success).
module spectral_map
    use types, only: dp
    use spectral_data, only: spectral_data_t, find_region, EREF_FERMI, EREF_ABSOLUTE
    implicit none
    private

    public :: parse_map_expression, write_map_file, write_kticks_file, map_file_tag

contains

    ! ==================================================================
    ! Map expressions
    ! ==================================================================

    !> Parse a map expression: region names joined by '+' or '-'; a leading
    !> '-' is allowed; spaces are ignored. Returns the distinct region indices
    !> (in order of first appearance) and their coefficients (sum of the +/-1
    !> of every occurrence; a repeated region is combined, e.g. 'top+top' -> 2).
    !> Errors: empty expression, empty term ('top+', '++', '+top'), unknown
    !> region (the message lists the available names).
    !> Note: a region whose name itself contains '+' or '-' cannot be addressed.
    !> On error (ierr = 1) region_idx and coeff are returned with size 0.
    !> The caller combines per-region intensities as sum(coeff(i)*I(region_idx(i))).
    subroutine parse_map_expression(expr, sd, region_idx, coeff, ierr, errmsg)
        character(len=*),      intent(in)  :: expr
        type(spectral_data_t), intent(in)  :: sd
        integer,  allocatable, intent(out) :: region_idx(:)
        real(dp), allocatable, intent(out) :: coeff(:)
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        character(len=:), allocatable :: s, name
        integer  :: i, n, nterm, idx, j
        real(dp) :: sgn
        integer,  allocatable :: tidx(:)
        real(dp), allocatable :: tcoef(:)
        logical :: found, at_end

        ierr = 0
        errmsg = ''
        allocate(region_idx(0), coeff(0))

        s = strip_spaces(expr)
        n = len(s)
        if (n == 0) then
            call fail("empty map expression")
            return
        end if

        ! Single left-to-right scan. 'name' accumulates the current term;
        ! '+'/'-' (or the end of the string, i = n+1) closes it, and the
        ! operator that closed it sets the sign of the NEXT term. At most n
        ! terms can exist, hence the size of the scratch arrays.
        allocate(tidx(n), tcoef(n))
        nterm = 0
        sgn = 1.0_dp
        name = ''
        i = 1
        ! A leading '-' negates the first term ('+' at the start is an
        ! empty-term error, caught below).
        if (s(1:1) == '-') then
            sgn = -1.0_dp
            i = 2
        end if

        do while (i <= n + 1)
            if (i > n) then
                at_end = .true.
            else
                at_end = (s(i:i) == '+' .or. s(i:i) == '-')
            end if
            if (at_end) then
                ! End of a term.
                if (len(name) == 0) then
                    call fail("empty term in map expression '" // trim(expr) // "'")
                    return
                end if
                idx = find_region(sd, name)
                if (idx == 0) then
                    call fail("unknown region '" // name // "' in map expression '" // &
                              trim(expr) // "'; available: " // available_names(sd))
                    return
                end if
                ! Merge repeated regions by adding their signed coefficients
                ! (so 'top-top' yields coefficient 0, not an error).
                found = .false.
                do j = 1, nterm
                    if (tidx(j) == idx) then
                        tcoef(j) = tcoef(j) + sgn
                        found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    nterm = nterm + 1
                    tidx(nterm) = idx
                    tcoef(nterm) = sgn
                end if
                ! Sign for the following term (none left at end of string).
                if (i <= n) then
                    if (s(i:i) == '+') then
                        sgn = 1.0_dp
                    else
                        sgn = -1.0_dp
                    end if
                end if
                name = ''
            else
                name = name // s(i:i)
            end if
            i = i + 1
        end do

        region_idx = tidx(1:nterm)
        coeff = tcoef(1:nterm)

    contains

        subroutine fail(msg)
            character(len=*), intent(in) :: msg
            ierr = 1
            errmsg = msg
            if (allocated(region_idx)) deallocate(region_idx)
            if (allocated(coeff)) deallocate(coeff)
            allocate(region_idx(0), coeff(0))
        end subroutine fail

    end subroutine parse_map_expression

    !> File-name tag of a map expression: spaces removed, a leading '-'
    !> becomes 'minus_', then '+' -> '_plus_' and '-' -> '_minus_'.
    !> Used to build '<prefix>.<tag>.map.dat'. The expression is not
    !> validated here.
    function map_file_tag(expr) result(tag)
        character(len=*), intent(in) :: expr
        character(len=:), allocatable :: tag
        character(len=:), allocatable :: s
        integer :: i, i0

        s = strip_spaces(expr)
        tag = ''
        i0 = 1
        if (len(s) > 0) then
            if (s(1:1) == '-') then
                tag = 'minus_'
                i0 = 2
            end if
        end if
        do i = i0, len(s)
            select case (s(i:i))
            case ('+'); tag = tag // '_plus_'
            case ('-'); tag = tag // '_minus_'
            case default; tag = tag // s(i:i)
            end select
        end do
    end function map_file_tag

    ! ==================================================================
    ! Output files
    ! ==================================================================

    !> Write a map A(ne,nk) on energy grid(ne) along sd%kdist(nk).
    !> grid in eV, kdist in 1/Angstrom; A is the intensity (arbitrary units).
    !> Each header line is written prefixed with '# ' (unless it already
    !> starts with '#'). Then, for every k, ne rows
    !> 'k_distance energy_eV intensity' followed by one blank line, so the
    !> file contains nk*ne data lines and exactly nk blank lines (the last
    !> one terminates the final block; harmless for gnuplot pm3d).
    !> ierr: 0 ok; 1 size mismatch or no kdist; 2 cannot open; 3 write/close error.
    subroutine write_map_file(filename, sd, A, grid, header_lines, ierr, errmsg)
        character(len=*),      intent(in)  :: filename
        type(spectral_data_t), intent(in)  :: sd
        real(dp),              intent(in)  :: A(:,:), grid(:)
        character(len=*),      intent(in)  :: header_lines(:)
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        integer :: u, ios, ios2, ik, ie, i
        character(len=*), parameter :: ROWFMT = '(ES17.9E3, 1X, ES17.9E3, 1X, ES18.10E3)'

        ierr = 0
        errmsg = ''
        if (size(A, 1) /= size(grid)) then
            ierr = 1
            errmsg = "write_map_file: size(A,1) /= size(grid)"
            return
        end if
        if (.not. allocated(sd%kdist)) then
            ierr = 1
            errmsg = "write_map_file: spectral data has no kdist"
            return
        end if
        if (size(A, 2) /= size(sd%kdist)) then
            ierr = 1
            errmsg = "write_map_file: size(A,2) /= nk"
            return
        end if

        open(newunit=u, file=filename, status='replace', action='write', &
             form='formatted', iostat=ios)
        if (ios /= 0) then
            ierr = 2
            errmsg = "cannot open map file for writing: " // trim(filename)
            return
        end if

        ! Header: every line is made a gnuplot comment.
        do i = 1, size(header_lines)
            if (len_trim(header_lines(i)) > 0) then
                if (header_lines(i)(1:1) == '#') then
                    write(u, '(A)', iostat=ios) trim(header_lines(i))
                else
                    write(u, '(A)', iostat=ios) '# ' // trim(header_lines(i))
                end if
            else
                write(u, '(A)', iostat=ios) '#'
            end if
            if (ios /= 0) exit
        end do
        if (ios == 0) write(u, '(A)', iostat=ios) '# columns: k_distance energy_eV intensity'

        ! Data: k-major blocks (energy runs fastest) separated by blank
        ! lines, the layout gnuplot's 'splot ... with pm3d' expects.
        do ik = 1, size(A, 2)
            if (ios /= 0) exit
            do ie = 1, size(grid)
                write(u, ROWFMT, iostat=ios) sd%kdist(ik), grid(ie), A(ie, ik)
                if (ios /= 0) exit
            end do
            if (ios == 0) write(u, '(A)', iostat=ios) ''
        end do

        if (ios /= 0) then
            ierr = 3
            errmsg = "error writing map file: " // trim(filename)
            close(u, iostat=ios2)
            return
        end if
        close(u, iostat=ios)
        if (ios /= 0) then
            ierr = 3
            errmsg = "error closing map file: " // trim(filename)
        end if
    end subroutine write_map_file

    !> Write a gnuplot include file with:
    !>   set xrange [0:kmax]
    !>   set xtics ("label" kdist, ...)     ('G'/'Gamma' -> {/Symbol G})
    !>   vertical dashed lines at interior vertices (colour VLINE_COLOR,
    !>   default "white"; a script may define VLINE_COLOR before loading)
    !>   EF = 0 (fermi reference) | E_F (absolute & known) | NaN
    !>   ELABEL = "E - E_F (eV)" | "E (eV)" | "E (eV, reference unknown)"
    !> The EF/ELABEL choice follows sd%energy_ref, so plots can mark the
    !> Fermi level correctly whichever zero the energies use. Vertex ticks
    !> are skipped if sd has no vertices.
    !> ierr: 0 ok; 1 no k-points; 2 cannot open; 3 error on close.
    subroutine write_kticks_file(filename, sd, ierr, errmsg)
        character(len=*),      intent(in)  :: filename
        type(spectral_data_t), intent(in)  :: sd
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        integer  :: u, ios, ios2, iv, ik
        logical  :: first
        real(dp) :: kmax, kd
        character(len=:), allocatable :: line
        character(len=32) :: num

        ierr = 0
        errmsg = ''
        if (.not. allocated(sd%kdist) .or. sd%nk < 1) then
            ierr = 1
            errmsg = "write_kticks_file: spectral data has no k-points"
            return
        end if
        kmax = sd%kdist(sd%nk)

        open(newunit=u, file=filename, status='replace', action='write', &
             form='formatted', iostat=ios)
        if (ios /= 0) then
            ierr = 2
            errmsg = "cannot open kticks file for writing: " // trim(filename)
            return
        end if

        if (ios == 0) write(u, '(A)', iostat=ios) '# gnuplot include written by specmap: k-path ticks, Fermi level, energy label'
        if (ios == 0) write(u, '(A)', iostat=ios) 'set xrange [0:' // fmt_real(kmax) // ']'

        if (sd%nvertex > 0 .and. allocated(sd%vertex_k) .and. allocated(sd%vertex_label)) then
            line = 'set xtics ('
            first = .true.
            do iv = 1, sd%nvertex
                ik = sd%vertex_k(iv)
                if (ik < 1 .or. ik > sd%nk) cycle
                if (.not. first) line = line // ', '
                first = .false.
                line = line // '"' // gp_label(sd%vertex_label(iv)) // '" ' // fmt_real(sd%kdist(ik))
            end do
            line = line // ')'
            if (ios == 0) write(u, '(A)', iostat=ios) line
            if (ios == 0) write(u, '(A)', iostat=ios) 'if (!exists("VLINE_COLOR")) VLINE_COLOR = "white"'
            do iv = 1, sd%nvertex
                ik = sd%vertex_k(iv)
                if (ik < 1 .or. ik > sd%nk) cycle
                kd = sd%kdist(ik)
                if (kd <= 0.0_dp .or. kd >= kmax) cycle
                num = fmt_real(kd)
                if (ios == 0) write(u, '(A)', iostat=ios) 'set arrow from ' // trim(num) // ', graph 0 to ' // trim(num) // &
                    ', graph 1 nohead dt 2 lw 1 lc rgb VLINE_COLOR front'
            end do
        end if

        ! Energy reference: with EREF_FERMI the energies are already
        ! relative to E_F, so the Fermi level sits at 0.
        if (sd%energy_ref == EREF_FERMI) then
            if (ios == 0) write(u, '(A)', iostat=ios) 'EF = 0.0'
            if (ios == 0) write(u, '(A)', iostat=ios) 'ELABEL = "E - E_F (eV)"'
        else if (sd%energy_ref == EREF_ABSOLUTE) then
            if (sd%fermi_known) then
                if (ios == 0) write(u, '(A)', iostat=ios) 'EF = ' // fmt_real(sd%fermi_energy_ev)
            else
                if (ios == 0) write(u, '(A)', iostat=ios) 'EF = NaN'
            end if
            if (ios == 0) write(u, '(A)', iostat=ios) 'ELABEL = "E (eV)"'
        else
            if (ios == 0) write(u, '(A)', iostat=ios) 'EF = NaN'
            if (ios == 0) write(u, '(A)', iostat=ios) 'ELABEL = "E (eV, reference unknown)"'
        end if

        if (ios /= 0) then
            ierr = 3
            errmsg = "error writing kticks file: " // trim(filename)
            close(u, iostat=ios2)
            return
        end if

        close(u, iostat=ios)
        if (ios /= 0) then
            ierr = 3
            errmsg = "error closing kticks file: " // trim(filename)
        end if
    end subroutine write_kticks_file

    ! ==================================================================
    ! Private helpers
    ! ==================================================================

    !> Remove blanks and tabs.
    pure function strip_spaces(s) result(r)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: r
        integer :: i
        r = ''
        do i = 1, len(s)
            if (s(i:i) /= ' ' .and. s(i:i) /= achar(9)) r = r // s(i:i)
        end do
    end function strip_spaces

    !> Comma-separated region names, for "unknown region" error messages.
    function available_names(sd) result(r)
        type(spectral_data_t), intent(in) :: sd
        character(len=:), allocatable :: r
        integer :: i
        r = ''
        if (.not. allocated(sd%region_name)) return
        do i = 1, size(sd%region_name)
            if (i > 1) r = r // ', '
            r = r // trim(sd%region_name(i))
        end do
    end function available_names

    !> Gnuplot (enhanced text) label: 'G'/'Gamma' -> '{/Symbol G}';
    !> double quotes are removed.
    function gp_label(lab) result(r)
        character(len=*), intent(in) :: lab
        character(len=:), allocatable :: r
        character(len=:), allocatable :: t
        integer :: i
        t = trim(adjustl(lab))
        if (t == 'G' .or. t == 'g' .or. t == 'Gamma' .or. t == 'GAMMA' .or. t == 'gamma') then
            r = '{/Symbol G}'
            return
        end if
        r = ''
        do i = 1, len(t)
            if (t(i:i) /= '"') r = r // t(i:i)
        end do
    end function gp_label

    !> Compact scientific notation without padding, for gnuplot numbers.
    function fmt_real(x) result(r)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: r
        character(len=32) :: buf
        write(buf, '(ES16.8E3)') x
        r = trim(adjustl(buf))
    end function fmt_real

end module spectral_map
