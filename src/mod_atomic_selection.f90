! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Atom selection helpers. Currently only numeric index lists
!> ("1-12", "15, 17, 22", "1-10, 15-20"). See docs/design/contracts.md §2.
!> Geometric selection and layer routines (§6) are added later (F8).
module atomic_selection
    implicit none
    private

    public :: parse_index_list

    character(len=*), parameter :: DIGITS = '0123456789'
    integer, parameter :: MAX_DIGITS = 9     ! keeps values well inside default integer

contains

    !> Parse a comma-separated list of indices and inclusive ranges a-b.
    !> Result is sorted ascending. Never stops; ierr /= 0 on error.
    subroutine parse_index_list(text, nmax, list, ierr, errmsg)
        character(len=*),     intent(in)  :: text
        integer,              intent(in)  :: nmax
        integer, allocatable, intent(out) :: list(:)
        integer,              intent(out) :: ierr
        character(len=*),     intent(out) :: errmsg

        logical, allocatable :: seen(:)
        integer :: pos, next, ia, ib, k, i, dash
        character(len=:), allocatable :: item

        ierr = 0
        errmsg = ''
        allocate(list(0))

        if (nmax < 1) then
            call fail('index list: nmax must be >= 1')
            return
        end if
        if (len_trim(text) == 0) then
            call fail('index list: empty text')
            return
        end if

        allocate(seen(nmax), source=.false.)
        pos = 1
        k = 0
        do
            next = index(text(pos:), ',')
            if (next == 0) then
                item = trim(adjustl(text(pos:)))
            else
                item = trim(adjustl(text(pos:pos+next-2)))
            end if
            k = k + 1

            if (len(item) == 0) then
                call fail('index list: empty item at position ' // itoa(k))
                return
            end if

            dash = index(item, '-')
            if (dash == 0) then
                call to_int(item, ia)
                if (ierr /= 0) return
                ib = ia
            else
                call to_int(trim(adjustl(item(1:dash-1))), ia)
                if (ierr /= 0) return
                call to_int(trim(adjustl(item(dash+1:))), ib)
                if (ierr /= 0) return
                if (ia > ib) then
                    call fail('index list: reversed range "' // item // '"')
                    return
                end if
            end if

            if (ia < 1 .or. ib > nmax) then
                if (ia < 1) then
                    call fail('index list: index ' // itoa(ia) // ' < 1')
                else
                    call fail('index list: index ' // itoa(ib) // ' > ' // itoa(nmax))
                end if
                return
            end if

            do i = ia, ib
                if (seen(i)) then
                    call fail('index list: duplicated index ' // itoa(i))
                    return
                end if
                seen(i) = .true.
            end do

            if (next == 0) exit
            pos = pos + next
            if (pos > len(text)) then
                call fail('index list: empty item at position ' // itoa(k + 1))
                return
            end if
        end do

        list = pack([(i, i = 1, nmax)], seen)

    contains

        subroutine fail(msg)
            character(len=*), intent(in) :: msg
            ierr = 1
            errmsg = msg
            if (allocated(list)) deallocate(list)
            allocate(list(0))
        end subroutine fail

        subroutine to_int(tok, val)
            character(len=*), intent(in)  :: tok
            integer,          intent(out) :: val
            integer :: ios
            val = 0
            if (len(tok) == 0 .or. len(tok) > MAX_DIGITS .or. verify(tok, DIGITS) /= 0) then
                call fail('index list: not a non-negative integer: "' // tok // '"')
                return
            end if
            read(tok, *, iostat=ios) val
            if (ios /= 0) call fail('index list: not an integer: "' // tok // '"')
        end subroutine to_int

    end subroutine parse_index_list

    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=16) :: buf
        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

end module atomic_selection
