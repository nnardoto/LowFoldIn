! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Minimal test harness shared by all fpm test programs.
!> Usage: call check(cond, "description"); at the end call test_summary("name")
!> which stops with exit code 1 if any check failed.
module testing
    use types, only: dp
    implicit none
    private
    public :: check, check_close, test_summary

    integer :: n_pass = 0, n_fail = 0

contains

    subroutine check(cond, msg)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: msg
        if (cond) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            print '(A, A)', "  FAIL: ", msg
        end if
    end subroutine check

    subroutine check_close(a, b, tol, msg)
        real(dp),         intent(in) :: a, b, tol
        character(len=*), intent(in) :: msg
        if (abs(a - b) <= tol) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            print '(A, A, A, ES12.4, A, ES12.4, A, ES10.2)', "  FAIL: ", msg, ": got ", a, &
                " expected ", b, " tol ", tol
        end if
    end subroutine check_close

    subroutine test_summary(name)
        character(len=*), intent(in) :: name
        print '(A, A, A, I0, A, I0, A)', "[", name, "] ", n_pass, " passed, ", n_fail, " failed"
        if (n_fail > 0) stop 1
    end subroutine test_summary

end module testing
