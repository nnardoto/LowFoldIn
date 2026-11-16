! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Piecewise-linear k-path in fractional reciprocal coordinates with
!> accumulated distance in Cartesian 1/Å (reciprocal metric,
!> docs/design/contracts.md §3.1).
!>
!> Pipeline position: the spectral engine builds the path once, evaluates
!> the Hamiltonian at every kfrac point, and copies kfrac/kdist/vertex data
!> into spectral_data_t for plotting.
!>
!> Conventions
!>   - kfrac are fractional coordinates (multiples of the reciprocal
!>     lattice vectors b_i); kdist is Cartesian, in 1/Å, kdist(1) = 0.
!>   - recip(3,3) holds b_i as ROWS and already includes the 2*pi factor.
!>   - vertex_k indices are 1-based positions in the k list.
!> Library routines never stop: errors are returned through ierr/errmsg.
module kpath
    use types, only: dp
    implicit none
    private

    integer, parameter, public :: KLABEL_LEN = 16

    !> A generated k-path. All arrays are allocated only on success.
    type, public :: kpath_t
        integer :: nk = 0, nvertex = 0
        real(dp), allocatable :: kfrac(:,:)     ! (3,nk)
        real(dp), allocatable :: kdist(:)       ! (nk) 1/Å, kdist(1)=0
        integer,  allocatable :: vertex_k(:)    ! (nvertex) index of vertex in k list
        character(len=KLABEL_LEN), allocatable :: label(:)   ! (nvertex)
    end type kpath_t

    public :: build_kpath

contains

    ! ==================================================================
    ! Public API
    ! ==================================================================

    !> Build the k-path.
    !> vertex(3,nvertex) fractional; npoints(nvertex-1) >= 1 intervals per
    !> segment; label(nvertex); recip(3,3) with rows b_i (1/Å, includes 2π).
    !> Segment s: vertex_s + j/npoints(s)*(vertex_{s+1}-vertex_s),
    !> j = 0..npoints(s)-1; the final vertex is appended once.
    !> nk = sum(npoints) + 1. kdist increments = norm2(matmul(dk_frac, recip)).
    !> Points are uniformly spaced within each segment and junction vertices
    !> are not duplicated, so kdist is continuous across segments.
    !> ierr codes: 1 vertex not (3,n); 2 nvertex < 2; 3 npoints size;
    !> 4 label size; 5 npoints < 1; 6 zero-length segment; 7 label longer
    !> than KLABEL_LEN. On error errmsg is set and path is left empty.
    subroutine build_kpath(vertex, npoints, label, recip, path, ierr, errmsg)
        real(dp),         intent(in)  :: vertex(:,:)
        integer,          intent(in)  :: npoints(:)
        character(len=*), intent(in)  :: label(:)
        real(dp),         intent(in)  :: recip(3,3)
        type(kpath_t),    intent(out) :: path
        integer,          intent(out) :: ierr
        character(len=*), intent(out) :: errmsg

        integer  :: nvertex, nk, s, j, k
        real(dp) :: dv(3), seglen, step, dist0

        ierr = 0
        errmsg = ''
        nvertex = size(vertex, 2)

        ! 1. Validate all inputs before allocating anything.
        if (size(vertex, 1) /= 3) then
            ierr = 1
            errmsg = 'build_kpath: vertex must have shape (3,nvertex)'
            return
        end if
        if (nvertex < 2) then
            ierr = 2
            errmsg = 'build_kpath: at least 2 vertices are required'
            return
        end if
        if (size(npoints) /= nvertex - 1) then
            ierr = 3
            write(errmsg, '(A,I0,A,I0)') 'build_kpath: npoints must have nvertex-1 = ', &
                nvertex - 1, ' values, got ', size(npoints)
            return
        end if
        if (size(label) /= nvertex) then
            ierr = 4
            write(errmsg, '(A,I0,A,I0)') 'build_kpath: label must have nvertex = ', &
                nvertex, ' values, got ', size(label)
            return
        end if
        do s = 1, nvertex - 1
            if (npoints(s) < 1) then
                ierr = 5
                write(errmsg, '(A,I0,A,I0)') 'build_kpath: npoints(', s, ') must be >= 1, got ', &
                    npoints(s)
                return
            end if
            if (all(vertex(:, s + 1) == vertex(:, s))) then
                ierr = 6
                write(errmsg, '(A,I0,A,I0)') 'build_kpath: zero-length segment between vertices ', &
                    s, ' and ', s + 1
                return
            end if
        end do
        do s = 1, nvertex
            if (len_trim(label(s)) > KLABEL_LEN) then
                ierr = 7
                write(errmsg, '(A,I0,A,I0,A)') 'build_kpath: label(', s, ') longer than ', &
                    KLABEL_LEN, ' characters'
                return
            end if
        end do

        ! 2. Allocate: each segment contributes npoints(s) points (its start
        !    included, its end excluded); the last vertex adds one more.
        nk = sum(npoints) + 1
        path%nk = nk
        path%nvertex = nvertex
        allocate(path%kfrac(3, nk), path%kdist(nk), path%vertex_k(nvertex), path%label(nvertex))
        do s = 1, nvertex
            path%label(s) = label(s)
        end do

        ! 3. Fill points segment by segment. dist0 is the accumulated
        !    distance at the start of segment s. The length is measured in
        !    Cartesian space (dv_frac * recip, rows = b_i) since fractional
        !    steps along different b_i are not equivalent.
        k = 0
        dist0 = 0.0_dp
        do s = 1, nvertex - 1
            dv = vertex(:, s + 1) - vertex(:, s)
            seglen = norm2(matmul(dv, recip))
            ! Guard for a degenerate recip (e.g. zero matrix) that the exact
            ! vertex-equality check above cannot detect.
            if (seglen <= 0.0_dp) then
                ierr = 6
                write(errmsg, '(A,I0,A,I0)') 'build_kpath: zero-length segment between vertices ', &
                    s, ' and ', s + 1
                call reset_path(path)
                return
            end if
            step = seglen / real(npoints(s), dp)
            path%vertex_k(s) = k + 1
            do j = 0, npoints(s) - 1
                k = k + 1
                path%kfrac(:, k) = vertex(:, s) + (real(j, dp) / real(npoints(s), dp)) * dv
                path%kdist(k) = dist0 + real(j, dp) * step
            end do
            dist0 = dist0 + seglen
        end do
        ! 4. Final vertex, appended once (closes the last segment).
        k = k + 1
        path%kfrac(:, k) = vertex(:, nvertex)
        path%kdist(k) = dist0
        path%vertex_k(nvertex) = k
    end subroutine build_kpath

    ! ==================================================================
    ! Internals
    ! ==================================================================

    !> Deallocate everything and zero the counts (used on late errors).
    subroutine reset_path(path)
        type(kpath_t), intent(inout) :: path
        if (allocated(path%kfrac))    deallocate(path%kfrac)
        if (allocated(path%kdist))    deallocate(path%kdist)
        if (allocated(path%vertex_k)) deallocate(path%vertex_k)
        if (allocated(path%label))    deallocate(path%label)
        path%nk = 0
        path%nvertex = 0
    end subroutine reset_path

end module kpath
