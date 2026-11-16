! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Atom selection helpers (docs/design/contracts.md §2 for parse_index_list,
!> §6 and §6.1 for the geometric rules).
!>
!> Pipeline role: used by the atomselect CLI to turn user requests (index
!> lists, slabs along a normal, boxes, top/bottom layers) into the atom lists
!> stored as regions in regions.nml (see atomic_regions).
!>
!> Conventions: atom indices are 1-based, in structure-file order; positions
!> are Cartesian Å (structure_t%cart); returned index lists are ascending.
!> Library routines never stop: errors come back as ierr /= 0 plus errmsg.
!>
!>  - numeric index lists ("1-12", "15, 17, 22", "1-10, 15-20");
!>  - geometric selection: projection on an arbitrary normal, s_i = r_i . n_hat
!>    (Cartesian Å; the normal is given by its CARTESIAN components and is
!>    never assumed to be z), range and Cartesian box selections;
!>  - layer detection along the normal (single-linkage on sorted s), manual
!>    reference layers and top/bottom surface regions.
!>
!> Periodicity: projected coordinates are NOT wrapped. Layer detection assumes
!> the slab is not cut by the cell boundary along n_hat; `check_periodic_split`
!> detects the opposite case so callers can warn (nothing is re-wrapped).
module atomic_selection
    use types,        only: dp
    use structure_io, only: structure_t
    implicit none
    private

    public :: parse_index_list
    public :: project_positions, select_by_range, select_by_box
    public :: detect_layers, check_periodic_split, assign_to_references
    public :: surface_layers

    real(dp), parameter :: NORMAL_EPS = 1.0e-12_dp   ! |normal| below this is treated as zero

    character(len=*), parameter :: DIGITS = '0123456789'
    integer, parameter :: MAX_DIGITS = 9     ! keeps values well inside default integer

contains

    ! ==================================================================
    ! Index lists
    ! ==================================================================

    !> Parses a comma-separated list of 1-based indices and inclusive ranges
    !> a-b (e.g. "1-10, 15, 17-20"; blanks around items are ignored) into a
    !> sorted ascending `list`. Indices must lie in 1..nmax and may not repeat
    !> (also across items); duplicates are reported, not removed. Only
    !> non-negative decimal integers are accepted (no signs).
    !> Never stops: on error ierr = 1, errmsg says why and `list` is size 0.
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

        ! Walk the text item by item (items are separated by commas). `seen`
        ! marks the indices already taken, which both detects duplicates and
        ! yields the sorted result at the end.
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

            ! single index "n" or range "a-b"
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
            ! a trailing comma leaves nothing after it
            if (pos > len(text)) then
                call fail('index list: empty item at position ' // itoa(k + 1))
                return
            end if
        end do

        list = pack([(i, i = 1, nmax)], seen)

    contains

        !> Sets the error state and resets `list` to size 0.
        subroutine fail(msg)
            character(len=*), intent(in) :: msg
            ierr = 1
            errmsg = msg
            if (allocated(list)) deallocate(list)
            allocate(list(0))
        end subroutine fail

        !> Converts a digits-only token; calls fail otherwise (so ierr is set).
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

    ! ==================================================================
    ! Geometric selection
    ! ==================================================================

    !> Projects every atom on the normal: s(i) = r_i . n_hat (Å, signed, measured
    !> from the Cartesian origin), n_hat = normal/|normal|. `normal` holds
    !> CARTESIAN components (not fractional) and need not be unit length.
    !> s has size st%nat. Error (ierr = 1) if |normal| < 1e-12; s is then zero.
    subroutine project_positions(st, normal, s, ierr, errmsg)
        type(structure_t),     intent(in)  :: st
        real(dp),              intent(in)  :: normal(3)
        real(dp), allocatable, intent(out) :: s(:)
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg
        real(dp) :: nhat(3), nrm
        integer :: i

        ierr = 0
        errmsg = ''
        allocate(s(st%nat))
        s = 0.0_dp
        nrm = norm2(normal)
        if (nrm < NORMAL_EPS) then
            ierr = 1
            errmsg = 'projection: normal vector has zero length'
            return
        end if
        nhat = normal / nrm
        do i = 1, st%nat
            s(i) = dot_product(st%cart(:, i), nhat)
        end do
    end subroutine project_positions

    !> Atoms with smin <= s(i) <= smax (inclusive, same units as s),
    !> ascending indices, possibly empty. Never fails.
    subroutine select_by_range(s, smin, smax, list)
        real(dp),             intent(in)  :: s(:)
        real(dp),             intent(in)  :: smin, smax
        integer, allocatable, intent(out) :: list(:)
        integer :: i
        list = pack([(i, i = 1, size(s))], s >= smin .and. s <= smax)
    end subroutine select_by_range

    !> Atoms with lo(k) <= r_i(k) <= hi(k) for k = x,y,z (Cartesian Å);
    !> use -huge/huge for unbounded limits. Ascending indices (possibly empty).
    subroutine select_by_box(st, lo, hi, list)
        type(structure_t),    intent(in)  :: st
        real(dp),             intent(in)  :: lo(3), hi(3)
        integer, allocatable, intent(out) :: list(:)
        logical, allocatable :: inside(:)
        integer :: i
        allocate(inside(st%nat))
        do i = 1, st%nat
            inside(i) = all(st%cart(:, i) >= lo) .and. all(st%cart(:, i) <= hi)
        end do
        list = pack([(i, i = 1, st%nat)], inside)
    end subroutine select_by_box

    ! ==================================================================
    ! Layers
    ! ==================================================================

    !> Single-linkage layer detection on projected coordinates (Å): atoms are
    !> sorted by s and a new layer starts when s(next) - s(prev) > tol, so a
    !> layer may be wider than tol if its atoms chain together. Layers are
    !> numbered 1..nlayer in ASCENDING s; layer_of(i) is the layer of atom i
    !> and layer_center(l) = mean s of layer l. tol > 0 required (ierr = 1
    !> otherwise; outputs are then empty/zero). Empty s gives nlayer = 0.
    !> s is not wrapped (see module notes).
    subroutine detect_layers(s, tol, layer_of, nlayer, layer_center, ierr, errmsg)
        real(dp),              intent(in)  :: s(:)
        real(dp),              intent(in)  :: tol
        integer,  allocatable, intent(out) :: layer_of(:)
        integer,               intent(out) :: nlayer
        real(dp), allocatable, intent(out) :: layer_center(:)
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg
        integer, allocatable :: idx(:), cnt(:)
        integer :: n, k, l

        ierr = 0
        errmsg = ''
        n = size(s)
        nlayer = 0
        allocate(layer_of(n))
        layer_of = 0
        allocate(layer_center(0))
        if (.not. (tol > 0.0_dp)) then
            ierr = 1
            errmsg = 'layers: tolerance must be > 0'
            return
        end if
        if (n == 0) return

        ! 1. Sort atoms by s and cut into layers at gaps larger than tol
        idx = sort_index(s)
        nlayer = 1
        layer_of(idx(1)) = 1
        do k = 2, n
            if (s(idx(k)) - s(idx(k - 1)) > tol) nlayer = nlayer + 1
            layer_of(idx(k)) = nlayer
        end do

        ! 2. Layer centers = mean s of the members
        deallocate(layer_center)
        allocate(layer_center(nlayer), cnt(nlayer))
        layer_center = 0.0_dp
        cnt = 0
        do k = 1, n
            l = layer_of(k)
            layer_center(l) = layer_center(l) + s(k)
            cnt(l) = cnt(l) + 1
        end do
        layer_center = layer_center / real(cnt, dp)
    end subroutine detect_layers

    !> warn = .true. if the slab appears continuous across a periodic boundary
    !> along n_hat: for every lattice vector a_i with p = |a_i . n_hat| > tol,
    !> the periodic image of the lowest atom lies at min(s) + p; if
    !> min(s) + p - max(s) <= tol the image joins the top of the slab (or
    !> interleaves with it) and layer detection on unwrapped s is unreliable.
    !> Lattice vectors with p <= tol (in-plane) shift layers onto themselves
    !> within tolerance and are ignored. Zero normal or empty s -> .false.
    subroutine check_periodic_split(st, normal, s, tol, warn)
        type(structure_t), intent(in)  :: st
        real(dp),          intent(in)  :: normal(3)
        real(dp),          intent(in)  :: s(:)
        real(dp),          intent(in)  :: tol
        logical,           intent(out) :: warn
        real(dp) :: nhat(3), nrm, p, smin, smax
        integer :: i

        warn = .false.
        nrm = norm2(normal)
        if (nrm < NORMAL_EPS .or. size(s) == 0) return
        nhat = normal / nrm
        smin = minval(s)
        smax = maxval(s)
        do i = 1, 3
            p = abs(dot_product(st%lattice(i, :), nhat))
            if (p <= tol) cycle
            if (smin + p - smax <= tol) warn = .true.
        end do
    end subroutine check_periodic_split

    !> Assigns every atom to the nearest manual reference center (same units
    !> as s, Å). Layers are numbered by ASCENDING center (layer_of(i) = rank
    !> of the chosen center), so numbering matches detect_layers; a center
    !> with no atoms still occupies its rank. Errors (ierr = 1, layer_of = 0): tol <= 0, no centers,
    !> duplicated centers, or an atom farther than tol from every center
    !> (ambiguous geometry is never forced).
    subroutine assign_to_references(s, ref_centers, tol, layer_of, ierr, errmsg)
        real(dp),             intent(in)  :: s(:)
        real(dp),             intent(in)  :: ref_centers(:)
        real(dp),             intent(in)  :: tol
        integer, allocatable, intent(out) :: layer_of(:)
        integer,              intent(out) :: ierr
        character(len=*),     intent(out) :: errmsg
        integer, allocatable :: idx(:)
        real(dp), allocatable :: c(:)
        real(dp) :: d, dbest
        integer :: i, k, kbest, nref
        character(len=32) :: b1, b2

        ierr = 0
        errmsg = ''
        allocate(layer_of(size(s)))
        layer_of = 0
        nref = size(ref_centers)
        if (.not. (tol > 0.0_dp)) then
            ierr = 1
            errmsg = 'layers: tolerance must be > 0'
            return
        end if
        if (nref == 0) then
            ierr = 1
            errmsg = 'layers: no reference centers given'
            return
        end if
        idx = sort_index(ref_centers)
        c = ref_centers(idx)
        do k = 2, nref
            if (c(k) == c(k - 1)) then
                write(b1, '(F12.6)') c(k)
                ierr = 1
                errmsg = 'layers: duplicated reference center ' // trim(adjustl(b1))
                return
            end if
        end do

        ! Nearest center for each atom, rejecting atoms too far from all of them
        do i = 1, size(s)
            kbest = 0
            dbest = huge(1.0_dp)
            do k = 1, nref
                d = abs(s(i) - c(k))
                if (d < dbest) then
                    dbest = d
                    kbest = k
                end if
            end do
            if (dbest > tol) then
                write(b1, '(F12.6)') s(i)
                write(b2, '(F12.6)') tol
                ierr = 1
                errmsg = 'layers: atom ' // itoa(i) // ' (projected ' // trim(adjustl(b1)) // &
                         ') is farther than tol = ' // trim(adjustl(b2)) // ' from every reference center'
                layer_of = 0
                return
            end if
            layer_of(i) = kbest
        end do
    end subroutine assign_to_references

    !> Surface regions from a layer assignment (layers numbered by ascending s):
    !> top_list = atoms of the ntop highest layers (nlayer-ntop+1..nlayer),
    !> bottom_list = atoms of the nbottom lowest layers (1..nbottom). ntop or
    !> nbottom may be 0 (not requested) but not both; negative values and
    !> ntop + nbottom > nlayer are errors (ierr = 1, lists empty), so top and
    !> bottom never overlap. Lists sorted ascending; atoms with layer_of < 1
    !> (unassigned) are ignored.
    subroutine surface_layers(layer_of, nlayer, ntop, nbottom, top_list, bottom_list, ierr, errmsg)
        integer,              intent(in)  :: layer_of(:)
        integer,              intent(in)  :: nlayer, ntop, nbottom
        integer, allocatable, intent(out) :: top_list(:), bottom_list(:)
        integer,              intent(out) :: ierr
        character(len=*),     intent(out) :: errmsg
        integer :: i

        ierr = 0
        errmsg = ''
        allocate(top_list(0), bottom_list(0))
        if (ntop < 0 .or. nbottom < 0) then
            ierr = 1
            errmsg = 'surface: number of layers must be >= 1'
            return
        end if
        if (ntop == 0 .and. nbottom == 0) then
            ierr = 1
            errmsg = 'surface: request at least one top or bottom layer'
            return
        end if
        if (ntop + nbottom > nlayer) then
            ierr = 1
            errmsg = 'surface: top (' // itoa(ntop) // ') + bottom (' // itoa(nbottom) // &
                     ') layers exceed the ' // itoa(nlayer) // ' layers found'
            return
        end if
        top_list = pack([(i, i = 1, size(layer_of))], layer_of > nlayer - ntop .and. layer_of >= 1)
        bottom_list = pack([(i, i = 1, size(layer_of))], layer_of >= 1 .and. layer_of <= nbottom)
    end subroutine surface_layers

    ! ==================================================================
    ! PRIVATE HELPERS
    ! ==================================================================

    !> Stable ascending sort permutation (bottom-up merge sort): idx such that
    !> v(idx) is ascending; ties keep their original order. Used instead of an
    !> intrinsic because Fortran has none, and stability keeps results
    !> deterministic for equal coordinates.
    function sort_index(v) result(idx)
        real(dp), intent(in) :: v(:)
        integer, allocatable :: idx(:)
        integer, allocatable :: tmp(:)
        integer :: n, i, width, lo, mid, hi, a, b, k

        n = size(v)
        idx = [(i, i = 1, n)]
        allocate(tmp(n))
        width = 1
        do while (width < n)
            lo = 1
            do while (lo <= n)
                mid = min(lo + width - 1, n)
                hi = min(lo + 2 * width - 1, n)
                a = lo
                b = mid + 1
                k = lo
                do while (k <= hi)
                    if (a <= mid .and. (b > hi .or. v(idx(a)) <= v(idx(min(b, n))))) then
                        tmp(k) = idx(a)
                        a = a + 1
                    else
                        tmp(k) = idx(b)
                        b = b + 1
                    end if
                    k = k + 1
                end do
                lo = lo + 2 * width
            end do
            idx = tmp
            width = 2 * width
        end do
    end function sort_index

    !> Integer to trimmed string (for error messages).
    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=16) :: buf
        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

end module atomic_selection
