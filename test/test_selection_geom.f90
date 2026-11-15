! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Tests for the geometric part of atomic_selection: projection on an
!> arbitrary normal, range/box selection, layer detection, reference layers,
!> surface layers and the periodic-split warning (docs/design/contracts.md §6).
program test_selection_geom
    use types, only: dp
    use testing, only: check, check_close, test_summary
    use structure_io, only: structure_t
    use atomic_selection, only: project_positions, select_by_range, select_by_box, &
                                detect_layers, check_periodic_split, &
                                assign_to_references, surface_layers
    implicit none

    real(dp), parameter :: TOL = 1.0e-10_dp
    real(dp), parameter :: R2 = 1.0_dp / sqrt(2.0_dp)

    call test_z_slab()
    call test_x_slab()
    call test_tilted_hex()
    call test_tolerance()
    call test_references()
    call test_range_box()
    call test_zero_normal()
    call test_periodic_split()
    call test_surface_errors()

    call test_summary('test_selection_geom')

contains

    ! ------------------------------------------------------------------
    subroutine make(st, lattice, pos)
        type(structure_t), intent(out) :: st
        real(dp),          intent(in)  :: lattice(3,3)   ! rows = a_i
        real(dp),          intent(in)  :: pos(:,:)       ! (3,nat)
        integer :: i
        st%nat = size(pos, 2)
        st%lattice = lattice
        st%cart = pos
        allocate(st%label(st%nat))
        do i = 1, st%nat
            if (mod(i, 2) == 1) then
                st%label(i) = 'Ga'
            else
                st%label(i) = 'N'
            end if
        end do
    end subroutine make

    function ortho(a, b, c) result(lat)
        real(dp), intent(in) :: a, b, c
        real(dp) :: lat(3,3)
        lat = 0.0_dp
        lat(1,1) = a; lat(2,2) = b; lat(3,3) = c
    end function ortho

    logical function same(a, b)
        integer, intent(in) :: a(:), b(:)
        same = size(a) == size(b)
        if (same) same = all(a == b)
    end function same

    ! ------------------------------------------------------------------
    !> (a) orthogonal slab along z, 4 layers, top layer first in file order.
    subroutine test_z_slab()
        type(structure_t) :: st
        real(dp), allocatable :: s(:), c(:)
        integer, allocatable :: lay(:), top(:), bot(:)
        integer :: nl, ierr
        character(len=256) :: msg

        call make(st, ortho(4.0_dp, 4.0_dp, 30.0_dp), reshape([ &
            0.0_dp, 0.0_dp, 9.5_dp,   2.0_dp, 2.0_dp, 9.5_dp, &
            0.0_dp, 2.0_dp, 7.0_dp,   2.0_dp, 0.0_dp, 7.0_dp, &
            0.0_dp, 0.0_dp, 4.5_dp,   2.0_dp, 2.0_dp, 4.5_dp, &
            0.0_dp, 2.0_dp, 2.0_dp,   2.0_dp, 0.0_dp, 2.0_dp], [3, 8]))
        call project_positions(st, [0.0_dp, 0.0_dp, 2.0_dp], s, ierr, msg)
        call check(ierr == 0, 'z slab: projection ok')
        call check_close(s(1), 9.5_dp, TOL, 'z slab: s normalised')
        call detect_layers(s, 0.5_dp, lay, nl, c, ierr, msg)
        call check(ierr == 0 .and. nl == 4, 'z slab: 4 layers')
        call check(same(lay, [4, 4, 3, 3, 2, 2, 1, 1]), 'z slab: ascending layer numbering')
        call check_close(c(1), 2.0_dp, TOL, 'z slab: center layer 1')
        call check_close(c(4), 9.5_dp, TOL, 'z slab: center layer 4')
        call surface_layers(lay, nl, 1, 1, top, bot, ierr, msg)
        call check(ierr == 0, 'z slab: surface ok')
        call check(same(top, [1, 2]), 'z slab: top = highest s')
        call check(same(bot, [7, 8]), 'z slab: bottom = lowest s')
        call surface_layers(lay, nl, 2, 0, top, bot, ierr, msg)
        call check(ierr == 0 .and. same(top, [1, 2, 3, 4]) .and. size(bot) == 0, 'z slab: top 2 only')
        ! reversed normal: top becomes the former bottom
        call project_positions(st, [0.0_dp, 0.0_dp, -1.0_dp], s, ierr, msg)
        call detect_layers(s, 0.5_dp, lay, nl, c, ierr, msg)
        call surface_layers(lay, nl, 1, 1, top, bot, ierr, msg)
        call check(same(top, [7, 8]) .and. same(bot, [1, 2]), 'z slab: -z normal swaps top/bottom')
    end subroutine test_z_slab

    ! ------------------------------------------------------------------
    !> (b) layers along x: detected with (1,0,0), not with z.
    subroutine test_x_slab()
        type(structure_t) :: st
        real(dp), allocatable :: s(:), c(:)
        integer, allocatable :: lay(:)
        integer :: nl, ierr
        character(len=256) :: msg

        call make(st, ortho(30.0_dp, 4.0_dp, 4.0_dp), reshape([ &
            1.0_dp, 0.0_dp, 0.0_dp,   1.0_dp, 2.0_dp, 0.3_dp, &
            3.0_dp, 0.0_dp, 0.6_dp,   3.0_dp, 2.0_dp, 0.9_dp, &
            5.0_dp, 0.0_dp, 1.2_dp,   5.0_dp, 2.0_dp, 1.5_dp], [3, 6]))
        call project_positions(st, [1.0_dp, 0.0_dp, 0.0_dp], s, ierr, msg)
        call detect_layers(s, 0.5_dp, lay, nl, c, ierr, msg)
        call check(ierr == 0 .and. nl == 3, 'x slab: 3 layers along x')
        call check(same(lay, [1, 1, 2, 2, 3, 3]), 'x slab: grouping along x')
        call project_positions(st, [0.0_dp, 0.0_dp, 1.0_dp], s, ierr, msg)
        call detect_layers(s, 0.5_dp, lay, nl, c, ierr, msg)
        call check(ierr == 0 .and. nl == 1, 'x slab: no layering along z')
    end subroutine test_x_slab

    ! ------------------------------------------------------------------
    !> (c) hexagonal (non-orthogonal) in-plane cell perpendicular to the tilted
    !> normal n = (1,1,0)/sqrt(2); tilted a3; buckled layers.
    subroutine test_tilted_hex()
        type(structure_t) :: st
        real(dp) :: n(3), u(3), w(3), lat(3,3), pos(3,8), sl(4)
        real(dp), allocatable :: s(:), c(:)
        integer, allocatable :: lay(:), top(:), bot(:)
        integer :: nl, ierr, k
        logical :: warn
        character(len=256) :: msg

        n = [R2, R2, 0.0_dp]
        u = [R2, -R2, 0.0_dp]
        w = [0.0_dp, 0.0_dp, 1.0_dp]
        lat(1,:) = 3.2_dp * u
        lat(2,:) = -1.6_dp * u + 2.771281_dp * w
        lat(3,:) = 20.0_dp * n + 1.5_dp * u
        sl = [9.5_dp, 7.0_dp, 4.5_dp, 2.0_dp]
        do k = 1, 4
            pos(:, 2*k-1) = sl(k) * n + 0.4_dp * (k - 1) * u
            pos(:, 2*k)   = (sl(k) + 0.1_dp) * n + (1.6_dp + 0.4_dp * (k - 1)) * u + 0.92376_dp * w
        end do
        call make(st, lat, pos)

        call project_positions(st, [1.0_dp, 1.0_dp, 0.0_dp], s, ierr, msg)
        call check_close(s(2), 9.6_dp, 1.0e-9_dp, 'tilted: projection of atom 2')
        call detect_layers(s, 0.5_dp, lay, nl, c, ierr, msg)
        call check(ierr == 0 .and. nl == 4, 'tilted: 4 layers along (1,1,0)')
        call check(same(lay, [4, 4, 3, 3, 2, 2, 1, 1]), 'tilted: grouping')
        call check_close(c(1), 2.05_dp, 1.0e-9_dp, 'tilted: center = mean s')
        call surface_layers(lay, nl, 1, 1, top, bot, ierr, msg)
        call check(same(top, [1, 2]) .and. same(bot, [7, 8]), 'tilted: top/bottom')
        call check_periodic_split(st, [1.0_dp, 1.0_dp, 0.0_dp], s, 0.5_dp, warn)
        call check(.not. warn, 'tilted: centered slab, no split warning')

        call project_positions(st, [0.0_dp, 0.0_dp, 1.0_dp], s, ierr, msg)
        call detect_layers(s, 0.5_dp, lay, nl, c, ierr, msg)
        call check(nl == 2 .and. same(lay, [1, 2, 1, 2, 1, 2, 1, 2]), 'tilted: z gives a different grouping')
    end subroutine test_tilted_hex

    ! ------------------------------------------------------------------
    !> (d) tolerance: too large merges layers, too small splits a buckled layer.
    subroutine test_tolerance()
        real(dp), allocatable :: c(:)
        integer, allocatable :: lay(:)
        integer :: nl, ierr
        character(len=256) :: msg
        real(dp) :: s(6)

        s = [0.0_dp, 0.3_dp, 2.0_dp, 2.2_dp, 4.0_dp, 4.0_dp]   ! layer 1 buckled by 0.3
        call detect_layers(s, 0.5_dp, lay, nl, c, ierr, msg)
        call check(nl == 3, 'tol 0.5: 3 layers')
        call detect_layers(s, 2.5_dp, lay, nl, c, ierr, msg)
        call check(nl == 1, 'tol 2.5: merged into 1 layer')
        call detect_layers(s, 0.25_dp, lay, nl, c, ierr, msg)
        call check(nl == 4 .and. same(lay, [1, 2, 3, 3, 4, 4]), 'tol 0.25: buckled layer split')
        call detect_layers(s, 0.0_dp, lay, nl, c, ierr, msg)
        call check(ierr /= 0, 'tol 0: error')
    end subroutine test_tolerance

    ! ------------------------------------------------------------------
    !> (e) manual references, numbering by ascending center, ambiguous atom error.
    subroutine test_references()
        integer, allocatable :: lay(:)
        integer :: ierr
        character(len=256) :: msg
        real(dp) :: s(5)

        s = [9.4_dp, 2.1_dp, 5.0_dp, 1.9_dp, 9.6_dp]
        call assign_to_references(s, [9.5_dp, 2.0_dp, 5.0_dp], 0.5_dp, lay, ierr, msg)
        call check(ierr == 0, 'refs: ok')
        call check(same(lay, [3, 1, 2, 1, 3]), 'refs: nearest center, ascending numbering')
        s(3) = 3.5_dp
        call assign_to_references(s, [9.5_dp, 2.0_dp, 5.0_dp], 0.5_dp, lay, ierr, msg)
        call check(ierr /= 0, 'refs: atom far from all centers -> error')
        call check(index(msg, 'atom 3') > 0, 'refs: message names the atom')
        call assign_to_references(s, [2.0_dp, 2.0_dp], 0.5_dp, lay, ierr, msg)
        call check(ierr /= 0, 'refs: duplicated centers -> error')
        call assign_to_references(s, [real(dp) ::], 0.5_dp, lay, ierr, msg)
        call check(ierr /= 0, 'refs: empty centers -> error')
    end subroutine test_references

    ! ------------------------------------------------------------------
    !> (f) select_by_range / select_by_box.
    subroutine test_range_box()
        type(structure_t) :: st
        integer, allocatable :: list(:)
        real(dp), parameter :: BIG = huge(1.0_dp)

        call select_by_range([3.0_dp, 1.0_dp, 2.0_dp, 5.0_dp], 1.0_dp, 3.0_dp, list)
        call check(same(list, [1, 2, 3]), 'range: inclusive limits, sorted')
        call select_by_range([3.0_dp, 1.0_dp], 4.0_dp, 5.0_dp, list)
        call check(size(list) == 0, 'range: empty')

        call make(st, ortho(10.0_dp, 10.0_dp, 10.0_dp), reshape([ &
            1.0_dp, 1.0_dp, 1.0_dp,   5.0_dp, 1.0_dp, 1.0_dp, &
            1.0_dp, 5.0_dp, 9.0_dp,   5.0_dp, 5.0_dp, 5.0_dp], [3, 4]))
        call select_by_box(st, [0.0_dp, -BIG, -BIG], [2.0_dp, BIG, BIG], list)
        call check(same(list, [1, 3]), 'box: x slab, y/z unbounded')
        call select_by_box(st, [-BIG, -BIG, 4.0_dp], [BIG, BIG, 6.0_dp], list)
        call check(same(list, [4]), 'box: z window')
        call select_by_box(st, [0.0_dp, 0.0_dp, 0.0_dp], [5.0_dp, 5.0_dp, 5.0_dp], list)
        call check(same(list, [1, 2, 4]), 'box: closed limits')
    end subroutine test_range_box

    ! ------------------------------------------------------------------
    !> (g) zero normal.
    subroutine test_zero_normal()
        type(structure_t) :: st
        real(dp), allocatable :: s(:)
        integer :: ierr
        character(len=256) :: msg
        call make(st, ortho(1.0_dp, 1.0_dp, 1.0_dp), reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]))
        call project_positions(st, [0.0_dp, 0.0_dp, 1.0e-14_dp], s, ierr, msg)
        call check(ierr /= 0, 'zero normal: error')
        call check(index(msg, 'zero') > 0, 'zero normal: message')
    end subroutine test_zero_normal

    ! ------------------------------------------------------------------
    !> (h) periodic split: slab across the z boundary warns; centered slab does not.
    subroutine test_periodic_split()
        type(structure_t) :: st
        real(dp), allocatable :: s(:)
        integer :: ierr
        logical :: warn
        character(len=256) :: msg

        call make(st, ortho(4.0_dp, 4.0_dp, 20.0_dp), reshape([ &
            0.0_dp, 0.0_dp, 0.5_dp,   2.0_dp, 2.0_dp, 19.5_dp, &
            0.0_dp, 2.0_dp, 17.0_dp,  2.0_dp, 0.0_dp, 3.0_dp], [3, 4]))
        call project_positions(st, [0.0_dp, 0.0_dp, 1.0_dp], s, ierr, msg)
        call check_periodic_split(st, [0.0_dp, 0.0_dp, 1.0_dp], s, 1.5_dp, warn)
        call check(warn, 'split: slab across boundary warns')

        call make(st, ortho(4.0_dp, 4.0_dp, 20.0_dp), reshape([ &
            0.0_dp, 0.0_dp, 8.0_dp,   2.0_dp, 2.0_dp, 10.0_dp, &
            0.0_dp, 2.0_dp, 12.0_dp,  2.0_dp, 0.0_dp, 9.0_dp], [3, 4]))
        call project_positions(st, [0.0_dp, 0.0_dp, 1.0_dp], s, ierr, msg)
        call check_periodic_split(st, [0.0_dp, 0.0_dp, 1.0_dp], s, 1.5_dp, warn)
        call check(.not. warn, 'split: centered slab does not warn')
        ! along x (a = 4 A) atoms at x = 0, 2: image gap 2 A <= tol 2 -> warning
        call project_positions(st, [1.0_dp, 0.0_dp, 0.0_dp], s, ierr, msg)
        call check_periodic_split(st, [1.0_dp, 0.0_dp, 0.0_dp], s, 2.0_dp, warn)
        call check(warn, 'split: bulk-like direction warns')
    end subroutine test_periodic_split

    ! ------------------------------------------------------------------
    subroutine test_surface_errors()
        integer, allocatable :: top(:), bot(:)
        integer :: ierr
        character(len=256) :: msg
        call surface_layers([1, 2, 3], 3, 2, 2, top, bot, ierr, msg)
        call check(ierr /= 0, 'surface: top+bottom > nlayer -> error')
        call surface_layers([1, 2, 3], 3, 0, 0, top, bot, ierr, msg)
        call check(ierr /= 0, 'surface: nothing requested -> error')
        call surface_layers([1, 2, 3], 3, -1, 1, top, bot, ierr, msg)
        call check(ierr /= 0, 'surface: negative -> error')
        call surface_layers([1, 2, 3, 3], 3, 2, 1, top, bot, ierr, msg)
        call check(ierr == 0 .and. same(top, [2, 3, 4]) .and. same(bot, [1]), 'surface: 2N == nlayer allowed')
    end subroutine test_surface_errors

end program test_selection_geom
