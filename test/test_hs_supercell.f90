! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Supercélula de .hsio (hs_supercell), num modelo aleatório não ortogonal
!> (2 átomos com 2 e 3 orbitais, vizinhos até |n_a| = 1):
!>  - tile 3 x 2 x 1 e dobra: autovalores de (H, S) da supercélula em K = união
!>    dos do primitivo em k = (K + m)/n, a 1e-10;
!>  - eixo cortado: nenhum vetor R com componente no eixo, e H(k) do filme =
!>    bloco da supercélula periódica sem os termos que cruzam a borda;
!>  - remoção: átomos e blocos que tocam o removido saem.
program test_hs_supercell
    use, intrinsic :: iso_fortran_env, only: int64
    use types, only: dp
    use testing, only: check, check_close, test_summary
    use HS_io, only: HSBlockSparse
    use hs_supercell, only: hs_supercell_build, hs_supercell_index, bs_hk_dense
    implicit none

    interface
        subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
            import :: dp
            integer, intent(in) :: itype, n, lda, ldb, lwork
            character(len=1), intent(in) :: jobz, uplo
            complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: w(*), rwork(*)
            complex(dp), intent(out) :: work(*)
            integer, intent(out) :: info
        end subroutine zhegv
    end interface

    type(HSBlockSparse) :: hp, sp

    call random_model(hp, sp)
    call fold_test()
    call cut_test()
    call remove_test()
    call test_summary('test_hs_supercell')

contains

    !> Modelo: blocos (i, j, R) para R em {-1,0,1}^3 com |R|_1 <= 2; hermiticidade
    !> H(j,i,-R) = H(i,j,R)^dag; S = I + 0.05 * (aleatório hermitiano).
    subroutine random_model(h, s)
        type(HSBlockSparse), intent(out) :: h, s
        integer :: orbs(2) = [2, 3], r(3), nr, i, j, ib, k, n1, n2, n3, ir, irm, nb
        integer, allocatable :: rl(:, :), bi(:, :)
        complex(dp), allocatable :: hb(:, :, :), sb(:, :, :)
        real(dp) :: x(3, 3, 2), y(3, 3, 2)
        integer :: seed(64)
        seed = 12345
        call random_seed(put=seed(1:size_seed()))
        nr = 0
        allocate(rl(3, 27))
        do n3 = -1, 1
            do n2 = -1, 1
                do n1 = -1, 1
                    if (abs(n1) + abs(n2) + abs(n3) > 2) cycle
                    nr = nr + 1
                    rl(:, nr) = [n1, n2, n3]
                end do
            end do
        end do
        ! blocos para R "positivos" (lexicográfico) e o conjugado em -R
        allocate(bi(3, 4 * nr), hb(3, 3, 4 * nr), sb(3, 3, 4 * nr))
        nb = 0
        do ir = 1, nr
            r = rl(:, ir)
            do irm = 1, nr
                if (all(rl(:, irm) == -r)) exit
            end do
            do i = 1, 2
                do j = 1, 2
                    if (ir > irm) cycle
                    if (ir == irm .and. j < i) cycle
                    call random_number(x); call random_number(y)
                    nb = nb + 1
                    bi(:, nb) = [i, j, ir]
                    hb(:, :, nb) = cmplx(x(:, :, 1) - 0.5_dp, y(:, :, 1) - 0.5_dp, dp)
                    sb(:, :, nb) = 0.05_dp * cmplx(x(:, :, 2) - 0.5_dp, y(:, :, 2) - 0.5_dp, dp)
                    if (ir == irm .and. i == j) then
                        hb(:, :, nb) = 0.5_dp * (hb(:, :, nb) + conjg(transpose(hb(:, :, nb))))
                        sb(:, :, nb) = 0.5_dp * (sb(:, :, nb) + conjg(transpose(sb(:, :, nb))))
                        do k = 1, 3
                            sb(k, k, nb) = sb(k, k, nb) + (1.0_dp, 0.0_dp)
                        end do
                    else
                        nb = nb + 1
                        bi(:, nb) = [j, i, irm]
                        hb(:, :, nb) = conjg(transpose(hb(:, :, nb - 1)))
                        sb(:, :, nb) = conjg(transpose(sb(:, :, nb - 1)))
                    end if
                end do
            end do
        end do
        call pack_bs(h, orbs, rl(:, 1:nr), bi(:, 1:nb), hb(:, :, 1:nb))
        call pack_bs(s, orbs, rl(:, 1:nr), bi(:, 1:nb), sb(:, :, 1:nb))
        ib = 0
    end subroutine random_model

    integer function size_seed()
        call random_seed(size=size_seed)
        size_seed = min(size_seed, 64)
    end function size_seed

    subroutine pack_bs(bs, orbs, rl, bi, blk)
        type(HSBlockSparse), intent(out) :: bs
        integer, intent(in) :: orbs(:), rl(:, :), bi(:, :)
        complex(dp), intent(in) :: blk(:, :, :)
        integer :: k, ni, nj
        bs%nAtoms = size(orbs)
        bs%Orbs = orbs
        bs%NOrbs = sum(orbs)
        bs%Neighbors = size(rl, 2)
        bs%IndexList = rl
        allocate(bs%Degenerecency(bs%Neighbors))
        bs%Degenerecency = 1
        bs%nBlocks = int(size(bi, 2), int64)
        bs%BlockIndex = bi
        allocate(bs%BlockOffset(size(bi, 2) + 1))
        bs%BlockOffset(1) = 0
        do k = 1, size(bi, 2)
            bs%BlockOffset(k + 1) = bs%BlockOffset(k) + orbs(bi(1, k)) * orbs(bi(2, k))
        end do
        allocate(bs%BlockData(bs%BlockOffset(size(bi, 2) + 1)))
        do k = 1, size(bi, 2)
            ni = orbs(bi(1, k)); nj = orbs(bi(2, k))
            bs%BlockData(bs%BlockOffset(k) + 1:bs%BlockOffset(k + 1)) = reshape(blk(1:ni, 1:nj, k), [ni * nj])
        end do
    end subroutine pack_bs

    subroutine geig(h, s, w)
        complex(dp), intent(in) :: h(:, :), s(:, :)
        real(dp), intent(out) :: w(:)
        complex(dp), allocatable :: a(:, :), b(:, :), work(:)
        real(dp), allocatable :: rwork(:)
        integer :: n, info
        n = size(h, 1)
        a = h; b = s
        allocate(work(4 * n), rwork(3 * n))
        call zhegv(1, 'N', 'U', n, a, n, b, n, w, work, 4 * n, rwork, info)
        call check(info == 0, 'zhegv')
    end subroutine geig

    subroutine fold_test()
        type(HSBlockSparse) :: hs, ss
        integer, allocatable :: idx(:), idx2(:)
        logical, allocatable :: keep(:)
        complex(dp), allocatable :: hk(:, :), sk(:, :)
        real(dp), allocatable :: wsup(:), wprim(:), wall(:)
        integer :: n(3) = [3, 2, 1], ierr, m1, m2, cnt, ik
        character(len=:), allocatable :: errmsg
        real(dp) :: kk(3), dmax, kpts(3, 3)
        allocate(keep(2 * product(n)))
        keep = .true.
        call hs_supercell_build(hp, n, [.true., .true., .true.], keep, hs, idx, ierr, errmsg)
        call check(ierr == 0, 'tile H: ' // errmsg)
        call hs_supercell_build(sp, n, [.true., .true., .true.], keep, ss, idx2, ierr, errmsg)
        call check(ierr == 0 .and. all(hs%IndexList == ss%IndexList), 'tile S, same cell list')
        call check(hs%nAtoms == 12 .and. hs%NOrbs == 30, 'tile 3x2x1: 12 atoms, 30 orbitals')
        kpts(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        kpts(:, 2) = [0.17_dp, -0.31_dp, 0.29_dp]
        kpts(:, 3) = [0.5_dp, 0.5_dp, 0.5_dp]
        dmax = 0.0_dp
        allocate(wsup(30), wprim(5))
        do ik = 1, 3
            call bs_hk_dense(hs, kpts(:, ik), hk)
            call bs_hk_dense(ss, kpts(:, ik), sk)
            call geig(hk, sk, wsup)
            allocate(wall(0))
            do m1 = 0, n(1) - 1
                do m2 = 0, n(2) - 1
                    kk = [(kpts(1, ik) + m1) / n(1), (kpts(2, ik) + m2) / n(2), kpts(3, ik)]
                    call bs_hk_dense(hp, kk, hk)
                    call bs_hk_dense(sp, kk, sk)
                    call geig(hk, sk, wprim)
                    wall = [wall, wprim]
                end do
            end do
            call sort(wall)
            cnt = size(wall)
            call check(cnt == 30, 'fold: 30 eigenvalues')
            dmax = max(dmax, maxval(abs(wall - wsup)))
            deallocate(wall)
        end do
        call check_close(dmax, 0.0_dp, 1.0e-10_dp, 'tile 3x2x1 then fold = primitive bands (1e-10)')
    end subroutine fold_test

    subroutine cut_test()
        type(HSBlockSparse) :: hs, hper
        integer, allocatable :: idx(:)
        logical, allocatable :: keep(:)
        complex(dp), allocatable :: hk(:, :), hpk(:, :)
        integer :: n(3) = [1, 3, 1], ierr, i, j, ib, a, b
        character(len=:), allocatable :: errmsg
        real(dp) :: k(3), dmax
        allocate(keep(6))
        keep = .true.
        call hs_supercell_build(hp, n, [.true., .false., .true.], keep, hs, idx, ierr, errmsg)
        call check(ierr == 0, 'cut axis 2: ' // errmsg)
        call check(all(hs%IndexList(2, :) == 0), 'cut axis: no cell along the cut axis')
        call hs_supercell_build(hp, n, [.true., .true., .true.], keep, hper, idx, ierr, errmsg)
        ! filme em k = H periódico em k2 = 0 menos os termos com R'_2 /= 0; com k2 arbitrário
        ! o filme não depende de k2
        k = [0.21_dp, 0.37_dp, -0.13_dp]
        call bs_hk_dense(hs, k, hk)
        call bs_hk_dense(hs, [k(1), 0.0_dp, k(3)], hpk)
        call check_close(maxval(abs(hk - hpk)), 0.0_dp, 1.0e-14_dp, 'film H(k) independent of k along the cut axis')
        ! blocos do filme = blocos da supercélula periódica com R'_2 = 0
        dmax = 0.0_dp
        do ib = 1, int(hper%nBlocks)
            if (hper%IndexList(2, hper%BlockIndex(3, ib)) /= 0) cycle
            i = hper%BlockIndex(1, ib); j = hper%BlockIndex(2, ib)
            do a = 1, int(hs%nBlocks)
                if (hs%BlockIndex(1, a) /= i .or. hs%BlockIndex(2, a) /= j) cycle
                if (any(hs%IndexList(:, hs%BlockIndex(3, a)) /= hper%IndexList(:, hper%BlockIndex(3, ib)))) cycle
                b = int(hper%BlockOffset(ib + 1) - hper%BlockOffset(ib))
                dmax = max(dmax, maxval(abs(hs%BlockData(hs%BlockOffset(a) + 1:hs%BlockOffset(a) + b) - &
                                            hper%BlockData(hper%BlockOffset(ib) + 1:hper%BlockOffset(ib + 1)))))
            end do
        end do
        call check(count(hper%IndexList(2, hper%BlockIndex(3, :)) == 0) == int(hs%nBlocks), &
                   'film keeps exactly the blocks with R2 = 0')
        call check_close(dmax, 0.0_dp, 0.0_dp, 'film blocks identical to the bulk ones')
    end subroutine cut_test

    subroutine remove_test()
        type(HSBlockSparse) :: hs
        integer, allocatable :: idx(:)
        logical, allocatable :: keep(:)
        integer :: n(3) = [2, 1, 1], ierr, rem
        character(len=:), allocatable :: errmsg
        allocate(keep(4))
        keep = .true.
        rem = hs_supercell_index(2, n, 2, [1, 0, 0])
        keep(rem) = .false.
        call hs_supercell_build(hp, n, [.true., .true., .true.], keep, hs, idx, ierr, errmsg)
        call check(ierr == 0 .and. hs%nAtoms == 3 .and. hs%NOrbs == 7, 'remove one atom: 3 atoms, 7 orbitals')
        call check(rem == 4 .and. idx(rem) == 0 .and. idx(3) == 3, 'renumbering after removal')
        call check(all(hs%BlockIndex(1:2, :) >= 1) .and. all(hs%BlockIndex(1:2, :) <= 3), 'no block on the removed atom')
    end subroutine remove_test

    subroutine sort(x)
        real(dp), intent(inout) :: x(:)
        integer :: i, j
        real(dp) :: t
        do i = 2, size(x)
            t = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= t) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = t
        end do
    end subroutine sort

end program test_hs_supercell
