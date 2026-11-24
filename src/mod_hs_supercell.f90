! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Supercélula diagonal de um H(R)/S(R) bloco-esparso (.hsio), com eixos
!> cortados (filme/slab) e remoção de átomos.
!>
!> Supercélula: A_a = n_a a_a (a = 1, 2, 3). O átomo b da cópia t = (t1, t2,
!> t3), 0 <= t_a < n_a, fica em r_b + sum_a t_a a_a; índice na supercélula
!>   I(b, t) = b + nat * (t1 + n1 * (t2 + n2 * t3))
!> (átomo mais rápido, depois t1, t2, t3), antes da remoção; os removidos
!> saem e os demais mantêm a ordem.
!> Blocos: cada bloco (b, b', R) do primitivo gera, para cada cópia t de b, o
!> bloco (I(b, t), I(b', t'), R') com u = t + R, t'_a = mod(u_a, n_a) e
!> R'_a = floor(u_a / n_a) nos eixos periódicos. Num eixo cortado
!> (periodic(a) = .false.) só entram os blocos com 0 <= u_a < n_a (R'_a = 0):
!> o filme/slab mantém H e S do volume entre os átomos que restam, inclusive
!> os termos de três centros no on-site (truncamento simples, sem relaxação).
!> Blocos que tocam um átomo removido são descartados.
!>
!> Com todos os eixos periódicos e sem remoção a dobra é exata: os autovalores
!> de (H, S) da supercélula em K são os do primitivo nos k com n_a k_a = K_a
!> (mod 1) (gauge da rede, H(k) = sum_R H(R) e^{2 pi i k.n}).
!>
!> Rotinas sem stop; erros por ierr/errmsg.
module hs_supercell
    use, intrinsic :: iso_fortran_env, only: int64
    use types, only: dp
    use HS_io, only: HSBlockSparse
    implicit none
    private

    public :: hs_supercell_build, hs_supercell_index, bs_hk_dense

contains

    pure integer function hs_supercell_index(nat, n, b, t) result(idx)
        integer, intent(in) :: nat, n(3), b, t(3)
        idx = b + nat * (t(1) + n(1) * (t(2) + n(2) * t(3)))
    end function hs_supercell_index

    !> Monta a supercélula de `prim` (H ou S). keep(I) = .false. remove o átomo I
    !> (índice antes da remoção, ver hs_supercell_index); newidx devolve o
    !> índice final (0 se removido).
    subroutine hs_supercell_build(prim, n, periodic, keep, sup, newidx, ierr, errmsg)
        type(HSBlockSparse), intent(in)  :: prim
        integer,             intent(in)  :: n(3)
        logical,             intent(in)  :: periodic(3)
        logical,             intent(in)  :: keep(:)
        type(HSBlockSparse), intent(out) :: sup
        integer, allocatable, intent(out) :: newidx(:)
        integer,             intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg

        integer :: nat, ntot, nkeep, i, ib, t1, t2, t3, t(3), u(3), tp(3), rp(3), a, ia, ja, nb, ir, nr
        integer :: rmin(3), rmax(3), span(3)
        integer, allocatable :: rindex(:), bidx(:, :)
        integer(int64), allocatable :: boff(:), src(:)
        integer(int64) :: tot, sz
        logical :: ok

        ierr = 0
        errmsg = ''
        nat = prim%nAtoms
        if (any(n < 1)) then
            ierr = 1; errmsg = 'hs_supercell_build: multiplicities must be >= 1'; return
        end if
        ntot = nat * n(1) * n(2) * n(3)
        if (size(keep) /= ntot) then
            ierr = 1; errmsg = 'hs_supercell_build: keep has the wrong size'; return
        end if
        allocate(newidx(ntot))
        nkeep = 0
        do i = 1, ntot
            if (keep(i)) then
                nkeep = nkeep + 1
                newidx(i) = nkeep
            else
                newidx(i) = 0
            end if
        end do
        if (nkeep == 0) then
            ierr = 1; errmsg = 'hs_supercell_build: every atom was removed'; return
        end if

        ! faixa dos novos vetores R'
        do a = 1, 3
            rmin(a) = floor(real(minval(prim%IndexList(a, :)), dp) / real(n(a), dp))
            rmax(a) = floor(real(maxval(prim%IndexList(a, :)) + n(a) - 1, dp) / real(n(a), dp))
            if (.not. periodic(a)) then
                rmin(a) = 0; rmax(a) = 0
            end if
        end do
        span = rmax - rmin + 1
        allocate(rindex(span(1) * span(2) * span(3)))
        rindex = 0

        ! contagem e lista de blocos
        nb = 0
        do ib = 1, int(prim%nBlocks)
            nb = nb + n(1) * n(2) * n(3)
        end do
        allocate(bidx(3, nb), src(nb))
        nb = 0
        nr = 0
        do t3 = 0, n(3) - 1
            do t2 = 0, n(2) - 1
                do t1 = 0, n(1) - 1
                    t = [t1, t2, t3]
                    do ib = 1, int(prim%nBlocks)
                        ia = hs_supercell_index(nat, n, prim%BlockIndex(1, ib), t)
                        if (newidx(ia) == 0) cycle
                        u = t + prim%IndexList(:, prim%BlockIndex(3, ib))
                        ok = .true.
                        do a = 1, 3
                            tp(a) = modulo(u(a), n(a))
                            rp(a) = (u(a) - tp(a)) / n(a)
                            if (.not. periodic(a) .and. rp(a) /= 0) ok = .false.
                        end do
                        if (.not. ok) cycle
                        ja = hs_supercell_index(nat, n, prim%BlockIndex(2, ib), tp)
                        if (newidx(ja) == 0) cycle
                        ir = 1 + (rp(1) - rmin(1)) + span(1) * ((rp(2) - rmin(2)) + span(2) * (rp(3) - rmin(3)))
                        if (rindex(ir) == 0) then
                            nr = nr + 1
                            rindex(ir) = nr
                        end if
                        nb = nb + 1
                        bidx(:, nb) = [newidx(ia), newidx(ja), rindex(ir)]
                        src(nb) = int(ib, int64)
                    end do
                end do
            end do
        end do

        sup%Title = 'supercell of: ' // trim(prim%Title)
        sup%Neighbors = nr
        sup%nAtoms = nkeep
        allocate(sup%Degenerecency(nr), sup%IndexList(3, nr), sup%Orbs(nkeep))
        sup%Degenerecency = 1
        do a = 1, size(rindex)
            if (rindex(a) == 0) cycle
            ir = a - 1
            sup%IndexList(1, rindex(a)) = rmin(1) + mod(ir, span(1))
            sup%IndexList(2, rindex(a)) = rmin(2) + mod(ir / span(1), span(2))
            sup%IndexList(3, rindex(a)) = rmin(3) + ir / (span(1) * span(2))
        end do
        do i = 1, ntot
            if (newidx(i) > 0) sup%Orbs(newidx(i)) = prim%Orbs(1 + mod(i - 1, nat))
        end do
        sup%NOrbs = sum(sup%Orbs)
        sup%nBlocks = int(nb, int64)
        allocate(boff(nb + 1))
        boff(1) = 0
        do i = 1, nb
            sz = prim%BlockOffset(src(i) + 1) - prim%BlockOffset(src(i))
            boff(i + 1) = boff(i) + sz
        end do
        tot = boff(nb + 1)
        allocate(sup%BlockIndex(3, nb), sup%BlockOffset(nb + 1), sup%BlockData(tot))
        sup%BlockIndex = bidx(:, 1:nb)
        sup%BlockOffset = boff
        do i = 1, nb
            sup%BlockData(boff(i) + 1:boff(i + 1)) = &
                prim%BlockData(prim%BlockOffset(src(i)) + 1:prim%BlockOffset(src(i) + 1))
        end do
        sup%has_meta = prim%has_meta
        sup%fermi_ev = prim%fermi_ev
        sup%fermi_shift_applied = prim%fermi_shift_applied
        sup%energy_reference = prim%energy_reference
    end subroutine hs_supercell_build

    !> H(k) denso, gauge da rede, k em frações da recíproca do próprio arquivo.
    subroutine bs_hk_dense(bs, kfrac, hk)
        type(HSBlockSparse), intent(in)  :: bs
        real(dp),            intent(in)  :: kfrac(3)
        complex(dp), allocatable, intent(out) :: hk(:, :)
        integer, allocatable :: aoff(:)
        integer :: i, j, ib, ni, nj, a
        complex(dp) :: ph
        real(dp), parameter :: TWO_PI = 6.283185307179586_dp
        allocate(aoff(bs%nAtoms + 1))
        aoff(1) = 0
        do a = 1, bs%nAtoms
            aoff(a + 1) = aoff(a) + bs%Orbs(a)
        end do
        allocate(hk(aoff(bs%nAtoms + 1), aoff(bs%nAtoms + 1)))
        hk = (0.0_dp, 0.0_dp)
        do ib = 1, int(bs%nBlocks)
            i = bs%BlockIndex(1, ib); j = bs%BlockIndex(2, ib)
            ni = bs%Orbs(i); nj = bs%Orbs(j)
            ph = exp(cmplx(0.0_dp, TWO_PI * dot_product(kfrac, real(bs%IndexList(:, bs%BlockIndex(3, ib)), dp)), dp))
            hk(aoff(i) + 1:aoff(i) + ni, aoff(j) + 1:aoff(j) + nj) = hk(aoff(i) + 1:aoff(i) + ni, aoff(j) + 1:aoff(j) + nj) &
                + ph * reshape(bs%BlockData(bs%BlockOffset(ib) + 1:bs%BlockOffset(ib + 1)), [ni, nj])
        end do
    end subroutine bs_hk_dense

end module hs_supercell
