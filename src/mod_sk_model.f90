! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Modelo Slater–Koster ortogonal (sem SOC) com on-site por l e tabela de
!> camadas de vizinhos, e montagem de H(R) bloco-esparso com S = I.
!>
!> Parâmetros (mesma convenção do sk_gen.py da trincheira de Te):
!>   p = [E_l para cada l da base, na ordem s, p, d] ++ [V_q para cada canal q]
!>   canal q = (la, lb, camada, c): integral de dois centros entre os orbitais
!>   la do átomo i e lb do átomo j na camada radial `camada` (0-based, índice
!>   em `shells`), c = 0 (sigma), 1 (pi), 2 (delta).
!> Bloco de uma ligação d = r_j + R - r_i (|d| < rc, camada = a de distância
!> mais próxima, a no máximo shell_tol):
!>   B[la, lb] += V_q E_{la lb}(d),   e, se la /= lb,   B[lb, la] += V_q E_{lb la}(d),
!> com E da tabela sk_table (menor l primeiro). On-site: diag(E_l).
!> Orbitais por átomo na ordem do OpenMX: [s]; px py pz; [d3z2-r2 dx2-y2 dxy dxz dyz].
!>
!> H(k) = sum_R H(R) e^{+i 2 pi k.n_R} (gauge da rede, como no resto do
!> LowFoldIn), H(R)_{ij} = <i,0|H|j,R>.
!>
!> Leitura dos parâmetros: JSON do ajuste em Python (te_skg_{pd,spd}.json),
!> chaves basis, channels, shells, rc, p, nso_per_atom (as demais são
!> ignoradas). O leitor só entende o subconjunto de JSON desses arquivos
!> (objeto plano; valores string, número ou listas de números aninhadas).
!>
!> Rotinas sem stop; erros por ierr/errmsg.
module sk_model
    use, intrinsic :: iso_fortran_env, only: int64
    use types, only: dp
    use sk_table, only: sk_channel_block
    use HS_io, only: HSBlockSparse
    implicit none
    private

    public :: sk_model_t, read_sk_model_json, write_sk_model_json, sk_onsite_block, sk_bond_block
    public :: sk_build_blocksparse, sk_identity_blocksparse, sk_hk_dense

    type :: sk_model_t
        character(len=16) :: basis = ''
        integer :: nl = 0
        integer :: ls(3) = -1                  !< valores de l presentes, crescentes
        integer :: loff(0:2) = -1              !< deslocamento (0-based) de cada l no átomo; -1 ausente
        integer :: norb = 0
        integer :: nch = 0
        integer, allocatable :: ch(:, :)       !< (4, nch): la, lb, camada (0-based), c
        real(dp), allocatable :: shells(:)     !< distâncias das camadas (Å)
        real(dp) :: rc = 0.0_dp                !< raio de corte (Å)
        real(dp) :: shell_tol = 0.03_dp        !< tolerância de casamento com as camadas (Å)
        real(dp), allocatable :: p(:)          !< nl + nch parâmetros (eV)
        integer :: nso_per_atom = 0            !< bandas ocupadas por átomo (informativo)
    end type sk_model_t

contains

    ! ==================================================================
    ! Leitura / escrita
    ! ==================================================================

    !> Lê o modelo de um JSON do ajuste (te_skg_*.json).
    subroutine read_sk_model_json(filename, m, ierr, errmsg)
        character(len=*),  intent(in)  :: filename
        type(sk_model_t),  intent(out) :: m
        integer,           intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=:), allocatable :: txt, sval
        real(dp), allocatable :: v(:)
        integer :: u, ios, nbytes, q

        ierr = 0
        errmsg = ''
        open(newunit=u, file=filename, status='old', action='read', access='stream', &
             form='unformatted', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = "cannot open SK model file '" // trim(filename) // "'"; return
        end if
        inquire(unit=u, size=nbytes)
        allocate(character(len=max(nbytes, 1)) :: txt)
        txt(:) = ' '
        if (nbytes > 0) read(u, iostat=ios) txt
        close(u)
        if (ios /= 0) then
            ierr = 1; errmsg = "cannot read SK model file '" // trim(filename) // "'"; return
        end if

        call json_string(txt, 'basis', sval, ierr)
        if (ierr /= 0) then
            errmsg = "SK model: key 'basis' missing or not a string"; return
        end if
        call set_basis(m, sval, ierr, errmsg)
        if (ierr /= 0) return

        call json_numbers(txt, 'channels', v, ierr)
        if (ierr /= 0 .or. mod(size(v), 4) /= 0) then
            ierr = 1; errmsg = "SK model: 'channels' must be a list of [la, lb, shell, c]"; return
        end if
        m%nch = size(v) / 4
        allocate(m%ch(4, m%nch))
        m%ch = reshape(nint(v), [4, m%nch])

        call json_numbers(txt, 'shells', v, ierr)
        if (ierr /= 0 .or. size(v) < 1) then
            ierr = 1; errmsg = "SK model: 'shells' missing"; return
        end if
        m%shells = v

        call json_numbers(txt, 'rc', v, ierr)
        if (ierr /= 0 .or. size(v) /= 1) then
            ierr = 1; errmsg = "SK model: 'rc' missing"; return
        end if
        m%rc = v(1)

        call json_numbers(txt, 'p', v, ierr)
        if (ierr /= 0) then
            ierr = 1; errmsg = "SK model: 'p' missing"; return
        end if
        if (size(v) /= m%nl + m%nch) then
            ierr = 1; errmsg = "SK model: size of 'p' differs from (number of l) + (number of channels)"; return
        end if
        m%p = v

        call json_numbers(txt, 'nso_per_atom', v, ierr)
        if (ierr == 0 .and. size(v) == 1) m%nso_per_atom = nint(v(1))
        ierr = 0

        do q = 1, m%nch
            if (any(m%ch(1:2, q) < 0) .or. any(m%ch(1:2, q) > 2) .or. m%ch(4, q) < 0 &
                .or. m%ch(4, q) > min(m%ch(1, q), m%ch(2, q)) .or. m%ch(3, q) < 0 &
                .or. m%ch(3, q) >= size(m%shells)) then
                ierr = 1; errmsg = 'SK model: invalid channel entry'; return
            end if
            if (m%loff(m%ch(1, q)) < 0 .or. m%loff(m%ch(2, q)) < 0) then
                ierr = 1; errmsg = 'SK model: channel uses an l absent from the basis'; return
            end if
        end do
    end subroutine read_sk_model_json

    !> Grava o modelo no mesmo JSON lido por read_sk_model_json.
    subroutine write_sk_model_json(filename, m, ierr, errmsg, note)
        character(len=*),  intent(in)  :: filename
        type(sk_model_t),  intent(in)  :: m
        integer,           intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=*), optional, intent(in) :: note
        integer :: u, ios, q
        character(len=*), parameter :: LSYM = 'spd', CSYM = 'spd'

        ierr = 0
        errmsg = ''
        open(newunit=u, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = "cannot write SK model file '" // trim(filename) // "'"; return
        end if
        write(u, '(a)') '{'
        write(u, '(a)') ' "basis": "' // trim(m%basis) // '",'
        write(u, '(a)', advance='no') ' "channels": ['
        do q = 1, m%nch
            write(u, '(a,3(i0,", "),i0,a)', advance='no') '[', m%ch(:, q), ']'
            if (q < m%nch) write(u, '(a)', advance='no') ', '
        end do
        write(u, '(a)') '],'
        write(u, '(a)', advance='no') ' "names": ['
        do q = 1, m%nl
            write(u, '(a)', advance='no') '"E_' // LSYM(m%ls(q) + 1:m%ls(q) + 1) // '"'
            if (q < m%nl + m%nch) write(u, '(a)', advance='no') ', '
        end do
        do q = 1, m%nch
            write(u, '(a,a,a,a,i0,a)', advance='no') '"', LSYM(m%ch(1, q) + 1:m%ch(1, q) + 1) // &
                LSYM(m%ch(2, q) + 1:m%ch(2, q) + 1), CSYM(m%ch(4, q) + 1:m%ch(4, q) + 1), '_', m%ch(3, q) + 1, '"'
            if (q < m%nch) write(u, '(a)', advance='no') ', '
        end do
        write(u, '(a)') '],'
        write(u, '(a)', advance='no') ' "shells": ['
        do q = 1, size(m%shells)
            write(u, '(es24.16e3)', advance='no') m%shells(q)
            if (q < size(m%shells)) write(u, '(a)', advance='no') ', '
        end do
        write(u, '(a)') '],'
        write(u, '(a,es24.16e3,a)') ' "rc": ', m%rc, ','
        write(u, '(a)', advance='no') ' "p": ['
        do q = 1, size(m%p)
            write(u, '(es24.16e3)', advance='no') m%p(q)
            if (q < size(m%p)) write(u, '(a)', advance='no') ', '
        end do
        write(u, '(a)') '],'
        write(u, '(a,i0,a)', advance='no') ' "nso_per_atom": ', m%nso_per_atom
        if (present(note)) then
            write(u, '(a)') ','
            write(u, '(a)') ' "fonte": "' // trim(note) // '"'
        else
            write(u, '(a)') ''
        end if
        write(u, '(a)') '}'
        close(u)
    end subroutine write_sk_model_json

    subroutine set_basis(m, basis, ierr, errmsg)
        type(sk_model_t), intent(inout) :: m
        character(len=*), intent(in)    :: basis
        integer,          intent(out)   :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: k, l
        ierr = 0
        errmsg = ''
        m%basis = basis
        m%nl = 0
        m%loff = -1
        m%norb = 0
        do l = 0, 2
            if (index(basis, 'spd'(l + 1:l + 1)) > 0) then
                m%nl = m%nl + 1
                m%ls(m%nl) = l
                m%loff(l) = m%norb
                m%norb = m%norb + 2 * l + 1
            end if
        end do
        do k = 1, len_trim(basis)
            if (index('spd', basis(k:k)) == 0) then
                ierr = 1; errmsg = "SK model: basis must be a combination of 's', 'p', 'd', got '" // &
                    trim(basis) // "'"; return
            end if
        end do
        if (m%nl == 0) then
            ierr = 1; errmsg = 'SK model: empty basis'
        end if
    end subroutine set_basis

    !> Posição logo após ':' da chave "key" (0 se ausente).
    pure integer function json_value_pos(txt, key) result(p)
        character(len=*), intent(in) :: txt, key
        integer :: i, j, n
        character(len=:), allocatable :: pat
        pat = '"' // key // '"'
        n = len(txt)
        p = 0
        i = 1
        do
            j = index(txt(i:), pat)
            if (j == 0) return
            j = i + j - 1 + len(pat)
            do while (j <= n)
                if (txt(j:j) /= ' ' .and. txt(j:j) /= achar(9) .and. txt(j:j) /= achar(10) &
                    .and. txt(j:j) /= achar(13)) exit
                j = j + 1
            end do
            if (j <= n) then
                if (txt(j:j) == ':') then
                    p = j + 1
                    return
                end if
            end if
            i = j
            if (i > n) return
        end do
    end function json_value_pos

    subroutine json_string(txt, key, val, ierr)
        character(len=*), intent(in) :: txt, key
        character(len=:), allocatable, intent(out) :: val
        integer, intent(out) :: ierr
        integer :: p, q
        val = ''
        ierr = 1
        p = json_value_pos(txt, key)
        if (p == 0) return
        q = index(txt(p:), '"')
        if (q == 0) return
        p = p + q
        if (verify(txt(p - q:p - 2), ' ' // achar(9) // achar(10) // achar(13)) /= 0) return
        q = index(txt(p:), '"')
        if (q == 0) return
        val = txt(p:p + q - 2)
        ierr = 0
    end subroutine json_string

    !> Todos os números do valor da chave (número ou listas aninhadas), achatados.
    subroutine json_numbers(txt, key, v, ierr)
        character(len=*), intent(in) :: txt, key
        real(dp), allocatable, intent(out) :: v(:)
        integer, intent(out) :: ierr
        integer :: p, n, depth, i, j, ios, cnt
        real(dp) :: x
        real(dp), allocatable :: buf(:)
        character(len=*), parameter :: NUMCH = '0123456789+-.eE'

        ierr = 1
        allocate(v(0))
        p = json_value_pos(txt, key)
        if (p == 0) return
        n = len(txt)
        allocate(buf(64))
        cnt = 0
        depth = 0
        i = p
        do while (i <= n)
            select case (txt(i:i))
            case ('[')
                depth = depth + 1
                i = i + 1
            case (']')
                depth = depth - 1
                i = i + 1
                if (depth <= 0) exit
            case (',')
                if (depth == 0) exit
                i = i + 1
            case ('}')
                exit
            case (' ', achar(9), achar(10), achar(13))
                i = i + 1
            case default
                if (index(NUMCH, txt(i:i)) == 0) return
                j = i
                do while (j <= n)
                    if (index(NUMCH, txt(j:j)) == 0) exit
                    j = j + 1
                end do
                read(txt(i:j - 1), *, iostat=ios) x
                if (ios /= 0) return
                cnt = cnt + 1
                if (cnt > size(buf)) buf = [buf, [(0.0_dp, ios = 1, size(buf))]]
                buf(cnt) = x
                i = j
                if (depth == 0) exit
            end select
        end do
        v = buf(1:cnt)
        ierr = 0
    end subroutine json_numbers

    ! ==================================================================
    ! Blocos
    ! ==================================================================

    !> Bloco on-site diag(E_l) (norb x norb).
    pure function sk_onsite_block(m) result(b)
        type(sk_model_t), intent(in) :: m
        real(dp) :: b(m%norb, m%norb)
        integer :: k, l, i
        b = 0.0_dp
        do k = 1, m%nl
            l = m%ls(k)
            do i = m%loff(l) + 1, m%loff(l) + 2 * l + 1
                b(i, i) = m%p(k)
            end do
        end do
    end function sk_onsite_block

    !> Bloco de uma ligação d (Å), de i para j. ierr = 1 se |d| não casa com
    !> nenhuma camada (a tolerância é shell_tol).
    subroutine sk_bond_block(m, d, b, ierr)
        type(sk_model_t), intent(in)  :: m
        real(dp),         intent(in)  :: d(3)
        real(dp),         intent(out) :: b(m%norb, m%norb)
        integer,          intent(out) :: ierr
        real(dp) :: r, v
        integer :: s, q, la, lb, c, oa, ob
        ierr = 0
        b = 0.0_dp
        r = norm2(d)
        s = minloc(abs(m%shells - r), dim=1)
        if (abs(m%shells(s) - r) > m%shell_tol) then
            ierr = 1; return
        end if
        do q = 1, m%nch
            if (m%ch(3, q) /= s - 1) cycle
            la = m%ch(1, q); lb = m%ch(2, q); c = m%ch(4, q)
            v = m%p(m%nl + q)
            oa = m%loff(la); ob = m%loff(lb)
            b(oa + 1:oa + 2 * la + 1, ob + 1:ob + 2 * lb + 1) = b(oa + 1:oa + 2 * la + 1, ob + 1:ob + 2 * lb + 1) &
                + v * sk_channel_block(la, lb, d, c)
            if (la /= lb) then
                b(ob + 1:ob + 2 * lb + 1, oa + 1:oa + 2 * la + 1) = b(ob + 1:ob + 2 * lb + 1, oa + 1:oa + 2 * la + 1) &
                    + v * sk_channel_block(lb, la, d, c)
            end if
        end do
    end subroutine sk_bond_block

    ! ==================================================================
    ! Montagem bloco-esparsa
    ! ==================================================================

    !> H(R) do modelo para os átomos `cart` (3, nat; Å) na rede `lattice`
    !> (linhas a_i, Å), periódica nos eixos com periodic(a) = .true. Todos os
    !> pares 1e-3 < |d| < rc entram; os blocos on-site ficam em R = 0, que é o
    !> primeiro vetor de IndexList. ierr = 1 se alguma distância não casa com
    !> as camadas do modelo.
    subroutine sk_build_blocksparse(m, lattice, cart, periodic, bs, ierr, errmsg)
        type(sk_model_t),    intent(in)  :: m
        real(dp),            intent(in)  :: lattice(3, 3), cart(:, :)
        logical,             intent(in)  :: periodic(3)
        type(HSBlockSparse), intent(out) :: bs
        integer,             intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg

        integer :: nat, nmax(3), a, n1, n2, n3, i, j, nR, ib, no, nb, nbmax, iR
        integer, allocatable :: rlist(:, :), bidx(:, :)
        real(dp) :: d(3), shift(3), vol, h, r
        real(dp), allocatable :: blk(:, :), bdat(:, :, :)
        logical :: any_in_R
        integer, allocatable :: tmpi(:, :)
        real(dp), allocatable :: tmpd(:, :, :)
        character(len=96) :: msg

        ierr = 0
        errmsg = ''
        nat = size(cart, 2)
        no = m%norb
        vol = abs(dot_product(lattice(1, :), cross(lattice(2, :), lattice(3, :))))
        if (vol < 1.0e-10_dp .or. nat < 1) then
            ierr = 1; errmsg = 'sk_build_blocksparse: singular lattice or no atoms'; return
        end if
        do a = 1, 3
            nmax(a) = 0
            if (.not. periodic(a)) cycle
            h = vol / norm2(cross(lattice(mod(a, 3) + 1, :), lattice(mod(a + 1, 3) + 1, :)))
            nmax(a) = ceiling(m%rc / h) + 1
        end do
        allocate(blk(no, no))
        nbmax = 64
        allocate(rlist(3, 0), bidx(3, nbmax), bdat(no, no, nbmax))
        nb = 0
        nR = 0
        ! R = 0 primeiro
        do iR = 0, 1
            do n1 = -nmax(1), nmax(1)
                do n2 = -nmax(2), nmax(2)
                    do n3 = -nmax(3), nmax(3)
                        if ((iR == 0) .neqv. (n1 == 0 .and. n2 == 0 .and. n3 == 0)) cycle
                        shift = n1 * lattice(1, :) + n2 * lattice(2, :) + n3 * lattice(3, :)
                        any_in_R = .false.
                        do i = 1, nat
                            do j = 1, nat
                                d = cart(:, j) + shift - cart(:, i)
                                r = norm2(d)
                                if (r >= m%rc) cycle
                                if (r <= 1.0e-3_dp) then
                                    if (i /= j) then
                                        ierr = 1; errmsg = 'sk_build_blocksparse: overlapping atoms'; return
                                    end if
                                    blk = sk_onsite_block(m)
                                else
                                    call sk_bond_block(m, d, blk, ierr)
                                    if (ierr /= 0) then
                                        write(msg, '(a,f9.4,a,i0,a,i0)') 'distance ', r, ' A between atoms ', i, &
                                            ' and ', j
                                        errmsg = 'sk_build_blocksparse: ' // trim(msg) // ' matches no shell'
                                        return
                                    end if
                                end if
                                if (.not. any_in_R) then
                                    any_in_R = .true.
                                    nR = nR + 1
                                    rlist = reshape([rlist, [n1, n2, n3]], [3, nR])
                                end if
                                nb = nb + 1
                                if (nb > nbmax) then
                                    tmpi = bidx; tmpd = bdat
                                    nbmax = 2 * nbmax
                                    deallocate(bidx, bdat)
                                    allocate(bidx(3, nbmax), bdat(no, no, nbmax))
                                    bidx(:, 1:nb - 1) = tmpi(:, 1:nb - 1)
                                    bdat(:, :, 1:nb - 1) = tmpd(:, :, 1:nb - 1)
                                end if
                                bidx(:, nb) = [i, j, nR]
                                bdat(:, :, nb) = blk
                            end do
                        end do
                    end do
                end do
            end do
        end do

        bs%Title = 'Slater-Koster H (LowFoldIn skgen)'
        bs%Neighbors = nR
        bs%nAtoms = nat
        bs%NOrbs = nat * no
        allocate(bs%Degenerecency(nR))
        bs%Degenerecency = 1
        allocate(bs%Orbs(nat))
        bs%Orbs = no
        bs%IndexList = rlist
        bs%nBlocks = int(nb, int64)
        allocate(bs%BlockIndex(3, nb), bs%BlockOffset(nb + 1), bs%BlockData(int(nb, int64) * no * no))
        bs%BlockOffset(1) = 0
        do ib = 1, nb
            bs%BlockIndex(:, ib) = bidx(:, ib)
            bs%BlockOffset(ib + 1) = bs%BlockOffset(ib) + int(no * no, int64)
            bs%BlockData(bs%BlockOffset(ib) + 1:bs%BlockOffset(ib + 1)) = &
                cmplx(reshape(bdat(:, :, ib), [no * no]), 0.0_dp, dp)
        end do
    end subroutine sk_build_blocksparse

    !> S = I com a mesma lista de vetores R de `h` (blocos só em R = 0, i = j).
    subroutine sk_identity_blocksparse(h, s)
        type(HSBlockSparse), intent(in)  :: h
        type(HSBlockSparse), intent(out) :: s
        integer :: i, k, no, r0, nat
        nat = h%nAtoms
        s%Title = 'Slater-Koster S = I (LowFoldIn skgen)'
        s%Neighbors = h%Neighbors
        s%nAtoms = nat
        s%NOrbs = h%NOrbs
        s%Degenerecency = h%Degenerecency
        s%Orbs = h%Orbs
        s%IndexList = h%IndexList
        r0 = 0
        do k = 1, h%Neighbors
            if (all(h%IndexList(:, k) == 0)) r0 = k
        end do
        s%nBlocks = int(nat, int64)
        allocate(s%BlockIndex(3, nat), s%BlockOffset(nat + 1))
        s%BlockOffset(1) = 0
        do i = 1, nat
            s%BlockIndex(:, i) = [i, i, r0]
            s%BlockOffset(i + 1) = s%BlockOffset(i) + int(h%Orbs(i), int64)**2
        end do
        allocate(s%BlockData(s%BlockOffset(nat + 1)))
        s%BlockData = (0.0_dp, 0.0_dp)
        do i = 1, nat
            no = h%Orbs(i)
            do k = 1, no
                s%BlockData(s%BlockOffset(i) + (k - 1) * no + k) = (1.0_dp, 0.0_dp)
            end do
        end do
    end subroutine sk_identity_blocksparse

    !> H(k) denso (gauge da rede) de um HSBlockSparse, k em frações da recíproca.
    subroutine sk_hk_dense(bs, kfrac, hk)
        type(HSBlockSparse), intent(in)  :: bs
        real(dp),            intent(in)  :: kfrac(3)
        complex(dp), allocatable, intent(out) :: hk(:, :)
        integer, allocatable :: aoff(:)
        integer :: i, j, ib, ni, nj, a
        complex(dp) :: ph
        real(dp), parameter :: TWO_PI = 6.283185307179586476925286766559_dp
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
    end subroutine sk_hk_dense

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
    end function cross

end module sk_model
