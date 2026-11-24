! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Componente G_par = 0 (média no plano) da eletrostática do ciclo SCC pela densidade completa (A19).
!>
!> Malha 1D na normal (altura h = r . nhat, vácuo em +h): nós y_k = y0 + (k-1) dy, funções chapéu g_k.
!> Partição por linhas (a de Mulliken das cargas q_i do SCC): a densidade de cada átomo central i é a soma
!> sobre os pares (i, j, d) da linha de i (d = r_j + L - r_i; j central, da costura ou, na referência, ausente):
!>   rho_i(r) = sum_{j, ab} D^{ij}_ab phi_ia(r) phi_jb(r),   int rho_i = q_i,
!> e a referência é a linha completa do mesmo átomo no volume (todos os pares até rmax, inclusive com os
!> vizinhos que faltam na superfície), transladada: int (rho_i - rho_i^ref) = dq_i (neutralidade do A5).
!> Para cada par o tensor
!>   W^k_ab = int phi_ia(r) phi_jb(r) g_k(h) d^3r            (sum_k W^k = S_ab da grade)
!> sai de perfis por fatia da grade cúbica das PAO (pao_planar_profile, um por vetor d canônico, como o
!> cache de dipolos), com a altura de cada fatia repartida linearmente entre os dois nós vizinhos.
!> Carga por nó (elétrons por célula 2D): c_k = s sum_q sum_ab D^q_ab W^{q,k}_ab (pares do sistema 1..npair,
!> pares da referência npair+1..npair+nref, com a D do volume).
!> Multipolos do SCC no mesmo esquema: c^mult_k = sum_i [dq_i <rho_s>_k(h_i) + dp_i.nhat <rho_s'>_k(h_i)],
!> <rho_s>(t) = 2 pi int_|t|^rc rho_s(r) r dr (média planar da forma rho_s), <rho_s>' = -2 pi t rho_s(|t|).
!> Potencial (energia do elétron, eV) com a referência G = 0 do Ewald do A5 (V = 0 e campo nulo no fundo):
!>   V_m = -(2 pi ke / A) sum_k (|y_m - y_k| + (y_m - y_k)) c_k   (sym: sem o termo linear).
!> Hamiltoniano (derivada de sum_k dV_k c_k na D simétrica): dH^q_ab = wh_q sum_k dV_k W^{q,k}_ab, wh = 1 para os
!> pares centrais (as duas ordens são linhas de átomos centrais) e 1/2 para os pares com a costura (a linha do
!> átomo da costura não entra no modelo; o bloco do plano recebe a metade de <ia|dV|jb>).
!> Rotinas sem stop; erros por ierr/errmsg. Ordem de soma fixa (resultado independe do número de threads).
module scc_planar
    use types, only: dp
    use pao_basis, only: pao_species_t, pao_planar_profile
    implicit none
    private

    public :: planar_t, planar_setup, planar_charge, planar_dh, planar_mult, planar_pot, planar_read_ref

    real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp
    real(dp), parameter :: HARTREE_EV = 27.211386245988_dp
    real(dp), parameter :: BOHR_ANG = 0.529177210903_dp
    real(dp), parameter :: KE_EVANG = HARTREE_EV * BOHR_ANG
    real(dp), parameter :: DTOL = 1.0e-6_dp           !< Bohr: vetores iguais

    type :: planar_t
        integer :: ny = 0, nn = 0, nprof = 0, npair = 0, nref = 0, no = 0, ncent = 0
        real(dp) :: y0 = 0.0_dp, dy = 0.1_dp, area = 1.0_dp
        real(dp) :: hb = 0.1_dp                           !< passo da grade das PAO (Bohr)
        real(dp) :: sgn = 1.0_dp                          !< nhat(ax) = sgn (normal alinhada ao eixo ax)
        real(dp), allocatable :: y(:)                     !< (ny) Å
        real(dp), allocatable :: prof(:, :, :, :)         !< (no, no, 2 nn + 1, nprof)
        real(dp), allocatable :: pdv(:, :)                !< (3, nprof) d canônico (Bohr)
        integer, allocatable  :: pidx(:)                  !< (npair) perfil do par
        logical, allocatable  :: ptr(:)                   !< (npair) .true.: W = perfil transposto (d não canônico)
        real(dp), allocatable :: href(:), wh(:)           !< (npair + nref) altura do centro da grade (Å); peso em H
        real(dp), allocatable :: fq(:, :), fp(:, :)       !< (ny, ncent) forma rho_s e rho_s' (por Bohr) nos nós
    end type planar_t

contains

    pure logical function canonical(d) result(yes)
        real(dp), intent(in) :: d(3)
        integer :: c
        yes = .true.
        do c = 1, 3
            if (abs(d(c)) > DTOL) then
                yes = d(c) > 0.0_dp
                return
            end if
        end do
    end function canonical

    !> Monta perfis, pares e malha. cent_h(n): altura dos átomos centrais (Å); pair_ic/pair_vec(Å)/pair_kind
    !> da lista de &RHO seguida (nref > 0) dos pares da referência (kind ignorado);
    !> ax: eixo cartesiano da normal, sgn = nhat(ax) = +-1; h_bohr: passo da grade das PAO;
    !> dy (Å), vac (Å acima do átomo mais alto); area (Å^2); rho(0:nr), dr: forma rho_s tabelada (Bohr).
    subroutine planar_setup(pl, sp, h_bohr, ax, sgn, cent_h, pair_ic, pair_vec, pair_kind, nref, dy, vac, area, &
                            rho, dr, ierr, errmsg)
        type(planar_t), intent(out) :: pl
        type(pao_species_t), intent(in) :: sp
        real(dp), intent(in) :: h_bohr, sgn, cent_h(:), pair_vec(:, :), dy, vac, area, rho(0:), dr
        integer, intent(in) :: ax, pair_ic(:), pair_kind(:), nref
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        real(dp), allocatable :: nd(:, :), rbar(:)
        real(dp) :: dc(3), hang, lo, hi, t, rc, f
        integer :: q, k, m, nr, i, nall
        ierr = 0
        errmsg = ''
        if (ax < 1 .or. ax > 3 .or. abs(abs(sgn) - 1.0_dp) > 0.0_dp) then
            ierr = 2; errmsg = 'planar: normal must be a Cartesian axis'; return
        end if
        if (dy <= 0.0_dp .or. h_bohr <= 0.0_dp .or. area <= 0.0_dp) then
            ierr = 2; errmsg = 'planar: needs dy > 0, grid step > 0, area > 0'; return
        end if
        pl%no = sp%norb
        nall = size(pair_ic)
        pl%nref = nref
        pl%npair = nall - nref
        pl%ncent = size(cent_h)
        pl%hb = h_bohr
        pl%sgn = sgn
        pl%dy = dy
        pl%area = area
        rc = min(sp%rcut, sp%rv(sp%nmesh))
        pl%nn = floor(rc / h_bohr)
        hang = h_bohr * BOHR_ANG
        ! perfis por d canônico
        allocate(nd(3, nall), pl%pidx(nall), pl%ptr(nall), pl%href(nall), pl%wh(nall))
        m = 0
        do q = 1, nall
            if (pair_ic(q) < 1 .or. pair_ic(q) > pl%ncent) then
                ierr = 2; errmsg = 'planar: a &RHO pair is not mapped to a central atom'; return
            end if
            dc = pair_vec(:, q) / BOHR_ANG
            pl%ptr(q) = .not. canonical(dc)
            if (pl%ptr(q)) dc = -dc
            do k = 1, m
                if (maxval(abs(nd(:, k) - dc)) <= DTOL) exit
            end do
            if (k > m) then
                m = m + 1
                nd(:, m) = dc
                k = m
            end if
            pl%pidx(q) = k
            pl%href(q) = cent_h(pair_ic(q))
            if (pl%ptr(q)) pl%href(q) = pl%href(q) + sgn * pair_vec(ax, q)
            pl%wh(q) = 1.0_dp
            if (q <= pl%npair) pl%wh(q) = merge(0.5_dp, 1.0_dp, pair_kind(q) == 4 .or. pair_kind(q) == 5)
        end do
        pl%nprof = m
        pl%pdv = nd(:, 1:m)
        allocate(pl%prof(pl%no, pl%no, 2 * pl%nn + 1, m))
        !$omp parallel do default(shared) private(k) schedule(dynamic)
        do k = 1, m
            call pao_planar_profile(sp, pl%pdv(:, k), h_bohr, ax, pl%prof(:, :, :, k))
        end do
        !$omp end parallel do
        ! malha
        lo = huge(1.0_dp)
        hi = -huge(1.0_dp)
        do q = 1, nall
            lo = min(lo, pl%href(q) - real(pl%nn, dp) * hang)
            hi = max(hi, pl%href(q) + real(pl%nn, dp) * hang)
        end do
        hi = max(hi, maxval(cent_h) + vac)
        pl%y0 = dy * real(floor(lo / dy) - 2, dp)
        pl%ny = ceiling((hi - pl%y0) / dy) + 3
        allocate(pl%y(pl%ny))
        do k = 1, pl%ny
            pl%y(k) = pl%y0 + real(k - 1, dp) * dy
        end do
        ! formas planares de rho_s e rho_s' (trapézio em t, passo dr; t de -rc a rc)
        nr = size(rho) - 1
        allocate(rbar(0:nr))
        rbar(nr) = 0.0_dp
        do k = nr - 1, 0, -1
            rbar(k) = rbar(k + 1) + 0.5_dp * dr * 2.0_dp * PI * (rho(k) * real(k, dp) * dr + &
                                                                rho(k + 1) * real(k + 1, dp) * dr)
        end do
        allocate(pl%fq(pl%ny, pl%ncent), pl%fp(pl%ny, pl%ncent))
        pl%fq = 0.0_dp
        pl%fp = 0.0_dp
        do i = 1, pl%ncent
            do k = -nr, nr
                t = real(k, dp) * dr
                f = dr
                if (abs(k) == nr) f = 0.5_dp * dr
                call bin(pl, cent_h(i) + t * BOHR_ANG, rbar(abs(k)) * f, pl%fq(:, i), ierr)
                if (ierr /= 0) exit
                call bin(pl, cent_h(i) + t * BOHR_ANG, -2.0_dp * PI * t * rho(abs(k)) * f, pl%fp(:, i), ierr)
                if (ierr /= 0) exit
            end do
            if (ierr /= 0) then
                errmsg = 'planar: rho_s support outside the grid'; return
            end if
        end do
    end subroutine planar_setup

    !> Soma v na malha em y (Å) pela função chapéu.
    subroutine bin(pl, y, v, c, ierr)
        type(planar_t), intent(in) :: pl
        real(dp), intent(in) :: y, v
        real(dp), intent(inout) :: c(:)
        integer, intent(out) :: ierr
        real(dp) :: u, w
        integer :: k
        ierr = 0
        u = (y - pl%y0) / pl%dy
        k = floor(u)
        w = u - real(k, dp)
        if (k < 0 .or. k + 2 > pl%ny) then
            ierr = 1; return
        end if
        c(k + 1) = c(k + 1) + (1.0_dp - w) * v
        c(k + 2) = c(k + 2) + w * v
    end subroutine bin

    !> Nó e peso da fatia s (1..2 nn + 1) do par q.
    pure subroutine slab_node(pl, q, s, k, w)
        type(planar_t), intent(in) :: pl
        integer, intent(in) :: q, s
        integer, intent(out) :: k
        real(dp), intent(out) :: w
        real(dp) :: u
        u = (pl%href(q) + pl%sgn * real(s - pl%nn - 1, dp) * pl%hb * BOHR_ANG - pl%y0) / pl%dy
        k = floor(u)
        w = u - real(k, dp)
        k = k + 1
    end subroutine slab_node

    !> c_k = sum_q sum_ab D^q_ab W^{q,k}_ab (sem fator de spin; linhas: átomo central i). ref = .false.: pares do
    !> sistema, d(no, no, npair); ref = .true.: pares da referência, d(no, no, nref).
    subroutine planar_charge(pl, d, c, ref)
        type(planar_t), intent(in) :: pl
        real(dp), intent(in) :: d(:, :, :)
        real(dp), intent(out) :: c(:)
        logical, intent(in), optional :: ref
        real(dp), allocatable :: dt(:, :)
        real(dp) :: v, w
        integer :: q, s, k, ip, q0, q1
        q0 = 0
        q1 = pl%npair
        if (present(ref)) then
            if (ref) then
                q0 = pl%npair
                q1 = pl%nref
            end if
        end if
        c = 0.0_dp
        allocate(dt(pl%no, pl%no))
        do q = q0 + 1, q0 + q1
            ip = pl%pidx(q)
            if (pl%ptr(q)) then
                dt = transpose(d(1:pl%no, 1:pl%no, q - q0))
            else
                dt = d(1:pl%no, 1:pl%no, q - q0)
            end if
            do s = 1, 2 * pl%nn + 1
                v = sum(dt * pl%prof(:, :, s, ip))
                if (v == 0.0_dp) cycle
                call slab_node(pl, q, s, k, w)
                c(k) = c(k) + (1.0_dp - w) * v
                c(k + 1) = c(k + 1) + w * v
            end do
        end do
    end subroutine planar_charge

    !> dh^q_ab = wh_q sum_k v_k W^{q,k}_ab (eV) para cada par do sistema.
    subroutine planar_dh(pl, v, dh)
        type(planar_t), intent(in) :: pl
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: dh(:, :, :)
        real(dp), allocatable :: a(:, :)
        real(dp) :: w, vs
        integer :: q, s, k, ip
        allocate(a(pl%no, pl%no))
        dh = 0.0_dp
        do q = 1, pl%npair
            ip = pl%pidx(q)
            a = 0.0_dp
            do s = 1, 2 * pl%nn + 1
                call slab_node(pl, q, s, k, w)
                vs = pl%wh(q) * ((1.0_dp - w) * v(k) + w * v(k + 1))
                a = a + vs * pl%prof(:, :, s, ip)
            end do
            if (pl%ptr(q)) then
                dh(1:pl%no, 1:pl%no, q) = transpose(a)
            else
                dh(1:pl%no, 1:pl%no, q) = a
            end if
        end do
    end subroutine planar_dh

    !> Carga dos multipolos do SCC nos nós: dq (e), dpn = dp . nhat (e Bohr).
    subroutine planar_mult(pl, dq, dpn, c)
        type(planar_t), intent(in) :: pl
        real(dp), intent(in) :: dq(:), dpn(:)
        real(dp), intent(out) :: c(:)
        integer :: i
        c = 0.0_dp
        do i = 1, pl%ncent
            c = c + dq(i) * pl%fq(:, i) + dpn(i) * pl%fp(:, i)
        end do
    end subroutine planar_mult

    !> V_m (eV) das cargas c_k (e por célula): referência no fundo; sym: só o termo |y - y'|.
    subroutine planar_pot(pl, c, v, sym)
        type(planar_t), intent(in) :: pl
        real(dp), intent(in) :: c(:)
        real(dp), intent(out) :: v(:)
        logical, intent(in), optional :: sym
        real(dp) :: dd, f
        integer :: m, k
        logical :: sy
        sy = .false.
        if (present(sym)) sy = sym
        f = -2.0_dp * PI * KE_EVANG / pl%area
        do m = 1, pl%ny
            v(m) = 0.0_dp
            do k = 1, pl%ny
                dd = pl%y(m) - pl%y(k)
                if (sy) then
                    v(m) = v(m) + abs(dd) * c(k)
                else
                    v(m) = v(m) + (abs(dd) + dd) * c(k)
                end if
            end do
            v(m) = f * v(m)
        end do
    end subroutine planar_pot

    !> Linhas de referência a partir de um PREFIX.rho_pairs.dat do volume (janela 1): para cada átomo central i
    !> (átomo de base b = mod(cent_atom(i) - 1, nbasis) + 1), todos os pares (i0, j, d) de uma cópia central i0 de
    !> b no volume (a primeira do arquivo), com a D do volume: ref_ic, ref_vec (Å), dref(no, no, :).
    !> nrow(b): pares por linha.
    subroutine planar_read_ref(file, nbasis, no, cent_atom, ref_ic, ref_vec, dref, nrow, ierr, errmsg)
        character(len=*), intent(in) :: file
        integer, intent(in) :: nbasis, no, cent_atom(:)
        integer, allocatable, intent(out) :: ref_ic(:), nrow(:)
        real(dp), allocatable, intent(out) :: ref_vec(:, :), dref(:, :, :)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=4096) :: ln
        integer :: u, ios, nrp, q, kind, ia, icol, ja, ni, nj, a, np, b, i, m
        integer, allocatable :: kia(:), kcol(:), kni(:), knj(:), i0(:), c0(:)
        real(dp), allocatable :: kd(:, :), blk(:, :, :)
        real(dp) :: dv(3), dn
        logical :: inwin
        ierr = 0
        errmsg = ''
        open(newunit=u, file=file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = 'planar: cannot open planar_ref_file ' // trim(file); return
        end if
        nrp = 0
        do
            read(u, '(A)', iostat=ios) ln
            if (ios /= 0) exit
            if (ln(1:5) == 'pair ') nrp = nrp + 1
        end do
        if (nrp == 0) then
            close(u); ierr = 2; errmsg = 'planar: no pairs in planar_ref_file'; return
        end if
        allocate(kia(nrp), kcol(nrp), kd(3, nrp), kni(nrp), knj(nrp), blk(no, no, nrp))
        blk = 0.0_dp
        rewind(u)
        inwin = .false.
        np = 0
        do
            read(u, '(A)', iostat=ios) ln
            if (ios /= 0) exit
            if (ln(1:5) == 'pair ') then
                read(ln(6:), *, iostat=ios) q, kind, ia, icol, ja, ni, nj, dv, dn
                if (ios /= 0 .or. q < 1 .or. q > nrp .or. ni > no .or. nj > no) then
                    close(u); ierr = 2; errmsg = 'planar: bad pair line in planar_ref_file'; return
                end if
                kia(q) = ia
                kcol(q) = icol
                kd(:, q) = dv
                kni(q) = ni
                knj(q) = nj
            else if (ln(1:6) == 'window') then
                if (inwin) exit
                inwin = .true.
            else if (ln(1:5) == 'block' .and. inwin) then
                read(ln(6:), *, iostat=ios) q
                if (ios /= 0 .or. q < 1 .or. q > nrp) then
                    close(u); ierr = 2; errmsg = 'planar: bad block line in planar_ref_file'; return
                end if
                do a = 1, kni(q)
                    read(u, *, iostat=ios) blk(a, 1:knj(q), q)
                    if (ios /= 0) then
                        close(u); ierr = 2; errmsg = 'planar: truncated block in planar_ref_file'; return
                    end if
                end do
                np = np + 1
            end if
        end do
        close(u)
        if (np /= nrp) then
            ierr = 2; errmsg = 'planar: planar_ref_file has fewer blocks than pairs (window 1)'; return
        end if
        ! cópia de cada átomo de base: o primeiro átomo central do arquivo
        allocate(i0(nbasis), c0(nbasis), nrow(nbasis))
        i0 = 0
        c0 = 0
        nrow = 0
        do q = 1, nrp
            b = modulo(kia(q) - 1, nbasis) + 1
            if (i0(b) == 0) then
                i0(b) = kia(q); c0(b) = kcol(q)
            end if
            if (kia(q) == i0(b) .and. kcol(q) == c0(b)) nrow(b) = nrow(b) + 1
        end do
        if (any(nrow == 0)) then
            ierr = 2; errmsg = 'planar: planar_ref_file lacks some basis atoms'; return
        end if
        m = 0
        do i = 1, size(cent_atom)
            m = m + nrow(modulo(cent_atom(i) - 1, nbasis) + 1)
        end do
        allocate(ref_ic(m), ref_vec(3, m), dref(no, no, m))
        m = 0
        do i = 1, size(cent_atom)
            b = modulo(cent_atom(i) - 1, nbasis) + 1
            do q = 1, nrp
                if (kia(q) /= i0(b) .or. kcol(q) /= c0(b)) cycle
                m = m + 1
                ref_ic(m) = i
                ref_vec(:, m) = kd(:, q)
                dref(:, :, m) = blk(:, :, q)
            end do
        end do
    end subroutine planar_read_ref

end module scc_planar
