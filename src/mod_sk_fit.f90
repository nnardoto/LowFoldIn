! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Ajuste de parâmetros de Slater–Koster ao acoplamento entre camadas
!> principais de um bulk, e regra de combinação para pares de espécies que
!> nunca foram vizinhas (o contato A/B).
!>
!> Alvo (acoplamento "ancorado" em E0):
!>
!>   T(R_par) = H_{0,1}(R_par) - E0 S_{0,1}(R_par)
!>
!> Trocar (H01, S01) por (T, 0) reproduz K01(z) = z S01 - H01 exatamente em
!> z = E0 e erra por (z - E0) S01 longe dela. Com spin não colinear o alvo
!> espacial é Re[(T_up,up + T_dn,dn)/2]; a parte fora da diagonal em spin e a
!> parte imaginária são descartadas e relatadas.
!>
!> Modelo: para cada par de átomos (i na PL 0, j na PL 1 transladada por
!> R_par) a menos de cutoff, e camadas (a em i, b em j) da lista de ajuste,
!>
!>   T_ab ~ sum_c V_c(entry) exp(-beta (d - d0(entry))) E^c_ab(d),
!>
!> linear nos V_c para beta fixo (mínimos quadrados, LAPACK dgels). beta é
!> global, escolhido por busca em grade (grossa e depois fina) ou fixo; d0 de
!> uma linha é a menor distância entre os pares que ela descreve.
!>
!> Combinação (contato entre espécies X de A e Y de B, nunca vizinhas):
!>   V_XY(a, b) = sgn * sqrt(|V_XX(a, b) V_YY(a, b)|) no d0 médio, com cada
!>   valor levado ao d0 médio pela sua própria lei exponencial; sinais
!>   opostos dão 0 e são relatados. beta = média.
module sk_fit
    use types, only: dp
    use sk_table, only: sk_block, sk_channel_block, species_basis_t, basis_shell_index, l_letter
    use sk_contact, only: sk_entry_t, sk_model_t, sk_geom_t, sk_pairs_t, make_sk_geom, build_sk_pairs, &
                          sk_contact_at_k, sk_find_entry, parse_shell, SK_NAME_LEN
    use layer_blocks, only: pl_blocks_t, pl_blocks_at_k
    use spectral_projection, only: generalized_eigensolve
    implicit none
    private

    public :: fit_shell_set_t, sk_fit_report_t, parse_fit_shells, sk_fit_material, sk_band_check, sk_combine

    real(dp), parameter :: TWO_PI = 6.283185307179586476925286766559_dp

    interface
        subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: trans
            integer,  intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dgels
    end interface

    !> Camadas (l, zeta) a ajustar, por espécie.
    type :: fit_shell_set_t
        character(len=SK_NAME_LEN) :: name = ''
        integer :: n = 0
        integer :: l(16) = 0, zeta(16) = 0
    end type fit_shell_set_t

    type :: sk_fit_report_t
        integer  :: npairs = 0, nrows = 0, nparams = 0
        real(dp) :: beta = 0.0_dp
        real(dp) :: rel_res_fitted = 0.0_dp     !< ||A x - t|| / ||t|| nos sub-blocos ajustados
        real(dp) :: rel_err_full = 0.0_dp       !< ||T_SK - T|| / ||T|| no bloco H01 inteiro (todas as camadas, todos R)
        real(dp) :: dropped_spin = 0.0_dp       !< ||parte descartada (Im, fora da diagonal de spin)|| / ||T||
        real(dp) :: anchor = 0.0_dp
    end type sk_fit_report_t

contains

    !> 'Te s1 p1 d1' -> camadas de Te. Entradas vazias são ignoradas.
    subroutine parse_fit_shells(texts, sets, ierr, errmsg)
        character(len=*), intent(in) :: texts(:)
        type(fit_shell_set_t), allocatable, intent(out) :: sets(:)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=64) :: tok(17)
        integer :: k, n, ntok, ios, t, l, z, e1
        ierr = 0
        errmsg = ''
        n = count([(len_trim(texts(k)) > 0, k = 1, size(texts))])
        allocate(sets(n))
        n = 0
        do k = 1, size(texts)
            if (len_trim(texts(k)) == 0) cycle
            n = n + 1
            tok = ''
            ntok = 0
            do t = 1, 17
                read(texts(k), *, iostat=ios) tok(1:t)
                if (ios /= 0) exit
                ntok = t
            end do
            if (ntok < 2) then
                ierr = 1; errmsg = "fit_shells entry must be 'Species s1 p1 ...', got '" // trim(texts(k)) // "'"
                return
            end if
            sets(n)%name = tok(1)
            do t = 2, ntok
                call parse_shell(tok(t), l, z, e1)
                if (e1 /= 0) then
                    ierr = 1; errmsg = "invalid shell '" // trim(tok(t)) // "' in fit_shells"; return
                end if
                sets(n)%n = sets(n)%n + 1
                sets(n)%l(sets(n)%n) = l
                sets(n)%zeta(sets(n)%n) = z
            end do
        end do
    end subroutine parse_fit_shells

    !> Ajusta as linhas SK ao acoplamento entre PLs do material.
    !> model: bases e cutoff (as linhas são substituídas pelo resultado).
    !> sets: camadas a ajustar por espécie (vazio: primeira zeta de s, p, d).
    !> beta_fixed >= 0: beta fixo; senão busca em [beta_min, beta_max].
    subroutine sk_fit_material(pl, labels, model, sets, anchor, beta_fixed, beta_min, beta_max, rep, ierr, errmsg)
        type(pl_blocks_t),     intent(in)    :: pl
        character(len=*),      intent(in)    :: labels(:)
        type(sk_model_t),      intent(inout) :: model
        type(fit_shell_set_t), intent(in)    :: sets(:)
        real(dp),              intent(in)    :: anchor, beta_fixed, beta_min, beta_max
        type(sk_fit_report_t), intent(out)   :: rep
        integer,               intent(out)   :: ierr
        character(len=:), allocatable, intent(out) :: errmsg

        type(sk_geom_t) :: g
        type(sk_entry_t), allocatable :: keys(:)
        integer, allocatable :: nch(:), col0(:)
        ! dados de cada sub-bloco: par, camada i, camada j, entrada, orientação
        integer, allocatable :: sb_i(:), sb_j(:), sb_r(:), sb_li(:), sb_lj(:), sb_oi(:), sb_oj(:), sb_key(:)
        real(dp), allocatable :: sb_d(:, :), t(:), a0(:, :), x(:), dmin(:)
        integer :: i, j, ir, si, sj, li, lj, zi, zj, k, nsb, nrow, ncol, r, a, b, c, q, nb
        real(dp) :: d(3), dist, a1(3), a2(3), best_res, res, beta, bstep, tnorm
        real(dp), allocatable :: tsp(:, :, :)
        complex(dp) :: tu

        ierr = 0
        errmsg = ''
        if (.not. pl%has_geometry) then
            ierr = 1; errmsg = 'sk_fit: principal layers have no geometry'; return
        end if
        call make_sk_geom(model, pl%norb, pl%tau, pl%orb_atom, pl%orb_cell, labels, pl%orbs_per_atom, g, ierr, errmsg)
        if (ierr /= 0) return
        a1 = pl%lattice(pl%inplane(1), :)
        a2 = pl%lattice(pl%inplane(2), :)
        rep%anchor = anchor

        ! alvo espacial T(R) (norb_espacial x norb_espacial por R) e norma descartada
        call spatial_target(pl, g, anchor, tsp, tnorm, rep%dropped_spin)

        ! --- sub-blocos e chaves ---
        allocate(keys(0), sb_i(0), sb_j(0), sb_r(0), sb_li(0), sb_lj(0), sb_oi(0), sb_oj(0), sb_key(0))
        allocate(sb_d(3, 0), dmin(0))
        nsb = 0
        do ir = 1, pl%nr
            do i = 1, g%nat
                do j = 1, g%nat
                    d = pl%pl_shift + g%pos(:, j) + real(pl%rpar(1, ir), dp) * a1 + &
                        real(pl%rpar(2, ir), dp) * a2 - g%pos(:, i)
                    dist = norm2(d)
                    if (dist > model%cutoff .or. dist < 1.0e-6_dp) cycle
                    associate (bi => model%basis(g%ibasis(i)), bj => model%basis(g%ibasis(j)))
                        do si = 1, nset_of(g%label(i))
                            call set_shell(g%label(i), si, li, zi)
                            if (basis_shell_index(bi, li, zi) == 0) cycle
                            do sj = 1, nset_of(g%label(j))
                                call set_shell(g%label(j), sj, lj, zj)
                                if (basis_shell_index(bj, lj, zj) == 0) cycle
                                k = key_index(g%label(i), li, zi, g%label(j), lj, zj, dist)
                                nsb = nsb + 1
                                sb_i = [sb_i, i]; sb_j = [sb_j, j]; sb_r = [sb_r, ir]
                                sb_li = [sb_li, li]; sb_lj = [sb_lj, lj]
                                sb_oi = [sb_oi, bi%shell(basis_shell_index(bi, li, zi))%offset]
                                sb_oj = [sb_oj, bj%shell(basis_shell_index(bj, lj, zj))%offset]
                                sb_key = [sb_key, k]
                                sb_d = reshape([sb_d, d], [3, nsb])
                            end do
                        end do
                    end associate
                end do
            end do
        end do
        if (nsb == 0 .or. size(keys) == 0) then
            ierr = 2
            errmsg = 'sk_fit: no atom pairs between neighbouring principal layers within the cutoff ' // &
                     '(raise cutoff_ang or check fit_shells/basis)'
            return
        end if
        do k = 1, size(keys)
            keys(k)%d0 = dmin(k)
        end do
        allocate(nch(size(keys)), col0(size(keys)))
        ncol = 0
        do k = 1, size(keys)
            nch(k) = min(keys(k)%la, keys(k)%lb) + 1
            col0(k) = ncol
            ncol = ncol + nch(k)
        end do
        nrow = 0
        do q = 1, nsb
            nrow = nrow + (2 * sb_li(q) + 1) * (2 * sb_lj(q) + 1)
        end do
        if (nrow < ncol) then
            ierr = 2; errmsg = 'sk_fit: fewer data than parameters (raise cutoff_ang or reduce fit_shells)'; return
        end if
        rep%nrows = nrow
        rep%nparams = ncol
        allocate(t(nrow), a0(nrow, ncol), x(ncol))

        ! --- vetor alvo e matriz sem o fator de distância ---
        r = 0
        a0 = 0.0_dp
        do q = 1, nsb
            li = sb_li(q); lj = sb_lj(q)
            k = sb_key(q)
            do b = 1, 2 * lj + 1
                do a = 1, 2 * li + 1
                    r = r + 1
                    t(r) = tsp(g%off_sp(sb_i(q)) + sb_oi(q) + a, g%off_sp(sb_j(q)) + sb_oj(q) + b, sb_r(q))
                    do c = 0, nch(k) - 1
                        a0(r, col0(k) + c + 1) = channel_elem(li, lj, sb_d(:, q), c, a, b)
                    end do
                end do
            end do
        end do

        ! --- beta ---
        if (beta_fixed >= 0.0_dp) then
            beta = beta_fixed
        else
            best_res = huge(1.0_dp)
            beta = beta_min
            nb = 40
            bstep = (beta_max - beta_min) / real(nb, dp)
            do c = 0, nb
                res = solve_at(beta_min + bstep * real(c, dp))
                if (res < best_res) then
                    best_res = res; beta = beta_min + bstep * real(c, dp)
                end if
            end do
            bstep = bstep / 20.0_dp
            res = beta
            do c = -20, 20
                if (res + bstep * real(c, dp) < 0.0_dp) cycle
                if (solve_at(res + bstep * real(c, dp)) < best_res) then
                    best_res = solve_at(res + bstep * real(c, dp)); beta = res + bstep * real(c, dp)
                end if
            end do
        end if
        rep%rel_res_fitted = solve_at(beta)
        if (ierr /= 0) return
        rep%beta = beta
        rep%npairs = nsb
        do k = 1, size(keys)
            keys(k)%beta = beta
            keys(k)%v = 0.0_dp
            keys(k)%v(1:nch(k)) = x(col0(k) + 1:col0(k) + nch(k))
        end do
        model%nentry = size(keys)
        call move_alloc(keys, model%e)
        model%lambda = 1.0_dp
        model%scale = 1.0_dp
        model%beta_override = -1.0_dp

        ! --- erro no bloco inteiro (todas as camadas, todos os R) ---
        rep%rel_err_full = full_block_error()

    contains

        integer function nset_of(name)
            character(len=*), intent(in) :: name
            integer :: s
            nset_of = 0
            do s = 1, size(sets)
                if (trim(sets(s)%name) == trim(name)) nset_of = sets(s)%n
            end do
            if (size(sets) == 0 .or. all([(trim(sets(s)%name) /= trim(name), s = 1, size(sets))])) then
                ! padrão: primeira zeta de s, p, d presentes na base
                nset_of = 0
                do s = 0, 2
                    if (basis_shell_index(model%basis(basis_of(name)), s, 1) > 0) nset_of = nset_of + 1
                end do
            end if
        end function nset_of

        subroutine set_shell(name, idx, l, z)
            character(len=*), intent(in) :: name
            integer, intent(in) :: idx
            integer, intent(out) :: l, z
            integer :: s, cnt
            do s = 1, size(sets)
                if (trim(sets(s)%name) == trim(name)) then
                    l = sets(s)%l(idx); z = sets(s)%zeta(idx); return
                end if
            end do
            cnt = 0
            do s = 0, 2
                if (basis_shell_index(model%basis(basis_of(name)), s, 1) > 0) then
                    cnt = cnt + 1
                    if (cnt == idx) then
                        l = s; z = 1; return
                    end if
                end if
            end do
            l = 0; z = 1
        end subroutine set_shell

        integer function basis_of(name)
            character(len=*), intent(in) :: name
            integer :: s
            basis_of = 1
            do s = 1, model%nbasis
                if (trim(model%basis(s)%name) == trim(name)) basis_of = s
            end do
        end function basis_of

        !> Índice da chave (par não ordenado de camadas); cria se não existe.
        integer function key_index(sa, la, za, sb, lb, zb, dist)
            character(len=*), intent(in) :: sa, sb
            integer, intent(in) :: la, za, lb, zb
            real(dp), intent(in) :: dist
            type(sk_entry_t) :: e
            integer :: s
            do s = 1, size(keys)
                if ((trim(keys(s)%spec_a) == trim(sa) .and. keys(s)%la == la .and. keys(s)%za == za .and. &
                     trim(keys(s)%spec_b) == trim(sb) .and. keys(s)%lb == lb .and. keys(s)%zb == zb) .or. &
                    (trim(keys(s)%spec_a) == trim(sb) .and. keys(s)%la == lb .and. keys(s)%za == zb .and. &
                     trim(keys(s)%spec_b) == trim(sa) .and. keys(s)%lb == la .and. keys(s)%zb == za)) then
                    key_index = s
                    dmin(s) = min(dmin(s), dist)
                    return
                end if
            end do
            ! forma canônica: menor (l, zeta) primeiro; empate pelo nome
            if (la < lb .or. (la == lb .and. (za < zb .or. (za == zb .and. sa <= sb)))) then
                e = sk_entry_t(sa, sb, la, za, lb, zb, 0.0_dp, 0.0_dp, [0.0_dp, 0.0_dp, 0.0_dp])
            else
                e = sk_entry_t(sb, sa, lb, zb, la, za, 0.0_dp, 0.0_dp, [0.0_dp, 0.0_dp, 0.0_dp])
            end if
            keys = [keys, e]
            dmin = [dmin, dist]
            key_index = size(keys)
        end function key_index

        !> Resíduo relativo para beta dado; preenche x.
        real(dp) function solve_at(bt)
            real(dp), intent(in) :: bt
            real(dp), allocatable :: am(:, :), rhs(:), work(:)
            real(dp) :: f, wq(1)
            integer :: rr, info, lwork, kk, cc, qq, nr_q
            allocate(am(nrow, ncol), rhs(nrow))
            rr = 0
            do qq = 1, nsb
                kk = sb_key(qq)
                f = exp(-bt * (norm2(sb_d(:, qq)) - keys(kk)%d0))
                nr_q = (2 * sb_li(qq) + 1) * (2 * sb_lj(qq) + 1)
                am(rr + 1:rr + nr_q, :) = f * a0(rr + 1:rr + nr_q, :)
                rr = rr + nr_q
            end do
            rhs = t
            call dgels('N', nrow, ncol, 1, am, nrow, rhs, nrow, wq, -1, info)
            lwork = max(1, int(wq(1)))
            allocate(work(lwork))
            call dgels('N', nrow, ncol, 1, am, nrow, rhs, nrow, work, lwork, info)
            if (info /= 0) then
                ierr = 3; errmsg = 'sk_fit: least-squares solve failed (rank-deficient design?)'
                solve_at = huge(1.0_dp); return
            end if
            x = rhs(1:ncol)
            ! resíduo: recomputa A x - t
            rr = 0
            f = 0.0_dp
            do qq = 1, nsb
                kk = sb_key(qq)
                nr_q = (2 * sb_li(qq) + 1) * (2 * sb_lj(qq) + 1)
                do cc = 1, nr_q
                    f = f + (exp(-bt * (norm2(sb_d(:, qq)) - keys(kk)%d0)) * &
                             dot_product(a0(rr + cc, :), x) - t(rr + cc))**2
                end do
                rr = rr + nr_q
            end do
            solve_at = sqrt(f) / max(norm2(t), tiny(1.0_dp))
        end function solve_at

        real(dp) function full_block_error()
            type(sk_pairs_t) :: pairs
            complex(dp), allocatable :: hk(:, :)
            real(dp) :: num, den
            integer :: ii, jj, rr
            character(len=:), allocatable :: m2
            integer :: e2
            ! T_SK(R) no espaço real: soma das contribuições por R (via pares)
            call build_sk_pairs(model, g, g, pl%pl_shift, a1, a2, pairs, e2, m2)
            if (e2 /= 0) then
                full_block_error = huge(1.0_dp); return
            end if
            num = 0.0_dp
            den = 0.0_dp
            allocate(hk(pl%norb, pl%norb))
            do rr = 1, pl%nr
                call pairs_at_r(pairs, pl%rpar(:, rr), hk)
                do jj = 1, pl%norb
                    do ii = 1, pl%norb
                        tu = pl%h(ii, jj, rr, 1)
                        if (pl%has_overlap) tu = tu - anchor * pl%s(ii, jj, rr, 1)
                        num = num + abs(hk(ii, jj) - tu)**2
                        den = den + abs(tu)**2
                    end do
                end do
            end do
            full_block_error = sqrt(num / max(den, tiny(1.0_dp)))
        end function full_block_error

    end subroutine sk_fit_material

    !> Contribuição dos pares com R_par = rr (bloco real no espaço real).
    subroutine pairs_at_r(pairs, rr, h)
        type(sk_pairs_t), intent(in) :: pairs
        integer,          intent(in) :: rr(2)
        complex(dp),      intent(out) :: h(:, :)
        integer :: p, a, b, q
        h = (0.0_dp, 0.0_dp)
        do p = 1, pairs%npair
            if (any(pairs%n12(:, p) /= rr)) cycle
            q = pairs%boff(p)
            do b = 1, pairs%nj(p)
                do a = 1, pairs%ni(p)
                    h(pairs%ioff(p) + a, pairs%joff(p) + b) = h(pairs%ioff(p) + a, pairs%joff(p) + b) + &
                        pairs%blk(q + (b - 1) * pairs%ni(p) + a - 1)
                    if (pairs%spinor(p)) h(pairs%ioff(p) + pairs%ni(p) + a, pairs%joff(p) + pairs%nj(p) + b) = &
                        h(pairs%ioff(p) + pairs%ni(p) + a, pairs%joff(p) + pairs%nj(p) + b) + &
                        pairs%blk(q + (b - 1) * pairs%ni(p) + a - 1)
                end do
            end do
        end do
    end subroutine pairs_at_r

    pure real(dp) function channel_elem(li, lj, d, c, a, b)
        integer,  intent(in) :: li, lj, c, a, b
        real(dp), intent(in) :: d(3)
        real(dp) :: e(2 * li + 1, 2 * lj + 1)
        e = sk_channel_block(li, lj, d, c)
        channel_elem = e(a, b)
    end function channel_elem

    !> Alvo espacial Re[(T_uu + T_dd)/2] (ou Re T), indexado pelos orbitais
    !> espaciais da PL; norma relativa descartada (Im e fora da diagonal de spin).
    subroutine spatial_target(pl, g, anchor, tsp, tnorm, dropped)
        type(pl_blocks_t), intent(in)    :: pl
        type(sk_geom_t),   intent(inout) :: g
        real(dp),          intent(in)    :: anchor
        real(dp), allocatable, intent(out) :: tsp(:, :, :)
        real(dp),          intent(out)   :: tnorm, dropped
        integer :: nsp_tot, i, j, ir, a, b, ri, rj
        complex(dp) :: tij, tuu, tdd
        real(dp) :: keep, all2
        integer, allocatable :: off_sp(:)

        allocate(off_sp(g%nat))
        nsp_tot = 0
        do i = 1, g%nat
            off_sp(i) = nsp_tot
            nsp_tot = nsp_tot + g%nsp(i)
        end do
        call move_alloc(off_sp, g%off_sp)
        allocate(tsp(nsp_tot, nsp_tot, pl%nr))
        tsp = 0.0_dp
        all2 = 0.0_dp
        keep = 0.0_dp
        do ir = 1, pl%nr
            do j = 1, g%nat
                do i = 1, g%nat
                    do b = 1, g%nsp(j)
                        do a = 1, g%nsp(i)
                            ri = g%off(i) + a
                            rj = g%off(j) + b
                            tuu = t_of(ri, rj, ir)
                            if (g%spinor(i)) then
                                tdd = t_of(ri + g%nsp(i), rj + g%nsp(j), ir)
                                tsp(g%off_sp(i) + a, g%off_sp(j) + b, ir) = 0.5_dp * real(tuu + tdd, dp)
                                keep = keep + 2.0_dp * (0.5_dp * real(tuu + tdd, dp))**2
                            else
                                tsp(g%off_sp(i) + a, g%off_sp(j) + b, ir) = real(tuu, dp)
                                keep = keep + real(tuu, dp)**2
                            end if
                        end do
                    end do
                end do
            end do
            do j = 1, pl%norb
                do i = 1, pl%norb
                    tij = t_of(i, j, ir)
                    all2 = all2 + abs(tij)**2
                end do
            end do
        end do
        tnorm = sqrt(all2)
        dropped = sqrt(max(all2 - keep, 0.0_dp) / max(all2, tiny(1.0_dp)))
    contains
        complex(dp) function t_of(ii, jj, rr)
            integer, intent(in) :: ii, jj, rr
            t_of = pl%h(ii, jj, rr, 1)
            if (pl%has_overlap) t_of = t_of - anchor * pl%s(ii, jj, rr, 1)
        end function t_of
    end subroutine spatial_target

    !> Teste de bandas do bulk ao longo de k_s em k_par fixo. Compara os
    !> autovalores exatos na janela [anchor - window, anchor + window] com os de
    !> (a) acoplamento ancorado (H01 - E0 S01, S01 = 0) e (b) SK (S01 = 0).
    !> Para cada autovalor exato, distância ao autovalor mais próximo do modelo.
    subroutine sk_band_check(pl, labels, model, kpar, nks, anchor, window, max_anchor, rms_anchor, &
                             max_sk, rms_sk, nval, ierr, errmsg)
        type(pl_blocks_t), intent(in)  :: pl
        character(len=*),  intent(in)  :: labels(:)
        type(sk_model_t),  intent(in)  :: model
        real(dp),          intent(in)  :: kpar(2), anchor, window
        integer,           intent(in)  :: nks
        real(dp),          intent(out) :: max_anchor, rms_anchor, max_sk, rms_sk
        integer,           intent(out) :: nval
        integer,           intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        type(sk_geom_t) :: g
        type(sk_pairs_t) :: pairs
        complex(dp), allocatable :: h(:, :, :), s(:, :, :), hsk(:, :), hk(:, :), sk(:, :), c(:, :)
        real(dp), allocatable :: e0(:), e1(:), e2(:)
        real(dp) :: ks, sa, ss
        complex(dp) :: ph
        integer :: n, iks, i

        ierr = 0
        errmsg = ''
        n = pl%norb
        call make_sk_geom(model, n, pl%tau, pl%orb_atom, pl%orb_cell, labels, pl%orbs_per_atom, g, ierr, errmsg)
        if (ierr /= 0) return
        call build_sk_pairs(model, g, g, pl%pl_shift, pl%lattice(pl%inplane(1), :), pl%lattice(pl%inplane(2), :), &
                            pairs, ierr, errmsg)
        if (ierr /= 0) return
        allocate(h(n, n, -1:1), s(n, n, -1:1), hsk(n, n), hk(n, n), sk(n, n), c(n, n), e0(n), e1(n), e2(n))
        call pl_blocks_at_k(pl, kpar, h, s)
        call sk_contact_at_k(pairs, kpar, hsk)
        max_anchor = 0.0_dp; max_sk = 0.0_dp; sa = 0.0_dp; ss = 0.0_dp; nval = 0
        do iks = 0, nks - 1
            ks = 0.5_dp * real(iks, dp) / real(max(nks - 1, 1), dp)
            ph = exp(cmplx(0.0_dp, TWO_PI * ks, dp))
            ! exato
            hk = h(:, :, 0) + ph * h(:, :, 1) + conjg(ph) * h(:, :, -1)
            sk = s(:, :, 0) + ph * s(:, :, 1) + conjg(ph) * s(:, :, -1)
            call generalized_eigensolve(n, herm(hk), herm(sk), e0, c, ierr)
            if (ierr /= 0) goto 9
            ! ancorado
            hk = h(:, :, 0) + ph * (h(:, :, 1) - anchor * s(:, :, 1)) + conjg(ph) * (h(:, :, -1) - anchor * s(:, :, -1))
            call generalized_eigensolve(n, herm(hk), herm(s(:, :, 0)), e1, c, ierr)
            if (ierr /= 0) goto 9
            ! SK
            hk = h(:, :, 0) + ph * hsk + conjg(ph) * conjg(transpose(hsk))
            call generalized_eigensolve(n, herm(hk), herm(s(:, :, 0)), e2, c, ierr)
            if (ierr /= 0) goto 9
            do i = 1, n
                if (abs(e0(i) - anchor) > window) cycle
                nval = nval + 1
                max_anchor = max(max_anchor, minval(abs(e1 - e0(i))))
                max_sk = max(max_sk, minval(abs(e2 - e0(i))))
                sa = sa + minval(abs(e1 - e0(i)))**2
                ss = ss + minval(abs(e2 - e0(i)))**2
            end do
        end do
        rms_anchor = sqrt(sa / real(max(nval, 1), dp))
        rms_sk = sqrt(ss / real(max(nval, 1), dp))
        return
9       errmsg = 'sk_band_check: generalized eigenproblem failed (S not positive definite?)'
    contains
        function herm(a) result(o)
            complex(dp), intent(in) :: a(:, :)
            complex(dp) :: o(size(a, 1), size(a, 2))
            o = 0.5_dp * (a + conjg(transpose(a)))
        end function herm
    end subroutine sk_band_check

    !> Regra de combinação: para cada espécie X de xa e Y de yb, e cada par
    !> de camadas presente nas linhas homonucleares X-X (de ea) e Y-Y (de eb),
    !> cria a linha X-Y. nwarn conta canais zerados por sinais opostos.
    subroutine sk_combine(ea, eb, xa, yb, out, nwarn)
        type(sk_entry_t), intent(in) :: ea(:), eb(:)
        character(len=*), intent(in) :: xa(:), yb(:)
        type(sk_entry_t), allocatable, intent(out) :: out(:)
        integer, intent(out) :: nwarn
        integer :: ix, iy, p, q, c
        real(dp) :: d0m, vx, vy
        type(sk_entry_t) :: e
        allocate(out(0))
        nwarn = 0
        do ix = 1, size(xa)
            if (len_trim(xa(ix)) == 0) cycle
            do iy = 1, size(yb)
                if (len_trim(yb(iy)) == 0) cycle
                do p = 1, size(ea)
                    if (trim(ea(p)%spec_a) /= trim(xa(ix)) .or. trim(ea(p)%spec_b) /= trim(xa(ix))) cycle
                    do q = 1, size(eb)
                        if (trim(eb(q)%spec_a) /= trim(yb(iy)) .or. trim(eb(q)%spec_b) /= trim(yb(iy))) cycle
                        if (ea(p)%la /= eb(q)%la .or. ea(p)%za /= eb(q)%za .or. &
                            ea(p)%lb /= eb(q)%lb .or. ea(p)%zb /= eb(q)%zb) cycle
                        d0m = 0.5_dp * (ea(p)%d0 + eb(q)%d0)
                        e = ea(p)
                        e%d0 = d0m
                        e%beta = 0.5_dp * (ea(p)%beta + eb(q)%beta)
                        do c = 1, 3
                            vx = ea(p)%v(c) * exp(-ea(p)%beta * (d0m - ea(p)%d0))
                            vy = eb(q)%v(c) * exp(-eb(q)%beta * (d0m - eb(q)%d0))
                            if (vx * vy > 0.0_dp) then
                                e%v(c) = sign(sqrt(vx * vy), vx)
                            else
                                if (vx /= 0.0_dp .or. vy /= 0.0_dp) nwarn = nwarn + 1
                                e%v(c) = 0.0_dp
                            end if
                        end do
                        ! X-Y: a camada la de X e lb de Y, e (se la /= lb ou za /= zb) também a inversa
                        e%spec_a = xa(ix); e%spec_b = yb(iy)
                        out = [out, e]
                        if (e%la /= e%lb .or. e%za /= e%zb) then
                            e%spec_a = yb(iy); e%spec_b = xa(ix)
                            out = [out, e]
                        end if
                    end do
                end do
            end do
        end do
    end subroutine sk_combine

end module sk_fit
