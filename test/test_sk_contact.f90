! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Contato de Slater–Koster e varredura de parâmetros (sk_contact, sk_fit,
!> green_engine com model = 'sk' e &SCAN).
!>
!> Material sintético: duas espécies (Xa s1p1, Xb s1p1d1), rede hexagonal no
!> plano, vetor de empilhamento inclinado, S(R) não ortogonal. O acoplamento
!> entre camadas principais é, por construção, exatamente SK:
!>   H_{0,1}(R) = SK_true(R) + E0 S_{0,1}(R).
!> O gerador abaixo percorre (i, j, R) por conta própria; a tabela SK foi
!> validada contra uma referência independente em test_sk_table.
!>
!> Testes:
!>  1. skfit recupera beta e os V de todas as linhas (espaço s, p, d e
!>     heteronuclear) e o erro do bloco H01 inteiro é ~0; o teste de bandas do
!>     SK coincide com o do acoplamento ancorado.
!>  2. junção A|A, sanduíche A|A|A e filme com contato 'sk' (parâmetros
!>     verdadeiros) = contato 'scaled' com lambda = 1 (o H01 exato), com spin
!>     não colinear: geometria, deslocamentos, sentido B-abaixo-de-A, fases e
!>     spin do contato estão certos;
!>  3. lambda_h = 0 = contato 'none';
!>  4. &SCAN de lambda_h, film_shift_ev (meio B refeito por valor) e
!>     offset_normal = execuções separadas;
!>  5. causalidade com contato e deslocamento genéricos; erros de entrada.
program test_sk_contact
    use types, only: dp
    use testing, only: check, check_close, test_summary
    use HS_io, only: HSDataFormat
    use structure_io, only: structure_t, read_structure
    use layer_blocks, only: pl_blocks_t, build_pl_blocks, attach_pl_geometry
    use sk_table, only: sk_block, species_basis_t, parse_basis
    use sk_contact, only: sk_model_t, sk_entry_t, sk_set_bases, write_sk_file, sk_find_entry, read_sk_file
    use sk_fit, only: fit_shell_set_t, sk_fit_report_t, sk_fit_material, sk_band_check
    use green_engine, only: green_input_t, read_green_input, run_green
    use green_data, only: green_data_t
    implicit none

    real(dp), parameter :: LAT(3, 3) = reshape([3.2_dp, -1.6_dp, 0.35_dp, &
                                                 0.0_dp, 2.7712812921102_dp, 0.20_dp, &
                                                 0.0_dp, 0.0_dp, 4.6_dp], [3, 3])
    real(dp), parameter :: TAU(3, 2) = reshape([0.10_dp, 0.20_dp, 0.30_dp, 1.70_dp, 0.95_dp, 2.10_dp], [3, 2])
    character(len=2), parameter :: LAB(2) = ['Xa', 'Xb']
    character(len=16), parameter :: BASIS_TXT(2) = ['Xa s1p1         ', 'Xb s1p1d1       ']
    real(dp), parameter :: RC = 5.5_dp, BETA_TRUE = 1.1_dp, D0_TRUE = 3.0_dp

    type(sk_model_t) :: truth, intra, ovl
    integer :: ierr, nbad
    character(len=:), allocatable :: errmsg

    call make_models()
    call test_fit()
    call test_engine()
    call cleanup()
    call test_summary('test_sk_contact')

contains

    ! ------------------------------------------------------------------
    ! Modelos SK: verdadeiro (entre PLs), intra-PL e de overlap
    ! ------------------------------------------------------------------
    subroutine make_models()
        call fill(truth, 0.8_dp, BETA_TRUE, 11)
        call fill(intra, 0.6_dp, 0.9_dp, 23)
        call fill(ovl, 0.03_dp, 1.0_dp, 37)
    end subroutine make_models

    subroutine fill(m, amp, beta, seed)
        type(sk_model_t), intent(out) :: m
        real(dp), intent(in) :: amp, beta
        integer, intent(in) :: seed
        integer :: st, k, la, lb, sa, sb
        character(len=2) :: spa, spb
        type(sk_entry_t), allocatable :: e(:)
        type(sk_entry_t) :: x
        st = seed
        allocate(e(0))
        do sa = 1, 2
            do sb = sa, 2
                do la = 0, merge(1, 2, sa == 1)
                    do lb = 0, merge(1, 2, sb == 1)
                        if (sa == sb .and. lb < la) cycle
                        spa = LAB(sa); spb = LAB(sb)
                        x = sk_entry_t(spa, spb, la, 1, lb, 1, D0_TRUE, beta, [0.0_dp, 0.0_dp, 0.0_dp])
                        do k = 1, min(la, lb) + 1
                            x%v(k) = amp * rnd(st)
                        end do
                        e = [e, x]
                    end do
                end do
            end do
        end do
        m%nentry = size(e)
        call move_alloc(e, m%e)
        m%cutoff = RC
        call sk_set_bases(m, BASIS_TXT, ierr, errmsg)
    end subroutine fill

    real(dp) function rnd(st)
        integer, intent(inout) :: st
        st = modulo(1103515245 * st + 12345, 2147483647)
        rnd = 2.0_dp * real(st, dp) / 2147483647.0_dp - 1.0_dp
        if (abs(rnd) < 0.2_dp) rnd = sign(0.2_dp, rnd)
    end function rnd

    !> Bloco espacial (nsp_i x nsp_j) do modelo m entre átomos i, j à distância d.
    function block_of(m, i, j, d) result(b)
        type(sk_model_t), intent(in) :: m
        integer, intent(in) :: i, j
        real(dp), intent(in) :: d(3)
        real(dp), allocatable :: b(:, :)
        type(species_basis_t) :: bi, bj
        integer :: si, sj, k, oi, oj, li, lj
        real(dp) :: v(3)
        call parse_basis(BASIS_TXT(i), bi, ierr, errmsg)
        call parse_basis(BASIS_TXT(j), bj, ierr, errmsg)
        allocate(b(bi%nspatial, bj%nspatial))
        b = 0.0_dp
        do si = 1, bi%nshell
            do sj = 1, bj%nshell
                li = bi%shell(si)%l; lj = bj%shell(sj)%l
                oi = bi%shell(si)%offset; oj = bj%shell(sj)%offset
                k = sk_find_entry(m, LAB(i), li, 1, LAB(j), lj, 1)
                if (k == 0) cycle
                v = m%e(k)%v * exp(-m%e(k)%beta * (norm2(d) - m%e(k)%d0))
                b(oi + 1:oi + 2 * li + 1, oj + 1:oj + 2 * lj + 1) = sk_block(li, lj, d, v)
            end do
        end do
    end function block_of

    !> Escreve PREFIX_H.hsio, PREFIX_S.hsio e PREFIX.nml.
    !> H(R) = intra (n3 = 0) ou truth + anchor S (n3 = +-1); S(R) = I + ovl
    !> (n3 = 0) ou s01 * ovl (n3 = +-1). spinor: V x I_2 + SOC fictício no sítio.
    subroutine make_material(prefix, spinor, s01, anchor)
        character(len=*), intent(in) :: prefix
        logical, intent(in) :: spinor
        real(dp), intent(in) :: s01, anchor
        type(HSDataFormat) :: hh, hs
        integer :: nsp(2), off(2), nsc, ns, n, nr, n1, n2, n3, ir, i, j, a, st, u, ir0
        real(dp) :: d(3), r(3)
        real(dp), allocatable :: bh(:, :), bs(:, :)
        nsp = [4, 9]
        nsc = merge(2, 1, spinor)
        off = [0, nsc * nsp(1)]
        n = nsc * sum(nsp)
        ns = 4
        nr = (2 * ns + 1)**2 * 3
        hh%nAtoms = 2; hh%NOrbs = n; hh%Neighbors = nr; hh%Title = 'sk synthetic H'
        allocate(hh%Orbs(2), hh%IndexList(3, nr), hh%Degenerecency(nr), hh%Data(n, n, nr))
        hh%Orbs = nsc * nsp
        hh%Degenerecency = 1
        hh%Data = (0.0_dp, 0.0_dp)
        ir = 0
        ir0 = 0
        do n3 = -1, 1
            do n2 = -ns, ns
                do n1 = -ns, ns
                    ir = ir + 1
                    hh%IndexList(:, ir) = [n1, n2, n3]
                    if (n1 == 0 .and. n2 == 0 .and. n3 == 0) ir0 = ir
                end do
            end do
        end do
        hs = hh
        hs%Title = 'sk synthetic S'
        do ir = 1, nr
            r = matmul(real(hh%IndexList(:, ir), dp), LAT)
            n3 = hh%IndexList(3, ir)
            do i = 1, 2
                do j = 1, 2
                    d = TAU(:, j) + r - TAU(:, i)
                    if (norm2(d) > RC .or. norm2(d) < 1.0e-8_dp) cycle
                    if (n3 == 0) then
                        bh = block_of(intra, i, j, d)
                        bs = block_of(ovl, i, j, d)
                    else
                        bs = s01 * block_of(ovl, i, j, d)
                        bh = block_of(truth, i, j, d) + anchor * bs
                    end if
                    do a = 0, nsc - 1
                        hh%Data(off(i) + a * nsp(i) + 1:off(i) + (a + 1) * nsp(i), &
                                off(j) + a * nsp(j) + 1:off(j) + (a + 1) * nsp(j), ir) = bh
                        hs%Data(off(i) + a * nsp(i) + 1:off(i) + (a + 1) * nsp(i), &
                                off(j) + a * nsp(j) + 1:off(j) + (a + 1) * nsp(j), ir) = bs
                    end do
                end do
            end do
        end do
        ! no sítio: energias distintas, S = I, SOC fictício (hermitiano) nos espinores
        st = 5
        do i = 1, n
            hh%Data(i, i, ir0) = hh%Data(i, i, ir0) + cmplx(-1.5_dp + 3.0_dp * real(i - 1, dp) / real(n - 1, dp), &
                                                             0.0_dp, dp)
            hs%Data(i, i, ir0) = hs%Data(i, i, ir0) + (1.0_dp, 0.0_dp)
        end do
        if (spinor) then
            do a = 1, 2
                do i = 1, nsp(a)
                    do j = 1, nsp(a)
                        hh%Data(off(a) + i, off(a) + nsp(a) + j, ir0) = cmplx(0.05_dp * rnd(st), 0.05_dp * rnd(st), dp)
                        hh%Data(off(a) + nsp(a) + j, off(a) + i, ir0) = conjg(hh%Data(off(a) + i, off(a) + nsp(a) + j, ir0))
                    end do
                end do
            end do
        end if
        call hh%WriteFile(prefix // '_H.hsio')
        call hs%WriteFile(prefix // '_S.hsio')
        open(newunit=u, file=prefix // '.nml', status='replace', action='write')
        write(u, '(A)') '&STRUCTURE'
        do i = 1, 3
            write(u, '(A,I0,A,3(F20.14,:,","))') '  LatticeVectors(', i, ',:) = ', LAT(i, :)
        end do
        write(u, '(A)') '/'
        write(u, '(A)') '&ATOMS'
        write(u, '(I0)') 2
        do i = 1, 2
            write(u, '(A,3F20.14)') LAB(i), TAU(:, i)
        end do
        write(u, '(A)') '/'
        close(u)
    end subroutine make_material

    ! ------------------------------------------------------------------
    ! 1. Ajuste
    ! ------------------------------------------------------------------
    subroutine test_fit()
        type(HSDataFormat) :: hh, hs
        type(structure_t) :: st
        type(pl_blocks_t) :: pl
        type(sk_model_t) :: m
        type(fit_shell_set_t) :: sets(0)
        type(sk_fit_report_t) :: rep
        character(len=512) :: msg
        integer :: k, kt, nval
        real(dp) :: dv, ma, ra, ms, rs, vexp(3)

        call make_material('tsk_fit', .false., 0.05_dp, 0.4_dp)
        call hh%ReadFile('tsk_fit_H.hsio')
        call hs%ReadFile('tsk_fit_S.hsio')
        call build_pl_blocks(hh, 3, 0, pl, ierr, errmsg, hs=hs)
        call check(ierr == 0 .and. pl%ncell == 1, 'fit: one-cell principal layer')
        call read_structure('tsk_fit.nml', st, ierr, msg)
        call attach_pl_geometry(pl, st, ierr, errmsg)
        m%cutoff = RC
        call sk_set_bases(m, BASIS_TXT, ierr, errmsg)
        call sk_fit_material(pl, st%label, m, sets, 0.4_dp, -1.0_dp, 0.0_dp, 4.0_dp, rep, ierr, errmsg)
        call check(ierr == 0, 'fit: runs')
        if (ierr /= 0) then
            print *, errmsg; return
        end if
        call check(m%nentry == truth%nentry, 'fit: same number of lines as the true model')
        call check_close(rep%beta, BETA_TRUE, 1.0e-6_dp, 'fit: beta recovered')
        call check_close(rep%rel_res_fitted, 0.0_dp, 1.0e-9_dp, 'fit: residual ~ 0')
        call check_close(rep%rel_err_full, 0.0_dp, 1.0e-9_dp, 'fit: full H01 - E0 S01 reproduced')
        dv = 0.0_dp
        do k = 1, m%nentry
            kt = sk_find_entry(truth, m%e(k)%spec_a, m%e(k)%la, m%e(k)%za, m%e(k)%spec_b, m%e(k)%lb, m%e(k)%zb)
            call check(kt > 0, 'fit: line exists in the true model')
            if (kt == 0) cycle
            vexp = truth%e(kt)%v * exp(-BETA_TRUE * (m%e(k)%d0 - D0_TRUE))
            dv = max(dv, maxval(abs(m%e(k)%v - vexp)))
        end do
        call check_close(dv, 0.0_dp, 1.0e-8_dp, 'fit: all V_sigma/pi/delta recovered (s, p, d, hetero)')
        call sk_band_check(pl, st%label, m, [0.1_dp, 0.0_dp], 11, 0.4_dp, 2.0_dp, ma, ra, ms, rs, nval, ierr, errmsg)
        call check(ierr == 0 .and. nval > 0, 'band check: runs')
        call check_close(ms, ma, 1.0e-8_dp, 'band check: SK = anchored coupling (exact fit)')
        call check(ma > 1.0e-6_dp, 'band check: dropping S01 is visible away from E0')
    end subroutine test_fit

    ! ------------------------------------------------------------------
    ! 2-5. Motor
    ! ------------------------------------------------------------------
    subroutine test_engine()
        type(green_data_t) :: g1, g2, g3
        type(green_data_t), allocatable :: gs(:)
        type(green_input_t) :: inp
        character(len=160) :: hdr(1)

        call make_material('tsk_A', .true., 0.0_dp, 0.0_dp)
        hdr(1) = 'true model'
        call write_sk_file('tsk_true.sk', truth%e, hdr, ierr, errmsg)
        call check(ierr == 0, 'write true .sk')

        ! A|A junção: sk = scaled (lambda = 1)
        call write_input('tsk_j_sk.nml', 'junction', 'sk')
        call write_input('tsk_j_sc.nml', 'junction', 'scaled')
        call run('tsk_j_sk.nml', g1)
        call check(ierr == 0, 'A|A junction sk runs')
        if (ierr /= 0) print *, errmsg
        call run('tsk_j_sc.nml', g2)
        call check(ierr == 0, 'A|A junction scaled runs')
        call check_close(maxval(abs(g1%amap - g2%amap)), 0.0_dp, 1.0e-10_dp, &
                         'A|A junction (spinor, tilted a_s): sk contact = exact H01')
        call check_close(maxval(abs(g1%total - g2%total)), 0.0_dp, 1.0e-10_dp, 'A|A junction: total DOS')

        ! sanduíche: usa também o sentido B abaixo de A
        call write_input('tsk_w_sk.nml', 'sandwich', 'sk')
        call write_input('tsk_w_sc.nml', 'sandwich', 'scaled')
        call run('tsk_w_sk.nml', g1)
        call run('tsk_w_sc.nml', g2)
        call check(ierr == 0, 'sandwich runs')
        call check_close(maxval(abs(g1%amap - g2%amap)), 0.0_dp, 1.0e-10_dp, &
                         'A|A|A sandwich: sk contact = exact H01 (both orientations)')

        ! filme com regiões e projeção coerente
        call write_input('tsk_f_sk.nml', 'film', 'sk', coherent=.true.)
        call write_input('tsk_f_sc.nml', 'film', 'scaled', coherent=.true.)
        call run('tsk_f_sk.nml', g1)
        call run('tsk_f_sc.nml', g2)
        call check_close(maxval(abs(g1%amap - g2%amap)), 0.0_dp, 1.0e-10_dp, 'film: sk = exact')
        call check_close(maxval(abs(g1%icoh - g2%icoh)), 0.0_dp, 1.0e-10_dp, 'film: coherent map sk = exact')

        ! 3. lambda_h = 0 = none
        call write_input('tsk_f_l0.nml', 'film', 'sk', lambda=0.0_dp)
        call write_input('tsk_f_none.nml', 'film', 'none')
        call run('tsk_f_l0.nml', g1)
        call run('tsk_f_none.nml', g2)
        call check_close(maxval(abs(g1%amap - g2%amap)), 0.0_dp, 1.0e-12_dp, 'lambda_h = 0 equals none')

        ! 4a. &SCAN lambda_h
        call write_input('tsk_scan_l.nml', 'film', 'sk', scan="parameter = 'lambda_h', values = 0.0, 0.5, 1.0")
        call read_green_input('tsk_scan_l.nml', inp, ierr, errmsg)
        call check(ierr == 0 .and. inp%nscan == 3, 'scan input lambda_h (3 values)')
        call run_green(inp, g1, ierr, errmsg, nonconverged=nbad, gscan=gs)
        call check(ierr == 0 .and. size(gs) == 3, 'scan lambda_h runs (3 outputs)')
        call write_input('tsk_f_l5.nml', 'film', 'sk', lambda=0.5_dp)
        call run('tsk_f_l5.nml', g3)
        call check(maxval(abs(gs(1)%amap - g2%amap)) <= 1.0e-12_dp, 'scan(lambda=0) = none run')
        call check(maxval(abs(gs(2)%amap - g3%amap)) <= 1.0e-12_dp, 'scan(lambda=0.5) = separate run')
        call run('tsk_f_sc.nml', g2)
        call check(maxval(abs(gs(3)%amap - g2%amap)) <= 1.0e-10_dp, 'scan(lambda=1) = exact')
        call check(maxval(abs(gs(1)%amap - g1%amap)) == 0.0_dp, 'gd = first scan value')

        ! 4b. &SCAN film_shift_ev na junção (meio B refeito por valor)
        call write_input('tsk_scan_e.nml', 'junction', 'sk', scan="parameter = 'film_shift_ev', n = 2, vmin = 0.0, vmax = 0.3")
        call read_green_input('tsk_scan_e.nml', inp, ierr, errmsg)
        call run_green(inp, g1, ierr, errmsg, nonconverged=nbad, gscan=gs)
        call check(ierr == 0 .and. size(gs) == 2, 'scan film_shift_ev runs')
        call write_input('tsk_j_e3.nml', 'junction', 'sk', film_shift=0.3_dp)
        call run('tsk_j_e3.nml', g3)
        call check(maxval(abs(gs(2)%amap - g3%amap)) <= 1.0e-12_dp, 'scan(film_shift=0.3) = separate run')
        call check(maxval(abs(gs(2)%niter - g3%niter)) == 0, 'scan film shift: B lead recomputed')

        ! 4c. &SCAN offset_normal (geometria do contato e posições)
        call write_input('tsk_scan_o.nml', 'film', 'sk', coherent=.true., &
                         scan="parameter = 'offset_normal', values = 0.0, 0.4")
        call read_green_input('tsk_scan_o.nml', inp, ierr, errmsg)
        call run_green(inp, g1, ierr, errmsg, nonconverged=nbad, gscan=gs)
        call check(ierr == 0, 'scan offset_normal runs')
        call write_input('tsk_f_o4.nml', 'film', 'sk', coherent=.true., offset_z=0.4_dp)
        call run('tsk_f_o4.nml', g3)
        call check(maxval(abs(gs(2)%amap - g3%amap)) <= 1.0e-12_dp, 'scan(offset_normal=0.4) = separate run')
        call check(maxval(abs(gs(2)%icoh - g3%icoh)) <= 1.0e-12_dp, 'scan offset: coherent geometry follows')
        call check(maxval(abs(gs(2)%amap - gs(1)%amap)) > 1.0e-6_dp, 'scan offset: contact changes the spectrum')

        ! 5. causalidade com contato genérico
        call write_input('tsk_f_gen.nml', 'film', 'sk', lambda=0.7_dp, offset_z=0.25_dp, film_shift=-0.2_dp)
        call run('tsk_f_gen.nml', g1)
        call check(ierr == 0, 'generic sk film runs')
        call check(minval(g1%min_neg_im_gii) > -1.0e-10_dp, 'generic sk film: -Im G_ii >= 0')

        ! erros
        call write_input('tsk_bad1.nml', 'film', 'sk', basis="sk_basis = 'Xa s1p1'")
        call run('tsk_bad1.nml', g1)
        call check(ierr /= 0 .and. index(errmsg, 'Xb') > 0, 'missing species basis is an error')
        call write_input('tsk_bad2.nml', 'film', 'sk', basis="sk_basis = 'Xa s1p1', 'Xb s1p2d1'")
        call run('tsk_bad2.nml', g1)
        call check(ierr /= 0, 'wrong orbital count for a basis is an error')
        call write_input('tsk_bad3.nml', 'film', 'none', scan="parameter = 'sk_scale_pi', values = 1.0")
        call read_green_input('tsk_bad3.nml', inp, ierr, errmsg)
        call check(ierr /= 0, "scan of sk_scale_pi needs model = 'sk'")
        call write_input('tsk_bad4.nml', 'film', 'sk', scan="parameter = 'nonsense', values = 1.0")
        call read_green_input('tsk_bad4.nml', inp, ierr, errmsg)
        call check(ierr /= 0, 'unknown scan parameter is an error')
    end subroutine test_engine

    subroutine run(nml, gd)
        character(len=*), intent(in) :: nml
        type(green_data_t), intent(out) :: gd
        type(green_input_t) :: inp
        call read_green_input(nml, inp, ierr, errmsg)
        if (ierr /= 0) return
        call run_green(inp, gd, ierr, errmsg, nonconverged=nbad)
    end subroutine run

    subroutine write_input(fname, mode, model, coherent, lambda, scan, film_shift, offset_z, basis)
        character(len=*), intent(in) :: fname, mode, model
        logical, intent(in), optional :: coherent
        real(dp), intent(in), optional :: lambda, film_shift, offset_z
        character(len=*), intent(in), optional :: scan, basis
        integer :: u
        open(newunit=u, file=fname, status='replace', action='write')
        write(u, '(A)') '&GREEN'
        write(u, '(A,A,A)') "  mode = '", mode, "'"
        write(u, '(A)') '  eta_ev = 0.05, emin = -1.0, emax = 1.0, ne = 7'
        write(u, '(A)') '  decim_tol = 1e-12, decim_maxit = 300'
        write(u, '(A,A,A)') "  output_prefix = '", fname(1:len(fname) - 4), "'"
        write(u, '(A)') '/'
        write(u, '(A)') '&SUBSTRATE'
        write(u, '(A)') "  h_file = 'tsk_A_H.hsio', s_file = 'tsk_A_S.hsio', structure_file = 'tsk_A.nml'"
        write(u, '(A)') "  stacking_axis = 3, central_layers = 1, side = 'top'"
        write(u, '(A)') '/'
        write(u, '(A)') '&FILM'
        write(u, '(A)') "  h_file = 'tsk_A_H.hsio', s_file = 'tsk_A_S.hsio', structure_file = 'tsk_A.nml'"
        write(u, '(A)') '  nlayers = 2'
        if (present(film_shift)) write(u, '(A,F12.6)') '  shift_ev = ', film_shift
        if (present(offset_z)) write(u, '(A,F12.6)') '  interface_offset = 0.0, 0.0, ', offset_z
        write(u, '(A)') '/'
        write(u, '(A)') '&INTERFACE'
        write(u, '(A,A,A)') "  model = '", model, "', reference = 'substrate'"
        if (present(lambda)) write(u, '(A,F12.6)') '  lambda_h = ', lambda
        write(u, '(A,F8.3)') "  sk_file = 'tsk_true.sk', sk_cutoff_ang = ", RC
        if (present(basis)) then
            write(u, '(2X,A)') basis
        else
            write(u, '(A)') "  sk_basis = 'Xa s1p1', 'Xb s1p1d1'"
        end if
        write(u, '(A)') '/'
        if (present(scan)) then
            write(u, '(A)') '&SCAN'
            write(u, '(2X,A)') scan
            write(u, '(A)') '/'
        end if
        if (present(coherent)) then
            write(u, '(A)') '&COHERENT'
            write(u, '(A)') '  nqz = 3, qz_min = 0.0, qz_max = 1.2, escape_length_ang = 5.0'
            write(u, '(A)') '/'
        end if
        write(u, '(A)') '&KPATH'
        write(u, '(A)') '  nvertex = 3'
        write(u, '(A)') '  vertex(:,1) = 0.0, 0.0, 0.0'
        write(u, '(A)') '  vertex(:,2) = 0.5, 0.0, 0.0'
        write(u, '(A)') '  vertex(:,3) = 0.3333333333, 0.3333333333, 0.0'
        write(u, '(A)') "  label = 'G', 'M', 'K'"
        write(u, '(A)') '  npoints = 2, 2'
        write(u, '(A)') '/'
        close(u)
    end subroutine write_input

    subroutine cleanup()
        call execute_command_line('rm -f tsk_*.hsio tsk_*.nml tsk_*.sk tsk_*.green.bin', wait=.true.)
    end subroutine cleanup

end program test_sk_contact
