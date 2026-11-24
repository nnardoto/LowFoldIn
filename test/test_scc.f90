! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Testes do ciclo autoconsistente (mod_scc, A5):
!>  1. Gamma de duas esferas uniformes (raio a) = fórmula analítica, U = 6/(5a), 1/R fora;
!>  2. p zeta1 do Te7.0.pao: norma, U_H > 0, Gamma(R -> 2 rc) = 1/R, monótona;
!>  3. Ewald 2D: independência de alpha; referência no fundo (z -> -inf: 0; z -> +inf: -(4 pi/A) z);
!>     par neutro no plano = soma direta;
!>  4. Anderson: ponto fixo linear não contrativo converge.
!>  5. Dipolos (A6): forma do dipolo (momentos), g' e g'' contra diferenças finitas de Gamma, T on-site,
!>     T -> dipolo pontual, Ewald de dipolos (alpha, Laplace, soma direta), integrais <mu|r|nu> nas PAO
!>     (S(0) = 1, X(0) simétrica com diagonal nula, sentido oposto exato) e o termo de H como dE/drho.
!>  6. A7: Ewald 2D sem substrato (sym): alpha, simetria em z, diferença (2 pi/A) z para a referência no fundo;
!>     campo pontual nos núcleos (esfera uniforme, analítico) = soma direta, com termo próprio.
!>  7. A10 (cargas de ligação): termo de H do modo 'bond' (on-site V_i S + F_i.X^(i), par V_b S) = dE/drho num
!>     modelo de dois átomos e um par; Gamma de sítios de par (pontos médios, dois sítios no mesmo ponto) e
!>     núcleo carga-dipolo cruzado = soma direta na rede 2D.
!>  8. A14 (substrate = 'dielectric'): beta = 0 não altera os núcleos (bit a bit); núcleos de imagem gk, gkf, tk
!>     = diferenças finitas de gam_img em R_i e R_j; carga pontual diante do meio: potencial e campo da imagem
!>     (-beta q no espelho) = soma direta na rede 2D (método das imagens, analítico).
!>  9. A19 (planar_full): perfil planar das PAO soma a S de pao_two_center; carga planar conserva
!>     sum_q tr(D_q S_q); dH^q / wh_q = sum_k dV_k W^{q,k} = dE/dD^q por diferenças finitas (E = 1/2 c K c, núcleo 1D
!>     simétrico); D = D_ref -> c = 0, dV = 0; formas planares de rho_s: norma 1, dipolo, e potencial da folha.
program test_scc
    use types, only: dp
    use scc
    use pao_basis, only: pao_species_t, read_pao_species, pao_two_center, pao_dip_cache_t, pao_dip_ensure, &
                         pao_dip_get, pao_dipole_dh, pao_planar_profile
    use scc_planar, only: planar_t, planar_setup, planar_charge, planar_dh, planar_mult, planar_pot
    implicit none
    real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp
    integer :: nfail
    nfail = 0
    call t_sphere()
    call t_pao()
    call t_ewald()
    call t_anderson()
    call t_dipole_radial()
    call t_dipole_ewald()
    call t_dipole_integrals()
    call t_ewald_sym()
    call t_point_field()
    call t_bond_model()
    call t_bond_ewald()
    call t_dielectric()
    call t_planar()
    if (nfail > 0) then
        print '(A,I0,A)', 'test_scc: ', nfail, ' FAILED'
        error stop 1
    end if
    print '(A)', 'test_scc: all passed'
contains

    !> A19: tensor planar, carga, termo de H (dE/dD) e formas dos multipolos.
    subroutine t_planar()
        real(dp), parameter :: BA = 0.529177210903_dp, HA = 27.211386245988_dp
        type(pao_species_t) :: sp
        type(planar_t) :: pl
        type(radial_density_t) :: rd
        integer :: ierr, no, nn, q, k, a, b, t, ip(4), jp(2)
        character(len=:), allocatable :: msg
        real(dp), allocatable :: p(:, :, :), s1(:, :), x1(:, :, :), d(:, :, :), dh(:, :, :), c(:), v(:), c0(:), &
                                 dd(:, :, :), sg(:, :, :)
        real(dp) :: dv(3), ch(2), pv(3, 4), e0, ep, em, fd, emax, tr, h, eps, dpn(2), dq(2)
        integer :: kind(4), ic(4)
        call read_pao_species('test/data/pao/Te7.0.pao', 's3p3d2f1', 7.0_dp, sp, ierr, msg)
        if (ierr /= 0) then
            call check(.false., 'planar: read pao'); return
        end if
        no = sp%norb
        h = 0.3_dp
        nn = floor(min(sp%rcut, sp%rv(sp%nmesh)) / h)
        allocate(p(no, no, 2 * nn + 1), s1(no, no), x1(no, no, 3))
        dv = [2.1_dp, -3.4_dp, 1.7_dp]
        call pao_planar_profile(sp, dv, h, 3, p)
        call pao_two_center(sp, dv, h, s1, x1)
        e0 = maxval(abs(sum(p, dim=3) - s1))
        print '(A,ES10.3)', '    planar profile: max|sum_k P - S(two_center)| = ', e0
        call check(e0 < 1.0e-13_dp, 'planar: sum of slab profiles = S of pao_two_center (same grid)')
        ! modelo: dois átomos centrais (alturas 0 e 1,6 Å ao longo de y), pares (1,1,0), (2,2,0), (1,2,d), (2,1,-d)
        ! e um par com a costura (1, j, d2) (peso 2)
        ch = [0.0_dp, 1.6_dp]
        pv(:, 1) = 0.0_dp
        pv(:, 2) = 0.0_dp
        pv(:, 3) = [1.3_dp, 1.6_dp, -0.4_dp]
        pv(:, 4) = -pv(:, 3)
        ic = [1, 2, 1, 2]
        kind = [1, 1, 1, 1]
        call planar_setup(pl, sp, h, 2, 1.0_dp, ch, [ic, 1], reshape([pv, [-0.9_dp, -2.4_dp, 0.7_dp]], [3, 5]), &
                          [kind, 4], 0, 0.1_dp, 3.0_dp, 20.0_dp, rhof(), 0.01_dp, ierr, msg)
        call check(ierr == 0, 'planar: setup')
        if (ierr /= 0) return
        call check(pl%nprof == 3 .and. pl%ptr(4) .neqv. pl%ptr(3), 'planar: profiles by canonical d (transpose)')
        allocate(d(no, no, 5), dh(no, no, 5), c(pl%ny), v(pl%ny), c0(pl%ny), dd(no, no, 5), sg(no, no, 5))
        ! D sintética: on-site simétricas, D(2,1,-d) = D(1,2,d)^T
        do q = 1, 5
            do b = 1, no
                do a = 1, no
                    d(a, b, q) = 0.01_dp * sin(1.3_dp * a + 0.7_dp * b + 2.1_dp * q)
                end do
            end do
        end do
        d(:, :, 1) = d(:, :, 1) + transpose(d(:, :, 1))
        d(:, :, 2) = d(:, :, 2) + transpose(d(:, :, 2))
        d(:, :, 4) = transpose(d(:, :, 3))
        call planar_charge(pl, d, c)
        tr = 0.0_dp
        do q = 1, 5
            sg(:, :, q) = sum(pl%prof(:, :, :, pl%pidx(q)), dim=3)
            if (pl%ptr(q)) sg(:, :, q) = transpose(sg(:, :, q))
            tr = tr + sum(d(:, :, q) * sg(:, :, q))
        end do
        print '(A,ES10.3)', '    planar: |sum_k c_k - sum_q tr(D S_grid)| = ', abs(sum(c) - tr)
        call check(abs(sum(c) - tr) < 1.0e-12_dp, 'planar: charge conservation sum_k c_k = sum_q tr(D_q S_q)')
        ! pares (1,2,d) e (2,1,-d) com D transposta dão a mesma carga
        dd = 0.0_dp
        dd(:, :, 3) = d(:, :, 3)
        call planar_charge(pl, dd, c0)
        dd = 0.0_dp
        dd(:, :, 4) = d(:, :, 4)
        call planar_charge(pl, dd, v)
        print '(A,ES10.3)', '    planar: max|c(1,2,d) - c(2,1,-d)| = ', maxval(abs(c0 - v))
        call check(maxval(abs(c0 - v)) < 1.0e-6_dp * maxval(abs(c0)), 'planar: pair and reverse pair (grid error)')
        ! D = D_ref -> c = 0, dV = 0
        call planar_charge(pl, d, c0)
        call planar_pot(pl, c - c0, v)
        call check(maxval(abs(v)) == 0.0_dp, 'planar: D = D_ref gives dV = 0 (bit for bit)')
        ! dH = dE/dD (E = 1/2 c K_sym c), diferenças finitas em elementos de todos os pares
        call planar_pot(pl, c, v, sym=.true.)
        call planar_dh(pl, v, dh)
        e0 = energy(pl, d)
        emax = 0.0_dp
        eps = 1.0e-4_dp
        ip = [3, 7, 12, 25]
        jp = [5, 18]
        do q = 1, 5
            do a = 1, 4
                do b = 1, 2
                    dd = d
                    dd(ip(a), jp(b), q) = d(ip(a), jp(b), q) + eps
                    ep = energy(pl, dd)
                    dd(ip(a), jp(b), q) = d(ip(a), jp(b), q) - eps
                    em = energy(pl, dd)
                    fd = (ep - em) / (2.0_dp * eps)
                    emax = max(emax, abs(fd - dh(ip(a), jp(b), q) / pl%wh(q)) / max(abs(fd), 1.0e-3_dp))
                end do
            end do
        end do
        print '(A,ES10.3)', '    planar: max rel |dE/dD (FD) - dH / wh_q| = ', emax
        call check(pl%wh(5) == 0.5_dp .and. all(pl%wh(1:4) == 1.0_dp), 'planar: H weight 1/2 for seam pairs')
        call check(emax < 1.0e-7_dp, 'planar: dH^q = dE/dD^q (finite differences)')
        ! multipolos: norma, dipolo, potencial acima de uma folha
        print '(A,2ES11.3,A,ES11.3)', '    planar: sum fq - 1 = ', sum(pl%fq(:, 1)) - 1.0_dp, sum(pl%fq(:, 2)) - &
            1.0_dp, ' ; sum y fp + Bohr = ', sum(pl%y * pl%fp(:, 1)) + BA
        call check(abs(sum(pl%fq(:, 1)) - 1.0_dp) < 1.0e-6_dp .and. abs(sum(pl%fp(:, 2))) < 1.0e-12_dp .and. &
                   abs(sum(pl%y * pl%fq(:, 2)) - 1.6_dp) < 1.0e-6_dp .and. &
                   abs(sum(pl%y * pl%fp(:, 1)) + BA) < 1.0e-6_dp, 'planar: rho_s forms (norm, center, dipole)')
        dq = [1.0_dp, 0.0_dp]
        dpn = 0.0_dp
        call planar_mult(pl, dq, dpn, c)
        call planar_pot(pl, c, v)
        e0 = (v(pl%ny) - v(pl%ny - 10)) / (pl%y(pl%ny) - pl%y(pl%ny - 10))
        print '(A,ES12.4,A,ES12.4)', '    planar: slope above = ', e0, ' ; -4 pi ke/A = ', -4.0_dp * PI * HA * BA / 20.0_dp
        call check(abs(e0 + 4.0_dp * PI * HA * BA / 20.0_dp * sum(c)) < 1.0e-10_dp .and. abs(v(1)) < 1.0e-14_dp, &
                   'planar: sheet potential (zero below, -(4 pi ke/A) h above)')
        t = 0
    end subroutine t_planar

    real(dp) function energy(pl, dm)
        type(planar_t), intent(in) :: pl
        real(dp), intent(in) :: dm(:, :, :)
        real(dp), allocatable :: cc(:), vv(:)
        allocate(cc(pl%ny), vv(pl%ny))
        call planar_charge(pl, dm, cc)
        call planar_pot(pl, cc, vv, sym=.true.)
        energy = 0.5_dp * dot_product(cc, vv)
    end function energy

    function rhof() result(r)
        real(dp) :: r(0:300)
        integer :: kk
        do kk = 0, 300
            r(kk) = exp(-(0.01_dp * kk / 0.8_dp)**2)
        end do
        r = r / (sum([(4.0_dp * PI * (0.01_dp * kk)**2 * r(kk) * 0.01_dp, kk = 0, 300)]))
    end function rhof

    !> A14: núcleos de imagem do meio dielétrico (scc_image_kernels).
    subroutine t_dielectric()
        real(dp), parameter :: HA = 27.211386245988_dp, BA = 0.529177210903_dp, KE = HA * BA
        type(radial_density_t) :: rd
        integer :: ierr, n, a, n1, n2, nm
        character(len=:), allocatable :: msg
        real(dp) :: a1(3), a2(3), nh(3), pos(3, 2), p2(3, 2), q(3, 2), beta, zp, e, eg, et, ml(3, 3)
        real(dp) :: gam(2, 2), gk(3, 2, 2), gkf(3, 2, 2), tk(3, 3, 2, 2)
        real(dp) :: gp(2, 2), gm(2, 2), fp(3, 2, 2), fm(3, 2, 2), dum(3, 2, 2), dt(3, 3, 2, 2)
        real(dp) :: rv(3), r, vs, vh, gs(3), gh(3), m1(3), m2(3), x(3), v, ex(3), hh
        n = 3000
        rd%n = n; rd%dr = 3.0_dp / n; rd%rc = 3.0_dp         ! esfera uniforme de raio 3 Bohr
        allocate(rd%rho(0:n))
        rd%rho = 1.0_dp
        call finish_radial_density(rd, ierr, msg)
        a1 = [4.7672460267921517_dp, 0.0_dp, 0.0_dp]
        a2 = [0.0_dp, 0.0_dp, 5.9335974928231314_dp]
        nh = [0.0_dp, -1.0_dp, 0.0_dp]
        zp = 0.7_dp
        pos(:, 1) = [0.3_dp, -2.4_dp, 0.2_dp]      ! h = 2,4 Å (imagem a 3,4 Å: dentro do alcance de gamma)
        pos(:, 2) = [2.1_dp, -3.1_dp, 3.0_dp]
        beta = (10.0_dp - 1.0_dp) / (10.0_dp + 1.0_dp)
        ! (a) beta = 0: nada muda (bit a bit)
        gam = 1.25_dp; gk = -0.5_dp; gkf = 0.75_dp; tk = 2.0_dp
        call scc_image_kernels(rd, pos, a1, a2, nh, 0.5_dp, zp, 0.0_dp, gam, gk, gkf, tk)
        call check(all(gam == 1.25_dp) .and. all(gk == -0.5_dp) .and. all(gkf == 0.75_dp) .and. all(tk == 2.0_dp), &
                   'dielectric: beta = 0 leaves the kernels bit for bit')
        ! (b) gk, gkf, tk = diferenças finitas centrais de gam_img (h = 1e-3 Å)
        call img(rd, a1, a2, nh, zp, beta, pos, gam, gk, gkf, tk)
        eg = 0.0_dp; et = 0.0_dp; hh = 1.0e-3_dp
        do a = 1, 3
            ex = 0.0_dp
            ex(a) = hh
            q = pos; q(:, 2) = q(:, 2) + ex; call img(rd, a1, a2, nh, zp, beta, q, gp, dum, fp, dt)
            q = pos; q(:, 2) = q(:, 2) - ex; call img(rd, a1, a2, nh, zp, beta, q, gm, dum, fm, dt)
            eg = max(eg, abs(gk(a, 1, 2) + (gp(1, 2) - gm(1, 2)) / (2.0_dp * hh) * BA))   ! gk = -d gam/d R_j
            et = max(et, maxval(abs(tk(:, a, 1, 2) - (fp(:, 1, 2) - fm(:, 1, 2)) / (2.0_dp * hh) * BA)))
            q = pos; q(:, 1) = q(:, 1) + ex; call img(rd, a1, a2, nh, zp, beta, q, gp, dum, fp, dt)
            q = pos; q(:, 1) = q(:, 1) - ex; call img(rd, a1, a2, nh, zp, beta, q, gm, dum, fm, dt)
            eg = max(eg, abs(gkf(a, 1, 2) - (gp(1, 2) - gm(1, 2)) / (2.0_dp * hh) * BA))  ! gkf = d gam/d R_i
        end do
        print '(A,2ES10.3)', '    image kernels vs finite differences (grad, hess) ', eg, et
        call check(eg < 1.0e-5_dp .and. et < 1.0e-5_dp, 'dielectric: gk, gkf, tk of the image = d gam_img (FD)')
        ! (c) carga pontual diante do meio: +1 em R_1 e -1 em R_2 na mesma altura (2,4 Å acima do plano: imagens a
        ! >= 4,8 Å, fora do alcance 2 rc = 3,2 Å de gamma -> 1/r). V_img(R_1) = -beta sum_L [1/|R_1 - M R_1 + L| - 1/|R_1 - M R_2 + L|], campo idem.
        p2(:, 1) = [0.3_dp, -3.1_dp, 0.2_dp]
        p2(:, 2) = [2.1_dp, -3.1_dp, 3.0_dp]
        call img(rd, a1, a2, nh, zp, beta, p2, gam, gk, gkf, tk)
        v = gam(1, 1) - gam(1, 2)
        x = gkf(:, 1, 1) - gkf(:, 1, 2)
        m1 = p2(:, 1) - 2.0_dp * (dot_product(p2(:, 1), nh) - zp) * nh
        m2 = p2(:, 2) - 2.0_dp * (dot_product(p2(:, 2), nh) - zp) * nh
        nm = 800
        vs = 0.0_dp; vh = 0.0_dp; gs = 0.0_dp; gh = 0.0_dp
        do n1 = -nm, nm
            do n2 = -nm, nm
                rv = real(n1, dp) * a1 + real(n2, dp) * a2
                r = norm2(rv)
                if (r > 3600.0_dp) cycle
                call pt_sum(p2(:, 1) - m1 + rv, 1.0_dp, r <= 1800.0_dp, vs, vh, gs, gh)
                call pt_sum(p2(:, 1) - m2 + rv, -1.0_dp, r <= 1800.0_dp, vs, vh, gs, gh)
            end do
        end do
        vs = -beta * KE * (2.0_dp * vs - vh)
        gs = -beta * KE * BA * (2.0_dp * gs - gh)
        e = max(abs(v - vs), maxval(abs(x - gs)))
        print '(A,2F14.8,ES10.3,A,3F10.6)', '    point charge before the medium: V_img kernel, direct, |diff| ', v, vs, e, &
            ' ; F_img ', x
        call check(e < 1.0e-5_dp, 'dielectric: point charge before the medium = images (direct sum)')
        ! regra de sinal do dipolo imagem -beta Ml p: Ml troca só a componente normal
        ml = 0.0_dp
        do a = 1, 3
            ml(a, a) = 1.0_dp
            ml(:, a) = ml(:, a) - 2.0_dp * nh * nh(a)
        end do
        call check(abs(ml(2, 2) + 1.0_dp) < 1.0e-15_dp .and. abs(ml(1, 1) - 1.0_dp) < 1.0e-15_dp, &
                   'dielectric: mirror flips the normal component only')
    end subroutine t_dielectric

    subroutine img(rd, a1, a2, nh, zp, beta, pp, g, k, kf, t)
        type(radial_density_t), intent(in) :: rd
        real(dp), intent(in) :: a1(3), a2(3), nh(3), zp, beta, pp(3, 2)
        real(dp), intent(out) :: g(2, 2), k(3, 2, 2), kf(3, 2, 2), t(3, 3, 2, 2)
        g = 0.0_dp; k = 0.0_dp; kf = 0.0_dp; t = 0.0_dp
        call scc_image_kernels(rd, pp, a1, a2, nh, 0.5_dp, zp, beta, g, k, kf, t)
    end subroutine img

    subroutine pt_sum(dd, qq, inner, vs, vh, gs, gh)
        real(dp), intent(in) :: dd(3), qq
        logical, intent(in) :: inner
        real(dp), intent(inout) :: vs, vh, gs(3), gh(3)
        real(dp) :: rr
        rr = norm2(dd)
        vs = vs + qq / rr
        gs = gs - qq * dd / rr**3
        if (inner) then
            vh = vh + qq / rr
            gh = gh - qq * dd / rr**3
        end if
    end subroutine pt_sum

    subroutine check(cond, msg)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: msg
        if (cond) then
            print '(A,A)', '  ok   ', msg
        else
            print '(A,A)', '  FAIL ', msg
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine t_sphere()
        type(radial_density_t) :: rd
        integer :: ierr, n, i
        character(len=:), allocatable :: msg
        real(dp) :: a, x, ex, g, err
        a = 3.0_dp
        n = 6000
        rd%n = n; rd%dr = a / n; rd%rc = a
        allocate(rd%rho(0:n))
        rd%rho = 1.0_dp
        call finish_radial_density(rd, ierr, msg)
        call check(ierr == 0, 'sphere: finish')
        err = 0.0_dp
        do i = 0, 8
            x = 0.25_dp * i             ! R/a de 0 a 2
            ex = (1.0_dp / a) * (1.2_dp - 0.5_dp * x**2 + 3.0_dp / 16.0_dp * x**3 - x**5 / 160.0_dp)
            g = gamma_au(rd, x * a)
            err = max(err, abs(g - ex) / ex)
        end do
        print '(A,ES10.3)', '    sphere max rel err = ', err
        call check(err < 1.0e-5_dp, 'sphere: Gamma(R) analytic (0 <= R <= 2a), U = 6/(5a)')
        call check(abs(gamma_au(rd, 7.0_dp) - 1.0_dp / 7.0_dp) < 1.0e-14_dp, 'sphere: 1/R beyond 2a')
    end subroutine t_sphere

    subroutine t_pao()
        type(pao_species_t) :: sp
        type(radial_density_t) :: rd
        integer :: ierr, i
        character(len=:), allocatable :: msg
        real(dp) :: r, g, gp, uh
        logical :: mono
        call read_pao_species('test/data/pao/Te7.0.pao', 's3p3d2f1', 7.0_dp, sp, ierr, msg)
        call check(ierr == 0, 'pao: read Te7.0.pao s3p3d2f1')
        if (ierr /= 0) return
        call build_radial_density(sp, 1, 1, 5.0e-4_dp, rd, ierr, msg)
        call check(ierr == 0, 'pao: p zeta1 density')
        print '(A,F12.8)', '    norm before normalization = ', rd%norm
        call check(abs(rd%norm - 1.0_dp) < 1.0e-3_dp, 'pao: p zeta1 norm ~ 1')
        uh = gamma_au(rd, 0.0_dp)
        print '(A,F12.6,A,F12.6,A)', '    U_H = ', uh, ' Ha = ', uh * 27.211386245988_dp, ' eV'
        call check(uh > 0.2_dp .and. uh < 1.0_dp, 'pao: U_H in a plausible range (0.2..1 Ha)')
        r = 2.0_dp * rd%rc - 1.0e-6_dp
        g = gamma_au(rd, r)
        print '(A,ES10.3)', '    |Gamma(2rc-) R - 1| = ', abs(g * r - 1.0_dp)
        call check(abs(g * r - 1.0_dp) < 1.0e-6_dp, 'pao: Gamma(R -> 2 rc) = 1/R')
        mono = .true.
        gp = uh
        do i = 1, 140
            g = gamma_au(rd, 0.1_dp * i)
            if (g > gp + 1.0e-14_dp .or. g > 1.0_dp / (0.1_dp * i) + 1.0e-12_dp) mono = .false.
            gp = g
        end do
        call check(mono, 'pao: Gamma decreasing and <= 1/R')
    end subroutine t_pao

    subroutine t_ewald()
        real(dp) :: a1(3), a2(3), nh(3), d(3), e1, e2, area, z, ds, dh, t, rv(3), r1(3), r2(3), md
        integer :: n1, n2, nm
        a1 = [4.7672460267921517_dp, 0.0_dp, 0.0_dp]
        a2 = [0.0_dp, 0.0_dp, 5.9335974928231314_dp]
        nh = [0.0_dp, -1.0_dp, 0.0_dp]      ! a1 x a2 = -y; normal para o vácuo em -y
        area = a1(1) * a2(3)
        d = [1.3_dp, -2.2_dp, 0.7_dp]
        e1 = ewald2d_ref(d, a1, a2, nh, 0.3_dp, .false.)
        e2 = ewald2d_ref(d, a1, a2, nh, 0.7_dp, .false.)
        print '(A,ES10.3)', '    alpha 0.3 vs 0.7: ', abs(e1 - e2)
        call check(abs(e1 - e2) < 1.0e-12_dp, 'ewald: independent of alpha (pair)')
        e1 = ewald2d_ref([0.0_dp, 0.0_dp, 0.0_dp], a1, a2, nh, 0.3_dp, .true.)
        e2 = ewald2d_ref([0.0_dp, 0.0_dp, 0.0_dp], a1, a2, nh, 0.7_dp, .true.)
        call check(abs(e1 - e2) < 1.0e-12_dp, 'ewald: independent of alpha (self)')
        z = -40.0_dp       ! i bem abaixo de j (lado do substrato)
        e1 = ewald2d_ref([0.5_dp, 0.0_dp, 0.0_dp] + z * nh, a1, a2, nh, 0.5_dp, .false.)
        call check(abs(e1) < 1.0e-12_dp, 'ewald: zero deep in the substrate')
        z = 40.0_dp
        e1 = ewald2d_ref([0.5_dp, 0.0_dp, 0.0_dp] + z * nh, a1, a2, nh, 0.5_dp, .false.)
        call check(abs(e1 + 4.0_dp * PI / area * z) < 1.0e-10_dp, 'ewald: -(4 pi/A) z on the vacuum side')
        ! par neutro no plano (z1 = z2): Ewald = soma direta (cauda ~ 1/R^2)
        r1 = [0.3_dp, 0.0_dp, 0.2_dp]
        r2 = [2.1_dp, 0.0_dp, 3.0_dp]
        d = [1.0_dp, 1.7_dp, 0.4_dp]
        e1 = ewald2d_ref(d - r1, a1, a2, nh, 0.5_dp, .false.) - ewald2d_ref(d - r2, a1, a2, nh, 0.5_dp, .false.)
        ! soma direta num disco de raio R: erro ~ 1/R (quadrupolo) -> Richardson com R e 2R
        nm = 800
        ds = 0.0_dp
        dh = 0.0_dp
        do n1 = -nm, nm
            do n2 = -nm, nm
                rv = real(n1, dp) * a1 + real(n2, dp) * a2
                if (norm2(rv) > 3600.0_dp) cycle
                t = 1.0_dp / norm2(d - r1 + rv) - 1.0_dp / norm2(d - r2 + rv)
                ds = ds + t
                if (norm2(rv) <= 1800.0_dp) dh = dh + t
            end do
        end do
        ds = 2.0_dp * ds - dh
        md = abs(e1 - ds)
        print '(A,2F16.10,ES10.3)', '    neutral in-plane pair: ewald, direct, diff ', e1, ds, md
        call check(md < 1.0e-5_dp, 'ewald: neutral in-plane pair = direct sum')
    end subroutine t_ewald

    subroutine t_anderson()
        type(anderson_t) :: st
        real(dp) :: a(4, 4), b(4), x(4), f(4), xs(4)
        integer :: it, i
        a = 0.0_dp
        do i = 1, 4
            a(i, i) = -3.0_dp * i       ! g(x) = A x + b, mistura linear com beta=0.5 diverge
        end do
        a(1, 2) = 0.5_dp; a(3, 4) = -0.7_dp
        b = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp]
        st%hist = 6; st%beta = 0.5_dp
        x = 0.0_dp
        do it = 1, 30
            f = matmul(a, x) + b - x
            if (maxval(abs(f)) < 1.0e-12_dp) exit
            call anderson_step(st, x, f)
        end do
        ! solução exata: (I - A) x = b
        xs = x
        print '(A,I0,A,ES10.3)', '    anderson: ', it, ' iterations, residual ', maxval(abs(matmul(a, xs) + b - xs))
        call check(maxval(abs(matmul(a, xs) + b - xs)) < 1.0e-10_dp .and. it < 15, 'anderson: linear problem')
    end subroutine t_anderson

    subroutine t_dipole_radial()
        type(pao_species_t) :: sp
        type(radial_density_t) :: rd
        integer :: ierr, k, i
        character(len=:), allocatable :: msg
        real(dp) :: m1, r, dro, g1, g2, h, fd1, fd2, e1, e2, rr(4), t0
        call read_pao_species('test/data/pao/Te7.0.pao', 's3p3d2f1', 7.0_dp, sp, ierr, msg)
        if (ierr /= 0) then
            call check(.false., 'dipole: read pao'); return
        end if
        call build_radial_density(sp, 1, 1, 5.0e-4_dp, rd, ierr, msg)
        ! forma do dipolo -p.grad rho_s: momento zero nulo por simetria; primeiro momento
        ! int x (-p_x rho'(r) x/r) d^3r = -(4 pi/3) p_x int rho' r^3 dr = p_x
        m1 = 0.0_dp
        do k = 1, rd%n - 1
            r = real(k, dp) * rd%dr
            dro = (rd%rho(k + 1) - rd%rho(k - 1)) / (2.0_dp * rd%dr)
            m1 = m1 - 4.0_dp * PI / 3.0_dp * dro * r**3 * rd%dr
        end do
        print '(A,ES10.3)', '    dipole shape: first moment / p - 1 = ', m1 - 1.0_dp
        call check(abs(m1 - 1.0_dp) < 1.0e-4_dp, 'dipole: first moment of -p.grad rho_s = p')
        ! g', g'' contra diferenças finitas de gamma_au
        rr = [3.0_dp, 5.0_dp, 9.0_dp, 12.5_dp]
        e1 = 0.0_dp; e2 = 0.0_dp
        h = 0.04_dp
        do i = 1, 4
            call gamma_derivs_au(rd, rr(i), g1, g2)
            ! diferenças centrais com Richardson (h e h/2)
            fd1 = (4.0_dp * (gamma_au(rd, rr(i) + h / 2) - gamma_au(rd, rr(i) - h / 2)) / h - &
                   (gamma_au(rd, rr(i) + h) - gamma_au(rd, rr(i) - h)) / (2.0_dp * h)) / 3.0_dp
            fd2 = (16.0_dp * (gamma_au(rd, rr(i) + h / 2) - 2.0_dp * gamma_au(rd, rr(i)) + gamma_au(rd, rr(i) - h / 2)) &
                   / h**2 - (gamma_au(rd, rr(i) + h) - 2.0_dp * gamma_au(rd, rr(i)) + gamma_au(rd, rr(i) - h)) / h**2) / 3.0_dp
            e1 = max(e1, abs(g1 - fd1) / abs(fd1))
            e2 = max(e2, abs(g2 - fd2) / abs(fd2))
        end do
        print '(A,2ES10.3)', '    g1, g2 vs finite differences (max rel): ', e1, e2
        call check(e1 < 1.0e-5_dp .and. e2 < 1.0e-4_dp, 'dipole: grad and hessian of Gamma = finite differences')
        ! on-site: T0 = (4 pi/3) int rho^2 = -g2(0)
        t0 = onsite_dipole_au(rd)
        call gamma_derivs_au(rd, 0.05_dp, g1, g2)
        print '(A,F10.6,A,F10.6,A,F10.5,A)', '    T_ii = ', t0, ' Ha/Bohr^2 ; -g2(R=0.05) = ', -g2, ' ; T_ii = ', &
            t0 * 27.211386245988_dp, ' eV/Bohr^2'
        call check(abs(t0 + g2) / t0 < 2.0e-3_dp, 'dipole: on-site T = -g2(0)')
        ! longe: g1 = -1/R^2, g2 = 2/R^3 (dipolo pontual), já em R -> 2 rc
        r = 2.0_dp * rd%rc - 0.1_dp
        call gamma_derivs_au(rd, r, g1, g2)
        print '(A,2ES10.3)', '    R = 2rc-0.1: g1 R^2 + 1, g2 R^3/2 - 1 = ', g1 * r * r + 1.0_dp, g2 * r**3 / 2.0_dp - 1.0_dp
        call check(abs(g1 * r * r + 1.0_dp) < 1.0e-5_dp .and. abs(g2 * r**3 / 2.0_dp - 1.0_dp) < 1.0e-4_dp, &
                   'dipole: T R^3 -> point dipole')
    end subroutine t_dipole_radial

    subroutine t_dipole_ewald()
        real(dp) :: a1(3), a2(3), nh(3), d(3), g1(3), h1(3, 3), g2(3), h2(3, 3), rv(3), r1(3), r2(3), x(3)
        real(dp) :: gs(3), gh(3), hs(3, 3), hh(3, 3), e, r
        integer :: n1, n2, nm
        a1 = [4.7672460267921517_dp, 0.0_dp, 0.0_dp]
        a2 = [0.0_dp, 0.0_dp, 5.9335974928231314_dp]
        nh = [0.0_dp, -1.0_dp, 0.0_dp]
        d = [1.3_dp, -2.2_dp, 0.7_dp]
        call ewald2d_derivs(d, a1, a2, nh, 0.3_dp, .false., g1, h1)
        call ewald2d_derivs(d, a1, a2, nh, 0.7_dp, .false., g2, h2)
        e = max(maxval(abs(g1 - g2)), maxval(abs(h1 - h2)))
        print '(A,ES10.3,A,ES10.3)', '    pair: alpha 0.3 vs 0.7 ', e, ' ; trace(hess) ', h1(1, 1) + h1(2, 2) + h1(3, 3)
        call check(e < 1.0e-8_dp, 'dipole ewald: grad and hessian independent of alpha (pair)')
        call check(abs(h1(1, 1) + h1(2, 2) + h1(3, 3)) < 1.0e-7_dp, 'dipole ewald: Laplace (pair)')
        call ewald2d_derivs([0.0_dp, 0.0_dp, 0.0_dp], a1, a2, nh, 0.3_dp, .true., g1, h1)
        call ewald2d_derivs([0.0_dp, 0.0_dp, 0.0_dp], a1, a2, nh, 0.7_dp, .true., g2, h2)
        e = max(maxval(abs(g1 - g2)), maxval(abs(h1 - h2)))
        print '(A,ES10.3,A,3F10.6)', '    self: alpha 0.3 vs 0.7 ', e, ' ; grad (own sheet, reference) ', g1
        call check(e < 1.0e-8_dp, 'dipole ewald: independent of alpha (self)')
        call check(maxval(abs(g1 + 2.0_dp * PI / (a1(1) * a2(3)) * nh)) < 1.0e-8_dp, &
                   'dipole ewald: self gradient = -(2 pi/A) n (bottom reference)')
        ! soma direta: campo/hessiana de um par neutro no plano e hessiana do self (sum_{L/=0})
        r1 = [0.3_dp, 0.0_dp, 0.2_dp]
        r2 = [2.1_dp, 0.0_dp, 3.0_dp]
        x = [1.0_dp, 1.7_dp, 0.4_dp]
        call ewald2d_derivs(x - r1, a1, a2, nh, 0.5_dp, .false., g1, h1)
        call ewald2d_derivs(x - r2, a1, a2, nh, 0.5_dp, .false., g2, h2)
        g1 = g1 - g2
        h1 = h1 - h2
        nm = 800
        gs = 0.0_dp; gh = 0.0_dp; hs = 0.0_dp; hh = 0.0_dp
        do n1 = -nm, nm
            do n2 = -nm, nm
                rv = real(n1, dp) * a1 + real(n2, dp) * a2
                r = norm2(rv)
                if (r > 3600.0_dp) cycle
                call add_point(x - r1 + rv, 1.0_dp, r <= 1800.0_dp, gs, gh, hs, hh)
                call add_point(x - r2 + rv, -1.0_dp, r <= 1800.0_dp, gs, gh, hs, hh)
            end do
        end do
        gs = 2.0_dp * gs - gh
        hs = 2.0_dp * hs - hh
        e = max(maxval(abs(g1 - gs)), maxval(abs(h1 - hs)))
        print '(A,ES10.3)', '    neutral in-plane pair: |ewald - direct| (grad, hess) = ', e
        call check(e < 1.0e-6_dp, 'dipole ewald: charge-dipole and dipole-dipole kernels = direct sum')
        call ewald2d_derivs([0.0_dp, 0.0_dp, 0.0_dp], a1, a2, nh, 0.5_dp, .true., g1, h1)
        gs = 0.0_dp; gh = 0.0_dp; hs = 0.0_dp; hh = 0.0_dp
        do n1 = -nm, nm
            do n2 = -nm, nm
                if (n1 == 0 .and. n2 == 0) cycle
                rv = real(n1, dp) * a1 + real(n2, dp) * a2
                r = norm2(rv)
                if (r > 3600.0_dp) cycle
                call add_point(rv, 1.0_dp, r <= 1800.0_dp, gs, gh, hs, hh)
            end do
        end do
        hs = 2.0_dp * hs - hh
        e = maxval(abs(h1 - hs))
        print '(A,ES10.3,A,ES10.3)', '    self hessian: |ewald - direct| = ', e, ' ; trace ', h1(1, 1) + h1(2, 2) + h1(3, 3)
        call check(e < 1.0e-6_dp, 'dipole ewald: self dipole-dipole lattice sum = direct sum')
    end subroutine t_dipole_ewald

    subroutine t_dipole_integrals()
        type(pao_species_t) :: sp
        type(pao_dip_cache_t) :: dc
        integer :: ierr, no, i, j, c, mu, nu, k
        character(len=:), allocatable :: msg
        real(dp), allocatable :: s0(:, :), x0(:, :, :), s1(:, :), x1(:, :, :), s2(:, :), x2(:, :, :), dh(:, :), dh2(:, :)
        real(dp), allocatable :: sf(:, :), xf(:, :, :, :), rho(:, :)
        real(dp) :: d(3), e, gam(2, 2), gk(3, 2, 2), tk(3, 3, 2, 2), v(2), f(3, 2), eps, ep, em, fd, an, emax
        logical :: ok
        call read_pao_species('test/data/pao/Te7.0.pao', 's3p3d2f1', 7.0_dp, sp, ierr, msg)
        if (ierr /= 0) then
            call check(.false., 'dipole: read pao'); return
        end if
        no = sp%norb
        allocate(s0(no, no), x0(no, no, 3), s1(no, no), x1(no, no, 3), s2(no, no), x2(no, no, 3), dh(no, no), &
                 dh2(no, no))
        dc%h = 0.2_dp
        d = [2.1_dp, -3.4_dp, 1.7_dp]       ! Bohr
        call pao_dip_ensure(dc, sp, reshape([0.0_dp, 0.0_dp, 0.0_dp, d], [3, 2]))
        call pao_dip_get(dc, [0.0_dp, 0.0_dp, 0.0_dp], s0, x0, ok)
        e = 0.0_dp
        do i = 1, no
            e = max(e, maxval(abs(x0(i, i, :))))
        end do
        print '(A,ES10.3,A,ES10.3,A,ES10.3)', '    d = 0: max|S - 1| = ', maxval(abs(s0 - eye(no))), &
            ' ; max|X - X^T| = ', maxval(abs(x0(:, :, 1) - transpose(x0(:, :, 1)))), ' ; max|X_mumu| = ', e
        call check(maxval(abs(s0 - eye(no))) < 2.0e-3_dp, 'dipole integrals: S(0) = 1 (orthonormal PAO)')
        call check(maxval(abs(x0(:, :, 2) - transpose(x0(:, :, 2)))) == 0.0_dp .and. e < 1.0e-12_dp, &
                   'dipole integrals: X(0) symmetric, <mu|r - R|mu> = 0')
        call check(maxval(abs(x0)) > 0.1_dp, 'dipole integrals: on-site s-p, p-d couplings present')
        ! sentido oposto: derivado do cache x calculado direto na grade
        call pao_dip_get(dc, -d, s1, x1, ok)
        call pao_two_center(sp, -d, dc%h, s2, x2)
        e = max(maxval(abs(s1 - s2)), maxval(abs(x1 - x2)))
        print '(A,ES10.3)', '    X(-d) from cache vs direct grid: ', e
        call check(ok .and. e < 1.0e-3_dp, 'dipole integrals: X(-d)_ab = X(d)_ba + (-d) S_ba')
        ! termo de H: hermiticidade exata
        call pao_dip_get(dc, d, s1, x1, ok)
        call pao_dipole_dh([0.1_dp, -0.2_dp, 0.05_dp], [-0.03_dp, 0.07_dp, 0.11_dp], x1, s1, d, dh)
        call pao_dip_get(dc, -d, s2, x2, ok)
        call pao_dipole_dh([-0.03_dp, 0.07_dp, 0.11_dp], [0.1_dp, -0.2_dp, 0.05_dp], x2, s2, -d, dh2)
        print '(A,ES10.3)', '    dH(i->j) - dH(j->i)^T = ', maxval(abs(dh - transpose(dh2)))
        call check(maxval(abs(dh - transpose(dh2))) < 1.0e-14_dp, 'dipole H term: hermitian')
        ! dE/drho num modelo de dois sítios: S, X^(1), X^(2) completos (2 no x 2 no), rho simétrica
        allocate(sf(2 * no, 2 * no), xf(2 * no, 2 * no, 3, 2), rho(2 * no, 2 * no))
        sf(1:no, 1:no) = s0; sf(no + 1:, no + 1:) = s0
        sf(1:no, no + 1:) = s1; sf(no + 1:, 1:no) = transpose(s1)
        do c = 1, 3
            xf(1:no, 1:no, c, 1) = x0(:, :, c)
            xf(1:no, no + 1:, c, 1) = x1(:, :, c)
            xf(no + 1:, 1:no, c, 1) = transpose(x1(:, :, c))
            xf(no + 1:, no + 1:, c, 1) = x0(:, :, c) + d(c) * s0
            xf(:, :, c, 2) = xf(:, :, c, 1) - d(c) * sf
        end do
        do j = 1, 2 * no
            do i = 1, j
                rho(i, j) = 0.3_dp * sin(1.7_dp * i + 0.37_dp * j * j)
                rho(j, i) = rho(i, j)
            end do
        end do
        gam = reshape([0.4_dp, 0.12_dp, 0.12_dp, 0.35_dp], [2, 2])
        gk = 0.0_dp
        gk(:, 1, 2) = [0.02_dp, -0.05_dp, 0.03_dp]
        gk(:, 2, 1) = -gk(:, 1, 2)
        tk = 0.0_dp
        do i = 1, 2
            do c = 1, 3
                tk(c, c, i, i) = 0.3_dp
            end do
        end do
        tk(:, :, 1, 2) = reshape([0.01_dp, 0.004_dp, -0.02_dp, 0.004_dp, 0.03_dp, 0.0_dp, -0.02_dp, 0.0_dp, -0.04_dp], &
                                 [3, 3])
        tk(:, :, 2, 1) = transpose(tk(:, :, 1, 2))
        eps = 1.0e-5_dp
        emax = 0.0_dp
        do k = 1, 5
            mu = 1 + mod(7 * k, no)
            nu = no + 1 + mod(11 * k + 3, no)
            rho(mu, nu) = rho(mu, nu) + eps; rho(nu, mu) = rho(mu, nu)
            call model_vf(rho, no, sf, xf, gam, gk, tk, v, f, ep)
            rho(mu, nu) = rho(mu, nu) - 2.0_dp * eps; rho(nu, mu) = rho(mu, nu)
            call model_vf(rho, no, sf, xf, gam, gk, tk, v, f, em)
            rho(mu, nu) = rho(mu, nu) + eps; rho(nu, mu) = rho(mu, nu)
            fd = (ep - em) / (2.0_dp * eps)
            call model_vf(rho, no, sf, xf, gam, gk, tk, v, f, e)
            call pao_dipole_dh(f(:, 1), f(:, 2), x1, s1, d, dh)
            an = 2.0_dp * (0.5_dp * (v(1) + v(2)) * s1(mu, nu - no) + dh(mu, nu - no))
            emax = max(emax, abs(fd - an) / max(abs(fd), 1.0e-3_dp))
        end do
        print '(A,ES10.3)', '    dE/drho vs 2 dH (max rel) = ', emax
        call check(emax < 1.0e-6_dp, 'dipole H term = dE/drho with V = dE/dq, F = -dE/dp (sign +)')
    end subroutine t_dipole_integrals

    subroutine t_ewald_sym()
        real(dp) :: a1(3), a2(3), nh(3), d(3), e1, e2, area, dm(3), g1(3), h1(3, 3), g2(3), h2(3, 3)
        a1 = [4.7672460267921517_dp, 0.0_dp, 0.0_dp]
        a2 = [0.0_dp, 0.0_dp, 5.9335974928231314_dp]
        nh = [0.0_dp, -1.0_dp, 0.0_dp]
        area = a1(1) * a2(3)
        d = [1.3_dp, -2.2_dp, 0.7_dp]
        e1 = ewald2d_ref(d, a1, a2, nh, 0.3_dp, .false., .true.)
        e2 = ewald2d_ref(d, a1, a2, nh, 0.7_dp, .false., .true.)
        call check(abs(e1 - e2) < 1.0e-12_dp, 'ewald sym: independent of alpha (pair)')
        e2 = ewald2d_ref(d, a1, a2, nh, 0.3_dp, .false.)
        call check(abs(e1 - e2 - 2.0_dp * PI / area * dot_product(d, nh)) < 1.0e-12_dp, &
                   'ewald sym: = bottom reference + (2 pi/A) z')
        call check(ewald2d_ref(d, a1, a2, nh, 0.3_dp, .false., .false.) == e2, &
                   'ewald: sym = .false. is bit-identical to the default')
        dm = d - 2.0_dp * dot_product(d, nh) * nh
        e2 = ewald2d_ref(dm, a1, a2, nh, 0.3_dp, .false., .true.)
        call check(abs(e1 - e2) < 1.0e-12_dp, 'ewald sym: even in z')
        e1 = ewald2d_ref([0.5_dp, 0.0_dp, 0.0_dp] + 40.0_dp * nh, a1, a2, nh, 0.5_dp, .false., .true.)
        e2 = ewald2d_ref([0.5_dp, 0.0_dp, 0.0_dp] - 40.0_dp * nh, a1, a2, nh, 0.5_dp, .false., .true.)
        call check(abs(e1 + 2.0_dp * PI / area * 40.0_dp) < 1.0e-10_dp .and. abs(e1 - e2) < 1.0e-12_dp, &
                   'ewald sym: -(2 pi/A) |z| on both sides')
        call ewald2d_derivs([0.0_dp, 0.0_dp, 0.0_dp], a1, a2, nh, 0.3_dp, .true., g1, h1, .true.)
        call ewald2d_derivs([0.0_dp, 0.0_dp, 0.0_dp], a1, a2, nh, 0.7_dp, .true., g2, h2, .true.)
        call check(maxval(abs(g1)) < 1.0e-9_dp .and. maxval(abs(h1 - h2)) < 1.0e-8_dp, &
                   'ewald sym: self gradient zero, hessian independent of alpha')
    end subroutine t_ewald_sym

    !> Campo pontual (scc_point_field) com densidade de esfera uniforme (raio a): v' = -r/a^3, v'' = -1/a^3
    !> dentro, 1/r fora; termo próprio -p/a^3. Átomos no mesmo plano, par neutro, dipolos gerais:
    !> soma direta num disco (Richardson R, 2R) contra Ewald 2D simétrico e com referência no fundo.
    subroutine t_point_field()
        real(dp), parameter :: HA = 27.211386245988_dp, BA = 0.529177210903_dp
        type(radial_density_t) :: rd
        integer :: ierr, n, i, j, n1, n2, nm, a
        character(len=:), allocatable :: msg
        real(dp) :: aa, a1(3), a2(3), nh(3), pos(3, 2), qn(2), p(3, 2), ef(3, 2), es(3, 2), ef0(3, 2), es0(3, 2)
        real(dp) :: ed(3, 2), eh(3, 2), rv(3), d(3), rb, u(3), v1, v2, hm(3, 3), t(3), e, alpha
        aa = 3.0_dp
        n = 6000
        rd%n = n; rd%dr = aa / n; rd%rc = aa
        allocate(rd%rho(0:n))
        rd%rho = 1.0_dp
        call finish_radial_density(rd, ierr, msg)
        a1 = [4.7672460267921517_dp, 0.0_dp, 0.0_dp]
        a2 = [0.0_dp, 0.0_dp, 5.9335974928231314_dp]
        nh = [0.0_dp, -1.0_dp, 0.0_dp]
        alpha = sqrt(PI / (a1(1) * a2(3)))
        pos(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        pos(:, 2) = [1.0_dp, 0.0_dp, 0.6_dp]          ! 1,17 Å < a = 1,59 Å: curto alcance dentro da esfera
        qn = [0.3_dp, -0.3_dp]
        p(:, 1) = [0.2_dp, -0.1_dp, 0.15_dp]
        p(:, 2) = [-0.05_dp, 0.3_dp, 0.1_dp]
        call scc_point_field(rd, pos, a1, a2, nh, alpha, qn, p, ef, es, .true.)
        call scc_point_field(rd, pos, a1, a2, nh, alpha, qn, p, ef0, es0)
        call check(maxval(abs(es + p / aa**3 * HA / BA)) < 1.0e-6_dp, 'point field: self term -(4 pi/3) rho(0) p')
        ! soma direta (Hartree/Bohr), disco de raio R e R/2
        nm = 800
        ed = 0.0_dp
        eh = 0.0_dp
        do i = 1, 2
            do j = 1, 2
                do n1 = -nm, nm
                    do n2 = -nm, nm
                        rv = real(n1, dp) * a1 + real(n2, dp) * a2
                        if (norm2(rv) > 3600.0_dp) cycle
                        d = (pos(:, i) - pos(:, j) - rv) / BA
                        rb = norm2(d)
                        if (rb < 1.0e-8_dp) then
                            t = -p(:, i) / aa**3
                        else
                            u = d / rb
                            if (rb < aa) then
                                v1 = -rb / aa**3; v2 = -1.0_dp / aa**3
                            else
                                v1 = -1.0_dp / rb**2; v2 = 2.0_dp / rb**3
                            end if
                            do a = 1, 3
                                hm(:, a) = (v2 - v1 / rb) * u * u(a)
                                hm(a, a) = hm(a, a) + v1 / rb
                            end do
                            t = qn(j) * v1 * u + matmul(hm, p(:, j))
                        end if
                        ed(:, i) = ed(:, i) + t
                        if (norm2(rv) <= 1800.0_dp) eh(:, i) = eh(:, i) + t
                    end do
                end do
            end do
        end do
        ed = (2.0_dp * ed - eh) * HA / BA
        e = maxval(abs(ef - ed))
        print '(A,6F12.6)', '    point field (V/A), Ewald sym: ', ef
        print '(A,6F12.6)', '    point field (V/A), direct   : ', ed
        print '(A,ES10.3,A,ES10.3)', '    |ewald - direct| = ', e, ' ; |sym - bottom| (neutral) = ', maxval(abs(ef - ef0))
        call check(e < 1.0e-5_dp, 'point field: Ewald 2D + short range = direct sum')
        call check(maxval(abs(ef - ef0)) < 1.0e-8_dp, 'point field: neutral cell, same field for both G = 0 references')
    end subroutine t_point_field

    subroutine add_point(v, qq, inner, gs, gh, hs, hh)
        real(dp), intent(in) :: v(3), qq
        logical, intent(in) :: inner
        real(dp), intent(inout) :: gs(3), gh(3), hs(3, 3), hh(3, 3)
        real(dp) :: rr, gg(3), hq(3, 3)
        integer :: a, b
        rr = norm2(v)
        gg = -qq * v / rr**3
        do b = 1, 3
            do a = 1, 3
                hq(a, b) = qq * 3.0_dp * v(a) * v(b) / rr**5
            end do
            hq(b, b) = hq(b, b) - qq / rr**3
        end do
        gs = gs + gg; hs = hs + hq
        if (inner) then
            gh = gh + gg; hh = hh + hq
        end if
    end subroutine add_point

    function eye(n) result(a)
        integer, intent(in) :: n
        real(dp) :: a(n, n)
        integer :: ii
        a = 0.0_dp
        do ii = 1, n
            a(ii, ii) = 1.0_dp
        end do
    end function eye

    !> Modelo de dois sítios: q, p de rho; E, V = dE/dq, F = -dE/dp (núcleos simétricos).
    !> A10: dois átomos (on-site q, p) e um par (q_b no ponto médio). E = 1/2 q G q + q C p + 1/2 p T p
    !> (sítios 1, 2 = átomos, 3 = par; C(:, s, a) carga-dipolo, F_a = -dE/dp_a); termo de H do modo 'bond'
    !> conferido contra dE/drho por diferenças finitas (perturbação simétrica: dE = 2 dH drho).
    subroutine t_bond_model()
        type(pao_species_t) :: sp
        type(pao_dip_cache_t) :: dc
        integer :: ierr, no, i, j, k, mu, nu, kind
        character(len=:), allocatable :: msg
        real(dp), allocatable :: s0(:, :), x0(:, :, :), s1(:, :), x1(:, :, :), dh(:, :), rho(:, :)
        real(dp) :: d(3), gm(3, 3), cg(3, 3, 2), tk(3, 3, 2, 2), v(3), f(3, 2), e, ep, em, fd, an, emax, eps
        logical :: ok
        call read_pao_species('test/data/pao/Te7.0.pao', 's3p3d2f1', 7.0_dp, sp, ierr, msg)
        if (ierr /= 0) then
            call check(.false., 'bond: read pao'); return
        end if
        no = sp%norb
        allocate(s0(no, no), x0(no, no, 3), s1(no, no), x1(no, no, 3), dh(no, no), rho(2 * no, 2 * no))
        dc%h = 0.2_dp
        d = [2.1_dp, -3.4_dp, 1.7_dp]
        call pao_dip_ensure(dc, sp, reshape([0.0_dp, 0.0_dp, 0.0_dp, d], [3, 2]))
        call pao_dip_get(dc, [0.0_dp, 0.0_dp, 0.0_dp], s0, x0, ok)
        call pao_dip_get(dc, d, s1, x1, ok)
        do j = 1, 2 * no
            do i = 1, j
                rho(i, j) = 0.3_dp * sin(1.7_dp * i + 0.37_dp * j * j)
                rho(j, i) = rho(i, j)
            end do
        end do
        gm = reshape([0.40_dp, 0.12_dp, 0.20_dp, 0.12_dp, 0.35_dp, 0.18_dp, 0.20_dp, 0.18_dp, 0.38_dp], [3, 3])
        do k = 1, 3
            cg(:, k, 1) = [0.02_dp, -0.05_dp, 0.03_dp] * real(k, dp)
            cg(:, k, 2) = [-0.04_dp, 0.01_dp, 0.06_dp] / real(k, dp)
        end do
        tk = 0.0_dp
        do i = 1, 2
            do k = 1, 3
                tk(k, k, i, i) = 0.3_dp
            end do
        end do
        tk(:, :, 1, 2) = reshape([0.01_dp, 0.004_dp, -0.02_dp, 0.004_dp, 0.03_dp, 0.0_dp, -0.02_dp, 0.0_dp, -0.04_dp], &
                                 [3, 3])
        tk(:, :, 2, 1) = transpose(tk(:, :, 1, 2))
        eps = 1.0e-5_dp
        emax = 0.0_dp
        do kind = 1, 3
            do k = 1, 4
                select case (kind)
                case (1)                                       ! on-site do átomo 1
                    mu = 1 + mod(7 * k, no); nu = 1 + mod(11 * k + 3, no)
                case (2)                                       ! on-site do átomo 2
                    mu = no + 1 + mod(5 * k, no); nu = no + 1 + mod(13 * k + 1, no)
                case default                                   ! par 1-2
                    mu = 1 + mod(7 * k, no); nu = no + 1 + mod(11 * k + 3, no)
                end select
                if (mu == nu) cycle
                rho(mu, nu) = rho(mu, nu) + eps; rho(nu, mu) = rho(mu, nu)
                call model_bond(rho, no, s0, x0, s1, gm, cg, tk, v, f, ep)
                rho(mu, nu) = rho(mu, nu) - 2.0_dp * eps; rho(nu, mu) = rho(mu, nu)
                call model_bond(rho, no, s0, x0, s1, gm, cg, tk, v, f, em)
                rho(mu, nu) = rho(mu, nu) + eps; rho(nu, mu) = rho(mu, nu)
                fd = (ep - em) / (2.0_dp * eps)
                call model_bond(rho, no, s0, x0, s1, gm, cg, tk, v, f, e)
                select case (kind)
                case (1)
                    call pao_dipole_dh(f(:, 1), f(:, 1), x0, s0, [0.0_dp, 0.0_dp, 0.0_dp], dh)
                    an = 2.0_dp * (v(1) * s0(mu, nu) + dh(mu, nu))
                case (2)
                    call pao_dipole_dh(f(:, 2), f(:, 2), x0, s0, [0.0_dp, 0.0_dp, 0.0_dp], dh)
                    an = 2.0_dp * (v(2) * s0(mu - no, nu - no) + dh(mu - no, nu - no))
                case default
                    an = 2.0_dp * v(3) * s1(mu, nu - no)
                end select
                emax = max(emax, abs(fd - an) / max(abs(fd), 1.0e-3_dp))
            end do
        end do
        print '(A,ES10.3)', '    bond mode: dE/drho vs 2 dH (on-site V S + F.X, bond V_b S), max rel = ', emax
        call check(emax < 1.0e-6_dp, 'bond H term = dE/drho (on-site charge and dipole, bond charge)')
    end subroutine t_bond_model

    subroutine model_bond(rh, no, s0, x0, s1, gm, cg, tk, vv, ff, en)
        real(dp), intent(in) :: rh(:, :), s0(:, :), x0(:, :, :), s1(:, :), gm(3, 3), cg(3, 3, 2), tk(3, 3, 2, 2)
        integer, intent(in) :: no
        real(dp), intent(out) :: vv(3), ff(3, 2), en
        real(dp) :: qq(3), pp(3, 2)
        integer :: ia, o, cc, s, b
        do ia = 1, 2
            o = (ia - 1) * no
            qq(ia) = sum(rh(o + 1:o + no, o + 1:o + no) * transpose(s0))
            do cc = 1, 3
                pp(cc, ia) = -sum(rh(o + 1:o + no, o + 1:o + no) * transpose(x0(:, :, cc)))
            end do
        end do
        qq(3) = sum(rh(1:no, no + 1:) * s1) + sum(rh(no + 1:, 1:no) * transpose(s1))
        en = 0.5_dp * dot_product(qq, matmul(gm, qq))
        vv = matmul(gm, qq)
        ff = 0.0_dp
        do s = 1, 3
            do ia = 1, 2
                en = en + qq(s) * dot_product(cg(:, s, ia), pp(:, ia))
                vv(s) = vv(s) + dot_product(cg(:, s, ia), pp(:, ia))
                ff(:, ia) = ff(:, ia) - cg(:, s, ia) * qq(s)
            end do
        end do
        do ia = 1, 2
            do b = 1, 2
                en = en + 0.5_dp * dot_product(pp(:, ia), matmul(tk(:, :, ia, b), pp(:, b)))
                ff(:, ia) = ff(:, ia) - matmul(tk(:, :, ia, b), pp(:, b))
            end do
        end do
    end subroutine model_bond

    !> A10: Gamma (scc_gamma_matrix) de átomos + pontos médios (dois sítios no mesmo ponto) e o núcleo cruzado
    !> carga-dipolo (scc_cross_kernel) contra a soma direta de gamma(|d + L|) num disco (Richardson R, R/2); o sítio 4 é o 3 transladado
    !> de a1 (dois sítios no mesmo ponto módulo a rede); sítios no mesmo plano, conjunto neutro, corte por célula.
    subroutine t_bond_ewald()
        real(dp), parameter :: HA = 27.211386245988_dp, BA = 0.529177210903_dp
        type(pao_species_t) :: sp
        type(radial_density_t) :: rd
        integer :: ierr, i, j, n1, n2, nm
        character(len=:), allocatable :: msg
        real(dp) :: a1(3), a2(3), nh(3), pos(3, 4), qn(4), gam(4, 4), uh, vd(4), vh(4), rv(3), r, g, e, alpha
        real(dp) :: gk(3, 4, 2), gd(3, 4), gdh(3, 4), g1, g2, gx(3), rl
        call read_pao_species('test/data/pao/Te7.0.pao', 's3p3d2f1', 7.0_dp, sp, ierr, msg)
        if (ierr == 0) call build_radial_density(sp, 1, 1, 5.0e-4_dp, rd, ierr, msg)
        if (ierr /= 0) then
            call check(.false., 'bond ewald: pao'); return
        end if
        a1 = [4.7672460267921517_dp, 0.0_dp, 0.0_dp]
        a2 = [0.0_dp, 0.0_dp, 5.9335974928231314_dp]
        nh = [0.0_dp, -1.0_dp, 0.0_dp]
        alpha = sqrt(PI / (a1(1) * a2(3)))
        pos(:, 1) = [0.3_dp, 0.0_dp, 0.4_dp]
        pos(:, 2) = [2.2_dp, 0.0_dp, 2.1_dp]
        pos(:, 3) = 0.5_dp * (pos(:, 1) + pos(:, 2))
        pos(:, 4) = pos(:, 3) + a1                                     ! mesmo ponto módulo a rede
        qn = [1.0_dp, 1.0_dp, -1.5_dp, -0.5_dp]
        call scc_gamma_matrix(rd, pos, a1, a2, nh, alpha, gam, uh)
        call scc_cross_kernel(rd, pos, pos(:, 1:2), a1, a2, nh, alpha, gk)
        vd = 0.0_dp; vh = 0.0_dp; gd = 0.0_dp; gdh = 0.0_dp
        nm = ceiling(3600.0_dp / a1(1)) + 1
        do i = 1, 4
            do j = 1, 4
                do n1 = -nm, nm
                    do n2 = -nm, nm
                        rl = norm2(real(n1, dp) * a1 + real(n2, dp) * a2)      ! corte por célula (neutra)
                        if (rl > 3600.0_dp) cycle
                        rv = pos(:, i) - pos(:, j) + real(n1, dp) * a1 + real(n2, dp) * a2
                        r = norm2(rv)
                        if (r < 1.0e-8_dp) then
                            g = uh
                        else
                            g = HA * gamma_au(rd, r / BA)
                        end if
                        vd(i) = vd(i) + qn(j) * g
                        if (rl <= 1800.0_dp) vh(i) = vh(i) + qn(j) * g
                        if (j <= 2 .and. r > 1.0e-8_dp) then        ! grad gamma (eV/Bohr), dipolo nos átomos
                            call gamma_derivs_au(rd, r / BA, g1, g2)
                            gx = HA * g1 * rv / r
                            gd(:, i) = gd(:, i) + gx * merge(1.0_dp, -1.0_dp, j == 1)
                            if (rl <= 1800.0_dp) gdh(:, i) = gdh(:, i) + gx * merge(1.0_dp, -1.0_dp, j == 1)
                        end if
                    end do
                end do
            end do
        end do
        vd = 2.0_dp * vd - vh
        gd = 2.0_dp * gd - gdh
        e = maxval(abs(matmul(gam, qn) - vd))
        print '(A,ES10.3,A,ES10.3)', '    bond sites: |Gamma q - direct| = ', e, ' eV ; Gamma(3,4) - Gamma(3,3) = ', &
            gam(3, 4) - gam(3, 3)
        call check(e < 1.0e-5_dp, 'bond ewald: Gamma of atoms + midpoints (coincident sites) = direct sum')
        call check(abs(gam(3, 4) - gam(3, 3)) < 1.0e-12_dp .and. abs(gam(3, 4) - gam(4, 3)) < 1.0e-12_dp, &
                   'bond ewald: two sites at one point couple with the on-site value')
        e = maxval(abs(gk(:, :, 1) - gk(:, :, 2) - gd))
        print '(A,ES10.3,A)', '    cross kernel (neutral dipole pair): |ewald - direct| = ', e, ' eV/Bohr'
        call check(e < 1.0e-5_dp, 'bond ewald: charge-dipole cross kernel = direct sum')
    end subroutine t_bond_ewald

    subroutine model_vf(rh, no, sf, xf, gam, gk, tk, vv, ff, en)
        real(dp), intent(in) :: rh(:, :), sf(:, :), xf(:, :, :, :), gam(2, 2), gk(3, 2, 2), tk(3, 3, 2, 2)
        integer, intent(in) :: no
        real(dp), intent(out) :: vv(2), ff(3, 2), en
        real(dp) :: qq(2), pp(3, 2)
        integer :: ia, o, cc, jj
        do ia = 1, 2
            o = (ia - 1) * no
            qq(ia) = sum(rh(o + 1:o + no, :) * transpose(sf(:, o + 1:o + no)))
            do cc = 1, 3
                pp(cc, ia) = -sum(rh(o + 1:o + no, :) * transpose(xf(:, o + 1:o + no, cc, ia)))
            end do
        end do
        en = 0.5_dp * dot_product(qq, matmul(gam, qq))
        do ia = 1, 2
            vv(ia) = dot_product(gam(ia, :), qq)
            ff(:, ia) = 0.0_dp
            do jj = 1, 2
                en = en + qq(ia) * dot_product(gk(:, ia, jj), pp(:, jj)) + &
                     0.5_dp * dot_product(pp(:, ia), matmul(tk(:, :, ia, jj), pp(:, jj)))
                vv(ia) = vv(ia) + dot_product(gk(:, ia, jj), pp(:, jj))
                ff(:, ia) = ff(:, ia) + gk(:, ia, jj) * qq(jj) - matmul(tk(:, :, ia, jj), pp(:, jj))
            end do
        end do
    end subroutine model_vf

end program test_scc
