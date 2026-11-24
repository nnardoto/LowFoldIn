! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Funções de base PAO do OpenMX em espaço real.
!>
!> Lê as radiais R_{L,mul}(r) de um arquivo `.pao` (blocos
!> `<pseudo.atomic.orbitals.L=l`, linhas `x r R_0 ... R_{Mul-1}`, r em Bohr)
!> e avalia phi(r) = R(|r|) Y(r^) exatamente como o Get_Orbitals.c do
!> OpenMX 4.0: interpolação cúbica de Hermite na malha do arquivo,
!> extrapolação polinomial para r < r_1, zero para r > r_último, e
!> harmônicos reais normalizados na ordem do OpenMX
!>   p: x, y, z
!>   d: 3z^2-r^2, x^2-y^2, xy, xz, yz
!>   f: z^3, xz^2, yz^2, z(x^2-y^2), xyz, x^3, y^3 (fatores do Get_Orbitals.c)
!> A ordem dos orbitais de um átomo é a do OpenMX: L crescente, para cada L
!> os mul = 0..n_L-1, para cada mul as componentes m acima. A base é dada
!> como no OpenMX ('s3p2d2' = 3 s, 2 p, 2 d), sem contração. `rcut` (Bohr)
!> limita a esfera do átomo como o GridListAtom do OpenMX (r <= rcut).
!>
!> Rotinas sem stop; erros por ierr/errmsg.
module pao_basis
    use types, only: dp
    implicit none
    private

    public :: pao_species_t, read_pao_species, pao_eval, pao_eval_grad
    public :: pao_dip_cache_t, pao_two_center, pao_dip_ensure, pao_dip_get, pao_dipole_dh
    public :: pao_planar_profile

    type :: pao_species_t
        integer :: nmesh = 0, lmax = -1, norb = 0
        integer :: nmul(0:3) = 0
        real(dp) :: rcut = 0.0_dp                    !< Bohr
        real(dp), allocatable :: rv(:)               !< (nmesh) Bohr
        real(dp), allocatable :: rwf(:, :, :)        !< (nmesh, maxmul, 0:lmax)
    end type pao_species_t

    !> Integrais de dois centros da mesma espécie (Bohr) por vetor d = R_b - R_a (A6, dipolos):
    !>   S_ab = <a|b(d)>,  X_ab(:, c) = <a| (r - R_a)_c |b(d)>.
    !> Guardadas só para d canônico (primeira componente não nula > 0, ou d = 0); o sentido
    !> oposto sai exato de S(-d)_ab = S(d)_ba e X(-d)_ab = X(d)_ba - d S(d)_ba, o que torna
    !> hermitianos, bit a bit, os blocos montados com elas.
    type :: pao_dip_cache_t
        integer :: n = 0
        real(dp) :: h = 0.1_dp                       !< passo da grade cúbica (Bohr)
        real(dp), allocatable :: d(:, :)             !< (3, n) Bohr, canônicos
        real(dp), allocatable :: s(:, :, :)          !< (norb, norb, n)
        real(dp), allocatable :: x(:, :, :, :)       !< (norb, norb, 3, n)
    end type pao_dip_cache_t

    real(dp), parameter :: DTOL = 1.0e-6_dp          !< Bohr: chaves iguais

contains

    !> Lê `file` e monta a espécie com a base `basis` ('s3p2d2' etc.).
    subroutine read_pao_species(file, basis, rcut, sp, ierr, errmsg)
        character(len=*), intent(in) :: file, basis
        real(dp), intent(in) :: rcut
        type(pao_species_t), intent(out) :: sp
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: u, ios, l, i, k, pmul, plmax, nm, maxmul
        character(len=512) :: line
        real(dp), allocatable :: row(:)
        character(len=1), parameter :: lch(0:3) = ['s', 'p', 'd', 'f']
        ierr = 0
        errmsg = ''
        sp%rcut = rcut
        ! base: letras s/p/d/f seguidas do número de radiais
        i = 1
        do while (i <= len_trim(basis))
            l = index('spdf', basis(i:i)) - 1
            if (l < 0) then
                ierr = 2; errmsg = 'pao: bad basis spec ''' // trim(basis) // ''''; return
            end if
            k = i + 1
            do while (k <= len_trim(basis))
                if (index('0123456789', basis(k:k)) == 0) exit
                k = k + 1
            end do
            if (k == i + 1) then
                ierr = 2; errmsg = 'pao: bad basis spec ''' // trim(basis) // ''''; return
            end if
            read(basis(i + 1:k - 1), *) sp%nmul(l)
            sp%lmax = max(sp%lmax, l)
            i = k
        end do
        sp%norb = 0
        do l = 0, 3
            sp%norb = sp%norb + sp%nmul(l) * (2 * l + 1)
        end do
        open(newunit=u, file=file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = 'pao: cannot open ' // trim(file); return
        end if
        nm = 0; pmul = 0; plmax = -1
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (line(1:15) == 'grid.num.output') read(line(16:), *) nm
            if (line(1:8) == 'PAO.Lmax') read(line(9:), *) plmax
            if (line(1:7) == 'PAO.Mul') read(line(8:), *) pmul
        end do
        if (nm < 4 .or. pmul < 1 .or. plmax < sp%lmax) then
            close(u); ierr = 3; errmsg = 'pao: header (grid.num.output/PAO.Lmax/PAO.Mul) missing in ' // trim(file)
            return
        end if
        maxmul = maxval(sp%nmul)
        if (maxmul > pmul) then
            close(u); ierr = 3; errmsg = 'pao: basis asks more radials than PAO.Mul in ' // trim(file); return
        end if
        sp%nmesh = nm
        allocate(sp%rv(nm), sp%rwf(nm, max(maxmul, 1), 0:sp%lmax), row(pmul + 2))
        sp%rwf = 0.0_dp
        do l = 0, sp%lmax
            rewind(u)
            do
                read(u, '(A)', iostat=ios) line
                if (ios /= 0) then
                    close(u); ierr = 3
                    errmsg = 'pao: block <pseudo.atomic.orbitals.L=' // lch(l) // ' missing'; return
                end if
                line = adjustl(line)
                if (line(1:28) == '<pseudo.atomic.orbitals.L=' // achar(48 + l) // ' ' .or. &
                    trim(line) == '<pseudo.atomic.orbitals.L=' // achar(48 + l)) exit
            end do
            do i = 1, nm
                read(u, *, iostat=ios) row
                if (ios /= 0) then
                    close(u); ierr = 3; errmsg = 'pao: short radial block in ' // trim(file); return
                end if
                sp%rv(i) = row(2)
                if (sp%nmul(l) > 0) sp%rwf(i, 1:sp%nmul(l), l) = row(3:2 + sp%nmul(l))
            end do
        end do
        close(u)
    end subroutine read_pao_species

    !> phi(1:norb) dos orbitais da espécie no deslocamento d = r - R_atomo (Bohr).
    pure subroutine pao_eval(sp, d, phi)
        type(pao_species_t), intent(in) :: sp
        real(dp), intent(in) :: d(3)
        real(dp), intent(out) :: phi(:)
        real(dp) :: rf(8, 0:3), af(7), r, q, p, siq, coq, sip, cop, dum, dum1, dum2
        real(dp) :: h1, h2, h3, x1, x2, y1, y2, y12, y22, g1, g2, f1, f2, f3, f4, f, df, rm, a, b, c, dd
        integer :: m, lo, hi, l, mu, o, nm, k
        phi = 0.0_dp
        r = norm2(d)
        nm = sp%nmesh
        if (r > sp%rcut .or. r > sp%rv(nm)) return
        if (r < 1.0e-13_dp) then           ! xyz2spherical do OpenMX: ângulos nulos na origem
            q = 0.0_dp; p = 0.0_dp
        else
            q = acos(max(-1.0_dp, min(1.0_dp, d(3) / r)))
            p = atan2(d(2), d(1))
        end if
        if (r < sp%rv(1)) then
            m = 5                           ! m = 4 em base 0
            rm = sp%rv(m)
        else
            lo = 1; hi = nm
            do while (hi - lo /= 1)
                m = (lo + hi) / 2
                if (sp%rv(m) < r) then
                    lo = m
                else
                    hi = m
                end if
            end do
            m = hi
            rm = r
        end if
        h1 = sp%rv(m - 1) - sp%rv(max(m - 2, 1))
        h2 = sp%rv(m) - sp%rv(m - 1)
        h3 = sp%rv(min(m + 1, nm)) - sp%rv(m)
        x1 = rm - sp%rv(m - 1)
        x2 = rm - sp%rv(m)
        y1 = x1 / h2; y2 = x2 / h2; y12 = y1 * y1; y22 = y2 * y2
        do l = 0, sp%lmax
            do mu = 1, sp%nmul(l)
                f2 = sp%rwf(m - 1, mu, l); f3 = sp%rwf(m, mu, l)
                if (m == 2) then                ! m == 1 em base 0
                    h1 = -(h2 + h3)
                    f4 = sp%rwf(m + 1, mu, l); f1 = f4
                else if (m == nm) then          ! m == nmesh-1 em base 0
                    f1 = sp%rwf(m - 2, mu, l)
                    h3 = -(h1 + h2); f4 = f1
                else
                    f1 = sp%rwf(m - 2, mu, l); f4 = sp%rwf(m + 1, mu, l)
                end if
                dum = h1 + h2; dum1 = h1 / h2 / dum; dum2 = h2 / h1 / dum
                g1 = (f3 - f2) * dum1 + (f2 - f1) * dum2
                dum = h2 + h3; dum1 = h2 / h3 / dum; dum2 = h3 / h2 / dum
                g2 = (f4 - f3) * dum1 + (f3 - f2) * dum2
                f = y22 * (3.0_dp * f2 + h2 * g1 + (2.0_dp * f2 + h2 * g1) * y2) &
                  + y12 * (3.0_dp * f3 - h2 * g2 - (2.0_dp * f3 - h2 * g2) * y1)
                if (r < sp%rv(1)) then
                    df = 2.0_dp * y2 / h2 * (3.0_dp * f2 + h2 * g1 + (2.0_dp * f2 + h2 * g1) * y2) &
                       + y22 * (2.0_dp * f2 + h2 * g1) / h2 &
                       + 2.0_dp * y1 / h2 * (3.0_dp * f3 - h2 * g2 - (2.0_dp * f3 - h2 * g2) * y1) &
                       - y12 * (2.0_dp * f3 - h2 * g2) / h2
                    if (l == 0) then
                        a = 0.0_dp; b = 0.5_dp * df / rm; c = 0.0_dp; dd = f - b * rm * rm
                    else if (l == 1) then
                        a = (rm * df - f) / (2.0_dp * rm**3); b = 0.0_dp; c = df - 3.0_dp * a * rm * rm; dd = 0.0_dp
                    else
                        b = (3.0_dp * f - rm * df) / (rm * rm); a = (f - b * rm * rm) / rm**3; c = 0.0_dp; dd = 0.0_dp
                    end if
                    f = a * r**3 + b * r * r + c * r + dd
                end if
                rf(mu, l) = f
            end do
        end do
        siq = sin(q); coq = cos(q); sip = sin(p); cop = cos(p)
        o = 0
        do l = 0, sp%lmax
            select case (l)
            case (0)
                af(1) = 0.282094791773878_dp
            case (1)
                dum = 0.48860251190292_dp * siq
                af(1) = dum * cop; af(2) = dum * sip; af(3) = 0.48860251190292_dp * coq
            case (2)
                dum1 = siq * siq
                dum2 = 1.09254843059208_dp * siq * coq
                af(1) = 0.94617469575756_dp * coq * coq - 0.31539156525252_dp
                af(2) = 0.54627421529604_dp * dum1 * (1.0_dp - 2.0_dp * sip * sip)
                af(3) = 1.09254843059208_dp * dum1 * sip * cop
                af(4) = dum2 * cop
                af(5) = dum2 * sip
            case (3)
                af(1) = 0.373176332590116_dp * (5.0_dp * coq**3 - 3.0_dp * coq)
                af(2) = 0.457045799464466_dp * cop * siq * (5.0_dp * coq * coq - 1.0_dp)
                af(3) = 0.457045799464466_dp * sip * siq * (5.0_dp * coq * coq - 1.0_dp)
                af(4) = 1.44530572132028_dp * siq * siq * coq * (cop * cop - sip * sip)
                af(5) = 2.89061144264055_dp * siq * siq * coq * sip * cop
                af(6) = 0.590043589926644_dp * siq**3 * (4.0_dp * cop**3 - 3.0_dp * cop)
                af(7) = 0.590043589926644_dp * siq**3 * (3.0_dp * sip - 4.0_dp * sip**3)
            end select
            do mu = 1, sp%nmul(l)
                do k = 1, 2 * l + 1
                    o = o + 1
                    phi(o) = rf(mu, l) * af(k)
                end do
            end do
        end do
    end subroutine pao_eval

    !> phi(1:norb) e grad phi (3, norb) (Bohr^-1) no deslocamento d = r - R_atomo (Bohr) (A9, ELF).
    !> Mesma radial que pao_eval (Hermite do Get_Orbitals.c e extrapolação para r < r_1), com a
    !> derivada analítica R'(r) do mesmo polinômio; harmônicos reais escritos como polinômios
    !> P(u) no versor u = d/r (mesmos fatores do pao_eval), grad phi = R' u P + (R/r)(grad_u P - u (u.grad_u P)).
    !> Em r -> 0, R/r -> R'(0) (só l = 1 sobrevive). phi coincide com pao_eval a menos de arredondamento.
    pure subroutine pao_eval_grad(sp, d, phi, gphi)
        type(pao_species_t), intent(in) :: sp
        real(dp), intent(in) :: d(3)
        real(dp), intent(out) :: phi(:), gphi(:, :)
        real(dp) :: rf(8, 0:3), drf(8, 0:3), af(7), gaf(3, 7), u(3), r, x, y, z, ror, ud, tg(3)
        real(dp) :: dum, dum1, dum2
        real(dp) :: h1, h2, h3, x1, x2, y1, y2, y12, y22, g1, g2, f1, f2, f3, f4, f, df, rm, a, b, c, dd
        integer :: m, lo, hi, l, mu, o, nm, k
        phi = 0.0_dp
        gphi = 0.0_dp
        r = norm2(d)
        nm = sp%nmesh
        if (r > sp%rcut .or. r > sp%rv(nm)) return
        if (r < 1.0e-13_dp) then
            u = [0.0_dp, 0.0_dp, 1.0_dp]
        else
            u = d / r
        end if
        if (r < sp%rv(1)) then
            m = 5
            rm = sp%rv(m)
        else
            lo = 1; hi = nm
            do while (hi - lo /= 1)
                m = (lo + hi) / 2
                if (sp%rv(m) < r) then
                    lo = m
                else
                    hi = m
                end if
            end do
            m = hi
            rm = r
        end if
        h1 = sp%rv(m - 1) - sp%rv(max(m - 2, 1))
        h2 = sp%rv(m) - sp%rv(m - 1)
        h3 = sp%rv(min(m + 1, nm)) - sp%rv(m)
        x1 = rm - sp%rv(m - 1)
        x2 = rm - sp%rv(m)
        y1 = x1 / h2; y2 = x2 / h2; y12 = y1 * y1; y22 = y2 * y2
        do l = 0, sp%lmax
            do mu = 1, sp%nmul(l)
                f2 = sp%rwf(m - 1, mu, l); f3 = sp%rwf(m, mu, l)
                if (m == 2) then
                    h1 = -(h2 + h3)
                    f4 = sp%rwf(m + 1, mu, l); f1 = f4
                else if (m == nm) then
                    f1 = sp%rwf(m - 2, mu, l)
                    h3 = -(h1 + h2); f4 = f1
                else
                    f1 = sp%rwf(m - 2, mu, l); f4 = sp%rwf(m + 1, mu, l)
                end if
                dum = h1 + h2; dum1 = h1 / h2 / dum; dum2 = h2 / h1 / dum
                g1 = (f3 - f2) * dum1 + (f2 - f1) * dum2
                dum = h2 + h3; dum1 = h2 / h3 / dum; dum2 = h3 / h2 / dum
                g2 = (f4 - f3) * dum1 + (f3 - f2) * dum2
                f = y22 * (3.0_dp * f2 + h2 * g1 + (2.0_dp * f2 + h2 * g1) * y2) &
                  + y12 * (3.0_dp * f3 - h2 * g2 - (2.0_dp * f3 - h2 * g2) * y1)
                df = 2.0_dp * y2 / h2 * (3.0_dp * f2 + h2 * g1 + (2.0_dp * f2 + h2 * g1) * y2) &
                   + y22 * (2.0_dp * f2 + h2 * g1) / h2 &
                   + 2.0_dp * y1 / h2 * (3.0_dp * f3 - h2 * g2 - (2.0_dp * f3 - h2 * g2) * y1) &
                   - y12 * (2.0_dp * f3 - h2 * g2) / h2
                if (r < sp%rv(1)) then
                    if (l == 0) then
                        a = 0.0_dp; b = 0.5_dp * df / rm; c = 0.0_dp; dd = f - b * rm * rm
                    else if (l == 1) then
                        a = (rm * df - f) / (2.0_dp * rm**3); b = 0.0_dp; c = df - 3.0_dp * a * rm * rm; dd = 0.0_dp
                    else
                        b = (3.0_dp * f - rm * df) / (rm * rm); a = (f - b * rm * rm) / rm**3; c = 0.0_dp; dd = 0.0_dp
                    end if
                    f = a * r**3 + b * r * r + c * r + dd
                    df = 3.0_dp * a * r * r + 2.0_dp * b * r + c
                end if
                rf(mu, l) = f
                drf(mu, l) = df
            end do
        end do
        x = u(1); y = u(2); z = u(3)
        o = 0
        do l = 0, sp%lmax
            gaf = 0.0_dp
            select case (l)
            case (0)
                af(1) = 0.282094791773878_dp
            case (1)
                af(1:3) = 0.48860251190292_dp * u
                gaf(1, 1) = 0.48860251190292_dp; gaf(2, 2) = 0.48860251190292_dp; gaf(3, 3) = 0.48860251190292_dp
            case (2)
                af(1) = 0.94617469575756_dp * z * z - 0.31539156525252_dp
                gaf(3, 1) = 2.0_dp * 0.94617469575756_dp * z
                af(2) = 0.54627421529604_dp * (x * x - y * y)
                gaf(:, 2) = 0.54627421529604_dp * [2.0_dp * x, -2.0_dp * y, 0.0_dp]
                af(3) = 1.09254843059208_dp * x * y
                gaf(:, 3) = 1.09254843059208_dp * [y, x, 0.0_dp]
                af(4) = 1.09254843059208_dp * x * z
                gaf(:, 4) = 1.09254843059208_dp * [z, 0.0_dp, x]
                af(5) = 1.09254843059208_dp * y * z
                gaf(:, 5) = 1.09254843059208_dp * [0.0_dp, z, y]
            case (3)
                af(1) = 0.373176332590116_dp * (5.0_dp * z**3 - 3.0_dp * z)
                gaf(3, 1) = 0.373176332590116_dp * (15.0_dp * z * z - 3.0_dp)
                af(2) = 0.457045799464466_dp * x * (5.0_dp * z * z - 1.0_dp)
                gaf(:, 2) = 0.457045799464466_dp * [5.0_dp * z * z - 1.0_dp, 0.0_dp, 10.0_dp * x * z]
                af(3) = 0.457045799464466_dp * y * (5.0_dp * z * z - 1.0_dp)
                gaf(:, 3) = 0.457045799464466_dp * [0.0_dp, 5.0_dp * z * z - 1.0_dp, 10.0_dp * y * z]
                af(4) = 1.44530572132028_dp * z * (x * x - y * y)
                gaf(:, 4) = 1.44530572132028_dp * [2.0_dp * x * z, -2.0_dp * y * z, x * x - y * y]
                af(5) = 2.89061144264055_dp * x * y * z
                gaf(:, 5) = 2.89061144264055_dp * [y * z, x * z, x * y]
                af(6) = 0.590043589926644_dp * (x**3 - 3.0_dp * x * y * y)
                gaf(:, 6) = 0.590043589926644_dp * [3.0_dp * x * x - 3.0_dp * y * y, -6.0_dp * x * y, 0.0_dp]
                af(7) = 0.590043589926644_dp * (3.0_dp * x * x * y - y**3)
                gaf(:, 7) = 0.590043589926644_dp * [6.0_dp * x * y, 3.0_dp * x * x - 3.0_dp * y * y, 0.0_dp]
            end select
            do mu = 1, sp%nmul(l)
                if (r < 1.0e-10_dp) then
                    ror = drf(mu, l)
                else
                    ror = rf(mu, l) / r
                end if
                do k = 1, 2 * l + 1
                    o = o + 1
                    phi(o) = rf(mu, l) * af(k)
                    ud = dot_product(u, gaf(:, k))
                    tg = gaf(:, k) - u * ud
                    gphi(:, o) = drf(mu, l) * af(k) * u + ror * tg
                end do
            end do
        end do
    end subroutine pao_eval_grad

    ! ------------------------------------------------------------------
    ! Integrais de dois centros (A6)
    ! ------------------------------------------------------------------

    !> S_ab = <a|b(d)> e X_ab(:, c) = <a|(r - R_a)_c|b(d)> (Bohr), a e b da espécie sp, por quadratura
    !> de trapézio na grade cúbica de passo h centrada no átomo a (pontos na lente |r| <= rc, |r-d| <= rc).
    subroutine pao_two_center(sp, d, h, s, x)
        type(pao_species_t), intent(in) :: sp
        real(dp), intent(in) :: d(3), h
        real(dp), intent(out) :: s(:, :), x(:, :, :)
        real(dp), allocatable :: pts(:, :), pa(:, :), pb(:, :), tmp(:, :)
        integer :: nn, kx, ky, kz, np, ip, c, no
        real(dp) :: r(3), rc
        no = sp%norb
        rc = min(sp%rcut, sp%rv(sp%nmesh))
        nn = floor(rc / h)
        allocate(pts(3, (2 * nn + 1)**2), pa(no, (2 * nn + 1)**2), pb(no, (2 * nn + 1)**2), &
                 tmp(no, (2 * nn + 1)**2))
        s = 0.0_dp
        x = 0.0_dp
        do kz = -nn, nn
            r(3) = real(kz, dp) * h
            if (abs(r(3)) > rc .or. abs(r(3) - d(3)) > rc) cycle
            np = 0
            do ky = -nn, nn
                r(2) = real(ky, dp) * h
                do kx = -nn, nn
                    r(1) = real(kx, dp) * h
                    if (norm2(r) > rc .or. norm2(r - d) > rc) cycle
                    np = np + 1
                    pts(:, np) = r
                end do
            end do
            if (np == 0) cycle
            do ip = 1, np
                call pao_eval(sp, pts(:, ip), pa(:, ip))
                call pao_eval(sp, pts(:, ip) - d, pb(:, ip))
            end do
            s = s + matmul(pa(:, 1:np), transpose(pb(:, 1:np)))
            do c = 1, 3
                do ip = 1, np
                    tmp(:, ip) = pa(:, ip) * pts(c, ip)
                end do
                x(:, :, c) = x(:, :, c) + matmul(tmp(:, 1:np), transpose(pb(:, 1:np)))
            end do
        end do
        s = s * h**3
        x = x * h**3
    end subroutine pao_two_center

    !> A19: perfil planar do produto de dois orbitais ao longo do eixo cartesiano ax (1..3):
    !>   p(:, :, k + nn + 1) = h^3 sum_{r na fatia r_ax = k h} phi_a(r) phi_b(r - d),  k = -nn..nn, nn = floor(rc/h),
    !> mesma grade cúbica (centrada no átomo a) e mesma lente de pao_two_center; sum_k p(:, :, k) = S_ab da grade.
    !> d = 0: cada fatia simetrizada (exata por construção).
    subroutine pao_planar_profile(sp, d, h, ax, p)
        type(pao_species_t), intent(in) :: sp
        real(dp), intent(in) :: d(3), h
        integer, intent(in) :: ax
        real(dp), intent(out) :: p(:, :, :)
        real(dp), allocatable :: pa(:, :), pb(:, :)
        integer :: nn, k1, k2, kz, np, ip, no, b1, b2
        real(dp) :: r(3), rc
        no = sp%norb
        rc = min(sp%rcut, sp%rv(sp%nmesh))
        nn = floor(rc / h)
        b1 = modulo(ax, 3) + 1
        b2 = modulo(ax + 1, 3) + 1
        allocate(pa(no, (2 * nn + 1)**2), pb(no, (2 * nn + 1)**2))
        p = 0.0_dp
        do kz = -nn, nn
            r(ax) = real(kz, dp) * h
            if (abs(r(ax)) > rc .or. abs(r(ax) - d(ax)) > rc) cycle
            np = 0
            do k2 = -nn, nn
                r(b2) = real(k2, dp) * h
                do k1 = -nn, nn
                    r(b1) = real(k1, dp) * h
                    if (norm2(r) > rc .or. norm2(r - d) > rc) cycle
                    np = np + 1
                    call pao_eval(sp, r, pa(:, np))
                    call pao_eval(sp, r - d, pb(:, np))
                end do
            end do
            if (np == 0) cycle
            p(:, :, kz + nn + 1) = matmul(pa(:, 1:np), transpose(pb(:, 1:np))) * h**3
            if (maxval(abs(d)) <= DTOL) p(:, :, kz + nn + 1) = 0.5_dp * (p(:, :, kz + nn + 1) + &
                                                                     transpose(p(:, :, kz + nn + 1)))
        end do
    end subroutine pao_planar_profile

    pure logical function dip_canonical(d) result(yes)
        real(dp), intent(in) :: d(3)
        integer :: c
        yes = .true.
        do c = 1, 3
            if (abs(d(c)) > DTOL) then
                yes = d(c) > 0.0_dp
                return
            end if
        end do
    end function dip_canonical

    pure integer function dip_find(c, d) result(idx)
        type(pao_dip_cache_t), intent(in) :: c
        real(dp), intent(in) :: d(3)
        integer :: k
        idx = 0
        do k = 1, c%n
            if (maxval(abs(c%d(:, k) - d)) <= DTOL) then
                idx = k
                return
            end if
        end do
    end function dip_find

    !> Garante no cache as integrais de todos os vetores dl(3, m) (Bohr); as novas são calculadas
    !> em paralelo (OpenMP sobre vetores). d = 0: S e X simetrizados (exatos por paridade).
    subroutine pao_dip_ensure(c, sp, dl)
        type(pao_dip_cache_t), intent(inout) :: c
        type(pao_species_t), intent(in) :: sp
        real(dp), intent(in) :: dl(:, :)
        real(dp), allocatable :: nd(:, :), ts(:, :, :), tx(:, :, :, :), td(:, :)
        integer :: k, m, nnew, no, j
        real(dp) :: dc(3)
        no = sp%norb
        allocate(nd(3, size(dl, 2)))
        nnew = 0
        do k = 1, size(dl, 2)
            dc = dl(:, k)
            if (.not. dip_canonical(dc)) dc = -dc
            if (dip_find(c, dc) > 0) cycle
            m = 0
            do j = 1, nnew
                if (maxval(abs(nd(:, j) - dc)) <= DTOL) then
                    m = j; exit
                end if
            end do
            if (m > 0) cycle
            nnew = nnew + 1
            nd(:, nnew) = dc
        end do
        if (nnew == 0) return
        allocate(ts(no, no, nnew), tx(no, no, 3, nnew))
        !$omp parallel do default(shared) private(k) schedule(dynamic)
        do k = 1, nnew
            call pao_two_center(sp, nd(:, k), c%h, ts(:, :, k), tx(:, :, :, k))
        end do
        !$omp end parallel do
        do k = 1, nnew
            if (maxval(abs(nd(:, k))) <= DTOL) then
                ts(:, :, k) = 0.5_dp * (ts(:, :, k) + transpose(ts(:, :, k)))
                do j = 1, 3
                    tx(:, :, j, k) = 0.5_dp * (tx(:, :, j, k) + transpose(tx(:, :, j, k)))
                end do
            end if
        end do
        m = c%n
        if (m == 0) then
            c%d = nd(:, 1:nnew)
            c%s = ts
            c%x = tx
        else
            allocate(td(3, m + nnew))
            td(:, 1:m) = c%d(:, 1:m)
            td(:, m + 1:) = nd(:, 1:nnew)
            call move_alloc(td, c%d)
            block
                real(dp), allocatable :: s2(:, :, :), x2(:, :, :, :)
                allocate(s2(no, no, m + nnew), x2(no, no, 3, m + nnew))
                s2(:, :, 1:m) = c%s(:, :, 1:m)
                s2(:, :, m + 1:) = ts
                x2(:, :, :, 1:m) = c%x(:, :, :, 1:m)
                x2(:, :, :, m + 1:) = tx
                call move_alloc(s2, c%s)
                call move_alloc(x2, c%x)
            end block
        end if
        c%n = m + nnew
    end subroutine pao_dip_ensure

    !> S e X do vetor d (Bohr) a partir do cache; ok = .false. se ausente.
    pure subroutine pao_dip_get(c, d, s, x, ok)
        type(pao_dip_cache_t), intent(in) :: c
        real(dp), intent(in) :: d(3)
        real(dp), intent(out) :: s(:, :), x(:, :, :)
        logical, intent(out) :: ok
        integer :: k, j
        ok = .false.
        if (dip_canonical(d)) then
            k = dip_find(c, d)
            if (k == 0) return
            s = c%s(:, :, k)
            x = c%x(:, :, :, k)
        else
            k = dip_find(c, -d)
            if (k == 0) return
            s = transpose(c%s(:, :, k))
            do j = 1, 3
                x(:, :, j) = transpose(c%x(:, :, j, k)) + d(j) * s
            end do
        end if
        ok = .true.
    end subroutine pao_dip_get

    !> Bloco do termo de dipolo do Hamiltoniano para mu em i, nu em j (d = R_j + L - R_i, Bohr):
    !>   dH = 1/2 [F_i . X^(i) + F_j . X^(j)],  X^(i) = <mu|r - R_i|nu>, X^(j) = X^(i) - d S,
    !> F em eV/Bohr (campo elétrico no sítio; energia do elétron +F.(r - R)), dH em eV.
    pure subroutine pao_dipole_dh(fi, fj, x, s, d, dh)
        real(dp), intent(in) :: fi(3), fj(3), x(:, :, :), s(:, :), d(3)
        real(dp), intent(out) :: dh(:, :)
        integer :: c
        dh = 0.0_dp
        do c = 1, 3
            dh = dh + 0.5_dp * ((fi(c) + fj(c)) * x(:, :, c) - fj(c) * d(c) * s)
        end do
    end subroutine pao_dipole_dh

end module pao_basis
