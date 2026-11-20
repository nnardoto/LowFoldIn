! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Contato de Slater–Koster entre duas camadas principais (PLs) de materiais
!> diferentes: o acoplamento A/B "posto à mão", como um termo de SOC num
!> modelo tight-binding.
!>
!> Para cada par de átomos i (PL de baixo) e j (PL de cima, transladada por
!> R_par) a menos de cutoff, e cada par de camadas radiais (a em i, b em j)
!> com um parâmetro no arquivo .sk,
!>
!>   V_ab(d) = lambda * s_c * V_c(d0) * exp(-beta (|d| - d0))       (c = sigma, pi, delta)
!>   H_AB(k_par) = sum_{i,j,R_par} E_ab(d; V) exp(2 pi i (k1 n1 + k2 n2)),   S_AB = 0
!>
!> com E_ab a tabela de sk_table e d = r_j + R_par - r_i. O bloco é real no
!> espaço real e diagonal em spin (V x I_2 nos átomos não colineares, ordem
!> atom-major [up | down]). S_AB = 0: os orbitais dos dois lados são tratados
!> como ortogonais entre si; a positividade de S não muda.
!>
!> Arquivo .sk (texto, '#' inicia comentário, uma linha por par de camadas):
!>
!>   # speciesA shellA  speciesB shellB   d0(Å)  beta(1/Å)  V_sigma  V_pi  V_delta (eV)
!>   Te p1  Sb p1   3.10  1.20   1.35  -0.42  0.00
!>
!> shell = letra (s, p, d) + índice zeta (1 = primeira função radial). Uma
!> linha vale nos dois sentidos (a em A e b em B, ou b em B e a em A). Canais
!> além de min(la, lb) são ignorados. Camadas f nunca são acopladas.
module sk_contact
    use types, only: dp
    use sk_table, only: sk_block, species_basis_t, parse_basis, basis_shell_index, l_letter
    implicit none
    private

    public :: sk_entry_t, sk_model_t, sk_geom_t, sk_pairs_t
    public :: read_sk_file, write_sk_file, sk_set_bases, parse_shell, sk_find_entry, sk_entry_value
    public :: make_sk_geom, build_sk_pairs, sk_contact_at_k

    real(dp), parameter :: TWO_PI = 6.283185307179586476925286766559_dp
    integer, parameter, public :: SK_NAME_LEN = 16

    type :: sk_entry_t
        character(len=SK_NAME_LEN) :: spec_a = '', spec_b = ''
        integer  :: la = 0, za = 1, lb = 0, zb = 1
        real(dp) :: d0 = 0.0_dp
        real(dp) :: beta = 0.0_dp
        real(dp) :: v(3) = 0.0_dp      !< V_sigma, V_pi, V_delta em d0 (eV)
    end type sk_entry_t

    !> Parâmetros do contato (arquivo .sk + fatores de varredura) e bases.
    type :: sk_model_t
        integer :: nentry = 0
        type(sk_entry_t), allocatable :: e(:)
        real(dp) :: cutoff = 6.0_dp           !< Å
        real(dp) :: lambda = 1.0_dp           !< fator global
        real(dp) :: scale(3) = 1.0_dp         !< fatores por canal (sigma, pi, delta)
        real(dp) :: beta_override = -1.0_dp   !< >= 0 substitui beta de todas as linhas
        integer :: nbasis = 0
        type(species_basis_t), allocatable :: basis(:)
    end type sk_model_t

    !> Geometria de uma PL para o contato: átomos, posições (Å, no referencial
    !> da PL), deslocamento orbital e número de orbitais espaciais de cada
    !> átomo, base da espécie e se o átomo é espinorial.
    type :: sk_geom_t
        integer :: nat = 0
        integer :: norb = 0
        character(len=SK_NAME_LEN), allocatable :: label(:)
        real(dp), allocatable :: pos(:, :)        !< (3, nat)
        integer,  allocatable :: off(:)           !< (nat) primeiro orbital do átomo na PL (0-based)
        integer,  allocatable :: nsp(:)           !< (nat) orbitais espaciais
        integer,  allocatable :: ibasis(:)        !< (nat) índice em model%basis
        logical,  allocatable :: spinor(:)        !< (nat)
        integer,  allocatable :: off_sp(:)        !< (nat) deslocamento só espacial (uso de sk_fit)
    end type sk_geom_t

    !> Pares (i, j, R_par) com bloco real não nulo, para um deslocamento fixo
    !> entre as PLs.
    type :: sk_pairs_t
        integer :: npair = 0
        integer :: nrow = 0, ncol = 0
        integer,  allocatable :: ioff(:), joff(:), ni(:), nj(:)
        logical,  allocatable :: spinor(:)
        integer,  allocatable :: n12(:, :)        !< (2, npair)
        integer,  allocatable :: boff(:)          !< (npair+1) início do bloco em blk
        real(dp), allocatable :: blk(:)           !< blocos ni x nj (ordem de colunas)
        real(dp) :: dmin = huge(1.0_dp)           !< menor distância com acoplamento
    end type sk_pairs_t

contains

    !> 'p1' -> l = 1, zeta = 1.
    subroutine parse_shell(tok, l, zeta, ierr)
        character(len=*), intent(in)  :: tok
        integer,          intent(out) :: l, zeta, ierr
        character(len=:), allocatable :: t
        integer :: ios
        ierr = 0
        t = trim(adjustl(tok))
        l = index('spd', t(1:1)) - 1
        zeta = 1
        if (len(t) < 1 .or. l < 0) then
            ierr = 1; return
        end if
        if (len(t) > 1) then
            read(t(2:), *, iostat=ios) zeta
            if (ios /= 0 .or. zeta < 1) ierr = 1
        end if
    end subroutine parse_shell

    !> Lê o arquivo .sk (substitui as linhas de model; preserva bases e fatores).
    subroutine read_sk_file(filename, model, ierr, errmsg)
        character(len=*), intent(in)    :: filename
        type(sk_model_t), intent(inout) :: model
        integer,          intent(out)   :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=1024) :: line
        character(len=64) :: sa, sha, sb, shb
        type(sk_entry_t) :: en
        type(sk_entry_t), allocatable :: tmp(:)
        integer :: u, ios, nl, k, e1
        logical :: ex

        ierr = 0
        errmsg = ''
        inquire(file=filename, exist=ex)
        if (.not. ex) then
            ierr = 1; errmsg = "SK file not found: '" // trim(filename) // "'"; return
        end if
        open(newunit=u, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = "cannot open SK file '" // trim(filename) // "'"; return
        end if
        allocate(tmp(0))
        nl = 0
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            nl = nl + 1
            k = index(line, '#')
            if (k > 0) line(k:) = ''
            if (len_trim(line) == 0) cycle
            en%v = 0.0_dp
            read(line, *, iostat=ios) sa, sha, sb, shb, en%d0, en%beta, en%v(1)
            if (ios /= 0) then
                ierr = 2; errmsg = 'SK file ' // trim(filename) // ': cannot parse line ' // itoa(nl) // &
                                   " (expected 'A p1 B p1 d0 beta Vsigma [Vpi [Vdelta]]')"
                close(u); return
            end if
            read(line, *, iostat=ios) sa, sha, sb, shb, en%d0, en%beta, en%v(1), en%v(2)
            if (ios /= 0) en%v(2) = 0.0_dp
            read(line, *, iostat=ios) sa, sha, sb, shb, en%d0, en%beta, en%v(1), en%v(2), en%v(3)
            if (ios /= 0) en%v(3) = 0.0_dp
            en%spec_a = sa
            en%spec_b = sb
            call parse_shell(sha, en%la, en%za, e1)
            if (e1 == 0) call parse_shell(shb, en%lb, en%zb, e1)
            if (e1 /= 0) then
                ierr = 2; errmsg = 'SK file ' // trim(filename) // ': invalid shell on line ' // itoa(nl) // &
                                   " (use s1, p1, d1, p2, ...; f is not supported)"
                close(u); return
            end if
            if (en%d0 <= 0.0_dp .or. en%beta < 0.0_dp) then
                ierr = 2; errmsg = 'SK file ' // trim(filename) // ': d0 must be > 0 and beta >= 0 (line ' // &
                                   itoa(nl) // ')'
                close(u); return
            end if
            tmp = [tmp, en]
        end do
        close(u)
        model%nentry = size(tmp)
        call move_alloc(tmp, model%e)
        if (model%nentry == 0) then
            ierr = 2; errmsg = 'SK file ' // trim(filename) // ' has no parameter lines'
        end if
    end subroutine read_sk_file

    !> Grava linhas no formato .sk, com cabeçalho de comentários.
    subroutine write_sk_file(filename, entries, header, ierr, errmsg)
        character(len=*), intent(in) :: filename
        type(sk_entry_t), intent(in) :: entries(:)
        character(len=*), intent(in) :: header(:)
        integer,          intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: u, ios, k
        ierr = 0
        errmsg = ''
        open(newunit=u, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = "cannot write SK file '" // trim(filename) // "'"; return
        end if
        do k = 1, size(header)
            write(u, '(A,A)') '# ', trim(header(k))
        end do
        write(u, '(A)') '# speciesA shellA  speciesB shellB      d0(A)           beta(1/A)       ' // &
                        'V_sigma(eV)               V_pi(eV)                  V_delta(eV)'
        do k = 1, size(entries)
            write(u, '(A,1X,A,I0,2X,A,1X,A,I0,2X,2F17.12,3ES26.17)') trim(entries(k)%spec_a), &
                l_letter(entries(k)%la), entries(k)%za, trim(entries(k)%spec_b), l_letter(entries(k)%lb), &
                entries(k)%zb, entries(k)%d0, entries(k)%beta, entries(k)%v
        end do
        close(u)
    end subroutine write_sk_file

    !> Bases a partir de strings 'Te s3p2d2'.
    subroutine sk_set_bases(model, texts, ierr, errmsg)
        type(sk_model_t), intent(inout) :: model
        character(len=*), intent(in)    :: texts(:)
        integer,          intent(out)   :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: k, n
        n = 0
        do k = 1, size(texts)
            if (len_trim(texts(k)) > 0) n = n + 1
        end do
        ierr = 0
        errmsg = ''
        if (allocated(model%basis)) deallocate(model%basis)
        allocate(model%basis(n))
        model%nbasis = n
        n = 0
        do k = 1, size(texts)
            if (len_trim(texts(k)) == 0) cycle
            n = n + 1
            call parse_basis(texts(k), model%basis(n), ierr, errmsg)
            if (ierr /= 0) return
        end do
    end subroutine sk_set_bases

    !> Linha que acopla (spec_i, li, zi) a (spec_j, lj, zj) em qualquer
    !> ordem; 0 se não há.
    pure integer function sk_find_entry(model, spec_i, li, zi, spec_j, lj, zj)
        type(sk_model_t), intent(in) :: model
        character(len=*), intent(in) :: spec_i, spec_j
        integer, intent(in) :: li, zi, lj, zj
        integer :: k
        sk_find_entry = 0
        do k = 1, model%nentry
            associate (e => model%e(k))
                if ((trim(e%spec_a) == trim(spec_i) .and. e%la == li .and. e%za == zi .and. &
                     trim(e%spec_b) == trim(spec_j) .and. e%lb == lj .and. e%zb == zj) .or. &
                    (trim(e%spec_a) == trim(spec_j) .and. e%la == lj .and. e%za == zj .and. &
                     trim(e%spec_b) == trim(spec_i) .and. e%lb == li .and. e%zb == zi)) then
                    sk_find_entry = k
                    return
                end if
            end associate
        end do
    end function sk_find_entry

    !> (V_sigma, V_pi, V_delta) da linha k à distância d, com lambda, fatores
    !> por canal e beta_override.
    pure function sk_entry_value(model, k, d) result(v)
        type(sk_model_t), intent(in) :: model
        integer,  intent(in) :: k
        real(dp), intent(in) :: d
        real(dp) :: v(3), beta
        beta = model%e(k)%beta
        if (model%beta_override >= 0.0_dp) beta = model%beta_override
        v = model%lambda * model%scale * model%e(k)%v * exp(-beta * (d - model%e(k)%d0))
    end function sk_entry_value

    !> Geometria da PL: tau (3, norb) e orb_atom/orb_cell de layer_blocks,
    !> rótulos dos átomos da célula. orbs_per_atom: orbitais de cada átomo da
    !> célula (com spin, se houver).
    subroutine make_sk_geom(model, norb, tau, orb_atom, orb_cell, labels, orbs_per_atom, g, ierr, errmsg)
        type(sk_model_t), intent(in) :: model
        integer,          intent(in) :: norb
        real(dp),         intent(in) :: tau(:, :)
        integer,          intent(in) :: orb_atom(:), orb_cell(:), orbs_per_atom(:)
        character(len=*), intent(in) :: labels(:)
        type(sk_geom_t),  intent(out) :: g
        integer,          intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: i, ia, at, k, nat

        ierr = 0
        errmsg = ''
        ! um átomo da PL começa onde (orb_atom, orb_cell) muda
        nat = 0
        do i = 1, norb
            if (i == 1) then
                nat = nat + 1
            else if (orb_atom(i) /= orb_atom(i - 1) .or. orb_cell(i) /= orb_cell(i - 1)) then
                nat = nat + 1
            end if
        end do
        g%nat = nat
        g%norb = norb
        allocate(g%label(nat), g%pos(3, nat), g%off(nat), g%nsp(nat), g%ibasis(nat), g%spinor(nat))
        ia = 0
        do i = 1, norb
            if (i > 1) then
                if (orb_atom(i) == orb_atom(i - 1) .and. orb_cell(i) == orb_cell(i - 1)) cycle
            end if
            ia = ia + 1
            at = orb_atom(i)
            g%label(ia) = labels(at)
            g%pos(:, ia) = tau(:, i)
            g%off(ia) = i - 1
            g%ibasis(ia) = 0
            do k = 1, model%nbasis
                if (trim(model%basis(k)%name) == trim(labels(at))) g%ibasis(ia) = k
            end do
            if (g%ibasis(ia) == 0) then
                ierr = 1
                errmsg = "SK contact: no basis given for species '" // trim(labels(at)) // &
                         "' (add e.g. sk_basis = '" // trim(labels(at)) // " s3p2d2' to &INTERFACE)"
                return
            end if
            g%nsp(ia) = model%basis(g%ibasis(ia))%nspatial
            if (orbs_per_atom(at) == g%nsp(ia)) then
                g%spinor(ia) = .false.
            else if (orbs_per_atom(at) == 2 * g%nsp(ia)) then
                g%spinor(ia) = .true.
            else
                ierr = 1
                errmsg = "SK contact: species '" // trim(labels(at)) // "' has " // itoa(orbs_per_atom(at)) // &
                         " orbitals per atom, but basis '" // trim(model%basis(g%ibasis(ia))%spec) // &
                         "' gives " // itoa(g%nsp(ia)) // ' (or ' // itoa(2 * g%nsp(ia)) // ' with spin)'
                return
            end if
        end do
    end subroutine make_sk_geom

    !> Lista de pares entre a PL de baixo (glow) e a de cima (gup), com a PL de
    !> cima deslocada de disp (Å) e transladada por n1 a1 + n2 a2.
    subroutine build_sk_pairs(model, glow, gup, disp, a1, a2, pairs, ierr, errmsg)
        type(sk_model_t), intent(in)  :: model
        type(sk_geom_t),  intent(in)  :: glow, gup
        real(dp),         intent(in)  :: disp(3), a1(3), a2(3)
        type(sk_pairs_t), intent(out) :: pairs
        integer,          intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        real(dp) :: d(3), dist, span, b1(3), b2(3), area(3), aa, vv(3)
        real(dp), allocatable :: blk(:, :), store(:)
        integer :: i, j, n1, n2, nr1, nr2, si, sj, k, li, lj, oi, oj, np, nstore
        type(species_basis_t) :: bi, bj
        integer, allocatable :: ioff(:), joff(:), ni(:), nj(:), n12(:, :), boff(:)
        logical, allocatable :: spn(:)

        ierr = 0
        errmsg = ''
        pairs%nrow = glow%norb
        pairs%ncol = gup%norb
        ! duais no plano para limitar n1, n2
        area = cross(a1, a2)
        aa = dot_product(area, area)
        if (aa <= 1.0e-20_dp) then
            ierr = 1; errmsg = 'SK contact: in-plane lattice vectors are collinear'; return
        end if
        b1 = cross(a2, area) / aa
        b2 = cross(area, a1) / aa
        span = 0.0_dp
        do i = 1, glow%nat
            do j = 1, gup%nat
                span = max(span, norm2(disp + gup%pos(:, j) - glow%pos(:, i)))
            end do
        end do
        nr1 = ceiling((model%cutoff + span) * norm2(b1)) + 1
        nr2 = ceiling((model%cutoff + span) * norm2(b2)) + 1

        allocate(ioff(0), joff(0), ni(0), nj(0), n12(2, 0), boff(1), spn(0), store(0))
        boff(1) = 1
        np = 0
        nstore = 0
        do i = 1, glow%nat
            bi = model%basis(glow%ibasis(i))
            do j = 1, gup%nat
                bj = model%basis(gup%ibasis(j))
                if (glow%spinor(i) .neqv. gup%spinor(j)) then
                    ierr = 1
                    errmsg = 'SK contact: one material is spinorial (non-collinear) and the other is not'
                    return
                end if
                do n2 = -nr2, nr2
                    do n1 = -nr1, nr1
                        d = disp + gup%pos(:, j) + real(n1, dp) * a1 + real(n2, dp) * a2 - glow%pos(:, i)
                        dist = norm2(d)
                        if (dist > model%cutoff .or. dist < 1.0e-6_dp) cycle
                        allocate(blk(glow%nsp(i), gup%nsp(j)))
                        blk = 0.0_dp
                        do si = 1, bi%nshell
                            li = bi%shell(si)%l
                            if (li > 2) cycle
                            oi = bi%shell(si)%offset
                            do sj = 1, bj%nshell
                                lj = bj%shell(sj)%l
                                if (lj > 2) cycle
                                oj = bj%shell(sj)%offset
                                k = sk_find_entry(model, glow%label(i), li, bi%shell(si)%zeta, &
                                                  gup%label(j), lj, bj%shell(sj)%zeta)
                                if (k == 0) cycle
                                vv = sk_entry_value(model, k, dist)
                                blk(oi + 1:oi + 2 * li + 1, oj + 1:oj + 2 * lj + 1) = sk_block(li, lj, d, vv)
                            end do
                        end do
                        if (any(blk /= 0.0_dp)) then
                            np = np + 1
                            ioff = [ioff, glow%off(i)]
                            joff = [joff, gup%off(j)]
                            ni = [ni, glow%nsp(i)]
                            nj = [nj, gup%nsp(j)]
                            spn = [spn, glow%spinor(i)]
                            n12 = reshape([n12, [n1, n2]], [2, np])
                            store = [store, reshape(blk, [size(blk)])]
                            nstore = nstore + size(blk)
                            boff = [boff, nstore + 1]
                            pairs%dmin = min(pairs%dmin, dist)
                        end if
                        deallocate(blk)
                    end do
                end do
            end do
        end do
        pairs%npair = np
        call move_alloc(ioff, pairs%ioff)
        call move_alloc(joff, pairs%joff)
        call move_alloc(ni, pairs%ni)
        call move_alloc(nj, pairs%nj)
        call move_alloc(spn, pairs%spinor)
        call move_alloc(n12, pairs%n12)
        call move_alloc(boff, pairs%boff)
        call move_alloc(store, pairs%blk)
    end subroutine build_sk_pairs

    !> H_AB(k_par) (nrow x ncol) com fase exp(2 pi i (k1 n1 + k2 n2)).
    pure subroutine sk_contact_at_k(pairs, kpar, hc)
        type(sk_pairs_t), intent(in)  :: pairs
        real(dp),         intent(in)  :: kpar(2)
        complex(dp),      intent(out) :: hc(:, :)
        integer :: p, a, b, i0, j0, nsi, nsj, q
        complex(dp) :: ph
        hc = (0.0_dp, 0.0_dp)
        do p = 1, pairs%npair
            ph = exp(cmplx(0.0_dp, TWO_PI * (kpar(1) * pairs%n12(1, p) + kpar(2) * pairs%n12(2, p)), dp))
            i0 = pairs%ioff(p)
            j0 = pairs%joff(p)
            nsi = pairs%ni(p)
            nsj = pairs%nj(p)
            q = pairs%boff(p)
            do b = 1, nsj
                do a = 1, nsi
                    hc(i0 + a, j0 + b) = hc(i0 + a, j0 + b) + pairs%blk(q + (b - 1) * nsi + a - 1) * ph
                end do
            end do
            if (pairs%spinor(p)) then
                do b = 1, nsj
                    do a = 1, nsi
                        hc(i0 + nsi + a, j0 + nsj + b) = hc(i0 + nsi + a, j0 + nsj + b) + &
                                                         pairs%blk(q + (b - 1) * nsi + a - 1) * ph
                    end do
                end do
            end if
        end do
    end subroutine sk_contact_at_k

    pure function cross(u, v) result(w)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: w(3)
        w = [u(2) * v(3) - u(3) * v(2), u(3) * v(1) - u(1) * v(3), u(1) * v(2) - u(2) * v(1)]
    end function cross

    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=24) :: b
        write(b, '(I0)') i
        s = trim(b)
    end function itoa

end module sk_contact
