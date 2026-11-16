! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Testes do módulo spectral_projection (contracts §4), somente com matrizes sintéticas.
program test_projection
    use types,               only: dp
    use spectral_projection, only: orbital_groups_t, atom_orbital_offsets, build_orbital_groups, &
                                   overlap_times_vectors, mulliken_weights, generalized_eigensolve, &
                                   lowdin_weights
    use testing,             only: check, check_close, test_summary
    implicit none

    integer(8) :: lcg_state = 20260315_8   ! semente fixa

    call test_offsets_and_groups()
    call test_identity_overlap()
    call test_hpd_overlap()
    call test_overlapping_groups()
    call test_negative_mulliken()
    call test_edge_localization()
    call test_degeneracy()
    call test_solver_inputs_and_errors()
    call test_lowdin_identity()
    call test_lowdin_hpd()
    call test_lowdin_negative_mulliken_case()
    call test_lowdin_not_pd()
    call test_lowdin_degeneracy()

    call test_summary("test_projection")

contains

    ! ------------------------------------------------------------------
    ! Rotinas auxiliares
    ! ------------------------------------------------------------------

    !> Desvio uniforme determinístico em [-1, 1): gerador padrão mínimo de Park-Miller
    !> (x <- 48271 x mod (2^31-1)); a aritmética int64 nunca estoura.
    real(dp) function rnd()
        integer(8), parameter :: amul = 48271_8, mmod = 2147483647_8
        lcg_state = mod(amul * lcg_state, mmod)
        rnd = 2.0_dp * real(lcg_state, dp) / real(mmod, dp) - 1.0_dp
    end function rnd

    subroutine random_hermitian(n, A)
        integer,     intent(in)  :: n
        complex(dp), intent(out) :: A(n, n)
        integer                  :: i, j
        do j = 1, n
            do i = 1, j - 1
                A(i, j) = cmplx(rnd(), rnd(), dp)
                A(j, i) = conjg(A(i, j))
            end do
            A(j, j) = cmplx(rnd(), 0.0_dp, dp)
        end do
    end subroutine random_hermitian

    !> S HPD = (X^H X)/n + I, depois normalizada para diagonal unitária (como um overlap LCAO).
    subroutine random_hpd(n, S)
        integer,     intent(in)  :: n
        complex(dp), intent(out) :: S(n, n)
        complex(dp)              :: X(n, n)
        real(dp)                 :: d(n)
        integer                  :: i, j
        do j = 1, n
            do i = 1, n
                X(i, j) = cmplx(rnd(), rnd(), dp)
            end do
        end do
        S = matmul(conjg(transpose(X)), X) / real(n, dp)
        do i = 1, n
            S(i, i) = S(i, i) + 1.0_dp
        end do
        do i = 1, n
            d(i) = sqrt(real(S(i, i), dp))
        end do
        do j = 1, n
            do i = 1, n
                S(i, j) = S(i, j) / (d(i) * d(j))
            end do
            S(j, j) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
    end subroutine random_hpd

    !> Grupos formados por átomos individuais 1..nat (uma partição completa).
    subroutine atom_partition(orbs, groups)
        integer,                intent(in)  :: orbs(:)
        type(orbital_groups_t), intent(out) :: groups
        integer                             :: a, starts(size(orbs) + 1), atoms(size(orbs))
        do a = 1, size(orbs)
            starts(a) = a
            atoms(a) = a
        end do
        starts(size(orbs) + 1) = size(orbs) + 1
        call build_orbital_groups(orbs, size(orbs), starts, atoms, groups)
    end subroutine atom_partition

    ! ------------------------------------------------------------------
    ! 8. Offsets e grupos para orbs = [2,1,3]
    ! ------------------------------------------------------------------
    subroutine test_offsets_and_groups()
        integer                :: offset(4)
        type(orbital_groups_t) :: groups

        call atom_orbital_offsets([2, 1, 3], offset)
        call check(all(offset == [0, 2, 3, 6]), "offsets for orbs=[2,1,3]")

        ! g1 = {atom 3, atom 1} (ordem mantida), g2 = {atom 2}, g3 = {} , g4 = {1,2,3}
        call build_orbital_groups([2, 1, 3], 4, [1, 3, 4, 4, 7], [3, 1, 2, 1, 2, 3], groups)
        call check(groups%ngroup == 4, "groups: ngroup")
        call check(groups%norb == 6, "groups: norb")
        call check(all(groups%start == [1, 6, 7, 7, 13]), "groups: start")
        call check(all(groups%orb(1:5) == [4, 5, 6, 1, 2]), "groups: g1 = atoms 3,1")
        call check(all(groups%orb(6:6) == [3]), "groups: g2 = atom 2")
        call check(all(groups%orb(7:12) == [1, 2, 3, 4, 5, 6]), "groups: g4 = all")
    end subroutine test_offsets_and_groups

    ! ------------------------------------------------------------------
    ! 1. S = I: W = Σ |C|^2 sobre o grupo
    ! ------------------------------------------------------------------
    subroutine test_identity_overlap()
        integer, parameter     :: n = 7
        integer, parameter     :: orbs(4) = [2, 1, 3, 1]
        complex(dp)            :: H(n, n), S(n, n), C(n, n), SC(n, n)
        real(dp)               :: E(n), W(n, 4), expect
        type(orbital_groups_t) :: groups
        integer                :: ierr, i, b, g
        logical                :: ok

        call random_hermitian(n, H)
        S = (0.0_dp, 0.0_dp)
        do i = 1, n
            S(i, i) = (1.0_dp, 0.0_dp)
        end do
        call generalized_eigensolve(n, H, S, E, C, ierr)
        call check(ierr == 0, "S=I: solver ierr")
        call overlap_times_vectors(n, n, S, C, SC)
        call atom_partition(orbs, groups)
        call mulliken_weights(n, n, C, SC, groups, W)

        ok = .true.
        do g = 1, groups%ngroup
            do b = 1, n
                expect = 0.0_dp
                do i = groups%start(g), groups%start(g + 1) - 1
                    expect = expect + abs(C(groups%orb(i), b))**2
                end do
                ok = ok .and. abs(W(b, g) - expect) < 1.0e-12_dp
            end do
        end do
        call check(ok, "S=I: Mulliken == sum |C|^2")
    end subroutine test_identity_overlap

    ! ------------------------------------------------------------------
    ! 2. S HPD: ortonormalidade, resíduo, total = 1, partição soma 1
    ! ------------------------------------------------------------------
    subroutine test_hpd_overlap()
        integer, parameter     :: n = 7
        integer, parameter     :: orbs(4) = [2, 1, 3, 1]
        complex(dp)            :: H(n, n), S(n, n), C(n, n), SC(n, n), R(n, n), M(n, n)
        real(dp)               :: E(n), Wp(n, 4), Wt(n, 1), dev
        type(orbital_groups_t) :: part, tot
        integer                :: ierr, i, b

        call random_hermitian(n, H)
        call random_hpd(n, S)
        call generalized_eigensolve(n, H, S, E, C, ierr)
        call check(ierr == 0, "HPD: solver ierr")
        call overlap_times_vectors(n, n, S, C, SC)

        call check(maxval(abs(SC - matmul(S, C))) < 1.0e-12_dp, "HPD: zhemm SC == matmul(S,C)")

        M = matmul(conjg(transpose(C)), SC)
        do i = 1, n
            M(i, i) = M(i, i) - 1.0_dp
        end do
        call check(maxval(abs(M)) < 1.0e-10_dp, "HPD: C^H S C == I")

        do b = 1, n
            R(:, b) = matmul(H, C(:, b)) - E(b) * SC(:, b)
        end do
        call check(maxval(abs(R)) < 1.0e-10_dp, "HPD: residual H C - S C E")
        call check(all(E(2:) >= E(:n - 1)), "HPD: eigenvalues ascending")

        call build_orbital_groups(orbs, 1, [1, 5], [1, 2, 3, 4], tot)
        call mulliken_weights(n, n, C, SC, tot, Wt)
        call check(maxval(abs(Wt(:, 1) - 1.0_dp)) < 1.0e-10_dp, "HPD: total group W == 1")

        call atom_partition(orbs, part)
        call mulliken_weights(n, n, C, SC, part, Wp)
        dev = maxval(abs(sum(Wp, dim=2) - 1.0_dp))
        call check(dev < 1.0e-10_dp, "HPD: partition sums to 1")
    end subroutine test_hpd_overlap

    ! ------------------------------------------------------------------
    ! 3. Grupos sobrepostos: sem erro, aditivo sobre átomos individuais
    ! ------------------------------------------------------------------
    subroutine test_overlapping_groups()
        integer, parameter     :: n = 7
        integer, parameter     :: orbs(4) = [2, 1, 3, 1]
        complex(dp)            :: H(n, n), S(n, n), C(n, n), SC(n, n)
        real(dp)               :: E(n), Wo(n, 2), Wa(n, 4)
        type(orbital_groups_t) :: ov, part
        integer                :: ierr

        call random_hermitian(n, H)
        call random_hpd(n, S)
        call generalized_eigensolve(n, H, S, E, C, ierr)
        call overlap_times_vectors(n, n, S, C, SC)

        ! g1 = {1,2,3}, g2 = {3,4}: átomo 3 em ambos
        call build_orbital_groups(orbs, 2, [1, 4, 6], [1, 2, 3, 3, 4], ov)
        call check(ov%ngroup == 2 .and. size(ov%orb) == 6 + 4, "overlap: groups built")
        call mulliken_weights(n, n, C, SC, ov, Wo)
        call atom_partition(orbs, part)
        call mulliken_weights(n, n, C, SC, part, Wa)

        call check(maxval(abs(Wo(:, 1) - (Wa(:, 1) + Wa(:, 2) + Wa(:, 3)))) < 1.0e-12_dp, &
                   "overlap: g1 == sum of atoms 1,2,3")
        call check(maxval(abs(Wo(:, 2) - (Wa(:, 3) + Wa(:, 4)))) < 1.0e-12_dp, &
                   "overlap: g2 == sum of atoms 3,4")
        call check(maxval(abs(Wo(:, 1) + Wo(:, 2) - (sum(Wa, dim=2) + Wa(:, 3)))) < 1.0e-12_dp, &
                   "overlap: g1+g2 == total + atom 3 (double counted)")
    end subroutine test_overlapping_groups

    ! ------------------------------------------------------------------
    ! 4. População de Mulliken negativa, preservada (sem recorte)
    ! ------------------------------------------------------------------
    !> Dois átomos com um orbital real cada:
    !>   S = [[1, s],[s, 1]], s = 0.8;  H = [[e1, t],[t, e2]], e1 = 0, e2 = 1, t = 0.2.
    !> Equação secular generalizada det(H - E S) = 0:
    !>   (1-s^2) E^2 - (e1 + e2 - 2 t s) E + (e1 e2 - t^2) = 0
    !>   0.36 E^2 - 0.68 E - 0.04 = 0  ->  E = (0.68 ∓ sqrt(0.68^2 + 4*0.36*0.04)) / 0.72
    !>   E_1 = -0.0570975765..., E_2 = 1.9459864654...
    !> Autovetor a partir da linha 1: (e1 - E) c1 + (t - s E) c2 = 0, normalizado com
    !>   c1^2 + c2^2 + 2 s c1 c2 = 1.
    !> População de Mulliken no átomo 1: w_1 = |c1|^2 + s Re(c1* c2); no átomo 2:
    !>   w_2 = |c2|^2 + s Re(c1* c2), with w_1 + w_2 = 1.
    !> Para a banda 1 (tipo ligante, majoritariamente átomo 1): c ∝ (-1.2108, 0.2814) (sinais opostos,
    !>   |c2| pequeno) -> w_2 = 0.0792 - 0.2725 = -0.19337524528153...
    !> Para a banda 2: c ∝ (1.1454, -1.6427) -> w_1 = 1.3119 - 1.5052 = -0.19337524528153...
    !> Ambas as bandas carregam uma população negativa; o valor é mantido como está.
    subroutine test_negative_mulliken()
        integer, parameter     :: n = 2
        real(dp), parameter    :: sov = 0.8_dp, e1 = 0.0_dp, e2 = 1.0_dp, t = 0.2_dp
        real(dp), parameter    :: wneg_ref = -0.1933752452815365_dp
        complex(dp)            :: H(n, n), S(n, n), C(n, n), SC(n, n)
        real(dp)               :: E(n), W(n, 2), w_an(n, 2), Eref(2), disc
        type(orbital_groups_t) :: part
        integer                :: ierr, b

        H = reshape([cmplx(e1, 0, dp), cmplx(t, 0, dp), cmplx(t, 0, dp), cmplx(e2, 0, dp)], [2, 2])
        S = reshape([cmplx(1, 0, dp), cmplx(sov, 0, dp), cmplx(sov, 0, dp), cmplx(1, 0, dp)], [2, 2])
        call generalized_eigensolve(n, H, S, E, C, ierr)
        call check(ierr == 0, "neg: solver ierr")

        disc = sqrt((e1 + e2 - 2*t*sov)**2 - 4*(1 - sov*sov)*(e1*e2 - t*t))
        Eref(1) = ((e1 + e2 - 2*t*sov) - disc) / (2*(1 - sov*sov))
        Eref(2) = ((e1 + e2 - 2*t*sov) + disc) / (2*(1 - sov*sov))
        call check(maxval(abs(E - Eref)) < 1.0e-12_dp, "neg: eigenvalues match analytic")

        call overlap_times_vectors(n, n, S, C, SC)
        call atom_partition([1, 1], part)
        call mulliken_weights(n, n, C, SC, part, W)

        do b = 1, n
            w_an(b, 1) = abs(C(1, b))**2 + sov * real(conjg(C(1, b)) * C(2, b), dp)
            w_an(b, 2) = abs(C(2, b))**2 + sov * real(conjg(C(1, b)) * C(2, b), dp)
        end do
        call check(maxval(abs(W - w_an)) < 1.0e-12_dp, "neg: W == analytic |c|^2 + s Re(c1* c2)")
        call check(any(W < 0.0_dp), "neg: at least one negative weight exists")
        call check(minval(W) < -0.19_dp, "neg: negative weight is not clipped")
        call check_close(W(1, 2), wneg_ref, 1.0e-10_dp, "neg: W(band1, atom2) == closed-form value")
        call check_close(W(2, 1), wneg_ref, 1.0e-10_dp, "neg: W(band2, atom1) == closed-form value")
        call check(maxval(abs(sum(W, dim=2) - 1.0_dp)) < 1.0e-12_dp, "neg: populations still sum to 1")
    end subroutine test_negative_mulliken

    ! ------------------------------------------------------------------
    ! 5. Localização de borda em uma cadeia 1D com S não ortogonal
    ! ------------------------------------------------------------------
    subroutine test_edge_localization()
        integer, parameter     :: n = 12
        complex(dp)            :: H(n, n), S(n, n), C(n, n), SC(n, n)
        real(dp)               :: E(n), W(n, 2)
        integer                :: orbs(n), i, ierr, atoms(n)
        type(orbital_groups_t) :: groups

        H = (0.0_dp, 0.0_dp)
        S = (0.0_dp, 0.0_dp)
        do i = 1, n
            S(i, i) = 1.0_dp
            if (i < n) then
                H(i, i + 1) = -1.0_dp;  H(i + 1, i) = -1.0_dp
                S(i, i + 1) =  0.1_dp;  S(i + 1, i) =  0.1_dp
            end if
        end do
        H(1, 1) = -3.0_dp
        H(n, n) = -3.0_dp

        call generalized_eigensolve(n, H, S, E, C, ierr)
        call check(ierr == 0, "chain: solver ierr")
        call overlap_times_vectors(n, n, S, C, SC)

        orbs = 1
        ! g1 = edges {1,2,11,12}, g2 = centro {3..10}
        atoms(1:4) = [1, 2, 11, 12]
        atoms(5:12) = [(i, i = 3, 10)]
        call build_orbital_groups(orbs, 2, [1, 5, 13], atoms, groups)
        call mulliken_weights(n, n, C, SC, groups, W)

        call check(E(1) < -3.0_dp .and. E(2) < -3.0_dp, "chain: two edge states below band")
        call check(W(1, 1) > 0.9_dp, "chain: band 1 edge weight > 0.9")
        call check(W(2, 1) > 0.9_dp, "chain: band 2 edge weight > 0.9")
        call check(W(6, 2) > W(6, 1) .and. W(7, 2) > W(7, 1), "chain: mid-spectrum states central")
        call check(W(6, 2) > 0.5_dp .and. W(7, 2) > 0.5_dp, "chain: mid-spectrum central weight > 0.5")
        call check(maxval(abs(W(:, 1) + W(:, 2) - 1.0_dp)) < 1.0e-10_dp, "chain: partition sums to 1")
    end subroutine test_edge_localization

    ! ------------------------------------------------------------------
    ! 6. Degenerescência exata: pesos do subespaço invariantes sob rotação unitária
    ! ------------------------------------------------------------------
    subroutine test_degeneracy()
        integer, parameter     :: m = 4, n = 2 * m
        complex(dp)            :: H0(m, m), S0(m, m), C0(m, m), SC0(m, m)
        complex(dp)            :: H(n, n), S(n, n), C(n, n), SC(n, n), U(2, 2)
        complex(dp)            :: Cr(n, 2), SCr(n, 2)
        real(dp)               :: E0(m), E(n), W(n, 2), W0(m, 1), Wr(2, 2), th, ph
        type(orbital_groups_t) :: groups, g0
        integer                :: ierr, k, i, orbs(n)
        logical                :: okdeg, okrot, oktrace

        call random_hermitian(m, H0)
        call random_hpd(m, S0)
        H = (0.0_dp, 0.0_dp);  S = (0.0_dp, 0.0_dp)
        H(1:m, 1:m) = H0;  H(m+1:n, m+1:n) = H0
        S(1:m, 1:m) = S0;  S(m+1:n, m+1:n) = S0

        call generalized_eigensolve(n, H, S, E, C, ierr)
        call check(ierr == 0, "degen: solver ierr")
        call generalized_eigensolve(m, H0, S0, E0, C0, ierr)
        call overlap_times_vectors(n, n, S, C, SC)
        call overlap_times_vectors(m, m, S0, C0, SC0)

        ! g1 = orbitais {1,3} do bloco 1 mais {2} do bloco 2; g2 = orbitais {1,3} do bloco 1
        orbs = 1
        call build_orbital_groups(orbs, 2, [1, 4, 6], [1, 3, m + 2, 1, 3], groups)
        call mulliken_weights(n, n, C, SC, groups, W)
        call build_orbital_groups(orbs(1:m), 1, [1, 3], [1, 3], g0)
        call mulliken_weights(m, m, C0, SC0, g0, W0)

        ! Unitária 2x2 arbitrária
        th = 0.7_dp;  ph = 1.3_dp
        U(1, 1) = cos(th);                         U(1, 2) = -sin(th) * exp(cmplx(0, -ph, dp))
        U(2, 1) = sin(th) * exp(cmplx(0, ph, dp));  U(2, 2) = cos(th)
        U = U * exp(cmplx(0.0_dp, 0.4_dp, dp))

        okdeg = .true.;  okrot = .true.;  oktrace = .true.
        do k = 1, m
            i = 2 * k - 1
            okdeg = okdeg .and. abs(E(i) - E(i + 1)) < 1.0e-10_dp .and. abs(E(i) - E0(k)) < 1.0e-10_dp
            Cr = matmul(C(:, i:i+1), U)
            call overlap_times_vectors(n, 2, S, Cr, SCr)
            call mulliken_weights(n, 2, Cr, SCr, groups, Wr)
            okrot = okrot .and. abs(sum(Wr(:, 1)) - sum(W(i:i+1, 1))) < 1.0e-10_dp &
                          .and. abs(sum(Wr(:, 2)) - sum(W(i:i+1, 2))) < 1.0e-10_dp
            ! Expectativa baseada no traço: o par degenerado abrange uma cópia do
            ! autovetor de bloco em cada bloco, de modo que o grupo {1,3} do bloco 1 carrega
            ! exatamente o peso de bloco único W0(k).
            oktrace = oktrace .and. abs(sum(W(i:i+1, 2)) - W0(k, 1)) < 1.0e-10_dp
        end do
        call check(okdeg, "degen: eigenvalues doubly degenerate and equal to block spectrum")
        call check(okrot, "degen: subspace weights invariant under 2x2 unitary rotation")
        call check(oktrace, "degen: subspace weight equals single-block (trace) value")
    end subroutine test_degeneracy

    ! ------------------------------------------------------------------
    ! 7. Entradas intocadas; S não PD retorna ierr /= 0 sem parar
    ! ------------------------------------------------------------------
    subroutine test_solver_inputs_and_errors()
        integer, parameter :: n = 5
        complex(dp)        :: H(n, n), S(n, n), Hc(n, n), Sc(n, n), C(n, n)
        real(dp)           :: E(n)
        integer            :: ierr

        call random_hermitian(n, H)
        call random_hpd(n, S)
        Hc = H;  Sc = S
        call generalized_eigensolve(n, H, S, E, C, ierr)
        call check(ierr == 0, "inputs: solver ierr")
        call check(all(H == Hc), "inputs: H unchanged")
        call check(all(S == Sc), "inputs: S unchanged")

        ! S não PD: identidade com uma entrada diagonal negativa
        S = (0.0_dp, 0.0_dp)
        S(1, 1) = 1.0_dp;  S(2, 2) = 1.0_dp;  S(3, 3) = -0.5_dp;  S(4, 4) = 1.0_dp;  S(5, 5) = 1.0_dp
        Sc = S
        call generalized_eigensolve(n, H, S, E, C, ierr)
        call check(ierr /= 0, "non-PD S: ierr /= 0 (no stop)")
        call check(ierr > n, "non-PD S: ierr = info > n (zpotrf failure)")
        call check(all(H == Hc) .and. all(S == Sc), "non-PD S: inputs unchanged")
    end subroutine test_solver_inputs_and_errors

    ! ------------------------------------------------------------------
    ! L1. Löwdin com S = I é igual a |C|^2 (= Mulliken)
    ! ------------------------------------------------------------------
    subroutine test_lowdin_identity()
        integer, parameter     :: n = 5
        complex(dp)            :: H(n, n), S(n, n), C(n, n), SC(n, n)
        real(dp)               :: E(n), WL(n, n), WM(n, n), smin
        type(orbital_groups_t) :: part
        integer                :: ierr, i

        call random_hermitian(n, H)
        S = (0.0_dp, 0.0_dp)
        do i = 1, n
            S(i, i) = (1.0_dp, 0.0_dp)
        end do
        call generalized_eigensolve(n, H, S, E, C, ierr)
        call check(ierr == 0, "lowdin S=I: solver ierr")
        call atom_partition([(1, i = 1, n)], part)
        call lowdin_weights(n, n, S, C, part, WL, 1.0e-10_dp, smin, ierr)
        call check(ierr == 0, "lowdin S=I: ierr = 0")
        call check_close(smin, 1.0_dp, 1.0e-12_dp, "lowdin S=I: smin = 1")
        call check(maxval(abs(WL - transpose(abs(C)**2))) < 1.0e-12_dp, "lowdin S=I: W == |C|^2")
        call overlap_times_vectors(n, n, S, C, SC)
        call mulliken_weights(n, n, C, SC, part, WM)
        call check(maxval(abs(WL - WM)) < 1.0e-12_dp, "lowdin S=I: Lowdin == Mulliken")
    end subroutine test_lowdin_identity

    ! ------------------------------------------------------------------
    ! L2. S HPD: não negativo, total = 1, partição soma 1; S, C intocadas
    ! ------------------------------------------------------------------
    subroutine test_lowdin_hpd()
        integer, parameter     :: n = 7
        integer, parameter     :: orbs(3) = [2, 1, 4]
        complex(dp)            :: H(n, n), S(n, n), C(n, n), Sc(n, n), Cc(n, n)
        real(dp)               :: E(n), W(n, 3), Wt(n, 1), smin
        type(orbital_groups_t) :: part, tot
        integer                :: ierr

        call random_hermitian(n, H)
        call random_hpd(n, S)
        call generalized_eigensolve(n, H, S, E, C, ierr)
        call check(ierr == 0, "lowdin HPD: solver ierr")
        Sc = S;  Cc = C
        call atom_partition(orbs, part)
        call build_orbital_groups(orbs, 1, [1, 4], [1, 2, 3], tot)
        call lowdin_weights(n, n, S, C, part, W, 1.0e-10_dp, smin, ierr)
        call check(ierr == 0, "lowdin HPD: ierr = 0")
        call check(smin > 0.0_dp, "lowdin HPD: smin > 0")
        call check(minval(W) >= -1.0e-12_dp, "lowdin HPD: all weights >= -1e-12")
        call check(maxval(abs(sum(W, dim=2) - 1.0_dp)) < 1.0e-10_dp, "lowdin HPD: partition sums to 1")
        call lowdin_weights(n, n, S, C, tot, Wt, 1.0e-10_dp, smin, ierr)
        call check(maxval(abs(Wt(:, 1) - 1.0_dp)) < 1.0e-10_dp, "lowdin HPD: total group W == 1")
        call check(all(S == Sc) .and. all(C == Cc), "lowdin HPD: S and C unchanged")
    end subroutine test_lowdin_hpd

    ! ------------------------------------------------------------------
    ! L3. O caso de 2 orbitais com Mulliken negativo: Löwdin é não negativo
    ! ------------------------------------------------------------------
    subroutine test_lowdin_negative_mulliken_case()
        integer, parameter     :: n = 2
        real(dp), parameter    :: sov = 0.8_dp, e1 = 0.0_dp, e2 = 1.0_dp, t = 0.2_dp
        complex(dp)            :: H(n, n), S(n, n), C(n, n)
        real(dp)               :: E(n), W(n, 2), smin
        type(orbital_groups_t) :: part
        integer                :: ierr

        H = reshape([cmplx(e1, 0, dp), cmplx(t, 0, dp), cmplx(t, 0, dp), cmplx(e2, 0, dp)], [2, 2])
        S = reshape([cmplx(1, 0, dp), cmplx(sov, 0, dp), cmplx(sov, 0, dp), cmplx(1, 0, dp)], [2, 2])
        call generalized_eigensolve(n, H, S, E, C, ierr)
        call atom_partition([1, 1], part)
        call lowdin_weights(n, n, S, C, part, W, 1.0e-10_dp, smin, ierr)
        call check(ierr == 0, "lowdin neg-case: ierr = 0")
        call check_close(smin, 1.0_dp - sov, 1.0e-12_dp, "lowdin neg-case: smin = 1 - s")
        call check(minval(W) >= 0.0_dp, "lowdin neg-case: weights non-negative")
        call check(maxval(abs(sum(W, dim=2) - 1.0_dp)) < 1.0e-12_dp, "lowdin neg-case: weights sum to 1")
    end subroutine test_lowdin_negative_mulliken_case

    ! ------------------------------------------------------------------
    ! L4. S quase singular / indefinida: ierr = 1, smin <= tol, S inalterada
    ! ------------------------------------------------------------------
    subroutine test_lowdin_not_pd()
        integer, parameter     :: n = 2
        real(dp), parameter    :: tol = 1.0e-10_dp
        complex(dp)            :: S(n, n), Sc(n, n), C(n, n)
        real(dp)               :: W(n, 2), smin
        type(orbital_groups_t) :: part
        integer                :: ierr

        call atom_partition([1, 1], part)
        C = reshape([cmplx(1, 0, dp), cmplx(0, 0, dp), cmplx(0, 0, dp), cmplx(1, 0, dp)], [2, 2])
        ! autovalores 2 + 1e-14 e 1e-14 (<= tol)
        S = reshape([cmplx(1 + 1.0e-14_dp, 0, dp), cmplx(1, 0, dp), &
                     cmplx(1, 0, dp), cmplx(1 + 1.0e-14_dp, 0, dp)], [2, 2])
        Sc = S
        call lowdin_weights(n, n, S, C, part, W, tol, smin, ierr)
        call check(ierr == 1, "lowdin near-singular S: ierr = 1")
        call check(smin <= tol, "lowdin near-singular S: smin <= tol_s_min")
        call check(all(S == Sc), "lowdin near-singular S: S unchanged")
        ! indefinida: autovalores 1 +/- 1.5
        S = reshape([cmplx(1, 0, dp), cmplx(1.5_dp, 0, dp), cmplx(1.5_dp, 0, dp), cmplx(1, 0, dp)], [2, 2])
        call lowdin_weights(n, n, S, C, part, W, tol, smin, ierr)
        call check(ierr == 1, "lowdin indefinite S: ierr = 1")
        call check_close(smin, -0.5_dp, 1.0e-12_dp, "lowdin indefinite S: smin = -0.5 (not floored)")
    end subroutine test_lowdin_not_pd

    ! ------------------------------------------------------------------
    ! L5. Par degenerado: pesos de Löwdin somados invariantes sob rotação
    ! ------------------------------------------------------------------
    subroutine test_lowdin_degeneracy()
        integer, parameter     :: m = 4, n = 2 * m
        complex(dp)            :: H0(m, m), S0(m, m), H(n, n), S(n, n), C(n, n), U(2, 2), Cr(n, 2)
        real(dp)               :: E(n), W(n, 2), Wr(2, 2), smin, th, ph
        type(orbital_groups_t) :: groups
        integer                :: ierr, k, i, orbs(n)
        logical                :: okdeg, okrot

        call random_hermitian(m, H0)
        call random_hpd(m, S0)
        H = (0.0_dp, 0.0_dp);  S = (0.0_dp, 0.0_dp)
        H(1:m, 1:m) = H0;  H(m+1:n, m+1:n) = H0
        S(1:m, 1:m) = S0;  S(m+1:n, m+1:n) = S0
        call generalized_eigensolve(n, H, S, E, C, ierr)
        call check(ierr == 0, "lowdin degen: solver ierr")

        orbs = 1
        call build_orbital_groups(orbs, 2, [1, 4, 6], [1, 3, m + 2, 1, 3], groups)
        call lowdin_weights(n, n, S, C, groups, W, 1.0e-10_dp, smin, ierr)
        call check(ierr == 0, "lowdin degen: ierr = 0")

        th = 0.7_dp;  ph = 1.3_dp
        U(1, 1) = cos(th);                         U(1, 2) = -sin(th) * exp(cmplx(0, -ph, dp))
        U(2, 1) = sin(th) * exp(cmplx(0, ph, dp));  U(2, 2) = cos(th)
        U = U * exp(cmplx(0.0_dp, 0.4_dp, dp))

        okdeg = .true.;  okrot = .true.
        do k = 1, m
            i = 2 * k - 1
            okdeg = okdeg .and. abs(E(i) - E(i + 1)) < 1.0e-10_dp
            Cr = matmul(C(:, i:i+1), U)
            call lowdin_weights(n, 2, S, Cr, groups, Wr, 1.0e-10_dp, smin, ierr)
            okrot = okrot .and. ierr == 0 &
                          .and. abs(sum(Wr(:, 1)) - sum(W(i:i+1, 1))) < 1.0e-10_dp &
                          .and. abs(sum(Wr(:, 2)) - sum(W(i:i+1, 2))) < 1.0e-10_dp
        end do
        call check(okdeg, "lowdin degen: eigenvalues doubly degenerate")
        call check(okrot, "lowdin degen: subspace weights invariant under 2x2 unitary rotation")
    end subroutine test_lowdin_degeneracy

end program test_projection
