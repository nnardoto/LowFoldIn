! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Projeções a partir da função de Green da região central.
!>
!> DOS de Mulliken por orbital (base não ortogonal), com os termos cruzados:
!>   d_i = -1/pi Im sum_j G_ij S_ji,
!> em que j percorre todos os orbitais acoplados a i por S: o próprio bloco,
!> os blocos vizinhos m±1 e, na borda, a PL de superfície do meio
!> semi-infinito (G_{1,A} S_{A,1}, G_{N,B} S_{B,N}). A soma sobre uma camada m
!> é rho_m = -1/pi Im sum_{n=m-1}^{m+1} Tr[G_mn S_nm]; usar só Tr[G_mm S_mm]
!> quebra o teste A/A. Regiões por átomo somam d_i sobre os orbitais do grupo.
!>
!> Löwdin (convenção alternativa, restrita à região central):
!>   d_i = -1/pi Im [S_C^{1/2} G_C S_C^{1/2}]_ii, com G_C formada inteira.
!>
!> Mapa coerente para um vetor de detecção covariante w (w_i = <phi_i|chi>):
!>   I_coh = -1/pi Im w^dagger G w (via K_eff x = w, sem G inteira),
!>   I_incoh = -1/pi Im sum_m w_m^dagger G_mm w_m (só termos m = n).
!> I_coh >= 0 para qualquer w, pois A = i(G - G^dagger) é semidefinida positiva.
module green_projection
    use types, only: dp
    use zlinalg, only: lu_factor, lu_solve_identity, zmm, herm_sqrt, herm_eigvals, diag_of_product
    use central_green, only: central_t, rgf_t, solve_projection
    use spectral_projection, only: orbital_groups_t
    implicit none
    private

    public :: orbital_dos_mulliken, orbital_dos_lowdin, region_dos, coherent_map
    public :: min_neg_im_gdiag, dense_central, central_overlap_min

    real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp

contains

    !> d(ntot): DOS de Mulliken por orbital da região central.
    subroutine orbital_dos_mulliken(sys, res, d)
        type(central_t), intent(in)  :: sys
        type(rgf_t),     intent(in)  :: res
        real(dp),        intent(out) :: d(:)
        complex(dp), allocatable :: acc(:)
        integer :: m, n, i0, i1
        n = sys%nblock
        do m = 1, n
            i0 = sys%boff(m) + 1
            i1 = sys%boff(m + 1)
            allocate(acc(sys%bsize(m)))
            acc = (0.0_dp, 0.0_dp)
            call diag_of_product(res%gd(m)%a, sys%sd(m)%a, acc)
            if (m > 1) call diag_of_product(res%gl(m - 1)%a, sys%su(m - 1)%a, acc)
            if (m < n) call diag_of_product(res%gu(m)%a, sys%sl(m)%a, acc)
            if (m == 1 .and. sys%has_bot) call diag_of_product(res%glead_bot%a, sys%sxc_bot%a, acc)
            if (m == n .and. sys%has_top) call diag_of_product(res%glead_top%a, sys%sxc_top%a, acc)
            d(i0:i1) = -aimag(acc) / PI
            deallocate(acc)
        end do
    end subroutine orbital_dos_mulliken

    !> Monta K_eff = K_C - Sigma_bot - Sigma_top e S_C densas (ntot x ntot).
    subroutine dense_central(sys, keff, sc)
        type(central_t), intent(in) :: sys
        complex(dp), allocatable, intent(out) :: keff(:, :), sc(:, :)
        integer :: m, n, i0, i1, j0, j1
        n = sys%nblock
        allocate(keff(sys%ntot, sys%ntot), sc(sys%ntot, sys%ntot))
        keff = (0.0_dp, 0.0_dp)
        sc = (0.0_dp, 0.0_dp)
        do m = 1, n
            i0 = sys%boff(m) + 1; i1 = sys%boff(m + 1)
            keff(i0:i1, i0:i1) = sys%kd(m)%a
            sc(i0:i1, i0:i1) = sys%sd(m)%a
            if (m < n) then
                j0 = sys%boff(m + 1) + 1; j1 = sys%boff(m + 2)
                keff(i0:i1, j0:j1) = sys%ku(m)%a
                keff(j0:j1, i0:i1) = sys%kl(m)%a
                sc(i0:i1, j0:j1) = sys%su(m)%a
                sc(j0:j1, i0:i1) = sys%sl(m)%a
            end if
        end do
        if (sys%has_bot) then
            i0 = 1; i1 = sys%boff(2)
            keff(i0:i1, i0:i1) = keff(i0:i1, i0:i1) - sys%sig_bot%a
        end if
        if (sys%has_top) then
            i0 = sys%boff(n) + 1; i1 = sys%boff(n + 1)
            keff(i0:i1, i0:i1) = keff(i0:i1, i0:i1) - sys%sig_top%a
        end if
    end subroutine dense_central

    !> Menor autovalor de S_C (diagnóstico da base não ortogonal).
    subroutine central_overlap_min(sys, smin, ierr, errmsg)
        type(central_t), intent(in)  :: sys
        real(dp),        intent(out) :: smin
        integer,         intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        complex(dp), allocatable :: keff(:, :), sc(:, :)
        real(dp), allocatable :: w(:)
        call dense_central(sys, keff, sc)
        allocate(w(sys%ntot))
        call herm_eigvals(sc, w, ierr, errmsg)
        smin = 0.0_dp
        if (ierr == 0) smin = w(1)
    end subroutine central_overlap_min

    !> d(ntot): DOS de Löwdin por orbital (G_C formada inteira; custo O(ntot^3)).
    subroutine orbital_dos_lowdin(sys, d, smin, ierr, errmsg)
        type(central_t), intent(in)  :: sys
        real(dp),        intent(out) :: d(:)
        real(dp),        intent(out) :: smin
        integer,         intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        complex(dp), allocatable :: keff(:, :), sc(:, :), g(:, :), shalf(:, :), t(:, :), acc(:)
        integer, allocatable :: ipiv(:)
        integer :: nt

        nt = sys%ntot
        call dense_central(sys, keff, sc)
        allocate(shalf(nt, nt), g(nt, nt), t(nt, nt), acc(nt))
        call herm_sqrt(sc, shalf, smin, ierr, errmsg)
        if (ierr /= 0) return
        call lu_factor(keff, ipiv, ierr, errmsg)
        if (ierr /= 0) return
        call lu_solve_identity(keff, ipiv, g, ierr, errmsg)
        if (ierr /= 0) return
        call zmm('N', 'N', g, shalf, t)
        acc = (0.0_dp, 0.0_dp)
        call diag_of_product(shalf, t, acc)
        d = -aimag(acc) / PI
    end subroutine orbital_dos_lowdin

    !> out(g) = soma de d sobre os orbitais do grupo g.
    pure subroutine region_dos(d, groups, out)
        real(dp),               intent(in)  :: d(:)
        type(orbital_groups_t), intent(in)  :: groups
        real(dp),               intent(out) :: out(:)
        integer :: g, k
        do g = 1, groups%ngroup
            out(g) = 0.0_dp
            do k = groups%start(g), groups%start(g + 1) - 1
                out(g) = out(g) + d(groups%orb(k))
            end do
        end do
    end subroutine region_dos

    !> Mapas coerente e incoerente para os vetores de detecção w(ntot, nq).
    subroutine coherent_map(sys, res, w, icoh, iinc, ierr, errmsg)
        type(central_t), intent(in)  :: sys
        type(rgf_t),     intent(in)  :: res
        complex(dp),     intent(in)  :: w(:, :)
        real(dp),        intent(out) :: icoh(:), iinc(:)
        integer,         intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        complex(dp), allocatable :: x(:, :), gw(:, :)
        integer :: iq, m, i0, i1, nq

        nq = size(w, 2)
        x = w
        call solve_projection(sys, res, x, ierr, errmsg)
        if (ierr /= 0) return
        do iq = 1, nq
            icoh(iq) = -aimag(dot_product(w(:, iq), x(:, iq))) / PI
        end do
        iinc = 0.0_dp
        do m = 1, sys%nblock
            i0 = sys%boff(m) + 1; i1 = sys%boff(m + 1)
            allocate(gw(i1 - i0 + 1, nq))
            call zmm('N', 'N', res%gd(m)%a, w(i0:i1, :), gw)
            do iq = 1, nq
                iinc(iq) = iinc(iq) - aimag(dot_product(w(i0:i1, iq), gw(:, iq))) / PI
            end do
            deallocate(gw)
        end do
    end subroutine coherent_map

    !> min_i -Im G_ii sobre a região central (>= 0 para uma G causal).
    pure real(dp) function min_neg_im_gdiag(res) result(v)
        type(rgf_t), intent(in) :: res
        integer :: m, i
        v = huge(1.0_dp)
        do m = 1, res%nblock
            do i = 1, size(res%gd(m)%a, 1)
                v = min(v, -aimag(res%gd(m)%a(i, i)))
            end do
        end do
    end function min_neg_im_gdiag

end module green_projection
