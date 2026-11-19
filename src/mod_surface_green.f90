! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Green retartada de uma superficie semi-infinita em base ortonormal.
!>
!> `build_principal_layer` agrupa celulas consecutivas na direcao normal ate
!> que a Hamiltoniana tenha apenas acoplamento entre camadas principais
!> vizinhas. `surface_greens_function` usa a decimacao iterativa de Lopez
!> Sancho. Os lados `left` e `right` trocam H01 por H01^dagger e permitem
!> representar as duas terminacoes de um mesmo bulk.
!>
!> `surface_gf` é a versão geral (base não ortogonal) usada pelo surfgreen:
!> trabalha com a matriz resolvente K(z) = zS - H em blocos de PL,
!> D = K00, U = K01, L = K10 (L /= U^dagger para eta > 0), e devolve
!> g_s = D_s^{-1} da superfície, com contagem de iterações e resíduo
!> ||(K00 - U' g_s L') g_s - I||_max (U', L' dependem do lado).
!>
!> Convenção de lado: 'top' (ou 'right') = superfície livre voltada para +a_s,
!> bulk nas PLs <= 0; 'bottom' (ou 'left') = superfície voltada para -a_s, bulk
!> nas PLs >= 0.
module surface_green
    use types, only: dp
    use zlinalg, only: lu_factor, lu_solve, lu_solve_identity, zmm, identity_matrix
    implicit none
    private

    public :: build_principal_layer, surface_greens_function
    public :: spectral_matrix, trace_surface_dos
    public :: surface_gf, surface_gf_residual, side_code
    public :: SIDE_TOP, SIDE_BOTTOM

    integer, parameter :: SIDE_TOP = 1, SIDE_BOTTOM = -1

    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            use types, only: dp
            integer     :: n, nrhs, lda, ldb, info
            integer     :: ipiv(*)
            complex(dp) :: a(lda, *), b(ldb, *)
        end subroutine zgesv
    end interface

contains

    subroutine build_principal_layer(hr, rpt, kfrac, normal_axis, h00, h01, thickness, ierr, errmsg)
        complex(dp), intent(in)  :: hr(:, :, :)
        integer,     intent(in)  :: rpt(:, :)
        real(dp),    intent(in)  :: kfrac(3)
        integer,     intent(in)  :: normal_axis
        complex(dp), allocatable, intent(out) :: h00(:, :), h01(:, :)
        integer,     intent(out) :: thickness, ierr
        character(len=:), allocatable, intent(out) :: errmsg

        real(dp), parameter :: TWO_PI = 6.283185307179586476925286766559_dp
        integer :: n, nr, ir, a, b, rz, row0, col0, p
        real(dp) :: phase_arg
        complex(dp) :: phase

        ierr = 0
        errmsg = ''
        thickness = 0
        if (normal_axis < 1 .or. normal_axis > 3) then
            ierr = 1; errmsg = 'build_principal_layer: normal_axis must be 1, 2 or 3'; return
        end if
        n = size(hr, 1)
        nr = size(hr, 3)
        if (size(hr, 2) /= n .or. size(rpt, 1) /= 3 .or. size(rpt, 2) /= nr) then
            ierr = 1; errmsg = 'build_principal_layer: inconsistent H(R)/R dimensions'; return
        end if

        thickness = max(1, maxval(abs(rpt(normal_axis, :))))
        p = n * thickness
        allocate(h00(p, p), h01(p, p))
        h00 = cmplx(0.0_dp, 0.0_dp, dp)
        h01 = cmplx(0.0_dp, 0.0_dp, dp)

        do ir = 1, nr
            rz = rpt(normal_axis, ir)
            phase_arg = TWO_PI * (dot_product(kfrac, real(rpt(:, ir), dp)) - &
                                   kfrac(normal_axis) * real(rz, dp))
            phase = exp(cmplx(0.0_dp, phase_arg, dp))
            do a = 0, thickness - 1
                b = a + rz
                if (b >= 0 .and. b < thickness) then
                    row0 = a * n
                    col0 = b * n
                    h00(row0 + 1:row0 + n, col0 + 1:col0 + n) = &
                        h00(row0 + 1:row0 + n, col0 + 1:col0 + n) + phase * hr(:, :, ir)
                else if (b >= thickness .and. b < 2 * thickness) then
                    row0 = a * n
                    col0 = (b - thickness) * n
                    h01(row0 + 1:row0 + n, col0 + 1:col0 + n) = &
                        h01(row0 + 1:row0 + n, col0 + 1:col0 + n) + phase * hr(:, :, ir)
                end if
            end do
        end do
    end subroutine build_principal_layer

    subroutine surface_greens_function(h00, h01, energy, eta, side, tol, max_iter, gs, niter, ierr, errmsg)
        complex(dp), intent(in)  :: h00(:, :), h01(:, :)
        real(dp),    intent(in)  :: energy, eta, tol
        character(len=*), intent(in) :: side
        integer,     intent(in)  :: max_iter
        complex(dp), intent(out) :: gs(:, :)
        integer,     intent(out) :: niter, ierr
        character(len=:), allocatable, intent(out) :: errmsg

        integer :: n, info
        logical :: converged
        complex(dp) :: z
        complex(dp), allocatable :: eps(:, :), eps_s(:, :), alpha(:, :), beta(:, :)
        complex(dp), allocatable :: g(:, :), agb(:, :), bga(:, :), next_a(:, :), next_b(:, :)

        ierr = 0
        errmsg = ''
        niter = 0
        n = size(h00, 1)
        if (size(h00, 2) /= n .or. any(shape(h01) /= [n, n]) .or. any(shape(gs) /= [n, n])) then
            ierr = 1; errmsg = 'surface_greens_function: inconsistent matrix dimensions'; return
        end if
        if (eta <= 0.0_dp .or. tol <= 0.0_dp .or. max_iter < 1) then
            ierr = 1; errmsg = 'surface_greens_function: eta, tol and max_iter must be positive'; return
        end if

        allocate(eps(n,n), eps_s(n,n), alpha(n,n), beta(n,n), g(n,n))
        allocate(agb(n,n), bga(n,n), next_a(n,n), next_b(n,n))
        eps = h00
        eps_s = h00
        select case (trim(side))
        case ('left')
            alpha = h01
            beta = conjg(transpose(h01))
        case ('right')
            alpha = conjg(transpose(h01))
            beta = h01
        case default
            ierr = 1; errmsg = "surface_greens_function: side must be 'left' or 'right'"; return
        end select
        z = cmplx(energy, eta, dp)
        converged = .false.

        do niter = 1, max_iter
            call resolvent(z, eps, g, info)
            if (info /= 0) then
                ierr = 2; errmsg = 'surface_greens_function: singular decimation resolvent'; return
            end if
            agb = matmul(alpha, matmul(g, beta))
            bga = matmul(beta, matmul(g, alpha))
            next_a = matmul(alpha, matmul(g, alpha))
            next_b = matmul(beta, matmul(g, beta))
            eps_s = eps_s + agb
            eps = eps + agb + bga
            alpha = next_a
            beta = next_b
            if (max(maxval(abs(alpha)), maxval(abs(beta))) < tol) then
                converged = .true.
                exit
            end if
        end do
        if (.not. converged) then
            ierr = 3; errmsg = 'surface_greens_function: decimation did not converge'; return
        end if

        call resolvent(z, eps_s, gs, info)
        if (info /= 0) then
            ierr = 2; errmsg = 'surface_greens_function: singular surface resolvent'; return
        end if

    contains

        subroutine resolvent(zin, h, inverse, lapack_info)
            complex(dp), intent(in)  :: zin, h(:, :)
            complex(dp), intent(out) :: inverse(:, :)
            integer,     intent(out) :: lapack_info
            complex(dp), allocatable :: a(:, :)
            integer, allocatable :: ipiv(:)
            integer :: k, nn
            nn = size(h, 1)
            allocate(a(nn, nn), ipiv(nn))
            a = -h
            do k = 1, nn
                a(k, k) = a(k, k) + zin
            end do
            inverse = cmplx(0.0_dp, 0.0_dp, dp)
            do k = 1, nn
                inverse(k, k) = cmplx(1.0_dp, 0.0_dp, dp)
            end do
            call zgesv(nn, nn, a, nn, ipiv, inverse, nn, lapack_info)
        end subroutine resolvent
    end subroutine surface_greens_function

    pure subroutine spectral_matrix(gs, a)
        complex(dp), intent(in)  :: gs(:, :)
        complex(dp), intent(out) :: a(:, :)
        real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp
        a = (cmplx(0.0_dp, 1.0_dp, dp) / (2.0_dp * PI)) * &
            (gs - conjg(transpose(gs)))
    end subroutine spectral_matrix

    pure real(dp) function trace_surface_dos(gs) result(dos)
        complex(dp), intent(in) :: gs(:, :)
        real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp
        integer :: i
        dos = 0.0_dp
        do i = 1, min(size(gs, 1), size(gs, 2))
            dos = dos - aimag(gs(i, i)) / PI
        end do
    end function trace_surface_dos

    !> Converte o nome do lado em SIDE_TOP/SIDE_BOTTOM (0 se inválido).
    pure integer function side_code(side)
        character(len=*), intent(in) :: side
        select case (trim(adjustl(side)))
        case ('top', 'right')
            side_code = SIDE_TOP
        case ('bottom', 'left')
            side_code = SIDE_BOTTOM
        case default
            side_code = 0
        end select
    end function side_code

    !> Green de superfície de um meio semi-infinito por decimação de
    !> López-Sancho na forma K = zS - H (base não ortogonal).
    !>
    !> k00, k01, k10: blocos K_{m,m}, K_{m,m+1}, K_{m+1,m} no mesmo z.
    !> side = SIDE_TOP: bulk nas PLs m <= 0; a superfície acopla para baixo
    !>   via U = K_{0,-1} = k10 e L = K_{-1,0} = k01.
    !> side = SIDE_BOTTOM: bulk nas PLs m >= 0; U = k01, L = k10.
    !>
    !> Iteração (g_b = D_b^{-1} aplicado por LU):
    !>   D_s <- D_s - U g_b L ;  D_b <- D_b - U g_b L - L g_b U
    !>   U   <- -U g_b U      ;  L   <- -L g_b L
    !> até max(|U|,|L|) < tol. gs = D_s^{-1}; niter = iterações usadas;
    !> residual = ||(K00 - U0 gs L0) gs - I||_max com os blocos originais.
    !> ierr: 1 argumentos, 2 matriz singular, 3 sem convergência (gs ainda é
    !> devolvido com a última iteração, para diagnóstico).
    subroutine surface_gf(k00, k01, k10, side, tol, max_iter, gs, niter, residual, ierr, errmsg)
        complex(dp), intent(in)  :: k00(:, :), k01(:, :), k10(:, :)
        integer,     intent(in)  :: side, max_iter
        real(dp),    intent(in)  :: tol
        complex(dp), intent(out) :: gs(:, :)
        integer,     intent(out) :: niter, ierr
        real(dp),    intent(out) :: residual
        character(len=:), allocatable, intent(out) :: errmsg

        complex(dp), allocatable :: ds(:, :), db(:, :), u(:, :), l(:, :)
        complex(dp), allocatable :: lu(:, :), gbl(:, :), gbu(:, :), tmp(:, :)
        integer, allocatable :: ipiv(:)
        integer :: n, it
        logical :: converged
        complex(dp), parameter :: MONE = (-1.0_dp, 0.0_dp), ONE = (1.0_dp, 0.0_dp)

        ierr = 0
        errmsg = ''
        niter = 0
        residual = huge(1.0_dp)
        n = size(k00, 1)
        if (size(k00, 2) /= n .or. any(shape(k01) /= [n, n]) .or. any(shape(k10) /= [n, n]) &
            .or. any(shape(gs) /= [n, n])) then
            ierr = 1; errmsg = 'surface_gf: inconsistent block dimensions'; return
        end if
        if (tol <= 0.0_dp .or. max_iter < 1) then
            ierr = 1; errmsg = 'surface_gf: tol and max_iter must be positive'; return
        end if
        allocate(ds(n, n), db(n, n), u(n, n), l(n, n), lu(n, n), gbl(n, n), gbu(n, n), tmp(n, n))
        select case (side)
        case (SIDE_TOP)
            u = k10
            l = k01
        case (SIDE_BOTTOM)
            u = k01
            l = k10
        case default
            ierr = 1; errmsg = "surface_gf: side must be 'top' or 'bottom'"; return
        end select
        ds = k00
        db = k00
        converged = .false.

        do it = 1, max_iter
            niter = it
            lu = db
            call lu_factor(lu, ipiv, ierr, errmsg)
            if (ierr /= 0) then
                ierr = 2; errmsg = 'surface_gf: singular bulk block during decimation'; return
            end if
            gbl = l
            call lu_solve(lu, ipiv, gbl, ierr, errmsg)       ! g_b L
            if (ierr /= 0) return
            gbu = u
            call lu_solve(lu, ipiv, gbu, ierr, errmsg)       ! g_b U
            if (ierr /= 0) return
            call zmm('N', 'N', u, gbl, tmp)                  ! U g_b L
            ds = ds - tmp
            db = db - tmp
            call zmm('N', 'N', l, gbu, db, alpha=MONE, beta=ONE)   ! - L g_b U
            call zmm('N', 'N', u, gbu, tmp)                  ! U g_b U
            u = -tmp
            call zmm('N', 'N', l, gbl, tmp)                  ! L g_b L
            l = -tmp
            if (max(maxval(abs(u)), maxval(abs(l))) < tol) then
                converged = .true.
                exit
            end if
        end do

        lu = ds
        call lu_factor(lu, ipiv, ierr, errmsg)
        if (ierr /= 0) then
            ierr = 2; errmsg = 'surface_gf: singular surface block'; return
        end if
        call lu_solve_identity(lu, ipiv, gs, ierr, errmsg)
        if (ierr /= 0) return

        residual = surface_gf_residual(k00, k01, k10, side, gs)
        if (.not. converged) then
            ierr = 3; errmsg = 'surface_gf: decimation did not converge within max_iter'
        end if
    end subroutine surface_gf

    !> ||(K00 - U gs L) gs - I||_max, com U, L do lado indicado: mede quanto gs
    !> satisfaz a equação de Dyson do meio semi-infinito.
    function surface_gf_residual(k00, k01, k10, side, gs) result(res)
        complex(dp), intent(in) :: k00(:, :), k01(:, :), k10(:, :), gs(:, :)
        integer,     intent(in) :: side
        real(dp) :: res
        complex(dp), allocatable :: a(:, :), b(:, :), c(:, :)
        integer :: n, i
        n = size(k00, 1)
        allocate(a(n, n), b(n, n), c(n, n))
        if (side == SIDE_TOP) then
            call zmm('N', 'N', gs, k01, a)       ! gs L
            call zmm('N', 'N', k10, a, b)        ! U gs L
        else
            call zmm('N', 'N', gs, k10, a)
            call zmm('N', 'N', k01, a, b)
        end if
        b = k00 - b
        call zmm('N', 'N', b, gs, c)
        do i = 1, n
            c(i, i) = c(i, i) - (1.0_dp, 0.0_dp)
        end do
        res = maxval(abs(c))
    end function surface_gf_residual

end module surface_green
