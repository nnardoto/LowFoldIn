! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Região central bloco-tridiagonal acoplada a meios semi-infinitos.
!>
!> A região central tem N blocos (camadas principais) ordenados ao longo de
!> +a_s. Num ponto (k_par, z) ela é descrita pelos blocos de
!>   K = zS - H',   H'_{mn} = H_{mn} + (V_m + V_n)/2 S_{mn},
!> em que V_m é o potencial (deslocamento rígido ou perfil) da camada m. Os
!> meios semi-infinitos entram por autoenergias nas pontas:
!>   Sigma_bot = K_{1,A} g_A K_{A,1},   Sigma_top = K_{N,B} g_B K_{B,N},
!> com g_A, g_B as Greens de superfície (surface_green::surface_gf).
!>
!> rgf_solve calcula, sem inversas explícitas fora de blocos (só LU), os blocos
!> G_{mm}, G_{m,m+1}, G_{m+1,m} de G_C = [K_C - Sigma_bot - Sigma_top]^{-1}:
!>   Dt_1 = K_11 - Sigma_bot,  Dt_m = K_mm - K_{m,m-1} Dt_{m-1}^{-1} K_{m-1,m}
!>   Et_N = K_NN - Sigma_top,  Et_m = K_mm - K_{m,m+1} Et_{m+1}^{-1} K_{m+1,m}
!>   G_mm = [Dt_m - K_{m,m+1} Et_{m+1}^{-1} K_{m+1,m}]^{-1}   (G_NN = [Dt_N - Sigma_top]^{-1})
!>   G_{m,m+1} = -G_mm K_{m,m+1} Et_{m+1}^{-1},  G_{m+1,m} = -Et_{m+1}^{-1} K_{m+1,m} G_mm
!> e guarda as fatorações de Dt_m para solve_projection (K_eff x = w, vários
!> lados direitos), que dá w^dagger G w sem formar G inteira.
!>
!> Os blocos fora da região central acoplados por S (necessários para a DOS de
!> Mulliken na borda) são G_{1,A} = -G_11 K_{1,A} g_A e G_{N,B} = -G_NN K_{N,B} g_B.
module central_green
    use types, only: dp
    use zlinalg, only: lu_factor, lu_solve, lu_solve_identity, zmm
    implicit none
    private

    public :: zmat_t, central_t, rgf_t
    public :: central_init, assemble_central, attach_lead, rgf_solve, solve_projection
    public :: lead_block_g

    type :: zmat_t
        complex(dp), allocatable :: a(:, :)
    end type zmat_t

    !> Blocos de uma região central num ponto (k_par, z).
    type :: central_t
        integer :: nblock = 0, ntot = 0
        integer, allocatable :: bsize(:), boff(:)       !< tamanhos e offsets (nblock+1)
        type(zmat_t), allocatable :: kd(:), ku(:), kl(:)    !< K_mm, K_{m,m+1}, K_{m+1,m}
        type(zmat_t), allocatable :: sd(:), su(:), sl(:)    !< S_mm, S_{m,m+1}, S_{m+1,m}
        ! -- meios semi-infinitos (attach_lead) --
        logical :: has_bot = .false., has_top = .false.
        type(zmat_t) :: sig_bot, sig_top        !< autoenergias
        type(zmat_t) :: kcx_bot, gx_bot, sxc_bot   !< K_{1,A}, g_A, S_{A,1}
        type(zmat_t) :: kcx_top, gx_top, sxc_top   !< K_{N,B}, g_B, S_{B,N}
    end type central_t

    !> Resultado da RGF.
    type :: rgf_t
        integer :: nblock = 0
        type(zmat_t), allocatable :: gd(:), gu(:), gl(:)   !< G_mm, G_{m,m+1}, G_{m+1,m}
        type(zmat_t), allocatable :: dlu(:)                 !< LU de Dt_m
        type(zmat_t) :: glead_bot, glead_top                !< G_{1,A}, G_{N,B}
        integer, allocatable :: ipiv(:, :)                  !< (max bsize, nblock)
    end type rgf_t

contains

    !> Reserva os blocos para `nblock` camadas de tamanhos `bsize`.
    subroutine central_init(sys, bsize)
        type(central_t), intent(out) :: sys
        integer,         intent(in)  :: bsize(:)
        integer :: m, n
        n = size(bsize)
        sys%nblock = n
        sys%bsize = bsize
        allocate(sys%boff(n + 1))
        sys%boff(1) = 0
        do m = 1, n
            sys%boff(m + 1) = sys%boff(m) + bsize(m)
        end do
        sys%ntot = sys%boff(n + 1)
        allocate(sys%kd(n), sys%sd(n), sys%ku(max(n - 1, 0)), sys%kl(max(n - 1, 0)))
        allocate(sys%su(max(n - 1, 0)), sys%sl(max(n - 1, 0)))
        do m = 1, n
            allocate(sys%kd(m)%a(bsize(m), bsize(m)), sys%sd(m)%a(bsize(m), bsize(m)))
        end do
        do m = 1, n - 1
            allocate(sys%ku(m)%a(bsize(m), bsize(m + 1)), sys%kl(m)%a(bsize(m + 1), bsize(m)))
            allocate(sys%su(m)%a(bsize(m), bsize(m + 1)), sys%sl(m)%a(bsize(m + 1), bsize(m)))
        end do
    end subroutine central_init

    !> Preenche K e S de uma camada ou de um acoplamento a partir de H, S e dos
    !> potenciais. Para m = 1..nblock: hd/sd são H_mm, S_mm; para o acoplamento
    !> m -> m+1, hu/su são H_{m,m+1}, S_{m,m+1} e os blocos inferiores são os
    !> adjuntos (H hermitiana).
    subroutine assemble_central(sys, m, z, hd, sd, vm, hu, su, vnext)
        type(central_t), intent(inout) :: sys
        integer,         intent(in)    :: m
        complex(dp),     intent(in)    :: z
        complex(dp),     intent(in)    :: hd(:, :), sd(:, :)
        real(dp),        intent(in)    :: vm
        complex(dp),     intent(in), optional :: hu(:, :), su(:, :)
        real(dp),        intent(in), optional :: vnext
        real(dp) :: vbar
        sys%sd(m)%a = sd
        sys%kd(m)%a = (z - vm) * sd - hd
        if (present(hu) .and. m < sys%nblock) then
            vbar = 0.5_dp * (vm + vnext)
            sys%su(m)%a = su
            sys%sl(m)%a = conjg(transpose(su))
            sys%ku(m)%a = (z - vbar) * su - hu
            sys%kl(m)%a = (z - vbar) * sys%sl(m)%a - conjg(transpose(hu))
        end if
    end subroutine assemble_central

    !> Acopla um meio semi-infinito X a uma ponta. `where` = 'bot' (abaixo do
    !> bloco 1) ou 'top' (acima do bloco N). hxc, sxc: H_{X0,C} e S_{X0,C} na
    !> orientação "de baixo para cima" do empilhamento:
    !>   bot: H_{A,1} (linhas = orbitais de A, colunas = bloco 1)
    !>   top: H_{N,B} (linhas = bloco N, colunas = orbitais de B)
    !> vc, vx: potenciais da camada central e do meio. gx: Green de superfície.
    subroutine attach_lead(sys, where, z, hxc, sxc, vc, vx, gx)
        type(central_t),  intent(inout) :: sys
        character(len=*), intent(in)    :: where
        complex(dp),      intent(in)    :: z
        complex(dp),      intent(in)    :: hxc(:, :), sxc(:, :), gx(:, :)
        real(dp),         intent(in)    :: vc, vx
        complex(dp), allocatable :: kup(:, :), kdn(:, :), tmp(:, :)
        real(dp) :: vbar
        vbar = 0.5_dp * (vc + vx)
        ! kup = K na orientação de baixo para cima; kdn = de cima para baixo
        kup = (z - vbar) * sxc - hxc
        kdn = (z - vbar) * conjg(transpose(sxc)) - conjg(transpose(hxc))
        if (where == 'bot') then
            ! K_{A,1} = kup, K_{1,A} = kdn; Sigma = K_{1,A} g_A K_{A,1}
            sys%has_bot = .true.
            sys%kcx_bot%a = kdn
            sys%gx_bot%a = gx
            sys%sxc_bot%a = sxc                       ! S_{A,1}
            allocate(tmp(size(gx, 1), size(kup, 2)))
            call zmm('N', 'N', gx, kup, tmp)
            if (allocated(sys%sig_bot%a)) deallocate(sys%sig_bot%a)
            allocate(sys%sig_bot%a(size(kdn, 1), size(kup, 2)))
            call zmm('N', 'N', kdn, tmp, sys%sig_bot%a)
        else
            ! K_{N,B} = kup, K_{B,N} = kdn; Sigma = K_{N,B} g_B K_{B,N}
            sys%has_top = .true.
            sys%kcx_top%a = kup
            sys%gx_top%a = gx
            sys%sxc_top%a = conjg(transpose(sxc))     ! S_{B,N}
            allocate(tmp(size(gx, 1), size(kdn, 2)))
            call zmm('N', 'N', gx, kdn, tmp)
            if (allocated(sys%sig_top%a)) deallocate(sys%sig_top%a)
            allocate(sys%sig_top%a(size(kup, 1), size(kdn, 2)))
            call zmm('N', 'N', kup, tmp, sys%sig_top%a)
        end if
    end subroutine attach_lead

    !> RGF: blocos G_mm, G_{m,m+1}, G_{m+1,m} e blocos de borda G_{1,A}, G_{N,B}.
    subroutine rgf_solve(sys, res, ierr, errmsg)
        type(central_t), intent(in)  :: sys
        type(rgf_t),     intent(out) :: res
        integer,         intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg

        type(zmat_t), allocatable :: dt(:), sr(:), y(:), zz(:)
        complex(dp), allocatable :: a(:, :), x(:, :), b(:, :)
        integer, allocatable :: ipiv(:)
        integer :: n, m, nm, nmax

        ierr = 0
        errmsg = ''
        n = sys%nblock
        res%nblock = n
        if (n < 1) then
            ierr = 1; errmsg = 'rgf_solve: empty central region'; return
        end if
        nmax = maxval(sys%bsize)
        allocate(res%gd(n), res%gu(max(n - 1, 0)), res%gl(max(n - 1, 0)), res%dlu(n))
        allocate(res%ipiv(nmax, n))
        res%ipiv = 0
        allocate(dt(n), sr(n), y(n), zz(n))

        ! --- passo para a frente: Dt_m e suas fatorações ---
        do m = 1, n
            nm = sys%bsize(m)
            dt(m)%a = sys%kd(m)%a
            if (m == 1) then
                if (sys%has_bot) dt(m)%a = dt(m)%a - sys%sig_bot%a
            else
                ! X = Dt_{m-1}^{-1} K_{m-1,m}
                x = sys%ku(m - 1)%a
                call lu_solve(res%dlu(m - 1)%a, res%ipiv(:, m - 1), x, ierr, errmsg)
                if (ierr /= 0) return
                allocate(b(nm, nm))
                call zmm('N', 'N', sys%kl(m - 1)%a, x, b)
                dt(m)%a = dt(m)%a - b
                deallocate(b)
            end if
            res%dlu(m)%a = dt(m)%a
            call lu_factor(res%dlu(m)%a, ipiv, ierr, errmsg)
            if (ierr /= 0) then
                errmsg = 'rgf_solve: singular forward block'; return
            end if
            res%ipiv(1:nm, m) = ipiv(1:nm)
        end do

        ! --- passo para trás: autoenergias à direita SR_m, Y_{m+1}, Z_m ---
        ! sr(m) = K_{m,m+1} Et_{m+1}^{-1} K_{m+1,m} (sr(N) = Sigma_top)
        ! y(m+1) = Et_{m+1}^{-1} K_{m+1,m};  zz(m) = K_{m,m+1} Et_{m+1}^{-1}
        nm = sys%bsize(n)
        allocate(sr(n)%a(nm, nm))
        if (sys%has_top) then
            sr(n)%a = sys%sig_top%a
        else
            sr(n)%a = (0.0_dp, 0.0_dp)
        end if
        do m = n - 1, 1, -1
            ! Et_{m+1} = K_{m+1,m+1} - sr(m+1)
            a = sys%kd(m + 1)%a - sr(m + 1)%a
            call lu_factor(a, ipiv, ierr, errmsg)
            if (ierr /= 0) then
                errmsg = 'rgf_solve: singular backward block'; return
            end if
            y(m + 1)%a = sys%kl(m)%a
            call lu_solve(a, ipiv, y(m + 1)%a, ierr, errmsg)
            if (ierr /= 0) return
            ! zz(m)^T = Et^{-T} K_{m,m+1}^T
            b = transpose(sys%ku(m)%a)
            call lu_solve(a, ipiv, b, ierr, errmsg, trans='T')
            if (ierr /= 0) return
            zz(m)%a = transpose(b)
            allocate(sr(m)%a(sys%bsize(m), sys%bsize(m)))
            call zmm('N', 'N', sys%ku(m)%a, y(m + 1)%a, sr(m)%a)
        end do

        ! --- blocos diagonais e vizinhos ---
        do m = 1, n
            nm = sys%bsize(m)
            a = dt(m)%a - sr(m)%a
            call lu_factor(a, ipiv, ierr, errmsg)
            if (ierr /= 0) then
                errmsg = 'rgf_solve: singular diagonal block'; return
            end if
            allocate(res%gd(m)%a(nm, nm))
            call lu_solve_identity(a, ipiv, res%gd(m)%a, ierr, errmsg)
            if (ierr /= 0) return
        end do
        do m = 1, n - 1
            allocate(res%gu(m)%a(sys%bsize(m), sys%bsize(m + 1)))
            allocate(res%gl(m)%a(sys%bsize(m + 1), sys%bsize(m)))
            call zmm('N', 'N', res%gd(m)%a, zz(m)%a, res%gu(m)%a, alpha=(-1.0_dp, 0.0_dp))
            call zmm('N', 'N', y(m + 1)%a, res%gd(m)%a, res%gl(m)%a, alpha=(-1.0_dp, 0.0_dp))
        end do

        ! --- blocos de borda com os meios semi-infinitos ---
        if (sys%has_bot) call lead_block_g(res%gd(1)%a, sys%kcx_bot%a, sys%gx_bot%a, res%glead_bot%a)
        if (sys%has_top) call lead_block_g(res%gd(n)%a, sys%kcx_top%a, sys%gx_top%a, res%glead_top%a)
    end subroutine rgf_solve

    !> G_{C,X} = -G_CC K_{C,X} g_X.
    subroutine lead_block_g(gcc, kcx, gx, gcx)
        complex(dp), intent(in) :: gcc(:, :), kcx(:, :), gx(:, :)
        complex(dp), allocatable, intent(out) :: gcx(:, :)
        complex(dp), allocatable :: tmp(:, :)
        allocate(tmp(size(kcx, 1), size(gx, 2)), gcx(size(gcc, 1), size(gx, 2)))
        call zmm('N', 'N', kcx, gx, tmp)
        call zmm('N', 'N', gcc, tmp, gcx, alpha=(-1.0_dp, 0.0_dp))
    end subroutine lead_block_g

    !> Resolve K_eff x = w (K_eff = K_C - Sigma_bot - Sigma_top) para vários lados
    !> direitos, com as fatorações de rgf_solve. w(ntot, nrhs) é sobrescrito por x.
    subroutine solve_projection(sys, res, w, ierr, errmsg)
        type(central_t), intent(in)    :: sys
        type(rgf_t),     intent(in)    :: res
        complex(dp),     intent(inout) :: w(:, :)
        integer,         intent(out)   :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        complex(dp), allocatable :: t(:, :)
        integer :: n, m, i0, i1, j0, j1, nrhs

        ierr = 0
        errmsg = ''
        n = sys%nblock
        nrhs = size(w, 2)
        if (size(w, 1) /= sys%ntot) then
            ierr = 1; errmsg = 'solve_projection: wrong right-hand side length'; return
        end if
        ! eliminação para a frente: y_m = w_m - K_{m,m-1} Dt_{m-1}^{-1} y_{m-1}
        do m = 2, n
            j0 = sys%boff(m - 1) + 1; j1 = sys%boff(m)
            i0 = sys%boff(m) + 1;     i1 = sys%boff(m + 1)
            t = w(j0:j1, :)
            call lu_solve(res%dlu(m - 1)%a, res%ipiv(:, m - 1), t, ierr, errmsg)
            if (ierr /= 0) return
            block
                complex(dp), allocatable :: wm(:, :)
                wm = w(i0:i1, :)
                call zmm('N', 'N', sys%kl(m - 1)%a, t, wm, alpha=(-1.0_dp, 0.0_dp), beta=(1.0_dp, 0.0_dp))
                w(i0:i1, :) = wm
            end block
        end do
        ! retrossubstituição: x_N = G_NN y_N; x_m = Dt_m^{-1} (y_m - K_{m,m+1} x_{m+1})
        i0 = sys%boff(n) + 1; i1 = sys%boff(n + 1)
        if (allocated(t)) deallocate(t)
        allocate(t(i1 - i0 + 1, nrhs))
        call zmm('N', 'N', res%gd(n)%a, w(i0:i1, :), t)
        w(i0:i1, :) = t
        do m = n - 1, 1, -1
            i0 = sys%boff(m) + 1;     i1 = sys%boff(m + 1)
            j0 = sys%boff(m + 1) + 1; j1 = sys%boff(m + 2)
            t = w(i0:i1, :)
            call zmm('N', 'N', sys%ku(m)%a, w(j0:j1, :), t, alpha=(-1.0_dp, 0.0_dp), beta=(1.0_dp, 0.0_dp))
            call lu_solve(res%dlu(m)%a, res%ipiv(:, m), t, ierr, errmsg)
            if (ierr /= 0) return
            w(i0:i1, :) = t
        end do
    end subroutine solve_projection

end module central_green
