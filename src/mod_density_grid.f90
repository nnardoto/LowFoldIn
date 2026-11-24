! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Densidade eletrônica em espaço real a partir de uma matriz densidade
!> por pares de átomos numa base PAO (OpenMX):
!>
!>   n(r) = sum_T sum_(i,j,n) sum_{mu,nu} D^{ij,n}_{mu nu} phi_mu(r - r_i - T) phi_nu(r - r_j - L n - T)
!>
!> T percorre as translações dos eixos periódicos (colunas 1..nper de lat),
!> (i, j, n) a lista de pares (n: coeficientes inteiros nos eixos
!> periódicos). Um par ausente cujo inverso (j, i, -n) existe entra como
!> D^T (matriz densidade real simétrica: D^{ji,-n} = (D^{ij,n})^T).
!> Unidades: posições em Bohr; n em e/Bohr^3 (D em elétrons).
!>
!> Grade: r(iu, iv, iw) = origin + iu du + iv dv + iw dw (iu = 0..nu-1, ...).
!> integrate = .true.: map(iu, iv) = sum_iw n |dw| (e/Bohr^2; com dw = c/nw
!> é a integral ao longo de c, regra dos retângulos periódica);
!> integrate = .false.: map = n no plano iw = 0 (e/Bohr^3).
!> Paralelo em colunas (iu, iv): o resultado não depende do número de threads.
!>
!> Rotinas sem stop; erros por ierr/errmsg.
module density_grid
    use types, only: dp
    use pao_basis, only: pao_species_t, pao_eval, pao_eval_grad
    !$ use omp_lib
    implicit none
    private

    public :: dens_system_t, density_map, dens_set_pairs

    type :: dens_system_t
        integer :: nsite = 0, nper = 0, npair = 0
        real(dp), allocatable :: pos(:, :)          !< (3, nsite) Bohr
        integer, allocatable :: spec(:)             !< (nsite) índice em sp
        type(pao_species_t), allocatable :: sp(:)
        real(dp) :: lat(3, 3) = 0.0_dp              !< colunas 1..nper: vetores periódicos, Bohr
        integer, allocatable :: pi(:), pj(:), pn(:, :)   !< (npair), (3, npair)
        integer(8), allocatable :: poff(:)          !< offset (0-based) do bloco (ni, nj) em dmat
        real(dp), allocatable :: dmat(:)
        ! índice por átomo (montado por dens_set_pairs)
        integer, allocatable :: first(:), nxt(:)
    end type dens_system_t

contains

    !> Monta as listas encadeadas de pares por átomo i (chamar depois de preencher pi/pj/pn/poff/dmat).
    subroutine dens_set_pairs(sys)
        type(dens_system_t), intent(inout) :: sys
        integer :: p
        if (allocated(sys%first)) deallocate(sys%first, sys%nxt)
        allocate(sys%first(sys%nsite), sys%nxt(sys%npair))
        sys%first = 0
        do p = sys%npair, 1, -1
            sys%nxt(p) = sys%first(sys%pi(p))
            sys%first(sys%pi(p)) = p
        end do
    end subroutine dens_set_pairs

    !> Índice com sinal do bloco (a, b, n): +p se o par existe, -p se só o inverso existe, 0 se nenhum.
    pure integer function find_pair(sys, a, b, n) result(ip)
        type(dens_system_t), intent(in) :: sys
        integer, intent(in) :: a, b, n(3)
        integer :: p
        p = sys%first(a)
        do while (p > 0)
            if (sys%pj(p) == b .and. all(sys%pn(:, p) == n)) then
                ip = p; return
            end if
            p = sys%nxt(p)
        end do
        p = sys%first(b)
        do while (p > 0)
            if (sys%pj(p) == a .and. all(sys%pn(:, p) == -n)) then
                ip = -p; return
            end if
            p = sys%nxt(p)
        end do
        ip = 0
    end function find_pair

    !> cube (opcional, (nu, nv, nw), Bohr^-3): n em todos os pontos r0 + iw dw da coluna (A8, grade 3D).
    !> grad (opcional, (3, nu, nv, nw), Bohr^-4) e tau (opcional, (nu, nv, nw), Ha/Bohr^3) (A9, ELF):
    !>   grad n = sum D (grad phi_mu phi_nu + phi_mu grad phi_nu),  tau = 1/2 sum D grad phi_mu . grad phi_nu,
    !> com o gradiente analítico das PAO (pao_eval_grad). Com grad ou tau presentes, cube e map saem da
    !> mesma soma (iguais ao caminho sem gradiente a menos de arredondamento); sem eles nada muda.
    subroutine density_map(sys, origin, du, dv, dw, nu, nv, nw, integrate, map, ierr, errmsg, cube, grad, tau)
        type(dens_system_t), intent(in) :: sys
        real(dp), intent(in) :: origin(3), du(3), dv(3), dw(3)
        integer, intent(in) :: nu, nv, nw
        logical, intent(in) :: integrate
        real(dp), intent(out) :: map(nu, nv)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        real(dp), optional, intent(out) :: cube(:, :, :)
        real(dp), optional, intent(out) :: grad(:, :, :, :), tau(:, :, :)
        real(dp) :: bmin(3), bmax(3), cor(3), rmax, p(3), seg, dwl
        integer, allocatable :: isite(:), iimg(:, :)
        real(dp), allocatable :: ipos(:, :)
        integer :: nimg, s, a, b, c, n1, n2, n3, nr(3), pass, col, iu, iv, nmo, i
        logical :: kin
        ierr = 0
        errmsg = ''
        map = 0.0_dp
        kin = present(grad) .or. present(tau)
        if (present(grad)) then
            if (size(grad, 1) /= 3 .or. size(grad, 2) /= nu .or. size(grad, 3) /= nv .or. size(grad, 4) /= nw) then
                ierr = 3; errmsg = 'density_map: grad must be (3, nu, nv, nw)'; return
            end if
            grad = 0.0_dp
        end if
        if (present(tau)) then
            if (size(tau, 1) /= nu .or. size(tau, 2) /= nv .or. size(tau, 3) /= nw) then
                ierr = 3; errmsg = 'density_map: tau must be (nu, nv, nw)'; return
            end if
            tau = 0.0_dp
        end if
        if (present(cube)) then
            if (size(cube, 1) /= nu .or. size(cube, 2) /= nv .or. size(cube, 3) /= nw) then
                ierr = 3; errmsg = 'density_map: cube must be (nu, nv, nw)'; return
            end if
            cube = 0.0_dp
        end if
        if (sys%nsite < 1 .or. .not. allocated(sys%first)) then
            ierr = 2; errmsg = 'density_map: empty system or pairs not indexed'; return
        end if
        rmax = 0.0_dp
        nmo = 0
        do s = 1, size(sys%sp)
            rmax = max(rmax, sys%sp(s)%rcut)
            nmo = max(nmo, sys%sp(s)%norb)
        end do
        ! caixa da região
        bmin = huge(1.0_dp); bmax = -huge(1.0_dp)
        do a = 0, 1
            do b = 0, 1
                do c = 0, 1
                    cor = origin + a * (nu - 1) * du + b * (nv - 1) * dv + c * (nw - 1) * dw
                    bmin = min(bmin, cor); bmax = max(bmax, cor)
                end do
            end do
        end do
        bmin = bmin - rmax; bmax = bmax + rmax
        ! imagens periódicas dos átomos que tocam a caixa
        nr = 0
        do a = 1, sys%nper
            nr(a) = ceiling((norm2(bmax - bmin) + 2.0_dp * maxval(abs(sys%pos))) / norm2(sys%lat(:, a))) + 1
        end do
        do pass = 1, 2
            nimg = 0
            do s = 1, sys%nsite
                do n1 = -nr(1), nr(1)
                    do n2 = -nr(2), nr(2)
                        do n3 = -nr(3), nr(3)
                            p = sys%pos(:, s) + n1 * sys%lat(:, 1) + n2 * sys%lat(:, 2) + n3 * sys%lat(:, 3)
                            if (any(p < bmin) .or. any(p > bmax)) cycle
                            nimg = nimg + 1
                            if (pass == 2) then
                                isite(nimg) = s; iimg(:, nimg) = [n1, n2, n3]; ipos(:, nimg) = p
                            end if
                        end do
                    end do
                end do
            end do
            if (pass == 1) allocate(isite(nimg), iimg(3, nimg), ipos(3, nimg))
        end do
        dwl = norm2(dw)
        seg = real(nw - 1, dp)

        !$omp parallel default(shared) private(col, iu, iv, i)
        block
            integer, allocatable :: cand(:), tab(:, :), act(:)
            real(dp), allocatable :: phi(:, :), tmp(:), gph(:, :, :), tg(:, :), g3(:, :)
            real(dp) :: r0(3), r(3), d(3), t, acc, val, gv(3), tv
            integer :: nc, na, iw, k1, k2, sa, sb, p_, ni, nj, pidx, q
            allocate(cand(nimg), act(nimg), phi(nmo, nimg), tmp(nmo))
            if (kin) allocate(gph(nmo, 3, nimg), tg(nmo, 3), g3(3, nmo))
            allocate(tab(0, 0))
            !$omp do schedule(dynamic, 16)
            do col = 1, nu * nv
                iu = mod(col - 1, nu)
                iv = (col - 1) / nu
                r0 = origin + iu * du + iv * dv
                ! candidatos: imagens a até rcut do segmento r0 + [0, nw-1] dw
                nc = 0
                do i = 1, nimg
                    d = ipos(:, i) - r0
                    if (nw > 1) then
                        t = max(0.0_dp, min(seg, dot_product(d, dw) / (dwl * dwl)))
                        d = d - t * dw
                    end if
                    if (norm2(d) <= sys%sp(sys%spec(isite(i)))%rcut) then
                        nc = nc + 1
                        cand(nc) = i
                    end if
                end do
                if (nc == 0) cycle
                if (size(tab, 1) < nc) then
                    deallocate(tab); allocate(tab(nc, nc))
                end if
                do k1 = 1, nc
                    do k2 = 1, nc
                        tab(k1, k2) = find_pair(sys, isite(cand(k1)), isite(cand(k2)), &
                                                iimg(:, cand(k2)) - iimg(:, cand(k1)))
                    end do
                end do
                acc = 0.0_dp
                do iw = 0, max(nw, 1) - 1
                    r = r0 + iw * dw
                    na = 0
                    do k1 = 1, nc
                        sa = isite(cand(k1))
                        d = r - ipos(:, cand(k1))
                        if (norm2(d) > sys%sp(sys%spec(sa))%rcut) cycle
                        na = na + 1
                        act(na) = k1
                        if (kin) then
                            call pao_eval_grad(sys%sp(sys%spec(sa)), d, phi(:, na), g3)
                            gph(:, :, na) = transpose(g3)
                        else
                            call pao_eval(sys%sp(sys%spec(sa)), d, phi(:, na))
                        end if
                    end do
                    val = 0.0_dp
                    gv = 0.0_dp
                    tv = 0.0_dp
                    do k1 = 1, na
                        sa = isite(cand(act(k1)))
                        ni = sys%sp(sys%spec(sa))%norb
                        do k2 = 1, na
                            pidx = tab(act(k1), act(k2))
                            if (pidx == 0) cycle
                            sb = isite(cand(act(k2)))
                            nj = sys%sp(sys%spec(sb))%norb
                            p_ = abs(pidx)
                            if (pidx > 0) then
                                ! D (ni, nj) coluna-maior: tmp = D phi_b
                                call dgemv_lite(ni, nj, sys%dmat(sys%poff(p_) + 1:), phi(:, k2), tmp, .false.)
                            else
                                call dgemv_lite(ni, nj, sys%dmat(sys%poff(p_) + 1:), phi(:, k2), tmp, .true.)
                            end if
                            val = val + dot_product(phi(1:ni, k1), tmp(1:ni))
                            if (kin) then
                                ! tg(:, q) = D grad_q phi_b; grad n += grad phi_a . (D phi_b) + phi_a . tg
                                do q = 1, 3
                                    call dgemv_lite(ni, nj, sys%dmat(sys%poff(p_) + 1:), gph(:, q, k2), tg(:, q), &
                                                    pidx < 0)
                                    gv(q) = gv(q) + dot_product(gph(1:ni, q, k1), tmp(1:ni)) &
                                          + dot_product(phi(1:ni, k1), tg(1:ni, q))
                                    tv = tv + 0.5_dp * dot_product(gph(1:ni, q, k1), tg(1:ni, q))
                                end do
                            end if
                        end do
                    end do
                    if (present(cube)) cube(iu + 1, iv + 1, iw + 1) = val
                    if (present(grad)) grad(:, iu + 1, iv + 1, iw + 1) = gv
                    if (present(tau)) tau(iu + 1, iv + 1, iw + 1) = tv
                    if (integrate) then
                        acc = acc + val * dwl
                    else
                        acc = val
                    end if
                end do
                map(iu + 1, iv + 1) = acc
            end do
            !$omp end do
        end block
        !$omp end parallel
    end subroutine density_map

    !> y(1:m) = A x (trans = .false., A m x n) ou A^T x (trans = .true., A n x m), A coluna-maior.
    pure subroutine dgemv_lite(m, n, a, x, y, trans)
        integer, intent(in) :: m, n
        real(dp), intent(in) :: a(*), x(*)
        real(dp), intent(out) :: y(*)
        logical, intent(in) :: trans
        integer :: i, j
        if (.not. trans) then
            y(1:m) = 0.0_dp
            do j = 1, n
                y(1:m) = y(1:m) + a((j - 1) * m + 1:j * m) * x(j)
            end do
        else
            ! A é (n x m); y_i = sum_j A(j, i) x_j
            do i = 1, m
                y(i) = dot_product(a((i - 1) * n + 1:i * n), x(1:n))
            end do
        end if
    end subroutine dgemv_lite

end module density_grid
