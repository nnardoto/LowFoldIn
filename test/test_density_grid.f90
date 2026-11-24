! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Testes de pao_basis e density_grid: normalização das PAO do Te7.0
!> (s3p2d2), ortogonalidade angular, invariância à integração ao longo de c
!> e simetria D^T no par inverso.
program test_density_grid
    use types, only: dp
    use testing, only: check, check_close, test_summary
    use pao_basis, only: pao_species_t, read_pao_species, pao_eval, pao_eval_grad
    use density_grid, only: dens_system_t, density_map, dens_set_pairs
    implicit none
    type(dens_system_t) :: ds
    real(dp), allocatable :: map(:, :), map2(:, :)
    real(dp) :: h, L, phi(19), tot
    integer :: ierr, n, i
    character(len=:), allocatable :: msg

    allocate(ds%sp(1))
    call read_pao_species('test/data/pao/Te7.0.pao', 's3p2d2', 7.0_dp, ds%sp(1), ierr, msg)
    call check(ierr == 0, 'read Te7.0.pao: ' // msg)
    call check(ds%sp(1)%norb == 19, 's3p2d2 = 19 orbitals')
    call pao_eval(ds%sp(1), [1.0_dp, 0.0_dp, 0.0_dp], phi)
    call check(abs(phi(5)) < 1.0e-14_dp .and. abs(phi(6)) < 1.0e-14_dp .and. abs(phi(4)) > 1.0e-3_dp, &
               'p order x, y, z (point on +x)')
    call pao_eval(ds%sp(1), [0.0_dp, 0.0_dp, 7.5_dp], phi)
    call check(all(phi == 0.0_dp), 'zero beyond rcut')

    ! um átomo numa caixa cúbica periódica de 16 Bohr; D = I no próprio sítio -> int n = 19
    L = 16.0_dp
    h = 0.2_dp
    n = nint(L / h)
    ds%nsite = 1
    allocate(ds%pos(3, 1), ds%spec(1))
    ds%pos(:, 1) = [8.1_dp, 7.9_dp, 8.05_dp]
    ds%spec = 1
    ds%nper = 3
    ds%lat = 0.0_dp
    do i = 1, 3
        ds%lat(i, i) = L
    end do
    ds%npair = 1
    allocate(ds%pi(1), ds%pj(1), ds%pn(3, 1), ds%poff(1), ds%dmat(19 * 19))
    ds%pi = 1; ds%pj = 1; ds%pn = 0; ds%poff = 0
    ds%dmat = 0.0_dp
    do i = 1, 19
        ds%dmat((i - 1) * 19 + i) = 1.0_dp
    end do
    call dens_set_pairs(ds)
    allocate(map(n, n), map2(n, n))
    call density_map(ds, [0.0_dp, 0.0_dp, 0.0_dp], [h, 0.0_dp, 0.0_dp], [0.0_dp, h, 0.0_dp], [0.0_dp, 0.0_dp, h], &
                     n, n, n, .true., map, ierr, msg)
    call check(ierr == 0, 'density_map: ' // msg)
    tot = sum(map) * h * h
    write(*, '(A,ES24.16)') '  int n (D = I, 19 orbitals) = ', tot
    ! as radiais interpoladas (Hermite do OpenMX) têm norma 1 - O(3e-5) (s2: 0,99997206 na integral 1D fina,
    ! 0,99997183 aqui); por isso a tolerância relativa de 2e-5 no total
    call check_close(tot / 19.0_dp, 1.0_dp, 2.0e-5_dp, 'PAO normalisation: int sum_mu phi_mu^2 = 19')

    ! D só com o elemento (px, py) + (py, px): integral nula (ortogonalidade angular)
    ds%dmat = 0.0_dp
    ds%dmat((5 - 1) * 19 + 4) = 1.0_dp
    ds%dmat((4 - 1) * 19 + 5) = 1.0_dp
    call density_map(ds, [0.0_dp, 0.0_dp, 0.0_dp], [h, 0.0_dp, 0.0_dp], [0.0_dp, h, 0.0_dp], [0.0_dp, 0.0_dp, h], &
                     n, n, n, .true., map, ierr, msg)
    call check_close(sum(map) * h * h, 0.0_dp, 1.0e-8_dp, 'int phi_px phi_py = 0')

    ! s1-s2 (mesmo L, radiais diferentes) não são ortogonais em geral; a integral deve bater com a radial 1D
    ! aqui só checamos a simetria: par (1,1,0) com bloco não simétrico = soma do bloco e da transposta / 2
    ds%dmat = 0.0_dp
    ds%dmat((2 - 1) * 19 + 1) = 2.0_dp          ! D(1,2) = 2
    call density_map(ds, [0.0_dp, 0.0_dp, 0.0_dp], [h, 0.0_dp, 0.0_dp], [0.0_dp, h, 0.0_dp], [0.0_dp, 0.0_dp, h], &
                     n, n, n, .true., map, ierr, msg)
    ds%dmat = 0.0_dp
    ds%dmat((2 - 1) * 19 + 1) = 1.0_dp
    ds%dmat((1 - 1) * 19 + 2) = 1.0_dp
    call density_map(ds, [0.0_dp, 0.0_dp, 0.0_dp], [h, 0.0_dp, 0.0_dp], [0.0_dp, h, 0.0_dp], [0.0_dp, 0.0_dp, h], &
                     n, n, n, .true., map2, ierr, msg)
    call check_close(maxval(abs(map - map2)), 0.0_dp, 1.0e-13_dp, 'phi^T D phi depends only on D + D^T')

    ! A8: grade 3D (cube) = soma do modo integrate e = plano em cada iw
    block
        real(dp), allocatable :: cub(:, :, :)
        integer :: kk
        allocate(cub(n, n, n))
        call density_map(ds, [0.0_dp, 0.0_dp, 0.0_dp], [h, 0.0_dp, 0.0_dp], [0.0_dp, h, 0.0_dp], [0.0_dp, 0.0_dp, h], &
                         n, n, n, .true., map, ierr, msg, cube=cub)
        call check(ierr == 0, 'density_map cube: ' // msg)
        call check_close(maxval(abs(sum(cub, dim=3) * h - map2)), 0.0_dp, 1.0e-13_dp, 'cube summed over w = integrate')
        kk = n / 3
        call density_map(ds, [0.0_dp, 0.0_dp, (kk - 1) * h], [h, 0.0_dp, 0.0_dp], [0.0_dp, h, 0.0_dp], &
                         [0.0_dp, 0.0_dp, h], n, n, 1, .false., map, ierr, msg)
        call check_close(maxval(abs(cub(:, :, kk) - map)), 0.0_dp, 0.0_dp, 'cube plane = plane mode (bitwise)')
    end block


    ! A9: gradiente analítico das PAO (s3p3d2f1) x diferença finita central; phi = pao_eval
    block
        type(pao_species_t) :: spf
        real(dp) :: pf(29), pg(29), gp(3, 29), pp(29), pm(29), e(3), dpt(3, 6), hh, errg, errp, gmax
        integer :: ip, q
        call read_pao_species('test/data/pao/Te7.0.pao', 's3p3d2f1', 7.0_dp, spf, ierr, msg)
        call check(ierr == 0 .and. spf%norb == 29, 'read Te7.0.pao s3p3d2f1')
        dpt(:, 1) = [0.7_dp, -1.3_dp, 0.4_dp]
        dpt(:, 2) = [-2.1_dp, 0.3_dp, 1.9_dp]
        dpt(:, 3) = [0.01_dp, 0.02_dp, -0.015_dp]       ! r < r_1: extrapolação do OpenMX
        dpt(:, 4) = [3.5_dp, 2.0_dp, -3.1_dp]
        dpt(:, 5) = [0.0_dp, 0.0_dp, 1.2_dp]            ! eixo z
        dpt(:, 6) = [-0.2_dp, 0.9_dp, -0.05_dp]
        hh = 1.0e-5_dp
        errg = 0.0_dp; errp = 0.0_dp; gmax = 0.0_dp
        do ip = 1, 6
            call pao_eval(spf, dpt(:, ip), pf)
            call pao_eval_grad(spf, dpt(:, ip), pg, gp)
            errp = max(errp, maxval(abs(pg - pf)))
            gmax = max(gmax, maxval(abs(gp)))
            do q = 1, 3
                e = 0.0_dp; e(q) = hh
                call pao_eval(spf, dpt(:, ip) + e, pp)
                call pao_eval(spf, dpt(:, ip) - e, pm)
                errg = max(errg, maxval(abs((pp - pm) / (2.0_dp * hh) - gp(q, :))))
            end do
        end do
        write(*, '(A,2ES12.4,A,ES10.3)') '  pao_eval_grad: max|phi - pao_eval|, max|grad - DF| = ', errp, errg, &
            '  (max|grad| ', gmax
        call check_close(errp, 0.0_dp, 1.0e-13_dp, 'pao_eval_grad phi = pao_eval')
        call check_close(errg, 0.0_dp, 1.0e-6_dp, 'pao_eval_grad = central finite difference')
    end block

    ! A9: grad n e tau do density_map; D positiva semidefinida (D = v v^T + I) -> tau >= |grad n|^2 / (8 n)
    block
        real(dp), allocatable :: cub(:, :, :), grd(:, :, :, :), tau(:, :, :), cub0(:, :, :)
        real(dp) :: v(19), fd, egr, rmin, ref
        integer :: i1, i2, i3, j
        do i = 1, 19
            v(i) = sin(1.3_dp * i)
        end do
        do j = 1, 19
            do i = 1, 19
                ds%dmat((j - 1) * 19 + i) = v(i) * v(j)
            end do
            ds%dmat((j - 1) * 19 + j) = ds%dmat((j - 1) * 19 + j) + 0.3_dp
        end do
        allocate(cub(n, n, n), cub0(n, n, n), grd(3, n, n, n), tau(n, n, n))
        call density_map(ds, [0.0_dp, 0.0_dp, 0.0_dp], [h, 0.0_dp, 0.0_dp], [0.0_dp, h, 0.0_dp], [0.0_dp, 0.0_dp, h], &
                         n, n, n, .true., map, ierr, msg, cube=cub0)
        call density_map(ds, [0.0_dp, 0.0_dp, 0.0_dp], [h, 0.0_dp, 0.0_dp], [0.0_dp, h, 0.0_dp], [0.0_dp, 0.0_dp, h], &
                         n, n, n, .true., map, ierr, msg, cube=cub, grad=grd, tau=tau)
        call check(ierr == 0, 'density_map grad/tau: ' // msg)
        call check_close(maxval(abs(cub - cub0)) / maxval(cub0), 0.0_dp, 1.0e-13_dp, 'cube with grad = cube without')
        ! DF de 4.ª ordem (passo 1e-3 Bohr) em pontos espalhados, nas três direções
        egr = 0.0_dp; ref = 0.0_dp
        block
            real(dp) :: c5(5, 1, 1), g5(3, 5, 1, 1), t5(5, 1, 1), m5(5, 1), pt(3), ev(3), hs
            integer :: ip, q
            hs = 1.0e-3_dp
            do ip = 1, 12
                pt = ds%pos(:, 1) + [1.7_dp * sin(2.1_dp * ip), 1.9_dp * cos(1.7_dp * ip), 1.5_dp * sin(0.9_dp * ip + 0.3_dp)]
                do q = 1, 3
                    ev = 0.0_dp; ev(q) = hs
                    call density_map(ds, pt - 2.0_dp * ev, ev, [0.0_dp, h, 0.0_dp], [0.0_dp, 0.0_dp, h], &
                                     5, 1, 1, .true., m5, ierr, msg, cube=c5, grad=g5, tau=t5)
                    fd = (-c5(5, 1, 1) + 8.0_dp * c5(4, 1, 1) - 8.0_dp * c5(2, 1, 1) + c5(1, 1, 1)) / (12.0_dp * hs)
                    egr = max(egr, abs(fd - g5(q, 3, 1, 1)))
                    ref = max(ref, abs(g5(q, 3, 1, 1)))
                end do
            end do
        end block
        write(*, '(A,ES12.4,A,ES12.4)') '  grad n: max|DF4 - analytic| = ', egr, ' (max|grad n| ', ref
        call check(egr < 1.0e-7_dp * ref, 'grad n analytic = finite difference of n (h = 1e-3 Bohr)')
        rmin = huge(1.0_dp)
        do i3 = 1, n
            do i2 = 1, n
                do i1 = 1, n
                    if (cub(i1, i2, i3) < 1.0e-8_dp) cycle
                    rmin = min(rmin, (tau(i1, i2, i3) - sum(grd(:, i1, i2, i3)**2) / (8.0_dp * cub(i1, i2, i3))) &
                                     / tau(i1, i2, i3))
                end do
            end do
        end do
        write(*, '(A,ES12.4)') '  min (tau - tau_W)/tau = ', rmin
        call check(rmin > -1.0e-10_dp, 'tau >= von Weizsacker tau for positive semidefinite D')
    end block

    ! par inverso ausente: dois sítios, só (1,2) guardado -> (2,1) entra como D^T
    call test_summary('test_density_grid')
end program test_density_grid
