! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Testes do módulo kpath (contracts §3.1).
program test_kpath
    use types,   only: dp
    use kpath,   only: kpath_t, build_kpath, KLABEL_LEN, kgrid_t, build_kgrid
    use testing, only: check, check_close, test_summary
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: a = 3.0_dp
    real(dp)                  :: lat(3,3), recip(3,3), vertex(3,3), step, step0
    character(len=KLABEL_LEN) :: label(3)
    character(len=256)        :: errmsg
    type(kpath_t)             :: path
    type(kgrid_t)             :: kg
    integer                   :: ierr, k
    logical                   :: ok

    ! Vetores de rede como LINHAS: a1, a2, a3 (Å)
    lat(1, :) = [a, 0.0_dp, 0.0_dp]
    lat(2, :) = [-a / 2.0_dp, a * sqrt(3.0_dp) / 2.0_dp, 0.0_dp]
    lat(3, :) = [0.0_dp, 0.0_dp, 20.0_dp]
    ! Linhas b_i com b_i·a_j = 2π δ_ij  ->  B = 2π (A^{-1})^T
    recip = 2.0_dp * pi * transpose(inv3(lat))

    ok = .true.
    do k = 1, 3
        ok = ok .and. abs(dot_product(recip(k, :), lat(k, :)) - 2.0_dp * pi) < 1.0e-12_dp
    end do
    call check(ok, "reciprocal metric b_i.a_i = 2pi")

    ! ---- 1. K-G-M -------------------------------------------------------
    vertex(:, 1) = [1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp, 0.0_dp]   ! K
    vertex(:, 2) = [0.0_dp, 0.0_dp, 0.0_dp]                     ! G
    vertex(:, 3) = [0.5_dp, 0.0_dp, 0.0_dp]                     ! M
    label = [character(len=KLABEL_LEN) :: 'K', 'G', 'M']

    call build_kpath(vertex, [10, 8], label, recip, path, ierr, errmsg)
    call check(ierr == 0, "KGM: ierr == 0 ("//trim(errmsg)//")")
    if (ierr == 0) then
        call check(path%nk == 19, "KGM: nk == 19")
        call check(path%nvertex == 3, "KGM: nvertex == 3")
        call check(all(path%vertex_k == [1, 11, 19]), "KGM: vertex_k == [1,11,19]")
        call check(all(path%kfrac(:, 11) == 0.0_dp), "KGM: kfrac(:,11) == G exactly")
        call check(all(path%kfrac(:, 1) == vertex(:, 1)), "KGM: kfrac(:,1) == K exactly")
        call check(all(path%kfrac(:, 19) == vertex(:, 3)), "KGM: kfrac(:,19) == M exactly")
        call check(path%label(1) == 'K' .and. path%label(2) == 'G' .and. path%label(3) == 'M', &
                   "KGM: labels")
        call check_close(path%kdist(1), 0.0_dp, 0.0_dp, "KGM: kdist(1) == 0")

        ok = .true.
        do k = 2, path%nk
            ok = ok .and. any(path%kfrac(:, k) /= path%kfrac(:, k - 1))
        end do
        call check(ok, "KGM: no duplicated consecutive points")

        call check(all(path%kdist(2:) > path%kdist(:path%nk - 1)), "KGM: kdist strictly increasing")

        call check_close(path%kdist(11), 4.0_dp * pi / (3.0_dp * a), 1.0e-12_dp, "KGM: kdist(11) = |K|")
        call check_close(path%kdist(19) - path%kdist(11), 2.0_dp * pi / (a * sqrt(3.0_dp)), 1.0e-12_dp, &
                         "KGM: kdist(19)-kdist(11) = |M|")

        ! Passos uniformes dentro de cada segmento
        ok = .true.
        step0 = path%kdist(2) - path%kdist(1)
        do k = 2, 11
            step = path%kdist(k) - path%kdist(k - 1)
            ok = ok .and. abs(step - step0) < 1.0e-12_dp
        end do
        call check(ok, "KGM: uniform steps in segment K-G")
        ok = .true.
        step0 = path%kdist(12) - path%kdist(11)
        do k = 12, 19
            step = path%kdist(k) - path%kdist(k - 1)
            ok = ok .and. abs(step - step0) < 1.0e-12_dp
        end do
        call check(ok, "KGM: uniform steps in segment G-M")

        ! Pontos fracionários seguem a interpolação linear
        call check(maxval(abs(path%kfrac(:, 6) - 0.5_dp * vertex(:, 1))) < 1.0e-15_dp, &
                   "KGM: midpoint of K-G")
        call check(maxval(abs(path%kfrac(:, 15) - 0.5_dp * vertex(:, 3))) < 1.0e-15_dp, &
                   "KGM: midpoint of G-M")
    end if

    ! ---- 2. Erros -------------------------------------------------------
    call build_kpath(vertex(:, 1:1), [integer ::], label(1:1), recip, path, ierr, errmsg)
    call check(ierr /= 0, "error: nvertex = 1")
    call check(len_trim(errmsg) > 0, "error: nvertex = 1 has message")

    call build_kpath(vertex, [10, 0], label, recip, path, ierr, errmsg)
    call check(ierr /= 0, "error: npoints = 0")

    vertex(:, 2) = vertex(:, 1)
    call build_kpath(vertex, [10, 8], label, recip, path, ierr, errmsg)
    call check(ierr /= 0, "error: repeated vertex (zero-length segment)")
    vertex(:, 2) = 0.0_dp

    block
        character(len=KLABEL_LEN + 4) :: longlab(3)
        longlab = [character(len=KLABEL_LEN + 4) :: 'K', 'GAMMA_POINT_TOO_LONG', 'M']
        call build_kpath(vertex, [10, 8], longlab, recip, path, ierr, errmsg)
        call check(ierr /= 0, "error: label longer than KLABEL_LEN")
    end block

    call build_kpath(vertex, [10], label, recip, path, ierr, errmsg)
    call check(ierr /= 0, "error: wrong size(npoints)")

    ! ---- grade k 2D (mapas de energia constante) --------------------------
    ! Célula hexagonal: b1 e b2 a 60 graus; grade centrada em G cobrindo b1 x b2.
    call build_kgrid([-0.5_dp, -0.5_dp, 0.0_dp], [1.0_dp, 0.0_dp, 0.0_dp], &
                     [0.0_dp, 1.0_dp, 0.0_dp], 5, 3, recip, kg, ierr, errmsg)
    call check(ierr == 0, "kgrid: ierr == 0 ("//trim(errmsg)//")")
    if (ierr == 0) then
        call check(kg%nk == 15 .and. kg%n1 == 5 .and. kg%n2 == 3, "kgrid: nk = n1*n2 = 15")
        call check(all(kg%kfrac(:, 1) == [-0.5_dp, -0.5_dp, 0.0_dp]), "kgrid: first point = origin")
        call check(all(kg%kfrac(:, 15) == [0.5_dp, 0.5_dp, 0.0_dp]), "kgrid: last point = origin+span1+span2")
        call check(all(kg%kfrac(:, 5) == [0.5_dp, -0.5_dp, 0.0_dp]), "kgrid: i runs fastest (ik=n1 -> end of span1)")
        call check(all(abs(kg%kfrac(:, 8)) < 1.0e-15_dp), "kgrid: centre point (i=3,j=2) is G")
        call check(maxval(abs(kg%kcart(:, 15) - (0.5_dp * recip(1, :) + 0.5_dp * recip(2, :)))) &
                   < 1.0e-12_dp, "kgrid: kcart = kfrac . recip")
        call check_close(norm2(kg%e_u), 1.0_dp, 1.0e-14_dp, "kgrid: |e_u| = 1")
        call check_close(norm2(kg%e_v), 1.0_dp, 1.0e-14_dp, "kgrid: |e_v| = 1")
        call check_close(dot_product(kg%e_u, kg%e_v), 0.0_dp, 1.0e-14_dp, "kgrid: e_u . e_v = 0")
        ok = .true.
        do k = 1, kg%nk
            ok = ok .and. abs(norm2(kg%kuv(:, k)) - norm2(kg%kcart(:, k))) < 1.0e-12_dp
        end do
        call check(ok, "kgrid: in-plane |k_uv| = |k_cart| (grid plane through G)")
    end if

    call build_kgrid([0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 0.0_dp, 0.0_dp], &
                     [0.0_dp, 1.0_dp, 0.0_dp], 1, 3, recip, kg, ierr, errmsg)
    call check(ierr == 1 .and. kg%nk == 0, "kgrid error: n1 < 2")
    call build_kgrid([0.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp, 0.0_dp], &
                     [0.0_dp, 1.0_dp, 0.0_dp], 3, 3, recip, kg, ierr, errmsg)
    call check(ierr == 2, "kgrid error: zero-length span1")
    call build_kgrid([0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 0.0_dp], &
                     [-0.5_dp, -0.5_dp, 0.0_dp], 3, 3, recip, kg, ierr, errmsg)
    call check(ierr == 3, "kgrid error: collinear spans")

    call test_summary("test_kpath")

contains

    !> Inversa de uma matriz 3x3 via adjunta.
    pure function inv3(m) result(r)
        real(dp), intent(in) :: m(3,3)
        real(dp)             :: r(3,3), det
        r(1,1) = m(2,2)*m(3,3) - m(2,3)*m(3,2)
        r(1,2) = m(1,3)*m(3,2) - m(1,2)*m(3,3)
        r(1,3) = m(1,2)*m(2,3) - m(1,3)*m(2,2)
        r(2,1) = m(2,3)*m(3,1) - m(2,1)*m(3,3)
        r(2,2) = m(1,1)*m(3,3) - m(1,3)*m(3,1)
        r(2,3) = m(1,3)*m(2,1) - m(1,1)*m(2,3)
        r(3,1) = m(2,1)*m(3,2) - m(2,2)*m(3,1)
        r(3,2) = m(1,2)*m(3,1) - m(1,1)*m(3,2)
        r(3,3) = m(1,1)*m(2,2) - m(1,2)*m(2,1)
        det = m(1,1)*r(1,1) + m(1,2)*r(2,1) + m(1,3)*r(3,1)
        r = r / det
    end function inv3

end program test_kpath
