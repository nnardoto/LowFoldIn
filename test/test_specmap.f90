! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Testes de spectral_broadening e spectral_map (contracts.md §7).
!> Um spectral_data_t sintético é construído no código, escrito com
!> write_spectral_file em build/, relido, alargado e escrito como mapas.
!> Todos os arquivos temporários são removidos ao final.
program test_specmap
    use types,               only: dp
    use testing,             only: check, check_close, test_summary
    use spectral_data,       only: spectral_data_t, write_spectral_file, read_spectral_file, &
                                   SPIN_NONE, METHOD_MULLIKEN, EREF_FERMI
    use spectral_broadening, only: KERNEL_LORENTZIAN, KERNEL_GAUSSIAN, lorentzian, gaussian, &
                                   fermi_dirac, energy_grid, broaden
    use spectral_map,        only: parse_map_expression, write_map_file, write_kticks_file, &
                                   map_file_tag
    implicit none

    real(dp), parameter         :: PI = acos(-1.0_dp)
    character(len=*), parameter :: BIN   = 'build/test_specmap.spectral.bin'
    character(len=*), parameter :: TRUNC = 'build/test_specmap.truncated.bin'
    character(len=*), parameter :: MAPF  = 'build/test_specmap.map.dat'
    character(len=*), parameter :: TICKF = 'build/test_specmap.kticks.gp'

    type(spectral_data_t) :: sd0, sd
    integer               :: ierr
    character(len=512)    :: errmsg

    call build_synthetic(sd0)
    call write_spectral_file(BIN, sd0, ierr, errmsg)
    call check(ierr == 0, "write synthetic spectral file: " // trim(errmsg))
    call read_spectral_file(BIN, sd, ierr, errmsg)
    call check(ierr == 0, "read synthetic spectral file: " // trim(errmsg))

    call test_normalisation()
    call test_intensity_conservation()
    call test_width_change_same_file()
    call test_combinations()
    call test_negative_weight()
    call test_fermi_dirac()
    call test_write_map_file()
    call test_tags()
    call test_kticks()
    call test_truncated_file()

    call delete_file(BIN)
    call delete_file(TRUNC)
    call delete_file(MAPF)
    call delete_file(TICKF)
    call test_summary("test_specmap")

contains

    !> nk = 5, nband = 3, regiões total/top/bottom (top + bottom = total),
    !> caminho G-K-M, energias referidas a E_F.
    subroutine build_synthetic(s)
        type(spectral_data_t), intent(out) :: s
        integer                            :: ik
        s%nk = 5; s%nband = 3; s%nspin = 1; s%nregion = 3; s%nvertex = 3
        s%spin_mode = SPIN_NONE; s%method = METHOD_MULLIKEN; s%energy_ref = EREF_FERMI
        s%fermi_known = .true.; s%fermi_energy_ev = -4.2_dp
        s%fingerprint = 'synthetic0000000'
        s%source = 'test_specmap synthetic data'
        allocate(s%kfrac(3, 5), s%kdist(5), s%vertex_k(3), s%vertex_label(3))
        allocate(s%region_name(3), s%region_natoms(3))
        allocate(s%energy(3, 5, 1), s%weight(3, 5, 3, 1))
        s%kfrac = 0.0_dp
        do ik = 1, 5
            s%kfrac(1, ik) = 0.1_dp * real(ik - 1, dp)
            s%kdist(ik) = 0.25_dp * real(ik - 1, dp)
            s%energy(1, ik, 1) = -1.5_dp + 0.2_dp * real(ik, dp)
            s%energy(2, ik, 1) = -0.1_dp * real(ik, dp)
            s%energy(3, ik, 1) =  1.0_dp + 0.05_dp * real(ik, dp)
            ! top: dominada pela banda 1, bottom: dominada pela banda 3, banda 2 mista
            s%weight(1, ik, 2, 1) = 0.8_dp
            s%weight(2, ik, 2, 1) = 0.1_dp * real(ik, dp)
            s%weight(3, ik, 2, 1) = 0.1_dp
        end do
        s%weight(:, :, 3, 1) = 1.0_dp - s%weight(:, :, 2, 1)
        s%weight(:, :, 1, 1) = 1.0_dp
        s%vertex_k = [1, 3, 5]
        s%vertex_label(1) = 'G'; s%vertex_label(2) = 'K'; s%vertex_label(3) = 'M'
        s%region_name(1) = 'total'; s%region_name(2) = 'top'; s%region_name(3) = 'bottom'
        s%region_natoms = [4, 2, 2]
    end subroutine build_synthetic

    !> Integral trapezoidal de y sobre uma grade uniforme x.
    pure function trapz(x, y) result(r)
        real(dp), intent(in) :: x(:), y(:)
        real(dp)             :: r
        integer              :: n
        n = size(x)
        r = (x(2) - x(1)) * (sum(y) - 0.5_dp * (y(1) + y(n)))
    end function trapz

    !> 1. Normalização dos núcleos em uma janela finita [-X, X].
    !> Os núcleos são normalizados na reta real; em uma janela finita as
    !> integrais exatas são (2/pi) atan(X/eta) (Lorentziana) e
    !> erf(X/(sigma sqrt 2)) (Gaussiana). Comparar com esses valores truncados
    !> isola a quadratura: a cauda ausente é propriedade da janela finita,
    !> não um erro do método.
    subroutine test_normalisation()
        real(dp), allocatable :: x(:)
        real(dp)              :: eta, sigma, xw
        integer               :: n
        eta = 0.05_dp; sigma = 0.05_dp
        xw = 3.0_dp; n = 60001
        allocate(x(n))
        call energy_grid(-xw, xw, n, x)
        call check_close(trapz(x, lorentzian(x, eta)), (2.0_dp/PI) * atan(xw/eta), 1.0e-4_dp, &
                         "Lorentzian truncated normalisation X=3")
        call check_close(trapz(x, gaussian(x, sigma)), erf(xw/(sigma*sqrt(2.0_dp))), 1.0e-6_dp, &
                         "Gaussian truncated normalisation X=3")
        xw = 0.1_dp
        call energy_grid(-xw, xw, n, x)
        call check_close(trapz(x, lorentzian(x, eta)), (2.0_dp/PI) * atan(xw/eta), 1.0e-4_dp, &
                         "Lorentzian truncated normalisation X=0.1")
        call check_close(trapz(x, gaussian(x, sigma)), erf(xw/(sigma*sqrt(2.0_dp))), 1.0e-6_dp, &
                         "Gaussian truncated normalisation X=0.1")
        call check_close(lorentzian(0.0_dp, eta), 1.0_dp/(PI*eta), 1.0e-12_dp, "Lorentzian peak 1/(pi eta)")
        call check_close(lorentzian(eta, eta), 0.5_dp/(PI*eta), 1.0e-12_dp, "Lorentzian eta is HWHM")
        call check_close(gaussian(sigma, sigma)/gaussian(0.0_dp, sigma), exp(-0.5_dp), 1.0e-14_dp, &
                         "Gaussian width is sigma")
    end subroutine test_normalisation

    !> 2. Uma banda de peso w em eps = 0 -> integral de A = w * integral truncada.
    subroutine test_intensity_conservation()
        real(dp) :: en(1,1), w(1,1), x(20001), A(20001,1), wt, xw
        xw = 4.0_dp
        wt = 0.7_dp
        en = 0.0_dp; w = wt
        call energy_grid(-xw, xw, size(x), x)
        call broaden(en, w, x, KERNEL_LORENTZIAN, 0.05_dp, A)
        call check_close(trapz(x, A(:,1)), wt * (2.0_dp/PI) * atan(xw/0.05_dp), 1.0e-4_dp, &
                         "Lorentzian intensity conservation")
        call broaden(en, w, x, KERNEL_GAUSSIAN, 0.05_dp, A)
        call check_close(trapz(x, A(:,1)), wt * erf(xw/(0.05_dp*sqrt(2.0_dp))), 1.0e-6_dp, &
                         "Gaussian intensity conservation")
    end subroutine test_intensity_conservation

    !> 3. Mudanças de largura reutilizam o mesmo arquivo (já lido): sem rediagonalização.
    subroutine test_width_change_same_file()
        real(dp), allocatable :: x(:), A1(:,:), A2(:,:)
        real(dp)              :: emin, emax, ref, i1, i2
        integer               :: ik, ib, ne
        emin = -12.0_dp; emax = 12.0_dp; ne = 24001
        allocate(x(ne), A1(ne, sd%nk), A2(ne, sd%nk))
        call energy_grid(emin, emax, ne, x)

        ! Gaussiana: janela larga o bastante para que o truncamento seja desprezível.
        call broaden(sd%energy(:,:,1), sd%weight(:,:,2,1), x, KERNEL_GAUSSIAN, 0.05_dp, A1)
        call broaden(sd%energy(:,:,1), sd%weight(:,:,2,1), x, KERNEL_GAUSSIAN, 0.10_dp, A2)
        i1 = 0; i2 = 0
        do ik = 1, sd%nk
            i1 = i1 + trapz(x, A1(:,ik)); i2 = i2 + trapz(x, A2(:,ik))
        end do
        call check_close(i1, sum(sd%weight(:,:,2,1)), 1.0e-8_dp, "Gaussian sigma=0.05 integrates to weight sum")
        call check_close(i2, sum(sd%weight(:,:,2,1)), 1.0e-8_dp, "Gaussian sigma=0.10 integrates to weight sum")
        call check(maxval(abs(A1 - A2)) > 1.0e-2_dp, "different widths give different maps")

        ! Lorentziana: compara com a soma analítica de pesos truncada.
        call broaden(sd%energy(:,:,1), sd%weight(:,:,2,1), x, KERNEL_LORENTZIAN, 0.05_dp, A1)
        call broaden(sd%energy(:,:,1), sd%weight(:,:,2,1), x, KERNEL_LORENTZIAN, 0.10_dp, A2)
        i1 = 0; i2 = 0
        do ik = 1, sd%nk
            i1 = i1 + trapz(x, A1(:,ik)); i2 = i2 + trapz(x, A2(:,ik))
        end do
        ref = 0
        do ik = 1, sd%nk
            do ib = 1, sd%nband
                ref = ref + sd%weight(ib,ik,2,1) * (atan((emax - sd%energy(ib,ik,1))/0.05_dp) - &
                                                    atan((emin - sd%energy(ib,ik,1))/0.05_dp)) / PI
            end do
        end do
        call check_close(i1, ref, 1.0e-4_dp, "Lorentzian eta=0.05 integrates to truncated weight sum")
        ref = 0
        do ik = 1, sd%nk
            do ib = 1, sd%nband
                ref = ref + sd%weight(ib,ik,2,1) * (atan((emax - sd%energy(ib,ik,1))/0.10_dp) - &
                                                    atan((emin - sd%energy(ib,ik,1))/0.10_dp)) / PI
            end do
        end do
        call check_close(i2, ref, 1.0e-4_dp, "Lorentzian eta=0.10 integrates to truncated weight sum")
    end subroutine test_width_change_same_file

    !> Mapa alargado de uma expressão analisada (combina os pesos primeiro, depois alarga).
    subroutine expr_map(expr, x, A, ierr)
        character(len=*), intent(in)  :: expr
        real(dp),         intent(in)  :: x(:)
        real(dp),         intent(out) :: A(:,:)
        integer,          intent(out) :: ierr
        integer,  allocatable         :: idx(:)
        real(dp), allocatable         :: c(:), w(:,:)
        character(len=512)            :: msg
        integer                       :: i
        call parse_map_expression(expr, sd, idx, c, ierr, msg)
        A = 0
        if (ierr /= 0) return
        allocate(w(sd%nband, sd%nk))
        w = 0
        do i = 1, size(idx)
            w = w + c(i) * sd%weight(:,:,idx(i),1)
        end do
        call broaden(sd%energy(:,:,1), w, x, KERNEL_LORENTZIAN, 0.05_dp, A)
    end subroutine expr_map

    !> 4. Mapas de diferença/soma e erros de análise.
    subroutine test_combinations()
        real(dp)              :: x(601), At(601,5), Ab(601,5), Ad(601,5), As(601,5), Atot(601,5)
        integer               :: ierr
        integer,  allocatable :: idx(:)
        real(dp), allocatable :: c(:)
        character(len=512)    :: msg
        real(dp)              :: tol

        call energy_grid(-3.0_dp, 3.0_dp, 601, x)
        call broaden(sd%energy(:,:,1), sd%weight(:,:,2,1), x, KERNEL_LORENTZIAN, 0.05_dp, At)
        call broaden(sd%energy(:,:,1), sd%weight(:,:,3,1), x, KERNEL_LORENTZIAN, 0.05_dp, Ab)
        call broaden(sd%energy(:,:,1), sd%weight(:,:,1,1), x, KERNEL_LORENTZIAN, 0.05_dp, Atot)
        ! Tolerância pontual de 1e-14 escalada pela magnitude do mapa (combinar antes
        ! ou depois do alargamento difere apenas por arredondamento de ponto flutuante).
        tol = 1.0e-14_dp * max(1.0_dp, maxval(abs(Atot)))

        call expr_map('top-bottom', x, Ad, ierr)
        call check(ierr == 0, "parse 'top-bottom'")
        call check(maxval(abs(Ad - (At - Ab))) <= tol, "difference map = A_top - A_bottom pointwise")
        call check(minval(Ad) < -1.0_dp, "difference map keeps negative values")
        call check(maxval(Ad) > 1.0_dp, "difference map keeps positive values")

        call expr_map(' top + bottom ', x, As, ierr)
        call check(ierr == 0, "parse ' top + bottom '")
        call check(maxval(abs(As - (At + Ab))) <= tol, "sum map = A_top + A_bottom pointwise")
        call check(maxval(abs(As - Atot)) <= tol, "top+bottom = total for a partition")

        call parse_map_expression('-top+bottom', sd, idx, c, ierr, msg)
        call check(ierr == 0 .and. size(idx) == 2, "parse leading '-'")
        if (ierr == 0 .and. size(idx) == 2) then
            call check(idx(1) == 2 .and. idx(2) == 3, "leading '-' indices")
            call check(c(1) == -1.0_dp .and. c(2) == 1.0_dp, "leading '-' coefficients")
        end if
        call parse_map_expression('top+top', sd, idx, c, ierr, msg)
        call check(ierr == 0 .and. size(idx) == 1, "repeated region combined")
        if (ierr == 0 .and. size(idx) == 1) call check(c(1) == 2.0_dp, "repeated region coefficient 2")

        call parse_map_expression('top-nonexistent', sd, idx, c, ierr, msg)
        call check(ierr /= 0, "unknown region rejected")
        call check(index(msg, 'nonexistent') > 0 .and. index(msg, 'bottom') > 0, &
                   "unknown region message lists available names")
        call parse_map_expression('top+', sd, idx, c, ierr, msg)
        call check(ierr /= 0, "empty trailing term 'top+' rejected")
        call parse_map_expression('++', sd, idx, c, ierr, msg)
        call check(ierr /= 0, "'++' rejected")
        call parse_map_expression('top--bottom', sd, idx, c, ierr, msg)
        call check(ierr /= 0, "'top--bottom' rejected")
        call parse_map_expression('+top', sd, idx, c, ierr, msg)
        call check(ierr /= 0, "leading '+' (empty first term) rejected")
        call parse_map_expression('   ', sd, idx, c, ierr, msg)
        call check(ierr /= 0, "empty expression rejected")
    end subroutine test_combinations

    !> 5. Pesos negativos são preservados (sem abs, sem recorte).
    subroutine test_negative_weight()
        real(dp) :: en(1,1), w(1,1), x(201), A(201,1)
        en = 0.0_dp; w = -0.2_dp
        call energy_grid(-1.0_dp, 1.0_dp, 201, x)   ! x(101) = 0
        call broaden(en, w, x, KERNEL_LORENTZIAN, 0.05_dp, A)
        call check(A(101,1) < 0.0_dp, "negative weight -> negative Lorentzian peak")
        call check_close(A(101,1), -0.2_dp/(PI*0.05_dp), 1.0e-12_dp, "negative Lorentzian peak value")
        call broaden(en, w, x, KERNEL_GAUSSIAN, 0.05_dp, A)
        call check(A(101,1) < 0.0_dp, "negative weight -> negative Gaussian peak")
    end subroutine test_negative_weight

    !> 6. Ocupação de Fermi-Dirac.
    subroutine test_fermi_dirac()
        real(dp) :: f, x
        integer  :: i
        logical  :: ok
        call check_close(fermi_dirac(0.3_dp, 0.3_dp, 300.0_dp), 0.5_dp, 1.0e-15_dp, "f(mu) = 1/2 at T>0")
        call check_close(fermi_dirac(0.3_dp, 0.3_dp, 0.0_dp), 0.5_dp, 0.0_dp, "f(mu) = 1/2 at T=0")
        call check_close(fermi_dirac(0.2_dp, 0.3_dp, 0.0_dp), 1.0_dp, 0.0_dp, "T=0 step below mu")
        call check_close(fermi_dirac(0.4_dp, 0.3_dp, 0.0_dp), 0.0_dp, 0.0_dp, "T=0 step above mu")
        f = fermi_dirac(10.0_dp, 0.0_dp, 1.0_dp)
        call check(f == f .and. f >= 0.0_dp .and. f < 1.0e-300_dp, "huge positive z: f -> 0, no NaN")
        f = fermi_dirac(-10.0_dp, 0.0_dp, 1.0_dp)
        call check(f == f .and. f == 1.0_dp, "huge negative z: f -> 1, no NaN")
        ok = .true.
        do i = 1, 50
            x = 0.013_dp * real(i, dp)
            f = fermi_dirac(0.1_dp + x, 0.1_dp, 300.0_dp) + fermi_dirac(0.1_dp - x, 0.1_dp, 300.0_dp)
            if (abs(f - 1.0_dp) > 1.0e-14_dp) ok = .false.
        end do
        call check(ok, "f(mu+x) + f(mu-x) = 1")
        call check_close(fermi_dirac(0.025852_dp, 0.0_dp, 300.0_dp), 1.0_dp/(1.0_dp + exp(1.0_dp)), &
                         1.0e-4_dp, "f at kT ~ 1/(1+e)")
    end subroutine test_fermi_dirac

    !> 7. Layout de write_map_file: nk*ne linhas de dados, nk linhas em branco (uma após
    !> cada bloco, inclusive o último), linhas de cabeçalho começam com '#'.
    subroutine test_write_map_file()
        real(dp)           :: x(11), A(11,5)
        character(len=64)  :: hdr(2)
        character(len=256) :: line
        integer            :: u, ios, ndata, nblank, nhash, nother, ierr
        character(len=512) :: msg
        call energy_grid(-1.0_dp, 1.0_dp, 11, x)
        call broaden(sd%energy(:,:,1), sd%weight(:,:,1,1), x, KERNEL_GAUSSIAN, 0.1_dp, A)
        hdr(1) = 'test header line'
        hdr(2) = '# already commented'
        call write_map_file(MAPF, sd, A, x, hdr, ierr, msg)
        call check(ierr == 0, "write_map_file: " // trim(msg))
        ndata = 0; nblank = 0; nhash = 0; nother = 0
        open(newunit=u, file=MAPF, status='old', action='read')
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) then
                nblank = nblank + 1
            else if (line(1:1) == '#') then
                nhash = nhash + 1
            else if (ndata + nblank == 0 .and. nhash == 0) then
                nother = nother + 1
            else
                ndata = ndata + 1
            end if
        end do
        close(u)
        call check(ndata == 5*11, "map file has nk*ne data lines")
        call check(nblank == 5, "map file has nk blank lines")
        call check(nhash >= 2 .and. nother == 0, "map header lines start with '#'")
    end subroutine test_write_map_file

    !> 8. Marcadores (tags) de arquivo.
    subroutine test_tags()
        call check(map_file_tag('top-bottom') == 'top_minus_bottom', "tag top-bottom")
        call check(map_file_tag('a+b') == 'a_plus_b', "tag a+b")
        call check(map_file_tag(' a + b - c ') == 'a_plus_b_minus_c', "tag with spaces")
        call check(map_file_tag('-top') == 'minus_top', "tag leading '-'")
        call check(map_file_tag('total') == 'total', "tag plain name")
    end subroutine test_tags

    subroutine test_kticks()
        integer            :: ierr, u, ios
        character(len=512) :: msg, line
        logical            :: has_xtics, has_ef, has_gamma, has_range, has_label
        call write_kticks_file(TICKF, sd, ierr, msg)
        call check(ierr == 0, "write_kticks_file: " // trim(msg))
        has_xtics = .false.; has_ef = .false.; has_gamma = .false.; has_range = .false.
        has_label = .false.
        open(newunit=u, file=TICKF, status='old', action='read')
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'set xtics (') == 1) has_xtics = .true.
            if (index(line, '{/Symbol G}') > 0) has_gamma = .true.
            if (index(line, 'set xrange [0:') == 1) has_range = .true.
            if (trim(line) == 'EF = 0.0') has_ef = .true.
            if (trim(line) == 'ELABEL = "E - E_F (eV)"') has_label = .true.
        end do
        close(u)
        call check(has_xtics .and. has_gamma .and. has_range, "kticks: xtics with Gamma and xrange")
        call check(has_ef .and. has_label, "kticks: EF = 0 and ELABEL for fermi reference")
    end subroutine test_kticks

    !> 9. Um arquivo espectral truncado é rejeitado (caminho de erro de leitura de specmap).
    subroutine test_truncated_file()
        character(len=:), allocatable :: buf
        integer                       :: u, fsize, ierr
        character(len=512)            :: msg
        type(spectral_data_t)         :: st
        open(newunit=u, file=BIN, access='stream', form='unformatted', status='old', action='read')
        inquire(unit=u, size=fsize)
        allocate(character(len=fsize/2) :: buf)
        read(u) buf
        close(u)
        open(newunit=u, file=TRUNC, access='stream', form='unformatted', status='replace', action='write')
        write(u) buf
        close(u)
        call read_spectral_file(TRUNC, st, ierr, msg)
        call check(ierr /= 0, "truncated spectral file rejected")
        call check(.not. allocated(st%energy) .and. .not. allocated(st%kfrac) .and. st%nk == 0, &
                   "failed read leaves sd reset")
    end subroutine test_truncated_file

    subroutine delete_file(fn)
        character(len=*), intent(in) :: fn
        integer                      :: u, ios
        open(newunit=u, file=fn, status='old', iostat=ios)
        if (ios == 0) close(u, status='delete')
    end subroutine delete_file

end program test_specmap
