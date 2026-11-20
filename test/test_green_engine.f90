! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Portões 2 e 3 no nível do programa (green_engine + green_data), com modelos
!> sintéticos gravados em .hsio e arquivos de estrutura com a_3 inclinado:
!>  - modo surface (&KPATH): soma das regiões = total, causalidade, resíduo,
!>    ida e volta do .green.bin;
!>  - modo film com contato 'none': a camada do substrato coincide com o modo
!>    surface (montagem consistente) e o filme se desloca exatamente por shift_ev;
!>  - teste A/A (sandwich e junction, contato 'scaled', lambda = 1): todas as
!>    camadas centrais têm a mesma DOS (bulk) a 1e-8;
!>  - &KGRID com &COHERENT: I_coh >= 0, I_incoh gravado;
!>  - erros: modelo 'extracted' (fase 4), pl_cells insuficiente, rede no plano
!>    incompatível; E_F do dataset meta do .hsio vira energy_reference.
program test_green_engine
    use types, only: dp
    use HS_io, only: HSDataFormat
    use structure_io, only: structure_t
    use green_engine, only: green_input_t, read_green_input, run_green
    use green_data, only: green_data_t, write_green_file, read_green_file, meta_get
    use green_models, only: make_model, make_structure, slab_eigen
    use testing, only: check, check_close, test_summary
    implicit none

    type(HSDataFormat) :: hh, hs
    type(structure_t) :: st
    real(dp) :: lat(3, 3), latb(3, 3)
    type(green_data_t) :: g1, g2, g3
    integer :: ierr, nbad
    character(len=:), allocatable :: errmsg
    character(len=*), parameter :: KPATH_BLOCK = 'path', KGRID_BLOCK = 'grid'

    lat(1, :) = [3.0_dp, 0.0_dp, 0.0_dp]
    lat(2, :) = [-1.5_dp, 2.598076211353316_dp, 0.0_dp]
    lat(3, :) = [1.0_dp, 0.6_dp, 5.0_dp]          ! a_3 inclinado

    ! Material A (substrato) e B (filme), mesma rede no plano
    call make_model([2, 1], [1, 1, 1], 0.3_dp, 0.02_dp, .true., 101, hh, hs)
    hh%has_meta = .true.
    hh%fermi_ev = -4.25_dp
    hh%energy_reference = 'absolute'
    call hh%WriteFile('tge_A_H.hsio')
    call hs%WriteFile('tge_A_S.hsio')
    call make_structure(lat, 2, 'Sa', 3, st)
    call write_structure('tge_A.nml', st)
    call make_model([2, 1], [1, 1, 1], 0.25_dp, 0.02_dp, .true., 202, hh, hs)
    call hh%WriteFile('tge_B_H.hsio')
    call hs%WriteFile('tge_B_S.hsio')
    latb = lat
    latb(3, :) = [0.4_dp, -0.3_dp, 4.0_dp]
    call make_structure(latb, 2, 'Fb', 9, st)
    call write_structure('tge_B.nml', st)
    ! segundo vizinho ao longo de a_3 (para o teste de pl_cells)
    call make_model([1], [1, 0, 2], 0.3_dp, 0.0_dp, .true., 303, hh, hs)
    call hh%WriteFile('tge_C_H.hsio')
    call make_structure(lat, 1, 'Cc', 4, st)
    call write_structure('tge_C.nml', st)
    ! rede no plano diferente
    latb(1, 1) = 3.3_dp
    call make_structure(latb, 2, 'Fb', 9, st)
    call write_structure('tge_Bbad.nml', st)

    ! --- 1. modo surface com &KPATH ---
    call write_input('tge_surface.nml', mode='surface', central=2, kblock=KPATH_BLOCK)
    call run('tge_surface.nml', g1)
    call check(ierr == 0, 'surface run: ' // errmsg)
    if (ierr == 0) then
        call check(g1%nreg == 2 .and. g1%nk == 7 .and. g1%ne == 9, 'surface: dimensions')
        call check_close(maxval(abs(sum(g1%amap, 3) - g1%total)), 0.0_dp, 1.0e-12_dp, 'surface: regions sum to total')
        call check(minval(g1%min_neg_im_gii) >= 0.0_dp, 'surface: -Im G_ii >= 0')
        call check(maxval(g1%residual) < 1.0e-9_dp, 'surface: decimation residual')
        call check(minval(g1%niter) > 0, 'surface: decimation converged everywhere')
        call check(minval(g1%smin) > 0.0_dp, 'surface: S_C positive definite')
        call check(meta_get(g1, 'energy_reference') == 'absolute', 'energy reference from .hsio meta')
        call check(index(meta_get(g1, 'substrate_fermi_ev'), '-4.25') > 0, 'substrate E_F recorded')
        call check(g1%kdist(g1%nk) > 0.0_dp, 'surface: k-path distances')
        call write_green_file('tge_surface.green.bin', g1, ierr, errmsg)
        call check(ierr == 0, 'write_green_file')
        call read_green_file('tge_surface.green.bin', g2, ierr, errmsg)
        call check(ierr == 0, 'read_green_file')
        call check(all(g2%total == g1%total) .and. all(g2%amap == g1%amap) .and. all(g2%niter == g1%niter) &
                   .and. all(g2%energy == g1%energy) .and. all(g2%kfrac == g1%kfrac) .and. &
                   all(g2%region == g1%region) .and. g2%nmeta == g1%nmeta .and. &
                   all(g2%label == g1%label), 'green.bin round trip is exact')
    end if

    ! --- 2. modo film, contato none: substrato idêntico ao modo surface (1 PL) ---
    call write_input('tge_surf1.nml', mode='surface', central=1, kblock=KPATH_BLOCK)
    call run('tge_surf1.nml', g1)
    call write_input('tge_film.nml', mode='film', central=1, nlayers=2, model='none', kblock=KPATH_BLOCK)
    call run('tge_film.nml', g2)
    call check(ierr == 0, 'film run: ' // errmsg)
    if (ierr == 0) then
        call check(g2%nreg == 3, 'film: one region per layer')
        call check_close(maxval(abs(g2%amap(:, :, 1) - g1%total)), 0.0_dp, 1.0e-12_dp, &
                         'film/none: substrate layer = bare surface')
        ! filme deslocado: DOS(E; Delta) = DOS(E - Delta; 0)
        call write_input('tge_film_shift.nml', mode='film', central=1, nlayers=2, model='none', &
                         kblock=KPATH_BLOCK, film_shift=0.3_dp, eshift=0.3_dp)
        call run('tge_film_shift.nml', g3)
        call check_close(maxval(abs(g3%amap(:, :, 2:3) - g2%amap(:, :, 2:3))), 0.0_dp, 1.0e-11_dp, &
                         'film shift_ev moves the film spectrum by Delta')
    end if

    ! --- 3. teste A/A ---
    call write_input('tge_aa.nml', mode='sandwich', central=1, nlayers=2, model='scaled', &
                     kblock=KPATH_BLOCK, film_is_substrate=.true.)
    call run('tge_aa.nml', g1)
    call check(ierr == 0, 'A/A sandwich run: ' // errmsg)
    if (ierr == 0) then
        call check(g1%nreg == 4, 'A/A: four layers')
        call check_close(maxval(abs(g1%amap(:, :, 2) - g1%amap(:, :, 1))) + &
                         maxval(abs(g1%amap(:, :, 3) - g1%amap(:, :, 1))) + &
                         maxval(abs(g1%amap(:, :, 4) - g1%amap(:, :, 1))), 0.0_dp, 1.0e-8_dp, &
                         'A/A sandwich: every layer has the bulk DOS')
        call write_input('tge_aa_j.nml', mode='junction', central=1, nlayers=1, model='scaled', &
                         kblock=KPATH_BLOCK, film_is_substrate=.true.)
        call run('tge_aa_j.nml', g2)
        call check(ierr == 0, 'A/A junction run')
        call check_close(maxval(abs(g2%amap(:, :, 1) - g1%amap(:, :, 1))), 0.0_dp, 1.0e-8_dp, &
                         'A/A junction = sandwich bulk DOS')
    end if

    ! --- 4. &KGRID + &COHERENT ---
    call write_input('tge_cohe.nml', mode='film', central=1, nlayers=2, model='scaled', &
                     kblock=KGRID_BLOCK, coherent=.true.)
    call run('tge_cohe.nml', g1)
    call check(ierr == 0, 'k-grid coherent run: ' // errmsg)
    if (ierr == 0) then
        call check(g1%has_coherent .and. g1%nq == 3 .and. g1%nk == 9, 'coherent: dimensions')
        call check(minval(g1%icoh) >= -1.0e-12_dp, 'coherent: I_coh >= 0')
        call check(size(g1%iinc) == size(g1%icoh) .and. maxval(abs(g1%iinc)) > 0.0_dp, 'incoherent map recorded')
        call check(maxval(abs(g1%icoh(:, :, 1) - g1%icoh(:, :, 3))) > 1.0e-6_dp, 'coherent map depends on q_z')
        call write_green_file('tge_cohe.green.bin', g1, ierr, errmsg)
        call read_green_file('tge_cohe.green.bin', g2, ierr, errmsg)
        call check(ierr == 0 .and. g2%has_coherent .and. all(g2%icoh == g1%icoh), 'coherent round trip')
    end if

    ! --- 5. Löwdin ---
    call write_input('tge_low.nml', mode='film', central=1, nlayers=1, model='scaled', &
                     kblock=KPATH_BLOCK, projection='lowdin')
    call run('tge_low.nml', g1)
    call check(ierr == 0, 'Lowdin run')
    if (ierr == 0) call check(minval(g1%amap) >= 0.0_dp, 'Lowdin region DOS non-negative')

    ! --- 6. erros ---
    call write_input('tge_err1.nml', mode='film', central=1, nlayers=1, model='extracted', kblock=KPATH_BLOCK)
    call run('tge_err1.nml', g1)
    call check(ierr /= 0 .and. index(errmsg, 'phase 4') > 0, "model 'extracted' reports phase 4")
    call write_input('tge_err2.nml', mode='surface', central=1, kblock=KPATH_BLOCK, sub='C', pl_cells=1)
    call run('tge_err2.nml', g1)
    call check(ierr /= 0 .and. index(errmsg, 'pl_cells') > 0, 'insufficient pl_cells is an error')
    call write_input('tge_ok2.nml', mode='surface', central=1, kblock=KPATH_BLOCK, sub='C', pl_cells=0)
    call run('tge_ok2.nml', g1)
    call check(ierr == 0 .and. index(meta_get(g1, 'substrate_pl'), '2 x 1') > 0, 'automatic pl_cells = 2')
    call write_input('tge_err3.nml', mode='film', central=1, nlayers=1, model='none', kblock=KPATH_BLOCK, &
                     film_struct='tge_Bbad.nml')
    call run('tge_err3.nml', g1)
    call check(ierr /= 0 .and. index(errmsg, 'strained') > 0, 'in-plane lattice mismatch is an error')

    ! --- 7. filme com número de células que não é múltiplo da PL (&FILM ncells) ---
    call odd_film_tests()

    call cleanup()
    call test_summary('test_green_engine')

contains

    !> Filme D com segundo vizinho ao longo de a_3 (PL = 2 células) sobre A:
    !>  - ncells = 4 é idêntico a nlayers = 2;
    !>  - ncells = 3 com contato 'none': a DOS do filme em Γ é a soma de
    !>    lorentzianas nos autovalores do slab finito de 3 células, montado
    !>    diretamente de H(R), S(R);
    !>  - ncells fora do modo film / side top é rejeitado.
    subroutine odd_film_tests()
        type(HSDataFormat) :: dh, ds
        type(structure_t) :: dst
        type(green_data_t) :: ga, gb
        complex(dp), allocatable :: hsl(:, :), ssl(:, :), c(:, :)
        real(dp), allocatable :: ev(:)
        real(dp) :: dev, ref, e, eta
        integer :: nb, ir, a, b, ie, ig, j
        call make_model([1, 1], [1, 1, 2], 0.25_dp, 0.02_dp, .true., 404, dh, ds)
        call dh%WriteFile('tge_D_H.hsio')
        call ds%WriteFile('tge_D_S.hsio')
        call make_structure(lat, 2, 'Dd', 12, dst)
        call write_structure('tge_D.nml', dst)
        call write_odd('tge_odd_pl.nml', 'nlayers = 2', 'none')
        call run('tge_odd_pl.nml', ga)
        call check(ierr == 0, 'film nlayers = 2 with P = 2: ' // errmsg)
        call write_odd('tge_odd_c4.nml', 'ncells = 4', 'none')
        call run('tge_odd_c4.nml', gb)
        call check(ierr == 0, 'film ncells = 4: ' // errmsg)
        if (ierr == 0) then
            call check(all(gb%total == ga%total) .and. all(gb%amap == ga%amap), 'ncells = 4 == nlayers = 2 (P = 2)')
            call check(meta_get(gb, 'film_cells') == '4', 'film_cells metadata')
        end if
        call write_odd('tge_odd_c3.nml', 'ncells = 3', 'none')
        call run('tge_odd_c3.nml', gb)
        call check(ierr == 0, 'film ncells = 3: ' // errmsg)
        if (ierr == 0) then
            call check(gb%nreg == 3 .and. meta_get(gb, 'film_cells') == '3', 'ncells = 3: one full PL + one partial block')
            ! slab finito de 3 células em Γ
            nb = dh%NOrbs
            allocate(hsl(3 * nb, 3 * nb), ssl(3 * nb, 3 * nb), c(3 * nb, 3 * nb), ev(3 * nb))
            hsl = (0.0_dp, 0.0_dp); ssl = (0.0_dp, 0.0_dp)
            do ir = 1, dh%Neighbors
                do a = 0, 2
                    b = a + dh%IndexList(3, ir)
                    if (b < 0 .or. b > 2) cycle
                    hsl(a * nb + 1:(a + 1) * nb, b * nb + 1:(b + 1) * nb) = &
                        hsl(a * nb + 1:(a + 1) * nb, b * nb + 1:(b + 1) * nb) + dh%Data(:, :, ir)
                    ssl(a * nb + 1:(a + 1) * nb, b * nb + 1:(b + 1) * nb) = &
                        ssl(a * nb + 1:(a + 1) * nb, b * nb + 1:(b + 1) * nb) + ds%Data(:, :, ir)
                end do
            end do
            call slab_eigen(hsl, ssl, ev, c)
            ig = minloc(norm2(gb%kcart, dim=1), 1)
            eta = 0.05_dp
            dev = 0.0_dp
            do ie = 1, gb%ne
                e = gb%energy(ie)
                ref = 0.0_dp
                do j = 1, 3 * nb
                    ref = ref + (eta / 3.1415926535897932384626433832795_dp) / ((e - ev(j))**2 + eta**2)
                end do
                dev = max(dev, abs(gb%amap(ie, ig, 2) + gb%amap(ie, ig, 3) - ref) / ref)
            end do
            call check_close(dev, 0.0_dp, 1.0e-9_dp, 'ncells = 3, decoupled: film DOS = 3-cell slab eigenvalues')
        end if
        call write_odd('tge_odd_bad.nml', 'ncells = 3', 'none', mode='junction')
        call run('tge_odd_bad.nml', gb)
        call check(ierr /= 0 .and. index(errmsg, 'ncells') > 0, "ncells outside mode 'film' is rejected")
    end subroutine odd_film_tests

    subroutine write_odd(fname, thick, model, mode)
        character(len=*), intent(in) :: fname, thick, model
        character(len=*), intent(in), optional :: mode
        integer :: u
        character(len=16) :: md
        md = 'film'
        if (present(mode)) md = mode
        open(newunit=u, file=fname, status='replace', action='write')
        write(u, '(A)') '&GREEN'
        write(u, '(A,A,A)') "  mode = '", trim(md), "', eta_ev = 0.05, emin = -1.2, emax = 1.2, ne = 9"
        write(u, '(A,A,A)') "  decim_tol = 1e-12, decim_maxit = 300, output_prefix = '", fname(1:len(fname) - 4), "'"
        write(u, '(A)') '/'
        write(u, '(A)') '&SUBSTRATE'
        write(u, '(A)') "  h_file = 'tge_A_H.hsio', s_file = 'tge_A_S.hsio', structure_file = 'tge_A.nml'"
        write(u, '(A)') "  stacking_axis = 3, central_layers = 1, side = 'top'"
        write(u, '(A)') '/'
        write(u, '(A)') '&FILM'
        write(u, '(A)') "  h_file = 'tge_D_H.hsio', s_file = 'tge_D_S.hsio', structure_file = 'tge_D.nml'"
        write(u, '(A,A)') '  ', thick
        write(u, '(A)') '/'
        write(u, '(A)') '&INTERFACE'
        write(u, '(A,A,A)') "  model = '", model, "'"
        write(u, '(A)') '/'
        write(u, '(A)') '&KPATH'
        write(u, '(A)') '  nvertex = 2'
        write(u, '(A)') '  vertex(:,1) = 0.0, 0.0, 0.0'
        write(u, '(A)') '  vertex(:,2) = 0.2, 0.0, 0.0'
        write(u, '(A)') "  label = 'G', 'X'"
        write(u, '(A)') '  npoints = 2'
        write(u, '(A)') '/'
        close(u)
    end subroutine write_odd

    subroutine run(nml, gd)
        character(len=*), intent(in) :: nml
        type(green_data_t), intent(out) :: gd
        type(green_input_t) :: inp
        call read_green_input(nml, inp, ierr, errmsg)
        if (ierr /= 0) return
        call run_green(inp, gd, ierr, errmsg, nonconverged=nbad)
    end subroutine run

    subroutine write_structure(fname, s)
        character(len=*),  intent(in) :: fname
        type(structure_t), intent(in) :: s
        integer :: u, i
        open(newunit=u, file=fname, status='replace', action='write')
        write(u, '(A)') '&STRUCTURE'
        do i = 1, 3
            write(u, '(A,I0,A,3(F20.14,:,","))') '  LatticeVectors(', i, ',:) = ', s%lattice(i, :)
        end do
        write(u, '(A)') '/'
        write(u, '(A)') '&ATOMS'
        write(u, '(I0)') s%nat
        do i = 1, s%nat
            write(u, '(A,3F20.14)') trim(s%label(i)), s%cart(:, i)
        end do
        write(u, '(A)') '/'
        close(u)
    end subroutine write_structure

    subroutine write_input(fname, mode, central, kblock, nlayers, model, film_shift, eshift, &
                           film_is_substrate, coherent, projection, sub, pl_cells, film_struct)
        character(len=*), intent(in) :: fname, mode, kblock
        integer, intent(in) :: central
        integer, intent(in), optional :: nlayers, pl_cells
        character(len=*), intent(in), optional :: model, projection, sub, film_struct
        real(dp), intent(in), optional :: film_shift, eshift
        logical, intent(in), optional :: film_is_substrate, coherent
        integer :: u
        real(dp) :: e0
        character(len=1) :: sm, fm
        sm = 'A'; fm = 'B'
        if (present(sub)) sm = sub(1:1)
        if (present(film_is_substrate)) then
            if (film_is_substrate) fm = 'A'
        end if
        e0 = 0.0_dp
        if (present(eshift)) e0 = eshift
        open(newunit=u, file=fname, status='replace', action='write')
        write(u, '(A)') '&GREEN'
        write(u, '(A,A,A)') "  mode = '", mode, "'"
        write(u, '(A)') '  eta_ev = 0.05'
        write(u, '(A,F12.6,A,F12.6,A)') '  emin = ', -1.2_dp + e0, ', emax = ', 1.2_dp + e0, ', ne = 9'
        if (kblock == KGRID_BLOCK) write(u, '(A)') '  ne = 2'
        write(u, '(A)') '  decim_tol = 1e-12, decim_maxit = 300'
        if (present(projection)) write(u, '(A,A,A)') "  projection = '", projection, "'"
        write(u, '(A,A,A)') "  output_prefix = '", fname(1:len(fname) - 4), "'"
        write(u, '(A)') '/'
        write(u, '(A)') '&SUBSTRATE'
        if (sm == 'C') then
            write(u, '(A)') "  h_file = 'tge_C_H.hsio', structure_file = 'tge_C.nml'"
        else
            write(u, '(A)') "  h_file = 'tge_A_H.hsio', s_file = 'tge_A_S.hsio', structure_file = 'tge_A.nml'"
        end if
        write(u, '(A,I0)') '  stacking_axis = 3, central_layers = ', central
        if (present(pl_cells)) write(u, '(A,I0)') '  pl_cells = ', pl_cells
        write(u, '(A)') "  side = 'top'"
        write(u, '(A)') '/'
        if (present(nlayers)) then
            write(u, '(A)') '&FILM'
            if (fm == 'A') then
                write(u, '(A)') "  h_file = 'tge_A_H.hsio', s_file = 'tge_A_S.hsio', structure_file = 'tge_A.nml'"
            else if (present(film_struct)) then
                write(u, '(A,A,A)') "  h_file = 'tge_B_H.hsio', s_file = 'tge_B_S.hsio', structure_file = '", &
                                    film_struct, "'"
            else
                write(u, '(A)') "  h_file = 'tge_B_H.hsio', s_file = 'tge_B_S.hsio', structure_file = 'tge_B.nml'"
            end if
            write(u, '(A,I0)') '  nlayers = ', nlayers
            if (present(film_shift)) write(u, '(A,F12.6)') '  shift_ev = ', film_shift
            write(u, '(A)') '/'
        end if
        if (present(model)) then
            write(u, '(A)') '&INTERFACE'
            write(u, '(A,A,A)') "  model = '", model, "', reference = 'substrate', lambda_h = 1.0, lambda_s = 1.0"
            write(u, '(A)') '/'
        end if
        if (present(coherent)) then
            write(u, '(A)') '&COHERENT'
            write(u, '(A)') '  nqz = 3, qz_min = 0.0, qz_max = 1.2, escape_length_ang = 5.0'
            write(u, '(A)') '/'
        end if
        if (kblock == KPATH_BLOCK) then
            write(u, '(A)') '&KPATH'
            write(u, '(A)') '  nvertex = 3'
            write(u, '(A)') '  vertex(:,1) = 0.0, 0.0, 0.0'
            write(u, '(A)') '  vertex(:,2) = 0.5, 0.0, 0.0'
            write(u, '(A)') '  vertex(:,3) = 0.3333333333, 0.3333333333, 0.0'
            write(u, '(A)') "  label = 'G', 'M', 'K'"
            write(u, '(A)') '  npoints = 3, 3'
            write(u, '(A)') '/'
        else
            write(u, '(A)') '&KGRID'
            write(u, '(A)') '  origin = -0.1, -0.1, 0.0'
            write(u, '(A)') '  span1 = 0.2, 0.0, 0.0'
            write(u, '(A)') '  span2 = 0.0, 0.2, 0.0'
            write(u, '(A)') '  n1 = 3, n2 = 3'
            write(u, '(A)') '/'
        end if
        close(u)
    end subroutine write_input

    subroutine cleanup()
        call execute_command_line('rm -f tge_*.hsio tge_*.nml tge_*.green.bin', wait=.true.)
    end subroutine cleanup

end program test_green_engine
