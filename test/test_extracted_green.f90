! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Modo 'extracted' do surfgreen (região central tirada de uma supercélula),
!> com um modelo sintético não ortogonal, complexo, a_3 inclinado:
!>  1. supercélula = o próprio lead, lead à esquerda e corte à direita: igual ao
!>     modo surface do surfgreen (camada a camada, 1e-10), também com 1 e 3
!>     colunas por bloco de RGF;
!>  2. lead dos dois lados (A/A): todas as colunas com a mesma DOS (1e-8) com os
!>     termos cruzados; sem eles a coluna da borda difere (o teste os exige);
!>     &LEAD_RIGHT explícito (duas decimações) = um só material (decimação dupla);
!>  3. supercélula explícita (6 células do lead ao longo de a_3, periódica,
!>     transladada de um t uniforme, átomos em ordem trocada): igual ao caso 2
!>     (1e-10) com match_translation = 'auto'; 'none' falha com erro claro;
!>  4. supercélula com H + Delta0 S: delta_mode = 'fit' recupera Delta0
!>     (1e-12) e o mapa coincide com o caso 3 deslocado em energia;
!>  5. erros: acoplamento através do vácuo, coluna da costura com átomo extra,
!>     central_columns ausente.
!>  0. surface_gf_pair = surface_gf nos dois lados (1e-12);
!>  7. &SITEPOT: V = c em todos os átomos centrais e na costura (lead_ev = c)
!>     com Delta = c desloca o espectro rigidamente, A(E) = A(E - c) (1e-10);
!>     o mesmo V lido de um arquivo (casamento por posição) = uniform_ev;
!>  8. &KLIST: média uniforme em k_1 de um ponto da lista = média dos mesmos
!>     pontos calculados por &KPATH (1e-12);
!>  9. &RHO num modelo real (simetria de inversão temporal) e k simétricos:
!>     a carga de Mulliken por átomo sum_j tr[D_ij S_ji] (todos os pares,
!>     inclusive com a costura) = integral do DOS de Mulliken da janela (1e-10);
!>  6. .hsio bloco-esparso: read_hs_blocksparse/write_hs_blocksparse
!>     reproduzem o arquivo gravado por HSDataFormat%WriteFile byte a byte.
program test_extracted_green
    use, intrinsic :: iso_fortran_env, only: int64
    use types, only: dp
    use HS_io, only: HSDataFormat, HSBlockSparse, read_hs_blocksparse, write_hs_blocksparse
    use structure_io, only: structure_t
    use green_engine, only: green_input_t, read_green_input, run_green
    use extracted_green, only: extracted_input_t, read_extracted_input, run_extracted, input_is_extracted
    use green_data, only: green_data_t, meta_get
    use green_models, only: make_model, make_structure, rnd, set_seed
    use surface_green, only: surface_gf, surface_gf_pair, SIDE_TOP, SIDE_BOTTOM
    use testing, only: check, check_close, test_summary
    implicit none

    type(HSDataFormat) :: hh, hs, sh, ss
    type(structure_t) :: st, sst
    type(green_data_t) :: gs, gx, gy, gz
    real(dp) :: lat(3, 3), t(3), d0, scale
    integer :: ierr, nbad, i, l
    character(len=512) :: line9
    character(len=:), allocatable :: errmsg
    integer, parameter :: NSUP = 6
    character(len=*), parameter :: LEAD_NML = &
        "&LEAD h_file = 'tx_L_H.hsio', s_file = 'tx_L_S.hsio', structure_file = 'tx_L.nml' /"

    lat(1, :) = [3.0_dp, 0.0_dp, 0.0_dp]
    lat(2, :) = [-1.5_dp, 2.598076211353316_dp, 0.0_dp]
    lat(3, :) = [1.0_dp, 0.6_dp, 5.0_dp]          ! a_3 inclinado

    call make_model([2, 1, 2], [1, 1, 1], 0.3_dp, 0.03_dp, .true., 707, hh, hs)
    hh%has_meta = .true.
    hh%fermi_ev = -3.5_dp
    hh%energy_reference = 'absolute'
    call hh%WriteFile('tx_L_H.hsio')
    call hs%WriteFile('tx_L_S.hsio')
    call make_structure(lat, 3, 'X', 11, st)
    ! espalha os átomos ao longo de a_3 sem cruzar a fronteira da célula
    do i = 1, 3
        st%cart(:, i) = st%cart(:, i) + (0.1_dp * real(i - 2, dp)) * lat(3, :)
    end do
    call write_structure('tx_L.nml', st)
    call check(input_is_extracted('tx_surface_x.nml') .eqv. .false., 'missing file is not extracted')

    ! --- 0. decimação dos dois lados numa só iteração = duas decimações ---
    block
        integer, parameter :: n = 5
        complex(dp) :: k00(n, n), k01(n, n), k10(n, n), g1(n, n), g2(n, n), p1(n, n), p2(n, n), h(n, n), x(n, n)
        complex(dp) :: z
        real(dp) :: r1, r2, q1, q2
        integer :: it1, it2, itp, a, b
        call set_seed(99)
        do b = 1, n
            do a = 1, n
                h(a, b) = cmplx(rnd(), rnd(), dp)
                x(a, b) = 0.4_dp * cmplx(rnd(), rnd(), dp)
            end do
        end do
        h = 0.5_dp * (h + conjg(transpose(h)))
        z = (0.3_dp, 0.05_dp)
        k00 = -h
        do a = 1, n
            k00(a, a) = k00(a, a) + z
        end do
        k01 = -x
        k10 = -conjg(transpose(x))
        call surface_gf(k00, k01, k10, SIDE_TOP, 1.0e-12_dp, 500, g1, it1, r1, ierr, errmsg)
        call surface_gf(k00, k01, k10, SIDE_BOTTOM, 1.0e-12_dp, 500, g2, it2, r2, ierr, errmsg)
        call surface_gf_pair(k00, k01, k10, 1.0e-12_dp, 500, p1, p2, itp, q1, q2, ierr, errmsg)
        call check(ierr == 0, 'surface_gf_pair converged: ' // errmsg)
        call check_close(maxval(abs(p1 - g1)) / maxval(abs(g1)), 0.0_dp, 1.0e-12_dp, 'pair decimation: top side')
        call check_close(maxval(abs(p2 - g2)) / maxval(abs(g2)), 0.0_dp, 1.0e-12_dp, 'pair decimation: bottom side')
        call check(itp == max(it1, it2) .and. q1 < 1.0e-10_dp .and. q2 < 1.0e-10_dp, 'pair decimation: iterations, residuals')
    end block

    ! --- 6. ida e volta bloco-esparsa ---
    block
        type(HSBlockSparse) :: b
        call read_hs_blocksparse('tx_L_H.hsio', b, ierr, errmsg)
        call check(ierr == 0, 'read_hs_blocksparse: ' // errmsg)
        call check(b%has_meta .and. abs(b%fermi_ev + 3.5_dp) < 1.0e-14_dp, 'block-sparse meta')
        call write_hs_blocksparse('tx_L_H_copy.hsio', b, 'H', ierr, errmsg)
        call check(ierr == 0, 'write_hs_blocksparse: ' // errmsg)
        call check(same_bytes('tx_L_H.hsio', 'tx_L_H_copy.hsio'), 'block-sparse round trip is byte-identical')
        call read_hs_blocksparse('missing.hsio', b, ierr, errmsg)
        call check(ierr == 1, 'read_hs_blocksparse: missing file gives ierr = 1')
    end block

    ! --- 1. lead como supercélula, corte à direita = modo surface ---
    call write_surface('tx_surface.nml', 4)
    call run_sg('tx_surface.nml', gs)
    call check(ierr == 0, 'surfgreen surface: ' // errmsg)
    call write_x('tx_x1.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 3, 'lead', 'cut')
    call run_x('tx_x1.nml', gx)
    call check(ierr == 0, 'extracted (lead as supercell, cut): ' // errmsg)
    if (ierr == 0 .and. allocated(gs%amap)) then
        call check(gx%nreg == 4 .and. gs%nreg == 4, 'one region per column / per layer')
        scale = maxval(abs(gs%amap))
        call check_close(maxval(abs(gx%amap - gs%amap)) / scale, 0.0_dp, 1.0e-10_dp, &
                         'extracted = surfgreen surface, layer by layer')
        call check_close(maxval(abs(gx%total - gs%total)) / maxval(abs(gs%total)), 0.0_dp, 1.0e-10_dp, &
                         'extracted = surfgreen surface, total')
        call check(minval(gx%min_neg_im_gii) >= 0.0_dp, '-Im G_ii >= 0')
        call check(maxval(gx%residual) < 1.0e-9_dp, 'decimation residual')
        call check(meta_get(gx, 'mode') == 'extracted', 'meta mode')
        call check(meta_get(gx, 'energy_reference') == 'absolute', 'energy reference from .hsio meta')
    end if

    ! --- 2. A/A: colunas iguais ao bulk ---
    call write_x('tx_x2.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 2, 'lead', 'lead')
    call run_x('tx_x2.nml', gy)
    call check(ierr == 0, 'extracted A/A: ' // errmsg)
    if (ierr == 0) then
        scale = maxval(abs(gy%amap))
        do l = 2, 3
            call check_close(maxval(abs(gy%amap(:, :, l) - gy%amap(:, :, 1))) / scale, 0.0_dp, 1.0e-8_dp, &
                             'A/A: column equals the first one')
        end do
        call check(minval(gy%niter) > 0, 'A/A: decimation converged (paired decimation)')
    end if
    call write_x('tx_x2n.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 2, 'lead', 'lead', cross='.false.')
    call run_x('tx_x2n.nml', gz)
    if (ierr == 0) call check(maxval(abs(gz%amap(:, :, 1) - gz%amap(:, :, 2))) / maxval(abs(gz%amap)) > 1.0e-6_dp, &
                              'A/A without lead cross terms: edge column differs (cross terms are needed)')

    ! --- 2b. &LEAD_RIGHT explícito (decimações separadas) = um só material ---
    call write_x('tx_x2r.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 2, 'lead', 'lead', &
                 extra="/" // new_line('a') // "&LEAD_RIGHT h_file = 'tx_L_H.hsio', s_file = 'tx_L_S.hsio', " // &
                       "structure_file = 'tx_L.nml'")
    call run_x('tx_x2r.nml', gz)
    call check(ierr == 0, '&LEAD_RIGHT run: ' // errmsg)
    if (ierr == 0 .and. allocated(gy%amap)) &
        call check_close(maxval(abs(gz%amap - gy%amap)) / maxval(abs(gy%amap)), 0.0_dp, 1.0e-11_dp, &
                         '&LEAD_RIGHT (two decimations) = one lead material (paired decimation)')

    ! --- 3. supercélula explícita, transladada e reordenada ---
    t = [0.37_dp, -0.21_dp, 0.83_dp]
    call make_supercell(0.0_dp, t, sh, ss, sst)
    call sh%WriteFile('tx_S_H.hsio')
    call ss%WriteFile('tx_S_S.hsio')
    call write_structure('tx_S.nml', sst)
    call write_x('tx_x3.nml', 'tx_S_H.hsio', 'tx_S_S.hsio', 'tx_S.nml', 2, 4, 'lead', 'lead', &
                 extra="match_translation = 'auto'")
    call run_x('tx_x3.nml', gz)
    call check(ierr == 0, 'explicit supercell: ' // errmsg)
    if (ierr == 0 .and. allocated(gy%amap)) then
        call check_close(maxval(abs(gz%amap - gy%amap)) / maxval(abs(gy%amap)), 0.0_dp, 1.0e-10_dp, &
                         'explicit translated supercell = lead as supercell')
    end if
    call write_x('tx_x3n.nml', 'tx_S_H.hsio', 'tx_S_S.hsio', 'tx_S.nml', 2, 4, 'lead', 'lead')
    call run_x('tx_x3n.nml', gx)
    call check(ierr /= 0 .and. index(errmsg, 'do not match') > 0, &
               "translated supercell with match_translation = 'none' is an error")

    ! --- 4. Delta ajustado ---
    d0 = 0.137_dp
    call make_supercell(d0, t, sh, ss, sst)
    call sh%WriteFile('tx_D_H.hsio')
    call write_x('tx_x4.nml', 'tx_D_H.hsio', 'tx_S_S.hsio', 'tx_S.nml', 2, 4, 'lead', 'lead', &
                 extra="match_translation = 'auto', delta_mode = 'fit', fit_neighbors = 1, " // &
                       "fit_k(:,1) = 0.0, 0.0, 0.0, fit_k(:,2) = 0.2, 0.1, 0.0", eshift=d0)
    call run_x('tx_x4.nml', gx)
    call check(ierr == 0, 'Delta fit run: ' // errmsg)
    if (ierr == 0 .and. allocated(gz%amap)) then
        block
            real(dp) :: dfit
            character(len=64) :: buf
            buf = meta_get(gx, 'delta_left_ev')
            read(buf, *) dfit
            call check_close(dfit, d0, 1.0e-12_dp, 'fitted Delta recovers H + Delta0 S')
        end block
        call check_close(maxval(abs(gx%amap - gz%amap)) / maxval(abs(gz%amap)), 0.0_dp, 1.0e-9_dp, &
                         'shifted supercell with fitted Delta = unshifted map at E + Delta0')
    end if

    ! --- 5. erros ---
    call write_x('tx_e1.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 3, 'lead', 'cut', extra='vacuum_axes = 1')
    call run_x('tx_e1.nml', gx)
    call check(ierr /= 0 .and. index(errmsg, 'vacuum') > 0, 'coupling through a vacuum axis is an error')
    call write_x('tx_e2.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 3, 'lead', 'vacuum')
    call run_x('tx_e2.nml', gx)
    call check(ierr /= 0 .and. index(errmsg, 'neither central nor') > 0, &
               "right = 'vacuum' with atoms beyond the central region is an error")
    open(newunit=i, file='tx_e3.nml', status='replace', action='write')
    write(i, '(A)') "&GREEN mode = 'extracted' /"
    write(i, '(A)') "&EXTRACTED h_file = 'tx_L_H.hsio', structure_file = 'tx_L.nml', embed_axis = 3 /"
    write(i, '(A)') LEAD_NML
    write(i, '(A)') "&KPATH nvertex = 2, vertex(:,1) = 0,0,0, vertex(:,2) = 0.5,0,0, npoints = 2 /"
    close(i)
    call run_x('tx_e3.nml', gx)
    call check(ierr /= 0 .and. index(errmsg, 'central_columns') > 0, 'missing central_columns is an error')

    ! --- 7. &SITEPOT: deslocamento rígido ---
    d0 = 0.21_dp
    call write_x('tx_x7.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 2, 'lead', 'lead', &
                 extra="delta_ev = 0.21 /" // new_line('a') // "&SITEPOT uniform_ev = 0.21, lead_ev = 0.21", eshift=d0)
    call run_x('tx_x7.nml', gx)
    call check(ierr == 0, '&SITEPOT uniform run: ' // errmsg)
    if (ierr == 0 .and. allocated(gy%amap)) &
        call check_close(maxval(abs(gx%amap - gy%amap)) / maxval(abs(gy%amap)), 0.0_dp, 1.0e-10_dp, &
                         '&SITEPOT V = c everywhere (1/2{S,V}) + Delta = c: rigid shift A(E) = A(E - c)')
    ! mesmo V por arquivo (escalar por átomo, posições das colunas 0..2 ao longo de a_3)
    open(newunit=i, file='tx_v.sitepot', status='replace', action='write')
    write(i, '(A)') 'LFSITEPOT 1'
    write(i, '(A)') 'norb 2'
    write(i, '(A)') 'nentry 9'
    do l = 0, 2
        do nbad = 1, 3
            write(i, '(A,I0,1X,I0,3(1X,ES24.16),A)') 'entry ', l, nbad, st%cart(:, nbad) + real(l, dp) * lat(3, :), &
                ' scalar 0.21'
        end do
    end do
    close(i)
    call write_x('tx_x7f.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 2, 'lead', 'lead', &
                 extra="delta_ev = 0.21 /" // new_line('a') // "&SITEPOT file = 'tx_v.sitepot', lead_ev = 0.21", &
                 eshift=d0)
    call run_x('tx_x7f.nml', gz)
    call check(ierr == 0, '&SITEPOT file run: ' // errmsg)
    if (ierr == 0 .and. allocated(gx%amap)) &
        call check_close(maxval(abs(gz%amap - gx%amap)) / maxval(abs(gx%amap)), 0.0_dp, 1.0e-12_dp, &
                         '&SITEPOT file (matched by position) = uniform_ev')
    call write_x('tx_x7o.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 2, 'lead', 'lead', &
                 extra="delta_ev = 0.21 /" // new_line('a') // "&SITEPOT uniform_ev = 0.21, mode = 'onsite'", eshift=d0)
    call run_x('tx_x7o.nml', gz)
    if (ierr == 0 .and. allocated(gy%amap)) &
        call check(maxval(abs(gz%amap - gy%amap)) / maxval(abs(gy%amap)) > 1.0e-6_dp, &
                   "&SITEPOT mode = 'onsite' (no off-diagonal V S) is not a rigid shift")

    ! --- 7b. &SITEPOT com bloco por átomo (A12): bloco zero = sem &SITEPOT bit a bit; bloco não nulo atua ---
    call write_x('tx_x7r.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 2, 'lead', 'lead')
    call run_x('tx_x7r.nml', gy)
    call write_blk('tx_z.sitepot', 0.0_dp)
    call write_x('tx_x7z.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 2, 'lead', 'lead', &
                 extra="/" // new_line('a') // "&SITEPOT file = 'tx_z.sitepot'")
    call run_x('tx_x7z.nml', gz)
    call check(ierr == 0, '&SITEPOT zero-block run: ' // errmsg)
    if (ierr == 0 .and. allocated(gy%amap) .and. allocated(gz%amap)) &
        call check(all(gz%amap == gy%amap), '&SITEPOT zero 2x2 blocks on the central atoms: A(k,E) bit-identical')
    call write_blk('tx_z.sitepot', 0.13_dp)
    call run_x('tx_x7z.nml', gz)
    if (ierr == 0 .and. allocated(gy%amap) .and. allocated(gz%amap)) &
        call check(maxval(abs(gz%amap - gy%amap)) / maxval(abs(gy%amap)) > 1.0e-6_dp, &
                   '&SITEPOT non-zero block changes A(k,E)')

    ! --- 8. &KLIST ---
    call write_x('tx_x8.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 2, 'lead', 'lead', &
                 kblock="&KLIST avg_axis = 1, navg = 2, list_axis = 2, nlist = 2, list_k = 0.1, 0.3, " // &
                        "list_w = 1.0, 3.0 /")
    call run_x('tx_x8.nml', gx)
    call check(ierr == 0, '&KLIST run: ' // errmsg)
    call write_x('tx_x8p.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 2, 'lead', 'lead', &
                 kblock="&KPATH nvertex = 2, vertex(:,1) = -0.5, 0.1, 0.0, vertex(:,2) = 0.0, 0.1, 0.0, npoints = 1 /")
    call run_x('tx_x8p.nml', gz)
    if (ierr == 0 .and. allocated(gx%amap)) then
        call check(gx%nk == 2 .and. abs(gx%kfrac(2, 2) - 0.3_dp) < 1.0e-15_dp, '&KLIST: one output k per list point')
        call check_close(maxval(abs(gx%amap(:, 1, :) - 0.5_dp * (gz%amap(:, 1, :) + gz%amap(:, 2, :)))) / &
                         maxval(abs(gz%amap)), 0.0_dp, 1.0e-12_dp, '&KLIST average over k_1 = mean of &KPATH points')
        call check(index(meta_get(gx, 'klist_weights'), '2.500000E-01') > 0, '&KLIST normalised weights in meta')
    end if

    ! --- 9. &RHO: Mulliken por D S = DOS integrado (modelo real, k simétricos) ---
    call make_model([2, 1, 2], [1, 1, 1], 0.3_dp, 0.03_dp, .false., 909, hh, hs)
    call hh%WriteFile('tx_L_H.hsio')
    call hs%WriteFile('tx_L_S.hsio')
    call write_x('tx_x9.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 3, 'lead', 'lead', &
                 extra="/" // new_line('a') // "&RHO nwin = 2, name = 'low', 'high', lo = -2.0, -0.3, " // &
                       "hi = -0.3, 1.5, rmax = 100.0", &
                 kblock="&KLIST avg_axis = 1, navg = 2, list_axis = 2, nlist = 2, list_k = -0.2, 0.2 /")
    call run_x('tx_x9.nml', gx)
    call check(ierr == 0, '&RHO run: ' // errmsg)
    if (ierr == 0) then
        block
            integer :: u9, ios9, n9
            character(len=16) :: wn
            integer :: i9, a9, c9
            real(dp) :: q1, q2, dq, qs
            dq = 0.0_dp; qs = 0.0_dp; n9 = 0
            open(newunit=u9, file='tx_x.rho_mulliken.dat', status='old', action='read', iostat=ios9)
            call check(ios9 == 0, '&RHO: rho_mulliken.dat written')
            if (ios9 == 0) then
                do
                    read(u9, '(A)', iostat=ios9) line9
                    if (ios9 /= 0) exit
                    if (line9(1:1) == '#') cycle
                    read(line9, *) wn, i9, a9, c9, q1, q2
                    n9 = n9 + 1
                    dq = max(dq, abs(q1 - q2))
                    qs = max(qs, abs(q2))
                end do
                close(u9, status='delete')
            end if
            call check(n9 == 2 * 12, '&RHO: two windows x 12 central atoms')
            call check(qs > 1.0e-3_dp, '&RHO: non-trivial charges in the windows')
            call check_close(dq / qs, 0.0_dp, 1.0e-10_dp, &
                             '&RHO: Mulliken sum_j tr[D_ij S_ji] = integrated Mulliken DOS (with seam terms)')
        end block
    end if
    call execute_command_line('rm -f tx_x.rho_pairs.dat', wait=.true.)

    ! --- 10. &RHO por contorno (mesmo modelo do 9) ---
    ! (a) janela que contém todo o espectro: Mulliken por átomo = número de orbitais (sum_j (S^-1)_ij S_ji = 1)
    block
        real(dp) :: qa(64), qb(64), r1, r2, rc
        integer :: na, nb10
        call write_x('tx_x10.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 3, 'lead', 'lead', &
                     extra="/" // new_line('a') // "&RHO nwin = 1, name = 'all', lo = -60.0, hi = 60.0, " // &
                           "rmax = 100.0, contour = .true., npoints = 128", &
                     kblock="&KLIST avg_axis = 1, navg = 2, list_axis = 2, nlist = 2, list_k = -0.2, 0.2 /")
        call run_x('tx_x10.nml', gx)
        call check(ierr == 0, '&RHO contour run: ' // errmsg)
        call read_qds(qa, na)
        call check(na == 12 .and. all(nint(qa(1:na)) >= 1), '&RHO contour: 12 atoms with integer charges')
        call check_close(maxval(abs(qa(1:na) - nint(qa(1:na)))), 0.0_dp, 1.0e-8_dp, &
                         '&RHO contour, whole spectrum: Mulliken per atom = number of orbitals')
        ! (b) janela parcial [-0.3, 1.5]: o contorno deslocado de i eta é a integral no eixo real com o mesmo
        !     eta; eixo real pelo ponto médio com dE = 2e-3 e 1e-3 e Richardson em dE.
        call write_x('tx_x10.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 3, 'lead', 'lead', &
                     extra="/" // new_line('a') // "&RHO nwin = 1, name = 'w', lo = -0.3, hi = 1.5, " // &
                           "rmax = 100.0, contour = .true., npoints = 192, contour_eta = 0.05", &
                     kblock="&KLIST avg_axis = 1, navg = 2, list_axis = 2, nlist = 2, list_k = -0.2, 0.2 /")
        call run_x('tx_x10.nml', gx)
        call read_qds(qa, na)
        rc = sum(qa(1:na))
        call real_axis(900, r1)
        call real_axis(1800, r2)
        write(*, '(A,ES24.16,A,2ES12.4,A,ES12.4)') '  contour total ', rc, '; real axis - contour (dE 2e-3, 1e-3) ', &
            r1 - rc, r2 - rc, '; Richardson ', (4.0_dp * r2 - r1) / 3.0_dp - rc
        call check_close(abs((4.0_dp * r2 - r1) / 3.0_dp - rc) / abs(rc), 0.0_dp, 1.0e-8_dp, &
                         '&RHO contour (offset i eta) = real-axis window charge (same eta), partial window, total')
        nb10 = 0
        qb = 0.0_dp
    end block
    call execute_command_line('rm -f tx_x.rho_pairs.dat tx_x.rho_mulliken.dat', wait=.true.)

    call execute_command_line('rm -f tx_*.hsio tx_*.nml tx_v.sitepot tx_z.sitepot', wait=.true.)
    call test_summary('test_extracted_green')

contains

    subroutine run_sg(nml, gd)
        character(len=*), intent(in) :: nml
        type(green_data_t), intent(out) :: gd
        type(green_input_t) :: inp
        call read_green_input(nml, inp, ierr, errmsg)
        if (ierr /= 0) return
        call run_green(inp, gd, ierr, errmsg, nonconverged=nbad)
    end subroutine run_sg

    !> LFSITEPOT com um bloco 2x2 [[0, b], [b, -b]] nos átomos 1 e 3 (2 orbitais) das colunas 0..2 (ao longo de a_3).
    subroutine write_blk(fname, b)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: b
        integer :: u, c, a
        open(newunit=u, file=fname, status='replace', action='write')
        write(u, '(A)') 'LFSITEPOT 1'
        write(u, '(A)') 'norb 2'
        write(u, '(A)') 'nentry 6'
        do c = 0, 2
            do a = 1, 3, 2
                write(u, '(A,I0,1X,I0,3(1X,ES24.16),A)') 'entry ', c, a, st%cart(:, a) + real(c, dp) * lat(3, :), ' block'
                write(u, '(2ES24.16)') 0.0_dp, b
                write(u, '(2ES24.16)') b, -b
            end do
        end do
        close(u)
    end subroutine write_blk

    subroutine run_x(nml, gd)
        character(len=*), intent(in) :: nml
        type(green_data_t), intent(out) :: gd
        type(extracted_input_t) :: inp
        call check(input_is_extracted(nml), 'input recognised as extracted: ' // nml)
        call read_extracted_input(nml, inp, ierr, errmsg)
        if (ierr /= 0) return
        call run_extracted(inp, gd, ierr, errmsg, nonconverged=nbad)
    end subroutine run_x

    subroutine write_kpath(u)
        integer, intent(in) :: u
        write(u, '(A)') '&KPATH'
        write(u, '(A)') '  nvertex = 3'
        write(u, '(A)') '  vertex(:,1) = 0.0, 0.0, 0.0'
        write(u, '(A)') '  vertex(:,2) = 0.5, 0.0, 0.0'
        write(u, '(A)') '  vertex(:,3) = 0.333333333333, 0.333333333333, 0.0'
        write(u, '(A)') '  npoints = 3, 2'
        write(u, '(A)') '/'
    end subroutine write_kpath

    subroutine write_surface(fname, ncentral)
        character(len=*), intent(in) :: fname
        integer, intent(in) :: ncentral
        integer :: u
        open(newunit=u, file=fname, status='replace', action='write')
        write(u, '(A)') "&GREEN mode = 'surface', eta_ev = 0.05, emin = -2.0, emax = 2.0, ne = 9, " // &
                        "output_prefix = 'tx_sg' /"
        write(u, '(A,I0,A)') "&SUBSTRATE h_file = 'tx_L_H.hsio', s_file = 'tx_L_S.hsio', " // &
            "structure_file = 'tx_L.nml', stacking_axis = 3, side = 'top', central_layers = ", ncentral, ' /'
        call write_kpath(u)
        close(u)
    end subroutine write_surface

    !> Lê q_DS (coluna 5) de tx_x.rho_mulliken.dat.
    subroutine read_qds(q, n)
        real(dp), intent(out) :: q(:)
        integer, intent(out) :: n
        integer :: u, ios, i, a, c
        character(len=16) :: wn
        character(len=512) :: ln
        real(dp) :: q2
        n = 0
        q = 0.0_dp
        open(newunit=u, file='tx_x.rho_mulliken.dat', status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read(u, '(A)', iostat=ios) ln
            if (ios /= 0) exit
            if (ln(1:1) == '#') cycle
            n = n + 1
            read(ln, *) wn, i, a, c, q(n), q2
        end do
        close(u)
    end subroutine read_qds

    !> Total q_DOS no eixo real (eta = 0.05), regra do ponto médio em [-0.3, 1.5] com nst pontos.
    subroutine real_axis(nst, tot)
        integer, intent(in) :: nst
        real(dp), intent(out) :: tot
        integer :: u, ios, i, a, c
        character(len=16) :: wn
        character(len=512) :: ln
        real(dp) :: q1, q2
        character(len=64) :: cmd
        call write_x('tx_x10r.nml', 'tx_L_H.hsio', 'tx_L_S.hsio', 'tx_L.nml', 0, 3, 'lead', 'lead', &
                     extra="/" // new_line('a') // "&RHO nwin = 1, name = 'w', lo = -0.3, hi = 1.5, rmax = 100.0", &
                     kblock="&KLIST avg_axis = 1, navg = 2, list_axis = 2, nlist = 2, list_k = -0.2, 0.2 /")
        write(cmd, '(A,F10.7,A,F10.7,A,I0)') 'emin = ', -0.3_dp + 0.9_dp / nst, ', emax = ', 1.5_dp - 0.9_dp / nst, &
            ', ne = ', nst
        call execute_command_line("sed -i 's/emin = *-2.00000000, emax = *2.00000000, ne = 9/" // &
            trim(cmd) // "/' tx_x10r.nml", wait=.true.)
        call run_x('tx_x10r.nml', gz)
        tot = 0.0_dp
        open(newunit=u, file='tx_x.rho_mulliken.dat', status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read(u, '(A)', iostat=ios) ln
            if (ios /= 0) exit
            if (ln(1:1) == '#') cycle
            read(ln, *) wn, i, a, c, q1, q2
            tot = tot + q2
        end do
        close(u)
    end subroutine real_axis

    subroutine write_x(fname, hf, sf, stf, clo, chi, left, right, extra, cross, eshift, kblock)
        character(len=*), intent(in) :: fname, hf, sf, stf, left, right
        integer, intent(in) :: clo, chi
        character(len=*), intent(in), optional :: extra, cross, kblock
        real(dp), intent(in), optional :: eshift
        integer :: u
        real(dp) :: e0
        e0 = 0.0_dp
        if (present(eshift)) e0 = eshift
        open(newunit=u, file=fname, status='replace', action='write')
        write(u, '(A,F12.8,A,F12.8,A)') "&GREEN mode = 'extracted', eta_ev = 0.05, emin = ", -2.0_dp + e0, &
            ", emax = ", 2.0_dp + e0, ", ne = 9, output_prefix = 'tx_x' /"
        write(u, '(A)') '&EXTRACTED'
        write(u, '(A)') "  h_file = '" // hf // "', s_file = '" // sf // "', structure_file = '" // stf // "'"
        write(u, '(A,I0,A,I0)') "  embed_axis = 3, central_columns = ", clo, ', ', chi
        write(u, '(A)') "  left = '" // left // "', right = '" // right // "'"
        if (present(cross)) write(u, '(A)') '  lead_cross = ' // cross
        if (present(extra)) write(u, '(A)') '  ' // extra
        write(u, '(A)') '/'
        write(u, '(A)') LEAD_NML
        if (present(kblock)) then
            write(u, '(A)') kblock
        else
            call write_kpath(u)
        end if
        close(u)
    end subroutine write_x

    !> Supercélula de NSUP células do lead ao longo de a_3 (periódica), com os
    !> átomos em ordem trocada (célula mais rápida), transladada de tt e com
    !> H -> H + dd S.
    subroutine make_supercell(dd, tt, h2, s2, st2)
        real(dp), intent(in) :: dd, tt(3)
        type(HSDataFormat), intent(out) :: h2, s2
        type(structure_t), intent(out) :: st2
        integer :: nat, nb, ia, ic, ja, jc, ir, jr, r1, r2, r3, n3, k, nr2
        integer, allocatable :: aoff(:), orbs(:), smap(:, :), soff(:)
        nat = hh%nAtoms
        nb = hh%NOrbs
        allocate(aoff(nat + 1))
        aoff(1) = 0
        do ia = 1, nat
            aoff(ia + 1) = aoff(ia) + hh%Orbs(ia)
        end do
        ! átomo da supercélula s = (ic-1)*... com a célula variando mais rápido
        allocate(smap(nat, 0:NSUP - 1), orbs(nat * NSUP), soff(nat * NSUP + 1))
        k = 0
        do ia = 1, nat
            do ic = 0, NSUP - 1
                k = k + 1
                smap(ia, ic) = k
                orbs(k) = hh%Orbs(ia)
            end do
        end do
        soff(1) = 0
        do k = 1, nat * NSUP
            soff(k + 1) = soff(k) + orbs(k)
        end do
        nr2 = 27
        h2%nAtoms = nat * NSUP; h2%NOrbs = nb * NSUP; h2%Neighbors = nr2
        allocate(h2%Orbs(nat * NSUP), h2%IndexList(3, nr2), h2%Degenerecency(nr2))
        allocate(h2%Data(h2%NOrbs, h2%NOrbs, nr2))
        h2%Orbs = orbs
        h2%Degenerecency = 1
        h2%Data = (0.0_dp, 0.0_dp)
        h2%Title = 'synthetic supercell H'
        h2%has_meta = .true.
        h2%fermi_ev = -3.4_dp
        h2%energy_reference = 'absolute'
        s2 = h2
        s2%Title = 'synthetic supercell S'
        jr = 0
        do r3 = -1, 1
            do r2 = -1, 1
                do r1 = -1, 1
                    jr = jr + 1
                    h2%IndexList(:, jr) = [r1, r2, r3]
                end do
            end do
        end do
        s2%IndexList = h2%IndexList
        do ir = 1, hh%Neighbors
            do ic = 0, NSUP - 1
                n3 = ic + hh%IndexList(3, ir)          ! célula de destino (estendida)
                jc = modulo(n3, NSUP)
                r3 = (n3 - jc) / NSUP
                jr = 1 + (hh%IndexList(1, ir) + 1) + 3 * (hh%IndexList(2, ir) + 1) + 9 * (r3 + 1)
                do ia = 1, nat
                    do ja = 1, nat
                        associate(i0 => soff(smap(ia, ic)), j0 => soff(smap(ja, jc)), &
                                  a0 => aoff(ia), b0 => aoff(ja), ni => hh%Orbs(ia), nj => hh%Orbs(ja))
                            h2%Data(i0 + 1:i0 + ni, j0 + 1:j0 + nj, jr) = hh%Data(a0 + 1:a0 + ni, b0 + 1:b0 + nj, ir) &
                                + dd * hs%Data(a0 + 1:a0 + ni, b0 + 1:b0 + nj, ir)
                            s2%Data(i0 + 1:i0 + ni, j0 + 1:j0 + nj, jr) = hs%Data(a0 + 1:a0 + ni, b0 + 1:b0 + nj, ir)
                        end associate
                    end do
                end do
            end do
        end do
        st2%nat = nat * NSUP
        st2%lattice = lat
        st2%lattice(3, :) = real(NSUP, dp) * lat(3, :)
        allocate(st2%cart(3, st2%nat), st2%label(st2%nat))
        do ia = 1, nat
            do ic = 0, NSUP - 1
                st2%cart(:, smap(ia, ic)) = st%cart(:, ia) + real(ic, dp) * lat(3, :) + tt
                st2%label(smap(ia, ic)) = st%label(ia)
            end do
        end do
    end subroutine make_supercell

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

    logical function same_bytes(f1, f2)
        character(len=*), intent(in) :: f1, f2
        integer :: u1, u2, n1, n2
        character(len=:), allocatable :: b1, b2
        inquire(file=f1, size=n1)
        inquire(file=f2, size=n2)
        same_bytes = .false.
        if (n1 /= n2 .or. n1 <= 0) return
        allocate(character(len=n1) :: b1, b2)
        open(newunit=u1, file=f1, access='stream', form='unformatted', action='read')
        open(newunit=u2, file=f2, access='stream', form='unformatted', action='read')
        read(u1) b1
        read(u2) b2
        close(u1)
        close(u2)
        same_bytes = (b1 == b2)
    end function same_bytes

end program test_extracted_green
