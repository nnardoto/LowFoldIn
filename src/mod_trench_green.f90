! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Trincheira (ou qualquer defeito de linha) num slab de base ORTOGONAL,
!> embutida lateralmente entre dois meios semi-infinitos do próprio slab.
!>
!> Geometria. O "lead" é o slab perfeito com uma coluna por célula: H(R) do
!> .hsio (S = I, p.ex. do skgen) com R = n_e a_e + n_k a_k (+ 0 a_vac), em que
!> a_e é o eixo de embedding (colunas), a_k o eixo periódico restante (k_z) e
!> a_vac o eixo de vácuo; só |n_e| <= 1 (camada principal = 1 coluna). Em k_z:
!>   H_{0,dc}(k_z) = sum_{R: n_e = dc} H(R) e^{+2 pi i k_z n_k}   (gauge da rede).
!> A região central (região D) tem ncol colunas c = 1..ncol, cada uma com os
!> átomos do lead menos os removidos (lista (coluna, átomo do lead), ou os
!> remove_top átomos de maior coordenada normal na coluna da trincheira t_c),
!> sem átomos fantasmas. Leads em c <= 0 e c >= ncol+1.
!>
!> Potencial on-site (site_potential): no átomo s da coluna c entra o bloco
!> da entrada com rótulo c - t_c e posição r_s (do lead), H_ss += V.
!>
!> Cálculo, por (k_z, E), z = E_ref + E + i eta, K = z - H:
!>  - Greens de superfície dos dois leads numa decimação (surface_gf_pair);
!>  - RGF coluna a coluna (central_green): G_cc, G_{c,c+1}, G_{c+1,c};
!>  - A_i = -(1/pi) Im sum_orb (G_cc)_ii por átomo (S = I).
!> Bordas projetadas do lead: ev(k_z) = máx_kx da banda nocc, ec(k_z) = mín da
!> banda nocc+1, com edge_nkx pontos k_x = -1/2 + j/edge_nkx.
!>
!> Referência de energia (energy_zero): 'vbm_true' = máximo de ev numa grade
!> densa própria (vbm_nkz pontos em [0, 1/2] x edge_nkx), 'vbm_grid' = máximo
!> de ev na grade de k_z em uso, 'fermi' = vbm_<fermi_shift_from> +
!> fermi_shift_ev, 'absolute' = 0. As energias de E, ev, ec e das janelas são
!> relativas a ela.
!>
!> Janelas (por k_z) [lo, hi) (lo aberto se win_lo_strict), com lo = ref_lo +
!> win_lo e hi = ref_hi + win_hi, ref em {'zero', 'vbm_true', 'vbm_grid',
!> 'fermi', 'emid', 'ev', 'ec', 'mid'}: 'ev'/'ec'/'mid' são as bordas e o meio
!> do gap locais em k_z; 'emid' = (máx ev + mín ec)/2 na grade em uso.
!> Observáveis (soma de Riemann em E com dE da grade; média ponderada em k_z):
!>   w_i = sum_kz w_k sum_{E em jan} A_i dE
!>   rho_i[m,m'] = Re sum_kz w_k sum_{E em jan} -(G_ii - G_ii^dag)/(2 pi i) dE
!>   rho_{i0,jn} = Re sum_kz w_k e^{-2 pi i k_z n} sum_{E em jan} -(G_ij - G_ji^dag)/(2 pi i) dE
!> com G_{c,c+1} da RGF; pares com i nas colunas |c - t_c| <= pair_columns,
!> j na mesma coluna ou na seguinte, |n| <= pair_nz e 1e-3 < |R_j - R_i| < pair_rmax.
!> Pesos de k_z: linspace (pesos iguais) ou lista explícita (normalizados).
!>
!> Banda complexa do lead (opcional): max |lambda| < 1 de
!>   (E - H00) - H01 lambda - H10 / lambda = 0   (E real),
!> o decaimento evanescente mais lento.
!>
!> Rotinas sem stop; erros por ierr/errmsg.
module trench_green
    use, intrinsic :: iso_fortran_env, only: int64
    use types, only: dp
    use HS_io, only: HSBlockSparse, read_hs_blocksparse
    use structure_io, only: structure_t, read_structure, write_structure
    use surface_green, only: surface_gf_pair
    use central_green, only: central_t, rgf_t, central_init, assemble_central, attach_lead, rgf_solve
    use zlinalg, only: herm_eigvals
    use site_potential, only: site_potential_t, read_site_potential, site_potential_find
    use green_data, only: green_data_t, meta_add, write_green_file, GREEN_KGRID
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    implicit none
    private

    public :: trench_input_t, trench_result_t, trench_read_input, trench_run, trench_write_outputs
    public :: TR_MAXW, TR_MAXK, TR_MAXREM

    integer, parameter :: TR_MAXW = 16, TR_MAXK = 1024, TR_MAXREM = 4096, PATHL = 512
    real(dp), parameter :: PI = 3.14159265358979323846_dp
    real(dp), parameter :: TWO_PI = 2.0_dp * PI

    !> Espelho de &TRENCH, &KZ, &WINDOWS.
    type :: trench_input_t
        ! &TRENCH
        character(len=PATHL) :: h_file = '', s_file = '', structure_file = ''
        character(len=PATHL) :: site_potential_file = '', output_prefix = 'trench'
        integer  :: embed_axis = 1, kz_axis = 3
        integer  :: nocc = 0                       !< bandas ocupadas do lead (total)
        integer  :: ncol = 15, trench_col = 0      !< 0: ncol/2 + 1
        integer  :: remove_top = 0, normal_axis = 2, atoms_per_layer = 3
        integer  :: nremove = 0
        integer  :: remove_atoms(2, TR_MAXREM) = 0 !< (coluna, átomo do lead)
        real(dp) :: pot_match_tol = 1.0e-3_dp
        real(dp) :: eta_ev = 0.010_dp
        real(dp) :: emin = -0.6_dp, emax = 1.2_dp
        integer  :: ne = 361
        character(len=16) :: energy_zero = 'vbm_true'
        integer  :: edge_nkx = 48, vbm_nkz = 21
        real(dp) :: decim_tol = 1.0e-10_dp
        integer  :: decim_maxit = 300
        logical  :: complex_band = .false.
        logical  :: rho = .true.
        integer  :: pair_columns = 3, pair_nz = 1
        real(dp) :: pair_rmax = 3.9_dp
        ! &KZ
        character(len=16) :: kz_mode = 'linspace'
        integer  :: nkz = 21
        real(dp) :: kz_min = 0.0_dp, kz_max = 0.5_dp
        real(dp) :: kz_points(TR_MAXK) = 0.0_dp, kz_weight(TR_MAXK) = 1.0_dp
        ! &WINDOWS
        integer  :: nwin = 4
        character(len=16) :: win_name(TR_MAXW) = ''
        character(len=16) :: win_lo_ref(TR_MAXW) = '', win_hi_ref(TR_MAXW) = ''
        real(dp) :: win_lo(TR_MAXW) = 0.0_dp, win_hi(TR_MAXW) = 0.0_dp
        logical  :: win_lo_strict(TR_MAXW) = .false.
        real(dp) :: fermi_shift_ev = 0.0_dp
        character(len=16) :: fermi_shift_from = 'vbm_true'
    end type trench_input_t

    type :: trench_result_t
        integer :: nk = 0, ne = 0, ncol = 0, tc = 0, nlat = 0, natc = 0, nwin = 0, npair = 0, norb_atom = 0
        real(dp), allocatable :: kz(:), wk(:), energy(:)
        real(dp), allocatable :: ev(:), ec(:)              !< relativas a eref
        real(dp) :: eref = 0.0_dp, vbm_true = 0.0_dp, vbm_grid = 0.0_dp, emid = 0.0_dp, de = 0.0_dp
        ! átomos da região central (c externo, átomo do lead interno)
        integer, allocatable :: acol(:), asite(:), alayer(:)
        real(dp), allocatable :: apos(:, :)                 !< (3, natc) Å, com o deslocamento de coluna
        real(dp), allocatable :: avtr(:)                    !< traço do potencial aplicado (eV)
        real(dp), allocatable :: vapplied(:, :, :)          !< (no, no, natc)
        real(dp), allocatable :: a(:, :, :)                 !< (ne, nk, natc) eV^-1
        real(dp), allocatable :: lam(:, :)                  !< (ne, nk)
        real(dp), allocatable :: w(:, :)                    !< (natc, nwin)
        real(dp), allocatable :: ron(:, :, :, :)            !< (no, no, natc, nwin)
        integer, allocatable :: pi(:), pj(:), pnz(:)
        real(dp), allocatable :: pvec(:, :)
        real(dp), allocatable :: rpar(:, :, :, :)           !< (no, no, npair, nwin)
        integer, allocatable :: niter(:, :)
        real(dp), allocatable :: resid(:, :), minneg(:, :)
        character(len=16), allocatable :: wname(:)
        logical, allocatable :: inwin(:, :, :)              !< (ne, nk, nwin)
        real(dp) :: lattice(3, 3) = 0.0_dp
        character(len=256), allocatable :: info(:)
    end type trench_result_t

    !> Lead em forma bloco-esparsa com os índices de célula já separados.
    type :: lead_t
        integer :: nat = 0, norb = 0, nb = 0
        integer, allocatable :: orbs(:), aoff(:)
        integer, allocatable :: bi(:), bj(:), bne(:), bnk(:)
        integer(int64), allocatable :: boff(:)
        complex(dp), allocatable :: dat(:)
    end type lead_t

    interface
        subroutine zggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
            import :: dp
            character(len=1), intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldb, ldvl, ldvr, lwork
            complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
            complex(dp), intent(out) :: alpha(*), beta(*), vl(ldvl, *), vr(ldvr, *), work(*)
            real(dp), intent(out) :: rwork(*)
            integer, intent(out) :: info
        end subroutine zggev
    end interface

contains

    ! ==================================================================
    ! Entrada
    ! ==================================================================

    subroutine trench_read_input(filename, inp, ierr, errmsg)
        character(len=*),     intent(in)  :: filename
        type(trench_input_t), intent(out) :: inp
        integer,              intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg

        character(len=PATHL) :: h_file, s_file, structure_file, site_potential_file, output_prefix
        integer  :: embed_axis, kz_axis, nocc, ncol, trench_col, remove_top, normal_axis, atoms_per_layer, nremove
        integer  :: remove_atoms(2, TR_MAXREM)
        real(dp) :: pot_match_tol, eta_ev, emin, emax, decim_tol, pair_rmax
        integer  :: ne, edge_nkx, vbm_nkz, decim_maxit, pair_columns, pair_nz
        character(len=16) :: energy_zero, kz_mode, fermi_shift_from
        logical  :: complex_band, rho
        integer  :: nkz, nwin
        real(dp) :: kz_min, kz_max, kz_points(TR_MAXK), kz_weight(TR_MAXK), fermi_shift_ev
        character(len=16) :: win_name(TR_MAXW), win_lo_ref(TR_MAXW), win_hi_ref(TR_MAXW)
        real(dp) :: win_lo(TR_MAXW), win_hi(TR_MAXW)
        logical  :: win_lo_strict(TR_MAXW)
        namelist /TRENCH/ h_file, s_file, structure_file, site_potential_file, output_prefix, embed_axis, kz_axis, &
            nocc, ncol, trench_col, remove_top, normal_axis, atoms_per_layer, nremove, remove_atoms, pot_match_tol, &
            eta_ev, emin, emax, ne, energy_zero, edge_nkx, vbm_nkz, decim_tol, decim_maxit, complex_band, rho, &
            pair_columns, pair_nz, pair_rmax
        namelist /KZ/ kz_mode, nkz, kz_min, kz_max, kz_points, kz_weight
        namelist /WINDOWS/ nwin, win_name, win_lo_ref, win_lo, win_hi_ref, win_hi, win_lo_strict, &
            fermi_shift_ev, fermi_shift_from
        integer :: u, ios, k
        character(len=256) :: iomsg

        ierr = 0
        errmsg = ''
        call default_windows(inp)
        h_file = inp%h_file; s_file = inp%s_file; structure_file = inp%structure_file
        site_potential_file = inp%site_potential_file; output_prefix = inp%output_prefix
        embed_axis = inp%embed_axis; kz_axis = inp%kz_axis; nocc = inp%nocc; ncol = inp%ncol
        trench_col = inp%trench_col; remove_top = inp%remove_top; normal_axis = inp%normal_axis
        atoms_per_layer = inp%atoms_per_layer; nremove = inp%nremove; remove_atoms = inp%remove_atoms
        pot_match_tol = inp%pot_match_tol; eta_ev = inp%eta_ev; emin = inp%emin; emax = inp%emax; ne = inp%ne
        energy_zero = inp%energy_zero; edge_nkx = inp%edge_nkx; vbm_nkz = inp%vbm_nkz
        decim_tol = inp%decim_tol; decim_maxit = inp%decim_maxit; complex_band = inp%complex_band; rho = inp%rho
        pair_columns = inp%pair_columns; pair_nz = inp%pair_nz; pair_rmax = inp%pair_rmax
        kz_mode = inp%kz_mode; nkz = inp%nkz; kz_min = inp%kz_min; kz_max = inp%kz_max
        kz_points = inp%kz_points; kz_weight = inp%kz_weight
        nwin = inp%nwin; win_name = inp%win_name; win_lo_ref = inp%win_lo_ref; win_lo = inp%win_lo
        win_hi_ref = inp%win_hi_ref; win_hi = inp%win_hi; win_lo_strict = inp%win_lo_strict
        fermi_shift_ev = inp%fermi_shift_ev; fermi_shift_from = inp%fermi_shift_from

        open(newunit=u, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = 'cannot open ' // trim(filename); return
        end if
        read(u, nml=TRENCH, iostat=ios, iomsg=iomsg)
        if (ios /= 0) then
            ierr = 1; errmsg = 'error reading &TRENCH: ' // trim(iomsg); close(u); return
        end if
        rewind(u)
        read(u, nml=KZ, iostat=ios, iomsg=iomsg)
        if (ios > 0) then
            ierr = 1; errmsg = 'error reading &KZ: ' // trim(iomsg); close(u); return
        end if
        rewind(u)
        read(u, nml=WINDOWS, iostat=ios, iomsg=iomsg)
        if (ios > 0) then
            ierr = 1; errmsg = 'error reading &WINDOWS: ' // trim(iomsg); close(u); return
        end if
        close(u)

        inp%h_file = h_file; inp%s_file = s_file; inp%structure_file = structure_file
        inp%site_potential_file = site_potential_file; inp%output_prefix = output_prefix
        inp%embed_axis = embed_axis; inp%kz_axis = kz_axis; inp%nocc = nocc; inp%ncol = ncol
        inp%trench_col = trench_col; inp%remove_top = remove_top; inp%normal_axis = normal_axis
        inp%atoms_per_layer = atoms_per_layer; inp%nremove = nremove; inp%remove_atoms = remove_atoms
        inp%pot_match_tol = pot_match_tol; inp%eta_ev = eta_ev; inp%emin = emin; inp%emax = emax; inp%ne = ne
        inp%energy_zero = energy_zero; inp%edge_nkx = edge_nkx; inp%vbm_nkz = vbm_nkz
        inp%decim_tol = decim_tol; inp%decim_maxit = decim_maxit; inp%complex_band = complex_band; inp%rho = rho
        inp%pair_columns = pair_columns; inp%pair_nz = pair_nz; inp%pair_rmax = pair_rmax
        inp%kz_mode = kz_mode; inp%nkz = nkz; inp%kz_min = kz_min; inp%kz_max = kz_max
        inp%kz_points = kz_points; inp%kz_weight = kz_weight
        inp%nwin = nwin; inp%win_name = win_name; inp%win_lo_ref = win_lo_ref; inp%win_lo = win_lo
        inp%win_hi_ref = win_hi_ref; inp%win_hi = win_hi; inp%win_lo_strict = win_lo_strict
        inp%fermi_shift_ev = fermi_shift_ev; inp%fermi_shift_from = fermi_shift_from

        if (inp%trench_col == 0) inp%trench_col = inp%ncol / 2 + 1
        if (inp%ncol < 1 .or. inp%trench_col < 1 .or. inp%trench_col > inp%ncol) then
            ierr = 1; errmsg = '&TRENCH: invalid ncol/trench_col'; return
        end if
        if (inp%nocc < 1) then
            ierr = 1; errmsg = '&TRENCH: nocc (number of occupied lead bands) is required'; return
        end if
        if (inp%ne < 1 .or. inp%eta_ev <= 0.0_dp .or. inp%edge_nkx < 1 .or. inp%vbm_nkz < 1) then
            ierr = 1; errmsg = '&TRENCH: invalid ne/eta_ev/edge_nkx/vbm_nkz'; return
        end if
        if (inp%nremove < 0 .or. inp%nremove > TR_MAXREM) then
            ierr = 1; errmsg = '&TRENCH: invalid nremove'; return
        end if
        if (inp%nkz < 1 .or. inp%nkz > TR_MAXK) then
            ierr = 1; errmsg = '&KZ: invalid nkz'; return
        end if
        if (inp%kz_mode /= 'linspace' .and. inp%kz_mode /= 'list') then
            ierr = 1; errmsg = "&KZ: kz_mode must be 'linspace' or 'list'"; return
        end if
        if (inp%nwin < 0 .or. inp%nwin > TR_MAXW) then
            ierr = 1; errmsg = '&WINDOWS: invalid nwin'; return
        end if
        do k = 1, inp%nwin
            if (.not. valid_ref(inp%win_lo_ref(k)) .or. .not. valid_ref(inp%win_hi_ref(k))) then
                ierr = 1; errmsg = '&WINDOWS: invalid reference in window ' // trim(inp%win_name(k)); return
            end if
        end do
        if (.not. any(inp%energy_zero == [character(len=16) :: 'vbm_true', 'vbm_grid', 'fermi', 'absolute'])) then
            ierr = 1; errmsg = "&TRENCH: energy_zero must be 'vbm_true', 'vbm_grid', 'fermi' or 'absolute'"; return
        end if
        if (.not. any(inp%fermi_shift_from == [character(len=16) :: 'vbm_true', 'vbm_grid'])) then
            ierr = 1; errmsg = "&WINDOWS: fermi_shift_from must be 'vbm_true' or 'vbm_grid'"; return
        end if
    end subroutine trench_read_input

    !> As quatro janelas do atom_windows.py / trench_cross.py.
    subroutine default_windows(inp)
        type(trench_input_t), intent(inout) :: inp
        inp%nwin = 4
        inp%win_name(1:4) = [character(len=16) :: 'wv_gap', 'wc_gap', 'wv_05', 'wc_05']
        inp%win_lo_ref(1:4) = [character(len=16) :: 'ev', 'mid', 'emid', 'emid']
        inp%win_lo(1:4) = [0.015_dp, 0.0_dp, -0.5_dp, 0.0_dp]
        inp%win_hi_ref(1:4) = [character(len=16) :: 'mid', 'ec', 'emid', 'emid']
        inp%win_hi(1:4) = [0.0_dp, -0.015_dp, 0.0_dp, 0.5_dp]
        inp%win_lo_strict(1:4) = [.true., .false., .false., .false.]
    end subroutine default_windows

    pure logical function valid_ref(r)
        character(len=*), intent(in) :: r
        valid_ref = any(r == [character(len=16) :: 'zero', 'vbm_true', 'vbm_grid', 'fermi', 'emid', 'ev', 'ec', 'mid'])
    end function valid_ref

    ! ==================================================================
    ! Cálculo
    ! ==================================================================

    subroutine trench_run(inp, res, ierr, errmsg)
        type(trench_input_t),  intent(in)  :: inp
        type(trench_result_t), intent(out) :: res
        integer,               intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg

        type(structure_t) :: st
        type(HSBlockSparse) :: hb, sb
        type(lead_t) :: ld
        type(site_potential_t) :: sp
        character(len=256) :: emsg
        integer :: nl, nk, ne, ncol, tc, ik, ie, c, s, g, k, iw, nthr, no, vac_axis, nrem
        logical, allocatable :: present_(:, :)
        integer, allocatable :: corb(:, :), cn(:), gfirst(:), gloc(:), order(:), layer(:)
        complex(dp), allocatable :: h00k(:, :, :), h01k(:, :, :), h10k(:, :, :)
        real(dp), allocatable :: vat(:, :, :), eband(:)
        real(dp) :: trans(3, 1), dkz, de, lo, hi
        integer :: nerr_shared

        ierr = 0
        errmsg = ''

        ! ---------------- lead ----------------
        call read_structure(trim(inp%structure_file), st, ierr, emsg)
        if (ierr /= 0) then
            errmsg = trim(emsg); return
        end if
        call read_hs_blocksparse(trim(inp%h_file), hb, ierr, errmsg)
        if (ierr /= 0) return
        call read_hs_blocksparse(trim(inp%s_file), sb, ierr, errmsg)
        if (ierr /= 0) return
        if (hb%nAtoms /= st%nat .or. sb%nAtoms /= st%nat) then
            ierr = 1; errmsg = 'trench: structure and .hsio have different numbers of atoms'; return
        end if
        vac_axis = 6 - inp%embed_axis - inp%kz_axis
        if (inp%embed_axis == inp%kz_axis .or. vac_axis < 1 .or. vac_axis > 3) then
            ierr = 1; errmsg = 'trench: embed_axis and kz_axis must be two different lattice axes'; return
        end if
        call check_identity_s(sb, ierr, errmsg)
        if (ierr /= 0) return
        call lead_from_blocksparse(hb, inp%embed_axis, inp%kz_axis, vac_axis, ld, ierr, errmsg)
        if (ierr /= 0) return
        nl = ld%norb
        if (inp%nocc >= nl) then
            ierr = 1; errmsg = 'trench: nocc must be smaller than the number of lead orbitals'; return
        end if
        no = ld%orbs(1)
        if (inp%rho .and. any(ld%orbs /= no)) then
            ierr = 1; errmsg = 'trench: rho output needs the same number of orbitals on every atom'; return
        end if
        res%norb_atom = no
        res%lattice = st%lattice

        ! ---------------- grades ----------------
        nk = inp%nkz
        ne = inp%ne
        ncol = inp%ncol
        tc = inp%trench_col
        res%nk = nk; res%ne = ne; res%ncol = ncol; res%tc = tc; res%nlat = st%nat
        allocate(res%kz(nk), res%wk(nk), res%energy(ne))
        if (inp%kz_mode == 'linspace') then
            res%kz = linspace(inp%kz_min, inp%kz_max, nk)
            res%wk = 1.0_dp / real(nk, dp)
        else
            res%kz = inp%kz_points(1:nk)
            if (any(inp%kz_weight(1:nk) < 0.0_dp) .or. sum(inp%kz_weight(1:nk)) <= 0.0_dp) then
                ierr = 1; errmsg = '&KZ: weights must be non-negative with positive sum'; return
            end if
            res%wk = inp%kz_weight(1:nk) / sum(inp%kz_weight(1:nk))
        end if
        res%energy = linspace(inp%emin, inp%emax, ne)
        de = 0.0_dp
        if (ne > 1) de = res%energy(2) - res%energy(1)
        res%de = de

        ! ---------------- bordas e referência ----------------
        allocate(h00k(nl, nl, nk), h01k(nl, nl, nk), h10k(nl, nl, nk), res%ev(nk), res%ec(nk), eband(nl))
        do ik = 1, nk
            call lead_kblocks(ld, res%kz(ik), h00k(:, :, ik), h01k(:, :, ik), h10k(:, :, ik))
            call edges_at(h00k(:, :, ik), h01k(:, :, ik), h10k(:, :, ik), inp%edge_nkx, inp%nocc, &
                          res%ev(ik), res%ec(ik), ierr, errmsg)
            if (ierr /= 0) return
        end do
        res%vbm_grid = maxval(res%ev)
        block
            complex(dp), allocatable :: a0(:, :), a1(:, :), a2(:, :)
            real(dp) :: v, cdum, kzd
            allocate(a0(nl, nl), a1(nl, nl), a2(nl, nl))
            res%vbm_true = -huge(1.0_dp)
            do k = 1, inp%vbm_nkz
                kzd = 0.0_dp
                if (inp%vbm_nkz > 1) then
                    dkz = 0.5_dp / real(inp%vbm_nkz - 1, dp)
                    kzd = real(k - 1, dp) * dkz
                    if (k == inp%vbm_nkz) kzd = 0.5_dp
                end if
                call lead_kblocks(ld, kzd, a0, a1, a2)
                call edges_at(a0, a1, a2, inp%edge_nkx, inp%nocc, v, cdum, ierr, errmsg)
                if (ierr /= 0) return
                res%vbm_true = max(res%vbm_true, v)
            end do
        end block
        select case (trim(inp%energy_zero))
        case ('vbm_true')
            res%eref = res%vbm_true
        case ('vbm_grid')
            res%eref = res%vbm_grid
        case ('fermi')
            res%eref = merge(res%vbm_true, res%vbm_grid, inp%fermi_shift_from == 'vbm_true') + inp%fermi_shift_ev
        case default
            res%eref = 0.0_dp
        end select
        res%ev = res%ev - res%eref
        res%ec = res%ec - res%eref
        res%emid = 0.5_dp * (maxval(res%ev) + minval(res%ec))

        ! janelas
        res%nwin = inp%nwin
        allocate(res%wname(inp%nwin), res%inwin(ne, nk, inp%nwin))
        res%inwin = .false.
        do iw = 1, inp%nwin
            res%wname(iw) = inp%win_name(iw)
            do ik = 1, nk
                lo = ref_value(inp%win_lo_ref(iw), ik) + inp%win_lo(iw)
                hi = ref_value(inp%win_hi_ref(iw), ik) + inp%win_hi(iw)
                do ie = 1, ne
                    if (inp%win_lo_strict(iw)) then
                        res%inwin(ie, ik, iw) = res%energy(ie) > lo .and. res%energy(ie) < hi
                    else
                        res%inwin(ie, ik, iw) = res%energy(ie) >= lo .and. res%energy(ie) < hi
                    end if
                end do
            end do
        end do

        ! ---------------- região central ----------------
        allocate(present_(st%nat, ncol))
        present_ = .true.
        allocate(order(st%nat), layer(st%nat))
        order = argsort_desc(st%cart(inp%normal_axis, :))
        do k = 1, st%nat
            layer(order(k)) = (k - 1) / max(inp%atoms_per_layer, 1)
        end do
        if (inp%remove_top > 0) then
            if (inp%remove_top > st%nat) then
                ierr = 1; errmsg = '&TRENCH: remove_top larger than the number of lead atoms'; return
            end if
            present_(order(1:inp%remove_top), tc) = .false.
        end if
        do k = 1, inp%nremove
            c = inp%remove_atoms(1, k); s = inp%remove_atoms(2, k)
            if (c < 1 .or. c > ncol .or. s < 1 .or. s > st%nat) then
                ierr = 1; errmsg = '&TRENCH: remove_atoms entry out of range'; return
            end if
            present_(s, c) = .false.
        end do
        nrem = count(.not. present_)
        res%natc = count(present_)
        if (res%natc < 1 .or. any(count(present_, dim=1) < 1)) then
            ierr = 1; errmsg = 'trench: every column must keep at least one atom'; return
        end if
        allocate(res%acol(res%natc), res%asite(res%natc), res%alayer(res%natc), res%apos(3, res%natc))
        allocate(cn(ncol), corb(nl, ncol), gfirst(ncol + 1), gloc(res%natc))
        g = 0
        do c = 1, ncol
            gfirst(c) = g + 1
            cn(c) = 0
            do s = 1, st%nat
                if (.not. present_(s, c)) cycle
                g = g + 1
                res%acol(g) = c; res%asite(g) = s; res%alayer(g) = layer(s)
                res%apos(:, g) = st%cart(:, s) + real(c - tc, dp) * st%lattice(inp%embed_axis, :)
                gloc(g) = cn(c)                          ! deslocamento (0-based) no bloco da coluna
                corb(cn(c) + 1:cn(c) + ld%orbs(s), c) = [(ld%aoff(s) + k, k = 1, ld%orbs(s))]
                cn(c) = cn(c) + ld%orbs(s)
            end do
        end do
        gfirst(ncol + 1) = g + 1

        ! potencial on-site
        allocate(vat(no, no, res%natc), res%avtr(res%natc))
        vat = 0.0_dp
        if (len_trim(inp%site_potential_file) > 0) then
            call read_site_potential(trim(inp%site_potential_file), sp, ierr, errmsg)
            if (ierr /= 0) return
            if (sp%norb /= no) then
                ierr = 1; errmsg = 'trench: site potential norb differs from the lead basis'; return
            end if
            trans(:, 1) = st%lattice(inp%kz_axis, :)
            do g = 1, res%natc
                k = site_potential_find(sp, res%acol(g) - tc, st%cart(:, res%asite(g)), trans, inp%pot_match_tol)
                if (k > 0) vat(:, :, g) = 0.5_dp * (sp%v(:, :, k) + transpose(sp%v(:, :, k)))
            end do
        end if
        do g = 1, res%natc
            res%avtr(g) = sum([(vat(k, k, g), k = 1, no)])
        end do
        res%vapplied = vat

        ! pares
        call build_pairs()

        ! ---------------- laço (k_z, E) ----------------
        allocate(res%a(ne, nk, res%natc), res%niter(ne, nk), res%resid(ne, nk), res%minneg(ne, nk))
        res%a = 0.0_dp; res%niter = 0; res%resid = 0.0_dp; res%minneg = 0.0_dp
        if (inp%complex_band) then
            allocate(res%lam(ne, nk))
            res%lam = 0.0_dp
        end if
        nthr = 1
        !$ nthr = omp_get_max_threads()
        block
            complex(dp), allocatable :: ronb(:, :, :, :, :), rparb(:, :, :, :, :)
            integer :: nwr, npr
            nwr = 0; npr = 0
            if (inp%rho) then
                nwr = inp%nwin; npr = res%npair
            end if
            allocate(ronb(no, no, res%natc, nwr, nthr), rparb(no, no, npr, nwr, nthr))
            ronb = (0.0_dp, 0.0_dp); rparb = (0.0_dp, 0.0_dp)
            nerr_shared = 0
            !$omp parallel default(shared) private(ik, ie)
            block
                integer :: it, tid, ierr_l
                character(len=:), allocatable :: emsg_l
                tid = 1
                !$ tid = omp_get_thread_num() + 1
                !$omp do schedule(dynamic, 1)
                do it = 1, nk * ne
                    ik = (it - 1) / ne + 1
                    ie = mod(it - 1, ne) + 1
                    call point(ik, ie, ronb(:, :, :, :, tid), rparb(:, :, :, :, tid), ierr_l, emsg_l)
                    if (ierr_l /= 0) then
                        !$omp critical (trench_err)
                        if (nerr_shared == 0) errmsg = emsg_l
                        nerr_shared = nerr_shared + 1
                        !$omp end critical (trench_err)
                    end if
                end do
                !$omp end do
            end block
            !$omp end parallel
            if (nerr_shared > 0) then
                ierr = 2; return
            end if
            if (inp%rho) then
                allocate(res%ron(no, no, res%natc, nwr), res%rpar(no, no, npr, nwr))
                res%ron = real(sum(ronb, dim=5), dp)
                res%rpar = real(sum(rparb, dim=5), dp)
            end if
        end block

        ! pesos por átomo (determinísticos, a partir de A)
        allocate(res%w(res%natc, inp%nwin))
        res%w = 0.0_dp
        do iw = 1, inp%nwin
            do g = 1, res%natc
                do ik = 1, nk
                    res%w(g, iw) = res%w(g, iw) + res%wk(ik) * sum(res%a(:, ik, g), mask=res%inwin(:, ik, iw)) * de
                end do
            end do
        end do

        allocate(res%info(8))
        write(res%info(1), '(a,i0,a,i0,a,i0,a,i0)') 'lead atoms = ', st%nat, ', lead orbitals = ', nl, &
            ', central atoms = ', res%natc, ', removed = ', nrem
        write(res%info(2), '(a,f12.6,a,f12.6,a,f12.6)') 'vbm_true = ', res%vbm_true, ' eV, vbm_grid = ', res%vbm_grid, &
            ' eV, eref = ', res%eref
        write(res%info(3), '(a,f10.6,a,f10.6,a)') 'projected gap of the lead: ', minval(res%ec) - maxval(res%ev), &
            ' eV; emid = ', res%emid, ' eV (relative to eref)'
        write(res%info(4), '(a,i0,a,i0)') 'decimation iterations (max) = ', maxval(res%niter), ', npairs = ', res%npair
        write(res%info(5), '(a,es10.3,a,es10.3)') 'max surface-gf residual = ', maxval(res%resid), &
            ', min(-Im G_ii) = ', minval(res%minneg)
        write(res%info(6), '(a,i0,a,f8.4)') 'atoms with a site potential: ', count(maxval(maxval(abs(vat), dim=1), dim=1) > 0.0_dp), &
            '; max |V| = ', maxval(abs(vat))
        write(res%info(7), '(a,i0,a,i0)') 'kz points = ', nk, ', energies = ', ne
        write(res%info(8), '(a,a)') 'energy zero = ', trim(inp%energy_zero)

    contains

        real(dp) function ref_value(r, ik) result(v)
            character(len=*), intent(in) :: r
            integer, intent(in) :: ik
            select case (trim(r))
            case ('zero')
                v = 0.0_dp
            case ('vbm_true')
                v = res%vbm_true - res%eref
            case ('vbm_grid')
                v = res%vbm_grid - res%eref
            case ('fermi')
                v = merge(res%vbm_true, res%vbm_grid, inp%fermi_shift_from == 'vbm_true') + inp%fermi_shift_ev - res%eref
            case ('emid')
                v = res%emid
            case ('ev')
                v = res%ev(ik)
            case ('ec')
                v = res%ec(ik)
            case default        ! 'mid'
                v = 0.5_dp * (res%ev(ik) + res%ec(ik))
            end select
        end function ref_value

        subroutine build_pairs()
            integer :: c1, dc, c2, n3, gi, gj, np
            real(dp) :: d(3), r
            integer, allocatable :: ti(:), tj(:), tn(:)
            real(dp), allocatable :: tv(:, :)
            np = 0
            allocate(ti(0), tj(0), tn(0), tv(3, 0))
            if (inp%rho) then
                do c1 = max(1, tc - inp%pair_columns), min(ncol, tc + inp%pair_columns)
                    do dc = 0, 1
                        c2 = c1 + dc
                        if (c2 > ncol) cycle
                        do n3 = -inp%pair_nz, inp%pair_nz
                            do gi = gfirst(c1), gfirst(c1 + 1) - 1
                                do gj = gfirst(c2), gfirst(c2 + 1) - 1
                                    d = st%cart(:, res%asite(gj)) + real(dc, dp) * st%lattice(inp%embed_axis, :) &
                                        + real(n3, dp) * st%lattice(inp%kz_axis, :) - st%cart(:, res%asite(gi))
                                    r = norm2(d)
                                    if (r > 1.0e-3_dp .and. r < inp%pair_rmax) then
                                        np = np + 1
                                        ti = [ti, gi]; tj = [tj, gj]; tn = [tn, n3]
                                        tv = reshape([tv, d], [3, np])
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
            end if
            res%npair = np
            res%pi = ti; res%pj = tj; res%pnz = tn; res%pvec = tv
        end subroutine build_pairs

        !> Um ponto (k_z, E): decimação, RGF, A, rho.
        subroutine point(ik, ie, ron_t, rpar_t, ierr_p, emsg_p)
            integer, intent(in) :: ik, ie
            complex(dp), intent(inout) :: ron_t(:, :, :, :), rpar_t(:, :, :, :)
            integer, intent(out) :: ierr_p
            character(len=:), allocatable, intent(out) :: emsg_p
            complex(dp), allocatable :: k00(:, :), k01(:, :), k10(:, :), gtop(:, :), gbot(:, :)
            complex(dp), allocatable :: hd(:, :), sd(:, :), hu(:, :), su(:, :), hl(:, :), sl(:, :), zero_xc(:, :)
            complex(dp), allocatable :: gij(:, :), gji(:, :)
            type(central_t) :: sys
            type(rgf_t) :: rg
            complex(dp) :: z, fac, ph
            real(dp) :: rt, rb, mn, acc
            integer :: niter, c1, i, gg, m, n1, iw2, p, gi, gj, oi, oj, ci, cj
            logical :: anywin

            ierr_p = 0
            emsg_p = ''
            z = cmplx(res%energy(ie) + res%eref, inp%eta_ev, dp)
            allocate(k00(nl, nl), k01(nl, nl), k10(nl, nl), gtop(nl, nl), gbot(nl, nl))
            k00 = -h00k(:, :, ik)
            do i = 1, nl
                k00(i, i) = k00(i, i) + z
            end do
            k01 = -h01k(:, :, ik)
            k10 = -h10k(:, :, ik)
            call surface_gf_pair(k00, k01, k10, inp%decim_tol, inp%decim_maxit, gtop, gbot, niter, rt, rb, &
                                 ierr_p, emsg_p)
            if (ierr_p /= 0) return
            res%niter(ie, ik) = niter
            res%resid(ie, ik) = max(rt, rb)

            call central_init(sys, cn)
            do c1 = 1, ncol
                n1 = cn(c1)
                hd = h00k(corb(1:n1, c1), corb(1:n1, c1), ik)
                do gg = gfirst(c1), gfirst(c1 + 1) - 1
                    m = ld%orbs(res%asite(gg))
                    oi = gloc(gg)
                    hd(oi + 1:oi + m, oi + 1:oi + m) = hd(oi + 1:oi + m, oi + 1:oi + m) + &
                        cmplx(vat(1:m, 1:m, gg), 0.0_dp, dp)
                end do
                if (allocated(sd)) deallocate(sd)
                allocate(sd(n1, n1))
                sd = (0.0_dp, 0.0_dp)
                do i = 1, n1
                    sd(i, i) = (1.0_dp, 0.0_dp)
                end do
                if (c1 < ncol) then
                    hu = h01k(corb(1:n1, c1), corb(1:cn(c1 + 1), c1 + 1), ik)
                    hl = h10k(corb(1:cn(c1 + 1), c1 + 1), corb(1:n1, c1), ik)
                    if (allocated(su)) deallocate(su)
                    if (allocated(sl)) deallocate(sl)
                    allocate(su(n1, cn(c1 + 1)), sl(cn(c1 + 1), n1))
                    su = (0.0_dp, 0.0_dp); sl = (0.0_dp, 0.0_dp)
                    call assemble_central(sys, c1, z, hd, sd, 0.0_dp, hu, su, 0.0_dp, hl, sl)
                else
                    call assemble_central(sys, c1, z, hd, sd, 0.0_dp)
                end if
            end do
            ! lead esquerdo (c = 0) -> coluna 1: H_{A,1} = H01(:, coluna 1)
            hu = h01k(:, corb(1:cn(1), 1), ik)
            allocate(zero_xc(nl, cn(1)))
            zero_xc = (0.0_dp, 0.0_dp)
            call attach_lead(sys, 'bot', z, hu, zero_xc, 0.0_dp, 0.0_dp, gtop)
            ! coluna ncol -> lead direito: H_{N,B} = H01(coluna N, :)
            hu = h01k(corb(1:cn(ncol), ncol), :, ik)
            deallocate(zero_xc)
            allocate(zero_xc(cn(ncol), nl))
            zero_xc = (0.0_dp, 0.0_dp)
            call attach_lead(sys, 'top', z, hu, zero_xc, 0.0_dp, 0.0_dp, gbot)
            call rgf_solve(sys, rg, ierr_p, emsg_p)
            if (ierr_p /= 0) return

            mn = huge(1.0_dp)
            do gg = 1, res%natc
                c1 = res%acol(gg); m = ld%orbs(res%asite(gg)); oi = gloc(gg)
                acc = 0.0_dp
                do i = oi + 1, oi + m
                    acc = acc + aimag(rg%gd(c1)%a(i, i))
                    mn = min(mn, -aimag(rg%gd(c1)%a(i, i)))
                end do
                res%a(ie, ik, gg) = -acc / PI
            end do
            res%minneg(ie, ik) = mn

            if (inp%rho .and. inp%nwin > 0) then
                anywin = any(res%inwin(ie, ik, :))
                if (anywin) then
                    ! -(G - G^dag)/(2 pi i) dE w_k = (i/2pi)(G - G^dag) dE w_k
                    fac = cmplx(0.0_dp, res%de * res%wk(ik) / TWO_PI, dp)
                    do gg = 1, res%natc
                        c1 = res%acol(gg); oi = gloc(gg)
                        gij = rg%gd(c1)%a(oi + 1:oi + no, oi + 1:oi + no)
                        gij = fac * (gij - conjg(transpose(gij)))
                        do iw2 = 1, inp%nwin
                            if (res%inwin(ie, ik, iw2)) ron_t(:, :, gg, iw2) = ron_t(:, :, gg, iw2) + gij
                        end do
                    end do
                    do p = 1, res%npair
                        gi = res%pi(p); gj = res%pj(p)
                        ci = res%acol(gi); cj = res%acol(gj); oi = gloc(gi); oj = gloc(gj)
                        if (ci == cj) then
                            gij = rg%gd(ci)%a(oi + 1:oi + no, oj + 1:oj + no)
                            gji = rg%gd(ci)%a(oj + 1:oj + no, oi + 1:oi + no)
                        else
                            gij = rg%gu(ci)%a(oi + 1:oi + no, oj + 1:oj + no)
                            gji = rg%gl(ci)%a(oj + 1:oj + no, oi + 1:oi + no)
                        end if
                        ph = exp(cmplx(0.0_dp, -TWO_PI * res%kz(ik) * real(res%pnz(p), dp), dp))
                        gij = (fac * ph) * (gij - conjg(transpose(gji)))
                        do iw2 = 1, inp%nwin
                            if (res%inwin(ie, ik, iw2)) rpar_t(:, :, p, iw2) = rpar_t(:, :, p, iw2) + gij
                        end do
                    end do
                end if
            end if

            if (inp%complex_band) then
                res%lam(ie, ik) = slowest_decay(h00k(:, :, ik), h01k(:, :, ik), h10k(:, :, ik), &
                                                res%energy(ie) + res%eref, ierr_p)
                if (ierr_p /= 0) emsg_p = 'trench: zggev failed in the complex band'
            end if
        end subroutine point

    end subroutine trench_run

    !> S deve ser a identidade (só blocos on-site em R = 0).
    subroutine check_identity_s(sb, ierr, errmsg)
        type(HSBlockSparse), intent(in) :: sb
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: ib, i, j, n
        integer(int64) :: k
        complex(dp) :: x
        logical :: onsite
        ierr = 0
        errmsg = ''
        do ib = 1, int(sb%nBlocks)
            i = sb%BlockIndex(1, ib); j = sb%BlockIndex(2, ib)
            onsite = (i == j) .and. all(sb%IndexList(:, sb%BlockIndex(3, ib)) == 0)
            n = sb%Orbs(i)
            do k = 1_int64, sb%BlockOffset(ib + 1) - sb%BlockOffset(ib)
                x = sb%BlockData(sb%BlockOffset(ib) + k)
                if (onsite .and. mod(k - 1, int(n, int64)) == (k - 1) / n) x = x - (1.0_dp, 0.0_dp)
                if (abs(x) > 1.0e-10_dp) then
                    ierr = 1; errmsg = 'trench: the lead basis must be orthogonal (S = I); use an SK .hsio from skgen'
                    return
                end if
            end do
        end do
    end subroutine check_identity_s

    subroutine lead_from_blocksparse(hb, ea, ka, va, ld, ierr, errmsg)
        type(HSBlockSparse), intent(in) :: hb
        integer, intent(in) :: ea, ka, va
        type(lead_t), intent(out) :: ld
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: ib, a, r(3)
        ierr = 0
        errmsg = ''
        ld%nat = hb%nAtoms
        ld%orbs = hb%Orbs
        allocate(ld%aoff(ld%nat + 1))
        ld%aoff(1) = 0
        do a = 1, ld%nat
            ld%aoff(a + 1) = ld%aoff(a) + ld%orbs(a)
        end do
        ld%norb = ld%aoff(ld%nat + 1)
        ld%nb = int(hb%nBlocks)
        allocate(ld%bi(ld%nb), ld%bj(ld%nb), ld%bne(ld%nb), ld%bnk(ld%nb))
        do ib = 1, ld%nb
            r = hb%IndexList(:, hb%BlockIndex(3, ib))
            if (r(va) /= 0) then
                if (any(abs(hb%BlockData(hb%BlockOffset(ib) + 1:hb%BlockOffset(ib + 1))) > 0.0_dp)) then
                    ierr = 1; errmsg = 'trench: the lead has couplings along the vacuum axis'; return
                end if
            end if
            if (abs(r(ea)) > 1) then
                if (any(abs(hb%BlockData(hb%BlockOffset(ib) + 1:hb%BlockOffset(ib + 1))) > 0.0_dp)) then
                    ierr = 1; errmsg = 'trench: coupling beyond the neighbouring column (principal layer > 1 cell)'
                    return
                end if
            end if
            ld%bi(ib) = hb%BlockIndex(1, ib); ld%bj(ib) = hb%BlockIndex(2, ib)
            ld%bne(ib) = r(ea); ld%bnk(ib) = r(ka)
            if (r(va) /= 0) ld%bne(ib) = 99
        end do
        ld%boff = hb%BlockOffset
        ld%dat = hb%BlockData
    end subroutine lead_from_blocksparse

    !> H00, H01 = H_{0,+1}, H10 = H_{0,-1} em k_z (gauge da rede).
    subroutine lead_kblocks(ld, kz, h00, h01, h10)
        type(lead_t), intent(in) :: ld
        real(dp), intent(in) :: kz
        complex(dp), intent(out) :: h00(:, :), h01(:, :), h10(:, :)
        integer :: ib, i, j, ni, nj
        complex(dp) :: ph
        complex(dp), allocatable :: b(:, :)
        h00 = (0.0_dp, 0.0_dp); h01 = (0.0_dp, 0.0_dp); h10 = (0.0_dp, 0.0_dp)
        do ib = 1, ld%nb
            if (abs(ld%bne(ib)) > 1) cycle
            i = ld%bi(ib); j = ld%bj(ib)
            ni = ld%orbs(i); nj = ld%orbs(j)
            ph = exp(cmplx(0.0_dp, TWO_PI * kz * real(ld%bnk(ib), dp), dp))
            b = ph * reshape(ld%dat(ld%boff(ib) + 1:ld%boff(ib + 1)), [ni, nj])
            select case (ld%bne(ib))
            case (0)
                h00(ld%aoff(i) + 1:ld%aoff(i) + ni, ld%aoff(j) + 1:ld%aoff(j) + nj) = &
                    h00(ld%aoff(i) + 1:ld%aoff(i) + ni, ld%aoff(j) + 1:ld%aoff(j) + nj) + b
            case (1)
                h01(ld%aoff(i) + 1:ld%aoff(i) + ni, ld%aoff(j) + 1:ld%aoff(j) + nj) = &
                    h01(ld%aoff(i) + 1:ld%aoff(i) + ni, ld%aoff(j) + 1:ld%aoff(j) + nj) + b
            case (-1)
                h10(ld%aoff(i) + 1:ld%aoff(i) + ni, ld%aoff(j) + 1:ld%aoff(j) + nj) = &
                    h10(ld%aoff(i) + 1:ld%aoff(i) + ni, ld%aoff(j) + 1:ld%aoff(j) + nj) + b
            end select
        end do
    end subroutine lead_kblocks

    !> ev = máx_kx E_nocc(kx), ec = mín_kx E_{nocc+1}(kx), kx = -1/2 + j/nkx.
    subroutine edges_at(h00, h01, h10, nkx, nocc, ev, ec, ierr, errmsg)
        complex(dp), intent(in) :: h00(:, :), h01(:, :), h10(:, :)
        integer, intent(in) :: nkx, nocc
        real(dp), intent(out) :: ev, ec
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        complex(dp), allocatable :: hk(:, :)
        real(dp), allocatable :: w(:)
        real(dp) :: kx
        integer :: j, n
        n = size(h00, 1)
        allocate(hk(n, n), w(n))
        ev = -huge(1.0_dp); ec = huge(1.0_dp)
        ierr = 0
        errmsg = ''
        do j = 0, nkx - 1
            kx = -0.5_dp + real(j, dp) / real(nkx, dp)
            hk = h00 + exp(cmplx(0.0_dp, TWO_PI * kx, dp)) * h01 + exp(cmplx(0.0_dp, -TWO_PI * kx, dp)) * h10
            call herm_eigvals(hk, w, ierr, errmsg)
            if (ierr /= 0) return
            ev = max(ev, w(nocc))
            ec = min(ec, w(nocc + 1))
        end do
    end subroutine edges_at

    !> max |lambda| < 1 - 1e-6 do problema quadrático do lead em E real.
    function slowest_decay(h00, h01, h10, e, info) result(lmax)
        complex(dp), intent(in) :: h00(:, :), h01(:, :), h10(:, :)
        real(dp), intent(in) :: e
        integer, intent(out) :: info
        real(dp) :: lmax
        complex(dp), allocatable :: a(:, :), b(:, :), alpha(:), beta(:), work(:), vdl(:, :), vdr(:, :)
        real(dp), allocatable :: rwork(:)
        complex(dp) :: lw(1)
        integer :: n, i, lwork
        real(dp) :: x
        n = size(h00, 1)
        allocate(a(2 * n, 2 * n), b(2 * n, 2 * n), alpha(2 * n), beta(2 * n), rwork(16 * n), vdl(1, 1), vdr(1, 1))
        a = (0.0_dp, 0.0_dp); b = (0.0_dp, 0.0_dp)
        do i = 1, n
            a(i, n + i) = (1.0_dp, 0.0_dp)
            b(i, i) = (1.0_dp, 0.0_dp)
        end do
        a(n + 1:2 * n, 1:n) = -h10
        a(n + 1:2 * n, n + 1:2 * n) = -h00
        do i = 1, n
            a(n + i, n + i) = a(n + i, n + i) + cmplx(e, 0.0_dp, dp)
        end do
        b(n + 1:2 * n, n + 1:2 * n) = h01
        call zggev('N', 'N', 2 * n, a, 2 * n, b, 2 * n, alpha, beta, vdl, 1, vdr, 1, lw, -1, rwork, info)
        lwork = max(1, int(real(lw(1), dp)))
        allocate(work(lwork))
        call zggev('N', 'N', 2 * n, a, 2 * n, b, 2 * n, alpha, beta, vdl, 1, vdr, 1, work, lwork, rwork, info)
        lmax = 0.0_dp
        if (info /= 0) return
        do i = 1, 2 * n
            if (.not. (abs(beta(i)) > 0.0_dp)) cycle
            x = abs(alpha(i) / beta(i))
            if (x < 1.0_dp - 1.0e-6_dp) lmax = max(lmax, x)
        end do
    end function slowest_decay

    ! ==================================================================
    ! Saídas
    ! ==================================================================

    subroutine trench_write_outputs(inp, res, ierr, errmsg)
        type(trench_input_t),  intent(in)  :: inp
        type(trench_result_t), intent(in)  :: res
        integer,               intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        type(green_data_t) :: gd
        type(structure_t) :: cs
        character(len=PATHL) :: pre
        character(len=256) :: emsg
        integer :: u, g, ik, iw, i, p, ie
        character(len=64) :: nm

        ierr = 0
        errmsg = ''
        pre = inp%output_prefix

        ! ---- green.bin: A por átomo como regiões ----
        gd%kmode = GREEN_KGRID
        gd%nk = res%nk; gd%ne = res%ne; gd%nq = 0; gd%n1 = res%nk; gd%n2 = 1; gd%nvertex = 0
        allocate(gd%kfrac(3, res%nk), gd%kcart(3, res%nk), gd%kdist(res%nk))
        gd%kfrac = 0.0_dp; gd%kcart = 0.0_dp
        gd%kfrac(inp%kz_axis, :) = res%kz
        do ik = 1, res%nk
            gd%kcart(:, ik) = res%kz(ik) * TWO_PI * res%lattice(inp%kz_axis, :) / &
                              dot_product(res%lattice(inp%kz_axis, :), res%lattice(inp%kz_axis, :))
            gd%kdist(ik) = norm2(gd%kcart(:, ik))
        end do
        gd%energy = res%energy
        gd%nreg = res%natc
        allocate(gd%region(res%natc))
        do g = 1, res%natc
            write(nm, '(a,i0,a,i0)') 'c', res%acol(g), '_s', res%asite(g)
            gd%region(g) = nm
        end do
        gd%total = sum(res%a, dim=3)
        gd%amap = res%a
        gd%niter = res%niter
        gd%residual = res%resid
        gd%min_neg_im_gii = res%minneg
        allocate(gd%smin(res%nk))
        gd%smin = 1.0_dp
        call meta_add(gd, 'program', 'sktrench')
        call meta_add(gd, 'energy_zero', trim(inp%energy_zero))
        write(nm, '(es24.16e3)') res%eref
        call meta_add(gd, 'eref_ev', trim(adjustl(nm)))
        write(nm, '(es24.16e3)') inp%eta_ev
        call meta_add(gd, 'eta_ev', trim(adjustl(nm)))
        write(nm, '(i0)') res%ncol
        call meta_add(gd, 'ncol', trim(nm))
        write(nm, '(i0)') res%tc
        call meta_add(gd, 'trench_col', trim(nm))
        call meta_add(gd, 'regions', 'one per central atom: c<column>_s<lead atom>; A = -(1/pi) Im tr G_ii')
        call write_green_file(trim(pre) // '.green.bin', gd, ierr, errmsg)
        if (ierr /= 0) return

        ! ---- bordas ----
        open(newunit=u, file=trim(pre) // '.edges.dat', status='replace', action='write')
        write(u, '(a)') '# sktrench: projected lead edges per kz (eV, relative to eref)'
        write(u, '(a,a)') '# energy_zero = ', trim(inp%energy_zero)
        write(u, '(a,3(1x,es24.16e3))') '# eref vbm_true vbm_grid (absolute, eV) =', res%eref, res%vbm_true, res%vbm_grid
        write(u, '(a,es24.16e3)') '# emid (relative) = ', res%emid
        write(u, '(a)') '# ik kz weight ev ec'
        do ik = 1, res%nk
            write(u, '(i6,4(1x,es24.16e3))') ik, res%kz(ik), res%wk(ik), res%ev(ik), res%ec(ik)
        end do
        close(u)

        if (allocated(res%lam)) then
            open(newunit=u, file=trim(pre) // '.lambda.dat', status='replace', action='write')
            write(u, '(a)') '# sktrench: slowest evanescent decay of the lead, max |lambda| < 1; ik ie kz E lambda'
            do ik = 1, res%nk
                do ie = 1, res%ne
                    write(u, '(2(i6),3(1x,es24.16e3))') ik, ie, res%kz(ik), res%energy(ie), res%lam(ie, ik)
                end do
            end do
            close(u)
        end if

        ! ---- átomos e pesos por janela ----
        open(newunit=u, file=trim(pre) // '.atoms.dat', status='replace', action='write')
        write(u, '(a)') '# sktrench: central atoms (column-major); weights per window = sum_kz w_k sum_E A_i dE'
        do iw = 1, res%nwin
            write(u, '(a,i0,a,a,a,a,a,f9.5,a,a,a,f9.5,a,l1)') '# window ', iw, ' ', trim(res%wname(iw)), &
                ': lo = ', trim(inp%win_lo_ref(iw)), ' + ', inp%win_lo(iw), ', hi = ', trim(inp%win_hi_ref(iw)), &
                ' + ', inp%win_hi(iw), ', lo strict = ', inp%win_lo_strict(iw)
        end do
        write(u, '(a)', advance='no') '# g col off site layer x y z trV'
        do iw = 1, res%nwin
            write(u, '(1x,a)', advance='no') trim(res%wname(iw))
        end do
        write(u, '(a)') ''
        do g = 1, res%natc
            write(u, '(i6,1x,i4,1x,i4,1x,i4,1x,i3,4(1x,es24.16e3))', advance='no') g, res%acol(g), &
                res%acol(g) - res%tc, res%asite(g), res%alayer(g), res%apos(:, g), res%avtr(g)
            write(u, '(*(1x,es24.16e3))') res%w(g, :)
        end do
        close(u)

        ! ---- blocos rho ----
        if (allocated(res%ron)) then
            open(newunit=u, file=trim(pre) // '.rho_onsite.dat', status='replace', action='write')
            write(u, '(a)') '# sktrench: rho_i = Re sum_kz w_k sum_E -(G_ii - G_ii^dag)/(2 pi i) dE (norb x norb)'
            write(u, '(a,i0,a,i0,a,i0)') '# nwin ', res%nwin, ' natc ', res%natc, ' norb ', res%norb_atom
            do iw = 1, res%nwin
                do g = 1, res%natc
                    write(u, '(a,1x,a,1x,i0)') 'block', trim(res%wname(iw)), g
                    do i = 1, res%norb_atom
                        write(u, '(*(1x,es24.16e3))') res%ron(i, :, g, iw)
                    end do
                end do
            end do
            close(u)
            open(newunit=u, file=trim(pre) // '.rho_pairs.dat', status='replace', action='write')
            write(u, '(a)') '# sktrench: rho_{i0,jn} = Re sum_kz w_k e^{-2 pi i kz n} sum_E -(G_ij - G_ji^dag)/(2 pi i) dE'
            write(u, '(a)') '# lattice gauge; the inverse pair is the transpose. Pair list: p gi gj nz Rx Ry Rz dist'
            write(u, '(a,i0,a,i0,a,i0)') '# nwin ', res%nwin, ' npair ', res%npair, ' norb ', res%norb_atom
            do p = 1, res%npair
                write(u, '(a,3(1x,i0),4(1x,es24.16e3))') 'pair', p, res%pi(p), res%pj(p), res%pnz(p), &
                    res%pvec(:, p), norm2(res%pvec(:, p))
            end do
            do iw = 1, res%nwin
                do p = 1, res%npair
                    write(u, '(a,1x,a,1x,i0)') 'block', trim(res%wname(iw)), p
                    do i = 1, res%norb_atom
                        write(u, '(*(1x,es24.16e3))') res%rpar(i, :, p, iw)
                    end do
                end do
            end do
            close(u)
        end if

        ! ---- estrutura da região central (para atomselect e figuras) ----
        cs%nat = res%natc
        cs%lattice = res%lattice
        cs%lattice(inp%embed_axis, :) = real(res%ncol, dp) * res%lattice(inp%embed_axis, :)
        cs%cart = res%apos
        allocate(cs%label(res%natc))
        cs%label = 'X'
        call write_structure(trim(pre) // '.central.structure.nml', cs, ierr, emsg)
        if (ierr /= 0) then
            errmsg = trim(emsg); return
        end if
    end subroutine trench_write_outputs

    ! ==================================================================
    ! Utilitários
    ! ==================================================================

    !> Como numpy.linspace (último ponto exato).
    pure function linspace(a, b, n) result(x)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: n
        real(dp) :: x(n)
        real(dp) :: step
        integer :: i
        if (n == 1) then
            x(1) = a
            return
        end if
        step = (b - a) / real(n - 1, dp)
        do i = 1, n
            x(i) = a + real(i - 1, dp) * step
        end do
        x(n) = b
    end function linspace

    function argsort_desc(x) result(idx)
        real(dp), intent(in) :: x(:)
        integer :: idx(size(x))
        integer :: i, j, t
        idx = [(i, i = 1, size(x))]
        do i = 2, size(x)
            t = idx(i)
            j = i - 1
            do while (j >= 1)
                if (x(idx(j)) >= x(t)) exit
                idx(j + 1) = idx(j)
                j = j - 1
            end do
            idx(j + 1) = t
        end do
    end function argsort_desc

end module trench_green
