! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Região central extraída de um DFT (fase 4 do surfgreen, modo 'extracted').
!>
!> A região central C são as colunas (camadas) c_lo..c_hi, ao longo do eixo de
!> embedding, de uma supercélula DFT (H, S bloco-esparsos do .hsio). Os meios
!> semi-infinitos (leads) são um material periódico ao longo do mesmo eixo, com
!> camada principal (PL) de P células; o acoplamento C <-> lead usa os blocos
!> da PRÓPRIA supercélula entre C e as P colunas vizinhas de cada lado
!> (colunas de costura), casadas por posição com a PL do lead, na ordem dos
!> átomos do lead. O lead é alinhado à supercélula por H + Delta S (Delta dado
!> ou ajustado por mínimos quadrados em colunas da supercélula).
!>
!> Colunas. Seja B a rede da supercélula com o vetor do eixo de embedding e
!> trocado pelo vetor a_e do lead; f = r B^{-1}. A coluna do átomo i é
!> c(i) = floor(f_e(i) - s), com s = column_shift (automático: meio da maior
!> lacuna da distribuição de f_e mod 1). A imagem de j por R = sum n_a a_a da
!> supercélula cai na coluna estendida X = c(j) + ncx n_e, ncx = a_e^sup / a_e
!> (ncx = 0 se o eixo e é de vácuo na supercélula: blocos com n_e /= 0 devem
!> ser nulos). Um "átomo estendido" é o par (j, X).
!>
!> k. Os pontos k são frações da recíproca de B nos eixos periódicos que
!> restam (componentes nos eixos de embedding e de vácuo nulas), de modo que
!> k . a_e = 0 para o lead e para a supercélula. Todos os blocos são montados
!> no gauge atômico, direto da forma bloco-esparsa (sem H(R) denso):
!>   H_{(i,Xi),(j,Xj)}(k) = sum_{pares (i,j,R) com X(j,R) - c(i) = Xj - Xi} H_ij(R) e^{i k.(r_j + R - r_i)}.
!> O gauge atômico usa a posição absoluta; uma translação UNIFORME t entre o
!> referencial do lead e o da supercélula multiplica as funções de Bloch da
!> costura por e^{i k.t}, que se cancela em Sigma = K_CL g K_LC; por isso o
!> casamento aceita, opcionalmente (match_translation = 'auto'), um t único por
!> grupo de colunas casadas, e sempre identifica posições módulo vetores de rede
!> dos eixos de k (a função de Bloch atômica não depende da imagem escolhida).
!>
!> Solução. Os leads: decimação de López-Sancho na forma K (surface_green); com
!> o mesmo material e o mesmo Delta dos dois lados, as duas Greens de superfície
!> saem de uma só decimação (surface_gf_pair). A região central é dividida em
!> blocos de colunas de largura >= alcance do acoplamento e resolvida por RGF
!> (central_green). DOS de Mulliken por orbital com os termos cruzados com os
!> blocos vizinhos e com as PLs de costura (G_{C,L} S_{L,C}); lead_cross =
!> .false. os omite só para comparação com referências antigas.
!>
!> Rotinas sem stop; erros por ierr/errmsg.
module extracted_green
    use, intrinsic :: iso_fortran_env, only: iostat_end, int64
    use types, only: dp
    use HS_io, only: HSBlockSparse, read_hs_blocksparse
    use structure_io, only: structure_t, read_structure, structure_fingerprint, LABEL_LEN
    use kpath, only: kpath_t, kgrid_t, build_kpath, build_kgrid
    use atomic_regions, only: region_set_t, read_regions_file, validate_regions
    use spectral_projection, only: orbital_groups_t, build_orbital_groups
    use surface_green, only: surface_gf, surface_gf_pair, SIDE_TOP, SIDE_BOTTOM
    use central_green, only: central_t, rgf_t, central_init, assemble_central, attach_lead, rgf_solve
    use green_projection, only: orbital_dos_mulliken, orbital_dos_lowdin, region_dos, &
                                min_neg_im_gdiag, central_overlap_min
    use green_data, only: green_data_t, meta_add, GREEN_KPATH, GREEN_KGRID
    use site_potential, only: site_potential_t, read_site_potential, site_potential_find
    use pao_basis, only: pao_species_t, pao_dip_cache_t, pao_dip_ensure, pao_dip_get, pao_dipole_dh
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    implicit none
    private

    integer, parameter :: XPATH_LEN = 512, MAX_XVERTEX = 64, MAX_FITK = 32, MAX_FITCOL = 64, MAX_KLIST = 1024
    integer, parameter :: NONE_INT = -huge(0)
    real(dp), parameter :: AUTO_SHIFT = huge(1.0_dp)
    real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp
    real(dp), parameter :: BOHR_ANG_X = 0.529177210903_dp

    !> Espelho de &GREEN (mode = 'extracted'), &EXTRACTED, &LEAD, &LEAD_RIGHT,
    !> &KPATH ou &KGRID.
    type, public :: extracted_input_t
        ! &GREEN
        real(dp)                  :: eta_ev = 0.002_dp
        real(dp)                  :: emin = -1.0_dp, emax = 1.0_dp
        integer                   :: ne = 101
        real(dp)                  :: decim_tol = 1.0e-10_dp
        integer                   :: decim_maxit = 200
        character(len=16)         :: projection = 'mulliken'
        character(len=XPATH_LEN)  :: output_prefix = 'surfgreen'
        character(len=XPATH_LEN)  :: regions_file = ''
        logical                   :: allow_no_fingerprint = .true.
        real(dp)                  :: vacuum_ang = 20.0_dp
        real(dp)                  :: lattice_tol_ang = 0.05_dp
        ! &EXTRACTED (supercélula)
        character(len=XPATH_LEN)  :: h_file = '', s_file = '', structure_file = ''
        integer                   :: embed_axis = 1
        integer                   :: vacuum_axes(2) = 0
        integer                   :: central_columns(2) = NONE_INT
        real(dp)                  :: column_shift = AUTO_SHIFT
        character(len=8)          :: left = 'lead', right = 'lead'
        real(dp)                  :: match_tol_ang = 0.05_dp
        character(len=8)          :: match_translation = 'none'
        character(len=8)          :: delta_mode = 'given'
        real(dp)                  :: delta_ev = 0.0_dp, delta_right_ev = 0.0_dp
        integer                   :: nfit_col = 0, nfit_col_right = 0
        integer                   :: fit_columns(MAX_FITCOL) = NONE_INT
        integer                   :: fit_columns_right(MAX_FITCOL) = NONE_INT
        integer                   :: fit_neighbors = 0
        integer                   :: nfit_k = 0
        real(dp)                  :: fit_k(3, MAX_FITK) = 0.0_dp
        character(len=16)         :: energy_origin = 'absolute'
        logical                   :: lead_cross = .true.
        character(len=8)          :: region_mode = 'columns'
        real(dp)                  :: vacuum_tol = 1.0e-8_dp
        logical                   :: overlap_diagnostic = .true.
        character(len=12)         :: k_basis = 'lead'     !< 'lead' | 'supercell' (só com eixo e de vácuo)
        ! &SITEPOT (opcional): potencial on-site por átomo da região central
        character(len=XPATH_LEN)  :: sitepot_file = ''
        character(len=16)         :: sitepot_mode = 'anticommutator'   !< ou 'onsite'
        logical                   :: sitepot_use_tag = .false.        !< .true.: tag = coluna estendida X
        real(dp)                  :: sitepot_tol = 1.0e-3_dp           !< Å
        real(dp)                  :: sitepot_shift(3) = 0.0_dp         !< somado à posição do átomo antes da busca
        real(dp)                  :: sitepot_uniform_ev = 0.0_dp       !< escalar somado a todos os átomos centrais
        real(dp)                  :: sitepot_lead_ev = 0.0_dp          !< V dos átomos da costura (lado do lead)
        logical                   :: have_sitepot = .false.
        ! &RHO (opcional): matriz densidade por pares e Mulliken por janela
        logical                   :: have_rho = .false.
        integer                   :: rho_nwin = 0
        character(len=16)         :: rho_name(16) = ''
        real(dp)                  :: rho_lo(16) = 0.0_dp, rho_hi(16) = 0.0_dp
        real(dp)                  :: rho_rmax = 3.9_dp
        logical                   :: rho_contour = .false.   ! janelas por contorno no semiplano superior
        integer                   :: rho_ncontour = 48       ! pontos de Gauss-Legendre por janela
        real(dp)                  :: rho_contour_eta = 0.0_dp ! deslocamento imaginário do contorno (eV; 0 = eixo real exato)
        ! &LEAD (1) e &LEAD_RIGHT (2, opcional: sem ele o lado direito usa &LEAD)
        logical                   :: have_lead_right = .false.
        character(len=XPATH_LEN)  :: lead_h(2) = '', lead_s(2) = '', lead_structure(2) = ''
        integer                   :: lead_pl_cells(2) = 0
        real(dp)                  :: lead_column_shift(2) = AUTO_SHIFT
        ! &KPATH / &KGRID
        character(len=8)          :: kmode = 'path'
        integer                   :: nvertex = 0
        real(dp)                  :: vertex(3, MAX_XVERTEX) = 0.0_dp
        character(len=16)         :: label(MAX_XVERTEX) = ''
        integer                   :: npoints(MAX_XVERTEX - 1) = 0
        real(dp)                  :: grid_origin(3) = 0.0_dp, grid_span1(3) = 0.0_dp, grid_span2(3) = 0.0_dp
        integer                   :: grid_n1 = 0, grid_n2 = 0
        ! &KLIST: média uniforme num eixo e lista com pesos no outro
        integer                   :: klist_avg_axis = 0, klist_navg = 0, klist_axis = 0, klist_n = 0
        real(dp)                  :: klist_k(MAX_KLIST) = 0.0_dp, klist_w(MAX_KLIST) = 1.0_dp
    end type extracted_input_t

    !> Sistema bloco-esparso com geometria e colunas (supercélula ou lead).
    type :: xsys_t
        integer :: nat = 0, norb = 0
        integer, allocatable :: orbs(:), aoff(:)              !< (nat), (nat+1)
        character(len=LABEL_LEN), allocatable :: label(:)
        real(dp) :: lattice(3, 3) = 0.0_dp                    !< linhas a_i (Å)
        real(dp) :: basis(3, 3) = 0.0_dp, binv(3, 3) = 0.0_dp !< B (linhas) e B^{-1}
        real(dp) :: recip(3, 3) = 0.0_dp                      !< recíproca de B (linhas, com 2 pi)
        real(dp), allocatable :: cart(:, :), fe(:)
        integer, allocatable :: col(:)
        real(dp) :: shift = 0.0_dp                            !< column_shift usado
        integer :: ncx = 0                                    !< 0: eixo e sem periodicidade
        integer :: axis = 1
        integer :: npair = 0
        integer, allocatable :: pi(:), pj(:), pdx(:), pne(:)  !< átomos, desl. de coluna, n_e
        real(dp), allocatable :: pvec(:, :)                   !< (3,npair) r_j + R - r_i
        integer(int64), allocatable :: hoff(:), soff(:)       !< -1: bloco ausente (nulo)
        complex(dp), allocatable :: hdat(:), sdat(:)
        integer, allocatable :: pstart(:)                     !< CSR por átomo da linha (nat+1)
        integer :: range = 0                                  !< max |pdx|
        logical :: has_meta = .false., shifted = .false.
        real(dp) :: fermi_ev = 0.0_dp
        type(structure_t) :: st
    end type xsys_t

    !> Lista ordenada de átomos estendidos (j, X) e offsets de orbital.
    !> Potencial on-site (&SITEPOT) de um bloco da RGF: escalar por orbital (v),
    !> por átomo (av, com aoff/aorb) e blocos orbitais somados ao on-site (a).
    type :: vblock_t
        real(dp), allocatable :: v(:), av(:), a(:, :)
        integer, allocatable :: aoff(:), aorb(:)
        logical :: has_block = .false.
    end type vblock_t

    type :: xlist_t
        integer :: n = 0, norb = 0
        integer, allocatable :: atom(:), x(:), off(:)          !< off(n+1)
    end type xlist_t

    !> Plano de montagem de um bloco H(k), S(k) entre duas listas.
    type :: plan_t
        integer :: nr = 0, nc = 0, nent = 0
        integer, allocatable :: pair(:), ro(:), co(:)
        integer, allocatable :: ri(:), ci(:)                   !< item da linha/coluna nas listas (A6)
    end type plan_t

    !> Termo de dipolo (A6) de um plano: blocos reais por entrada, mesmo leiaute de hdat.
    type :: dplan_t
        integer, allocatable :: off(:)
        real(dp), allocatable :: dat(:)
    end type dplan_t

    !> Um lead resolvido: sistema, PL, Delta.
    type :: lead_t
        type(xsys_t) :: sys
        integer :: P = 0, c0 = 0, norb = 0
        type(xlist_t) :: pl0, pl1
        type(plan_t) :: p00, p01
        real(dp) :: delta = 0.0_dp
    end type lead_t

    public :: input_is_extracted, read_extracted_input, run_extracted

contains

    ! ==================================================================
    ! Entrada
    ! ==================================================================

    !> .true. se o &GREEN do arquivo tem mode = 'extracted'.
    logical function input_is_extracted(filename) result(yes)
        character(len=*), intent(in) :: filename
        integer :: iu, ios
        character(len=16) :: mode, projection
        real(dp) :: eta_ev, emin, emax, decim_tol, vacuum_ang, lattice_tol_ang
        integer :: ne, decim_maxit
        character(len=XPATH_LEN) :: output_prefix, regions_file
        logical :: allow_no_fingerprint
        namelist /GREEN/ mode, eta_ev, emin, emax, ne, decim_tol, decim_maxit, projection, &
                         output_prefix, regions_file, allow_no_fingerprint, vacuum_ang, lattice_tol_ang
        yes = .false.
        mode = ''
        open(newunit=iu, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        read(iu, nml=GREEN, iostat=ios)
        close(iu)
        if (ios == 0) yes = trim(adjustl(mode)) == 'extracted'
    end function input_is_extracted

    !> Lê os namelists do modo 'extracted'. &GREEN, &EXTRACTED e &LEAD (com um
    !> lado 'lead') são obrigatórios; exatamente um entre &KPATH e &KGRID.
    subroutine read_extracted_input(filename, inp, ierr, errmsg)
        character(len=*),        intent(in)  :: filename
        type(extracted_input_t), intent(out) :: inp
        integer,                 intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: iu, ios
        logical :: found, have_path, have_grid, have_list
        character(len=512) :: iomsg

        ierr = 0
        errmsg = ''
        open(newunit=iu, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = "cannot open input file: '" // trim(filename) // "'"; return
        end if
        call read_green_nml(found)
        if (ierr /= 0) goto 900
        if (.not. found) then
            ierr = 2; errmsg = 'namelist &GREEN not found'; goto 900
        end if
        call read_extracted_nml(found)
        if (ierr /= 0) goto 900
        if (.not. found) then
            ierr = 2; errmsg = "namelist &EXTRACTED not found (required with mode = 'extracted')"; goto 900
        end if
        call read_lead_nml(1, found)
        if (ierr /= 0) goto 900
        if (.not. found .and. (trim(inp%left) == 'lead' .or. trim(inp%right) == 'lead')) then
            ierr = 2; errmsg = 'namelist &LEAD not found'; goto 900
        end if
        call read_lead_nml(2, inp%have_lead_right)
        if (ierr /= 0) goto 900
        call read_kpath_nml(have_path)
        if (ierr /= 0) goto 900
        call read_kgrid_nml(have_grid)
        if (ierr /= 0) goto 900
        call read_klist_nml(have_list)
        if (ierr /= 0) goto 900
        if (count([have_path, have_grid, have_list]) /= 1) then
            ierr = 2; errmsg = 'exactly one of &KPATH, &KGRID and &KLIST must be present'; goto 900
        end if
        inp%kmode = merge('grid', 'path', have_grid)
        if (have_list) inp%kmode = 'list'
        call read_sitepot_nml(inp%have_sitepot)
        if (ierr /= 0) goto 900
        call read_rho_nml(inp%have_rho)
        if (ierr /= 0) goto 900
        close(iu)
        call validate_input(inp, ierr, errmsg)
        return
900     close(iu)

    contains

        subroutine nml_error(name)
            character(len=*), intent(in) :: name
            ierr = 2
            errmsg = 'error reading &' // name // ': ' // trim(iomsg)
        end subroutine nml_error

        subroutine read_green_nml(found)
            logical, intent(out) :: found
            character(len=16) :: mode, projection
            real(dp) :: eta_ev, emin, emax, decim_tol, vacuum_ang, lattice_tol_ang
            integer :: ne, decim_maxit
            character(len=XPATH_LEN) :: output_prefix, regions_file
            logical :: allow_no_fingerprint
            namelist /GREEN/ mode, eta_ev, emin, emax, ne, decim_tol, decim_maxit, projection, &
                             output_prefix, regions_file, allow_no_fingerprint, vacuum_ang, lattice_tol_ang
            mode = ''; eta_ev = inp%eta_ev; emin = inp%emin; emax = inp%emax; ne = inp%ne
            decim_tol = inp%decim_tol; decim_maxit = inp%decim_maxit; projection = inp%projection
            output_prefix = inp%output_prefix; regions_file = inp%regions_file
            allow_no_fingerprint = inp%allow_no_fingerprint; vacuum_ang = inp%vacuum_ang
            lattice_tol_ang = inp%lattice_tol_ang
            rewind(iu)
            read(iu, nml=GREEN, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('GREEN'); return
            end if
            if (found .and. trim(adjustl(mode)) /= 'extracted') then
                ierr = 2; errmsg = "&GREEN mode must be 'extracted' for this input"; return
            end if
            inp%eta_ev = eta_ev; inp%emin = emin; inp%emax = emax; inp%ne = ne
            inp%decim_tol = decim_tol; inp%decim_maxit = decim_maxit; inp%projection = projection
            inp%output_prefix = output_prefix; inp%regions_file = regions_file
            inp%allow_no_fingerprint = allow_no_fingerprint; inp%vacuum_ang = vacuum_ang
            inp%lattice_tol_ang = lattice_tol_ang
        end subroutine read_green_nml

        subroutine read_extracted_nml(found)
            logical, intent(out) :: found
            character(len=XPATH_LEN) :: h_file, s_file, structure_file
            integer :: embed_axis, vacuum_axes(2), central_columns(2), fit_neighbors
            integer :: fit_columns(MAX_FITCOL), fit_columns_right(MAX_FITCOL)
            real(dp) :: column_shift, match_tol_ang, delta_ev, delta_right_ev, vacuum_tol
            real(dp) :: fit_k(3, MAX_FITK)
            character(len=8) :: left, right, match_translation, delta_mode, region_mode
            character(len=16) :: energy_origin
            logical :: lead_cross, overlap_diagnostic
            character(len=12) :: k_basis
            integer :: k
            namelist /EXTRACTED/ h_file, s_file, structure_file, embed_axis, vacuum_axes, central_columns, &
                                 column_shift, left, right, match_tol_ang, match_translation, delta_mode, &
                                 delta_ev, delta_right_ev, fit_columns, fit_columns_right, fit_neighbors, &
                                 fit_k, energy_origin, lead_cross, region_mode, vacuum_tol, overlap_diagnostic, &
                                 k_basis
            h_file = ''; s_file = ''; structure_file = ''
            embed_axis = inp%embed_axis; vacuum_axes = 0; central_columns = NONE_INT
            column_shift = AUTO_SHIFT; left = inp%left; right = inp%right
            match_tol_ang = inp%match_tol_ang; match_translation = inp%match_translation
            delta_mode = inp%delta_mode; delta_ev = 0.0_dp; delta_right_ev = huge(1.0_dp)
            fit_columns = NONE_INT; fit_columns_right = NONE_INT; fit_neighbors = 0
            fit_k = huge(1.0_dp); energy_origin = inp%energy_origin; lead_cross = .true.
            region_mode = inp%region_mode; vacuum_tol = inp%vacuum_tol; overlap_diagnostic = .true.
            k_basis = inp%k_basis
            rewind(iu)
            read(iu, nml=EXTRACTED, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('EXTRACTED'); return
            end if
            if (.not. found) return
            inp%h_file = h_file; inp%s_file = s_file; inp%structure_file = structure_file
            inp%embed_axis = embed_axis; inp%vacuum_axes = vacuum_axes; inp%central_columns = central_columns
            inp%column_shift = column_shift; inp%left = adjustl(left); inp%right = adjustl(right)
            inp%match_tol_ang = match_tol_ang; inp%match_translation = adjustl(match_translation)
            inp%delta_mode = adjustl(delta_mode); inp%delta_ev = delta_ev
            inp%delta_right_ev = merge(delta_ev, delta_right_ev, delta_right_ev == huge(1.0_dp))
            inp%nfit_col = count(fit_columns /= NONE_INT)
            inp%fit_columns(1:inp%nfit_col) = pack(fit_columns, fit_columns /= NONE_INT)
            inp%nfit_col_right = count(fit_columns_right /= NONE_INT)
            inp%fit_columns_right(1:inp%nfit_col_right) = pack(fit_columns_right, fit_columns_right /= NONE_INT)
            inp%fit_neighbors = fit_neighbors
            inp%nfit_k = 0
            do k = 1, MAX_FITK
                if (fit_k(1, k) == huge(1.0_dp)) exit
                inp%nfit_k = k
            end do
            where (fit_k == huge(1.0_dp)) fit_k = 0.0_dp
            inp%fit_k = fit_k
            inp%energy_origin = adjustl(energy_origin); inp%lead_cross = lead_cross
            inp%region_mode = adjustl(region_mode); inp%vacuum_tol = vacuum_tol
            inp%overlap_diagnostic = overlap_diagnostic
            inp%k_basis = adjustl(k_basis)
        end subroutine read_extracted_nml

        subroutine read_lead_nml(il, found)
            integer, intent(in)  :: il
            logical, intent(out) :: found
            character(len=XPATH_LEN) :: h_file, s_file, structure_file
            integer :: pl_cells
            real(dp) :: column_shift
            namelist /LEAD/ h_file, s_file, structure_file, pl_cells, column_shift
            namelist /LEAD_RIGHT/ h_file, s_file, structure_file, pl_cells, column_shift
            h_file = ''; s_file = ''; structure_file = ''; pl_cells = 0; column_shift = AUTO_SHIFT
            rewind(iu)
            if (il == 1) then
                read(iu, nml=LEAD, iostat=ios, iomsg=iomsg)
            else
                read(iu, nml=LEAD_RIGHT, iostat=ios, iomsg=iomsg)
            end if
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error(merge('LEAD      ', 'LEAD_RIGHT', il == 1)); return
            end if
            if (.not. found) return
            inp%lead_h(il) = h_file; inp%lead_s(il) = s_file; inp%lead_structure(il) = structure_file
            inp%lead_pl_cells(il) = pl_cells; inp%lead_column_shift(il) = column_shift
        end subroutine read_lead_nml

        subroutine read_kpath_nml(found)
            logical, intent(out) :: found
            integer :: nvertex
            real(dp) :: vertex(3, MAX_XVERTEX)
            character(len=16) :: label(MAX_XVERTEX)
            integer :: npoints(MAX_XVERTEX - 1)
            namelist /KPATH/ nvertex, vertex, label, npoints
            nvertex = 0; vertex = 0.0_dp; label = ''; npoints = 0
            rewind(iu)
            read(iu, nml=KPATH, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('KPATH'); return
            end if
            inp%nvertex = nvertex; inp%vertex = vertex; inp%label = label; inp%npoints = npoints
        end subroutine read_kpath_nml

        subroutine read_kgrid_nml(found)
            logical, intent(out) :: found
            real(dp) :: origin(3), span1(3), span2(3)
            integer :: n1, n2
            namelist /KGRID/ origin, span1, span2, n1, n2
            origin = 0.0_dp; span1 = 0.0_dp; span2 = 0.0_dp; n1 = 0; n2 = 0
            rewind(iu)
            read(iu, nml=KGRID, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('KGRID'); return
            end if
            inp%grid_origin = origin; inp%grid_span1 = span1; inp%grid_span2 = span2
            inp%grid_n1 = n1; inp%grid_n2 = n2
        end subroutine read_kgrid_nml

        subroutine read_klist_nml(found)
            logical, intent(out) :: found
            integer :: avg_axis, navg, list_axis, nlist
            real(dp) :: list_k(MAX_KLIST), list_w(MAX_KLIST)
            namelist /KLIST/ avg_axis, navg, list_axis, nlist, list_k, list_w
            avg_axis = 0; navg = 1; list_axis = 0; nlist = 0; list_k = 0.0_dp; list_w = 1.0_dp
            rewind(iu)
            read(iu, nml=KLIST, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('KLIST'); return
            end if
            inp%klist_avg_axis = avg_axis; inp%klist_navg = navg; inp%klist_axis = list_axis
            inp%klist_n = nlist; inp%klist_k = list_k; inp%klist_w = list_w
        end subroutine read_klist_nml

        subroutine read_rho_nml(found)
            logical, intent(out) :: found
            integer :: nwin
            character(len=16) :: name(16)
            real(dp) :: lo(16), hi(16), rmax
            logical :: contour
            integer :: npoints
            real(dp) :: contour_eta
            namelist /RHO/ nwin, name, lo, hi, rmax, contour, npoints, contour_eta
            nwin = 0; name = ''; lo = 0.0_dp; hi = 0.0_dp; rmax = 3.9_dp; contour = .false.; npoints = 48
            contour_eta = 0.0_dp
            rewind(iu)
            read(iu, nml=RHO, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('RHO'); return
            end if
            if (.not. found) return
            if (nwin < 1 .or. nwin > 16) then
                ierr = 2; errmsg = '&RHO nwin must be in 1..16'; return
            end if
            inp%rho_nwin = nwin; inp%rho_name = name; inp%rho_lo = lo; inp%rho_hi = hi; inp%rho_rmax = rmax
            inp%rho_contour = contour; inp%rho_ncontour = npoints; inp%rho_contour_eta = contour_eta
            if (contour .and. npoints < 2) then
                ierr = 2; errmsg = '&RHO npoints must be >= 2'; return
            end if
        end subroutine read_rho_nml

        subroutine read_sitepot_nml(found)
            logical, intent(out) :: found
            character(len=XPATH_LEN) :: file
            character(len=16) :: mode
            logical :: use_tag
            real(dp) :: tol, shift(3), uniform_ev, lead_ev
            namelist /SITEPOT/ file, mode, use_tag, tol, shift, uniform_ev, lead_ev
            file = ''; mode = 'anticommutator'; use_tag = .false.; tol = 1.0e-3_dp; shift = 0.0_dp
            uniform_ev = 0.0_dp; lead_ev = 0.0_dp
            rewind(iu)
            read(iu, nml=SITEPOT, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('SITEPOT'); return
            end if
            if (.not. found) return
            inp%sitepot_file = file; inp%sitepot_mode = adjustl(mode); inp%sitepot_use_tag = use_tag
            inp%sitepot_tol = tol; inp%sitepot_shift = shift; inp%sitepot_uniform_ev = uniform_ev
            inp%sitepot_lead_ev = lead_ev
        end subroutine read_sitepot_nml

    end subroutine read_extracted_input

    subroutine validate_input(inp, ierr, errmsg)
        type(extracted_input_t), intent(inout) :: inp
        integer,                 intent(out)   :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        logical :: kaxis(3)
        integer :: a, s
        ierr = 2
        errmsg = ''
        if (inp%eta_ev <= 0.0_dp) then
            errmsg = '&GREEN eta_ev must be > 0'; return
        end if
        if (inp%ne < 1 .or. (inp%ne > 1 .and. inp%emax <= inp%emin)) then
            errmsg = '&GREEN needs ne >= 1 and emax > emin'; return
        end if
        if (inp%decim_tol <= 0.0_dp .or. inp%decim_maxit < 1) then
            errmsg = '&GREEN decim_tol and decim_maxit must be positive'; return
        end if
        if (trim(inp%projection) /= 'mulliken' .and. trim(inp%projection) /= 'lowdin') then
            errmsg = "&GREEN projection must be 'mulliken' or 'lowdin'"; return
        end if
        if (len_trim(inp%h_file) == 0 .or. len_trim(inp%structure_file) == 0) then
            errmsg = '&EXTRACTED h_file and structure_file are required'; return
        end if
        if (inp%embed_axis < 1 .or. inp%embed_axis > 3) then
            errmsg = '&EXTRACTED embed_axis must be 1, 2 or 3'; return
        end if
        if (any(inp%vacuum_axes < 0) .or. any(inp%vacuum_axes > 3)) then
            errmsg = '&EXTRACTED vacuum_axes must be in 1..3 (0 = unused)'; return
        end if
        if (inp%vacuum_axes(1) > 0 .and. inp%vacuum_axes(1) == inp%vacuum_axes(2)) then
            errmsg = '&EXTRACTED vacuum_axes repeated'; return
        end if
        kaxis = .true.
        kaxis(inp%embed_axis) = .false.
        do a = 1, 2
            if (inp%vacuum_axes(a) > 0) kaxis(inp%vacuum_axes(a)) = .false.
        end do
        if (count(inp%vacuum_axes > 0 .and. inp%vacuum_axes /= inp%embed_axis) == 2) then
            errmsg = '&EXTRACTED: no periodic axis left for k'; return
        end if
        if (any(inp%central_columns == NONE_INT)) then
            errmsg = '&EXTRACTED central_columns = c_lo, c_hi is required'; return
        end if
        if (inp%central_columns(2) < inp%central_columns(1)) then
            errmsg = '&EXTRACTED central_columns: c_hi < c_lo'; return
        end if
        do s = 1, 2
            select case (trim(merge(inp%left, inp%right, s == 1)))
            case ('lead', 'vacuum', 'cut')
            case default
                errmsg = "&EXTRACTED left/right must be 'lead', 'vacuum' or 'cut'"; return
            end select
        end do
        ! A7: os dois lados sem lead (vacuum/cut) = região central isolada (slab livre); o &LEAD continua
        ! obrigatório (vetores da rede, E_F e Delta), mas não entra em autoenergia nem costura.
        if (len_trim(inp%lead_h(1)) == 0 .or. len_trim(inp%lead_structure(1)) == 0) then
            errmsg = '&LEAD h_file and structure_file are required'; return
        end if
        if (inp%have_lead_right) then
            if (len_trim(inp%lead_h(2)) == 0 .or. len_trim(inp%lead_structure(2)) == 0) then
                errmsg = '&LEAD_RIGHT h_file and structure_file are required'; return
            end if
        end if
        if (any(inp%lead_pl_cells < 0)) then
            errmsg = 'lead pl_cells must be >= 0 (0 = automatic)'; return
        end if
        if (trim(inp%match_translation) /= 'none' .and. trim(inp%match_translation) /= 'auto') then
            errmsg = "&EXTRACTED match_translation must be 'none' or 'auto'"; return
        end if
        if (inp%match_tol_ang <= 0.0_dp) then
            errmsg = '&EXTRACTED match_tol_ang must be > 0'; return
        end if
        select case (trim(inp%delta_mode))
        case ('given')
        case ('fit')
            if (inp%nfit_k < 1) then
                inp%nfit_k = 1
                inp%fit_k(:, 1) = 0.0_dp
            end if
            if (inp%fit_neighbors < 0) then
                errmsg = '&EXTRACTED fit_neighbors must be >= 0'; return
            end if
        case default
            errmsg = "&EXTRACTED delta_mode must be 'given' or 'fit'"; return
        end select
        do s = 1, inp%nfit_k
            if (any(abs(inp%fit_k(:, s)) > 0.0_dp .and. .not. kaxis)) then
                errmsg = '&EXTRACTED fit_k must have zero components along the embedding and vacuum axes'; return
            end if
        end do
        select case (trim(inp%energy_origin))
        case ('absolute', 'lead_fermi', 'supercell_fermi')
        case default
            errmsg = "&EXTRACTED energy_origin must be 'absolute', 'lead_fermi' or 'supercell_fermi'"; return
        end select
        if (trim(inp%k_basis) /= 'lead' .and. trim(inp%k_basis) /= 'supercell') then
            errmsg = "&EXTRACTED k_basis must be 'lead' or 'supercell'"; return
        end if
        if (trim(inp%k_basis) == 'supercell' .and. .not. any(inp%vacuum_axes == inp%embed_axis)) then
            errmsg = "&EXTRACTED k_basis = 'supercell' needs the embedding axis in vacuum_axes"; return
        end if
        if (trim(inp%region_mode) /= 'columns' .and. trim(inp%region_mode) /= 'atoms') then
            errmsg = "&EXTRACTED region_mode must be 'columns' or 'atoms'"; return
        end if
        if (trim(inp%sitepot_mode) /= 'anticommutator' .and. trim(inp%sitepot_mode) /= 'onsite') then
            errmsg = "&SITEPOT mode must be 'anticommutator' or 'onsite'"; return
        end if
        if (inp%kmode == 'path') then
            if (inp%nvertex < 2 .or. inp%nvertex > MAX_XVERTEX) then
                errmsg = '&KPATH nvertex must be in 2..64'; return
            end if
            if (any(inp%npoints(1:inp%nvertex - 1) < 1)) then
                errmsg = '&KPATH npoints must be >= 1 for every segment'; return
            end if
            do a = 1, 3
                if (.not. kaxis(a) .and. any(abs(inp%vertex(a, 1:inp%nvertex)) > 0.0_dp)) then
                    errmsg = '&KPATH vertices must have zero components along the embedding and vacuum axes'
                    return
                end if
            end do
        else if (inp%kmode == 'list') then
            if (inp%klist_n < 1 .or. inp%klist_n > MAX_KLIST .or. inp%klist_navg < 1) then
                errmsg = '&KLIST needs 1 <= nlist <= 1024 and navg >= 1'; return
            end if
            if (inp%klist_axis < 1 .or. inp%klist_axis > 3) then
                errmsg = '&KLIST list_axis must be 1, 2 or 3'; return
            end if
            if (.not. kaxis(inp%klist_axis)) then
                errmsg = '&KLIST list_axis must be a periodic k axis (not embedding/vacuum)'; return
            end if
            if (inp%klist_avg_axis /= 0) then
                if (inp%klist_avg_axis < 1 .or. inp%klist_avg_axis > 3 .or. inp%klist_avg_axis == inp%klist_axis) then
                    errmsg = '&KLIST avg_axis must be 0 or another axis'; return
                end if
                if (.not. kaxis(inp%klist_avg_axis)) then
                    errmsg = '&KLIST avg_axis must be a periodic k axis (not embedding/vacuum)'; return
                end if
            end if
            if (any(inp%klist_w(1:inp%klist_n) < 0.0_dp) .or. sum(inp%klist_w(1:inp%klist_n)) <= 0.0_dp) then
                errmsg = '&KLIST weights must be non-negative with a positive sum'; return
            end if
        else
            if (inp%grid_n1 < 2 .or. inp%grid_n2 < 2) then
                errmsg = '&KGRID n1 and n2 must be >= 2'; return
            end if
            do a = 1, 3
                if (.not. kaxis(a) .and. (abs(inp%grid_origin(a)) > 0.0_dp .or. abs(inp%grid_span1(a)) > 0.0_dp &
                    .or. abs(inp%grid_span2(a)) > 0.0_dp)) then
                    errmsg = '&KGRID origin/span must have zero components along the embedding and vacuum axes'
                    return
                end if
            end do
        end if
        ierr = 0
    end subroutine validate_input

    ! ==================================================================
    ! Sistemas bloco-esparsos
    ! ==================================================================

    !> Carrega H, S (s_file = '' => S = I) e a estrutura; colunas ao longo de
    !> axis com B = rede com a linha axis trocada por ae; ncx = 0 se o eixo é de
    !> vácuo (blocos com n_axis /= 0 devem ser nulos). vac: eixos de vácuo.
    subroutine load_xsys(hfile, sfile, stfile, axis, ae, vac, embed_is_vacuum, shift_in, vtol, sys, ierr, errmsg)
        character(len=*), intent(in)  :: hfile, sfile, stfile
        integer,          intent(in)  :: axis, vac(:)
        real(dp),         intent(in)  :: ae(3), shift_in, vtol
        logical,          intent(in)  :: embed_is_vacuum
        type(xsys_t),     intent(out) :: sys
        integer,          intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        type(HSBlockSparse) :: bh, bs
        character(len=512) :: msg
        integer(int64), allocatable :: keyh(:), keys(:), key(:)
        integer, allocatable :: oh(:), os(:), ord(:), srcH(:), srcS(:), rn(:)
        integer(int64) :: nh, ns, k, ih, is, n
        integer :: i, j, r, a, nat, p, zero_r
        real(dp) :: v(3), rv(3), amax
        logical :: have_s, drop
        logical, allocatable :: keep(:)

        ierr = 0
        errmsg = ''
        call read_hs_blocksparse(trim(hfile), bh, ierr, errmsg)
        if (ierr /= 0) return
        have_s = len_trim(sfile) > 0
        if (have_s) then
            call read_hs_blocksparse(trim(sfile), bs, ierr, errmsg)
            if (ierr /= 0) return
            if (bs%nAtoms /= bh%nAtoms .or. any(bs%Orbs /= bh%Orbs)) then
                ierr = 2; errmsg = "H and S files differ in atoms/orbitals: '" // trim(hfile) // "'"; return
            end if
            if (bs%Neighbors /= bh%Neighbors) then
                ierr = 2; errmsg = 'H and S files have different cell lists'; return
            end if
            if (any(bs%IndexList /= bh%IndexList)) then
                ierr = 2; errmsg = 'H and S files have different cell lists'; return
            end if
        end if
        call read_structure(trim(stfile), sys%st, ierr, msg)
        if (ierr /= 0) then
            errmsg = 'structure: ' // trim(msg); return
        end if
        nat = bh%nAtoms
        if (sys%st%nat /= nat) then
            ierr = 2
            errmsg = "structure '" // trim(stfile) // "' has " // itoa(sys%st%nat) // ' atoms but ' // &
                     trim(hfile) // ' has ' // itoa(nat)
            return
        end if
        sys%nat = nat
        sys%orbs = bh%Orbs
        allocate(sys%aoff(nat + 1))
        sys%aoff(1) = 0
        do i = 1, nat
            sys%aoff(i + 1) = sys%aoff(i) + sys%orbs(i)
        end do
        sys%norb = sys%aoff(nat + 1)
        sys%label = sys%st%label
        sys%cart = sys%st%cart
        sys%lattice = sys%st%lattice
        sys%axis = axis
        sys%has_meta = bh%has_meta
        sys%shifted = bh%fermi_shift_applied
        sys%fermi_ev = bh%fermi_ev
        sys%basis = sys%lattice
        sys%basis(axis, :) = ae
        sys%binv = inv3(sys%basis)
        sys%recip = recip_of(sys%basis)
        if (abs(det3(sys%basis)) < 1.0e-8_dp) then
            ierr = 2; errmsg = 'embedding vector is coplanar with the other lattice vectors'; return
        end if

        ! colunas
        allocate(sys%fe(nat), sys%col(nat))
        do i = 1, nat
            v = matmul(sys%cart(:, i), sys%binv)
            sys%fe(i) = v(axis)
        end do
        if (shift_in == AUTO_SHIFT) then
            sys%shift = largest_gap_shift(sys%fe)
        else
            sys%shift = shift_in
        end if
        do i = 1, nat
            sys%col(i) = floor(sys%fe(i) - sys%shift)
        end do

        ! periodicidade ao longo do eixo
        if (embed_is_vacuum) then
            sys%ncx = 0
        else
            v = matmul(sys%lattice(axis, :), sys%binv)
            sys%ncx = nint(v(axis))
            if (sys%ncx < 1 .or. abs(v(axis) - real(sys%ncx, dp)) > 1.0e-4_dp .or. &
                any(abs(v(other_axes(axis))) > 1.0e-4_dp)) then
                ierr = 2
                errmsg = "lattice vector a" // itoa(axis) // " of '" // trim(stfile) // &
                         "' is not an integer multiple of the lead embedding vector (list the axis in vacuum_axes " // &
                         'if it has vacuum)'
                return
            end if
        end if

        ! pares: união dos blocos de H e S por chave (célula, i, j)
        nh = bh%nBlocks
        ns = 0
        if (have_s) ns = bs%nBlocks
        allocate(keyh(nh))
        do k = 1, nh
            keyh(k) = blk_key(bh%BlockIndex(:, k), nat)
        end do
        call stable_order64(keyh, oh)
        if (have_s) then
            allocate(keys(ns))
            do k = 1, ns
                keys(k) = blk_key(bs%BlockIndex(:, k), nat)
            end do
            call stable_order64(keys, os)
        else
            ! S = I: blocos identidade (i, i, R = 0)
            zero_r = 0
            do r = 1, bh%Neighbors
                if (all(bh%IndexList(:, r) == 0)) zero_r = r
            end do
            if (zero_r == 0) then
                ierr = 2; errmsg = 'H file has no R = 0 cell'; return
            end if
            ns = nat
            allocate(keys(ns), os(ns))
            do i = 1, nat
                keys(i) = blk_key([i, i, zero_r], nat)
                os(i) = i
            end do
        end if
        ! fusão das duas listas ordenadas
        allocate(key(nh + ns), srcH(nh + ns), srcS(nh + ns))
        n = 0; ih = 1; is = 1
        do while (ih <= nh .or. is <= ns)
            n = n + 1
            srcH(n) = 0; srcS(n) = 0
            if (is > ns) then
                key(n) = keyh(oh(ih)); srcH(n) = oh(ih); ih = ih + 1
            else if (ih > nh) then
                key(n) = keys(os(is)); srcS(n) = os(is); is = is + 1
            else if (keyh(oh(ih)) == keys(os(is))) then
                key(n) = keyh(oh(ih)); srcH(n) = oh(ih); srcS(n) = os(is); ih = ih + 1; is = is + 1
            else if (keyh(oh(ih)) < keys(os(is))) then
                key(n) = keyh(oh(ih)); srcH(n) = oh(ih); ih = ih + 1
            else
                key(n) = keys(os(is)); srcS(n) = os(is); is = is + 1
            end if
            if (n > 1) then
                if (key(n) == key(n - 1)) then
                    ierr = 2; errmsg = 'repeated (cell, atom, atom) block in ' // trim(hfile); return
                end if
            end if
        end do

        ! dados (cópia dos blocos usados) e filtro de vácuo
        allocate(keep(n), rn(n))
        do k = 1, n
            if (srcH(k) > 0) then
                rn(k) = bh%BlockIndex(3, srcH(k))
            else if (have_s) then
                rn(k) = bs%BlockIndex(3, srcS(k))
            else
                rn(k) = zero_r
            end if
            drop = .false.
            do a = 1, size(vac)
                if (vac(a) > 0) then
                    if (bh%IndexList(vac(a), rn(k)) /= 0) drop = .true.
                end if
            end do
            if (sys%ncx == 0 .and. bh%IndexList(axis, rn(k)) /= 0) drop = .true.
            keep(k) = .not. drop
            if (drop) then
                amax = 0.0_dp
                if (srcH(k) > 0) amax = maxval(abs(bh%BlockData(bh%BlockOffset(srcH(k)) + 1: &
                                                                 bh%BlockOffset(srcH(k) + 1))))
                if (srcS(k) > 0 .and. have_s) amax = max(amax, maxval(abs(bs%BlockData( &
                    bs%BlockOffset(srcS(k)) + 1:bs%BlockOffset(srcS(k) + 1)))))
                if (amax > vtol) then
                    ierr = 2
                    write(msg, '(A,ES9.2,A,3(1X,I0),A,I0,A,I0)') 'coupling through vacuum (max ', amax, &
                        ') at R =', bh%IndexList(:, rn(k)), ', atoms ', &
                        int(modulo(key(k) / nat, int(nat, int64))) + 1, ' -> ', int(modulo(key(k), int(nat, int64))) + 1
                    errmsg = trim(msg) // '; check vacuum_axes/embed_axis'
                    return
                end if
            end if
        end do
        sys%npair = count(keep)
        allocate(sys%pi(sys%npair), sys%pj(sys%npair), sys%pdx(sys%npair), sys%pne(sys%npair))
        allocate(sys%pvec(3, sys%npair), sys%hoff(sys%npair), sys%soff(sys%npair))
        ! ordenação por átomo da linha (estável: dentro dele, ordem das chaves)
        block
            integer(int64), allocatable :: rkey(:)
            integer, allocatable :: kidx(:)
            allocate(rkey(sys%npair), kidx(sys%npair))
            p = 0
            do k = 1, n
                if (.not. keep(k)) cycle
                p = p + 1
                kidx(p) = int(k)
                rkey(p) = int(modulo(key(k) / nat, int(nat, int64)), int64)   ! átomo i - 1
            end do
            call stable_order64(rkey, ord)
            allocate(sys%hdat(0), sys%sdat(0))
            block
                integer(int64) :: hsz, ssz, hpos, spos
                hsz = 0; ssz = 0
                do p = 1, sys%npair
                    k = kidx(ord(p))
                    if (srcH(k) > 0) hsz = hsz + bh%BlockOffset(srcH(k) + 1) - bh%BlockOffset(srcH(k))
                    if (srcS(k) > 0) then
                        if (have_s) then
                            ssz = ssz + bs%BlockOffset(srcS(k) + 1) - bs%BlockOffset(srcS(k))
                        else
                            ssz = ssz + int(bh%Orbs(srcS(k)), int64)**2
                        end if
                    end if
                end do
                deallocate(sys%hdat, sys%sdat)
                allocate(sys%hdat(hsz), sys%sdat(ssz))
                hpos = 0; spos = 0
                do p = 1, sys%npair
                    k = kidx(ord(p))
                    i = int(modulo(key(k) / nat, int(nat, int64))) + 1
                    j = int(modulo(key(k), int(nat, int64))) + 1
                    sys%pi(p) = i
                    sys%pj(p) = j
                    sys%pne(p) = bh%IndexList(axis, rn(k))
                    sys%pdx(p) = sys%col(j) + sys%ncx * sys%pne(p) - sys%col(i)
                    rv = matmul(real(bh%IndexList(:, rn(k)), dp), sys%lattice)
                    sys%pvec(:, p) = sys%cart(:, j) + rv - sys%cart(:, i)
                    sys%hoff(p) = -1; sys%soff(p) = -1
                    if (srcH(k) > 0) then
                        sys%hoff(p) = hpos
                        associate(b0 => bh%BlockOffset(srcH(k)), b1 => bh%BlockOffset(srcH(k) + 1))
                            sys%hdat(hpos + 1:hpos + b1 - b0) = bh%BlockData(b0 + 1:b1)
                            hpos = hpos + b1 - b0
                        end associate
                    end if
                    if (srcS(k) > 0) then
                        sys%soff(p) = spos
                        if (have_s) then
                            associate(b0 => bs%BlockOffset(srcS(k)), b1 => bs%BlockOffset(srcS(k) + 1))
                                sys%sdat(spos + 1:spos + b1 - b0) = bs%BlockData(b0 + 1:b1)
                                spos = spos + b1 - b0
                            end associate
                        else
                            sys%sdat(spos + 1:spos + int(sys%orbs(i), int64)**2) = (0.0_dp, 0.0_dp)
                            do a = 1, sys%orbs(i)
                                sys%sdat(spos + int(a - 1, int64) * sys%orbs(i) + a) = (1.0_dp, 0.0_dp)
                            end do
                            spos = spos + int(sys%orbs(i), int64)**2
                        end if
                    end if
                end do
            end block
        end block
        allocate(sys%pstart(nat + 1))
        sys%pstart = 0
        do p = 1, sys%npair
            sys%pstart(sys%pi(p) + 1) = sys%pstart(sys%pi(p) + 1) + 1
        end do
        sys%pstart(1) = 1
        do i = 1, nat
            sys%pstart(i + 1) = sys%pstart(i + 1) + sys%pstart(i)
        end do
        sys%range = 0
        if (sys%npair > 0) sys%range = maxval(abs(sys%pdx))
    end subroutine load_xsys

    !> Deslocamento de coluna automático: meio da maior lacuna circular de f mod 1.
    function largest_gap_shift(f) result(s)
        real(dp), intent(in) :: f(:)
        real(dp) :: s
        real(dp), allocatable :: g(:)
        real(dp) :: gap, best, t
        integer :: i, j, n
        n = size(f)
        g = f - floor(f)
        ! ordenação por inserção (n pequeno a moderado)
        do i = 2, n
            t = g(i)
            j = i - 1
            do while (j >= 1)
                if (g(j) <= t) exit
                g(j + 1) = g(j)
                j = j - 1
            end do
            g(j + 1) = t
        end do
        best = -1.0_dp
        s = 0.0_dp
        do i = 1, n
            if (i < n) then
                gap = g(i + 1) - g(i)
                t = 0.5_dp * (g(i + 1) + g(i))
            else
                gap = g(1) + 1.0_dp - g(n)
                t = 0.5_dp * (g(1) + 1.0_dp + g(n))
            end if
            if (gap > best + 1.0e-9_dp) then
                best = gap
                s = t
            end if
        end do
        s = s - floor(s)
    end function largest_gap_shift

    !> O átomo j tem imagem na coluna estendida X? (n = índice da imagem ao longo do eixo)
    logical function ext_ok(sys, j, x, n)
        type(xsys_t), intent(in)  :: sys
        integer,      intent(in)  :: j, x
        integer,      intent(out) :: n
        n = 0
        if (sys%ncx == 0) then
            ext_ok = (sys%col(j) == x)
        else
            ext_ok = (modulo(x - sys%col(j), sys%ncx) == 0)
            if (ext_ok) n = (x - sys%col(j)) / sys%ncx
        end if
    end function ext_ok

    !> Posição cartesiana do átomo estendido (j, X).
    function ext_pos(sys, j, x) result(r)
        type(xsys_t), intent(in) :: sys
        integer,      intent(in) :: j, x
        real(dp) :: r(3)
        integer :: n
        if (ext_ok(sys, j, x, n)) then
            r = sys%cart(:, j) + real(n, dp) * sys%lattice(sys%axis, :)
        else
            r = huge(1.0_dp)
        end if
    end function ext_pos

    !> Todas as imagens de átomos nas colunas x0..x1 (ordem: coluna, átomo).
    subroutine columns_list(sys, x0, x1, lst)
        type(xsys_t),  intent(in)  :: sys
        integer,       intent(in)  :: x0, x1
        type(xlist_t), intent(out) :: lst
        integer :: x, j, n, cnt
        integer, allocatable :: at(:), xx(:)
        allocate(at(sys%nat * max(x1 - x0 + 1, 0)), xx(sys%nat * max(x1 - x0 + 1, 0)))
        cnt = 0
        do x = x0, x1
            do j = 1, sys%nat
                if (ext_ok(sys, j, x, n)) then
                    cnt = cnt + 1
                    at(cnt) = j
                    xx(cnt) = x
                end if
            end do
        end do
        call make_list(sys, at(1:cnt), xx(1:cnt), lst)
    end subroutine columns_list

    subroutine make_list(sys, at, xx, lst)
        type(xsys_t),  intent(in)  :: sys
        integer,       intent(in)  :: at(:), xx(:)
        type(xlist_t), intent(out) :: lst
        integer :: k
        lst%n = size(at)
        lst%atom = at
        lst%x = xx
        allocate(lst%off(lst%n + 1))
        lst%off(1) = 0
        do k = 1, lst%n
            lst%off(k + 1) = lst%off(k) + sys%orbs(at(k))
        end do
        lst%norb = lst%off(lst%n + 1)
    end subroutine make_list

    !> Sublista dos itens k0..k1.
    subroutine sub_list(sys, lst, k0, k1, out)
        type(xsys_t),  intent(in)  :: sys
        type(xlist_t), intent(in)  :: lst
        integer,       intent(in)  :: k0, k1
        type(xlist_t), intent(out) :: out
        call make_list(sys, lst%atom(k0:k1), lst%x(k0:k1), out)
    end subroutine sub_list

    !> Plano de montagem de H(rows, cols): todos os pares (i, j, R) com
    !> (i, Xi) em rows e (j, Xi + dX) em cols.
    subroutine build_plan(sys, rows, cols, plan)
        type(xsys_t),  intent(in)  :: sys
        type(xlist_t), intent(in)  :: rows, cols
        type(plan_t),  intent(out) :: plan
        integer, allocatable :: idx(:, :), tp(:), tr(:), tc(:), tri(:), tci(:)
        integer :: xmin, xmax, k, p, xj, c, cnt
        plan%nr = rows%norb
        plan%nc = cols%norb
        if (rows%n == 0 .or. cols%n == 0) then
            allocate(plan%pair(0), plan%ro(0), plan%co(0), plan%ri(0), plan%ci(0))
            return
        end if
        xmin = minval(cols%x)
        xmax = maxval(cols%x)
        allocate(idx(sys%nat, xmin:xmax))
        idx = 0
        do k = 1, cols%n
            idx(cols%atom(k), cols%x(k)) = k
        end do
        cnt = 0
        do k = 1, rows%n
            cnt = cnt + sys%pstart(rows%atom(k) + 1) - sys%pstart(rows%atom(k))
        end do
        allocate(tp(cnt), tr(cnt), tc(cnt), tri(cnt), tci(cnt))
        cnt = 0
        do k = 1, rows%n
            do p = sys%pstart(rows%atom(k)), sys%pstart(rows%atom(k) + 1) - 1
                xj = rows%x(k) + sys%pdx(p)
                if (xj < xmin .or. xj > xmax) cycle
                c = idx(sys%pj(p), xj)
                if (c == 0) cycle
                cnt = cnt + 1
                tp(cnt) = p
                tr(cnt) = rows%off(k)
                tc(cnt) = cols%off(c)
                tri(cnt) = k
                tci(cnt) = c
            end do
        end do
        plan%nent = cnt
        plan%pair = tp(1:cnt)
        plan%ro = tr(1:cnt)
        plan%co = tc(1:cnt)
        plan%ri = tri(1:cnt)
        plan%ci = tci(1:cnt)
    end subroutine build_plan

    !> h, s = blocos do plano em k, gauge atômico. kf: frações da recíproca de
    !> B; cada sistema usa a SUA recíproca (kc = kf . b_sys), de modo que a fase
    !> de rede é exatamente 2 pi k.n mesmo se as redes do lead e da supercélula
    !> diferem por arredondamento.
    subroutine fill_plan(sys, plan, kf, h, s)
        type(xsys_t), intent(in)  :: sys
        type(plan_t), intent(in)  :: plan
        real(dp),     intent(in)  :: kf(3)
        complex(dp),  intent(out) :: h(:, :), s(:, :)
        integer :: e, p, ni, nj, a, b
        integer(int64) :: o
        complex(dp) :: ph
        real(dp) :: kc(3)
        kc = matmul(kf, sys%recip)
        h = (0.0_dp, 0.0_dp)
        s = (0.0_dp, 0.0_dp)
        do e = 1, plan%nent
            p = plan%pair(e)
            ni = sys%orbs(sys%pi(p))
            nj = sys%orbs(sys%pj(p))
            ph = exp(cmplx(0.0_dp, dot_product(kc, sys%pvec(:, p)), dp))
            if (sys%hoff(p) >= 0) then
                o = sys%hoff(p)
                do b = 1, nj
                    do a = 1, ni
                        h(plan%ro(e) + a, plan%co(e) + b) = h(plan%ro(e) + a, plan%co(e) + b) + &
                            ph * sys%hdat(o + int(b - 1, int64) * ni + a)
                    end do
                end do
            end if
            if (sys%soff(p) >= 0) then
                o = sys%soff(p)
                do b = 1, nj
                    do a = 1, ni
                        s(plan%ro(e) + a, plan%co(e) + b) = s(plan%ro(e) + a, plan%co(e) + b) + &
                            ph * sys%sdat(o + int(b - 1, int64) * ni + a)
                    end do
                end do
            end if
        end do
    end subroutine fill_plan

    !> Casa as colunas x0..x0+ncol-1 da supercélula com as células 0..ncol-1 do
    !> lead (posição do átomo m na célula a: r_m + a a_e). Devolve a lista da
    !> supercélula na ordem do lead (célula, átomo), a translação t e a maior
    !> distância residual. Posições comparadas módulo vetores de rede dos eixos
    !> de k (kaxis); em mode 'none' t só pode ser um múltiplo inteiro de a_e
    !> (o rótulo de coluna é arbitrário, e k.a_e = 0); em 'auto', qualquer t
    !> uniforme (um por chamada).
    subroutine match_columns(sup, lead, x0, ncol, kaxis, tol, mode, lst, t, dmax, ierr, errmsg)
        type(xsys_t),     intent(in)  :: sup
        type(lead_t),     intent(in)  :: lead
        integer,          intent(in)  :: x0, ncol
        logical,          intent(in)  :: kaxis(3)
        real(dp),         intent(in)  :: tol
        character(len=*), intent(in)  :: mode
        type(xlist_t),    intent(out) :: lst
        real(dp),         intent(out) :: t(3), dmax
        integer,          intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        type(xlist_t) :: cand
        integer, allocatable :: at(:), xx(:)
        logical, allocatable :: used(:)
        real(dp) :: tc(3), d, rl(3), f(3)
        integer :: c, nl, a, m, k, best
        logical :: ok

        ierr = 0
        errmsg = ''
        nl = lead%sys%nat
        call columns_list(sup, x0, x0 + ncol - 1, cand)
        if (cand%n /= nl * ncol) then
            ierr = 2
            errmsg = 'columns ' // itoa(x0) // '..' // itoa(x0 + ncol - 1) // ' of the supercell have ' // &
                     itoa(cand%n) // ' atoms, the lead has ' // itoa(nl) // ' per cell (' // &
                     itoa(nl * ncol) // ' expected); check central_columns/column_shift'
            return
        end if
        allocate(at(nl * ncol), xx(nl * ncol), used(cand%n))
        ok = .false.
        do c = 1, cand%n
            if (cand%x(c) /= x0) cycle
            tc = ext_pos(sup, cand%atom(c), cand%x(c)) - lead%sys%cart(:, 1)
            if (trim(mode) == 'none') then
                ! só translações inteiras ao longo do eixo de embedding (k.a_e = 0)
                f = matmul(tc, sup%binv)
                f(sup%axis) = real(nint(f(sup%axis)), dp)
                do a = 1, 3
                    if (kaxis(a)) f(a) = real(nint(f(a)), dp)
                    if (a /= sup%axis .and. .not. kaxis(a)) f(a) = 0.0_dp
                end do
                if (norm2(tc - matmul(f, sup%basis)) > tol) cycle
                f = 0.0_dp
                f(sup%axis) = real(nint(dot_product(tc, sup%binv(:, sup%axis))), dp)
                tc = matmul(f, sup%basis)
            end if
            used = .false.
            dmax = 0.0_dp
            ok = .true.
            do a = 0, ncol - 1
                do m = 1, nl
                    rl = lead%sys%cart(:, m) + real(a, dp) * lead%sys%lattice(lead%sys%axis, :) + tc
                    best = 0
                    do k = 1, cand%n
                        if (used(k) .or. cand%x(k) /= x0 + a) cycle
                        d = reduced_dist(sup, ext_pos(sup, cand%atom(k), cand%x(k)) - rl, kaxis)
                        if (d < tol) then
                            best = k
                            exit
                        end if
                    end do
                    if (best == 0) then
                        ok = .false.
                        exit
                    end if
                    d = reduced_dist(sup, ext_pos(sup, cand%atom(best), cand%x(best)) - rl, kaxis)
                    dmax = max(dmax, d)
                    if (sup%orbs(cand%atom(best)) /= lead%sys%orbs(m)) then
                        ierr = 2
                        errmsg = 'supercell atom ' // itoa(cand%atom(best)) // ' matched to lead atom ' // &
                                 itoa(m) // ' has a different number of orbitals'
                        return
                    end if
                    used(best) = .true.
                    at(a * nl + m) = cand%atom(best)
                    xx(a * nl + m) = cand%x(best)
                end do
                if (.not. ok) exit
            end do
            if (ok) then
                t = tc
                exit
            end if
        end do
        if (.not. ok) then
            ierr = 2
            errmsg = 'supercell columns ' // itoa(x0) // '..' // itoa(x0 + ncol - 1) // &
                     ' do not match the lead principal layer within match_tol_ang'
            if (trim(mode) == 'none') errmsg = errmsg // " (match_translation = 'none'; try 'auto' " // &
                                               'if the two frames differ by a uniform translation)'
            return
        end if
        call make_list(sup, at, xx, lst)
    end subroutine match_columns

    !> |d| após remover vetores de rede dos eixos de k (base B do sistema).
    real(dp) function reduced_dist(sys, d, kaxis) result(r)
        type(xsys_t), intent(in) :: sys
        real(dp),     intent(in) :: d(3)
        logical,      intent(in) :: kaxis(3)
        real(dp) :: f(3)
        integer :: a
        f = matmul(d, sys%binv)
        do a = 1, 3
            if (kaxis(a)) f(a) = f(a) - real(nint(f(a)), dp)
        end do
        r = norm2(matmul(f, sys%basis))
    end function reduced_dist

    !> Carrega o lead e monta a PL (P células) e os planos H00, H01.
    subroutine setup_lead(inp, il, vac, kaxis, lead, ierr, errmsg)
        type(extracted_input_t), intent(in)  :: inp
        integer,                 intent(in)  :: il
        integer,                 intent(in)  :: vac(:)
        logical,                 intent(in)  :: kaxis(3)
        type(lead_t),            intent(out) :: lead
        integer,                 intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        type(structure_t) :: st
        character(len=512) :: msg
        integer :: a, e

        e = inp%embed_axis
        call read_structure(trim(inp%lead_structure(il)), st, ierr, msg)
        if (ierr /= 0) then
            errmsg = 'lead structure: ' // trim(msg); return
        end if
        call load_xsys(inp%lead_h(il), inp%lead_s(il), inp%lead_structure(il), e, st%lattice(e, :), &
                       pack(vac, vac /= e), .false., inp%lead_column_shift(il), inp%vacuum_tol, lead%sys, ierr, errmsg)
        if (ierr /= 0) then
            errmsg = 'lead: ' // errmsg; return
        end if
        if (any(lead%sys%col /= lead%sys%col(1))) then
            ierr = 2
            errmsg = 'lead atoms fall in more than one column along the embedding axis; set &LEAD column_shift'
            return
        end if
        lead%c0 = lead%sys%col(1)
        if (inp%lead_pl_cells(il) == 0) then
            lead%P = max(1, lead%sys%range)
        else
            lead%P = inp%lead_pl_cells(il)
            if (lead%sys%range > lead%P) then
                ierr = 2
                errmsg = 'lead couples cells up to ' // itoa(lead%sys%range) // ' apart along the embedding axis, ' // &
                         'beyond pl_cells = ' // itoa(lead%P)
                return
            end if
        end if
        call columns_list(lead%sys, lead%c0, lead%c0 + lead%P - 1, lead%pl0)
        call columns_list(lead%sys, lead%c0 + lead%P, lead%c0 + 2 * lead%P - 1, lead%pl1)
        lead%norb = lead%pl0%norb
        call build_plan(lead%sys, lead%pl0, lead%pl0, lead%p00)
        call build_plan(lead%sys, lead%pl0, lead%pl1, lead%p01)
        ierr = 0
        errmsg = ''
    end subroutine setup_lead

    ! ==================================================================
    ! Execução
    ! ==================================================================

    !> Executa o modo 'extracted'. structure_only: só grava <prefix>_central.nml.
    !> notes: linhas de resumo (colunas, casamento, Delta) para o terminal.
    !> Ganchos do ciclo autoconsistente (mod_scc, A5), todos opcionais:
    !>   scc_v(ncent)      V escalar por átomo central (eV, ordem de cent), somado ao &SITEPOT na forma do sitepot_mode;
    !>   cent_pos(3,ncent), cent_col(ncent), cent_atom(ncent): posição estendida (Å), coluna e átomo da supercélula;
    !>   rho_q(ncent,nwin) q_DS de &RHO por átomo central; rho_files = .false. não grava rho_pairs/rho_mulliken.
    !> Dipolos (A6), também opcionais (exigem scc_pao e scc_dip):
    !>   scc_f(3,ncent)    campo elétrico por átomo central (eV/Bohr); H_{mu nu} += 1/2 [F_i.X^(i) + F_j.X^(j)],
    !>                     X^(i) = <mu|r - R_i|nu> (integrais nas PAO de scc_pao, cache scc_dip, Bohr); F = 0 fora;
    !>   rho_p(3,ncent,nwin) p_i = -sum_{mu in i, nu} D_{mu nu} X^(i)_{nu mu} (e.Bohr por spin; elétron = -e).
    subroutine run_extracted(inp, gd, ierr, errmsg, structure_only, nonconverged, nthreads, notes, &
                             scc_v, cent_pos, cent_col, cent_atom, rho_q, rho_files, scc_f, scc_pao, scc_dip, rho_p, &
                             scc_vb, bond_ij, bond_d, bond_pos, rho_qb, rho_qon, rho_pon, &
                             rho_pair_ic, rho_pair_ja, rho_pair_kind, rho_pair_vec, rho_pair_s, rho_d, scc_dhp)
        type(extracted_input_t), intent(in)  :: inp
        type(green_data_t),      intent(out) :: gd
        integer,                 intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        logical, optional,       intent(in)  :: structure_only
        integer, optional,       intent(out) :: nonconverged
        integer, optional,       intent(out) :: nthreads
        character(len=256), allocatable, optional, intent(out) :: notes(:)
        real(dp), optional,      intent(in)  :: scc_v(:)
        real(dp), allocatable, optional, intent(out) :: cent_pos(:, :), rho_q(:, :)
        integer, allocatable, optional, intent(out)  :: cent_col(:), cent_atom(:)
        logical, optional,       intent(in)  :: rho_files
        real(dp), optional,      intent(in)  :: scc_f(:, :)
        type(pao_species_t), optional, intent(in) :: scc_pao
        type(pao_dip_cache_t), optional, intent(inout) :: scc_dip
        real(dp), allocatable, optional, intent(out) :: rho_p(:, :, :)
        ! cargas de ligação (A10)
        real(dp), optional,      intent(in)  :: scc_vb(:)
        integer, allocatable, optional, intent(out)  :: bond_ij(:, :)
        real(dp), allocatable, optional, intent(out) :: bond_d(:, :), bond_pos(:, :), rho_qb(:, :), rho_qon(:, :)
        real(dp), allocatable, optional, intent(out) :: rho_pon(:, :, :)
        ! A19 (componente planar pela densidade completa): lista de pares de &RHO (i central, j, d = r_j + L - r_i,
        ! S_ij do arquivo), D_ij real por par (janela 1, por spin) e termo de H real por par somado ao bloco
        ! correspondente (exige scc_f; pares com a costura: o bloco do plano é a transposta)
        integer, allocatable, optional, intent(out)  :: rho_pair_ic(:), rho_pair_ja(:), rho_pair_kind(:)
        real(dp), allocatable, optional, intent(out) :: rho_pair_vec(:, :), rho_pair_s(:, :, :), rho_d(:, :, :)
        real(dp), optional,      intent(in)  :: scc_dhp(:, :, :)

        type(xsys_t) :: sup
        type(lead_t), target :: leads(2)
        type(xlist_t) :: cent, seam(2)
        type(xlist_t), allocatable :: blk(:)
        type(plan_t), allocatable :: pd(:), pu(:), pl(:)
        type(plan_t) :: ps(2)
        type(structure_t) :: stc
        type(region_set_t) :: rs
        type(orbital_groups_t) :: groups
        type(kpath_t) :: path
        type(kgrid_t) :: kg
        character(len=512) :: msg
        character(len=256) :: line
        character(len=256), allocatable :: nt(:)
        character(len=:), allocatable :: central_file, perr_msg
        integer :: e, vac(2), nlead, il, side, nb, ncol, clo, chi, w, b, k, nk, ne, nreg, g, ik, ie, ip
        integer :: nthr, nbad, perr_ip, lmat(2), nfitc
        integer, allocatable :: bsize(:), bwid(:), cblock(:), list_start(:), list_atoms(:)
        logical :: kaxis(3), use_side(2), pair_decim
        real(dp) :: recip(3, 3), t(3), dmax, eorigin, fit_res(2, 2), match_d(2), match_t(3, 2)
        real(dp), allocatable :: fitk_c(:, :)
        integer, allocatable :: fitc(:)
        type(vblock_t), allocatable :: vb(:)
        logical :: have_v, have_dip
        type(dplan_t), allocatable :: dpd(:), dpu(:), dpl(:)
        type(dplan_t) :: dps(2)
        ! &RHO
        integer :: nrp, nmo
        integer, allocatable :: rp_kind(:), rp_b(:), rp_e(:), rp_ic(:), rp_ni(:), rp_nj(:), rp_ja(:), rp_jx(:)
        real(dp), allocatable :: rp_vec(:, :), rp_s(:, :, :), kwt(:)
        complex(dp), allocatable :: rbuf(:, :, :, :, :), qbuf(:, :, :)
        real(dp), allocatable :: glx(:), glw(:)
        real(dp) :: thc
        integer :: w2c, icn
        complex(dp) :: zc, wz
        integer :: nkfull
        real(dp), allocatable :: kw(:)
        ! A10: sítios de ligação (pares de &RHO com d /= 0, as duas ordens num só sítio)
        logical :: bmode
        integer :: nbd
        integer, allocatable :: rp_bond(:), rp_bw(:), rp_jc(:), bd_q(:), qlook(:, :, :)

        ierr = 0
        errmsg = ''
        allocate(nt(0))
        if (present(nonconverged)) nonconverged = 0
        nthr = 1
        !$ nthr = omp_get_max_threads()
        if (present(nthreads)) nthreads = nthr
        e = inp%embed_axis
        vac = inp%vacuum_axes
        kaxis = .true.
        kaxis(e) = .false.
        do k = 1, 2
            if (vac(k) > 0) kaxis(vac(k)) = .false.
        end do
        use_side = [trim(inp%left) == 'lead', trim(inp%right) == 'lead']
        clo = inp%central_columns(1)
        chi = inp%central_columns(2)
        ncol = chi - clo + 1

        ! --- leads ---
        nlead = merge(2, 1, inp%have_lead_right)
        do il = 1, nlead
            call setup_lead(inp, il, vac, kaxis, leads(il), ierr, errmsg)
            if (ierr /= 0) return
        end do
        lmat = [1, nlead]          ! material do lado esquerdo / direito

        ! --- supercélula ---
        call load_xsys(inp%h_file, inp%s_file, inp%structure_file, e, leads(1)%sys%lattice(e, :), vac, &
                       any(vac == e), inp%column_shift, inp%vacuum_tol, sup, ierr, errmsg)
        if (ierr /= 0) then
            errmsg = 'supercell: ' // errmsg; return
        end if
        if (trim(inp%k_basis) == 'supercell') sup%recip = recip_of(sup%lattice)
        do il = 1, nlead
            do k = 1, 3
                if (.not. kaxis(k)) cycle
                if (norm2(sup%lattice(k, :) - leads(il)%sys%lattice(k, :)) > inp%lattice_tol_ang) then
                    ierr = 2
                    errmsg = 'lead and supercell lattice vectors a' // itoa(k) // ' (periodic axis) differ by ' // &
                             'more than lattice_tol_ang'
                    return
                end if
            end do
            if (norm2(leads(il)%sys%lattice(e, :) - leads(1)%sys%lattice(e, :)) > inp%lattice_tol_ang) then
                ierr = 2; errmsg = 'left and right leads have different embedding vectors'; return
            end if
        end do
        line = 'supercell: ' // itoa(sup%nat) // ' atoms, ' // itoa(sup%norb) // ' orbitals, ' // &
               itoa(sup%npair) // ' blocks; columns ' // itoa(minval(sup%col)) // '..' // itoa(maxval(sup%col))
        if (sup%ncx > 0) then
            call note(trim(line) // ' (periodic along the embedding axis, ncx = ' // itoa(sup%ncx) // ')')
        else
            call note(trim(line) // ' (vacuum along the embedding axis)')
        end if
        write(line, '(A,F8.4)') '  column_shift = ', sup%shift
        call note(trim(line) // '; atoms per column:' // column_counts(sup))
        do il = 1, nlead
            call note('lead ' // itoa(il) // ': ' // itoa(leads(il)%sys%nat) // ' atoms/cell, PL = ' // &
                      itoa(leads(il)%P) // ' cells (' // itoa(leads(il)%norb) // ' orbitals), coupling range ' // &
                      itoa(leads(il)%sys%range))
        end do

        ! --- costuras ---
        match_d = 0.0_dp
        match_t = 0.0_dp
        do side = 1, 2
            if (.not. use_side(side)) cycle
            associate(ld => leads(lmat(side)))
                if (side == 1) then
                    call match_columns(sup, ld, clo - ld%P, ld%P, kaxis, inp%match_tol_ang, inp%match_translation, &
                                       seam(side), t, dmax, ierr, errmsg)
                else
                    call match_columns(sup, ld, chi + 1, ld%P, kaxis, inp%match_tol_ang, inp%match_translation, &
                                       seam(side), t, dmax, ierr, errmsg)
                end if
            end associate
            if (ierr /= 0) then
                errmsg = merge('left  seam: ', 'right seam: ', side == 1) // errmsg; return
            end if
            match_d(side) = dmax
            match_t(:, side) = t
            write(line, '(A,A,3F10.5,A,ES9.2,A)') merge('left  seam', 'right seam', side == 1), &
                ': translation (A) ', t, ', max distance ', dmax, ' A'
            call note(trim(line))
        end do

        ! --- Delta ---
        fit_res = 0.0_dp
        if (trim(inp%delta_mode) == 'fit') then
            fitk_c = inp%fit_k(:, 1:inp%nfit_k)
            do il = 1, nlead
                if (il == 1) then
                    if (inp%nfit_col > 0) then
                        fitc = inp%fit_columns(1:inp%nfit_col)
                    else
                        fitc = default_fit_columns(1)
                    end if
                else
                    if (inp%nfit_col_right > 0) then
                        fitc = inp%fit_columns_right(1:inp%nfit_col_right)
                    else
                        fitc = default_fit_columns(2)
                    end if
                end if
                nfitc = size(fitc)
                if (nfitc == 0) then
                    ierr = 2; errmsg = 'delta_mode = fit: no fit columns'; return
                end if
                call fit_delta(sup, leads(il), fitc, inp%fit_neighbors, fitk_c, kaxis, inp%match_tol_ang, &
                               inp%match_translation, leads(il)%delta, fit_res(1, il), fit_res(2, il), ierr, errmsg)
                if (ierr /= 0) then
                    errmsg = 'Delta fit: ' // errmsg; return
                end if
                write(line, '(A,I0,A,ES16.9,A,F9.5,A,ES9.2,A)') 'Delta (lead ', il, ') = ', leads(il)%delta, &
                    ' eV (fit); max onsite |H_sup - H_lead - D S| = ', fit_res(1, il), ' eV, max |dS| = ', &
                    fit_res(2, il), ''
                call note(trim(line))
            end do
        else
            leads(1)%delta = inp%delta_ev
            if (nlead == 2) leads(2)%delta = inp%delta_right_ev
            write(line, '(A,ES16.9,A)') 'Delta = ', leads(1)%delta, ' eV (given)'
            call note(trim(line))
        end if
        pair_decim = use_side(1) .and. use_side(2) .and. nlead == 1

        ! --- origem de energia ---
        select case (trim(inp%energy_origin))
        case ('absolute')
            eorigin = 0.0_dp
        case ('lead_fermi')
            if (.not. leads(1)%sys%has_meta) then
                ierr = 2; errmsg = "energy_origin = 'lead_fermi' needs E_F in the lead .hsio (meta)"; return
            end if
            eorigin = merge(0.0_dp, leads(1)%sys%fermi_ev, leads(1)%sys%shifted) + leads(1)%delta
        case ('supercell_fermi')
            if (.not. sup%has_meta) then
                ierr = 2; errmsg = "energy_origin = 'supercell_fermi' needs E_F in the supercell .hsio (meta)"; return
            end if
            eorigin = merge(0.0_dp, sup%fermi_ev, sup%shifted)
        end select

        ! --- região central e blocos da RGF ---
        call columns_list(sup, clo, chi, cent)
        if (cent%n == 0) then
            ierr = 2; errmsg = 'central region is empty (check central_columns)'; return
        end if
        w = 1
        do side = 1, 2
            if (use_side(side)) w = max(w, leads(lmat(side))%P)
        end do
        block
            integer :: p, kk
            do kk = 1, cent%n
                do p = sup%pstart(cent%atom(kk)), sup%pstart(cent%atom(kk) + 1) - 1
                    if (cent%x(kk) + sup%pdx(p) >= clo .and. cent%x(kk) + sup%pdx(p) <= chi) &
                        w = max(w, abs(sup%pdx(p)))
                end do
            end do
        end block
        nb = max(1, ncol / w)
        allocate(bwid(nb))
        bwid = ncol / nb
        do b = nb - modulo(ncol, nb) + 1, nb
            bwid(b) = bwid(b) + 1
        end do
        allocate(blk(nb), bsize(nb), cblock(cent%n))
        k = clo
        do b = 1, nb
            call columns_list(sup, k, k + bwid(b) - 1, blk(b))
            bsize(b) = blk(b)%norb
            k = k + bwid(b)
        end do
        do k = 1, cent%n
            cblock(k) = block_of(cent%x(k))
        end do
        call check_central_couplings(ierr, errmsg)
        if (ierr /= 0) return
        write(line, '(A,I0,A,I0,A,I0,A,I0,A)') 'central region: columns ', clo, '..', chi, ' (', cent%n, &
            ' atoms, ', cent%norb, ' orbitals)'
        call note(trim(line))
        line = '  RGF blocks (columns):'
        do b = 1, nb
            line = trim(line) // ' ' // itoa(bwid(b))
        end do
        call note(trim(line) // '; orbitals per block:' // ints(bsize))
        if (present(cent_pos)) then
            allocate(cent_pos(3, cent%n))
            do k = 1, cent%n
                cent_pos(:, k) = ext_pos(sup, cent%atom(k), cent%x(k))
            end do
        end if
        if (present(cent_col)) cent_col = cent%x(1:cent%n)
        if (present(cent_atom)) cent_atom = cent%atom(1:cent%n)
        if (present(scc_v)) then
            if (size(scc_v) /= cent%n) then
                ierr = 2; errmsg = 'scc_v: size differs from the number of central atoms'; return
            end if
        end if
        allocate(pd(nb), pu(max(nb - 1, 1)), pl(max(nb - 1, 1)))
        do b = 1, nb
            call build_plan(sup, blk(b), blk(b), pd(b))
            if (b < nb) then
                call build_plan(sup, blk(b), blk(b + 1), pu(b))
                call build_plan(sup, blk(b + 1), blk(b), pl(b))
            end if
        end do
        if (use_side(1)) call build_plan(sup, seam(1), blk(1), ps(1))
        if (use_side(2)) call build_plan(sup, blk(nb), seam(2), ps(2))

        ! --- potencial on-site por átomo (&SITEPOT) ---
        have_v = inp%have_sitepot .or. present(scc_v)
        if (have_v) then
            call setup_sitepot(ierr, errmsg)
            if (ierr /= 0) return
        end if
        nrp = 0
        ! --- cargas de ligação (A10): os sítios saem da lista de pares de &RHO (rho_rmax), mesmo sem &RHO ---
        bmode = present(scc_vb) .or. present(bond_ij) .or. present(rho_qb) .or. present(rho_qon) .or. present(rho_pon)
        if (inp%have_rho .or. bmode .or. present(scc_dhp) .or. present(rho_pair_ic)) call setup_rho_pairs()
        if (present(rho_pair_ic)) rho_pair_ic = rp_ic(1:nrp)
        if (present(rho_pair_ja)) rho_pair_ja = rp_ja(1:nrp)
        if (present(rho_pair_kind)) rho_pair_kind = rp_kind(1:nrp)
        if (present(rho_pair_vec)) rho_pair_vec = rp_vec(:, 1:nrp)
        if (present(rho_pair_s)) rho_pair_s = rp_s(:, :, 1:nrp)
        if (present(scc_dhp)) then
            if (.not. present(scc_f)) then
                ierr = 2; errmsg = 'scc_dhp: needs scc_f (dipole planes)'; return
            end if
            if (size(scc_dhp, 3) /= nrp .or. size(scc_dhp, 1) < nmo .or. size(scc_dhp, 2) < nmo) then
                ierr = 2; errmsg = 'scc_dhp: shape differs from (nmo, nmo, &RHO pairs)'; return
            end if
        end if
        nbd = 0
        if (bmode) then
            call setup_bonds(ierr, errmsg)
            if (ierr /= 0) return
        end if
        ! --- dipolos (A6) ---
        have_dip = present(scc_f)
        if (have_dip .or. present(rho_p) .or. present(rho_pon)) then
            call setup_dipole(ierr, errmsg)
            if (ierr /= 0) return
        end if

        ! --- estrutura central ---
        central_file = trim(inp%output_prefix) // '_central.nml'
        call write_xcentral(central_file, ierr, errmsg)
        if (ierr /= 0) return
        if (present(notes)) notes = nt
        if (present(structure_only)) then
            if (structure_only) return
        end if
        call read_structure(central_file, stc, ierr, msg)
        if (ierr /= 0) then
            errmsg = 'central structure: ' // trim(msg); return
        end if

        ! --- regiões ---
        if (len_trim(inp%regions_file) > 0) then
            call read_regions_file(trim(inp%regions_file), rs, ierr, msg)
            if (ierr == 0) call validate_regions(rs, stc, inp%allow_no_fingerprint, ierr, msg)
            if (ierr /= 0) then
                errmsg = 'regions: ' // trim(msg); return
            end if
            nreg = rs%nregion
            allocate(list_start(nreg + 1))
            list_start(1) = 1
            do g = 1, nreg
                list_start(g + 1) = list_start(g) + size(rs%regions(g)%atoms)
            end do
            allocate(list_atoms(list_start(nreg + 1) - 1))
            do g = 1, nreg
                list_atoms(list_start(g):list_start(g + 1) - 1) = rs%regions(g)%atoms
            end do
        else if (trim(inp%region_mode) == 'atoms') then
            nreg = cent%n
            allocate(list_start(nreg + 1), list_atoms(nreg))
            list_start = [(g, g = 1, nreg + 1)]
            list_atoms = [(g, g = 1, nreg)]
        else
            nreg = ncol
            allocate(list_start(nreg + 1), list_atoms(cent%n))
            list_start(1) = 1
            g = 0
            do k = 1, ncol
                do ik = 1, cent%n
                    if (cent%x(ik) == clo + k - 1) then
                        g = g + 1
                        list_atoms(g) = ik
                    end if
                end do
                list_start(k + 1) = g + 1
            end do
        end if
        call build_orbital_groups(sup%orbs(cent%atom), nreg, list_start, list_atoms, groups)

        ! --- pontos k (frações da recíproca de B) ---
        recip = sup%recip
        gd%nvertex = 0
        if (inp%kmode == 'path') then
            call build_kpath(inp%vertex(:, 1:inp%nvertex), inp%npoints(1:inp%nvertex - 1), &
                             inp%label(1:inp%nvertex), recip, path, ierr, msg)
            if (ierr /= 0) then
                errmsg = 'k-path: ' // trim(msg); return
            end if
            gd%kmode = GREEN_KPATH
            nk = path%nk
            gd%kfrac = path%kfrac
            gd%nvertex = path%nvertex
            gd%vertex_k = path%vertex_k
            gd%label = path%label
            gd%kdist = path%kdist
        else if (inp%kmode == 'list') then
            ! pontos completos: lista (lento) x média uniforme (rápido); reduzidos no fim
            gd%kmode = GREEN_KGRID
            nk = inp%klist_n * inp%klist_navg
            allocate(gd%kfrac(3, nk))
            gd%kfrac = 0.0_dp
            do ik = 1, inp%klist_n
                do k = 1, inp%klist_navg
                    ip = (ik - 1) * inp%klist_navg + k
                    gd%kfrac(inp%klist_axis, ip) = inp%klist_k(ik)
                    if (inp%klist_avg_axis > 0) gd%kfrac(inp%klist_avg_axis, ip) = &
                        -0.5_dp + real(k - 1, dp) / real(inp%klist_navg, dp)
                end do
            end do
            gd%n1 = inp%klist_n
            gd%n2 = 1
            allocate(gd%vertex_k(0), gd%label(0), gd%kdist(nk))
            gd%kdist = 0.0_dp
        else
            call build_kgrid(inp%grid_origin, inp%grid_span1, inp%grid_span2, inp%grid_n1, inp%grid_n2, &
                             recip, kg, ierr, msg)
            if (ierr /= 0) then
                errmsg = 'k-grid: ' // trim(msg); return
            end if
            gd%kmode = GREEN_KGRID
            nk = kg%nk
            gd%kfrac = kg%kfrac
            gd%n1 = kg%n1
            gd%n2 = kg%n2
            allocate(gd%vertex_k(0), gd%label(0), gd%kdist(nk))
            gd%kdist = 0.0_dp
        end if
        gd%nk = nk
        allocate(gd%kcart(3, nk))
        do ik = 1, nk
            gd%kcart(:, ik) = matmul(gd%kfrac(:, ik), recip)
        end do

        ! --- energias e saídas ---
        ne = inp%ne
        gd%ne = ne
        allocate(gd%energy(ne))
        do ie = 1, ne
            if (ne == 1) then
                gd%energy(ie) = inp%emin
            else
                gd%energy(ie) = inp%emin + real(ie - 1, dp) * (inp%emax - inp%emin) / real(ne - 1, dp)
            end if
        end do
        gd%nq = 0
        allocate(gd%qz(0))
        gd%nreg = nreg
        allocate(gd%region(nreg))
        if (len_trim(inp%regions_file) > 0) then
            do g = 1, nreg
                gd%region(g) = rs%regions(g)%name
            end do
        else if (trim(inp%region_mode) == 'atoms') then
            do g = 1, nreg
                gd%region(g) = 'a' // itoa(cent%atom(g)) // '_c' // itoa(cent%x(g))
            end do
        else
            do g = 1, nreg
                gd%region(g) = 'col' // itoa(clo + g - 1)
            end do
        end if
        allocate(gd%total(ne, nk), gd%amap(ne, nk, nreg))
        allocate(gd%niter(ne, nk), gd%residual(ne, nk), gd%min_neg_im_gii(ne, nk), gd%smin(nk))
        gd%smin = 0.0_dp
        gd%has_coherent = .false.
        allocate(gd%icoh(0, 0, 0), gd%iinc(0, 0, 0))

        if (inp%have_rho) then
            allocate(kwt(nk))
            if (inp%kmode == 'list') then
                do ik = 1, nk
                    kwt(ik) = inp%klist_w((ik - 1) / inp%klist_navg + 1) / sum(inp%klist_w(1:inp%klist_n)) / &
                              real(inp%klist_navg, dp)
                end do
            else
                kwt = 1.0_dp / real(nk, dp)
            end if
            allocate(rbuf(nmo, nmo, nrp, inp%rho_nwin, nthr), qbuf(cent%norb, inp%rho_nwin, nthr))
            rbuf = (0.0_dp, 0.0_dp)
            qbuf = (0.0_dp, 0.0_dp)
            write(line, '(A,I0,A,F10.1,A)') 'rho: ', nrp, ' pairs, buffers ', &
                16.0_dp * real(size(rbuf, kind=int64), dp) / 1048576.0_dp, ' MiB'
            call note(trim(line))
        end if

        ! --- laço principal sobre (k, E) ---
        perr_ip = huge(0)
        perr_msg = ''
        nbad = 0
        !$omp parallel do default(shared) schedule(dynamic) private(ip, ik, ie) reduction(+:nbad)
        do ip = 1, nk * ne
            ik = (ip - 1) / ne + 1
            ie = ip - (ik - 1) * ne
            call one_point(ik, ie, nbad)
        end do
        !$omp end parallel do
        if (perr_ip /= huge(0)) then
            ierr = 3
            errmsg = perr_msg
            return
        end if
        if (inp%have_rho .and. inp%rho_contour) then
            ! janelas por contorno: semicírculo de lo a hi no semiplano superior, Gauss-Legendre em theta
            allocate(glx(inp%rho_ncontour), glw(inp%rho_ncontour))
            call gauss_legendre(inp%rho_ncontour, glx, glw)
            !$omp parallel do default(shared) schedule(dynamic) private(ip, ik, ie, w2c, icn, zc, wz, thc) &
            !$omp reduction(+:nbad)
            do ip = 1, nk * inp%rho_nwin * inp%rho_ncontour
                ik = (ip - 1) / (inp%rho_nwin * inp%rho_ncontour) + 1
                ie = ip - (ik - 1) * inp%rho_nwin * inp%rho_ncontour
                w2c = (ie - 1) / inp%rho_ncontour + 1
                icn = ie - (w2c - 1) * inp%rho_ncontour
                thc = 0.5_dp * PI * (1.0_dp + glx(icn))
                zc = cmplx(0.5_dp * (inp%rho_lo(w2c) + inp%rho_hi(w2c)) + eorigin, inp%rho_contour_eta, dp) + &
                     0.5_dp * (inp%rho_hi(w2c) - inp%rho_lo(w2c)) * exp(cmplx(0.0_dp, thc, dp))
                ! dz = -i R e^{i theta} dtheta (theta de pi a 0)
                wz = -cmplx(0.0_dp, 1.0_dp, dp) * 0.5_dp * (inp%rho_hi(w2c) - inp%rho_lo(w2c)) * &
                     exp(cmplx(0.0_dp, thc, dp)) * 0.5_dp * PI * glw(icn)
                call one_point(ik, 1, nbad, zc, wz, w2c)
            end do
            !$omp end parallel do
            if (perr_ip /= huge(0)) then
                ierr = 3
                errmsg = perr_msg
                return
            end if
        end if
        if (present(nonconverged)) nonconverged = nbad
        if (inp%have_rho) then
            call write_rho(ierr, errmsg)
            if (ierr /= 0) return
        end if
        if (inp%kmode == 'list') call reduce_klist()
        call fill_meta()
        if (present(notes)) notes = nt

    contains

        subroutine note(s)
            character(len=*), intent(in) :: s
            nt = [character(len=256) :: nt, s]
        end subroutine note

        integer function block_of(x)
            integer, intent(in) :: x
            integer :: bb, x0
            x0 = clo
            block_of = 0
            do bb = 1, nb
                if (x >= x0 .and. x < x0 + bwid(bb)) then
                    block_of = bb
                    return
                end if
                x0 = x0 + bwid(bb)
            end do
        end function block_of

        function default_fit_columns(il) result(c)
            integer, intent(in) :: il
            integer, allocatable :: c(:)
            integer :: q
            allocate(c(0))
            do side = 1, 2
                if (.not. use_side(side) .or. lmat(side) /= il) cycle
                if (side == 1) then
                    c = [c, [(clo - leads(il)%P + q, q = 0, leads(il)%P - 1)]]
                else
                    c = [c, [(chi + 1 + q, q = 0, leads(il)%P - 1)]]
                end if
            end do
        end function default_fit_columns

        !> Todo acoplamento de um átomo central vai para a própria região (blocos
        !> vizinhos), para a PL de costura do mesmo lado (só no bloco da ponta) ou,
        !> num lado 'cut', para além do corte. Qualquer outro é erro.
        subroutine check_central_couplings(lerr, lmsg)
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            integer, allocatable :: cidx(:, :), sidx(:, :)
            integer :: xmin, xmax, kk, p, xj, j, tgt, sd
            lerr = 0
            lmsg = ''
            xmin = clo - maxval(leads(:)%P) - sup%range - 1
            xmax = chi + maxval(leads(:)%P) + sup%range + 1
            allocate(cidx(sup%nat, xmin:xmax), sidx(sup%nat, xmin:xmax))
            cidx = 0
            sidx = 0
            do kk = 1, cent%n
                cidx(cent%atom(kk), cent%x(kk)) = kk
            end do
            do sd = 1, 2
                if (.not. use_side(sd)) cycle
                do kk = 1, seam(sd)%n
                    sidx(seam(sd)%atom(kk), seam(sd)%x(kk)) = sd
                end do
            end do
            do kk = 1, cent%n
                do p = sup%pstart(cent%atom(kk)), sup%pstart(cent%atom(kk) + 1) - 1
                    xj = cent%x(kk) + sup%pdx(p)
                    j = sup%pj(p)
                    if (xj >= clo .and. xj <= chi) then
                        tgt = cidx(j, xj)
                        if (abs(cblock(tgt) - cblock(kk)) > 1) then
                            lerr = 2; lmsg = 'central coupling skips an RGF block (internal error)'; return
                        end if
                        cycle
                    end if
                    if (sidx(j, xj) == 1 .and. xj < clo) then
                        if (cblock(kk) /= 1) then
                            lerr = 2; lmsg = 'left seam couples beyond the first RGF block'; return
                        end if
                        cycle
                    end if
                    if (sidx(j, xj) == 2 .and. xj > chi) then
                        if (cblock(kk) /= nb) then
                            lerr = 2; lmsg = 'right seam couples beyond the last RGF block'; return
                        end if
                        cycle
                    end if
                    if (xj < clo .and. trim(inp%left) == 'cut') cycle
                    if (xj > chi .and. trim(inp%right) == 'cut') cycle
                    lerr = 2
                    lmsg = 'central atom ' // itoa(cent%atom(kk)) // ' (column ' // itoa(cent%x(kk)) // &
                           ') couples to supercell atom ' // itoa(j) // ' in column ' // itoa(xj) // &
                           ', which is neither central nor in a lead seam'
                    if (xj < clo) then
                        lmsg = lmsg // " (left = '" // trim(inp%left) // "')"
                    else
                        lmsg = lmsg // " (right = '" // trim(inp%right) // "')"
                    end if
                    return
                end do
            end do
        end subroutine check_central_couplings

        subroutine write_xcentral(fname, lerr, lmsg)
            character(len=*), intent(in) :: fname
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            integer :: u, ios, i, s, ns
            character(len=LABEL_LEN), allocatable :: names(:)
            integer, allocatable :: norb(:)
            real(dp) :: lat(3, 3), r(3), ue(3)
            logical :: found
            lerr = 0
            lmsg = ''
            allocate(names(cent%n), norb(cent%n))
            ns = 0
            do i = 1, cent%n
                found = .false.
                do s = 1, ns
                    if (names(s) == sup%label(cent%atom(i))) found = .true.
                end do
                if (.not. found) then
                    ns = ns + 1
                    names(ns) = sup%label(cent%atom(i))
                    norb(ns) = sup%orbs(cent%atom(i))
                end if
            end do
            lat = sup%basis
            ue = sup%basis(e, :) / norm2(sup%basis(e, :))
            lat(e, :) = real(ncol, dp) * sup%basis(e, :) + inp%vacuum_ang * ue
            open(newunit=u, file=fname, status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                lerr = 1; lmsg = 'cannot write central structure: ' // fname; return
            end if
            write(u, '(A)') '! Central region written by surfgreen (mode extracted): columns ' // itoa(clo) // &
                            '..' // itoa(chi) // ' of ' // trim(inp%structure_file)
            write(u, '(A)') '&STRUCTURE'
            do i = 1, 3
                write(u, '(A,I0,A,3(F16.8,:,","))') '    LatticeVectors(', i, ',:) = ', lat(i, :)
            end do
            write(u, '(A)') '/'
            write(u, '(A)') ''
            write(u, '(A)') '&SPECIES'
            write(u, '(I0)') ns
            do s = 1, ns
                write(u, '(A,1X,I0)') trim(names(s)), norb(s)
            end do
            write(u, '(A)') '/'
            write(u, '(A)') ''
            write(u, '(A)') '&ATOMS'
            write(u, '(I0)') cent%n
            do i = 1, cent%n
                r = ext_pos(sup, cent%atom(i), cent%x(i))
                write(u, '(A,3(1X,F16.8),A,I0,A,I0)') sup%label(cent%atom(i)), r, '   ! column ', cent%x(i), &
                    ', supercell atom ', cent%atom(i)
            end do
            write(u, '(A)', iostat=ios) '/'
            close(u)
            if (ios /= 0) then
                lerr = 2; lmsg = 'error writing central structure: ' // fname
            end if
        end subroutine write_xcentral

        subroutine record_error(ipt, msg_in)
            integer, intent(in) :: ipt
            character(len=*), intent(in) :: msg_in
            !$omp critical (xgreen_error)
            if (ipt < perr_ip) then
                perr_ip = ipt
                perr_msg = msg_in
            end if
            !$omp end critical (xgreen_error)
        end subroutine record_error

        subroutine one_point(ik, ie, nbad, zc, wz, iwin)
            integer, intent(in) :: ik, ie
            integer, intent(inout) :: nbad
            complex(dp), intent(in), optional :: zc, wz     ! contorno: z complexo e peso dz
            integer, intent(in), optional :: iwin           ! contorno: janela de &RHO
            complex(dp), allocatable :: h00(:, :), s00(:, :), h01(:, :), s01(:, :), k00(:, :), k01(:, :), k10(:, :)
            complex(dp), allocatable :: gl(:, :, :), hd(:, :), sd(:, :), hu(:, :), su(:, :), hlo(:, :), slo(:, :)
            type(central_t) :: sys
            type(rgf_t) :: res
            real(dp), allocatable :: d(:), rd(:)
            real(dp) :: kc(3), rr(2), smin, dl
            integer :: lerr, it(2), il, n, sd_, bb, iti
            character(len=:), allocatable :: lmsg
            complex(dp) :: z, zl
            logical :: done(2), anti
            complex(dp), allocatable :: ka1(:, :), kbn(:, :)

            anti = trim(inp%sitepot_mode) == 'anticommutator'

            kc = gd%kfrac(:, ik)
            z = cmplx(gd%energy(ie) + eorigin, inp%eta_ev, dp)
            if (present(zc)) z = zc
            it = 0
            rr = 0.0_dp
            done = .false.
            n = maxval(leads(1:nlead)%norb)
            allocate(gl(n, n, 2))
            do sd_ = 1, 2
                if (.not. use_side(sd_) .or. done(sd_)) cycle
                il = lmat(sd_)
                n = leads(il)%norb
                dl = leads(il)%delta
                allocate(h00(n, n), s00(n, n), h01(n, n), s01(n, n))
                call fill_plan(leads(il)%sys, leads(il)%p00, kc, h00, s00)
                call fill_plan(leads(il)%sys, leads(il)%p01, kc, h01, s01)
                zl = z - dl
                k00 = zl * s00 - h00
                k01 = zl * s01 - h01
                k10 = zl * conjg(transpose(s01)) - conjg(transpose(h01))
                if (pair_decim) then
                    call surface_gf_pair(k00, k01, k10, inp%decim_tol, inp%decim_maxit, gl(1:n, 1:n, 1), &
                                         gl(1:n, 1:n, 2), iti, rr(1), rr(2), lerr, lmsg)
                    it = iti
                    done = .true.
                else if (sd_ == 1) then
                    call surface_gf(k00, k01, k10, SIDE_TOP, inp%decim_tol, inp%decim_maxit, gl(1:n, 1:n, 1), &
                                    iti, rr(1), lerr, lmsg)
                    it(1) = iti
                    done(1) = .true.
                else
                    call surface_gf(k00, k01, k10, SIDE_BOTTOM, inp%decim_tol, inp%decim_maxit, gl(1:n, 1:n, 2), &
                                    iti, rr(2), lerr, lmsg)
                    it(2) = iti
                    done(2) = .true.
                end if
                if (lerr == 3) then
                    nbad = nbad + 1
                    it = -abs(it)
                else if (lerr /= 0) then
                    call record_error(ip_of(ik, ie), lmsg); return
                end if
                deallocate(h00, s00, h01, s01)
            end do

            call central_init(sys, bsize)
            do bb = 1, nb
                allocate(hd(bsize(bb), bsize(bb)), sd(bsize(bb), bsize(bb)))
                call fill_plan(sup, pd(bb), kc, hd, sd)
                if (have_v) call vpot_diag(vb(bb), hd, sd)
                if (have_dip) call dip_fill(pd(bb), dpd(bb), kc, hd)
                if (bb < nb) then
                    allocate(hu(bsize(bb), bsize(bb + 1)), su(bsize(bb), bsize(bb + 1)))
                    allocate(hlo(bsize(bb + 1), bsize(bb)), slo(bsize(bb + 1), bsize(bb)))
                    call fill_plan(sup, pu(bb), kc, hu, su)
                    call fill_plan(sup, pl(bb), kc, hlo, slo)
                    if (have_v .and. anti) then
                        call vpot_cross(hu, su, vb(bb)%v, vb(bb + 1)%v)
                        call vpot_cross(hlo, slo, vb(bb + 1)%v, vb(bb)%v)
                    end if
                    if (have_dip) then
                        call dip_fill(pu(bb), dpu(bb), kc, hu)
                        call dip_fill(pl(bb), dpl(bb), kc, hlo)
                    end if
                    call assemble_central(sys, bb, z, hd, sd, 0.0_dp, hu, su, 0.0_dp, hl=hlo, sl=slo)
                    deallocate(hu, su, hlo, slo)
                else
                    call assemble_central(sys, bb, z, hd, sd, 0.0_dp)
                end if
                deallocate(hd, sd)
            end do
            if (use_side(1)) then
                n = leads(lmat(1))%norb
                allocate(hu(n, bsize(1)), su(n, bsize(1)))
                call fill_plan(sup, ps(1), kc, hu, su)
                if (have_v .and. anti) call vpot_cross(hu, su, [(inp%sitepot_lead_ev, iti = 1, n)], vb(1)%v)
                if (have_dip) call dip_fill(ps(1), dps(1), kc, hu)
                if (inp%have_rho) ka1 = z * su - hu
                call attach_lead(sys, 'bot', z, hu, su, 0.0_dp, 0.0_dp, gl(1:n, 1:n, 1))
                deallocate(hu, su)
            end if
            if (use_side(2)) then
                n = leads(lmat(2))%norb
                allocate(hu(bsize(nb), n), su(bsize(nb), n))
                call fill_plan(sup, ps(2), kc, hu, su)
                if (have_v .and. anti) call vpot_cross(hu, su, vb(nb)%v, [(inp%sitepot_lead_ev, iti = 1, n)])
                if (have_dip) call dip_fill(ps(2), dps(2), kc, hu)
                if (inp%have_rho) kbn = z * conjg(transpose(su)) - conjg(transpose(hu))
                call attach_lead(sys, 'top', z, hu, su, 0.0_dp, 0.0_dp, gl(1:n, 1:n, 2))
                deallocate(hu, su)
            end if
            call rgf_solve(sys, res, lerr, lmsg)
            if (lerr /= 0) then
                call record_error(ip_of(ik, ie), lmsg); return
            end if
            if (present(zc)) then
                block
                    logical :: inwc(16)
                    inwc = .false.
                    inwc(iwin) = .true.
                    call rho_core(ik, wz, inwc, res, ka1, kbn, gl)
                end block
                return
            end if
            allocate(d(cent%norb), rd(nreg))
            if (trim(inp%projection) == 'lowdin') then
                call orbital_dos_lowdin(sys, d, smin, lerr, lmsg)
                if (lerr /= 0) then
                    call record_error(ip_of(ik, ie), 'Lowdin projection: ' // lmsg); return
                end if
            else
                call orbital_dos_mulliken(sys, res, d, lead_cross=inp%lead_cross)
            end if
            call region_dos(d, groups, rd)
            if (inp%have_rho .and. .not. inp%rho_contour) call accumulate_rho(ik, ie, res, d, ka1, kbn, gl)
            gd%total(ie, ik) = sum(d)
            gd%amap(ie, ik, :) = rd
            if (any(it < 0)) then
                gd%niter(ie, ik) = minval(it)
            else
                gd%niter(ie, ik) = maxval(it)
            end if
            gd%residual(ie, ik) = maxval(rr)
            gd%min_neg_im_gii(ie, ik) = min_neg_im_gdiag(res)
            if (ie == 1 .and. inp%overlap_diagnostic) then
                call central_overlap_min(sys, smin, lerr, lmsg)
                gd%smin(ik) = smin
            end if
        end subroutine one_point

        !> Dipolos (A6): garante as integrais S, X de todos os pares dos planos (e de &RHO) no cache
        !> e, com scc_f, monta os blocos reais do termo de dipolo por entrada de plano.
        subroutine setup_dipole(lerr, lmsg)
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            real(dp), allocatable :: dl(:, :)
            integer :: kind, bb, e2, cnt, p, no
            real(dp) :: smax
            character(len=256) :: ln
            lerr = 0
            lmsg = ''
            if (.not. (present(scc_pao) .and. present(scc_dip))) then
                lerr = 2; lmsg = 'dipoles: scc_pao and scc_dip are required'; return
            end if
            if (present(rho_p) .and. .not. inp%have_rho) then
                lerr = 2; lmsg = 'dipoles: rho_p needs &RHO'; return
            end if
            if (have_dip) then
                if (size(scc_f, 1) /= 3 .or. size(scc_f, 2) /= cent%n) then
                    lerr = 2; lmsg = 'scc_f: shape differs from (3, central atoms)'; return
                end if
            end if
            no = scc_pao%norb
            cnt = 0
            do kind = 1, 5
                do bb = 1, nb
                    if (.not. dplan_used(kind, bb)) cycle
                    associate(pp => plan_of(kind, bb))
                        do e2 = 1, pp%nent
                            p = pp%pair(e2)
                            if (sup%orbs(sup%pi(p)) /= no .or. sup%orbs(sup%pj(p)) /= no) then
                                lerr = 2
                                lmsg = 'dipoles: atom with ' // itoa(sup%orbs(sup%pi(p))) // ' orbitals; the PAO ' // &
                                       'basis has ' // itoa(no) // ' (one species only)'
                                return
                            end if
                            cnt = cnt + 1
                        end do
                    end associate
                end do
            end do
            allocate(dl(3, cnt + nrp))
            cnt = 0
            smax = 0.0_dp
            do kind = 1, 5
                do bb = 1, nb
                    if (.not. dplan_used(kind, bb)) cycle
                    associate(pp => plan_of(kind, bb))
                        do e2 = 1, pp%nent
                            cnt = cnt + 1
                            dl(:, cnt) = sup%pvec(:, pp%pair(e2)) / BOHR_ANG_X
                        end do
                    end associate
                end do
            end do
            if (nrp > 0) dl(:, cnt + 1:cnt + nrp) = rp_vec(:, 1:nrp) / BOHR_ANG_X
            call pao_dip_ensure(scc_dip, scc_pao, dl)
            if (have_dip .and. present(scc_vb)) then
                allocate(dpd(nb), dpu(max(nb - 1, 1)), dpl(max(nb - 1, 1)))
                no = 0
                do bb = 1, nb
                    call make_bplan(1, bb, pd(bb), dpd(bb), no)
                    if (bb < nb) then
                        call make_bplan(2, bb, pu(bb), dpu(bb), no)
                        call make_bplan(3, bb, pl(bb), dpl(bb), no)
                    end if
                end do
                if (use_side(1)) call make_bplan(4, 1, ps(1), dps(1), no)
                if (use_side(2)) call make_bplan(5, nb, ps(2), dps(2), no)
                write(ln, '(A,I0,A,I0,A,I0,A)') 'bond charges: ', cent%n, ' atoms + ', nbd, ' bond sites; ', no, &
                    ' coupling blocks with S /= 0 beyond &RHO rmax (no potential there)'
                call note(trim(ln))
            else if (have_dip) then
                allocate(dpd(nb), dpu(max(nb - 1, 1)), dpl(max(nb - 1, 1)))
                do bb = 1, nb
                    call make_dplan(pd(bb), blk(bb), blk(bb), .true., .true., dpd(bb), smax)
                    if (bb < nb) then
                        call make_dplan(pu(bb), blk(bb), blk(bb + 1), .true., .true., dpu(bb), smax)
                        call make_dplan(pl(bb), blk(bb + 1), blk(bb), .true., .true., dpl(bb), smax)
                    end if
                end do
                if (use_side(1)) call make_dplan(ps(1), seam(1), blk(1), .false., .true., dps(1), smax)
                if (use_side(2)) call make_dplan(ps(2), blk(nb), seam(2), .true., .false., dps(2), smax)
                if (present(scc_dhp)) call add_dhp()
                write(ln, '(A,I0,A,F6.3,A,ES10.3)') 'dipoles: ', scc_dip%n, ' two-center vectors (grid h = ', &
                    scc_dip%h, ' Bohr); max|S_grid - S_file| = ', smax
                call note(trim(ln))
            end if
        end subroutine setup_dipole

        !> A10: sítios de ligação. Cada par de &RHO com d /= 0 é uma ordem (i -> j, d) de um sítio; o par
        !> (j -> i, -d) entre dois átomos centrais é o mesmo sítio (peso 1 em cada ordem); o par com a costura
        !> (j fora da região central) entra uma vez com peso 2 (a outra ordem é a transposta, D_ji(-d) = D_ij(d)^T).
        !> Monta rp_bond/rp_bw/rp_jc, bd_q (ordem canônica de cada sítio) e qlook (entrada de plano -> par).
        subroutine setup_bonds(lerr, lmsg)
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            integer :: q, q2, bb, nmx, nsolo, ib, ic
            integer, allocatable :: first(:), nxt(:)
            lerr = 0
            lmsg = ''
            if (inp%rho_rmax <= 0.0_dp) then
                lerr = 2; lmsg = 'bond charges need rho_rmax > 0 (pair range)'; return
            end if
            if ((present(rho_qb) .or. present(rho_qon) .or. present(rho_pon)) .and. .not. inp%have_rho) then
                lerr = 2; lmsg = 'bond charges: rho_qb, rho_qon, rho_pon need &RHO'; return
            end if
            if (present(scc_vb)) then
                if (present(scc_v) .or. .not. present(scc_f)) then
                    lerr = 2; lmsg = 'scc_vb: needs scc_f and excludes scc_v'; return
                end if
            end if
            allocate(rp_bond(nrp), rp_bw(nrp), rp_jc(nrp), bd_q(nrp), first(cent%n), nxt(nrp))
            do q = 1, nrp
                bb = rp_b(q)
                select case (rp_kind(q))
                case (1)
                    rp_jc(q) = cent_of_offset(bb, pd(bb)%co(rp_e(q)))
                case (2)
                    rp_jc(q) = cent_of_offset(bb + 1, pu(bb)%co(rp_e(q)))
                case (3)
                    rp_jc(q) = cent_of_offset(bb, pl(bb)%co(rp_e(q)))
                case default
                    rp_jc(q) = -1
                end select
            end do
            first = 0
            nxt = 0
            do q = nrp, 1, -1
                if (rp_ic(q) < 1) then
                    lerr = 2; lmsg = 'bond charges: a &RHO pair is not mapped to a central atom'; return
                end if
                nxt(q) = first(rp_ic(q))
                first(rp_ic(q)) = q
            end do
            rp_bond = 0
            rp_bw = 0
            nbd = 0
            nsolo = 0
            do q = 1, nrp
                if (norm2(rp_vec(:, q)) < 1.0e-8_dp) then
                    if (rp_jc(q) /= rp_ic(q)) then
                        lerr = 2; lmsg = 'bond charges: d = 0 pair between different atoms'; return
                    end if
                    cycle
                end if
                if (rp_bond(q) > 0) cycle
                nbd = nbd + 1
                bd_q(nbd) = q
                rp_bond(q) = nbd
                rp_bw(q) = 2
                if (rp_jc(q) < 1) cycle
                q2 = first(rp_jc(q))
                do while (q2 > 0)
                    if (q2 /= q .and. rp_jc(q2) == rp_ic(q) .and. rp_bond(q2) == 0) then
                        if (norm2(rp_vec(:, q2) + rp_vec(:, q)) < 1.0e-6_dp) exit
                    end if
                    q2 = nxt(q2)
                end do
                if (q2 > 0) then
                    rp_bond(q2) = nbd
                    rp_bw(q) = 1
                    rp_bw(q2) = 1
                else
                    nsolo = nsolo + 1
                end if
            end do
            nmx = 1
            do bb = 1, nb
                nmx = max(nmx, pd(bb)%nent)
                if (bb < nb) nmx = max(nmx, pu(bb)%nent, pl(bb)%nent)
            end do
            if (use_side(1)) nmx = max(nmx, ps(1)%nent)
            if (use_side(2)) nmx = max(nmx, ps(2)%nent)
            allocate(qlook(nmx, 5, nb))
            qlook = 0
            do q = 1, nrp
                qlook(rp_e(q), rp_kind(q), rp_b(q)) = q
            end do
            call note('bond sites: ' // itoa(nbd) // ' (' // itoa(count(rp_bw == 2 .and. rp_jc < 1)) // &
                      ' with the seam, ' // itoa(nsolo) // ' central pairs without the reverse order)')
            if (present(bond_ij)) then
                allocate(bond_ij(4, nbd))
                do ib = 1, nbd
                    q = bd_q(ib)
                    bond_ij(:, ib) = [rp_ic(q), rp_jc(q), cent%atom(rp_ic(q)), rp_ja(q)]
                end do
            end if
            if (present(bond_d)) bond_d = rp_vec(:, bd_q(1:nbd))
            if (present(bond_pos)) then
                allocate(bond_pos(3, nbd))
                do ib = 1, nbd
                    q = bd_q(ib)
                    ic = rp_ic(q)
                    bond_pos(:, ib) = ext_pos(sup, cent%atom(ic), cent%x(ic)) + 0.5_dp * rp_vec(:, q)
                end do
            end if
            if (present(scc_vb)) then
                if (size(scc_vb) /= cent%n + nbd) then
                    lerr = 2; lmsg = 'scc_vb: size differs from central atoms + bond sites'; return
                end if
            end if
        end subroutine setup_bonds

        !> A10: termo de H por entrada de plano no modo de ligações (S do arquivo):
        !> on-site (i, i, d = 0): V_i S + F_i . X^(i); par de &RHO (sítio b): V_b S; fora de &RHO: nada.
        subroutine make_bplan(kind, bb, pp, dq, nmiss)
            integer, intent(in) :: kind, bb
            type(plan_t), intent(in) :: pp
            type(dplan_t), intent(out) :: dq
            integer, intent(inout) :: nmiss
            real(dp), allocatable :: sb(:, :), xb(:, :, :), dh(:, :), sf(:, :)
            integer :: e2, p, ni, nj, tot, q, ic
            logical :: okd
            allocate(dq%off(pp%nent))
            tot = 0
            do e2 = 1, pp%nent
                p = pp%pair(e2)
                dq%off(e2) = tot
                tot = tot + sup%orbs(sup%pi(p)) * sup%orbs(sup%pj(p))
            end do
            allocate(dq%dat(tot), sb(nmo_d(), nmo_d()), xb(nmo_d(), nmo_d(), 3), dh(nmo_d(), nmo_d()), &
                     sf(nmo_d(), nmo_d()))
            do e2 = 1, pp%nent
                p = pp%pair(e2)
                ni = sup%orbs(sup%pi(p)); nj = sup%orbs(sup%pj(p))
                sf(1:ni, 1:nj) = 0.0_dp
                if (sup%soff(p) >= 0) sf(1:ni, 1:nj) = real(reshape(sup%sdat(sup%soff(p) + 1:sup%soff(p) + &
                    int(ni * nj, int64)), [ni, nj]), dp)
                q = qlook(e2, kind, bb)
                if (q == 0) then
                    dh(1:ni, 1:nj) = 0.0_dp
                    if (any(sf(1:ni, 1:nj) /= 0.0_dp)) nmiss = nmiss + 1
                else if (rp_bond(q) == 0) then
                    ic = rp_ic(q)
                    call pao_dip_get(scc_dip, [0.0_dp, 0.0_dp, 0.0_dp], sb(1:ni, 1:nj), xb(1:ni, 1:nj, :), okd)
                    if (.not. okd) xb = 0.0_dp
                    if (.not. okd) nmiss = nmiss + 1
                    call pao_dipole_dh(scc_f(:, ic), scc_f(:, ic), xb(1:ni, 1:nj, :), sb(1:ni, 1:nj), &
                                       [0.0_dp, 0.0_dp, 0.0_dp], dh(1:ni, 1:nj))
                    dh(1:ni, 1:nj) = dh(1:ni, 1:nj) + scc_vb(ic) * sf(1:ni, 1:nj)
                else
                    dh(1:ni, 1:nj) = scc_vb(cent%n + rp_bond(q)) * sf(1:ni, 1:nj)
                end if
                dq%dat(dq%off(e2) + 1:dq%off(e2) + ni * nj) = reshape(dh(1:ni, 1:nj), [ni * nj])
            end do
        end subroutine make_bplan

        logical function dplan_used(kind, bb) result(yes)
            integer, intent(in) :: kind, bb
            yes = .true.
            if ((kind == 2 .or. kind == 3) .and. bb == nb) yes = .false.
            if (kind == 4 .and. (bb /= 1 .or. .not. use_side(1))) yes = .false.
            if (kind == 5 .and. (bb /= nb .or. .not. use_side(2))) yes = .false.
        end function dplan_used

        subroutine make_dplan(pp, rl, cl, rcen, ccen, dq, smax)
            type(plan_t), intent(in) :: pp
            type(xlist_t), intent(in) :: rl, cl
            logical, intent(in) :: rcen, ccen
            type(dplan_t), intent(out) :: dq
            real(dp), intent(inout) :: smax
            real(dp), allocatable :: sb(:, :), xb(:, :, :), dh(:, :)
            real(dp) :: fi(3), fj(3), dv(3)
            integer :: e2, p, ni, nj, tot, ic
            logical :: okd
            allocate(dq%off(pp%nent))
            tot = 0
            do e2 = 1, pp%nent
                p = pp%pair(e2)
                dq%off(e2) = tot
                tot = tot + sup%orbs(sup%pi(p)) * sup%orbs(sup%pj(p))
            end do
            allocate(dq%dat(tot), sb(nmo_d(), nmo_d()), xb(nmo_d(), nmo_d(), 3), dh(nmo_d(), nmo_d()))
            do e2 = 1, pp%nent
                p = pp%pair(e2)
                ni = sup%orbs(sup%pi(p)); nj = sup%orbs(sup%pj(p))
                dv = sup%pvec(:, p) / BOHR_ANG_X
                fi = 0.0_dp; fj = 0.0_dp
                if (rcen) then
                    ic = cent_index(rl%atom(pp%ri(e2)), rl%x(pp%ri(e2)))
                    if (ic > 0) fi = scc_f(:, ic)
                end if
                if (ccen) then
                    ic = cent_index(cl%atom(pp%ci(e2)), cl%x(pp%ci(e2)))
                    if (ic > 0) fj = scc_f(:, ic)
                end if
                call pao_dip_get(scc_dip, dv, sb(1:ni, 1:nj), xb(1:ni, 1:nj, :), okd)
                if (.not. okd) then
                    dq%dat(dq%off(e2) + 1:dq%off(e2) + ni * nj) = 0.0_dp
                    smax = huge(1.0_dp)
                    cycle
                end if
                if (sup%soff(p) >= 0) smax = max(smax, maxval(abs(sb(1:ni, 1:nj) - real(reshape(sup%sdat( &
                    sup%soff(p) + 1:sup%soff(p) + int(ni * nj, int64)), [ni, nj]), dp))))
                call pao_dipole_dh(fi, fj, xb(1:ni, 1:nj, :), sb(1:ni, 1:nj), dv, dh(1:ni, 1:nj))
                dq%dat(dq%off(e2) + 1:dq%off(e2) + ni * nj) = reshape(dh(1:ni, 1:nj), [ni * nj])
            end do
        end subroutine make_dplan

        !> A19: soma scc_dhp (bloco real por par de &RHO, linhas no átomo central) às entradas de plano.
        subroutine add_dhp()
            integer :: q
            do q = 1, nrp
                select case (rp_kind(q))
                case (1)
                    call add_dhp_block(dpd(rp_b(q))%dat, dpd(rp_b(q))%off(rp_e(q)), q, .false.)
                case (2)
                    call add_dhp_block(dpu(rp_b(q))%dat, dpu(rp_b(q))%off(rp_e(q)), q, .false.)
                case (3)
                    call add_dhp_block(dpl(rp_b(q))%dat, dpl(rp_b(q))%off(rp_e(q)), q, .false.)
                case (4)
                    call add_dhp_block(dps(1)%dat, dps(1)%off(rp_e(q)), q, .true.)
                case default
                    call add_dhp_block(dps(2)%dat, dps(2)%off(rp_e(q)), q, .false.)
                end select
            end do
        end subroutine add_dhp

        !> tr: bloco do plano (costura, central) = dh^T (linhas: os nj orbitais da costura).
        subroutine add_dhp_block(dat, o, q, tr)
            real(dp), intent(inout) :: dat(:)
            integer, intent(in) :: o, q
            logical, intent(in) :: tr
            integer :: a, b2, ni, nj
            ni = rp_ni(q); nj = rp_nj(q)
            if (tr) then
                do a = 1, ni
                    do b2 = 1, nj
                        dat(o + (a - 1) * nj + b2) = dat(o + (a - 1) * nj + b2) + scc_dhp(a, b2, q)
                    end do
                end do
            else
                do b2 = 1, nj
                    do a = 1, ni
                        dat(o + (b2 - 1) * ni + a) = dat(o + (b2 - 1) * ni + a) + scc_dhp(a, b2, q)
                    end do
                end do
            end if
        end subroutine add_dhp_block

        integer function nmo_d()
            nmo_d = maxval(sup%orbs)
        end function nmo_d

        !> h += sum_e e^{i k.d} dH_e (mesma fase e leiaute de fill_plan).
        subroutine dip_fill(pp, dq, kf, h)
            type(plan_t), intent(in) :: pp
            type(dplan_t), intent(in) :: dq
            real(dp), intent(in) :: kf(3)
            complex(dp), intent(inout) :: h(:, :)
            integer :: e2, p, ni, nj, a, b2
            complex(dp) :: ph
            real(dp) :: kcart(3)
            kcart = matmul(kf, sup%recip)
            do e2 = 1, pp%nent
                p = pp%pair(e2)
                ni = sup%orbs(sup%pi(p)); nj = sup%orbs(sup%pj(p))
                ph = exp(cmplx(0.0_dp, dot_product(kcart, sup%pvec(:, p)), dp))
                do b2 = 1, nj
                    do a = 1, ni
                        h(pp%ro(e2) + a, pp%co(e2) + b2) = h(pp%ro(e2) + a, pp%co(e2) + b2) + &
                            ph * dq%dat(dq%off(e2) + (b2 - 1) * ni + a)
                    end do
                end do
            end do
        end subroutine dip_fill

        !> Monta vb(:) a partir de &SITEPOT: casamento por posição (módulo os
        !> vetores de rede dos eixos de k) e, opcionalmente, pela coluna X.
        subroutine setup_sitepot(lerr, lmsg)
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            type(site_potential_t) :: sp
            real(dp), allocatable :: trans(:, :)
            real(dp) :: pos(3), vs
            integer :: bb, kk, j, no, o, idx, nmatch, nblk, a
            logical :: have_file
            lerr = 0
            lmsg = ''
            have_file = len_trim(inp%sitepot_file) > 0
            if (have_file) then
                call read_site_potential(trim(inp%sitepot_file), sp, lerr, lmsg)
                if (lerr /= 0) return
            end if
            allocate(trans(3, 0))
            do a = 1, 3
                if (kaxis(a)) trans = reshape([trans, sup%lattice(a, :)], [3, size(trans, 2) + 1])
            end do
            allocate(vb(nb))
            nmatch = 0
            nblk = 0
            do bb = 1, nb
                allocate(vb(bb)%v(bsize(bb)), vb(bb)%av(blk(bb)%n), vb(bb)%aoff(blk(bb)%n), vb(bb)%aorb(blk(bb)%n))
                vb(bb)%v = 0.0_dp
                do kk = 1, blk(bb)%n
                    j = blk(bb)%atom(kk)
                    no = sup%orbs(j)
                    o = blk(bb)%off(kk)
                    vs = inp%sitepot_uniform_ev
                    if (have_file) then
                        pos = ext_pos(sup, j, blk(bb)%x(kk)) + inp%sitepot_shift
                        idx = site_potential_find(sp, blk(bb)%x(kk), pos, trans, inp%sitepot_tol, &
                                                  match_tag=inp%sitepot_use_tag)
                        if (idx > 0) then
                            nmatch = nmatch + 1
                            if (sp%is_scalar(idx)) then
                                vs = vs + sp%v(1, 1, idx)
                            else
                                if (sp%norb /= no) then
                                    lerr = 2
                                    lmsg = '&SITEPOT: block size differs from the orbitals of atom ' // itoa(j)
                                    return
                                end if
                                if (.not. vb(bb)%has_block) then
                                    allocate(vb(bb)%a(bsize(bb), bsize(bb)))
                                    vb(bb)%a = 0.0_dp
                                    vb(bb)%has_block = .true.
                                end if
                                vb(bb)%a(o + 1:o + no, o + 1:o + no) = 0.5_dp * (sp%v(:, :, idx) + &
                                                                                 transpose(sp%v(:, :, idx)))
                                nblk = nblk + 1
                            end if
                        end if
                    end if
                    if (present(scc_v)) vs = vs + scc_v(cent_index(j, blk(bb)%x(kk)))
                    vb(bb)%v(o + 1:o + no) = vs
                    vb(bb)%av(kk) = vs
                    vb(bb)%aoff(kk) = o
                    vb(bb)%aorb(kk) = no
                end do
            end do
            call note('site potential (&SITEPOT, ' // trim(inp%sitepot_mode) // '): ' // itoa(nmatch) // &
                      ' of ' // itoa(cent%n) // ' central atoms matched in ''' // trim(inp%sitepot_file) // &
                      ''' (' // itoa(nblk) // ' orbital blocks)')
        end subroutine setup_sitepot

        !> On-site de um bloco: 'anticommutator' H += 1/2 {S, V} (V escalar por
        !> orbital); 'onsite' H_ii += V_i S_ii por átomo; blocos orbitais somados.
        subroutine vpot_diag(v, h, s)
            type(vblock_t), intent(in) :: v
            complex(dp), intent(inout) :: h(:, :)
            complex(dp), intent(in) :: s(:, :)
            integer :: kk, o, no
            if (trim(inp%sitepot_mode) == 'anticommutator') then
                call vpot_cross(h, s, v%v, v%v)
            else
                do kk = 1, size(v%av)
                    o = v%aoff(kk); no = v%aorb(kk)
                    h(o + 1:o + no, o + 1:o + no) = h(o + 1:o + no, o + 1:o + no) + &
                        cmplx(v%av(kk), 0.0_dp, dp) * s(o + 1:o + no, o + 1:o + no)
                end do
            end if
            if (v%has_block) h = h + cmplx(v%a, 0.0_dp, dp)
        end subroutine vpot_diag

        pure integer function ip_of(ik, ie)
            integer, intent(in) :: ik, ie
            ip_of = (ik - 1) * ne + ie
        end function ip_of

        !> Pares (i, j) de &RHO: i central, j central ou da costura, |r_j + R - r_i| < rmax.
        !> Tipos: 1 = mesmo bloco, 2 = bloco b -> b+1, 3 = b+1 -> b, 4 = central -> costura
        !> esquerda (inverso de ps(1)), 5 = central -> costura direita.
        subroutine setup_rho_pairs()
            integer :: pass, bb, e2, p, kind, ni, nj, cnt
            real(dp) :: v(3)
            nmo = maxval(sup%orbs)
            do pass = 1, 2
                cnt = 0
                do kind = 1, 5
                    do bb = 1, nb
                        if ((kind == 2 .or. kind == 3) .and. bb == nb) cycle
                        if (kind == 4 .and. (bb /= 1 .or. .not. use_side(1))) cycle
                        if (kind == 5 .and. (bb /= nb .or. .not. use_side(2))) cycle
                        associate(pp => plan_of(kind, bb))
                            do e2 = 1, pp%nent
                                p = pp%pair(e2)
                                v = sup%pvec(:, p)
                                if (kind == 4) v = -v
                                if (norm2(v) >= inp%rho_rmax .and. norm2(v) > 1.0e-8_dp) cycle
                                cnt = cnt + 1
                                if (pass == 1) cycle
                                rp_kind(cnt) = kind; rp_b(cnt) = bb; rp_e(cnt) = e2; rp_vec(:, cnt) = v
                                if (kind == 4) then
                                    ni = sup%orbs(sup%pj(p)); nj = sup%orbs(sup%pi(p))
                                    rp_ic(cnt) = cent_of_offset(1, pp%co(e2))
                                    rp_ja(cnt) = sup%pi(p)
                                else
                                    ni = sup%orbs(sup%pi(p)); nj = sup%orbs(sup%pj(p))
                                    rp_ic(cnt) = -1
                                    rp_ja(cnt) = sup%pj(p)
                                end if
                                rp_ni(cnt) = ni; rp_nj(cnt) = nj
                                rp_s(:, :, cnt) = 0.0_dp
                                if (sup%soff(p) >= 0) then
                                    if (kind == 4) then
                                        rp_s(1:ni, 1:nj, cnt) = transpose(real(reshape(sup%sdat(sup%soff(p) + 1: &
                                            sup%soff(p) + int(ni * nj, int64)), [nj, ni]), dp))
                                    else
                                        rp_s(1:ni, 1:nj, cnt) = real(reshape(sup%sdat(sup%soff(p) + 1: &
                                            sup%soff(p) + int(ni * nj, int64)), [ni, nj]), dp)
                                    end if
                                end if
                            end do
                        end associate
                    end do
                end do
                if (pass == 1) then
                    nrp = cnt
                    allocate(rp_kind(nrp), rp_b(nrp), rp_e(nrp), rp_ic(nrp), rp_ni(nrp), rp_nj(nrp), rp_ja(nrp))
                    allocate(rp_jx(nrp), rp_vec(3, nrp), rp_s(nmo, nmo, nrp))
                    rp_jx = 0
                end if
            end do
            ! átomo central de cada par (linha): índice em cent pela posição da linha no bloco
            do cnt = 1, nrp
                if (rp_kind(cnt) == 4) cycle
                associate(pp => plan_of(rp_kind(cnt), rp_b(cnt)))
                    bb = rp_b(cnt)
                    if (rp_kind(cnt) == 3) bb = bb + 1
                    rp_ic(cnt) = cent_of_offset(bb, pp%ro(rp_e(cnt)))
                end associate
            end do
            if (any(rp_ic < 1)) call note('rho: WARNING some pairs could not be mapped to a central atom')
        end subroutine setup_rho_pairs

        function plan_of(kind, bb) result(pp)
            integer, intent(in) :: kind, bb
            type(plan_t) :: pp
            select case (kind)
            case (1)
                pp = pd(bb)
            case (2)
                pp = pu(bb)
            case (3)
                pp = pl(bb)
            case (4)
                pp = ps(1)
            case default
                pp = ps(2)
            end select
        end function plan_of

        !> Índice em cent do átomo central cujo bloco RGF bb começa em ro.
        integer function cent_of_offset(bb, ro) result(ic)
            integer, intent(in) :: bb, ro
            integer :: kk
            ic = -1
            do kk = 1, blk(bb)%n
                if (blk(bb)%off(kk) == ro) then
                    ic = cent_index(blk(bb)%atom(kk), blk(bb)%x(kk))
                    return
                end if
            end do
        end function cent_of_offset

        integer function cent_index(at, x) result(ic)
            integer, intent(in) :: at, x
            integer :: kk
            ic = -1
            do kk = 1, cent%n
                if (cent%atom(kk) == at .and. cent%x(kk) == x) then
                    ic = kk
                    return
                end if
            end do
        end function cent_index

        !> Soma, no ponto (k, E), (i/2pi)(G_ij - G_ji^dag) dE w_k e^{-i k.d_ij} de cada
        !> par e o DOS de Mulliken orbital (d) das janelas de &RHO que contêm E.
        subroutine accumulate_rho(ik, ie, rg, dd, ka1, kbn, glead)
            integer, intent(in) :: ik, ie
            type(rgf_t), intent(in) :: rg
            real(dp), intent(in) :: dd(:)
            complex(dp), allocatable, intent(in) :: ka1(:, :), kbn(:, :)
            complex(dp), intent(in) :: glead(:, :, :)
            real(dp) :: de_, en
            integer :: tid, w2
            logical :: inw(16)
            en = gd%energy(ie)
            inw = .false.
            do w2 = 1, inp%rho_nwin
                inw(w2) = en >= inp%rho_lo(w2) .and. en < inp%rho_hi(w2)
            end do
            if (.not. any(inw)) return
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            de_ = 0.0_dp
            if (ne > 1) de_ = (inp%emax - inp%emin) / real(ne - 1, dp)
            do w2 = 1, inp%rho_nwin
                if (inw(w2)) qbuf(:, w2, tid) = qbuf(:, w2, tid) + cmplx(de_ * kwt(ik) * dd, 0.0_dp, dp)
            end do
            call rho_core(ik, cmplx(de_, 0.0_dp, dp), inw, rg, ka1, kbn, glead)
        end subroutine accumulate_rho

        !> Soma (i/2pi) w_k e^{-i k.d_ij} [wz G_ij(z) - conj(wz) G_ji(z)^dag] de cada par nas janelas inw.
        !> Eixo real: wz = dE, z = E + i eta. Contorno: wz = dz no semicírculo (G(z*) = G(z)^dag).
        subroutine rho_core(ik, wz, inw, rg, ka1, kbn, glead)
            integer, intent(in) :: ik
            complex(dp), intent(in) :: wz
            logical, intent(in) :: inw(16)
            type(rgf_t), intent(in) :: rg
            complex(dp), allocatable, intent(in) :: ka1(:, :), kbn(:, :)
            complex(dp), intent(in) :: glead(:, :, :)
            complex(dp), allocatable :: ga1(:, :), gbn(:, :), tmp(:, :), gij(:, :), gji(:, :)
            complex(dp) :: fac, ph
            real(dp) :: kc(3)
            integer :: tid, w2, q, ni, nj, ro, co, bb, nl
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            kc = matmul(gd%kfrac(:, ik), sup%recip)
            ! G_{A,1} = -g_A K_{A,1} G_11, G_{B,N} = -g_B K_{B,N} G_NN
            if (use_side(1)) then
                nl = size(ka1, 1)
                allocate(tmp(nl, bsize(1)), ga1(nl, bsize(1)))
                tmp = matmul(ka1, rg%gd(1)%a)
                ga1 = -matmul(glead(1:nl, 1:nl, 1), tmp)
                deallocate(tmp)
            end if
            if (use_side(2)) then
                nl = size(kbn, 1)
                allocate(tmp(nl, bsize(nb)), gbn(nl, bsize(nb)))
                tmp = matmul(kbn, rg%gd(nb)%a)
                gbn = -matmul(glead(1:nl, 1:nl, 2), tmp)
                deallocate(tmp)
            end if
            fac = cmplx(0.0_dp, kwt(ik) / (2.0_dp * PI), dp)
            do q = 1, nrp
                ni = rp_ni(q); nj = rp_nj(q); bb = rp_b(q)
                select case (rp_kind(q))
                case (1)
                    ro = pd(bb)%ro(rp_e(q)); co = pd(bb)%co(rp_e(q))
                    gij = rg%gd(bb)%a(ro + 1:ro + ni, co + 1:co + nj)
                    gji = rg%gd(bb)%a(co + 1:co + nj, ro + 1:ro + ni)
                case (2)
                    ro = pu(bb)%ro(rp_e(q)); co = pu(bb)%co(rp_e(q))
                    gij = rg%gu(bb)%a(ro + 1:ro + ni, co + 1:co + nj)
                    gji = rg%gl(bb)%a(co + 1:co + nj, ro + 1:ro + ni)
                case (3)
                    ro = pl(bb)%ro(rp_e(q)); co = pl(bb)%co(rp_e(q))
                    gij = rg%gl(bb)%a(ro + 1:ro + ni, co + 1:co + nj)
                    gji = rg%gu(bb)%a(co + 1:co + nj, ro + 1:ro + ni)
                case (4)
                    ro = ps(1)%ro(rp_e(q)); co = ps(1)%co(rp_e(q))    ! ro: costura (j), co: central (i)
                    gij = rg%glead_bot%a(co + 1:co + ni, ro + 1:ro + nj)
                    gji = ga1(ro + 1:ro + nj, co + 1:co + ni)
                case default
                    ro = ps(2)%ro(rp_e(q)); co = ps(2)%co(rp_e(q))
                    gij = rg%glead_top%a(ro + 1:ro + ni, co + 1:co + nj)
                    gji = gbn(co + 1:co + nj, ro + 1:ro + ni)
                end select
                ph = fac * exp(cmplx(0.0_dp, -dot_product(kc, rp_vec(:, q)), dp))
                gij = ph * (wz * gij - conjg(wz) * conjg(transpose(gji)))
                do w2 = 1, inp%rho_nwin
                    if (inw(w2)) rbuf(1:ni, 1:nj, q, w2, tid) = rbuf(1:ni, 1:nj, q, w2, tid) + gij
                end do
            end do
        end subroutine rho_core

        !> PREFIX.rho_pairs.dat (blocos D_ij reais por janela) e PREFIX.rho_mulliken.dat
        !> (carga de Mulliken por átomo central: sum_j tr[D_ij S_ji] e a integral do DOS).
        subroutine write_rho(lerr, lmsg)
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            real(dp), allocatable :: dsum(:, :, :, :), qds(:, :), qdos(:, :)
            integer :: u, w2, q, a, ic, o, ni, nj, p
            real(dp) :: dmax, tot1, tot2
            character(len=256) :: ln
            lerr = 0
            lmsg = ''
            allocate(dsum(nmo, nmo, nrp, inp%rho_nwin), qds(cent%norb, inp%rho_nwin), qdos(cent%norb, inp%rho_nwin))
            dsum = real(sum(rbuf, dim=5), dp)
            qdos = real(sum(qbuf, dim=3), dp)
            qds = 0.0_dp
            do w2 = 1, inp%rho_nwin
                do q = 1, nrp
                    ic = rp_ic(q)
                    if (ic < 1) cycle
                    ni = rp_ni(q); nj = rp_nj(q)
                    do a = 1, ni
                        o = cent%off(ic) + a
                        qds(o, w2) = qds(o, w2) + sum(dsum(a, 1:nj, q, w2) * rp_s(a, 1:nj, q))
                    end do
                end do
            end do
            if (present(rho_d)) rho_d = dsum(:, :, :, 1)
            if (present(rho_q)) then
                allocate(rho_q(cent%n, inp%rho_nwin))
                do w2 = 1, inp%rho_nwin
                    do ic = 1, cent%n
                        o = cent%off(ic)
                        rho_q(ic, w2) = sum(qds(o + 1:o + sup%orbs(cent%atom(ic)), w2))
                    end do
                end do
            end if
            if (present(rho_p)) then
                allocate(rho_p(3, cent%n, inp%rho_nwin))
                rho_p = 0.0_dp
                block
                    real(dp), allocatable :: sb(:, :), xb(:, :, :)
                    logical :: okd
                    integer :: c3
                    allocate(sb(nmo, nmo), xb(nmo, nmo, 3))
                    do q = 1, nrp
                        ic = rp_ic(q)
                        if (ic < 1) cycle
                        ni = rp_ni(q); nj = rp_nj(q)
                        call pao_dip_get(scc_dip, rp_vec(:, q) / BOHR_ANG_X, sb(1:ni, 1:nj), xb(1:ni, 1:nj, :), okd)
                        if (.not. okd) then
                            lerr = 4; lmsg = 'rho_p: dipole integrals missing for a pair'; return
                        end if
                        do w2 = 1, inp%rho_nwin
                            do c3 = 1, 3
                                rho_p(c3, ic, w2) = rho_p(c3, ic, w2) - sum(dsum(1:ni, 1:nj, q, w2) * xb(1:ni, 1:nj, c3))
                            end do
                        end do
                    end do
                end block
            end if
            if (bmode) then
                block
                    real(dp), allocatable :: sb(:, :), xb(:, :, :), qpr(:, :), b1(:), b2(:)
                    real(dp) :: tq(3), asym
                    logical :: okd
                    integer :: c3, ib
                    allocate(qpr(nrp, inp%rho_nwin), sb(nmo, nmo), xb(nmo, nmo, 3), b1(nbd), b2(nbd))
                    do w2 = 1, inp%rho_nwin
                        do q = 1, nrp
                            qpr(q, w2) = sum(dsum(1:rp_ni(q), 1:rp_nj(q), q, w2) * rp_s(1:rp_ni(q), 1:rp_nj(q), q))
                        end do
                    end do
                    if (present(rho_qon)) then
                        allocate(rho_qon(cent%n, inp%rho_nwin))
                        rho_qon = 0.0_dp
                        do q = 1, nrp
                            if (rp_bond(q) == 0) rho_qon(rp_ic(q), :) = rho_qon(rp_ic(q), :) + qpr(q, :)
                        end do
                    end if
                    if (present(rho_qb)) then
                        allocate(rho_qb(nbd, inp%rho_nwin))
                        rho_qb = 0.0_dp
                        do q = 1, nrp
                            if (rp_bond(q) > 0) rho_qb(rp_bond(q), :) = rho_qb(rp_bond(q), :) + &
                                real(rp_bw(q), dp) * qpr(q, :)
                        end do
                    end if
                    if (present(rho_pon)) then
                        allocate(rho_pon(3, cent%n, inp%rho_nwin))
                        rho_pon = 0.0_dp
                        do q = 1, nrp
                            if (rp_bond(q) /= 0) cycle
                            ic = rp_ic(q)
                            ni = rp_ni(q); nj = rp_nj(q)
                            call pao_dip_get(scc_dip, [0.0_dp, 0.0_dp, 0.0_dp], sb(1:ni, 1:nj), xb(1:ni, 1:nj, :), okd)
                            if (.not. okd) then
                                lerr = 4; lmsg = 'rho_pon: on-site dipole integrals missing'; return
                            end if
                            do w2 = 1, inp%rho_nwin
                                do c3 = 1, 3
                                    rho_pon(c3, ic, w2) = rho_pon(c3, ic, w2) - sum(dsum(1:ni, 1:nj, q, w2) * &
                                                          xb(1:ni, 1:nj, c3))
                                end do
                            end do
                        end do
                    end if
                    ! conferência (janela 1): sum q_on + sum q_b - metade das ligações com a costura = sum q_DS;
                    ! assimetria entre as duas ordens de cada ligação central
                    tq = 0.0_dp
                    b1 = 0.0_dp
                    b2 = 0.0_dp
                    do q = 1, nrp
                        if (rp_bond(q) == 0) then
                            tq(1) = tq(1) + qpr(q, 1)
                        else
                            tq(2) = tq(2) + real(rp_bw(q), dp) * qpr(q, 1)
                            if (rp_jc(q) < 1) tq(3) = tq(3) + qpr(q, 1)
                            if (rp_bw(q) == 1) b2(rp_bond(q)) = b2(rp_bond(q)) + qpr(q, 1)
                        end if
                    end do
                    asym = 0.0_dp
                    do ib = 1, nbd
                        if (rp_bw(bd_q(ib)) /= 1) cycle
                        b1(ib) = qpr(bd_q(ib), 1)
                        asym = max(asym, abs(2.0_dp * b1(ib) - b2(ib)))
                    end do
                    write(ln, '(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') 'bond: sum q_on + sum q_b - seam halves - ' // &
                        'sum q_DS = ', tq(1) + tq(2) - tq(3) - sum(qds(:, 1)), ' (sum q_DS ', sum(qds(:, 1)), &
                        '); max |q(i->j) - q(j->i)| = ', asym, '; sum q_b = ', tq(2)
                    call note(trim(ln))
                end block
            end if
            if (present(rho_files)) then
                if (.not. rho_files) return
            end if
            open(newunit=u, file=trim(inp%output_prefix) // '.rho_pairs.dat', status='replace', action='write')
            write(u, '(A)') '# surfgreen extracted &RHO: D_ij = Re sum_k w_k e^{-i k.d} sum_{E in window} ' // &
                '(i/2pi)(G_ij - G_ji^dag) dE (atomic gauge; d = r_j + R - r_i)'
            write(u, '(A)') '# pair: q kind i_atom i_col j_atom ni nj dx dy dz |d| xi yi zi label_i label_j ' // &
                '(Angstrom; r_i extended position; kind 4/5: j in the left/right seam)'
            write(u, '(A)') '# lattice: a x y z periodic(1 = k axis)'
            if (inp%rho_contour) then
                write(u, '(A,I0)') '# windows by contour (Gauss-Legendre in theta), npoints ', inp%rho_ncontour
            else
                write(u, '(A)') '# windows by real-axis sum'
            end if
            write(u, '(A,I0,A,I0,A,I0)') '# nwin ', inp%rho_nwin, ' npair ', nrp, ' nmo ', nmo
            do a = 1, 3
                write(u, '(A,I0,3(1X,ES24.16),1X,I0)') 'lattice ', a, sup%lattice(a, :), merge(1, 0, kaxis(a))
            end do
            do q = 1, nrp
                ic = max(rp_ic(q), 1)
                write(u, '(A,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,7(1X,ES24.16),2(1X,A))') 'pair ', q, &
                    rp_kind(q), cent%atom(ic), cent%x(ic), rp_ja(q), rp_ni(q), rp_nj(q), rp_vec(:, q), &
                    norm2(rp_vec(:, q)), ext_pos(sup, cent%atom(ic), cent%x(ic)), trim(sup%label(cent%atom(ic))), &
                    trim(sup%label(rp_ja(q)))
            end do
            do w2 = 1, inp%rho_nwin
                write(u, '(A,1X,A,2(1X,ES16.8))') 'window', trim(inp%rho_name(w2)), inp%rho_lo(w2), inp%rho_hi(w2)
                do q = 1, nrp
                    write(u, '(A,1X,I0)') 'block', q
                    do a = 1, rp_ni(q)
                        write(u, '(*(1X,ES24.16))') dsum(a, 1:rp_nj(q), q, w2)
                    end do
                end do
            end do
            close(u)
            open(newunit=u, file=trim(inp%output_prefix) // '.rho_mulliken.dat', status='replace', action='write')
            write(u, '(A)') '# surfgreen extracted &RHO: Mulliken charge per central atom and window'
            write(u, '(A)') '# q_DS = sum_j tr[D_ij S_ji] (pairs within rmax) ; q_DOS = window integral of the ' // &
                'Mulliken DOS (with seam cross terms)'
            write(u, '(A)') '# window ic atom col q_DS q_DOS'
            dmax = 0.0_dp
            do w2 = 1, inp%rho_nwin
                do ic = 1, cent%n
                    o = cent%off(ic)
                    p = sup%orbs(cent%atom(ic))
                    write(u, '(A,1X,I0,1X,I0,1X,I0,2(1X,ES24.16))') trim(inp%rho_name(w2)), ic, cent%atom(ic), &
                        cent%x(ic), sum(qds(o + 1:o + p, w2)), sum(qdos(o + 1:o + p, w2))
                    dmax = max(dmax, abs(sum(qds(o + 1:o + p, w2)) - sum(qdos(o + 1:o + p, w2))))
                end do
            end do
            close(u)
            tot1 = sum(qds); tot2 = sum(qdos)
            write(ln, '(A,ES10.3,A,2ES16.8)') 'rho: max |q_DS - q_DOS| per atom = ', dmax, '; totals ', tot1, tot2
            call note(trim(ln))
        end subroutine write_rho

        !> &KLIST: média uniforme no eixo avg para cada ponto da lista.
        subroutine reduce_klist()
            integer :: nl_, na, il_, ia, ipp
            real(dp), allocatable :: tot(:, :), am(:, :, :), res_(:, :), mn(:, :), sm(:), kf(:, :)
            integer, allocatable :: ni(:, :)
            nl_ = inp%klist_n
            na = inp%klist_navg
            allocate(tot(ne, nl_), am(ne, nl_, nreg), res_(ne, nl_), mn(ne, nl_), sm(nl_), ni(ne, nl_), kf(3, nl_))
            tot = 0.0_dp; am = 0.0_dp; res_ = 0.0_dp; mn = huge(1.0_dp); sm = huge(1.0_dp); ni = 0
            do il_ = 1, nl_
                do ia = 1, na
                    ipp = (il_ - 1) * na + ia
                    tot(:, il_) = tot(:, il_) + gd%total(:, ipp) / real(na, dp)
                    am(:, il_, :) = am(:, il_, :) + gd%amap(:, ipp, :) / real(na, dp)
                    res_(:, il_) = max(res_(:, il_), gd%residual(:, ipp))
                    mn(:, il_) = min(mn(:, il_), gd%min_neg_im_gii(:, ipp))
                    sm(il_) = min(sm(il_), gd%smin(ipp))
                    where (abs(gd%niter(:, ipp)) > abs(ni(:, il_)) .or. gd%niter(:, ipp) < 0) ni(:, il_) = gd%niter(:, ipp)
                end do
                kf(:, il_) = 0.0_dp
                kf(inp%klist_axis, il_) = inp%klist_k(il_)
            end do
            call move_alloc(tot, gd%total); call move_alloc(am, gd%amap); call move_alloc(res_, gd%residual)
            call move_alloc(mn, gd%min_neg_im_gii); call move_alloc(sm, gd%smin); call move_alloc(ni, gd%niter)
            call move_alloc(kf, gd%kfrac)
            gd%nk = nl_
            nk = nl_
            deallocate(gd%kcart, gd%kdist)
            allocate(gd%kcart(3, nl_), gd%kdist(nl_))
            do il_ = 1, nl_
                gd%kcart(:, il_) = matmul(gd%kfrac(:, il_), recip)
            end do
            gd%kdist = 0.0_dp
            kw = inp%klist_w(1:nl_) / sum(inp%klist_w(1:nl_))
        end subroutine reduce_klist

        subroutine fill_meta()
            character(len=128) :: buf
            integer :: sd_
            call meta_add(gd, 'program', 'surfgreen')
            call meta_add(gd, 'format', 'green.bin v1')
            call meta_add(gd, 'mode', 'extracted')
            write(buf, '(ES16.8)') inp%eta_ev
            call meta_add(gd, 'eta_ev', adjustl(buf))
            write(buf, '(ES12.4)') inp%decim_tol
            call meta_add(gd, 'decim_tol', adjustl(buf))
            write(buf, '(I0)') inp%decim_maxit
            call meta_add(gd, 'decim_maxit', adjustl(buf))
            call meta_add(gd, 'projection', trim(inp%projection))
            call meta_add(gd, 'mulliken_lead_cross_terms', trim(merge('yes', 'no ', inp%lead_cross)))
            call meta_add(gd, 'map_gauge', 'atomic (region DOS is gauge invariant)')
            if (inp%kmode == 'list') then
                write(buf, '(A,I0,A,I0,A,I0)') 'list axis ', inp%klist_axis, ', average axis ', inp%klist_avg_axis, &
                    ', navg ', inp%klist_navg
                call meta_add(gd, 'klist', trim(buf) // ' (uniform k = -1/2 + j/navg; maps averaged)')
                call meta_add(gd, 'klist_weights', reals(kw))
            end if
            if (have_v) then
                call meta_add(gd, 'site_potential_file', trim(inp%sitepot_file))
                call meta_add(gd, 'site_potential_mode', trim(inp%sitepot_mode) // &
                    ' (anticommutator: H_ij += (V_i+V_j)/2 S_ij in all central and seam blocks; onsite: H_ii += V_i S_ii)')
                write(buf, '(2ES16.8)') inp%sitepot_uniform_ev, inp%sitepot_lead_ev
                call meta_add(gd, 'site_potential_uniform_lead_ev', adjustl(buf))
            end if
            call meta_add(gd, 'supercell_h', trim(inp%h_file))
            call meta_add(gd, 'supercell_s', trim(inp%s_file))
            call meta_add(gd, 'supercell_structure_fingerprint', structure_fingerprint(sup%st))
            write(buf, '(I0,A,I0,A,I0)') e, ' (vacuum axes ', vac(1), ' ', vac(2)
            call meta_add(gd, 'embed_axis', trim(buf) // ')')
            write(buf, '(I0,A,I0,A,F10.6,A,I0)') clo, ' .. ', chi, '  (column_shift ', sup%shift, ', ncx ', sup%ncx
            call meta_add(gd, 'central_columns', trim(buf) // ')')
            call meta_add(gd, 'rgf_block_columns', ints(bwid))
            call meta_add(gd, 'leads_left_right', trim(inp%left) // ' ' // trim(inp%right))
            do sd_ = 1, 2
                if (.not. use_side(sd_)) cycle
                associate(ld => leads(lmat(sd_)))
                    call meta_add(gd, merge('lead_left_h ', 'lead_right_h', sd_ == 1), &
                                  trim(inp%lead_h(lmat(sd_))))
                    write(buf, '(I0,A,I0,A)') ld%P, ' cells x ', ld%sys%norb, ' orbitals'
                    call meta_add(gd, merge('lead_left_pl ', 'lead_right_pl', sd_ == 1), trim(buf))
                    write(buf, '(ES20.12)') ld%delta
                    call meta_add(gd, merge('delta_left_ev ', 'delta_right_ev', sd_ == 1), adjustl(buf))
                    write(buf, '(3F12.6,A,ES10.3)') match_t(:, sd_), '  max_dist ', match_d(sd_)
                    call meta_add(gd, merge('seam_left_translation_ang ', 'seam_right_translation_ang', sd_ == 1), &
                                  adjustl(buf))
                end associate
            end do
            call meta_add(gd, 'delta_mode', trim(inp%delta_mode))
            if (trim(inp%delta_mode) == 'fit') then
                write(buf, '(2ES12.4)') maxval(fit_res(1, :)), maxval(fit_res(2, :))
                call meta_add(gd, 'delta_fit_residual_onsite_h_s', adjustl(buf))
            end if
            call meta_add(gd, 'match_translation', trim(inp%match_translation))
            call meta_add(gd, 'energy_origin', trim(inp%energy_origin))
            write(buf, '(ES20.12)') eorigin
            call meta_add(gd, 'energy_origin_ev', adjustl(buf))
            if (trim(inp%energy_origin) == 'absolute') then
                if (sup%has_meta .and. sup%shifted) then
                    call meta_add(gd, 'energy_reference', 'E-E_F')
                else if (sup%has_meta) then
                    call meta_add(gd, 'energy_reference', 'absolute')
                else
                    call meta_add(gd, 'energy_reference', 'unknown')
                end if
            else
                call meta_add(gd, 'energy_reference', 'E-E_F')
            end if
            if (sup%has_meta) then
                write(buf, '(ES20.12)') sup%fermi_ev
                call meta_add(gd, 'supercell_fermi_ev', adjustl(buf))
            end if
            if (leads(1)%sys%has_meta) then
                write(buf, '(ES20.12)') leads(1)%sys%fermi_ev
                call meta_add(gd, 'lead_fermi_ev', adjustl(buf))
            end if
            call meta_add(gd, 'central_structure', central_file)
            call meta_add(gd, 'central_structure_fingerprint', structure_fingerprint(stc))
            if (len_trim(inp%regions_file) > 0) then
                call meta_add(gd, 'regions_file', trim(inp%regions_file))
            else
                call meta_add(gd, 'regions_file', '(default: one region per ' // &
                              trim(merge('column', 'atom  ', trim(inp%region_mode) == 'columns')) // ')')
            end if
        end subroutine fill_meta

    end subroutine run_extracted

    !> Ajuste de Delta: H_lead + Delta S_lead ~ H_sup nos blocos (c, c+p),
    !> |p| <= r, das colunas c dadas, nos pontos kc. Delta = sum Re<S_l, H_s - H_l>
    !> / sum |S_l|^2 (kcs: k fracionários). res_h, res_s: max |H_s - H_l - Delta S_l| e max |S_s - S_l|
    !> nos blocos de sítio (p = 0), como nos protótipos.
    subroutine fit_delta(sup, lead, cols, r, kcs, kaxis, tol, mode, delta, res_h, res_s, ierr, errmsg)
        type(xsys_t),     intent(in)  :: sup
        type(lead_t),     intent(in)  :: lead
        integer,          intent(in)  :: cols(:), r
        real(dp),         intent(in)  :: kcs(:, :), tol
        logical,          intent(in)  :: kaxis(3)
        character(len=*), intent(in)  :: mode
        real(dp),         intent(out) :: delta, res_h, res_s
        integer,          intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        type(xlist_t) :: grp, ra, rb, la, lb
        type(plan_t), allocatable :: psup(:, :), plead(:)
        complex(dp), allocatable :: hs(:, :), ss(:, :), hl(:, :), sl(:, :)
        real(dp) :: num, den, t(3), dmax
        integer :: ic, p, ik, nl, n

        ierr = 0
        errmsg = ''
        nl = lead%sys%nat
        n = lead%pl0%norb / lead%P
        allocate(psup(-r:r, size(cols)), plead(-r:r))
        call columns_list(lead%sys, lead%c0, lead%c0, la)
        do p = -r, r
            call columns_list(lead%sys, lead%c0 + p, lead%c0 + p, lb)
            call build_plan(lead%sys, la, lb, plead(p))
        end do
        do ic = 1, size(cols)
            call match_columns(sup, lead, cols(ic) - r, 2 * r + 1, kaxis, tol, mode, grp, t, dmax, ierr, errmsg)
            if (ierr /= 0) return
            call sub_list(sup, grp, r * nl + 1, (r + 1) * nl, ra)
            do p = -r, r
                call sub_list(sup, grp, (r + p) * nl + 1, (r + p + 1) * nl, rb)
                call build_plan(sup, ra, rb, psup(p, ic))
            end do
        end do
        allocate(hs(n, n), ss(n, n), hl(n, n), sl(n, n))
        num = 0.0_dp
        den = 0.0_dp
        do ik = 1, size(kcs, 2)
            do p = -r, r
                call fill_plan(lead%sys, plead(p), kcs(:, ik), hl, sl)
                do ic = 1, size(cols)
                    call fill_plan(sup, psup(p, ic), kcs(:, ik), hs, ss)
                    num = num + real(sum(conjg(sl) * (hs - hl)), dp)
                    den = den + real(sum(conjg(sl) * sl), dp)
                end do
            end do
        end do
        if (den <= 0.0_dp) then
            ierr = 2; errmsg = 'zero overlap in the fit blocks'; return
        end if
        delta = num / den
        res_h = 0.0_dp
        res_s = 0.0_dp
        do ik = 1, size(kcs, 2)
            call fill_plan(lead%sys, plead(0), kcs(:, ik), hl, sl)
            do ic = 1, size(cols)
                call fill_plan(sup, psup(0, ic), kcs(:, ik), hs, ss)
                res_h = max(res_h, maxval(abs(hs - hl - delta * sl)))
                res_s = max(res_s, maxval(abs(ss - sl)))
            end do
        end do
    end subroutine fit_delta

    ! ==================================================================
    ! Auxiliares
    ! ==================================================================

    function column_counts(sys) result(s)
        type(xsys_t), intent(in) :: sys
        character(len=:), allocatable :: s
        integer :: c
        s = ''
        do c = minval(sys%col), maxval(sys%col)
            s = s // ' ' // itoa(count(sys%col == c))
        end do
    end function column_counts

    function ints(v) result(s)
        integer, intent(in) :: v(:)
        character(len=:), allocatable :: s
        integer :: i
        s = ''
        do i = 1, size(v)
            s = s // ' ' // itoa(v(i))
        end do
    end function ints

    pure function other_axes(a) result(o)
        integer, intent(in) :: a
        integer :: o(2)
        select case (a)
        case (1);     o = [2, 3]
        case (2);     o = [1, 3]
        case default; o = [1, 2]
        end select
    end function other_axes

    pure integer(int64) function blk_key(bi, nat)
        integer, intent(in) :: bi(3), nat
        blk_key = (int(bi(3) - 1, int64) * nat + (bi(1) - 1)) * nat + (bi(2) - 1)
    end function blk_key

    !> Recíproca (linhas, com 2 pi) de uma base B (linhas).
    pure function recip_of(b) result(r)
        real(dp), intent(in) :: b(3, 3)
        real(dp) :: r(3, 3)
        r = 6.283185307179586476925286766559_dp * transpose(inv3(b))
    end function recip_of

    pure function det3(a) result(d)
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: d
        d = a(1, 1) * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)) - a(1, 2) * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1)) &
            + a(1, 3) * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
    end function det3

    pure function inv3(a) result(b)
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: b(3, 3), d
        d = det3(a)
        b(1, 1) = a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)
        b(1, 2) = a(1, 3) * a(3, 2) - a(1, 2) * a(3, 3)
        b(1, 3) = a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)
        b(2, 1) = a(2, 3) * a(3, 1) - a(2, 1) * a(3, 3)
        b(2, 2) = a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1)
        b(2, 3) = a(1, 3) * a(2, 1) - a(1, 1) * a(2, 3)
        b(3, 1) = a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1)
        b(3, 2) = a(1, 2) * a(3, 1) - a(1, 1) * a(3, 2)
        b(3, 3) = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
        b = b / d
    end function inv3

    !> Permutação estável que ordena key (merge sort de baixo para cima).
    subroutine stable_order64(key, order)
        integer(int64),       intent(in)  :: key(:)
        integer, allocatable, intent(out) :: order(:)
        integer, allocatable :: tmp(:)
        integer :: n, width, lo, mid, hi, i, j, k
        n = size(key)
        allocate(order(n), tmp(n))
        order = [(i, i = 1, n)]
        width = 1
        do while (width < n)
            lo = 1
            do while (lo <= n)
                mid = min(lo + width - 1, n)
                hi = min(lo + 2 * width - 1, n)
                i = lo; j = mid + 1; k = lo
                do while (i <= mid .and. j <= hi)
                    if (key(order(j)) < key(order(i))) then
                        tmp(k) = order(j); j = j + 1
                    else
                        tmp(k) = order(i); i = i + 1
                    end if
                    k = k + 1
                end do
                do while (i <= mid)
                    tmp(k) = order(i); i = i + 1; k = k + 1
                end do
                do while (j <= hi)
                    tmp(k) = order(j); j = j + 1; k = k + 1
                end do
                lo = lo + 2 * width
            end do
            order = tmp
            width = 2 * width
        end do
    end subroutine stable_order64

    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=16) :: buf
        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

    !> h(i,j) += (ur(i) + uc(j))/2 s(i,j): potencial escalar por orbital na
    !> forma 1/2 {S, V} (deslocamento rígido H + c S quando ur = uc = c).
    subroutine vpot_cross(h, s, ur, uc)
        complex(dp), intent(inout) :: h(:, :)
        complex(dp), intent(in) :: s(:, :)
        real(dp), intent(in) :: ur(:), uc(:)
        integer :: i, j
        do j = 1, size(h, 2)
            do i = 1, size(h, 1)
                h(i, j) = h(i, j) + cmplx(0.5_dp * (ur(i) + uc(j)), 0.0_dp, dp) * s(i, j)
            end do
        end do
    end subroutine vpot_cross

    !> Lista de reais para os metadados.
    function reals(x) result(txt)
        real(dp), intent(in) :: x(:)
        character(len=:), allocatable :: txt
        character(len=32) :: b
        integer :: i
        txt = ''
        do i = 1, size(x)
            write(b, '(ES14.6)') x(i)
            txt = txt // ' ' // trim(adjustl(b))
        end do
    end function reals

    !> Nós e pesos de Gauss-Legendre em [-1, 1] (Newton na recorrência de Legendre).
    subroutine gauss_legendre(n, x, w)
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        integer :: i, j, it
        real(dp) :: z, p1, p2, p3, pp, z1
        do i = 1, (n + 1) / 2
            z = cos(PI * (real(i, dp) - 0.25_dp) / (real(n, dp) + 0.5_dp))
            do it = 1, 100
                p1 = 1.0_dp; p2 = 0.0_dp
                do j = 1, n
                    p3 = p2; p2 = p1
                    p1 = ((2.0_dp * j - 1.0_dp) * z * p2 - (j - 1.0_dp) * p3) / real(j, dp)
                end do
                pp = n * (z * p1 - p2) / (z * z - 1.0_dp)
                z1 = z
                z = z1 - p1 / pp
                if (abs(z - z1) < 1.0e-15_dp) exit
            end do
            x(i) = -z; x(n + 1 - i) = z
            w(i) = 2.0_dp / ((1.0_dp - z * z) * pp * pp); w(n + 1 - i) = w(i)
        end do
    end subroutine gauss_legendre

end module extracted_green
