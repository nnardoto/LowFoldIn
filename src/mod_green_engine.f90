! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Motor do programa surfgreen: leitura da entrada, montagem do empilhamento,
!> laço sobre (k_par, E) e preenchimento do green_data_t.
!>
!> Empilhamento (de baixo para cima ao longo de +a_s), material A = substrato,
!> B = filme:
!>   surface:  side='top'    [A x nsub]                      meio A abaixo, vácuo acima
!>             side='bottom' [A x nsub]                      vácuo abaixo, meio A acima
!>   film:     side='top'    [A x nsub][B x N]               meio A abaixo, vácuo acima
!>             side='bottom' [B x N][A x nsub]               vácuo abaixo, meio A acima
!>   junction: side='top'    [A x nsub][B x N]               meio A abaixo, meio B acima
!>             side='bottom' [B x N][A x nsub]               meio B abaixo, meio A acima
!>   sandwich:               [A x nsub][B x N][A x nsub]     meio A nos dois lados
!> Entre PLs do mesmo material o acoplamento é o H_{m,m+1} do próprio
!> material; entre A e B é o bloco de contato de interface_model.
!>
!> Energias: V de cada camada = shift_ev do material (+ layer_potential_ev(m)
!> nas camadas do filme); H_mn -> H_mn + (V_m+V_n)/2 S_mn em todos os blocos,
!> inclusive nos meios semi-infinitos (deslocamento rígido H + Delta S).
module green_engine
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use types, only: dp
    use HS_io, only: HSDataFormat
    use structure_io, only: structure_t, read_structure, structure_fingerprint, reciprocal_lattice, LABEL_LEN
    use kpath, only: kpath_t, kgrid_t, build_kpath, build_kgrid
    use atomic_regions, only: region_set_t, read_regions_file, validate_regions, REGION_NAME_LEN
    use spectral_projection, only: orbital_groups_t, build_orbital_groups
    use layer_blocks, only: pl_blocks_t, build_pl_blocks, pl_blocks_at_k, attach_pl_geometry, kpar_cartesian
    use surface_green, only: surface_gf, side_code, SIDE_TOP, SIDE_BOTTOM
    use interface_model, only: contact_t, contact_validate, build_contact, REF_SUBSTRATE, REF_FILM, CONTACT_SK
    use sk_contact, only: sk_model_t, sk_geom_t, sk_pairs_t, read_sk_file, sk_set_bases, make_sk_geom, &
                          build_sk_pairs, sk_contact_at_k
    use central_green, only: central_t, rgf_t, central_init, assemble_central, attach_lead, rgf_solve
    use green_projection, only: orbital_dos_mulliken, orbital_dos_lowdin, region_dos, coherent_map, &
                                min_neg_im_gdiag, central_overlap_min
    use green_data, only: green_data_t, meta_add, GREEN_KPATH, GREEN_KGRID
    !$ use omp_lib, only: omp_get_max_threads
    implicit none
    private

    integer, parameter, public :: GPATH_LEN = 512, MAX_GKVERTEX = 64, MAX_GLAYERS = 512
    integer, parameter, public :: MAX_SK_BASIS = 32, MAX_SCAN = 256

    !> Parâmetros que &SCAN pode varrer.
    character(len=*), parameter :: SCAN_NAMES = &
        'lambda_h lambda_s sk_scale_sigma sk_scale_pi sk_scale_delta sk_beta film_shift_ev ' // &
        'offset_x offset_y offset_z offset_normal'

    !> Espelho dos namelists &GREEN, &SUBSTRATE, &FILM, &INTERFACE, &COHERENT,
    !> &KPATH ou &KGRID.
    type, public :: green_input_t
        ! &GREEN
        character(len=16)        :: mode = 'surface'
        real(dp)                 :: eta_ev = 0.002_dp
        real(dp)                 :: emin = -1.0_dp, emax = 1.0_dp
        integer                  :: ne = 101
        real(dp)                 :: decim_tol = 1.0e-10_dp
        integer                  :: decim_maxit = 200
        character(len=16)        :: projection = 'mulliken'
        character(len=GPATH_LEN) :: output_prefix = 'surfgreen'
        character(len=GPATH_LEN) :: regions_file = ''
        logical                  :: allow_no_fingerprint = .true.
        real(dp)                 :: vacuum_ang = 20.0_dp
        real(dp)                 :: lattice_tol_ang = 0.05_dp
        ! &SUBSTRATE
        character(len=GPATH_LEN) :: sub_h = '', sub_s = '', sub_structure = ''
        integer                  :: stacking_axis = 3
        integer                  :: sub_pl_cells = 0
        character(len=16)        :: side = 'top'
        real(dp)                 :: sub_shift_ev = 0.0_dp
        integer                  :: sub_layers = -1        !< PLs de A na região central (-1: 1 em surface, 0 nos demais)
        ! &FILM
        logical                  :: have_film = .false.
        character(len=GPATH_LEN) :: film_h = '', film_s = '', film_structure = ''
        integer                  :: nlayers = 0
        integer                  :: film_ncells = 0       !< > 0: espessura do filme em células (BL); sobrepõe nlayers
        integer                  :: film_pl_cells = 0
        real(dp)                 :: film_shift_ev = 0.0_dp
        real(dp)                 :: layer_potential_ev(MAX_GLAYERS) = 0.0_dp
        real(dp)                 :: interface_offset(3) = 0.0_dp
        ! &INTERFACE
        character(len=16)        :: contact_model = 'none'
        character(len=16)        :: contact_reference = 'substrate'
        real(dp)                 :: lambda_h = 1.0_dp, lambda_s = 1.0_dp
        character(len=GPATH_LEN) :: sk_file = ''
        real(dp)                 :: sk_cutoff_ang = 6.0_dp
        character(len=64)        :: sk_basis(MAX_SK_BASIS) = ''
        real(dp)                 :: sk_scale(3) = 1.0_dp      !< sigma, pi, delta
        real(dp)                 :: sk_beta = -1.0_dp         !< >= 0 substitui o beta do arquivo
        ! &SCAN
        character(len=32)        :: scan_parameter = ''
        integer                  :: nscan = 0
        real(dp)                 :: scan_values(MAX_SCAN) = 0.0_dp
        ! &COHERENT
        integer                  :: nqz = 0
        real(dp)                 :: qz_min = 0.0_dp, qz_max = 0.0_dp
        real(dp)                 :: escape_length_ang = 0.0_dp   !< 0: sem atenuação
        character(len=REGION_NAME_LEN) :: coherent_region = ''
        ! &KPATH / &KGRID
        character(len=8)         :: kmode = 'path'
        integer                  :: nvertex = 0
        real(dp)                 :: vertex(3, MAX_GKVERTEX) = 0.0_dp
        character(len=16)        :: label(MAX_GKVERTEX) = ''
        integer                  :: npoints(MAX_GKVERTEX - 1) = 0
        real(dp)                 :: grid_origin(3) = 0.0_dp, grid_span1(3) = 0.0_dp, grid_span2(3) = 0.0_dp
        integer                  :: grid_n1 = 0, grid_n2 = 0
    end type green_input_t

    !> Um material: blocos de PL, estrutura, deslocamento.
    type :: material_t
        type(pl_blocks_t) :: pl
        type(structure_t) :: st
        type(HSDataFormat) :: hh
        real(dp) :: shift = 0.0_dp
        character(len=1) :: tag = 'A'
    end type material_t

    !> Empilhamento resolvido.
    type, public :: green_stack_t
        integer :: nlayer = 0
        integer, allocatable :: mat(:)          !< (nlayer) 1 = A, 2 = B
        integer, allocatable :: ncell(:)        !< (nlayer) células no bloco (= P, ou menos no bloco parcial do topo)
        integer, allocatable :: film_index(:)   !< (nlayer) índice da camada no filme (0 se A)
        real(dp), allocatable :: v(:)           !< (nlayer) potencial
        integer :: bot_lead = 0, top_lead = 0   !< 0 = vácuo
        real(dp), allocatable :: origin(:, :)   !< (3, nlayer) origem cartesiana da camada
        ! estrutura da região central
        type(structure_t) :: st
        integer, allocatable :: atom_orbs(:)    !< (nat) orbitais por átomo central
        integer, allocatable :: atom_layer(:)   !< (nat)
    end type green_stack_t

    public :: read_green_input, run_green, build_stack, write_central_structure

    real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp

contains

    ! ==================================================================
    ! Entrada
    ! ==================================================================

    !> Lê os namelists. &GREEN e &SUBSTRATE são obrigatórios; &FILM é
    !> obrigatório fora do modo surface; &INTERFACE e &COHERENT são opcionais;
    !> exatamente um entre &KPATH e &KGRID.
    subroutine read_green_input(filename, inp, ierr, errmsg)
        character(len=*),    intent(in)  :: filename
        type(green_input_t), intent(out) :: inp
        integer,             intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: iu, ios
        logical :: exists, have_path, have_grid, found
        character(len=512) :: iomsg

        ierr = 0
        errmsg = ''
        inquire(file=filename, exist=exists)
        if (.not. exists) then
            ierr = 1; errmsg = "input file not found: '" // trim(filename) // "'"; return
        end if
        open(newunit=iu, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = "cannot open input file: '" // trim(filename) // "'"; return
        end if

        call read_green_nml(found)
        if (ierr /= 0) goto 900
        if (.not. found) then
            ierr = 2; errmsg = 'namelist &GREEN not found'; goto 900
        end if
        call read_substrate_nml(found)
        if (ierr /= 0) goto 900
        if (.not. found) then
            ierr = 2; errmsg = 'namelist &SUBSTRATE not found'; goto 900
        end if
        call read_film_nml(inp%have_film)
        if (ierr /= 0) goto 900
        call read_interface_nml(found)
        if (ierr /= 0) goto 900
        call read_scan_nml(found)
        if (ierr /= 0) goto 900
        call read_coherent_nml(found)
        if (ierr /= 0) goto 900
        call read_kpath_nml(have_path)
        if (ierr /= 0) goto 900
        call read_kgrid_nml(have_grid)
        if (ierr /= 0) goto 900
        if (have_path .eqv. have_grid) then
            ierr = 2
            errmsg = 'exactly one of &KPATH and &KGRID must be present'
            goto 900
        end if
        inp%kmode = merge('grid', 'path', have_grid)
        close(iu)
        call validate_input(inp, ierr, errmsg)
        return
900     close(iu)

    contains

        subroutine rew()
            rewind(iu)
        end subroutine rew

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
            character(len=GPATH_LEN) :: output_prefix, regions_file
            logical :: allow_no_fingerprint
            namelist /GREEN/ mode, eta_ev, emin, emax, ne, decim_tol, decim_maxit, projection, &
                             output_prefix, regions_file, allow_no_fingerprint, vacuum_ang, lattice_tol_ang
            mode = inp%mode; eta_ev = inp%eta_ev; emin = inp%emin; emax = inp%emax; ne = inp%ne
            decim_tol = inp%decim_tol; decim_maxit = inp%decim_maxit; projection = inp%projection
            output_prefix = inp%output_prefix; regions_file = inp%regions_file
            allow_no_fingerprint = inp%allow_no_fingerprint; vacuum_ang = inp%vacuum_ang
            lattice_tol_ang = inp%lattice_tol_ang
            call rew()
            read(iu, nml=GREEN, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('GREEN'); return
            end if
            inp%mode = mode; inp%eta_ev = eta_ev; inp%emin = emin; inp%emax = emax; inp%ne = ne
            inp%decim_tol = decim_tol; inp%decim_maxit = decim_maxit; inp%projection = projection
            inp%output_prefix = output_prefix; inp%regions_file = regions_file
            inp%allow_no_fingerprint = allow_no_fingerprint; inp%vacuum_ang = vacuum_ang
            inp%lattice_tol_ang = lattice_tol_ang
        end subroutine read_green_nml

        subroutine read_substrate_nml(found)
            logical, intent(out) :: found
            character(len=GPATH_LEN) :: h_file, s_file, structure_file
            integer :: stacking_axis, pl_cells, central_layers
            character(len=16) :: side
            real(dp) :: shift_ev
            namelist /SUBSTRATE/ h_file, s_file, structure_file, stacking_axis, pl_cells, side, &
                                 shift_ev, central_layers
            h_file = ''; s_file = ''; structure_file = ''
            stacking_axis = inp%stacking_axis; pl_cells = inp%sub_pl_cells; side = inp%side
            shift_ev = inp%sub_shift_ev; central_layers = inp%sub_layers
            call rew()
            read(iu, nml=SUBSTRATE, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('SUBSTRATE'); return
            end if
            inp%sub_h = h_file; inp%sub_s = s_file; inp%sub_structure = structure_file
            inp%stacking_axis = stacking_axis; inp%sub_pl_cells = pl_cells; inp%side = side
            inp%sub_shift_ev = shift_ev; inp%sub_layers = central_layers
        end subroutine read_substrate_nml

        subroutine read_film_nml(found)
            logical, intent(out) :: found
            character(len=GPATH_LEN) :: h_file, s_file, structure_file
            integer :: nlayers, pl_cells, ncells
            real(dp) :: shift_ev, layer_potential_ev(MAX_GLAYERS), interface_offset(3)
            namelist /FILM/ h_file, s_file, structure_file, nlayers, pl_cells, shift_ev, &
                            layer_potential_ev, interface_offset, ncells
            h_file = ''; s_file = ''; structure_file = ''
            nlayers = 0; pl_cells = 0; ncells = 0; shift_ev = 0.0_dp; layer_potential_ev = 0.0_dp
            interface_offset = 0.0_dp
            call rew()
            read(iu, nml=FILM, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('FILM'); return
            end if
            inp%film_h = h_file; inp%film_s = s_file; inp%film_structure = structure_file
            inp%nlayers = nlayers; inp%film_pl_cells = pl_cells; inp%film_shift_ev = shift_ev
            inp%film_ncells = ncells
            inp%layer_potential_ev = layer_potential_ev; inp%interface_offset = interface_offset
        end subroutine read_film_nml

        subroutine read_interface_nml(found)
            logical, intent(out) :: found
            character(len=16) :: model, reference
            real(dp) :: lambda_h, lambda_s
            character(len=GPATH_LEN) :: junction_h, junction_s, junction_structure, sk_file
            real(dp) :: sk_cutoff_ang, sk_scale_sigma, sk_scale_pi, sk_scale_delta, sk_beta
            character(len=64) :: sk_basis(MAX_SK_BASIS)
            namelist /INTERFACE/ model, reference, lambda_h, lambda_s, junction_h, junction_s, &
                                 junction_structure, sk_file, sk_cutoff_ang, sk_basis, sk_scale_sigma, &
                                 sk_scale_pi, sk_scale_delta, sk_beta
            model = inp%contact_model; reference = inp%contact_reference
            lambda_h = inp%lambda_h; lambda_s = inp%lambda_s
            junction_h = ''; junction_s = ''; junction_structure = ''
            sk_file = inp%sk_file; sk_cutoff_ang = inp%sk_cutoff_ang; sk_basis = inp%sk_basis
            sk_scale_sigma = inp%sk_scale(1); sk_scale_pi = inp%sk_scale(2); sk_scale_delta = inp%sk_scale(3)
            sk_beta = inp%sk_beta
            call rew()
            read(iu, nml=INTERFACE, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('INTERFACE'); return
            end if
            inp%contact_model = model; inp%contact_reference = reference
            inp%lambda_h = lambda_h; inp%lambda_s = lambda_s
            inp%sk_file = sk_file; inp%sk_cutoff_ang = sk_cutoff_ang; inp%sk_basis = sk_basis
            inp%sk_scale = [sk_scale_sigma, sk_scale_pi, sk_scale_delta]; inp%sk_beta = sk_beta
        end subroutine read_interface_nml

        !> &SCAN: parameter + values(:) explícitos, ou n valores igualmente
        !> espaçados de vmin a vmax.
        subroutine read_scan_nml(found)
            logical, intent(out) :: found
            character(len=32) :: parameter
            real(dp) :: values(MAX_SCAN), vmin, vmax
            integer :: n, k
            namelist /SCAN/ parameter, values, n, vmin, vmax
            parameter = ''; values = huge(1.0_dp); n = 0; vmin = 0.0_dp; vmax = 0.0_dp
            call rew()
            read(iu, nml=SCAN, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('SCAN'); return
            end if
            if (.not. found) return
            inp%scan_parameter = adjustl(parameter)
            if (n > 0) then
                if (n > MAX_SCAN) then
                    ierr = 2; errmsg = '&SCAN n must be <= 256'; return
                end if
                inp%nscan = n
                do k = 1, n
                    if (n == 1) then
                        inp%scan_values(k) = vmin
                    else
                        inp%scan_values(k) = vmin + real(k - 1, dp) * (vmax - vmin) / real(n - 1, dp)
                    end if
                end do
            else
                inp%nscan = count(values < huge(1.0_dp))
                if (any(values(1:inp%nscan) == huge(1.0_dp))) then
                    ierr = 2; errmsg = '&SCAN values must be contiguous (values(1), values(2), ...)'; return
                end if
                inp%scan_values(1:inp%nscan) = values(1:inp%nscan)
            end if
            if (inp%nscan < 1) then
                ierr = 2; errmsg = '&SCAN needs values(:) or n >= 1 with vmin, vmax'; return
            end if
        end subroutine read_scan_nml

        subroutine read_coherent_nml(found)
            logical, intent(out) :: found
            integer :: nqz
            real(dp) :: qz_min, qz_max, escape_length_ang
            character(len=REGION_NAME_LEN) :: region
            namelist /COHERENT/ nqz, qz_min, qz_max, escape_length_ang, region
            nqz = 0; qz_min = 0.0_dp; qz_max = 0.0_dp; escape_length_ang = 0.0_dp; region = ''
            call rew()
            read(iu, nml=COHERENT, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('COHERENT'); return
            end if
            inp%nqz = nqz; inp%qz_min = qz_min; inp%qz_max = qz_max
            inp%escape_length_ang = escape_length_ang; inp%coherent_region = region
        end subroutine read_coherent_nml

        subroutine read_kpath_nml(found)
            logical, intent(out) :: found
            integer :: nvertex
            real(dp) :: vertex(3, MAX_GKVERTEX)
            character(len=16) :: label(MAX_GKVERTEX)
            integer :: npoints(MAX_GKVERTEX - 1)
            namelist /KPATH/ nvertex, vertex, label, npoints
            nvertex = 0; vertex = 0.0_dp; label = ''; npoints = 0
            call rew()
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
            call rew()
            read(iu, nml=KGRID, iostat=ios, iomsg=iomsg)
            found = (ios /= iostat_end)
            if (found .and. ios /= 0) then
                call nml_error('KGRID'); return
            end if
            inp%grid_origin = origin; inp%grid_span1 = span1; inp%grid_span2 = span2
            inp%grid_n1 = n1; inp%grid_n2 = n2
        end subroutine read_kgrid_nml

    end subroutine read_green_input

    subroutine validate_input(inp, ierr, errmsg)
        type(green_input_t), intent(inout) :: inp
        integer,             intent(out)   :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        ierr = 2
        errmsg = ''
        select case (trim(inp%mode))
        case ('surface', 'film', 'junction', 'sandwich')
        case default
            errmsg = "&GREEN mode must be 'surface', 'film', 'junction' or 'sandwich'"; return
        end select
        if (inp%eta_ev <= 0.0_dp) then
            errmsg = '&GREEN eta_ev must be > 0'; return
        end if
        if (inp%ne < 1) then
            errmsg = '&GREEN ne must be >= 1'; return
        end if
        if (inp%ne > 1 .and. inp%emax <= inp%emin) then
            errmsg = '&GREEN emax must be > emin'; return
        end if
        if (inp%decim_tol <= 0.0_dp .or. inp%decim_maxit < 1) then
            errmsg = '&GREEN decim_tol and decim_maxit must be positive'; return
        end if
        select case (trim(inp%projection))
        case ('mulliken', 'lowdin')
        case default
            errmsg = "&GREEN projection must be 'mulliken' or 'lowdin'"; return
        end select
        if (len_trim(inp%sub_h) == 0 .or. len_trim(inp%sub_structure) == 0) then
            errmsg = '&SUBSTRATE h_file and structure_file are required'; return
        end if
        if (inp%stacking_axis < 1 .or. inp%stacking_axis > 3) then
            errmsg = '&SUBSTRATE stacking_axis must be 1, 2 or 3'; return
        end if
        if (side_code(inp%side) == 0) then
            errmsg = "&SUBSTRATE side must be 'top' or 'bottom'"; return
        end if
        if (inp%sub_layers < 0) inp%sub_layers = merge(1, 0, trim(inp%mode) == 'surface')
        if (trim(inp%mode) == 'surface') then
            if (inp%sub_layers < 1) then
                errmsg = "&SUBSTRATE central_layers must be >= 1 in mode 'surface'"; return
            end if
        else
            if (.not. inp%have_film) then
                errmsg = "&FILM is required in mode '" // trim(inp%mode) // "'"; return
            end if
            if (len_trim(inp%film_h) == 0 .or. len_trim(inp%film_structure) == 0) then
                errmsg = '&FILM h_file and structure_file are required'; return
            end if
            if (inp%film_ncells < 0 .or. inp%film_ncells > MAX_GLAYERS) then
                errmsg = '&FILM ncells must be in 0..512'; return
            end if
            if (inp%film_ncells > 0) then
                ! espessura em células: o bloco parcial fica no topo do filme, do lado do vácuo
                if (trim(inp%mode) /= 'film' .or. side_code(inp%side) /= SIDE_TOP) then
                    errmsg = "&FILM ncells is supported in mode 'film' with side = 'top' " // &
                             '(the partial principal layer sits at the vacuum side)'; return
                end if
                if (inp%nlayers /= 0) then
                    errmsg = '&FILM: give either nlayers (principal layers) or ncells (cells), not both'; return
                end if
            else if (inp%nlayers < 1 .or. inp%nlayers > MAX_GLAYERS) then
                errmsg = '&FILM nlayers must be in 1..512 (or use ncells)'; return
            end if
        end if
        if (trim(inp%contact_model) == 'sk') then
            if (len_trim(inp%sk_file) == 0) then
                errmsg = "&INTERFACE model = 'sk' needs sk_file"; return
            end if
            if (count(len_trim(inp%sk_basis) > 0) == 0) then
                errmsg = "&INTERFACE model = 'sk' needs sk_basis (e.g. 'Te s3p2d2', 'Sb s3p2d2')"; return
            end if
            if (inp%sk_cutoff_ang <= 0.0_dp) then
                errmsg = '&INTERFACE sk_cutoff_ang must be > 0'; return
            end if
        end if
        if (inp%nscan > 0) then
            if (index(' ' // SCAN_NAMES // ' ', ' ' // trim(inp%scan_parameter) // ' ') == 0 .or. &
                len_trim(inp%scan_parameter) == 0) then
                errmsg = '&SCAN parameter must be one of: ' // SCAN_NAMES; return
            end if
            if (trim(inp%mode) == 'surface') then
                errmsg = "&SCAN needs a film (not available in mode 'surface')"; return
            end if
            select case (trim(inp%scan_parameter))
            case ('lambda_s')
                if (trim(inp%contact_model) /= 'scaled') then
                    errmsg = "&SCAN lambda_s needs &INTERFACE model = 'scaled'"; return
                end if
            case ('sk_scale_sigma', 'sk_scale_pi', 'sk_scale_delta', 'sk_beta')
                if (trim(inp%contact_model) /= 'sk') then
                    errmsg = '&SCAN ' // trim(inp%scan_parameter) // " needs &INTERFACE model = 'sk'"; return
                end if
            case ('lambda_h')
                if (trim(inp%contact_model) /= 'scaled' .and. trim(inp%contact_model) /= 'sk') then
                    errmsg = "&SCAN lambda_h needs &INTERFACE model = 'scaled' or 'sk'"; return
                end if
            end select
        end if
        if (inp%nqz < 0) then
            errmsg = '&COHERENT nqz must be >= 0'; return
        end if
        if (inp%escape_length_ang < 0.0_dp) then
            errmsg = '&COHERENT escape_length_ang must be >= 0'; return
        end if
        if (inp%kmode == 'path') then
            if (inp%nvertex < 2 .or. inp%nvertex > MAX_GKVERTEX) then
                errmsg = '&KPATH nvertex must be in 2..64'; return
            end if
            if (any(inp%npoints(1:inp%nvertex - 1) < 1)) then
                errmsg = '&KPATH npoints must be >= 1 for every segment'; return
            end if
            if (any(abs(inp%vertex(inp%stacking_axis, 1:inp%nvertex)) > 0.0_dp)) then
                errmsg = '&KPATH vertices must have zero component along the stacking axis'; return
            end if
        else
            if (inp%grid_n1 < 2 .or. inp%grid_n2 < 2) then
                errmsg = '&KGRID n1 and n2 must be >= 2'; return
            end if
            if (abs(inp%grid_origin(inp%stacking_axis)) > 0.0_dp .or. &
                abs(inp%grid_span1(inp%stacking_axis)) > 0.0_dp .or. &
                abs(inp%grid_span2(inp%stacking_axis)) > 0.0_dp) then
                errmsg = '&KGRID origin/span must have zero component along the stacking axis'; return
            end if
        end if
        ierr = 0
    end subroutine validate_input

    ! ==================================================================
    ! Materiais e empilhamento
    ! ==================================================================

    subroutine load_material(hfile, sfile, stfile, axis, pl_cells, shift, tag, mat, ierr, errmsg)
        character(len=*), intent(in) :: hfile, sfile, stfile, tag
        integer,          intent(in) :: axis, pl_cells
        real(dp),         intent(in) :: shift
        type(material_t), intent(out) :: mat
        integer,          intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        type(HSDataFormat) :: hs
        character(len=512) :: msg
        logical :: ex

        ierr = 0
        errmsg = ''
        inquire(file=trim(hfile), exist=ex)
        if (.not. ex) then
            ierr = 1; errmsg = "H file not found: '" // trim(hfile) // "'"; return
        end if
        call mat%hh%ReadFile(trim(hfile))
        if (len_trim(sfile) > 0) then
            inquire(file=trim(sfile), exist=ex)
            if (.not. ex) then
                ierr = 1; errmsg = "S file not found: '" // trim(sfile) // "'"; return
            end if
            call hs%ReadFile(trim(sfile))
            call build_pl_blocks(mat%hh, axis, pl_cells, mat%pl, ierr, errmsg, hs=hs)
        else
            call build_pl_blocks(mat%hh, axis, pl_cells, mat%pl, ierr, errmsg)
        end if
        if (ierr /= 0) then
            errmsg = trim(tag) // ': ' // errmsg; return
        end if
        deallocate(mat%hh%Data)   ! só os metadados são usados depois
        call read_structure(trim(stfile), mat%st, ierr, msg)
        if (ierr /= 0) then
            errmsg = trim(tag) // ' structure: ' // trim(msg); return
        end if
        call attach_pl_geometry(mat%pl, mat%st, ierr, errmsg)
        if (ierr /= 0) then
            errmsg = trim(tag) // ': ' // errmsg; return
        end if
        mat%shift = shift
        mat%tag = tag(1:1)
    end subroutine load_material

    !> Resolve a sequência de camadas, os meios e a estrutura da região central.
    subroutine build_stack(inp, mats, stk, ierr, errmsg)
        type(green_input_t), intent(in)  :: inp
        type(material_t),    intent(in)  :: mats(:)
        type(green_stack_t), intent(out) :: stk
        integer,             intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: nA, nB, l, nat, a, at, m, ia, P, k
        real(dp) :: o(3), anorm(3), top
        integer, allocatable :: seq(:)
        logical :: top_side

        ierr = 0
        errmsg = ''
        nA = inp%sub_layers
        nB = inp%nlayers
        if (inp%film_ncells > 0) nB = (inp%film_ncells + mats(2)%pl%ncell - 1) / mats(2)%pl%ncell
        top_side = side_code(inp%side) == SIDE_TOP
        select case (trim(inp%mode))
        case ('surface')
            seq = [(1, l = 1, nA)]
            if (top_side) then
                stk%bot_lead = 1; stk%top_lead = 0
            else
                stk%bot_lead = 0; stk%top_lead = 1
            end if
        case ('film')
            if (top_side) then
                seq = [[(1, l = 1, nA)], [(2, l = 1, nB)]]
                stk%bot_lead = 1; stk%top_lead = 0
            else
                seq = [[(2, l = 1, nB)], [(1, l = 1, nA)]]
                stk%bot_lead = 0; stk%top_lead = 1
            end if
        case ('junction')
            if (top_side) then
                seq = [[(1, l = 1, nA)], [(2, l = 1, nB)]]
                stk%bot_lead = 1; stk%top_lead = 2
            else
                seq = [[(2, l = 1, nB)], [(1, l = 1, nA)]]
                stk%bot_lead = 2; stk%top_lead = 1
            end if
        case ('sandwich')
            seq = [[(1, l = 1, nA)], [(2, l = 1, nB)], [(1, l = 1, nA)]]
            stk%bot_lead = 1; stk%top_lead = 1
        end select
        stk%nlayer = size(seq)
        stk%mat = seq
        allocate(stk%ncell(stk%nlayer))
        do l = 1, stk%nlayer
            stk%ncell(l) = mats(stk%mat(l))%pl%ncell
        end do
        if (inp%film_ncells > 0) then
            ! último bloco do filme (topo, lado do vácuo): só as primeiras r células da PL
            stk%ncell(stk%nlayer) = inp%film_ncells - (nB - 1) * mats(2)%pl%ncell
        end if
        allocate(stk%film_index(stk%nlayer), stk%v(stk%nlayer), stk%origin(3, stk%nlayer))
        k = 0
        do l = 1, stk%nlayer
            if (stk%mat(l) == 2) then
                k = k + 1
                stk%film_index(l) = k
                stk%v(l) = mats(2)%shift + inp%layer_potential_ev(k)
            else
                stk%film_index(l) = 0
                stk%v(l) = mats(1)%shift
            end if
        end do

        ! Origens das camadas: O_1 = 0; O_{l+1} = O_l + P*a_s(mat_l) [+ offset na troca de material]
        o = 0.0_dp
        do l = 1, stk%nlayer
            stk%origin(:, l) = o
            o = o + real(stk%ncell(l), dp) * mats(stk%mat(l))%st%lattice(inp%stacking_axis, :)
            if (l < stk%nlayer) then
                if (stk%mat(l + 1) /= stk%mat(l)) then
                    if (stk%mat(l) == 1) then
                        o = o + inp%interface_offset
                    else
                        o = o - inp%interface_offset
                    end if
                end if
            end if
        end do

        ! Estrutura da região central (átomos na ordem dos orbitais)
        nat = 0
        do l = 1, stk%nlayer
            nat = nat + stk%ncell(l) * mats(stk%mat(l))%pl%nat_cell
        end do
        stk%st%nat = nat
        allocate(stk%st%cart(3, nat), stk%st%label(nat), stk%atom_orbs(nat), stk%atom_layer(nat))
        ia = 0
        do l = 1, stk%nlayer
            m = stk%mat(l)
            P = stk%ncell(l)
            do a = 0, P - 1
                do at = 1, mats(m)%pl%nat_cell
                    ia = ia + 1
                    stk%st%cart(:, ia) = stk%origin(:, l) + mats(m)%st%cart(:, at) + &
                                         real(a, dp) * mats(m)%st%lattice(inp%stacking_axis, :)
                    stk%st%label(ia) = mats(m)%st%label(at)
                    stk%atom_orbs(ia) = mats(m)%pl%orbs_per_atom(at)
                    stk%atom_layer(ia) = l
                end do
            end do
        end do
        ! Rede: a1, a2 do substrato; a_s = espessura total + vácuo ao longo da normal
        stk%st%lattice = mats(1)%st%lattice
        anorm = mats(1)%pl%normal
        top = maxval(matmul(anorm, stk%st%cart)) - minval(matmul(anorm, stk%st%cart))
        stk%st%lattice(inp%stacking_axis, :) = o - dot_product(o, anorm) * anorm + &
                                               (max(dot_product(o, anorm), top) + inp%vacuum_ang) * anorm
        call species_table(mats, stk, ierr, errmsg)
    end subroutine build_stack

    subroutine species_table(mats, stk, ierr, errmsg)
        type(material_t),    intent(in)    :: mats(:)
        type(green_stack_t), intent(inout) :: stk
        integer,             intent(out)   :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=LABEL_LEN) :: names(stk%st%nat)
        integer :: norb(stk%st%nat), ns, i, j
        logical :: found

        ierr = 0
        errmsg = ''
        if (size(mats) < 1) return
        ns = 0
        do i = 1, stk%st%nat
            found = .false.
            do j = 1, ns
                if (names(j) == stk%st%label(i)) then
                    found = .true.
                    if (norb(j) /= stk%atom_orbs(i)) then
                        ierr = 1
                        errmsg = "species '" // trim(stk%st%label(i)) // &
                                 "' has different orbital counts in substrate and film"
                        return
                    end if
                    exit
                end if
            end do
            if (.not. found) then
                ns = ns + 1
                names(ns) = stk%st%label(i)
                norb(ns) = stk%atom_orbs(i)
            end if
        end do
        stk%st%has_species = .true.
        stk%st%nspecies = ns
        stk%st%species_name = names(1:ns)
        stk%st%species_norb = norb(1:ns)
    end subroutine species_table

    !> Grava a estrutura da região central no formato de structure_io.
    subroutine write_central_structure(filename, stk, ierr, errmsg)
        character(len=*),    intent(in) :: filename
        type(green_stack_t), intent(in) :: stk
        integer,             intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: u, ios, i
        ierr = 0
        errmsg = ''
        open(newunit=u, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = 'cannot write central structure: ' // trim(filename); return
        end if
        write(u, '(A)') '! Central region written by surfgreen (layers bottom -> top).'
        write(u, '(A)') '&STRUCTURE'
        do i = 1, 3
            write(u, '(A,I0,A,3(F16.8,:,","))') '    LatticeVectors(', i, ',:) = ', stk%st%lattice(i, :)
        end do
        write(u, '(A)') '/'
        write(u, '(A)') ''
        write(u, '(A)') '&SPECIES'
        write(u, '(I0)') stk%st%nspecies
        do i = 1, stk%st%nspecies
            write(u, '(A,1X,I0)') trim(stk%st%species_name(i)), stk%st%species_norb(i)
        end do
        write(u, '(A)') '/'
        write(u, '(A)') ''
        write(u, '(A)') '&ATOMS'
        write(u, '(I0)') stk%st%nat
        do i = 1, stk%st%nat
            write(u, '(A,3(1X,F16.8),A,I0)') stk%st%label(i), stk%st%cart(:, i), '   ! layer ', stk%atom_layer(i)
        end do
        write(u, '(A)', iostat=ios) '/'
        close(u)
        if (ios /= 0) then
            ierr = 2; errmsg = 'error writing central structure: ' // trim(filename)
        end if
    end subroutine write_central_structure

    ! ==================================================================
    ! Execução
    ! ==================================================================

    !> Executa o cálculo. Se structure_only, só grava a estrutura central.
    !> nonconverged: número de pontos (k,E) em que a decimação não convergiu.
    !> Com &SCAN, gd é o primeiro valor e gscan (se presente) recebe todos.
    subroutine run_green(inp, gd, ierr, errmsg, structure_only, nonconverged, nthreads, gscan)
        type(green_input_t), intent(in)  :: inp
        type(green_data_t),  intent(out) :: gd
        integer,             intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        logical, optional,   intent(in)  :: structure_only
        integer, optional,   intent(out) :: nonconverged
        integer, optional,   intent(out) :: nthreads
        type(green_data_t), allocatable, optional, intent(out) :: gscan(:)

        type(material_t), allocatable :: mats(:)
        type(green_stack_t) :: stk
        type(contact_t) :: contact
        type(structure_t) :: stc
        type(region_set_t) :: rs
        type(orbital_groups_t) :: groups
        type(kpath_t) :: path
        type(kgrid_t) :: kg
        character(len=512) :: msg
        character(len=GPATH_LEN) :: central_file
        integer :: nmat, nk, ne, nq, nreg, ik, ie, ip, l, g, ntot, cohreg, nbad
        integer :: perr_ip, nthr
        integer, allocatable :: list_start(:), list_atoms(:), orb_layer(:)
        real(dp), allocatable :: depth(:), rpos(:, :), coh_mask(:)
        character(len=:), allocatable :: perr_msg
        real(dp) :: recip(3, 3), ntop
        ! varredura e contato SK
        integer :: nscan, is
        type(green_input_t) :: inp_s
        type(green_stack_t), allocatable :: stks(:)
        type(contact_t), allocatable :: contacts(:)
        type(sk_model_t) :: skm, skm_s
        type(sk_geom_t) :: skg(2)
        type(sk_pairs_t), allocatable :: pab(:), pba(:)
        real(dp), allocatable :: shiftv(:, :), depth_s(:, :), rpos_s(:, :, :)
        type(green_data_t), allocatable :: gall(:)
        real(dp) :: base_film_shift

        ierr = 0
        errmsg = ''
        if (present(nonconverged)) nonconverged = 0
        nthr = 1
        !$ nthr = omp_get_max_threads()
        if (present(nthreads)) nthreads = nthr

        ! --- materiais ---
        nmat = merge(1, 2, trim(inp%mode) == 'surface')
        allocate(mats(nmat))
        call load_material(inp%sub_h, inp%sub_s, inp%sub_structure, inp%stacking_axis, &
                           inp%sub_pl_cells, inp%sub_shift_ev, 'A', mats(1), ierr, errmsg)
        if (ierr /= 0) then
            errmsg = 'substrate ' // errmsg; return
        end if
        if (nmat == 2) then
            call load_material(inp%film_h, inp%film_s, inp%film_structure, inp%stacking_axis, &
                               inp%film_pl_cells, inp%film_shift_ev, 'B', mats(2), ierr, errmsg)
            if (ierr /= 0) then
                errmsg = 'film ' // errmsg; return
            end if
            do l = 1, 2
                g = mats(1)%pl%inplane(l)
                if (norm2(mats(1)%st%lattice(g, :) - mats(2)%st%lattice(g, :)) > inp%lattice_tol_ang) then
                    ierr = 2
                    errmsg = 'film and substrate in-plane lattice vectors differ by more than lattice_tol_ang; ' // &
                             'the film structure must be strained to the substrate'
                    return
                end if
            end do
            call contact_validate(inp%contact_model, inp%contact_reference, inp%lambda_h, inp%lambda_s, &
                                  contact, ierr, errmsg)
            if (ierr /= 0) return
            if (contact%model == CONTACT_SK) then
                skm%cutoff = inp%sk_cutoff_ang
                call read_sk_file(trim(inp%sk_file), skm, ierr, errmsg)
                if (ierr /= 0) return
                call sk_set_bases(skm, inp%sk_basis, ierr, errmsg)
                if (ierr /= 0) return
                do l = 1, 2
                    call make_sk_geom(skm, mats(l)%pl%norb, mats(l)%pl%tau, mats(l)%pl%orb_atom, &
                                      mats(l)%pl%orb_cell, mats(l)%st%label, mats(l)%pl%orbs_per_atom, &
                                      skg(l), ierr, errmsg)
                    if (ierr /= 0) then
                        errmsg = merge('substrate: ', 'film:      ', l == 1) // errmsg; return
                    end if
                end do
            end if
        end if

        ! --- empilhamento (um por valor de varredura) e estrutura central ---
        nscan = max(1, inp%nscan)
        allocate(stks(nscan), contacts(nscan), pab(nscan), pba(nscan), shiftv(nmat, nscan))
        base_film_shift = 0.0_dp
        if (nmat == 2) base_film_shift = mats(2)%shift
        do is = 1, nscan
            inp_s = scan_input(inp, is, mats(1)%pl%normal)
            if (nmat == 2) mats(2)%shift = inp_s%film_shift_ev
            call build_stack(inp_s, mats, stks(is), ierr, errmsg)
            if (nmat == 2) mats(2)%shift = base_film_shift
            if (ierr /= 0) return
            shiftv(1, is) = mats(1)%shift
            if (nmat == 2) then
                shiftv(2, is) = inp_s%film_shift_ev
                contacts(is) = contact
                contacts(is)%lambda_h = inp_s%lambda_h
                contacts(is)%lambda_s = inp_s%lambda_s
                if (contact%model == CONTACT_SK) then
                    skm_s = skm
                    skm_s%lambda = inp_s%lambda_h
                    skm_s%scale = inp_s%sk_scale
                    skm_s%beta_override = inp_s%sk_beta
                    call build_sk_pairs(skm_s, skg(1), skg(2), mats(1)%pl%pl_shift + inp_s%interface_offset, &
                                        mats(1)%pl%lattice(mats(1)%pl%inplane(1), :), &
                                        mats(1)%pl%lattice(mats(1)%pl%inplane(2), :), pab(is), ierr, errmsg)
                    if (ierr /= 0) return
                    call build_sk_pairs(skm_s, skg(2), skg(1), mats(2)%pl%pl_shift - inp_s%interface_offset, &
                                        mats(1)%pl%lattice(mats(1)%pl%inplane(1), :), &
                                        mats(1)%pl%lattice(mats(1)%pl%inplane(2), :), pba(is), ierr, errmsg)
                    if (ierr /= 0) return
                end if
            end if
        end do
        stk = stks(1)
        central_file = trim(inp%output_prefix) // '_central.nml'
        call write_central_structure(central_file, stk, ierr, errmsg)
        if (ierr /= 0) return
        if (present(structure_only)) then
            if (structure_only) return
        end if
        call read_structure(trim(central_file), stc, ierr, msg)
        if (ierr /= 0) then
            errmsg = 'central structure: ' // trim(msg); return
        end if

        ! --- orbitais da região central (posições por valor de varredura) ---
        ntot = sum(stk%atom_orbs)
        allocate(orb_layer(ntot), depth(ntot), rpos(3, ntot), depth_s(ntot, nscan), rpos_s(3, ntot, nscan))
        block
            integer :: ia, o0, j
            do is = 1, nscan
                ntop = maxval(matmul(mats(1)%pl%normal, stks(is)%st%cart))
                o0 = 0
                do ia = 1, stks(is)%st%nat
                    do j = 1, stk%atom_orbs(ia)
                        orb_layer(o0 + j) = stk%atom_layer(ia)
                        rpos_s(:, o0 + j, is) = stks(is)%st%cart(:, ia)
                        depth_s(o0 + j, is) = ntop - dot_product(mats(1)%pl%normal, stks(is)%st%cart(:, ia))
                    end do
                    o0 = o0 + stk%atom_orbs(ia)
                end do
            end do
            depth = depth_s(:, 1)
            rpos = rpos_s(:, :, 1)
        end block

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
        else
            ! padrão: uma região por camada
            nreg = stk%nlayer
            allocate(list_start(nreg + 1), list_atoms(stc%nat))
            list_start(1) = 1
            g = 0
            do l = 1, nreg
                do ik = 1, stc%nat
                    if (stk%atom_layer(ik) == l) then
                        g = g + 1
                        list_atoms(g) = ik
                    end if
                end do
                list_start(l + 1) = g + 1
            end do
        end if
        call build_orbital_groups(stk%atom_orbs, nreg, list_start, list_atoms, groups)

        ! --- pontos k ---
        recip = reciprocal_lattice(mats(1)%st)
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
            allocate(gd%vertex_k(0), gd%label(0))
        end if
        gd%nk = nk
        allocate(gd%kcart(3, nk), gd%kdist(nk))
        do ik = 1, nk
            gd%kcart(:, ik) = kpar_cartesian(mats(1)%pl, kpar_of(gd%kfrac(:, ik)))
        end do
        gd%kdist = 0.0_dp
        if (gd%kmode == GREEN_KPATH) then
            do ik = 2, nk
                gd%kdist(ik) = gd%kdist(ik - 1) + norm2(gd%kcart(:, ik) - gd%kcart(:, ik - 1))
            end do
        end if

        ! --- energias, q_z ---
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
        nq = inp%nqz
        gd%nq = nq
        allocate(gd%qz(nq))
        do ie = 1, nq
            if (nq == 1) then
                gd%qz(ie) = inp%qz_min
            else
                gd%qz(ie) = inp%qz_min + real(ie - 1, dp) * (inp%qz_max - inp%qz_min) / real(nq - 1, dp)
            end if
        end do
        allocate(coh_mask(ntot))
        if (nq > 0) then
            coh_mask = 1.0_dp
            if (len_trim(inp%coherent_region) > 0) then
                cohreg = 0
                if (len_trim(inp%regions_file) > 0) then
                    do g = 1, nreg
                        if (trim(rs%regions(g)%name) == trim(inp%coherent_region)) cohreg = g
                    end do
                end if
                if (cohreg == 0) then
                    ierr = 2
                    errmsg = "&COHERENT region '" // trim(inp%coherent_region) // "' not found in regions_file"
                    return
                end if
                coh_mask = 0.0_dp
                coh_mask(groups%orb(groups%start(cohreg):groups%start(cohreg + 1) - 1)) = 1.0_dp
            end if
            ! a atenuação exp(-d/2L) entra em one_point (depende da geometria de cada valor de varredura)
        end if

        ! --- saídas ---
        gd%nreg = nreg
        allocate(gd%region(nreg))
        if (len_trim(inp%regions_file) > 0) then
            do g = 1, nreg
                gd%region(g) = rs%regions(g)%name
            end do
        else
            do l = 1, nreg
                write(gd%region(l), '(A,I3.3,A,A)') 'L', l, '_', merge('sub ', 'film', stk%mat(l) == 1)
            end do
        end if
        allocate(gd%total(ne, nk), gd%amap(ne, nk, nreg))
        allocate(gd%niter(ne, nk), gd%residual(ne, nk), gd%min_neg_im_gii(ne, nk), gd%smin(nk))
        gd%has_coherent = nq > 0
        if (gd%has_coherent) then
            allocate(gd%icoh(ne, nk, nq), gd%iinc(ne, nk, nq))
        else
            allocate(gd%icoh(0, 0, 0), gd%iinc(0, 0, 0))
        end if

        ! --- uma saída por valor de varredura ---
        allocate(gall(nscan))
        do is = 1, nscan
            gall(is) = gd
        end do

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
        if (present(nonconverged)) nonconverged = nbad

        do is = 1, nscan
            call fill_meta(gall(is), is)
        end do
        gd = gall(1)
        if (present(gscan)) call move_alloc(gall, gscan)

    contains

        pure function kpar_of(kf) result(kp)
            real(dp), intent(in) :: kf(3)
            real(dp) :: kp(2)
            kp = [kf(mats(1)%pl%inplane(1)), kf(mats(1)%pl%inplane(2))]
        end function kpar_of

        subroutine record_error(ipt, msg_in)
            integer, intent(in) :: ipt
            character(len=*), intent(in) :: msg_in
            !$omp critical (green_error)
            if (ipt < perr_ip) then
                perr_ip = ipt
                perr_msg = msg_in
            end if
            !$omp end critical (green_error)
        end subroutine record_error

        !> Bloco de acoplamento H_{l,l+1}, S_{l,l+1} entre materiais ml (abaixo) e mu (acima).
        subroutine coupling(ml, mu, is, kp, h, s, hc, sc, lerr, lmsg)
            integer, intent(in) :: ml, mu, is
            real(dp), intent(in) :: kp(2)
            complex(dp), intent(in) :: h(:, :, -1:, :), s(:, :, -1:, :)
            complex(dp), allocatable, intent(out) :: hc(:, :), sc(:, :)
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            integer :: rm
            lerr = 0
            lmsg = ''
            if (ml == mu) then
                hc = h(1:mats(ml)%pl%norb, 1:mats(ml)%pl%norb, 1, ml)
                sc = s(1:mats(ml)%pl%norb, 1:mats(ml)%pl%norb, 1, ml)
            else if (contact%model == CONTACT_SK) then
                allocate(hc(mats(ml)%pl%norb, mats(mu)%pl%norb), sc(mats(ml)%pl%norb, mats(mu)%pl%norb))
                sc = (0.0_dp, 0.0_dp)
                if (ml == 1) then
                    call sk_contact_at_k(pab(is), kp, hc)
                else
                    call sk_contact_at_k(pba(is), kp, hc)
                end if
            else
                rm = merge(1, 2, contact%reference == REF_SUBSTRATE)
                call build_contact(contacts(is), mats(ml)%pl%norb, mats(mu)%pl%norb, &
                                   h(1:mats(rm)%pl%norb, 1:mats(rm)%pl%norb, 1, rm), &
                                   s(1:mats(rm)%pl%norb, 1:mats(rm)%pl%norb, 1, rm), hc, sc, lerr, lmsg)
            end if
        end subroutine coupling

        subroutine one_point(ik, ie, nbad)
            integer, intent(in) :: ik, ie
            integer, intent(inout) :: nbad
            complex(dp), allocatable :: h(:, :, :, :), s(:, :, :, :), gbot(:, :), gtop(:, :)
            complex(dp), allocatable :: hc(:, :), sc(:, :), w(:, :)
            type(central_t) :: sys
            type(rgf_t) :: res
            real(dp), allocatable :: d(:), rd(:), ic(:), ii(:)
            integer, allocatable :: bsize(:)
            integer :: lerr, nmax, m, it_b, it_t, iq, i, l, isc
            character(len=:), allocatable :: lmsg
            real(dp) :: res_b, res_t, smin, kc(3), sh_b, sh_t
            complex(dp) :: z
            real(dp) :: kp(2)
            logical :: have_b, have_t

            kp = kpar_of(gd%kfrac(:, ik))
            z = cmplx(gd%energy(ie), inp%eta_ev, dp)
            nmax = maxval([(mats(m)%pl%norb, m = 1, nmat)])
            allocate(h(nmax, nmax, -1:1, nmat), s(nmax, nmax, -1:1, nmat))
            h = (0.0_dp, 0.0_dp)
            s = (0.0_dp, 0.0_dp)
            do m = 1, nmat
                call pl_blocks_at_k(mats(m)%pl, kp, h(1:mats(m)%pl%norb, 1:mats(m)%pl%norb, :, m), &
                                    s(1:mats(m)%pl%norb, 1:mats(m)%pl%norb, :, m))
            end do
            allocate(bsize(stk%nlayer))
            do l = 1, stk%nlayer
                bsize(l) = stk%ncell(l) * mats(stk%mat(l))%pl%nbasis
            end do

            it_b = 0; it_t = 0; res_b = 0.0_dp; res_t = 0.0_dp
            have_b = .false.; have_t = .false.
            sh_b = 0.0_dp; sh_t = 0.0_dp
            do isc = 1, nscan
                ! meios semi-infinitos: decimação só quando o deslocamento do meio muda
                if (stk%bot_lead > 0) then
                    if (.not. have_b .or. shiftv(stk%bot_lead, isc) /= sh_b) then
                        sh_b = shiftv(stk%bot_lead, isc)
                        call lead_gf_impl(stk%bot_lead, sh_b, SIDE_TOP, h(:, :, :, stk%bot_lead), &
                                          s(:, :, :, stk%bot_lead), z, gbot, it_b, res_b, lerr, lmsg)
                        if (lerr == 3) then
                            nbad = nbad + 1
                            it_b = -it_b
                        else if (lerr /= 0) then
                            call record_error(ip_of(ik, ie), lmsg); return
                        end if
                        have_b = .true.
                    end if
                end if
                if (stk%top_lead > 0) then
                    if (.not. have_t .or. shiftv(stk%top_lead, isc) /= sh_t) then
                        sh_t = shiftv(stk%top_lead, isc)
                        call lead_gf_impl(stk%top_lead, sh_t, SIDE_BOTTOM, h(:, :, :, stk%top_lead), &
                                          s(:, :, :, stk%top_lead), z, gtop, it_t, res_t, lerr, lmsg)
                        if (lerr == 3) then
                            nbad = nbad + 1
                            it_t = -it_t
                        else if (lerr /= 0) then
                            call record_error(ip_of(ik, ie), lmsg); return
                        end if
                        have_t = .true.
                    end if
                end if

                ! região central
                call central_init(sys, bsize)
                do l = 1, stk%nlayer
                    m = stk%mat(l)
                    if (l < stk%nlayer) then
                        call coupling(m, stk%mat(l + 1), isc, kp, h, s, hc, sc, lerr, lmsg)
                        if (lerr /= 0) then
                            call record_error(ip_of(ik, ie), lmsg); return
                        end if
                        ! bloco parcial: primeiras células da PL seguinte (ordem ao longo de a_s)
                        call assemble_central(sys, l, z, h(1:bsize(l), 1:bsize(l), 0, m), &
                                              s(1:bsize(l), 1:bsize(l), 0, m), stks(isc)%v(l), &
                                              hu=hc(1:bsize(l), 1:bsize(l + 1)), su=sc(1:bsize(l), 1:bsize(l + 1)), &
                                              vnext=stks(isc)%v(l + 1))
                    else
                        call assemble_central(sys, l, z, h(1:bsize(l), 1:bsize(l), 0, m), &
                                              s(1:bsize(l), 1:bsize(l), 0, m), stks(isc)%v(l))
                    end if
                end do
                if (stk%bot_lead > 0) then
                    call coupling(stk%bot_lead, stk%mat(1), isc, kp, h, s, hc, sc, lerr, lmsg)
                    if (lerr /= 0) then
                        call record_error(ip_of(ik, ie), lmsg); return
                    end if
                    call attach_lead(sys, 'bot', z, hc(:, 1:bsize(1)), sc(:, 1:bsize(1)), stks(isc)%v(1), sh_b, gbot)
                end if
                if (stk%top_lead > 0) then
                    call coupling(stk%mat(stk%nlayer), stk%top_lead, isc, kp, h, s, hc, sc, lerr, lmsg)
                    if (lerr /= 0) then
                        call record_error(ip_of(ik, ie), lmsg); return
                    end if
                    call attach_lead(sys, 'top', z, hc, sc, stks(isc)%v(stk%nlayer), sh_t, gtop)
                end if

                call rgf_solve(sys, res, lerr, lmsg)
                if (lerr /= 0) then
                    call record_error(ip_of(ik, ie), lmsg); return
                end if

                if (allocated(d)) deallocate(d, rd)
                allocate(d(ntot), rd(nreg))
                if (trim(inp%projection) == 'lowdin') then
                    call orbital_dos_lowdin(sys, d, smin, lerr, lmsg)
                    if (lerr /= 0) then
                        call record_error(ip_of(ik, ie), 'Lowdin projection: ' // lmsg); return
                    end if
                else
                    call orbital_dos_mulliken(sys, res, d)
                end if
                call region_dos(d, groups, rd)
                gall(isc)%total(ie, ik) = sum(d)
                gall(isc)%amap(ie, ik, :) = rd
                gall(isc)%niter(ie, ik) = merge(min(it_b, it_t), max(it_b, it_t), it_b < 0 .or. it_t < 0)
                gall(isc)%residual(ie, ik) = max(res_b, res_t)
                gall(isc)%min_neg_im_gii(ie, ik) = min_neg_im_gdiag(res)
                if (ie == 1) then
                    call central_overlap_min(sys, smin, lerr, lmsg)
                    gall(isc)%smin(ik) = smin
                end if

                if (nq > 0) then
                    if (allocated(w)) deallocate(w, ic, ii)
                    allocate(w(ntot, nq), ic(nq), ii(nq))
                    kc = gd%kcart(:, ik)
                    do iq = 1, nq
                        do i = 1, ntot
                            ! gauge atômico: w' = U w, U_ii = exp(i k_par . r_i)
                            w(i, iq) = coh_mask(i) * exp(cmplx(0.0_dp, gd%qz(iq) * depth_s(i, isc) + &
                                                               dot_product(kc, rpos_s(:, i, isc)), dp))
                            if (inp%escape_length_ang > 0.0_dp) &
                                w(i, iq) = w(i, iq) * exp(-depth_s(i, isc) / (2.0_dp * inp%escape_length_ang))
                        end do
                    end do
                    call coherent_map(sys, res, w, ic, ii, lerr, lmsg)
                    if (lerr /= 0) then
                        call record_error(ip_of(ik, ie), lmsg); return
                    end if
                    gall(isc)%icoh(ie, ik, :) = ic
                    gall(isc)%iinc(ie, ik, :) = ii
                end if
            end do
        end subroutine one_point

        pure integer function ip_of(ik, ie)
            integer, intent(in) :: ik, ie
            ip_of = (ik - 1) * ne + ie
        end function ip_of

        !> Green de superfície do meio do material mm, no lado indicado.
        subroutine lead_gf_impl(mm, shift, side, hh, ss, z, gs, it, rr, lerr, lmsg)
            integer, intent(in) :: mm, side
            real(dp), intent(in) :: shift
            complex(dp), intent(in) :: hh(:, :, -1:), ss(:, :, -1:)
            complex(dp), intent(in) :: z
            complex(dp), allocatable, intent(out) :: gs(:, :)
            integer, intent(out) :: it, lerr
            real(dp), intent(out) :: rr
            character(len=:), allocatable, intent(out) :: lmsg
            complex(dp), allocatable :: k00(:, :), k01(:, :), k10(:, :)
            integer :: n
            n = mats(mm)%pl%norb
            k00 = (z - shift) * ss(1:n, 1:n, 0) - hh(1:n, 1:n, 0)
            k01 = (z - shift) * ss(1:n, 1:n, 1) - hh(1:n, 1:n, 1)
            k10 = (z - shift) * ss(1:n, 1:n, -1) - hh(1:n, 1:n, -1)
            allocate(gs(n, n))
            call surface_gf(k00, k01, k10, side, inp%decim_tol, inp%decim_maxit, gs, it, rr, lerr, lmsg)
        end subroutine lead_gf_impl

        subroutine fill_meta(gx, is)
            type(green_data_t), intent(inout) :: gx
            integer, intent(in) :: is
            character(len=64) :: buf
            type(green_input_t) :: ins
            ins = scan_input(inp, is, mats(1)%pl%normal)
            call meta_add(gx, 'program', 'surfgreen')
            call meta_add(gx, 'format', 'green.bin v1')
            call meta_add(gx, 'mode', trim(inp%mode))
            call meta_add(gx, 'side', trim(inp%side))
            write(buf, '(ES16.8)') inp%eta_ev
            call meta_add(gx, 'eta_ev', adjustl(buf))
            write(buf, '(ES12.4)') inp%decim_tol
            call meta_add(gx, 'decim_tol', adjustl(buf))
            write(buf, '(I0)') inp%decim_maxit
            call meta_add(gx, 'decim_maxit', adjustl(buf))
            call meta_add(gx, 'projection', trim(inp%projection))
            call meta_add(gx, 'map_gauge', 'lattice (region DOS is gauge invariant)')
            if (nq > 0) then
                call meta_add(gx, 'coherent_gauge', 'atomic: w_i = a_i exp(i q_z d_i) exp(-d_i/2L) exp(i k.r_i)')
                write(buf, '(ES16.8)') inp%escape_length_ang
                call meta_add(gx, 'escape_length_ang', adjustl(buf))
                call meta_add(gx, 'coherent_region', trim(inp%coherent_region))
            end if
            call meta_add(gx, 'substrate_h', trim(inp%sub_h))
            call meta_add(gx, 'substrate_s', trim(inp%sub_s))
            call meta_add(gx, 'substrate_structure_fingerprint', structure_fingerprint(mats(1)%st))
            write(buf, '(ES16.8)') inp%sub_shift_ev
            call meta_add(gx, 'substrate_shift_ev', adjustl(buf))
            write(buf, '(I0,A,I0)') mats(1)%pl%ncell, ' x ', mats(1)%pl%nbasis
            call meta_add(gx, 'substrate_pl', trim(buf) // ' (cells x orbitals)')
            if (nmat == 2) then
                call meta_add(gx, 'film_h', trim(inp%film_h))
                call meta_add(gx, 'film_s', trim(inp%film_s))
                call meta_add(gx, 'film_structure_fingerprint', structure_fingerprint(mats(2)%st))
                write(buf, '(ES16.8)') ins%film_shift_ev
                call meta_add(gx, 'film_shift_ev', adjustl(buf))
                write(buf, '(3F12.6)') ins%interface_offset
                call meta_add(gx, 'interface_offset_ang', adjustl(buf))
                write(buf, '(I0,A,I0)') mats(2)%pl%ncell, ' x ', mats(2)%pl%nbasis
                call meta_add(gx, 'film_pl', trim(buf) // ' (cells x orbitals)')
                write(buf, '(I0)') sum(stk%ncell, mask=stk%mat == 2)
                call meta_add(gx, 'film_cells', adjustl(buf))
                call meta_add(gx, 'contact_model', trim(inp%contact_model))
                call meta_add(gx, 'contact_reference', trim(inp%contact_reference))
                write(buf, '(2ES16.8)') ins%lambda_h, ins%lambda_s
                call meta_add(gx, 'lambda_h_s', adjustl(buf))
                if (contact%model == CONTACT_SK) then
                    call meta_add(gx, 'sk_file', trim(inp%sk_file))
                    write(buf, '(F10.4)') inp%sk_cutoff_ang
                    call meta_add(gx, 'sk_cutoff_ang', adjustl(buf))
                    write(buf, '(3ES14.6)') ins%sk_scale
                    call meta_add(gx, 'sk_scale_sigma_pi_delta', adjustl(buf))
                    if (ins%sk_beta >= 0.0_dp) then
                        write(buf, '(ES14.6)') ins%sk_beta
                        call meta_add(gx, 'sk_beta_override', adjustl(buf))
                    end if
                    write(buf, '(I0,A,I0)') pab(is)%npair, ' / ', pba(is)%npair
                    call meta_add(gx, 'sk_pairs_ab_ba', adjustl(buf))
                    if (pab(is)%npair > 0) then
                        write(buf, '(F10.4)') pab(is)%dmin
                        call meta_add(gx, 'sk_min_distance_ang', adjustl(buf))
                    end if
                    call meta_add(gx, 'contact_overlap', 'S_AB = 0')
                end if
                if (any(inp%layer_potential_ev(1:count(stk%mat == 2)) /= 0.0_dp)) &
                    call meta_add(gx, 'layer_potential', 'H_mn + (V_m+V_n)/2 S_mn')
            end if
            write(buf, '(I0)') stk%nlayer
            call meta_add(gx, 'central_layers', adjustl(buf))
            buf = ''
            do l = 1, min(stk%nlayer, 64)
                buf(l:l) = merge('A', 'B', stk%mat(l) == 1)
            end do
            call meta_add(gx, 'central_stack_bottom_to_top', trim(buf))
            call meta_add(gx, 'leads_bottom_top', merge('A', merge('B', '-', stk%bot_lead == 2), &
                          stk%bot_lead == 1) // ' ' // merge('A', merge('B', '-', stk%top_lead == 2), &
                          stk%top_lead == 1))
            call meta_add(gx, 'central_structure', trim(central_file))
            call meta_add(gx, 'central_structure_fingerprint', structure_fingerprint(stc))
            if (len_trim(inp%regions_file) > 0) then
                call meta_add(gx, 'regions_file', trim(inp%regions_file))
            else
                call meta_add(gx, 'regions_file', '(default: one region per principal layer)')
            end if
            call energy_reference_meta(gx)
            if (inp%nscan > 0) then
                call meta_add(gx, 'scan_parameter', trim(inp%scan_parameter))
                write(buf, '(ES20.12)') inp%scan_values(is)
                call meta_add(gx, 'scan_value', adjustl(buf))
                write(buf, '(I0,A,I0)') is, ' / ', inp%nscan
                call meta_add(gx, 'scan_index', adjustl(buf))
            end if
            ! (o número de threads não entra: o arquivo é idêntico bit a bit para qualquer contagem)
        end subroutine fill_meta

        subroutine energy_reference_meta(gx)
            type(green_data_t), intent(inout) :: gx
            character(len=64) :: buf
            if (mats(1)%hh%has_meta) then
                if (mats(1)%hh%fermi_shift_applied) then
                    call meta_add(gx, 'energy_reference', 'E-E_F')
                else
                    call meta_add(gx, 'energy_reference', 'absolute')
                end if
                write(buf, '(ES20.12)') mats(1)%hh%fermi_ev
                call meta_add(gx, 'substrate_fermi_ev', adjustl(buf))
            else
                call meta_add(gx, 'energy_reference', 'unknown')
            end if
            if (nmat == 2) then
                if (mats(2)%hh%has_meta) then
                    write(buf, '(ES20.12)') mats(2)%hh%fermi_ev
                    call meta_add(gx, 'film_fermi_ev', adjustl(buf))
                    if (mats(1)%hh%has_meta .and. (mats(1)%hh%fermi_shift_applied .or. &
                        mats(2)%hh%fermi_shift_applied)) then
                        call meta_add(gx, 'warning', 'Fermi-shifted .hsio in a heterostructure: ' // &
                                      'align with shift_ev from a junction calculation, not by equal E_F')
                    end if
                end if
            end if
        end subroutine energy_reference_meta

    end subroutine run_green

    !> Entrada com o valor is de &SCAN aplicado (is = 1 sem &SCAN: a própria entrada).
    function scan_input(inp, is, normal) result(o)
        type(green_input_t), intent(in) :: inp
        integer,             intent(in) :: is
        real(dp),            intent(in) :: normal(3)   !< normal unitária do substrato (sentido de a_s)
        type(green_input_t) :: o
        real(dp) :: v
        o = inp
        if (inp%nscan < 1) return
        v = inp%scan_values(is)
        select case (trim(inp%scan_parameter))
        case ('lambda_h');       o%lambda_h = v
        case ('lambda_s');       o%lambda_s = v
        case ('sk_scale_sigma'); o%sk_scale(1) = v
        case ('sk_scale_pi');    o%sk_scale(2) = v
        case ('sk_scale_delta'); o%sk_scale(3) = v
        case ('sk_beta');        o%sk_beta = v
        case ('film_shift_ev');  o%film_shift_ev = v
        case ('offset_x');       o%interface_offset(1) = v
        case ('offset_y');       o%interface_offset(2) = v
        case ('offset_z');       o%interface_offset(3) = v
        case ('offset_normal')
            ! deslocamento extra ao longo da normal do substrato (afasta o filme)
            o%interface_offset = inp%interface_offset + v * normal
        end select
    end function scan_input

end module green_engine
