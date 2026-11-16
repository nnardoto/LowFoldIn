! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Motor spectral (docs/design/contracts.md §3, §3.2; usa §4 e §5): lê
!> spectral.nml, os arquivos HSIO de H/S, o arquivo de estrutura e regions.nml,
!> diagonaliza H(k) ao longo de um caminho k (&KPATH) ou numa grade k 2D (&KGRID,
!> para mapas de energia constante) e calcula os pesos de projeção por
!> região (Mulliken ou Löwdin, escolhido por `method`).
!>
!> Posição no pipeline: backend do executável `spectral`. A saída é um
!> spectral_data_t (módulo spectral_data), escrito em disco pelo driver e
!> depois alargado/plotado pelo specmap.
!>
!> Convenções:
!>   - Energias em eV, pontos k em coordenadas recíprocas fracionárias, todos os
!>     índices base 1; ordem dos átomos = arquivo de estrutura = ordem do HSIO.
!>   - O grupo de região 1 é sempre 'total' (todos os átomos); as regiões do usuário
!>     são os grupos 2..nregion+1 (logo sd%nregion = 1 + número de regiões em regions.nml).
!>   - Layout dos arrays: sd%energy(nband, nk, nspin),
!>     sd%weight(nband, nk, group, nspin). nspin = 2 apenas para colinear.
!>
!> O trabalho por ponto k fica em `process_kpoint`, que só toca os arrays e os
!> espaços de trabalho que recebe. O laço (spin, k) de `run_spectral` é
!> paralelizado com OpenMP (quando compilado com -fopenmp): cada thread possui
!> um conjunto de espaços de trabalho, alocado uma única vez, e cada iteração escreve
!> apenas suas próprias fatias de sd%energy / sd%weight, de modo que o resultado não
!> depende do número de threads. Sem -fopenmp as linhas `!$` somem e o laço roda em série.
!> As rotinas da biblioteca nunca param o programa: os erros são devolvidos por ierr/errmsg.
!> (Exceções: HSDataFormat%ReadFile para o programa diante de um arquivo HSIO corrompido;
!> arquivos ausentes são detectados aqui antes de chamá-lo. process_kpoint executa
!> `error stop` diante de um erro de programação: Löwdin solicitado sem os argumentos
!> opcionais obrigatórios.)
module spectral_engine
    use, intrinsic :: iso_fortran_env, only: int64, iostat_end, output_unit
    use types, only: dp
    use hs_io, only: HSDataFormat
    use core, only: compute_Hk_at_kpoint
    use structure_io, only: structure_t, read_structure, structure_fingerprint, &
                            reciprocal_lattice, species_norb_of
    use atomic_regions, only: region_set_t, read_regions_file, validate_regions, &
                              regions_are_partition
    use spectral_projection, only: orbital_groups_t, build_orbital_groups, &
                                   overlap_times_vectors, mulliken_weights, lowdin_weights, &
                                   generalized_eigensolve
    use kpath, only: kpath_t, build_kpath, kgrid_t, build_kgrid
    !$ use omp_lib, only: omp_get_num_threads
    use spectral_data, only: spectral_data_t, SPIN_NONE, SPIN_COLLINEAR, SPIN_NONCOLLINEAR, &
                             METHOD_MULLIKEN, METHOD_LOWDIN, EREF_UNKNOWN, EREF_FERMI, EREF_ABSOLUTE
    implicit none
    private

    ! ==================================================================
    ! Parâmetros e tipos
    ! ==================================================================

    integer, parameter, public :: MAX_KVERTEX   = 100
    integer, parameter, public :: SPEC_PATH_LEN = 512
    !> Limite de avisos impressos por categoria em check_weights.
    integer, parameter :: MAX_WARN_LINES = 20

    !> Espelho dos namelists &SPECTRAL e &KPATH ou &KGRID (contracts §3). Os valores
    !> padrão aqui são os padrões do arquivo de entrada. Os campos são
    !> strings/números simples; a validação ocorre em read_spectral_input.
    !> vertex(:, i) = i-ésimo vértice do caminho k (fracionário); npoints(i) = pontos
    !> do segmento do vértice i ao i+1 (significado exato: veja o módulo kpath).
    type, public :: spectral_input_t
        character(len=SPEC_PATH_LEN) :: h_file               = ''
        character(len=SPEC_PATH_LEN) :: h_up_file            = ''
        character(len=SPEC_PATH_LEN) :: h_dn_file            = ''
        character(len=SPEC_PATH_LEN) :: s_file               = ''
        character(len=SPEC_PATH_LEN) :: structure_file       = ''
        character(len=SPEC_PATH_LEN) :: regions_file         = ''
        logical                      :: allow_no_fingerprint = .false.  ! aceita regions.nml sem fingerprint
        character(len=16)            :: spin_mode            = 'none'
        character(len=16)            :: method               = 'mulliken'
        character(len=16)            :: energy_reference     = 'unknown'
        real(dp)                     :: fermi_energy_ev      = 0.0_dp
        logical                      :: fermi_known          = .false.
        character(len=SPEC_PATH_LEN) :: output_prefix        = 'slab'
        logical                      :: write_text           = .false.
        real(dp)                     :: tol_norm             = 1.0e-8_dp   ! tolerância de normalização (check_weights)
        real(dp)                     :: tol_negative_warn    = 1.0e-3_dp   ! limiar de aviso de peso negativo (Mulliken)
        real(dp)                     :: tol_s_min            = 1.0e-10_dp  ! Löwdin: min autovalor de S deve exceder este valor
        ! &KPATH
        integer                      :: nvertex              = 0
        real(dp)                     :: vertex(3, MAX_KVERTEX) = 0.0_dp
        character(len=16)            :: label(MAX_KVERTEX)   = ''
        integer                      :: npoints(MAX_KVERTEX - 1) = 0
        ! &KGRID (alternativa a &KPATH; veja kpath::build_kgrid)
        character(len=8)             :: kmode                = 'path'   ! 'path' (&KPATH) ou 'grid' (&KGRID)
        real(dp)                     :: grid_origin(3)       = 0.0_dp
        real(dp)                     :: grid_span1(3)        = 0.0_dp
        real(dp)                     :: grid_span2(3)        = 0.0_dp
        integer                      :: grid_n1              = 0
        integer                      :: grid_n2              = 0
    end type spectral_input_t

    !> Todos os valores são medidos com system_clock (tempo de parede), em segundos.
    !> build_s, diag_s e proj_s são ACUMULADOS sobre todas as threads do laço
    !> em k (thread-segundos: soma do tempo de parede de cada thread na fase);
    !> loop_wall_s é o tempo de parede do laço em k inteiro. Com uma thread,
    !> build_s + diag_s + proj_s <= loop_wall_s.
    type, public :: spectral_timing_t
        real(dp) :: read_s      = 0.0_dp  ! estrutura, regiões, HSIO (parede)
        real(dp) :: build_s     = 0.0_dp  ! construção de H(k), S(k) (thread-s acumulados)
        real(dp) :: diag_s      = 0.0_dp  ! autoresolvedor generalizado (thread-s acumulados)
        real(dp) :: proj_s      = 0.0_dp  ! S*C + pesos (thread-s acumulados)
        real(dp) :: loop_wall_s = 0.0_dp  ! tempo de parede do laço em k
        real(dp) :: write_s     = 0.0_dp  ! preenchido pelo driver (parede)
        real(dp) :: total_s     = 0.0_dp  ! parede
        integer  :: nthreads    = 1       ! threads que executaram o laço em k (1 sem OpenMP)
    end type spectral_timing_t

    public :: read_spectral_input, run_spectral, process_kpoint

contains

    ! ==================================================================
    ! Entrada
    ! ==================================================================

    !> Lê e valida os namelists &SPECTRAL e &KPATH ou &KGRID de `filename`
    !> (os grupos podem aparecer em qualquer ordem no arquivo). Exatamente um
    !> entre &KPATH e &KGRID deve estar presente; inp%kmode = 'path' ou 'grid'.
    !> inp é inicializado com seus padrões primeiro, de modo que as chaves ausentes
    !> do arquivo os mantêm. spin_mode / method / energy_reference são convertidos para minúsculas.
    !> energy_reference='fermi' força fermi_known=.true. e E_F = 0.
    !> Nunca para o programa: em caso de falha ierr = 1 e errmsg descreve o primeiro problema
    !> encontrado (inp fica então apenas parcialmente preenchido).
    subroutine read_spectral_input(filename, inp, ierr, errmsg)
        character(len=*),       intent(in)  :: filename
        type(spectral_input_t), intent(out) :: inp
        integer,                intent(out) :: ierr
        character(len=*),       intent(out) :: errmsg

        ! --- variáveis do namelist &SPECTRAL ---
        ! arquivos de entrada
        character(len=SPEC_PATH_LEN) :: h_file, h_up_file, h_dn_file    ! hamiltonianos: sem spin ou não colinear | up | down
        character(len=SPEC_PATH_LEN) :: s_file                          ! matriz de overlap
        character(len=SPEC_PATH_LEN) :: structure_file, regions_file    ! estrutura e regiões atômicas
        logical                      :: allow_no_fingerprint            ! aceita regiões sem fingerprint
        ! modo / método
        character(len=16)            :: spin_mode, method               ! 'none'|'collinear'|'noncollinear'; 'mulliken'|'lowdin'
        ! nível de Fermi
        character(len=16)            :: energy_reference                ! 'fermi'|'absolute'|'unknown'
        real(dp)                     :: fermi_energy_ev                 ! E_F (eV)
        logical                      :: fermi_known
        ! saída
        character(len=SPEC_PATH_LEN) :: output_prefix
        logical                      :: write_text
        ! tolerâncias
        real(dp)                     :: tol_norm, tol_negative_warn, tol_s_min
        ! --- variáveis do namelist &KPATH ---
        integer                      :: nvertex
        real(dp)                     :: vertex(3, MAX_KVERTEX)
        character(len=64)            :: label(MAX_KVERTEX)   ! mais largo que o armazenado (16) para detectar truncamento
        integer                      :: npoints(MAX_KVERTEX - 1)
        ! --- variáveis do namelist &KGRID ---
        real(dp)                     :: origin(3), span1(3), span2(3)
        integer                      :: n1, n2

        namelist /SPECTRAL/ h_file, h_up_file, h_dn_file, s_file, structure_file, regions_file, &
                            allow_no_fingerprint, &
                            spin_mode, method, &
                            energy_reference, fermi_energy_ev, fermi_known, &
                            output_prefix, write_text, &
                            tol_norm, tol_negative_warn, tol_s_min
        namelist /KPATH/ nvertex, vertex, label, npoints
        namelist /KGRID/ origin, span1, span2, n1, n2

        integer            :: iu, ios, is
        character(len=512) :: iomsg
        logical            :: exists, have_path, have_grid

        ierr   = 0
        errmsg = ''
        iomsg  = ''   ! falhou() monta a mensagem mesmo sem erro: iomsg precisa estar definido

        ! 1. A leitura de namelist só sobrescreve as variáveis presentes no arquivo, então
        ! inicializa as variáveis locais do namelist com os padrões de spectral_input_t
        ! (inp é intent(out): sua inicialização padrão se aplica na entrada).
        h_file               = inp%h_file
        h_up_file            = inp%h_up_file
        h_dn_file            = inp%h_dn_file
        s_file               = inp%s_file
        structure_file       = inp%structure_file
        regions_file         = inp%regions_file
        allow_no_fingerprint = inp%allow_no_fingerprint

        spin_mode            = inp%spin_mode
        method               = inp%method

        energy_reference     = inp%energy_reference
        fermi_energy_ev      = inp%fermi_energy_ev
        fermi_known          = inp%fermi_known

        output_prefix        = inp%output_prefix
        write_text           = inp%write_text

        tol_norm             = inp%tol_norm
        tol_negative_warn    = inp%tol_negative_warn
        tol_s_min            = inp%tol_s_min

        nvertex              = inp%nvertex
        vertex               = inp%vertex
        label                = inp%label
        npoints              = inp%npoints

        origin               = inp%grid_origin
        span1                = inp%grid_span1
        span2                = inp%grid_span2
        n1                   = inp%grid_n1
        n2                   = inp%grid_n2
        have_path            = .false.
        have_grid            = .false.

        ! 2. Abre o arquivo e lê os dois namelists (rewind: independe da ordem).
        inquire(file=filename, exist=exists)
        if (falhou(.not. exists, "input file not found: '" // trim(filename) // "'")) return
        open(newunit=iu, file=filename, status='old', action='read', iostat=ios, iomsg=iomsg)
        if (falhou(ios /= 0, "cannot open '" // trim(filename) // "': " // trim(iomsg))) return

        leitura: block
            read(iu, nml=SPECTRAL, iostat=ios, iomsg=iomsg)
            if (falhou(ios == iostat_end, "namelist &SPECTRAL not found in '" // trim(filename) // "'")) exit leitura
            if (falhou(ios /= 0, "error reading &SPECTRAL: " // trim(iomsg))) exit leitura
            rewind(iu)
            read(iu, nml=KPATH, iostat=ios, iomsg=iomsg)
            have_path = (ios /= iostat_end)
            if (falhou(have_path .and. ios /= 0, "error reading &KPATH: " // trim(iomsg))) exit leitura
            rewind(iu)
            read(iu, nml=KGRID, iostat=ios, iomsg=iomsg)
            have_grid = (ios /= iostat_end)
            if (falhou(have_grid .and. ios /= 0, "error reading &KGRID: " // trim(iomsg))) exit leitura
            if (falhou(.not. have_path .and. .not. have_grid, &
                       "neither &KPATH nor &KGRID found in '" // trim(filename) // "'")) exit leitura
            if (falhou(have_path .and. have_grid, &
                       "&KPATH and &KGRID are mutually exclusive (found both in '" // &
                       trim(filename) // "')")) exit leitura
        end block leitura
        close(iu)
        if (ierr /= 0) return

        ! 3. Copia de volta para inp (normalizando as strings de enumeração) e valida.
        inp%h_file               = h_file
        inp%h_up_file            = h_up_file
        inp%h_dn_file            = h_dn_file
        inp%s_file               = s_file
        inp%structure_file       = structure_file
        inp%regions_file         = regions_file
        inp%allow_no_fingerprint = allow_no_fingerprint

        ! normalização das strings de enumeração
        inp%spin_mode            = lower(adjustl(spin_mode))
        inp%method               = lower(adjustl(method))
        inp%energy_reference     = lower(adjustl(energy_reference))

        inp%fermi_energy_ev      = fermi_energy_ev
        inp%fermi_known          = fermi_known

        inp%output_prefix        = output_prefix
        inp%write_text           = write_text

        inp%tol_norm             = tol_norm
        inp%tol_negative_warn    = tol_negative_warn
        inp%tol_s_min            = tol_s_min

        inp%nvertex              = nvertex
        inp%vertex               = vertex
        ! label é lido mais largo (64) que o armazenado (16): detecta truncamento
        if (falhou(any(len_trim(label) > len(inp%label)), '&KPATH label longer than 16 characters')) return
        inp%label                = label
        inp%npoints              = npoints

        inp%kmode = 'path'
        if (have_grid) inp%kmode = 'grid'
        inp%grid_origin          = origin
        inp%grid_span1           = span1
        inp%grid_span2           = span2
        inp%grid_n1              = n1
        inp%grid_n2              = n2

        ! --- enumerações -------------------------------------------------
        select case (trim(inp%spin_mode))
        case ('none', 'noncollinear')
            if (falhou(len_trim(inp%h_file) == 0, &
                       "h_file is required for spin_mode='" // trim(inp%spin_mode) // "'")) return
        case ('collinear')
            if (falhou(len_trim(inp%h_up_file) == 0 .or. len_trim(inp%h_dn_file) == 0, &
                       "h_up_file and h_dn_file are required for spin_mode='collinear'")) return
        case default
            call fail("invalid spin_mode '" // trim(inp%spin_mode) // &
                      "' (expected 'none', 'collinear' or 'noncollinear')")
            return
        end select

        select case (trim(inp%method))
        case ('mulliken', 'lowdin')
        case default
            call fail("invalid method '" // trim(inp%method) // "' (expected 'mulliken' or 'lowdin')")
            return
        end select

        select case (trim(inp%energy_reference))
        case ('fermi')
            ! As energias já são relativas a E_F: E_F = 0 por definição.
            inp%fermi_known = .true.
            inp%fermi_energy_ev = 0.0_dp
        case ('absolute', 'unknown')
        case default
            call fail("invalid energy_reference '" // trim(inp%energy_reference) // &
                      "' (expected 'fermi', 'absolute' or 'unknown')")
            return
        end select

        ! --- arquivos obrigatórios ---------------------------------------
        if (falhou(len_trim(inp%s_file) == 0, 's_file is required')) return
        if (falhou(len_trim(inp%structure_file) == 0, 'structure_file is required')) return
        if (falhou(len_trim(inp%output_prefix) == 0, 'output_prefix must not be empty')) return

        ! --- tolerâncias -------------------------------------------------
        if (falhou(.not. (inp%tol_norm > 0.0_dp), 'tol_norm must be > 0')) return
        if (falhou(.not. (inp%tol_negative_warn > 0.0_dp), 'tol_negative_warn must be > 0')) return
        if (falhou(.not. (inp%tol_s_min > 0.0_dp), 'tol_s_min must be > 0')) return

        ! --- grade k (a geometria completa é verificada em build_kgrid) ----
        if (have_grid) then
            if (falhou(inp%grid_n1 < 2 .or. inp%grid_n2 < 2, &
                       '&KGRID n1 and n2 must be >= 2 (got n1=' // itoa(inp%grid_n1) // &
                       ', n2=' // itoa(inp%grid_n2) // ')')) return
            if (falhou(all(inp%grid_span1 == 0.0_dp) .or. all(inp%grid_span2 == 0.0_dp), &
                       '&KGRID span1 and span2 must be non-zero')) return
            return
        end if

        ! --- caminho k ---------------------------------------------------
        if (falhou(inp%nvertex < 2 .or. inp%nvertex > MAX_KVERTEX, &
                   '&KPATH nvertex must be in 2..' // itoa(MAX_KVERTEX) // &
                   ' (got ' // itoa(inp%nvertex) // ')')) return
        do is = 1, inp%nvertex - 1
            if (falhou(inp%npoints(is) < 1, &
                       '&KPATH npoints(' // itoa(is) // ') must be >= 1 (got ' // &
                       itoa(inp%npoints(is)) // ')')) return
        end do

    contains

        !> Define o estado de erro.
        subroutine fail(m)
            character(len=*), intent(in) :: m
            ierr = 1
            errmsg = m
        end subroutine fail

        !> Se `cond` é verdadeira, registra a falha `m`; devolve `cond`.
        logical function falhou(cond, m)
            logical,          intent(in) :: cond
            character(len=*), intent(in) :: m
            falhou = cond
            if (cond) call fail(m)
        end function falhou

    end subroutine read_spectral_input

    ! ==================================================================
    ! Driver do cálculo
    ! ==================================================================

    !> Cálculo completo: lê as entradas, verifica a consistência, constrói os grupos de orbitais
    !> e o caminho k, diagonaliza H(k) em cada (spin, k) e preenche `sd`.
    !>   inp           entrada validada (veja read_spectral_input)
    !>   sd            dados espectrais de saída; sd%energy(nband,nk,nspin) em eV,
    !>                 sd%weight(nband,nk,group,nspin); grupo 1 = 'total'
    !>   timing        tempos de parede/thread; write_s fica para o driver
    !>   show_progress opcional (padrão .false.): imprime as linhas de etapa [1/3]..[3/3]
    !>                 e o progresso/ETA do laço em k no stdout
    !> Nunca para o programa: em caso de falha ierr = 1 e errmsg contém a mensagem (para falhas
    !> no laço em k ela nomeia ispin, ik e kfrac da menor iteração que falhou).
    !> Os avisos de check_weights são impressos, não devolvidos. Exceção:
    !> HSDataFormat%ReadFile pode parar o programa diante de um arquivo HSIO corrompido.
    subroutine run_spectral(inp, sd, timing, ierr, errmsg, show_progress)
        type(spectral_input_t),  intent(in)  :: inp
        type(spectral_data_t),   intent(out) :: sd
        type(spectral_timing_t), intent(out) :: timing
        integer,                 intent(out) :: ierr
        character(len=*),        intent(out) :: errmsg
        logical, optional,       intent(in)  :: show_progress

        type(structure_t)               :: st
        type(region_set_t)              :: rs
        type(HSDataFormat), allocatable :: hdat(:)
        type(HSDataFormat)              :: sdat
        type(orbital_groups_t)          :: groups
        type(kpath_t)                   :: path
        type(kgrid_t)                   :: kg
        real(dp), allocatable           :: kfrac(:,:)          ! (3,nk) pontos do caminho ou da grade
        integer, allocatable            :: list_start(:), list_atoms(:)
        complex(dp), allocatable        :: Hk(:,:), Sk(:,:), C(:,:), SC(:,:)
        real(dp), allocatable           :: E(:), W(:,:)
        character(len=SPEC_PATH_LEN)    :: hname(2)
        integer(int64)                  :: c_start, c0, c1, crate
        integer                         :: nspin, nreg, ngroup, nat, n, ispin, ik, ia, ir, g, pos, smode, kerr, perr
        logical                         :: have_regions, partition, use_lowdin
        real(dp)                        :: smin
        character(len=512)              :: msg
        ! Estado do laço em k (OpenMP)
        integer                         :: nk, nit, it, first_fail, ff, nthr
        integer, allocatable            :: it_kerr(:), it_perr(:)
        real(dp), allocatable           :: it_smin(:)
        real(dp)                        :: t_b, t_d, t_p
        integer(int64)                  :: cl0, cl1, cp
        integer                         :: completed, next_report, report_step
        logical                         :: report_progress
        real(dp)                        :: progress_elapsed, progress_eta

        ierr = 0
        errmsg = ''
        call system_clock(c_start, crate)
        use_lowdin = trim(inp%method) == 'lowdin'
        report_progress = .false.
        if (present(show_progress)) report_progress = show_progress

        if (report_progress) then
            print '(A)', '[1/3] Loading and validating input data'
            flush(output_unit)
        end if

        if (trim(inp%method) /= 'mulliken' .and. trim(inp%method) /= 'lowdin') then
            call fail("invalid method '" // trim(inp%method) // "'")
            return
        end if

        select case (trim(inp%spin_mode))
        case ('none')
            smode = SPIN_NONE;          nspin = 1; hname(1) = inp%h_file
        case ('noncollinear')
            smode = SPIN_NONCOLLINEAR;  nspin = 1; hname(1) = inp%h_file
        case ('collinear')
            smode = SPIN_COLLINEAR;     nspin = 2
            hname(1) = inp%h_up_file;   hname(2) = inp%h_dn_file
        case default
            call fail("invalid spin_mode '" // trim(inp%spin_mode) // "'")
            return
        end select

        ! ---------------- 1. leitura das entradas -----------------------
        call system_clock(c0)
        call read_structure(inp%structure_file, st, ierr, msg)
        if (ierr /= 0) then
            call fail('structure: ' // trim(msg))
            return
        end if

        have_regions = len_trim(inp%regions_file) > 0
        if (have_regions) then
            call read_regions_file(inp%regions_file, rs, ierr, msg)
            if (ierr /= 0) then
                call fail('regions: ' // trim(msg))
                return
            end if
            call validate_regions(rs, st, inp%allow_no_fingerprint, ierr, msg)
            if (ierr /= 0) then
                call fail('regions: ' // trim(msg))
                return
            end if
        end if

        if (.not. file_exists(inp%s_file)) then
            call fail("S file not found: '" // trim(inp%s_file) // "'")
            return
        end if
        do ispin = 1, nspin
            if (.not. file_exists(hname(ispin))) then
                call fail("H file not found: '" // trim(hname(ispin)) // "'")
                return
            end if
        end do
        allocate(hdat(nspin))
        call sdat%ReadFile(trim(inp%s_file))
        do ispin = 1, nspin
            call hdat(ispin)%ReadFile(trim(hname(ispin)))
        end do
        call system_clock(c1)
        timing%read_s = elapsed(c0, c1, crate)

        ! ---------------- 2. consistência -------------------------------
        do ispin = 1, nspin
            call compare_hs(hdat(ispin), sdat, trim(hname(ispin)), trim(inp%s_file), ierr, msg)
            if (ierr /= 0) then
                call fail(trim(msg))
                return
            end if
        end do
        if (.not. allocated(sdat%Orbs) .or. .not. allocated(sdat%Data) .or. &
            .not. allocated(sdat%IndexList)) then
            call fail("S file '" // trim(inp%s_file) // "' has no data")
            return
        end if
        nat = sdat%nAtoms
        n = sdat%NOrbs
        if (size(sdat%Orbs) /= nat .or. sum(sdat%Orbs) /= n) then
            call fail("HSIO Orbs inconsistent with nAtoms/NOrbs in '" // trim(inp%s_file) // "'")
            return
        end if
        if (nat /= st%nat) then
            call fail('HSIO nAtoms (' // itoa(nat) // ') differs from structure nat (' // &
                      itoa(st%nat) // ") in '" // trim(inp%structure_file) // "'")
            return
        end if
        if (st%has_species) then
            do ia = 1, nat
                if (species_norb_of(st, ia) /= sdat%Orbs(ia)) then
                    call fail('atom ' // itoa(ia) // " (species '" // trim(st%label(ia)) // &
                              "'): structure nOrb = " // itoa(species_norb_of(st, ia)) // &
                              ' but HSIO Orbs = ' // itoa(sdat%Orbs(ia)))
                    return
                end if
            end do
        end if

        ! ---------------- 3. grupos de orbitais -------------------------
        ! Grupo 1 = todos os átomos ('total'); grupos 2.. = regiões do usuário, na ordem
        ! do arquivo. As listas de átomos são concatenadas em formato CSR (list_start/list_atoms).
        nreg = 0
        if (have_regions) nreg = rs%nregion
        ngroup = 1 + nreg
        allocate(list_start(ngroup + 1))
        list_start(1) = 1
        list_start(2) = 1 + nat
        do ir = 1, nreg
            list_start(ir + 2) = list_start(ir + 1) + size(rs%regions(ir)%atoms)
        end do
        allocate(list_atoms(list_start(ngroup + 1) - 1))
        list_atoms(1:nat) = [(ia, ia = 1, nat)]
        pos = nat
        do ir = 1, nreg
            list_atoms(pos + 1:pos + size(rs%regions(ir)%atoms)) = rs%regions(ir)%atoms
            pos = pos + size(rs%regions(ir)%atoms)
        end do
        call build_orbital_groups(sdat%Orbs, ngroup, list_start, list_atoms, groups)
        ! partition: as regiões do usuário cobrem cada átomo exatamente uma vez, então seus
        ! pesos devem somar 1 (verificado em check_weights).
        partition = .false.
        if (have_regions) partition = regions_are_partition(rs)

        ! ---------------- 4. caminho k ou grade k -----------------------
        if (trim(inp%kmode) == 'grid') then
            call build_kgrid(inp%grid_origin, inp%grid_span1, inp%grid_span2, inp%grid_n1, &
                             inp%grid_n2, reciprocal_lattice(st), kg, ierr, msg)
            if (ierr /= 0) then
                call fail(trim(msg))
                return
            end if
            kfrac = kg%kfrac
        else
            call build_kpath(inp%vertex(:, 1:inp%nvertex), inp%npoints(1:inp%nvertex - 1), &
                             inp%label(1:inp%nvertex), reciprocal_lattice(st), path, ierr, msg)
            if (ierr /= 0) then
                call fail('k-path: ' // trim(msg))
                return
            end if
            kfrac = path%kfrac
        end if
        nk = size(kfrac, 2)

        ! ---------------- 5. metadados dos dados espectrais -------------
        ! Preenche tudo em sd exceto os resultados do laço em k (energy/weight são
        ! alocados zerados aqui e escritos na etapa 6).
        sd%nk        = nk
        sd%nband     = n
        sd%nspin     = nspin
        sd%nregion   = ngroup
        sd%nvertex   = path%nvertex     ! 0 para grade (path vazio)
        sd%spin_mode = smode
        sd%method    = METHOD_MULLIKEN
        if (use_lowdin) sd%method = METHOD_LOWDIN
        select case (trim(inp%energy_reference))
        case ('fermi');    sd%energy_ref = EREF_FERMI
        case ('absolute'); sd%energy_ref = EREF_ABSOLUTE
        case default;      sd%energy_ref = EREF_UNKNOWN
        end select
        sd%fermi_known = inp%fermi_known
        sd%fermi_energy_ev = 0.0_dp
        if (inp%fermi_known) sd%fermi_energy_ev = inp%fermi_energy_ev
        if (sd%energy_ref == EREF_FERMI) then
            sd%fermi_known = .true.
            sd%fermi_energy_ev = 0.0_dp
        end if
        sd%fingerprint = structure_fingerprint(st)
        if (smode == SPIN_COLLINEAR) then
            sd%source = 'H_up=' // trim(inp%h_up_file) // ' H_dn=' // trim(inp%h_dn_file) // &
                        ' S=' // trim(inp%s_file)
        else
            sd%source = 'H=' // trim(inp%h_file) // ' S=' // trim(inp%s_file)
        end if
        sd%kfrac = kfrac
        if (trim(inp%kmode) == 'grid') then
            ! Grade: kdist não tem significado (zeros); a geometria vai no bloco KGRID.
            allocate(sd%kdist(nk), source=0.0_dp)
            allocate(sd%vertex_k(0), sd%vertex_label(0))
            sd%is_grid     = .true.
            sd%grid_n1     = kg%n1
            sd%grid_n2     = kg%n2
            sd%grid_origin = kg%origin
            sd%grid_span1  = kg%span1
            sd%grid_span2  = kg%span2
            sd%recip       = reciprocal_lattice(st)
        else
            sd%kdist = path%kdist
            sd%vertex_k = path%vertex_k
            allocate(sd%vertex_label(path%nvertex))
            sd%vertex_label = path%label
        end if
        allocate(sd%region_name(ngroup), sd%region_natoms(ngroup))
        sd%region_name(1) = 'total'
        sd%region_natoms(1) = nat
        do ir = 1, nreg
            sd%region_name(ir + 1) = rs%regions(ir)%name
            sd%region_natoms(ir + 1) = size(rs%regions(ir)%atoms)
        end do
        allocate(sd%energy(n, nk, nspin), source=0.0_dp)
        allocate(sd%weight(n, nk, ngroup, nspin), source=0.0_dp)

        ! ---------------- 6. laço em k (OpenMP sobre o índice plano (ispin, ik)) --
        ! it = (ispin-1)*nk + ik (invertido abaixo como ispin = (it-1)/nk + 1).
        ! t_b/t_d/t_p são reduções OpenMP (thread-segundos). Cada iteração escreve apenas
        ! sd%energy(:,ik,ispin), sd%weight(:,ik,:,ispin) e sua própria entrada dos arrays
        ! de status, de modo que o resultado independe do número de threads e do
        ! escalonamento. Os espaços de trabalho são privados e alocados uma vez por thread.
        ! A E/S de progresso, quando habilitada, é serializada; não há stop dentro da
        ! região. A primeira iteração que falha (menor it) é reportada após a região paralela.
        ! Iterações acima da menor falha vista até então são puladas; as abaixo dela
        ! sempre executam, então a falha reportada é determinística.
        nit = nspin * nk
        allocate(it_kerr(nit), it_perr(nit), it_smin(nit))
        it_kerr = 0
        it_perr = 0
        it_smin = 0.0_dp
        first_fail = nit + 1
        nthr = 1
        completed = 0
        report_step = max(1, (nit + 9) / 10)
        next_report = report_step
        t_b = 0.0_dp; t_d = 0.0_dp; t_p = 0.0_dp
        if (report_progress) then
            print '(A,I0,A)', '[2/3] Solving ', nit, ' (spin, k) points'
        end if
        call system_clock(cl0)

        !$omp parallel default(none) &
        !$omp shared(n, ngroup, nk, nit, hdat, sdat, kfrac, groups, use_lowdin, inp, sd, &
        !$omp        it_kerr, it_perr, it_smin, first_fail, nthr, report_progress, completed, &
        !$omp        next_report, report_step, cl0, crate, cp, progress_elapsed, progress_eta) &
        !$omp private(it, ispin, ik, g, ff, kerr, perr, smin, Hk, Sk, C, SC, E, W) &
        !$omp reduction(+: t_b, t_d, t_p)
        allocate(Hk(n, n), Sk(n, n), C(n, n), SC(n, n), E(n), W(n, ngroup))
        !$omp single
        !$ nthr = omp_get_num_threads()
        !$omp end single
        !$omp do schedule(dynamic, 1)
        do it = 1, nit
            !$omp atomic read
            ff = first_fail
            if (it > ff) cycle
            ispin = (it - 1) / nk + 1
            ik = it - (ispin - 1) * nk
            call process_kpoint(n, ngroup, hdat(ispin)%Data, sdat%Data, sdat%IndexList, &
                                kfrac(:, ik), groups, Hk, Sk, C, SC, E, W, &
                                t_b, t_d, t_p, kerr, &
                                lowdin=use_lowdin, tol_s_min=inp%tol_s_min, smin=smin, &
                                proj_ierr=perr)
            if (kerr /= 0 .or. perr /= 0) then
                it_kerr(it) = kerr
                it_perr(it) = perr
                it_smin(it) = smin
                !$omp atomic update
                first_fail = min(first_fail, it)
                cycle
            end if
            sd%energy(:, ik, ispin) = E
            do g = 1, ngroup
                sd%weight(:, ik, g, ispin) = W(:, g)
            end do
            if (report_progress) then
                !$omp critical (spectral_progress)
                completed = completed + 1
                if (completed >= next_report .or. completed == nit) then
                    call system_clock(cp)
                    progress_elapsed = elapsed(cl0, cp, crate)
                    progress_eta = progress_elapsed * real(nit - completed, dp) / real(completed, dp)
                    write(*, '(A,I0,A,I0,A,F6.1,A,F9.2,A,F9.2,A)') '      ', completed, '/', nit, &
                        '  ', 100.0_dp * real(completed, dp) / real(nit, dp), &
                        '%  elapsed ', progress_elapsed, ' s  ETA ', progress_eta, ' s'
                    flush(output_unit)
                    next_report = min(nit, next_report + report_step)
                end if
                !$omp end critical (spectral_progress)
            end if
        end do
        !$omp end do
        deallocate(Hk, Sk, C, SC, E, W)
        !$omp end parallel

        call system_clock(cl1)
        timing%loop_wall_s = elapsed(cl0, cl1, crate)
        timing%build_s  = t_b
        timing%diag_s   = t_d
        timing%proj_s   = t_p
        timing%nthreads = nthr

        ! Reporta a menor iteração que falhou (determinístico para qualquer número de threads).
        ! perr == 1 é o código de S não definida positiva do Löwdin em lowdin_weights.
        if (first_fail <= nit) then
            it = first_fail
            ispin = (it - 1) / nk + 1
            ik = it - (ispin - 1) * nk
            perr = it_perr(it)
            kerr = it_kerr(it)
            smin = it_smin(it)
            if (perr == 1) then
                write(msg, '(A,I0,A,I0,A,3F12.7,A,ES12.4,A,ES10.3,A)') &
                    'overlap matrix not positive definite above tol_s_min at ispin=', ispin, &
                    ' ik=', ik, ' kfrac=(', kfrac(:, ik), '): min eigenvalue smin=', smin, &
                    ' (tol_s_min=', inp%tol_s_min, ')'
            else if (perr /= 0) then
                write(msg, '(A,I0,A,I0,A,I0,A,3F12.7,A)') 'Lowdin projection failed (ierr=', &
                    perr, ') at ispin=', ispin, ' ik=', ik, ' kfrac=(', kfrac(:, ik), ')'
            else if (kerr > n) then
                ! Convenção do zhegv (LAPACK): info > n significa que o menor principal líder
                ! de S(k) de ordem info-n não é definido positivo.
                write(msg, '(A,I0,A,I0,A,I0,A,3F12.7,A)') &
                    'overlap S(k) not positive definite (zhegv info=', &
                    kerr, ') at ispin=', ispin, ' ik=', ik, ' kfrac=(', kfrac(:, ik), ')'
            else
                write(msg, '(A,I0,A,I0,A,I0,A,3F12.7,A)') 'generalized eigensolver failed (info=', &
                    kerr, ') at ispin=', ispin, ' ik=', ik, ' kfrac=(', kfrac(:, ik), ')'
            end if
            call fail(trim(msg))
            return
        end if

        ! ---------------- 7. verificações -------------------------------
        if (report_progress) print '(A)', '[3/3] Checking projection weights'
        call check_weights(sd, partition, inp%tol_norm, inp%tol_negative_warn, use_lowdin)

        call system_clock(c1)
        timing%total_s = elapsed(c_start, c1, crate)

    contains

        subroutine fail(m)
            character(len=*), intent(in) :: m
            ierr = 1
            errmsg = m
        end subroutine fail

    end subroutine run_spectral

    ! ==================================================================
    ! Corpo por ponto k: toca apenas seus argumentos (thread-safe por construção).
    ! ==================================================================

    !> Um ponto (spin, k): constrói H(k), S(k), resolve H C = S C E e calcula os
    !> pesos de projeção W(banda, grupo).
    !>   Hr, Sr, Rpt  H(R), S(R) no espaço real e vetores de rede n_R (triplas
    !>                inteiras), como armazenados em HSDataFormat (Data / IndexList);
    !>                soma de Fourier H(k) = sum_R H(R) exp(+i 2 pi k.n_R)
    !>                (core::compute_Hk_at_kpoint, gauge de rede)
    !>   kfrac        ponto k em coordenadas recíprocas fracionárias
    !>   Hk,Sk,C,SC   espaços de trabalho (n,n) do chamador (um conjunto por thread);
    !>                no retorno Sk = S(k) original, C = autovetores (colunas)
    !>   E(n), W(n,ngroup)  autovalores em eV (crescentes) e pesos
    !>   t_build/t_diag/t_proj  segundos de parede SOMADOS aos valores passados
    !>   ierr         info do autoresolvedor generalizado (0 = ok); com ierr /= 0
    !>                E, C, W ficam indefinidos e t_proj não é atualizado
    !> Modo Löwdin opcional (lowdin=.true.): pesos de lowdin_weights sobre o Sk
    !> ORIGINAL (SC não é usado); tol_s_min, smin e proj_ierr passam a ser
    !> obrigatórios. proj_ierr = 1 significa min autovalor de Sk = smin <= tol_s_min.
    !> Sem os argumentos opcionais usa-se o caminho de Mulliken (inalterado).
    !> Violar a exigência de argumentos do Löwdin é um erro de programação e
    !> dispara `error stop` (o único stop na API de biblioteca deste módulo).
    subroutine process_kpoint(n, ngroup, Hr, Sr, Rpt, kfrac, groups, Hk, Sk, C, SC, E, W, &
                              t_build, t_diag, t_proj, ierr, lowdin, tol_s_min, smin, proj_ierr)
        integer,                 intent(in)    :: n, ngroup
        complex(dp), contiguous, intent(in)    :: Hr(:,:,:), Sr(:,:,:)
        integer,     contiguous, intent(in)    :: Rpt(:,:)
        real(dp),                intent(in)    :: kfrac(3)
        type(orbital_groups_t),  intent(in)    :: groups
        complex(dp), contiguous, intent(inout) :: Hk(:,:), Sk(:,:), C(:,:), SC(:,:)
        real(dp),                intent(inout) :: E(n), W(n, ngroup)
        real(dp),                intent(inout) :: t_build, t_diag, t_proj
        integer,                 intent(out)   :: ierr
        logical,     optional,   intent(in)    :: lowdin
        real(dp),    optional,   intent(in)    :: tol_s_min
        real(dp),    optional,   intent(out)   :: smin
        integer,     optional,   intent(out)   :: proj_ierr
        integer(int64)                         :: c0, c1, c2, c3, crate
        logical                                :: do_lowdin
        real(dp)                               :: smin_loc
        integer                                :: perr

        do_lowdin = .false.
        if (present(lowdin)) do_lowdin = lowdin
        if (do_lowdin .and. .not. (present(tol_s_min) .and. present(proj_ierr))) &
            error stop 'process_kpoint: lowdin requires tol_s_min and proj_ierr'
        if (present(proj_ierr)) proj_ierr = 0
        if (present(smin)) smin = 0.0_dp

        ! 1. H(k), S(k) por soma de Fourier sobre R.
        call system_clock(c0, crate)
        call compute_Hk_at_kpoint(Hr, Rpt, Hk, kfrac)
        call compute_Hk_at_kpoint(Sr, Rpt, Sk, kfrac)
        call system_clock(c1)
        ! 2. Autoproblema generalizado. generalized_eigensolve trabalha em cópias
        ! internas: Sk permanece o S(k) original, como as projeções exigem.
        call generalized_eigensolve(n, Hk, Sk, E, C, ierr)
        call system_clock(c2)
        if (ierr /= 0) then
            t_build = t_build + elapsed(c0, c1, crate)
            t_diag  = t_diag  + elapsed(c1, c2, crate)
            return
        end if
        ! 3. Pesos de projeção (Löwdin a partir de S^{1/2}C, Mulliken a partir de C e S·C).
        if (do_lowdin) then
            call lowdin_weights(n, n, Sk, C, groups, W, tol_s_min, smin_loc, perr)
            proj_ierr = perr
            if (present(smin)) smin = smin_loc
        else
            call overlap_times_vectors(n, n, Sk, C, SC)
            call mulliken_weights(n, n, C, SC, groups, W)
        end if
        call system_clock(c3)
        t_build = t_build + elapsed(c0, c1, crate)
        t_diag  = t_diag  + elapsed(c1, c2, crate)
        t_proj  = t_proj  + elapsed(c2, c3, crate)
    end subroutine process_kpoint

    ! ==================================================================
    ! Auxiliares
    ! ==================================================================

    !> Diagnósticos sobre os dados espectrais prontos; imprime avisos (limitados a
    !> MAX_WARN_LINES por categoria) e um resumo, nunca falha. Verifica:
    !>   - peso do grupo 1 ('total') = 1 por banda dentro de tol_norm
    !>   - se `partition`, soma sobre as regiões do usuário (grupos 2..) = 1 dentro de tol_norm
    !>   - nenhum peso abaixo de -tol_neg (populações de Mulliken podem ser levemente negativas)
    !> Para Löwdin (lowdin=.true.) os pesos são não negativos por construção; a
    !> checagem de negativos usa então -tol_norm em vez de -tol_negative_warn.
    subroutine check_weights(sd, partition, tol_norm, tol_neg_in, lowdin)
        type(spectral_data_t), intent(in) :: sd
        logical,               intent(in) :: partition
        real(dp),              intent(in) :: tol_norm, tol_neg_in
        logical,               intent(in) :: lowdin
        real(dp)               :: tol_neg
        character(len=8)       :: mname
        integer                :: ispin, ik, ib, g, n_tot, n_part, n_neg
        real(dp)               :: dev, max_tot, max_part, wmin, w

        if (lowdin) then
            tol_neg = tol_norm;   mname = 'Lowdin'
        else
            tol_neg = tol_neg_in; mname = 'Mulliken'
        end if
        n_tot = 0; n_part = 0; n_neg = 0
        max_tot = 0.0_dp; max_part = 0.0_dp; wmin = huge(1.0_dp)
        do ispin = 1, sd%nspin
            do ik = 1, sd%nk
                do ib = 1, sd%nband
                    dev = abs(sd%weight(ib, ik, 1, ispin) - 1.0_dp)
                    max_tot = max(max_tot, dev)
                    if (dev > tol_norm) then
                        n_tot = n_tot + 1
                        if (n_tot <= MAX_WARN_LINES) print '(A,I0,A,I0,A,I0,A,ES10.3)', &
                            '[spectral warning] total weight /= 1: spin=', ispin, ' k=', ik, &
                            ' band=', ib, ' |w-1|=', dev
                    end if
                    if (partition .and. sd%nregion > 1) then
                        dev = abs(sum(sd%weight(ib, ik, 2:sd%nregion, ispin)) - 1.0_dp)
                        max_part = max(max_part, dev)
                        if (dev > tol_norm) then
                            n_part = n_part + 1
                            if (n_part <= MAX_WARN_LINES) print '(A,I0,A,I0,A,I0,A,ES10.3)', &
                                '[spectral warning] sum over regions /= 1: spin=', ispin, ' k=', ik, &
                                ' band=', ib, ' |sum-1|=', dev
                        end if
                    end if
                    do g = 1, sd%nregion
                        w = sd%weight(ib, ik, g, ispin)
                        wmin = min(wmin, w)
                        if (w < -tol_neg) then
                            n_neg = n_neg + 1
                            if (n_neg <= MAX_WARN_LINES) print '(A,I0,A,I0,A,I0,A,A,A,ES10.3)', &
                                '[spectral warning] negative ' // trim(mname) // ' weight: spin=', ispin, ' k=', ik, &
                                ' band=', ib, ' region=', trim(sd%region_name(g)), ' w=', w
                        end if
                    end do
                end do
            end do
        end do
        if (n_tot > MAX_WARN_LINES) print '(A,I0,A)', '[spectral warning] ... ', &
            n_tot - MAX_WARN_LINES, ' more total-normalisation warnings'
        if (n_part > MAX_WARN_LINES) print '(A,I0,A)', '[spectral warning] ... ', &
            n_part - MAX_WARN_LINES, ' more region-sum warnings'
        if (n_neg > MAX_WARN_LINES) print '(A,I0,A)', '[spectral warning] ... ', &
            n_neg - MAX_WARN_LINES, ' more negative-weight warnings'

        print '(A)', 'Normalisation checks:'
        print '(A,ES10.3,A,I0,A,ES9.2,A)', '  max |total - 1|        = ', max_tot, &
            '  (', n_tot, ' above tol_norm=', tol_norm, ')'
        if (partition .and. sd%nregion > 1) then
            print '(A,ES10.3,A,I0,A)', '  max |sum regions - 1|  = ', max_part, &
                '  (', n_part, ' above tol_norm; regions form a partition)'
        else if (sd%nregion > 1) then
            print '(A)', '  sum over regions       : not checked (regions are not a partition)'
        end if
        if (lowdin) then
            print '(A,ES10.3,A,I0,A,ES9.2,A)', '  min Lowdin weight      = ', wmin, &
                '  (', n_neg, ' below -tol_norm=', -tol_neg, ')'
        else
            print '(A,ES10.3,A,I0,A,ES9.2,A)', '  min Mulliken weight    = ', wmin, &
                '  (', n_neg, ' below -tol_negative_warn=', -tol_neg, ')'
        end if
    end subroutine check_weights

    !> Verifica se dois arquivos HSIO (p. ex. H e S) descrevem o mesmo sistema:
    !> nAtoms, NOrbs, Neighbors, Orbs e IndexList (lista de vetores R) devem coincidir.
    !> ierr = 1 e errmsg nomeia o primeiro item divergente.
    subroutine compare_hs(A, B, nameA, nameB, ierr, errmsg)
        type(HSDataFormat), intent(in)  :: A, B
        character(len=*),   intent(in)  :: nameA, nameB
        integer,            intent(out) :: ierr
        character(len=*),   intent(out) :: errmsg
        character(len=:), allocatable   :: what

        ierr = 0
        errmsg = ''
        what = ''
        if (A%nAtoms /= B%nAtoms) then
            what = 'nAtoms'
        else if (A%NOrbs /= B%NOrbs) then
            what = 'NOrbs'
        else if (A%Neighbors /= B%Neighbors) then
            what = 'Neighbors'
        else if (allocated(A%Orbs) .neqv. allocated(B%Orbs)) then
            what = 'Orbs'
        else if (allocated(A%IndexList) .neqv. allocated(B%IndexList)) then
            what = 'IndexList'
        end if
        if (len(what) == 0 .and. allocated(A%Orbs)) then
            if (size(A%Orbs) /= size(B%Orbs)) then
                what = 'Orbs'
            else if (any(A%Orbs /= B%Orbs)) then
                what = 'Orbs'
            end if
        end if
        if (len(what) == 0 .and. allocated(A%IndexList)) then
            if (any(shape(A%IndexList) /= shape(B%IndexList))) then
                what = 'IndexList'
            else if (any(A%IndexList /= B%IndexList)) then
                what = 'IndexList'
            end if
        end if
        if (len(what) > 0) then
            ierr = 1
            errmsg = "HSIO files disagree in " // what // ": '" // nameA // "' vs '" // nameB // "'"
        end if
    end subroutine compare_hs

    logical function file_exists(fname)
        character(len=*), intent(in) :: fname
        inquire(file=trim(fname), exist=file_exists)
    end function file_exists

    !> Segundos entre duas contagens de system_clock (taxa protegida contra 0).
    pure real(dp) function elapsed(c0, c1, rate)
        integer(int64), intent(in) :: c0, c1, rate
        elapsed = real(c1 - c0, dp) / real(max(rate, 1_int64), dp)
    end function elapsed

    pure function lower(s) result(r)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: r
        integer               :: i, ic
        r = s
        do i = 1, len(s)
            ic = iachar(s(i:i))
            if (ic >= iachar('A') .and. ic <= iachar('Z')) r(i:i) = achar(ic + 32)
        end do
    end function lower

    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=24)             :: buf
        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

end module spectral_engine
