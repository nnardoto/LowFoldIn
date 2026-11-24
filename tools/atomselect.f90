! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! tools/atomselect.f90
!
! Ferramenta de linha de comando para regiões atômicas nomeadas (docs/design/contracts.md §6):
! inspeciona um arquivo de estrutura, imprime seu fingerprint, cria/estende um
! arquivo regions.nml com seleções de átomos numéricas ou geométricas, detecta
! camadas ao longo de uma normal (cartesiana) arbitrária, constrói regiões de
! superfície (topo/base) e valida um arquivo de regiões contra uma estrutura.
!
! Uso:
!   atomselect inspect     <structure.nml> [--xyz out.xyz] [--table out.txt]
!   atomselect fingerprint <structure.nml>
!   atomselect select      <structure.nml> --regions <regions.nml> --name NAME --atoms "LIST"
!   atomselect select      <structure.nml> --regions <regions.nml> --name NAME --normal nx ny nz [--smin A] [--smax B]
!   atomselect select      <structure.nml> --regions <regions.nml> --name NAME --box "xlo xhi ylo yhi zlo zhi"
!   atomselect layers      <structure.nml> --normal nx ny nz [--tol T] [--refs "c1,c2,..."] [--table out] [--regions file]
!   atomselect surface     <structure.nml> --regions <file> --normal nx ny nz [--tol T] [--refs "c1,..."]
!                          [--top N] [--bottom N] [--top-name top] [--bottom-name bottom]
!   atomselect validate    <structure.nml> <regions.nml> [--allow-no-fingerprint]
!   atomselect help
!
! Os arquivos de regiões (regions.nml) são criados no primeiro uso e estendidos
! depois; um arquivo existente é validado contra a estrutura antes de uma região
! ser adicionada. Os índices de átomos começam em 1; os comprimentos estão em Angstrom.
!
! Código de saída: 0 em sucesso (também para `help` ou sem argumentos); 1 em qualquer erro.
! Os erros são impressos como "[atomselect error] ..." na saída de erro padrão; os avisos
! ("[atomselect warning] ...") vão para a saída padrão.
!
! Documentação para o usuário: docs/atomselect.md
! ----------------------------------------------------------------------------------------------------------
program atomselect
    use, intrinsic :: iso_fortran_env, only: error_unit
    use types,            only: dp
    use structure_io,     only: structure_t, read_structure, structure_fingerprint
    use atomic_regions,   only: region_set_t, read_regions_file, write_regions_file, &
                                validate_regions, add_region, regions_are_partition
    use atomic_selection, only: parse_index_list, project_positions, select_by_range, &
                                select_by_box, detect_layers, check_periodic_split, &
                                assign_to_references, surface_layers
    implicit none

    type :: str_t
        character(len=:), allocatable :: s
    end type str_t

    integer,  parameter :: MSG_LEN = 1024
    integer,  parameter :: OPT_NAME_LEN = 32
    real(dp), parameter :: DEFAULT_LAYER_TOL = 0.5_dp   ! Å

    type(str_t),      allocatable :: args(:)     ! todos os argumentos da linha de comando
    type(str_t),      allocatable :: opt_key(:)  ! opções processadas (após os posicionais)
    type(str_t),      allocatable :: opt_val(:)
    integer                       :: nargs, nopt
    character(len=:), allocatable :: cmd

    ! 1. Coleta os argumentos (as opções são processadas por subcomando, ver parse_options)
    call collect_arguments()

    if (nargs < 1) then
        call print_usage()
        stop
    end if

    ! 2. Despacha conforme o subcomando
    cmd = args(1)%s
    select case (cmd)
    case ('help', '-h', '--help')
        call print_usage()
    case ('inspect')
        call cmd_inspect()
    case ('fingerprint')
        call cmd_fingerprint()
    case ('select')
        call cmd_select()
    case ('validate')
        call cmd_validate()
    case ('layers')
        call cmd_layers()
    case ('surface')
        call cmd_surface()
    case default
        call die("unknown subcommand '" // cmd // "' (run 'atomselect help')")
    end select

contains

    ! ==================================================================
    ! Subcomandos
    ! ==================================================================

    !> inspect: imprime o cabeçalho (célula, fingerprint) e a tabela de átomos; opcionalmente
    !> salva a tabela (--table) e um arquivo XYZ estendido (--xyz).
    subroutine cmd_inspect()
        type(structure_t) :: st
        character(len=16) :: fp
        integer           :: iu, ios

        call require_positionals(1, 'inspect <structure.nml>')
        call parse_options(3, [character(len=OPT_NAME_LEN) :: '--xyz', '--table'], &
                           [character(len=OPT_NAME_LEN) :: ''])
        call load_structure(args(2)%s, st)
        fp = structure_fingerprint(st)

        call write_header(6, args(2)%s, st, fp)
        call write_table(6, st)

        if (has_opt('--table')) then
            open(newunit=iu, file=opt_value('--table'), status='replace', action='write', iostat=ios)
            if (ios /= 0) call die("cannot open table output file '" // opt_value('--table') // "'")
            call write_header(iu, args(2)%s, st, fp)
            call write_table(iu, st)
            close(iu)
            print '(A)', 'table written to ' // opt_value('--table')
        end if

        if (has_opt('--xyz')) then
            open(newunit=iu, file=opt_value('--xyz'), status='replace', action='write', iostat=ios)
            if (ios /= 0) call die("cannot open xyz output file '" // opt_value('--xyz') // "'")
            call write_xyz(iu, st, fp)
            close(iu)
            print '(A)', 'xyz written to ' // opt_value('--xyz')
        end if
    end subroutine cmd_inspect

    ! ------------------------------------------------------------------
    !> fingerprint: imprime apenas o fingerprint da estrutura (vincula regiões/arquivos espectrais a uma estrutura).
    subroutine cmd_fingerprint()
        type(structure_t) :: st

        call require_positionals(1, 'fingerprint <structure.nml>')
        call parse_options(3, [character(len=OPT_NAME_LEN) :: ''], [character(len=OPT_NAME_LEN) :: ''])
        call load_structure(args(2)%s, st)
        print '(A)', structure_fingerprint(st)
    end subroutine cmd_fingerprint

    ! ------------------------------------------------------------------
    !> select: adiciona uma região nomeada a um arquivo de regiões. Exatamente um modo de seleção:
    !> lista explícita de índices (--atoms), intervalo de projeção ao longo de uma normal
    !> (--normal com --smin/--smax) ou caixa cartesiana (--box). Uma seleção
    !> geométrica vazia é um erro.
    subroutine cmd_select()
        type(structure_t)                   :: st
        type(region_set_t)                  :: rs
        integer,                allocatable :: list(:)
        real(dp),               allocatable :: s(:)
        real(dp)                            :: normal(3), smin, smax, lo(3), hi(3)
        character(len=:),       allocatable :: rfile, rname, how
        character(len=MSG_LEN)              :: errmsg
        integer                             :: ierr, nmode
        logical                             :: exists

        call require_positionals(1, 'select <structure.nml> --regions <regions.nml> --name NAME ' // &
                                    '(--atoms "LIST" | --normal nx ny nz [--smin A] [--smax B] | --box "...")')
        call parse_options(3, [character(len=OPT_NAME_LEN) :: '--regions', '--name', '--atoms', &
                                                              '--smin', '--smax', '--box'], &
                           [character(len=OPT_NAME_LEN) :: ''], &
                           [character(len=OPT_NAME_LEN) :: '--normal'])
        if (.not. has_opt('--regions')) call die('select: missing required option --regions')
        if (.not. has_opt('--name'))    call die('select: missing required option --name')
        nmode = count([has_opt('--atoms'), has_opt('--normal'), has_opt('--box')])
        if (nmode == 0) call die('select: one of --atoms, --normal or --box is required')
        if (nmode > 1)  call die('select: --atoms, --normal and --box are mutually exclusive')
        if ((has_opt('--smin') .or. has_opt('--smax')) .and. .not. has_opt('--normal')) &
            call die('select: --smin/--smax require --normal')
        rfile = opt_value('--regions')
        rname = opt_value('--name')

        call load_structure(args(2)%s, st)

        if (has_opt('--atoms')) then
            call parse_index_list(opt_value('--atoms'), st%nat, list, ierr, errmsg)
            if (ierr /= 0) call die('select: --atoms: ' // trim(errmsg))
            how = ''
        else if (has_opt('--normal')) then
            if (.not. (has_opt('--smin') .or. has_opt('--smax'))) &
                call die('select: --normal requires --smin and/or --smax')
            normal = opt_real3('--normal')
            smin = -huge(1.0_dp)
            smax = huge(1.0_dp)
            if (has_opt('--smin')) smin = opt_real('--smin')
            if (has_opt('--smax')) smax = opt_real('--smax')
            if (smin > smax) call die('select: --smin > --smax')
            call project_positions(st, normal, s, ierr, errmsg)
            if (ierr /= 0) call die('select: ' // trim(errmsg))
            call select_by_range(s, smin, smax, list)
            how = ' (projection range)'
        else
            call parse_box(opt_value('--box'), lo, hi)
            call select_by_box(st, lo, hi, list)
            how = ' (Cartesian box)'
        end if
        if (size(list) == 0) call die("select: geometric selection for region '" // rname // "' is empty")

        call open_regions(rfile, st, rs, exists)
        call add_region(rs, rname, list, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))
        call save_regions(rfile, rs)

        print '(A)', "region '" // rname // "': " // itoa(size(list)) // ' atoms [' // &
                     compress_list(list) // '] added to ' // rfile // status_word(exists) // how
    end subroutine cmd_select

    ! ------------------------------------------------------------------
    !> layers: detecta camadas ao longo de uma normal (ou atribui a --refs), imprime a
    !> tabela por átomo e o resumo das camadas; com --regions também armazena cada camada
    !> como uma região chamada layer_001, layer_002, ... (projeção crescente).
    subroutine cmd_layers()
        type(structure_t)                   :: st
        type(region_set_t)                  :: rs
        real(dp),               allocatable :: s(:), center(:)
        integer,                allocatable :: layer_of(:), list(:)
        real(dp)                            :: normal(3), tol
        character(len=MSG_LEN)              :: errmsg
        character(len=:),       allocatable :: rfile
        character(len=16)                   :: fp
        integer                             :: ierr, nlayer, l, iu, ios, i
        logical                             :: exists

        call require_positionals(1, 'layers <structure.nml> --normal nx ny nz [--tol T] [--refs "c1,..."] ' // &
                                    '[--table out] [--regions file]')
        call parse_options(3, [character(len=OPT_NAME_LEN) :: '--tol', '--refs', '--table', '--regions'], &
                           [character(len=OPT_NAME_LEN) :: ''], &
                           [character(len=OPT_NAME_LEN) :: '--normal'])
        if (.not. has_opt('--normal')) call die('layers: missing required option --normal')

        call load_structure(args(2)%s, st)
        call compute_layers(st, normal, tol, s, layer_of, nlayer, center)
        fp = structure_fingerprint(st)

        call write_header(6, args(2)%s, st, fp)
        call write_layer_info(6, normal, tol)
        call write_layer_table(6, st, s, layer_of)
        call write_layer_summary(6, layer_of, nlayer, center)

        if (has_opt('--table')) then
            open(newunit=iu, file=opt_value('--table'), status='replace', action='write', iostat=ios)
            if (ios /= 0) call die("cannot open table output file '" // opt_value('--table') // "'")
            call write_header(iu, args(2)%s, st, fp)
            call write_layer_info(iu, normal, tol)
            call write_layer_table(iu, st, s, layer_of)
            call write_layer_summary(iu, layer_of, nlayer, center)
            close(iu)
            print '(A)', 'table written to ' // opt_value('--table')
        end if

        if (has_opt('--regions')) then
            rfile = opt_value('--regions')
            call open_regions(rfile, st, rs, exists)
            do l = 1, nlayer
                list = pack([(i, i = 1, st%nat)], layer_of == l)
                if (size(list) == 0) call die('layers: layer ' // itoa(l) // &
                                              ' has no atoms; cannot create region ' // layer_name(l))
                call add_region(rs, layer_name(l), list, ierr, errmsg)
                if (ierr /= 0) call die(trim(errmsg))
            end do
            call save_regions(rfile, rs)
            print '(A)', itoa(nlayer) // ' regions ' // layer_name(1) // '..' // layer_name(nlayer) // &
                         ' added to ' // rfile // status_word(exists)
        end if
    end subroutine cmd_layers

    ! ------------------------------------------------------------------
    !> surface: constrói as regiões `top` (camadas de maior projeção) e/ou `bottom`
    !> (camadas de menor projeção) a partir das N camadas mais externas.
    subroutine cmd_surface()
        type(structure_t)                   :: st
        type(region_set_t)                  :: rs
        real(dp),               allocatable :: s(:), center(:)
        integer,                allocatable :: layer_of(:), top_list(:), bot_list(:)
        real(dp)                            :: normal(3), tol
        character(len=MSG_LEN)              :: errmsg
        character(len=:),       allocatable :: rfile, tname, bname
        integer                             :: ierr, nlayer, ntop, nbot
        logical                             :: exists

        call require_positionals(1, 'surface <structure.nml> --regions <file> --normal nx ny nz [--tol T] ' // &
                                    '[--refs "c1,..."] [--top N] [--bottom N] [--top-name top] [--bottom-name bottom]')
        call parse_options(3, [character(len=OPT_NAME_LEN) :: '--regions', '--tol', '--refs', '--top', &
                                                              '--bottom', '--top-name', '--bottom-name'], &
                           [character(len=OPT_NAME_LEN) :: ''], &
                           [character(len=OPT_NAME_LEN) :: '--normal'])
        if (.not. has_opt('--regions')) call die('surface: missing required option --regions')
        if (.not. has_opt('--normal'))  call die('surface: missing required option --normal')
        if (.not. (has_opt('--top') .or. has_opt('--bottom'))) &
            call die('surface: at least one of --top N or --bottom N is required')
        if (has_opt('--top-name') .and. .not. has_opt('--top')) call die('surface: --top-name requires --top')
        if (has_opt('--bottom-name') .and. .not. has_opt('--bottom')) &
            call die('surface: --bottom-name requires --bottom')
        ntop = 0
        nbot = 0
        if (has_opt('--top')) then
            ntop = opt_int('--top')
            if (ntop < 1) call die('surface: --top must be >= 1')
        end if
        if (has_opt('--bottom')) then
            nbot = opt_int('--bottom')
            if (nbot < 1) call die('surface: --bottom must be >= 1')
        end if
        tname = 'top'
        bname = 'bottom'
        if (has_opt('--top-name'))    tname = opt_value('--top-name')
        if (has_opt('--bottom-name')) bname = opt_value('--bottom-name')
        rfile = opt_value('--regions')

        call load_structure(args(2)%s, st)
        call compute_layers(st, normal, tol, s, layer_of, nlayer, center)
        call surface_layers(layer_of, nlayer, ntop, nbot, top_list, bot_list, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))
        if (ntop > 0 .and. size(top_list) == 0) call die('surface: top layers contain no atoms')
        if (nbot > 0 .and. size(bot_list) == 0) call die('surface: bottom layers contain no atoms')

        call open_regions(rfile, st, rs, exists)
        if (ntop > 0) then
            call add_region(rs, tname, top_list, ierr, errmsg)
            if (ierr /= 0) call die(trim(errmsg))
        end if
        if (nbot > 0) then
            call add_region(rs, bname, bot_list, ierr, errmsg)
            if (ierr /= 0) call die(trim(errmsg))
        end if
        call save_regions(rfile, rs)

        print '(A)', 'layers found: ' // itoa(nlayer)
        if (ntop > 0) print '(A)', "region '" // tname // "': layers " // itoa(nlayer - ntop + 1) // '-' // &
                                   itoa(nlayer) // ', ' // itoa(size(top_list)) // ' atoms [' // &
                                   compress_list(top_list) // '] added to ' // rfile // status_word(exists)
        if (nbot > 0) print '(A)', "region '" // bname // "': layers 1-" // itoa(nbot) // ', ' // &
                                   itoa(size(bot_list)) // ' atoms [' // compress_list(bot_list) // &
                                   '] added to ' // rfile // status_word(exists)
    end subroutine cmd_surface

    ! ------------------------------------------------------------------
    !> Compartilhada por layers/surface: lê --normal/--tol/--refs, projeta, avisa
    !> sobre divisão periódica e atribui as camadas (detectadas ou por referências).
    subroutine compute_layers(st, normal, tol, s, layer_of, nlayer, center)
        type(structure_t),      intent(in)               :: st
        real(dp),               intent(out)              :: normal(3), tol
        real(dp),               allocatable, intent(out) :: s(:), center(:)
        integer,                allocatable, intent(out) :: layer_of(:)
        integer,                intent(out)              :: nlayer
        real(dp),               allocatable              :: refs(:)
        character(len=MSG_LEN)                           :: errmsg
        integer                                          :: ierr
        logical                                          :: warn

        normal = opt_real3('--normal')
        tol = DEFAULT_LAYER_TOL
        if (has_opt('--tol')) tol = opt_real('--tol')
        if (.not. (tol > 0.0_dp)) call die('--tol must be > 0')

        call project_positions(st, normal, s, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))

        call check_periodic_split(st, normal, s, tol, warn)
        if (warn) then
            print '(A)', '[atomselect warning] the slab appears to cross the periodic cell boundary along ' // &
                         'the normal (projected coordinates are not wrapped); detected layers may be split.'
            print '(A)', '[atomselect warning] consider explicit layer centers with --refs, or translate ' // &
                         'the structure so the slab is not cut by the boundary.'
        end if

        if (has_opt('--refs')) then
            call parse_real_list(opt_value('--refs'), refs)
            call assign_to_references(s, refs, tol, layer_of, ierr, errmsg)
            if (ierr /= 0) call die(trim(errmsg))
            nlayer = size(refs)
            center = sorted_real(refs)       ! o resumo reporta os centros de referência
        else
            call detect_layers(s, tol, layer_of, nlayer, center, ierr, errmsg)
            if (ierr /= 0) call die(trim(errmsg))
        end if
    end subroutine compute_layers

    ! ------------------------------------------------------------------
    !> validate: confere um arquivo de regiões contra uma estrutura (fingerprint, contagens/índices
    !> de átomos) e reporta sobreposições, átomos não cobertos e se as regiões formam
    !> uma partição. Chegar à última linha significa "válido".
    subroutine cmd_validate()
        type(structure_t)                   :: st
        type(region_set_t)                  :: rs
        character(len=MSG_LEN)              :: errmsg
        integer,                allocatable :: cnt(:)
        logical,                allocatable :: in_i(:)
        integer                             :: ierr, ir, jr, a, nshared, noverlap, nuncovered

        call require_positionals(2, 'validate <structure.nml> <regions.nml>')
        call parse_options(4, [character(len=OPT_NAME_LEN) :: ''], &
                           [character(len=OPT_NAME_LEN) :: '--allow-no-fingerprint'])

        call load_structure(args(2)%s, st)
        call read_regions_file(args(3)%s, rs, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))
        call validate_regions(rs, st, has_opt('--allow-no-fingerprint'), ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))

        print '(A)', 'structure   : ' // args(2)%s
        print '(A)', 'regions     : ' // args(3)%s
        print '(A)', 'nat         : ' // itoa(st%nat)
        if (len_trim(rs%fingerprint) > 0) then
            print '(A)', 'fingerprint : ' // trim(rs%fingerprint) // ' (matches)'
        else
            print '(A)', 'fingerprint : none'
        end if
        print '(A)', 'nregion     : ' // itoa(rs%nregion)

        allocate(cnt(st%nat), source=0)
        allocate(in_i(st%nat))
        do ir = 1, rs%nregion
            print '(A)', "region '" // trim(rs%regions(ir)%name) // "': natoms = " // &
                         itoa(size(rs%regions(ir)%atoms)) // ' atoms = [' // &
                         compress_list(sorted(rs%regions(ir)%atoms)) // ']'
            do a = 1, size(rs%regions(ir)%atoms)
                cnt(rs%regions(ir)%atoms(a)) = cnt(rs%regions(ir)%atoms(a)) + 1
            end do
        end do

        noverlap = 0
        do ir = 1, rs%nregion
            in_i = .false.
            in_i(rs%regions(ir)%atoms) = .true.
            do jr = ir + 1, rs%nregion
                nshared = count(in_i(rs%regions(jr)%atoms))
                if (nshared > 0) then
                    noverlap = noverlap + 1
                    print '(A)', "overlap: '" // trim(rs%regions(ir)%name) // "' and '" // &
                                 trim(rs%regions(jr)%name) // "' share " // itoa(nshared) // ' atoms'
                end if
            end do
        end do
        if (noverlap == 0) print '(A)', 'overlaps    : none'

        nuncovered = count(cnt == 0)
        print '(A)', 'uncovered   : ' // itoa(nuncovered) // ' atoms'
        if (regions_are_partition(rs)) then
            print '(A)', 'partition   : yes'
        else
            print '(A)', 'partition   : no'
        end if
        print '(A)', 'valid       : yes'
    end subroutine cmd_validate

    ! ==================================================================
    ! Auxiliares do arquivo de regiões (lógica compartilhada de criar/acrescentar, contract §6)
    ! ==================================================================

    !> Lê e valida um arquivo de regiões existente, ou inicia um novo com o
    !> cabeçalho tirado da estrutura.
    subroutine open_regions(rfile, st, rs, exists)
        character(len=*),       intent(in)  :: rfile
        type(structure_t),      intent(in)  :: st
        type(region_set_t),     intent(out) :: rs
        logical,                intent(out) :: exists
        character(len=MSG_LEN)              :: errmsg
        integer                             :: ierr

        inquire(file=rfile, exist=exists)
        if (exists) then
            call read_regions_file(rfile, rs, ierr, errmsg)
            if (ierr /= 0) call die(trim(errmsg))
            call validate_regions(rs, st, .false., ierr, errmsg)
            if (ierr /= 0) call die(trim(errmsg))
        else
            rs%format_version = 1
            rs%natoms = st%nat
            rs%fingerprint = structure_fingerprint(st)
            rs%has_species = .true.
            rs%species = st%label
            rs%nregion = 0
            allocate(rs%regions(0))
        end if
    end subroutine open_regions

    !> Grava o conjunto de regiões em disco; encerra com erro em caso de falha de E/S.
    subroutine save_regions(rfile, rs)
        character(len=*),       intent(in) :: rfile
        type(region_set_t),     intent(in) :: rs
        character(len=MSG_LEN)             :: errmsg
        integer                            :: ierr
        call write_regions_file(rfile, rs, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))
    end subroutine save_regions

    !> Sufixo da linha de log: indica se o arquivo de regiões foi criado ou atualizado.
    function status_word(exists) result(w)
        logical,          intent(in)  :: exists
        character(len=:), allocatable :: w
        if (exists) then
            w = ' (updated)'
        else
            w = ' (created)'
        end if
    end function status_word

    !> Nome da região da camada l: layer_001, layer_002, ...
    function layer_name(l) result(name)
        integer,           intent(in)  :: l
        character(len=:),  allocatable :: name
        character(len=16)              :: buf
        write(buf, '(I0.3)') l
        name = 'layer_' // trim(buf)
    end function layer_name

    ! ==================================================================
    ! Saída das camadas
    ! ==================================================================

    !> Linhas de comentário prefixadas com '#' descrevendo como as camadas foram obtidas.
    subroutine write_layer_info(iu, normal, tol)
        integer,  intent(in) :: iu
        real(dp), intent(in) :: normal(3), tol
        write(iu, '(A,3F12.6)') '# normal (n)  : ', normal / norm2(normal)
        write(iu, '(A,F12.6)')  '# tol (A)     : ', tol
        if (has_opt('--refs')) then
            write(iu, '(A)') '# layers      : reference centers (--refs)'
        else
            write(iu, '(A)') '# layers      : detected (single-linkage, ascending projection)'
        end if
    end subroutine write_layer_info

    !> Tabela por átomo: índice, elemento, posição cartesiana, projeção s, camada.
    subroutine write_layer_table(iu, st, s, layer_of)
        integer,           intent(in) :: iu
        type(structure_t), intent(in) :: st
        real(dp),          intent(in) :: s(:)
        integer,           intent(in) :: layer_of(:)
        integer                       :: i

        write(iu, '(A6,2X,A8,4A12,A8)') 'Index', 'Element', 'X', 'Y', 'Z', 'Projected', 'Layer'
        do i = 1, st%nat
            write(iu, '(I6,2X,A8,4F12.6,I8)') i, adjustl(st%label(i)), st%cart(:, i), s(i), layer_of(i)
        end do
    end subroutine write_layer_table

    !> Centro, contagem de átomos e lista comprimida de átomos por camada (prefixada com '#').
    subroutine write_layer_summary(iu, layer_of, nlayer, center)
        integer,  intent(in)  :: iu
        integer,  intent(in)  :: layer_of(:), nlayer
        real(dp), intent(in)  :: center(:)
        integer,  allocatable :: list(:)
        integer               :: l, i

        write(iu, '(A)') ''
        write(iu, '(A)') '# nlayer      : ' // itoa(nlayer)
        write(iu, '(A1,A5,A12,A8,2X,A)') '#', 'Layer', 'Center', 'Count', 'Atoms'
        do l = 1, nlayer
            list = pack([(i, i = 1, size(layer_of))], layer_of == l)
            write(iu, '(A1,I5,F12.6,I8,2X,A)') '#', l, center(l), size(list), compress_list(list)
        end do
    end subroutine write_layer_summary

    ! ==================================================================
    ! Auxiliares de saída
    ! ==================================================================

    !> Cabeçalho prefixado com '#': nome do arquivo, número de átomos, vetores de rede, fingerprint.
    subroutine write_header(iu, filename, st, fp)
        integer,           intent(in) :: iu
        character(len=*),  intent(in) :: filename
        type(structure_t), intent(in) :: st
        character(len=16), intent(in) :: fp
        integer                       :: i

        write(iu, '(A)') '# file        : ' // filename
        write(iu, '(A)') '# nat         : ' // itoa(st%nat)
        do i = 1, 3
            write(iu, '(A,I1,A,3F12.6)') '# a', i, ' (A)      : ', st%lattice(i, :)
        end do
        write(iu, '(A)') '# fingerprint : ' // fp
    end subroutine write_header

    subroutine write_table(iu, st)
        integer,           intent(in) :: iu
        type(structure_t), intent(in) :: st
        integer                       :: i

        write(iu, '(A6,2X,A8,3A12)') 'Index', 'Element', 'X', 'Y', 'Z'
        do i = 1, st%nat
            write(iu, '(I6,2X,A8,3F12.6)') i, adjustl(st%label(i)), st%cart(:, i)
        end do
    end subroutine write_table

    !> XYZ estendido (rede + fingerprint na linha de comentário), para visualização.
    subroutine write_xyz(iu, st, fp)
        integer,           intent(in)  :: iu
        type(structure_t), intent(in)  :: st
        character(len=16), intent(in)  :: fp
        character(len=:),  allocatable :: lat
        character(len=24)              :: buf
        integer                        :: i, j

        lat = ''
        do i = 1, 3
            do j = 1, 3
                write(buf, '(F12.6)') st%lattice(i, j)
                if (len(lat) > 0) lat = lat // ' '
                lat = lat // trim(adjustl(buf))
            end do
        end do
        write(iu, '(I0)') st%nat
        write(iu, '(A)') 'Lattice="' // lat // '" fingerprint=' // fp
        do i = 1, st%nat
            write(iu, '(A,3(1X,F12.6))') trim(adjustl(st%label(i))), st%cart(:, i)
        end do
    end subroutine write_xyz

    !> Texto de uso; manter em sincronia com a sinopse do cabeçalho do arquivo.
    subroutine print_usage()
        print '(A)', 'atomselect - inspect structures and manage named atomic regions (regions.nml)'
        print '(A)', ''
        print '(A)', 'Usage:'
        print '(A)', '  atomselect inspect     <structure.nml> [--xyz out.xyz] [--table out.txt]'
        print '(A)', '  atomselect fingerprint <structure.nml>'
        print '(A)', '  atomselect select      <structure.nml> --regions <regions.nml> --name NAME --atoms "LIST"'
        print '(A)', '  atomselect select      <structure.nml> --regions <regions.nml> --name NAME ' // &
                     '--normal nx ny nz [--smin A] [--smax B]'
        print '(A)', '  atomselect select      <structure.nml> --regions <regions.nml> --name NAME ' // &
                     '--box "xlo xhi ylo yhi zlo zhi"'
        print '(A)', '  atomselect layers      <structure.nml> --normal nx ny nz [--tol T] [--refs "c1,c2,..."] ' // &
                     '[--table out] [--regions file]'
        print '(A)', '  atomselect surface     <structure.nml> --regions <file> --normal nx ny nz [--tol T] ' // &
                     '[--refs "c1,..."] [--top N] [--bottom N] [--top-name top] [--bottom-name bottom]'
        print '(A)', '  atomselect validate    <structure.nml> <regions.nml> [--allow-no-fingerprint]'
        print '(A)', '  atomselect help'
        print '(A)', ''
        print '(A)', 'LIST: comma-separated 1-based indices and inclusive ranges, e.g. "1-4, 7, 10-12".'
        print '(A)', 'Normal: three CARTESIAN components (any direction, not assumed z); projection s = r.n/|n| (A).'
        print '(A)', 'Layers: atoms sorted by s, new layer when the gap exceeds --tol (default 0.5 A); numbered'
        print '(A)', '        by ascending s ("top" = highest s). --refs assigns atoms to the nearest given center.'
        print '(A)', '        Projections are not wrapped: the slab must not be cut by the cell boundary.'
        print '(A)', 'Box: six limits in A; use * for an unbounded limit.'
    end subroutine print_usage

    ! ==================================================================
    ! Tratamento de argumentos
    ! ==================================================================

    !> Guarda todos os argumentos da linha de comando em args(:) e reinicia o armazenamento de opções.
    subroutine collect_arguments()
        integer :: i, n

        nargs = command_argument_count()
        allocate(args(max(nargs, 0)))
        do i = 1, nargs
            call get_command_argument(i, length=n)
            allocate(character(len=n) :: args(i)%s)
            if (n > 0) call get_command_argument(i, args(i)%s)
        end do
        nopt = 0
        allocate(opt_key(0), opt_val(0))
    end subroutine collect_arguments

    !> Exige `npos` argumentos posicionais após o subcomando (args(2:npos+1)),
    !> nenhum dos quais pode parecer uma opção.
    subroutine require_positionals(npos, synopsis)
        integer,          intent(in) :: npos
        character(len=*), intent(in) :: synopsis
        integer                      :: i

        if (nargs < npos + 1) call die('missing arguments; usage: atomselect ' // synopsis)
        do i = 2, npos + 1
            if (is_option(args(i)%s)) call die('missing arguments; usage: atomselect ' // synopsis)
        end do
    end subroutine require_positionals

    !> Processa args(istart:) como opções. `valued` lista as opções que recebem um valor,
    !> `flags` lista as chaves sem valor (entradas em branco são ignoradas). Qualquer outro
    !> argumento, ou uma opção repetida, é um erro.
    subroutine parse_options(istart, valued, flags, triples)
        integer,          intent(in)           :: istart
        character(len=*), intent(in)           :: valued(:), flags(:)
        character(len=*), intent(in), optional :: triples(:)   ! opções que recebem 3 tokens
        integer                                :: i
        character(len=:), allocatable          :: a

        i = istart
        do while (i <= nargs)
            a = args(i)%s
            if (.not. is_option(a)) call die("unexpected argument '" // a // "'")
            if (has_opt(a)) call die("option '" // a // "' given more than once")
            if (in_list(a, valued)) then
                if (i + 1 > nargs) call die("option '" // a // "' requires a value")
                call push_opt(a, args(i + 1)%s)
                i = i + 2
            else if (in_list(a, flags)) then
                call push_opt(a, '')
                i = i + 1
            else if (present(triples)) then
                if (.not. in_list(a, triples)) &
                    call die("unknown option '" // a // "' for subcommand '" // cmd // "'")
                if (i + 3 > nargs) call die("option '" // a // "' requires three values")
                call push_opt(a, args(i + 1)%s // ' ' // args(i + 2)%s // ' ' // args(i + 3)%s)
                i = i + 4
            else
                call die("unknown option '" // a // "' for subcommand '" // cmd // "'")
            end if
        end do
    end subroutine parse_options

    !> Acrescenta uma opção processada (chave/valor) ao armazenamento de opções.
    subroutine push_opt(key, val)
        character(len=*), intent(in)  :: key, val
        type(str_t),      allocatable :: tk(:), tv(:)

        allocate(tk(nopt + 1), tv(nopt + 1))
        tk(1:nopt) = opt_key
        tv(1:nopt) = opt_val
        tk(nopt + 1)%s = key
        tv(nopt + 1)%s = val
        call move_alloc(tk, opt_key)
        call move_alloc(tv, opt_val)
        nopt = nopt + 1
    end subroutine push_opt

    logical function has_opt(key)
        character(len=*), intent(in) :: key
        integer                      :: i
        has_opt = .false.
        do i = 1, nopt
            if (opt_key(i)%s == key) then
                has_opt = .true.
                return
            end if
        end do
    end function has_opt

    function opt_value(key) result(val)
        character(len=*), intent(in)  :: key
        character(len=:), allocatable :: val
        integer                       :: i
        val = ''
        do i = 1, nopt
            if (opt_key(i)%s == key) then
                val = opt_val(i)%s
                return
            end if
        end do
    end function opt_value

    logical function in_list(a, list)
        character(len=*), intent(in) :: a, list(:)
        integer                      :: i
        in_list = .false.
        do i = 1, size(list)
            if (len_trim(list(i)) == 0) cycle
            if (a == trim(list(i))) then
                in_list = .true.
                return
            end if
        end do
    end function in_list

    logical function is_option(a)
        character(len=*), intent(in) :: a
        is_option = .false.
        if (len(a) >= 2) is_option = (a(1:2) == '--')
    end function is_option

    ! ==================================================================
    ! Auxiliares diversos
    ! ==================================================================

    !> Converte um token em número real; encerra com uma mensagem que nomeia a opção.
    function to_real(tok, what) result(x)
        character(len=*), intent(in) :: tok, what
        real(dp)                     :: x
        integer                      :: ios
        x = 0.0_dp
        if (len_trim(tok) == 0 .or. scan(trim(adjustl(tok)), ' ,/') > 0) &
            call die(what // ": not a number: '" // tok // "'")
        read(tok, *, iostat=ios) x
        if (ios /= 0) call die(what // ": not a number: '" // tok // "'")
        if (x /= x .or. abs(x) > huge(1.0_dp)) call die(what // ": not a finite number: '" // tok // "'")
    end function to_real

    !> Valor da opção `key` como real finito.
    function opt_real(key) result(x)
        character(len=*), intent(in) :: key
        real(dp)                     :: x
        x = to_real(opt_value(key), key)
    end function opt_real

    !> Valor da opção `key` como inteiro não negativo (no máximo 9 dígitos).
    function opt_int(key) result(n)
        character(len=*), intent(in)  :: key
        integer                       :: n
        character(len=:), allocatable :: v
        integer                       :: ios
        v = trim(adjustl(opt_value(key)))
        n = 0
        if (len(v) == 0 .or. len(v) > 9 .or. verify(v, '0123456789') /= 0) &
            call die(key // ": not a non-negative integer: '" // v // "'")
        read(v, *, iostat=ios) n
        if (ios /= 0) call die(key // ": not an integer: '" // v // "'")
    end function opt_int

    !> Três números separados por espaço, guardados por parse_options para uma opção tripla.
    function opt_real3(key) result(x)
        character(len=*), intent(in) :: key
        real(dp)                     :: x(3)
        type(str_t)                  :: tok(6)
        integer                      :: ntok, k
        call split_ws(opt_value(key), tok, ntok)
        if (ntok /= 3) call die(key // ': expected three numbers')
        do k = 1, 3
            x(k) = to_real(tok(k)%s, key)
        end do
    end function opt_real3

    !> "c1,c2,..." -> lista de reais (ao menos um valor).
    subroutine parse_real_list(text, vals)
        character(len=*), intent(in)               :: text
        real(dp),         allocatable, intent(out) :: vals(:)
        integer                                    :: pos, next
        character(len=:), allocatable              :: item

        allocate(vals(0))
        if (len_trim(text) == 0) call die('--refs: empty list')
        pos = 1
        do
            next = index(text(pos:), ',')
            if (next == 0) then
                item = trim(adjustl(text(pos:)))
            else
                item = trim(adjustl(text(pos:pos + next - 2)))
            end if
            if (len(item) == 0) call die('--refs: empty item')
            vals = [vals, to_real(item, '--refs')]
            if (next == 0) exit
            pos = pos + next
            if (pos > len(text)) call die('--refs: empty item')
        end do
    end subroutine parse_real_list

    !> "xlo xhi ylo yhi zlo zhi" com '*' para limites sem restrição.
    subroutine parse_box(text, lo, hi)
        character(len=*), intent(in)  :: text
        real(dp),         intent(out) :: lo(3), hi(3)
        type(str_t)                   :: tok(7)
        integer                       :: ntok, k
        real(dp)                      :: v(6)

        call split_ws(text, tok, ntok)
        if (ntok /= 6) call die('--box: expected six limits "xlo xhi ylo yhi zlo zhi" (use * for unbounded)')
        do k = 1, 6
            if (tok(k)%s == '*') then
                if (mod(k, 2) == 1) then
                    v(k) = -huge(1.0_dp)
                else
                    v(k) = huge(1.0_dp)
                end if
            else
                v(k) = to_real(tok(k)%s, '--box')
            end if
        end do
        lo = v([1, 3, 5])
        hi = v([2, 4, 6])
        if (any(lo > hi)) call die('--box: lower limit greater than upper limit')
    end subroutine parse_box

    !> Divide em espaços; ntok pode exceder size(tok) (tokens extras são contados, não armazenados).
    subroutine split_ws(text, tok, ntok)
        character(len=*), intent(in)    :: text
        type(str_t),      intent(inout) :: tok(:)
        integer,          intent(out)   :: ntok
        integer                         :: i, j, n
        n = len(text)
        ntok = 0
        i = 1
        do while (i <= n)
            if (text(i:i) == ' ' .or. text(i:i) == achar(9)) then
                i = i + 1
                cycle
            end if
            j = i
            do while (j < n)
                if (text(j + 1:j + 1) == ' ' .or. text(j + 1:j + 1) == achar(9)) exit
                j = j + 1
            end do
            ntok = ntok + 1
            if (ntok <= size(tok)) tok(ntok)%s = text(i:j)
            i = j + 1
        end do
    end subroutine split_ws

    !> Cópia em ordem crescente (insertion sort; usada para a pequena lista --refs).
    function sorted_real(v) result(r)
        real(dp), intent(in)  :: v(:)
        real(dp), allocatable :: r(:)
        real(dp)              :: x
        integer               :: i, j
        r = v
        do i = 2, size(r)
            x = r(i)
            j = i - 1
            do while (j >= 1)
                if (r(j) <= x) exit
                r(j + 1) = r(j)
                j = j - 1
            end do
            r(j + 1) = x
        end do
    end function sorted_real

    !> Lê o arquivo de estrutura; em caso de falha, encerra com a mensagem do leitor.
    subroutine load_structure(filename, st)
        character(len=*),       intent(in)  :: filename
        type(structure_t),      intent(out) :: st
        character(len=MSG_LEN)              :: errmsg
        integer                             :: ierr

        call read_structure(filename, st, ierr, errmsg)
        if (ierr /= 0) call die(trim(errmsg))
    end subroutine load_structure

    !> Comprime uma lista crescente de índices em intervalos: [1,2,3,5] -> "1-3,5".
    function compress_list(list) result(s)
        integer,          intent(in)  :: list(:)
        character(len=:), allocatable :: s
        integer                       :: i, j

        s = ''
        i = 1
        do while (i <= size(list))
            j = i
            do while (j < size(list))
                if (list(j + 1) /= list(j) + 1) exit
                j = j + 1
            end do
            if (len(s) > 0) s = s // ','
            if (j > i) then
                s = s // itoa(list(i)) // '-' // itoa(list(j))
            else
                s = s // itoa(list(i))
            end if
            i = j + 1
        end do
    end function compress_list

    !> Retorna uma cópia em ordem crescente (insertion sort; as listas são pequenas e quase ordenadas).
    function sorted(list) result(r)
        integer, intent(in)  :: list(:)
        integer, allocatable :: r(:)
        integer              :: i, j, v

        r = list
        do i = 2, size(r)
            v = r(i)
            j = i - 1
            do while (j >= 1)
                if (r(j) <= v) exit
                r(j + 1) = r(j)
                j = j - 1
            end do
            r(j + 1) = v
        end do
    end function sorted

    pure function itoa(i) result(s)
        integer,           intent(in)  :: i
        character(len=:),  allocatable :: s
        character(len=16)              :: buf
        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

    !> Imprime o erro e encerra com código 1.
    subroutine die(msg)
        character(len=*), intent(in) :: msg
        write(error_unit, '(A)') '[atomselect error] ' // msg
        stop 1
    end subroutine die

end program atomselect
