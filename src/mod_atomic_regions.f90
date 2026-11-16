! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Regiões atômicas nomeadas (regions.nml): leitor, escritor canônico, checagens
!> intrínsecas, validação contra uma estrutura e teste de partição.
!>
!> Papel no pipeline: atomselect escreve regions.nml (grupos nomeados de átomos,
!> p. ex. 'top', 'bottom'); spectral/specmap o leem para projetar a função
!> espectral sobre cada região. O arquivo é um namelist: um cabeçalho
!> &REGION_FILE seguido de qualquer número de grupos &REGION.
!>
!> Convenções: índices de átomos começam em 1 e seguem a ordem do arquivo de
!> estrutura (= ordem HSIO). Regiões podem se sobrepor; as listas de átomos
!> guardadas mantêm a ordem do usuário. As checagens têm dois níveis:
!> read_regions_file / add_region verificam o que se pode checar só com o
!> arquivo (sintaxe, faixas, nomes, duplicatas) e validate_regions compara com
!> uma estrutura (nat, rótulos, fingerprint).
!>
!> As rotinas da biblioteca nunca param o programa: erros voltam como
!> ierr /= 0 mais errmsg. Ver docs/design/contracts.md §2.
module atomic_regions
    use structure_io, only: structure_t, structure_fingerprint, LABEL_LEN
    use, intrinsic :: iso_fortran_env, only: iostat_end
    implicit none
    private

    ! ==================================================================
    ! TIPOS
    ! ==================================================================
    !> REGION_NAME_LEN: comprimento máximo do nome da região. MAX_REGION_ATOMS:
    !> tamanho dos buffers de leitura do namelist, ou seja, a maior lista de
    !> átomos (ou espécies) que um arquivo pode conter.
    integer, parameter, public :: REGION_NAME_LEN = 64, MAX_REGION_ATOMS = 100000

    !> Uma região nomeada: `atoms` são índices da estrutura (base 1), sem duplicatas.
    type, public :: atom_region_t
        character(len=REGION_NAME_LEN) :: name = ''
        integer, allocatable :: atoms(:)
    end type atom_region_t

    !> Conteúdo de um arquivo de regiões. `natoms`, `fingerprint` e `species`
    !> descrevem a estrutura para a qual as regiões foram construídas (checada
    !> depois por validate_regions). fingerprint = '' significa "nenhum".
    !> regions(1:nregion) guarda os dados; nregion é mantido por add_region.
    type, public :: region_set_t
        integer :: format_version = 1
        integer :: natoms = 0
        character(len=16) :: fingerprint = ''
        logical :: has_species = .false.
        character(len=LABEL_LEN), allocatable :: species(:)   ! (natoms) se has_species
        integer :: nregion = 0
        type(atom_region_t), allocatable :: regions(:)
    end type region_set_t

    ! ==================================================================
    ! API PÚBLICA / CONSTANTES
    ! ==================================================================
    public :: read_regions_file, write_regions_file, validate_regions
    public :: add_region, regions_are_partition, is_valid_region_name

    ! Marca posições não usadas nos buffers do namelist (um índice real é sempre >= 1)
    integer, parameter :: ATOM_SENTINEL = -huge(0)
    ! Quebra de linha do escritor canônico
    integer, parameter :: SPECIES_PER_LINE = 10, ATOMS_PER_LINE = 20
    character(len=*), parameter :: NAME_CHARS = &
        'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.'

contains

    ! ==================================================================
    ! LEITURA / ESCRITA
    ! ==================================================================

    !> Lê um arquivo de regiões em `rs` e roda as checagens que não precisam de
    !> estrutura: &REGION_FILE primeiro e único, format_version = 1, natoms >= 1,
    !> fingerprint <= 16 caracteres, número de species = natoms (se dado) e, por
    !> região, as checagens de add_region (nome, duplicatas, faixa de índices).
    !> Nunca para o programa; em erro ierr = 1 e errmsg descreve o problema
    !> (inclui número do bloco / nome da região). rs é intent(out): o conteúdo
    !> não é confiável após um erro.
    subroutine read_regions_file(filename, rs, ierr, errmsg)
        character(len=*),   intent(in)  :: filename
        type(region_set_t), intent(out) :: rs
        integer,            intent(out) :: ierr
        character(len=*),   intent(out) :: errmsg

        ! variáveis do namelist (buffers alocados antes de cada leitura)
        integer                               :: format_version, natoms
        character(len=64)                     :: fingerprint
        character(len=LABEL_LEN), allocatable :: species(:)
        character(len=256)                    :: name
        integer,                  allocatable :: atoms(:)
        namelist /REGION_FILE/ format_version, natoms, fingerprint, species
        namelist /REGION/ name, atoms

        integer            :: iu, ios, n, iblock
        character(len=512) :: msg
        logical            :: ok

        ierr = 0
        errmsg = ''

        ! 1. Pré-varredura: a leitura de namelist não informa onde está o cabeçalho
        call check_header_position(filename, ok, ierr, errmsg)
        if (ierr /= 0) return

        open(newunit=iu, file=filename, status='old', action='read', iostat=ios, iomsg=msg)
        if (ios /= 0) then
            call fail('cannot open regions file "' // trim(filename) // '": ' // trim(msg))
            return
        end if

        allocate(species(MAX_REGION_ATOMS), atoms(MAX_REGION_ATOMS))

        ! 2. Cabeçalho &REGION_FILE. Os valores padrão são sentinelas, para que uma
        !    chave ausente seja detectada (natoms < 1, format_version /= 1).
        format_version = -huge(0)
        natoms = -huge(0)
        fingerprint = ''
        species = ''
        read(iu, nml=REGION_FILE, iostat=ios, iomsg=msg)
        if (ios /= 0) then
            call fail('regions file "' // trim(filename) // '": error reading &REGION_FILE: ' // trim(msg))
            close(iu)
            return
        end if
        if (format_version /= 1) then
            call fail('regions file: format_version must be 1 (got ' // itoa(format_version) // ')')
            close(iu)
            return
        end if
        if (natoms < 1) then
            call fail('regions file: natoms must be >= 1 (got ' // itoa(natoms) // ')')
            close(iu)
            return
        end if
        if (len_trim(fingerprint) > len(rs%fingerprint)) then
            call fail('regions file: fingerprint longer than 16 characters: "' // trim(fingerprint) // '"')
            close(iu)
            return
        end if
        rs%format_version = format_version
        rs%natoms = natoms
        rs%fingerprint = fingerprint(1:len(rs%fingerprint))

        n = count_leading_species(species)
        if (n < 0) then
            call fail('regions file: gap in species list')
            close(iu)
            return
        end if
        if (n > 0) then
            if (n /= natoms) then
                call fail('regions file: species has ' // itoa(n) // ' entries, natoms = ' // itoa(natoms))
                close(iu)
                return
            end if
            rs%has_species = .true.
            rs%species = species(1:n)
        end if

        ! 3. Blocos &REGION, lidos em sequência até o fim do arquivo. O buffer de
        !    átomos volta à sentinela antes de cada leitura, de modo que o tamanho
        !    da lista é o número de entradas iniciais não sentinela; add_region
        !    então aplica as checagens por região.
        allocate(rs%regions(0))
        rs%nregion = 0
        iblock = 0
        do
            iblock = iblock + 1
            name = ''
            atoms = ATOM_SENTINEL
            read(iu, nml=REGION, iostat=ios, iomsg=msg)
            if (ios == iostat_end) exit
            if (ios /= 0) then
                call fail('regions file: error reading &REGION block ' // itoa(iblock) // ': ' // trim(msg))
                close(iu)
                return
            end if
            if (len_trim(name) == 0) then
                call fail('regions file: &REGION block ' // itoa(iblock) // ' has no name')
                close(iu)
                return
            end if
            n = count_leading_atoms(atoms)
            if (n < 0) then
                call fail('regions file: region "' // trim(name) // '" (block ' // itoa(iblock) // &
                          '): gap in atoms list')
                close(iu)
                return
            end if
            call add_region(rs, trim(name), atoms(1:n), ierr, errmsg)
            if (ierr /= 0) then
                errmsg = 'regions file: block ' // itoa(iblock) // ': ' // trim(errmsg)
                close(iu)
                return
            end if
        end do
        close(iu)

    contains

        !> Define o estado de erro (o arquivo é fechado pelo chamador).
        subroutine fail(m)
            character(len=*), intent(in) :: m
            ierr = 1
            errmsg = m
        end subroutine fail

    end subroutine read_regions_file

    !> Pré-varredura de texto do arquivo: &REGION_FILE deve ser o primeiro grupo
    !> do namelist e aparecer exatamente uma vez (uma leitura de namelist sozinha
    !> pularia em silêncio até o primeiro grupo com o nome pedido). `ok` é
    !> verdadeiro em caso de sucesso; erros vão por ierr/errmsg. Grupos são
    !> reconhecidos por linhas cujo primeiro caractere não branco é '&'.
    subroutine check_header_position(filename, ok, ierr, errmsg)
        character(len=*), intent(in)  :: filename
        logical,          intent(out) :: ok
        integer,          intent(out) :: ierr
        character(len=*), intent(out) :: errmsg
        integer             :: iu, ios, nheader, ngroup
        character(len=1024) :: line, tok
        character(len=512)  :: msg

        ok = .false.
        ierr = 0
        errmsg = ''
        open(newunit=iu, file=filename, status='old', action='read', iostat=ios, iomsg=msg)
        if (ios /= 0) then
            ierr = 1
            errmsg = 'cannot open regions file "' // trim(filename) // '": ' // trim(msg)
            return
        end if
        nheader = 0
        ngroup = 0
        do
            read(iu, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (line(1:1) /= '&') cycle
            tok = upper(first_token(line))
            ngroup = ngroup + 1
            if (trim(tok) == '&REGION_FILE') then
                nheader = nheader + 1
                if (ngroup /= 1) then
                    ierr = 1
                    errmsg = 'regions file: &REGION_FILE must be the first namelist group'
                    exit
                end if
            end if
        end do
        close(iu)
        if (ierr /= 0) return
        if (nheader /= 1) then
            ierr = 1
            errmsg = 'regions file: expected exactly one &REGION_FILE group, found ' // itoa(nheader)
            return
        end if
        ok = .true.
    end subroutine check_header_position

    !> Escreve `rs` no layout canônico (cabeçalho, depois um &REGION por região,
    !> 20 átomos por linha, 10 species por linha). Sobrescreve um arquivo
    !> existente. Não valida `rs` (use add_region/validate_regions para isso).
    !> Nunca para o programa; ierr = 1 e errmsg se falhar ao abrir/fechar.
    subroutine write_regions_file(filename, rs, ierr, errmsg)
        character(len=*),   intent(in)  :: filename
        type(region_set_t), intent(in)  :: rs
        integer,            intent(out) :: ierr
        character(len=*),   intent(out) :: errmsg
        integer                       :: iu, ios, ir, i, n
        character(len=512)            :: msg
        character(len=:), allocatable :: line

        ierr = 0
        errmsg = ''
        open(newunit=iu, file=filename, status='replace', action='write', iostat=ios, iomsg=msg)
        if (ios /= 0) then
            ierr = 1
            errmsg = 'cannot open regions file "' // trim(filename) // '" for writing: ' // trim(msg)
            return
        end if

        if (ios == 0) write(iu, '(A)', iostat=ios) '&REGION_FILE'
        if (ios == 0) write(iu, '(A)', iostat=ios) '    format_version = ' // itoa(rs%format_version)
        if (ios == 0) write(iu, '(A)', iostat=ios) '    natoms         = ' // itoa(rs%natoms)
        if (ios == 0) write(iu, '(A)', iostat=ios) "    fingerprint    = '" // trim(rs%fingerprint) // "'"
        if (rs%has_species .and. allocated(rs%species)) then
            n = size(rs%species)
            line = '    species        = '
            do i = 1, n
                line = line // "'" // trim(rs%species(i)) // "'"
                if (i < n) line = line // ', '
                if (mod(i, SPECIES_PER_LINE) == 0 .or. i == n) then
                    if (ios == 0) write(iu, '(A)', iostat=ios) trim(line)
                    line = '                     '
                end if
            end do
        end if
        if (ios == 0) write(iu, '(A)', iostat=ios) '/'

        do ir = 1, rs%nregion
            if (ios == 0) write(iu, '(A)', iostat=ios) ''
            if (ios == 0) write(iu, '(A)', iostat=ios) '&REGION'
            if (ios == 0) write(iu, '(A)', iostat=ios) "    name  = '" // trim(rs%regions(ir)%name) // "'"
            n = size(rs%regions(ir)%atoms)
            line = '    atoms = '
            do i = 1, n
                line = line // itoa(rs%regions(ir)%atoms(i))
                if (i < n) line = line // ', '
                if (mod(i, ATOMS_PER_LINE) == 0 .or. i == n) then
                    if (ios == 0) write(iu, '(A)', iostat=ios) trim(line)
                    line = '            '
                end if
            end do
            if (ios == 0) write(iu, '(A)', iostat=ios) '/'
        end do

        if (ios /= 0) then
            ierr = 1
            errmsg = 'error writing regions file "' // trim(filename) // '"'
            close(iu, iostat=ios)
            return
        end if

        close(iu, iostat=ios, iomsg=msg)
        if (ios /= 0) then
            ierr = 1
            errmsg = 'error closing regions file "' // trim(filename) // '": ' // trim(msg)
        end if
    end subroutine write_regions_file

    ! ==================================================================
    ! VALIDAÇÃO E EDIÇÃO
    ! ==================================================================

    !> Checa um conjunto de regiões contra a estrutura à qual será aplicado, nesta
    !> ordem: natoms = st%nat; lista species (se presente) igual aos rótulos da
    !> estrutura átomo a átomo; fingerprint igual a structure_fingerprint(st).
    !> Fingerprint ausente é erro, a menos que allow_no_fingerprint, caso em que
    !> um aviso é impresso em stdout. Nunca para o programa; ierr = 1 e errmsg
    !> na primeira falha.
    subroutine validate_regions(rs, st, allow_no_fingerprint, ierr, errmsg)
        type(region_set_t), intent(in)  :: rs
        type(structure_t),  intent(in)  :: st
        logical,            intent(in)  :: allow_no_fingerprint
        integer,            intent(out) :: ierr
        character(len=*),   intent(out) :: errmsg
        integer :: i
        character(len=16) :: fp

        ierr = 0
        errmsg = ''
        if (rs%natoms /= st%nat) then
            ierr = 1
            errmsg = 'regions: natoms = ' // itoa(rs%natoms) // ' differs from structure nat = ' // itoa(st%nat)
            return
        end if

        if (rs%has_species) then
            if (.not. allocated(rs%species) .or. .not. allocated(st%label)) then
                ierr = 1
                errmsg = 'regions: species list present but labels unavailable'
                return
            end if
            if (size(rs%species) /= st%nat .or. size(st%label) /= st%nat) then
                ierr = 1
                errmsg = 'regions: species list size differs from structure nat'
                return
            end if
            do i = 1, st%nat
                if (trim(adjustl(rs%species(i))) /= trim(adjustl(st%label(i)))) then
                    ierr = 1
                    errmsg = 'regions: species of atom ' // itoa(i) // " is '" // trim(adjustl(rs%species(i))) // &
                             "' but structure label is '" // trim(adjustl(st%label(i))) // "'"
                    return
                end if
            end do
        end if

        if (len_trim(rs%fingerprint) > 0) then
            fp = structure_fingerprint(st)
            if (trim(rs%fingerprint) /= fp) then
                ierr = 1
                errmsg = "regions: fingerprint '" // trim(rs%fingerprint) // &
                         "' differs from structure fingerprint '" // fp // "'"
                return
            end if
        else
            if (.not. allow_no_fingerprint) then
                ierr = 1
                errmsg = 'regions: file has no fingerprint (structure identity cannot be verified)'
                return
            end if
            print '(A)', '[regions warning] regions file has no fingerprint; structure identity not verified'
        end if
    end subroutine validate_regions

    !> Acrescenta a região `name` com os `atoms` (base 1) a rs (ordem do usuário
    !> mantida). Exige rs%format_version = 1 e rs%natoms >= 1 já definidos.
    !> Erros (ierr = 1, rs inalterado): nome inválido ou duplicado, lista vazia,
    !> índice fora de 1..natoms, índice duplicado na região.
    !> Sobreposição com outras regiões é permitida.
    subroutine add_region(rs, name, atoms, ierr, errmsg)
        type(region_set_t), intent(inout) :: rs
        character(len=*),   intent(in)    :: name
        integer,            intent(in)    :: atoms(:)
        integer,            intent(out)   :: ierr
        character(len=*),   intent(out)   :: errmsg
        type(atom_region_t), allocatable :: tmp(:)
        logical,             allocatable :: seen(:)
        integer                          :: i, ir

        ierr = 0
        errmsg = ''
        if (rs%format_version /= 1) then
            call fail('format_version must be 1 (got ' // itoa(rs%format_version) // ')')
            return
        end if
        if (rs%natoms < 1) then
            call fail('natoms must be >= 1 before adding regions (got ' // itoa(rs%natoms) // ')')
            return
        end if
        if (.not. is_valid_region_name(name)) then
            call fail("invalid region name '" // trim(name) // &
                      "' (1-64 chars from [A-Za-z0-9_.], 'total' is reserved)")
            return
        end if
        do ir = 1, rs%nregion
            if (trim(rs%regions(ir)%name) == trim(name)) then
                call fail("duplicated region name '" // trim(name) // "'")
                return
            end if
        end do
        if (size(atoms) == 0) then
            call fail("region '" // trim(name) // "' is empty")
            return
        end if
        allocate(seen(rs%natoms), source=.false.)
        do i = 1, size(atoms)
            if (atoms(i) < 1 .or. atoms(i) > rs%natoms) then
                call fail("region '" // trim(name) // "': atom index " // itoa(atoms(i)) // &
                          ' out of range 1..' // itoa(rs%natoms))
                return
            end if
            if (seen(atoms(i))) then
                call fail("region '" // trim(name) // "': duplicated atom index " // itoa(atoms(i)))
                return
            end if
            seen(atoms(i)) = .true.
        end do

        ! Aumenta o array em um; move_alloc evita cópia profunda das listas de átomos
        if (.not. allocated(rs%regions)) allocate(rs%regions(0))
        allocate(tmp(rs%nregion + 1))
        do ir = 1, rs%nregion
            call move_alloc(rs%regions(ir)%atoms, tmp(ir)%atoms)
            tmp(ir)%name = rs%regions(ir)%name
        end do
        tmp(rs%nregion + 1)%name = trim(name)
        tmp(rs%nregion + 1)%atoms = atoms
        call move_alloc(tmp, rs%regions)
        rs%nregion = rs%nregion + 1

    contains

        subroutine fail(m)
            character(len=*), intent(in) :: m
            ierr = 1
            errmsg = m
        end subroutine fail

    end subroutine add_region

    !> Verdadeiro se todo átomo 1..natoms pertence a exatamente uma região (sem
    !> sobreposição e sem átomo de fora). Falso para conjunto vazio ou índices
    !> fora da faixa.
    logical function regions_are_partition(rs)
        type(region_set_t), intent(in) :: rs
        integer, allocatable :: cnt(:)
        integer              :: ir, i, a

        regions_are_partition = .false.
        if (rs%natoms < 1) return
        allocate(cnt(rs%natoms), source=0)
        do ir = 1, rs%nregion
            do i = 1, size(rs%regions(ir)%atoms)
                a = rs%regions(ir)%atoms(i)
                if (a < 1 .or. a > rs%natoms) return
                cnt(a) = cnt(a) + 1
            end do
        end do
        regions_are_partition = all(cnt == 1)
    end function regions_are_partition

    !> Verdadeiro se `name` (brancos finais ignorados) tem 1..REGION_NAME_LEN
    !> caracteres de [A-Za-z0-9_.] e não é a palavra reservada 'total'
    !> (diferencia maiúsculas).
    logical function is_valid_region_name(name)
        character(len=*), intent(in) :: name
        integer :: n
        n = len_trim(name)
        is_valid_region_name = .false.
        if (n < 1 .or. n > REGION_NAME_LEN) return
        if (verify(name(1:n), NAME_CHARS) /= 0) return
        if (name(1:n) == 'total') return
        is_valid_region_name = .true.
    end function is_valid_region_name

    ! ==================================================================
    ! AUXILIARES PRIVADAS
    ! ==================================================================

    !> Número de entradas iniciais não sentinela; -1 se um não sentinela vier após uma sentinela.
    integer function count_leading_atoms(buf) result(n)
        integer, intent(in) :: buf(:)
        integer :: i
        n = size(buf)
        do i = 1, size(buf)
            if (buf(i) == ATOM_SENTINEL) then
                n = i - 1
                exit
            end if
        end do
        if (n < size(buf)) then
            if (any(buf(n+1:) /= ATOM_SENTINEL)) n = -1
        end if
    end function count_leading_atoms

    !> Igual a count_leading_atoms, para o buffer de species (sentinela = branco).
    integer function count_leading_species(buf) result(n)
        character(len=*), intent(in) :: buf(:)
        integer :: i
        n = size(buf)
        do i = 1, size(buf)
            if (len_trim(buf(i)) == 0) then
                n = i - 1
                exit
            end if
        end do
        if (n < size(buf)) then
            if (any(len_trim(buf(n+1:)) /= 0)) n = -1
        end if
    end function count_leading_species

    !> Texto de `line` até o primeiro branco, tab, '!' ou '/' (p. ex. '&REGION').
    function first_token(line) result(tok)
        character(len=*), intent(in) :: line
        character(len=len(line)) :: tok
        integer :: k
        k = scan(line, ' ' // achar(9) // '!/')
        if (k == 0) then
            tok = line
        else
            tok = line(1:k-1)
        end if
    end function first_token

    !> Cópia de s em maiúsculas ASCII.
    pure function upper(s) result(u)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: u
        integer :: i, c
        u = s
        do i = 1, len(s)
            c = iachar(s(i:i))
            if (c >= iachar('a') .and. c <= iachar('z')) u(i:i) = achar(c - 32)
        end do
    end function upper

    !> Inteiro para string sem brancos (para mensagens de erro).
    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=16) :: buf
        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

end module atomic_regions
