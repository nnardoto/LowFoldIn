! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! mod_hs_io.f90 -- Formato .hsio v2 ("LFIOFMT")
!
! DESENHO GERAL
! -------------
! O arquivo em disco NAO guarda mais H(r)/S(r) como array denso
! (NOrbs,NOrbs,Neighbors). Em vez disso guarda apenas os blocos
! atomo-atomo (TNOa1 x TNOa2) que realmente existem para cada vizinho R
! -- "bloco-esparso", no espirito de BSR (Block Sparse Row) / do proprio
! formato nativo do OpenMX .scfout (natn/ncn/FNAN), do SIESTA (listh) e
! do DBCSR do CP2K. Isso e' o denominador comum de qualquer codigo DFT
! de base local, entao serve de base para conversores futuros alem do
! OpenMX sem reprojetar o formato.
!
! O arquivo tem uma estrutura de CONTAINER generico e versionado:
!
!   [ Header fixo ] [ payload dataset 1 ] [ payload dataset 2 ] ... [ TOC ]
!
! O Header aponta para o offset do TOC (tabela de conteudo), que lista,
! para cada dataset gravado, nome/tipo/versao-de-schema/offset/tamanho.
! Hoje so' gravamos um dataset por arquivo (kind=DATASET_BLOCKSPARSE_HS,
! nome "H" ou "S"), mas o container ja' suporta multiplos datasets no
! mesmo arquivo -- e' a base para, no futuro, um .hsio guardar H+S+
! geometria+resultados de unfolding juntos, sem quebrar leitores antigos
! (um leitor que nao reconhece um "kind" simplesmente pula via
! offset+nbytes).
!
! COMPATIBILIDADE DE API EM MEMORIA
! ----------------------------------
! O tipo publico continua se chamando HSDataFormat, com os MESMOS campos
! que o formato v1 denso expunha (Title, Neighbors, nAtoms, NOrbs, Orbs,
! IndexList, Data, Degenerecency) -- ReadFile/WriteFile tem a mesma
! assinatura de antes. Isso significa que mod_core.f90, mod_unfold.f90,
! tools/unfold.f90 e tools/openmx2hsio.f90 NAO precisam mudar uma linha:
! WriteFile converte o Data denso passado pelo chamador para bloco-
! esparso internamente (descartando blocos atomo-atomo inteiramente
! nulos) antes de gravar; ReadFile le' o bloco-esparso do disco e
! reconstroi Data denso em memoria (zero-preenchido) antes de devolver.
!
! O ganho desta etapa e' de ESPACO EM DISCO (arquivos MUITO menores para
! sistemas com raio de corte curto e muitos atomos) e de base para um
! formato unificado. O ganho de CPU/RAM em tempo de execucao (evitar
! materializar o array denso e somar so' blocos nao-nulos no somatorio
! de Fourier de compute_Hk_at_kpoint) fica para uma proxima etapa, feita
! separadamente dos modulos numericos (mod_core/mod_unfold) para nao
! arriscar a fisica ja validada sem testes dedicados.
!
! Nathanael N. Batista
! Agosto de 2026
! ----------------------------------------------------------------------------------------------------------

module HS_io
    use, intrinsic :: iso_fortran_env, only: int64
    use types, only: sp, dp, qp
    implicit none
    private

    public :: HSDataFormat

    ! ------------------------------------------------------------------
    ! Constantes do container
    ! ------------------------------------------------------------------
    character(len=8), parameter :: HSIO_MAGIC = "LFIOFMT" // achar(0)
    integer, parameter :: HSIO_VER_MAJOR = 2   ! bump: v1 era o formato denso antigo
    integer, parameter :: HSIO_VER_MINOR = 0
    integer, parameter :: HSIO_ENDIAN_CANARY = 168496141   ! = 0x0A0B0C0D, deteccao de endianness
    integer, parameter :: HSIO_HEADER_BYTES  = 8 + 4 + 4 + 4 + 4 + 8 + 4  ! ver header abaixo

    integer, parameter :: TOC_NAME_LEN = 64
    integer, parameter :: HSIO_TOC_ENTRY_BYTES = TOC_NAME_LEN + 4 + 4 + 8 + 8

    integer, parameter, public :: DATASET_BLOCKSPARSE_HS = 1
    integer, parameter :: SCHEMA_BLOCKSPARSE_HS_V1 = 1

    ! ------------------------------------------------------------------
    ! Entrada da tabela de conteudo (TOC). Tamanho fixo em disco =
    ! HSIO_TOC_ENTRY_BYTES, para permitir acesso aleatorio (seek) sem
    ! precisar ler o arquivo inteiro so' para descobrir onde estao os
    ! datasets.
    ! ------------------------------------------------------------------
    type :: TOCEntry
        character(len=TOC_NAME_LEN) :: name       = ""
        integer                     :: kind       = 0
        integer                     :: schema_ver = 0
        integer(int64)              :: offset     = 0
        integer(int64)              :: nbytes     = 0
    end type TOCEntry

    ! ------------------------------------------------------------------
    ! Tipo publico. Campos "de fora" (Title..Data) sao a MESMA interface
    ! do formato v1 -- ver nota de compatibilidade no cabecalho do
    ! modulo. Os campos de baixo de "-- representacao esparsa interna --"
    ! sao populados/consumidos so' internamente por ReadFile/WriteFile.
    ! ------------------------------------------------------------------
    type :: HSDataFormat
        character(len=256)         :: Title                 = ""
        integer                    :: Neighbors             = 0
        integer                    :: nAtoms                = 0
        integer                    :: NOrbs                 = 0
        integer,     allocatable   :: Degenerecency(:)
        integer,     allocatable   :: Orbs(:)
        integer,     allocatable   :: IndexList(:, :)
        complex(dp), allocatable   :: Data(:, :, :)

        contains
            procedure :: ReadFile
            procedure :: WriteFile
            procedure :: DisplayInfo
            procedure :: Clear
    end type HSDataFormat

contains

! ==================================================================
! API publica (compativel com o formato v1)
! ==================================================================

subroutine ReadFile(this, filename)
    class(HSDataFormat), intent(inout)        :: this
    character(len=*),  intent(in)             :: Filename
    integer                                   :: File, ierr
    character(len=8)                          :: magic_read
    integer                                   :: header_bytes, ver_major, ver_minor, canary
    integer(int64)                            :: toc_offset
    integer                                   :: toc_count
    type(TOCEntry), allocatable               :: toc(:)
    integer                                   :: i, ds_idx

    call this%Clear()

    open(newunit=File, file=Filename, status='old', action='read', iostat=ierr, &
         form='unformatted', access='stream')
    if (ierr /= 0) then
        print '(A,A)', "[HS_io] Erro: nao foi possivel abrir para leitura: ", trim(Filename)
        stop 1
    end if

    ! -- Header ------------------------------------------------------
    read(File, pos=1, iostat=ierr) magic_read
    if (ierr /= 0 .or. magic_read /= HSIO_MAGIC) then
        print '(A,A)', "[HS_io] Erro: '", trim(Filename)
        print '(A)',   "' nao comeca com a assinatura .hsio v2 ('LFIOFMT'). " // &
                        "Pode ser um .hsio do formato antigo (denso, v1) -- " // &
                        "reconverta com openmx2hsio para o formato atual."
        stop 1
    end if
    read(File, iostat=ierr) header_bytes, ver_major, ver_minor, canary, toc_offset, toc_count
    if (ierr /= 0) then
        print '(A,A)', "[HS_io] Erro: cabecalho truncado ou ilegivel em ", trim(Filename)
        stop 1
    end if

    if (canary /= HSIO_ENDIAN_CANARY) then
        print '(A)', "[HS_io] Erro: marca de endianness inesperada -- arquivo pode ter sido " // &
                      "gravado em uma arquitetura com ordem de bytes diferente (nao suportado)."
        stop 1
    end if
    if (ver_major /= HSIO_VER_MAJOR) then
        print '(A,I0,A,I0)', "[HS_io] Erro: versao de container .hsio incompativel: arquivo = ", &
                              ver_major, ", leitor suporta major = ", HSIO_VER_MAJOR
        stop 1
    end if

    ! -- TOC -----------------------------------------------------------
    ! (leitura estruturada campo-a-campo: character/int32/int32/int64/
    ! int64 nao formam um bloco homogeneo, entao le'-se entrada por
    ! entrada em vez de um unico READ do array de derived type)
    if (toc_count < 1) then
        print '(A,A)', "[HS_io] Erro: TOC vazio ou invalido em ", trim(Filename)
        stop 1
    end if
    allocate(toc(toc_count))
    call read_toc_entries(File, toc_offset, toc_count, toc)

    ds_idx = -1
    do i = 1, toc_count
        if (toc(i)%kind == DATASET_BLOCKSPARSE_HS) then
            ds_idx = i
            exit
        end if
    end do
    if (ds_idx < 0) then
        print '(A,A)', "[HS_io] Erro: nenhum dataset H/S bloco-esparso encontrado em ", trim(Filename)
        stop 1
    end if

    call read_blocksparse_hs_payload(File, toc(ds_idx)%offset, toc(ds_idx)%schema_ver, this)

    close(File)
end subroutine ReadFile

subroutine WriteFile(this, filename)
    class(HSDataFormat), intent(in)     :: this
    character(len=*), intent(in)        :: Filename
    integer                             :: File, ierr
    integer(int64)                      :: payload_offset, toc_offset, cur_pos
    type(TOCEntry)                      :: entry

    if (.not. allocated(this%Data) .or. .not. allocated(this%Orbs) &
        .or. .not. allocated(this%IndexList)) then
        print '(A)', "[HS_io] Erro: WriteFile chamado com Data/Orbs/IndexList nao alocados."
        stop 1
    end if

    open(newunit=File, file=Filename, status='replace', action='write', iostat=ierr, &
         form='unformatted', access='stream')
    if (ierr /= 0) then
        print '(A,A)', "[HS_io] Erro: nao foi possivel abrir para escrita: ", trim(Filename)
        stop 1
    end if

    ! -- Header (placeholder para toc_offset/toc_count, corrigido no fim) --
    write(File, pos=1) HSIO_MAGIC
    write(File) HSIO_HEADER_BYTES
    write(File) HSIO_VER_MAJOR
    write(File) HSIO_VER_MINOR
    write(File) HSIO_ENDIAN_CANARY
    write(File) 0_int64          ! toc_offset (placeholder)
    write(File) 0                ! toc_count  (placeholder)

    ! -- Payload: dataset unico, H ou S bloco-esparso -------------------
    inquire(File, pos=cur_pos)
    payload_offset = cur_pos
    call write_blocksparse_hs_payload(File, this)

    ! -- TOC --------------------------------------------------------
    inquire(File, pos=cur_pos)
    toc_offset = cur_pos

    entry%name       = merge("S", "H", is_overlap_like(this))
    entry%kind       = DATASET_BLOCKSPARSE_HS
    entry%schema_ver = SCHEMA_BLOCKSPARSE_HS_V1
    entry%offset     = payload_offset
    inquire(File, pos=cur_pos)
    entry%nbytes     = cur_pos - payload_offset   ! (nbytes ainda nao usado por este leitor, mas
                                                   !  fica disponivel para leitores futuros pularem
                                                   !  datasets desconhecidos sem parsear o payload)
    call write_toc_entry(File, entry)

    ! -- Corrige o header com o offset/contagem reais do TOC ------------
    write(File, pos=1 + 8 + 4 + 4 + 4 + 4) toc_offset
    write(File) 1   ! toc_count = 1 (um unico dataset por arquivo, por enquanto)

    close(File)
end subroutine WriteFile

subroutine Clear(this)
    class(HSDataFormat), intent(inout)   :: this
    if (allocated(this%Orbs))           deallocate(this%Orbs)
    if (allocated(this%IndexList))      deallocate(this%IndexList)
    if (allocated(this%Data))           deallocate(this%Data)
    if (allocated(this%Degenerecency))  deallocate(this%Degenerecency)

    this%Title     = ""
    this%Neighbors = 0
    this%nAtoms    = 0
    this%NOrbs     = 0
end subroutine Clear

! Sumário do conteúdo do arquivo (para depuração e verificação)
subroutine DisplayInfo(this)
    class(HSDataFormat), intent(in) :: this
    integer :: i, j
    integer(int64) :: nElemTotal, nElemNZ
    integer :: a1, a2, off1, off2, n1, n2
    integer, allocatable :: AtomOffset(:)
    logical :: any_nz
    integer(int64) :: nBlocksTotal, nBlocksNZ

    print '(A)', "================================================================="
    print '(A)', "                      INFORMAÇÕES DO ARQUIVO                     "
    print '(A)', "================================================================="

    print '(A, A)',       "  Título:         ", trim(adjustl(this%Title))
    print '(A, I8)',      "  Nº de Vizinhos: ", this%Neighbors
    print '(A, I8)',      "  Átomos / Célula:", this%nAtoms
    print '(A, I8)',      "  Total Orbitais: ", this%NOrbs
    print '(A)', "-----------------------------------------------------------------"

    if (allocated(this%Orbs)) then
        print '(A)', "  Orbitais por Átomo:"
        write(*, '("    ", 10(I0, 1X))') this%Orbs
    end if

    print '(A)', "-----------------------------------------------------------------"

    if (allocated(this%IndexList)) then
        print '(A)', "  Lista de Índices (h, k, l):"
        do i = 1, min(this%Neighbors, 5)
            write(*, '("    Vizinho ", I3, ": [", I4, ",", I4, ",", I4, " ]")') &
                i, this%IndexList(1, i), this%IndexList(2, i), this%IndexList(3, i)
        end do
        if (this%Neighbors > 5) then
            print '(A, I0, A)', "    ... e mais ", this%Neighbors - 5, " vizinhos."
        end if
    end if

    print '(A)', "-----------------------------------------------------------------"

    ! -- Estatística de esparsidade (elemento e bloco atomo-atomo) ------
    if (allocated(this%Data) .and. allocated(this%Orbs)) then
        allocate(AtomOffset(this%nAtoms + 1))
        AtomOffset(1) = 0
        do i = 1, this%nAtoms
            AtomOffset(i + 1) = AtomOffset(i) + this%Orbs(i)
        end do

        nElemTotal = int(this%NOrbs, int64) * int(this%NOrbs, int64) * int(this%Neighbors, int64)
        nElemNZ = count(this%Data /= cmplx(0.0_dp, 0.0_dp, kind=dp))

        nBlocksTotal = int(this%nAtoms, int64) * int(this%nAtoms, int64) * int(this%Neighbors, int64)
        nBlocksNZ = 0
        do i = 1, this%Neighbors
            do a1 = 1, this%nAtoms
                off1 = AtomOffset(a1); n1 = this%Orbs(a1)
                do a2 = 1, this%nAtoms
                    off2 = AtomOffset(a2); n2 = this%Orbs(a2)
                    any_nz = any(this%Data(off1+1:off1+n1, off2+1:off2+n2, i) /= cmplx(0.0_dp, 0.0_dp, kind=dp))
                    if (any_nz) nBlocksNZ = nBlocksNZ + 1
                end do
            end do
        end do
        deallocate(AtomOffset)

        print '(A)', "  Esparsidade:"
        write(*, '("    Elementos:          ", I0, " / ", I0, " nao-nulos (", F5.2, "%)")') &
            nElemNZ, nElemTotal, 100.0_dp * real(nElemNZ, dp) / real(max(nElemTotal, 1_int64), dp)
        write(*, '("    Blocos atomo-atomo: ", I0, " / ", I0, " nao-nulos (", F5.2, "%)")') &
            nBlocksNZ, nBlocksTotal, 100.0_dp * real(nBlocksNZ, dp) / real(max(nBlocksTotal, 1_int64), dp)
    end if

    print '(A)', "-----------------------------------------------------------------"

    if (allocated(this%Data)) then
        print '(A)', "  Dados Complexos Data(1, 1:min(NOrbs,3), 1:min(NOrbs,3)):"
        print '(A)', "  (Amostra do 1º Vizinho):"

        do i = 1, min(this%NOrbs, 3)
            write(*, '("    ")', advance='no')
            do j = 1, min(this%NOrbs, 3)
                if (aimag(this%Data(i, j, 1)) >= 0.0_8) then
                    write(*, '("(", F7.3, "+", F7.3, "i) ")', advance='no') &
                        real(this%Data(i, j, 1)), aimag(this%Data(i, j, 1))
                else
                    write(*, '("(", F7.3, "-", F7.3, "i) ")', advance='no') &
                        real(this%Data(i, j, 1)), abs(aimag(this%Data(i, j, 1)))
                end if
            end do
            print *
        end do

        if (this%NOrbs > 3) then
            print '(A)', "    ..."
        end if
    end if

    print '(A)', "================================================================="
end subroutine DisplayInfo

! ==================================================================
! Internals do container: TOC
! ==================================================================

subroutine write_toc_entry(File, entry)
    integer, intent(in)          :: File
    type(TOCEntry), intent(in)   :: entry

    write(File) entry%name
    write(File) entry%kind
    write(File) entry%schema_ver
    write(File) entry%offset
    write(File) entry%nbytes
end subroutine write_toc_entry

subroutine read_toc_entries(File, toc_offset, toc_count, toc)
    integer, intent(in)                       :: File
    integer(int64), intent(in)                :: toc_offset
    integer, intent(in)                       :: toc_count
    type(TOCEntry), intent(out)               :: toc(toc_count)
    integer :: i, ierr
    integer(int64) :: pos_i

    do i = 1, toc_count
        ! toc_offset ja' vem de INQUIRE(POS=) (1-based, convencao do
        ! proprio Fortran), entao NAO se soma +1 aqui -- ver nota em
        ! ReadFile/WriteFile sobre offsets em disco vs. POS 1-based.
        pos_i = toc_offset + int(i - 1, int64) * int(HSIO_TOC_ENTRY_BYTES, int64)
        read(File, pos=pos_i, iostat=ierr) toc(i)%name
        if (ierr == 0) read(File, iostat=ierr) toc(i)%kind
        if (ierr == 0) read(File, iostat=ierr) toc(i)%schema_ver
        if (ierr == 0) read(File, iostat=ierr) toc(i)%offset
        if (ierr == 0) read(File, iostat=ierr) toc(i)%nbytes
        if (ierr /= 0) then
            print '(A)', "[HS_io] Erro: TOC truncado ou ilegivel."
            stop 1
        end if
    end do
end subroutine read_toc_entries

! ==================================================================
! Internals do container: payload do dataset bloco-esparso H/S
!
! Layout em disco (schema_ver = SCHEMA_BLOCKSPARSE_HS_V1):
!
!   Title            char(256)
!   Neighbors        int32              ! = nCells
!   nAtoms           int32
!   NOrbs            int32
!   Degenerecency(Neighbors)  int32
!   Orbs(nAtoms)              int32
!   IndexList(3,Neighbors)    int32
!   nBlocks          int64
!   BlockIndex(3,nBlocks)     int32     ! (a1, a2, Rn) 1-based
!   BlockOffset(nBlocks+1)    int64     ! prefix-sum em elementos escalares
!   BlockData(sum TNOa1*TNOa2 sobre os blocos)  complex128
! ==================================================================

subroutine write_blocksparse_hs_payload(File, this)
    integer, intent(in)              :: File
    class(HSDataFormat), intent(in)  :: this

    integer, allocatable    :: AtomOffset(:)
    integer, allocatable    :: BlockIndex(:, :)
    integer(int64), allocatable :: BlockOffset(:)
    complex(dp), allocatable :: BlockData(:)
    integer(int64) :: nBlocks
    integer :: a1, a2, Rn, off1, off2, n1, n2, i, j
    integer(int64) :: bpos, k

    allocate(AtomOffset(this%nAtoms + 1))
    AtomOffset(1) = 0
    do a1 = 1, this%nAtoms
        AtomOffset(a1 + 1) = AtomOffset(a1) + this%Orbs(a1)
    end do

    ! Passo 1: conta blocos e escalares nao-nulos
    call count_nonzero_blocks(this, AtomOffset, nBlocks)

    allocate(BlockIndex(3, max(nBlocks, 1_int64)))
    allocate(BlockOffset(nBlocks + 1))
    BlockOffset(1) = 0

    ! Passo 2: preenche indices/offsets
    k = 0
    do Rn = 1, this%Neighbors
        do a1 = 1, this%nAtoms
            off1 = AtomOffset(a1); n1 = this%Orbs(a1)
            do a2 = 1, this%nAtoms
                off2 = AtomOffset(a2); n2 = this%Orbs(a2)
                if (any(this%Data(off1+1:off1+n1, off2+1:off2+n2, Rn) /= cmplx(0.0_dp, 0.0_dp, kind=dp))) then
                    k = k + 1
                    BlockIndex(:, k) = [a1, a2, Rn]
                    BlockOffset(k + 1) = BlockOffset(k) + int(n1, int64) * int(n2, int64)
                end if
            end do
        end do
    end do

    allocate(BlockData(BlockOffset(nBlocks + 1)))
    do k = 1, nBlocks
        a1 = BlockIndex(1, k); a2 = BlockIndex(2, k); Rn = BlockIndex(3, k)
        off1 = AtomOffset(a1); n1 = this%Orbs(a1)
        off2 = AtomOffset(a2); n2 = this%Orbs(a2)
        bpos = BlockOffset(k)
        do j = 1, n2
            do i = 1, n1
                ! ordem column-major dentro do bloco (i varia mais rapido) para casar com o
                ! layout usado no scfout/openmx2hsio (ver mod_openmx_convert.f90)
                BlockData(bpos + int(j - 1, int64) * int(n1, int64) + i) = this%Data(off1 + i, off2 + j, Rn)
            end do
        end do
    end do

    write(File) this%Title
    write(File) this%Neighbors
    write(File) this%nAtoms
    write(File) this%NOrbs
    if (allocated(this%Degenerecency)) then
        write(File) this%Degenerecency
    else
        write(File) [(1, i = 1, this%Neighbors)]
    end if
    write(File) this%Orbs
    write(File) this%IndexList
    write(File) nBlocks
    if (nBlocks > 0) write(File) BlockIndex(:, 1:nBlocks)
    write(File) BlockOffset
    write(File) BlockData
end subroutine write_blocksparse_hs_payload

subroutine count_nonzero_blocks(this, AtomOffset, nBlocks)
    class(HSDataFormat), intent(in) :: this
    integer, intent(in)             :: AtomOffset(:)
    integer(int64), intent(out)     :: nBlocks
    integer :: a1, a2, Rn, off1, off2, n1, n2

    nBlocks = 0
    do Rn = 1, this%Neighbors
        do a1 = 1, this%nAtoms
            off1 = AtomOffset(a1); n1 = this%Orbs(a1)
            do a2 = 1, this%nAtoms
                off2 = AtomOffset(a2); n2 = this%Orbs(a2)
                if (any(this%Data(off1+1:off1+n1, off2+1:off2+n2, Rn) /= cmplx(0.0_dp, 0.0_dp, kind=dp))) then
                    nBlocks = nBlocks + 1
                end if
            end do
        end do
    end do
end subroutine count_nonzero_blocks

subroutine read_blocksparse_hs_payload(File, offset, schema_ver, this)
    integer, intent(in)                :: File
    integer(int64), intent(in)         :: offset
    integer, intent(in)                :: schema_ver
    class(HSDataFormat), intent(inout) :: this

    integer, allocatable         :: AtomOffset(:)
    integer, allocatable         :: BlockIndex(:, :)
    integer(int64), allocatable  :: BlockOffset(:)
    complex(dp), allocatable     :: BlockData(:)
    integer(int64)               :: nBlocks, bpos, k
    integer :: a1, a2, Rn, off1, off2, n1, n2, i, j, ierr

    if (schema_ver /= SCHEMA_BLOCKSPARSE_HS_V1) then
        print '(A,I0)', "[HS_io] Erro: schema do dataset H/S bloco-esparso desconhecido: ", schema_ver
        stop 1
    end if

    ! offset ja' vem de INQUIRE(POS=) (1-based) -- mesma convencao acima
    read(File, pos=offset, iostat=ierr) this%Title, this%Neighbors, this%nAtoms, this%NOrbs
    if (ierr /= 0) call payload_error("cabecalho do dataset truncado")
    if (this%Neighbors < 1 .or. this%nAtoms < 1 .or. this%NOrbs < 1) &
        call payload_error("Neighbors, nAtoms e NOrbs devem ser >= 1")

    allocate(this%Degenerecency(this%Neighbors))
    allocate(this%Orbs(this%nAtoms))
    allocate(this%IndexList(3, this%Neighbors))

    read(File, iostat=ierr) this%Degenerecency
    if (ierr == 0) read(File, iostat=ierr) this%Orbs
    if (ierr == 0) read(File, iostat=ierr) this%IndexList
    if (ierr == 0) read(File, iostat=ierr) nBlocks
    if (ierr /= 0) call payload_error("dataset truncado (Degenerecency/Orbs/IndexList/nBlocks)")
    if (any(this%Orbs < 1)) call payload_error("Orbs deve ser >= 1 para todo atomo")
    if (sum(this%Orbs) /= this%NOrbs) call payload_error("sum(Orbs) /= NOrbs")
    if (nBlocks < 0) call payload_error("nBlocks negativo")

    allocate(BlockIndex(3, max(nBlocks, 1_int64)))
    if (nBlocks > 0) then
        read(File, iostat=ierr) BlockIndex(:, 1:nBlocks)
        if (ierr /= 0) call payload_error("BlockIndex truncado")
    end if
    allocate(BlockOffset(nBlocks + 1))
    read(File, iostat=ierr) BlockOffset
    if (ierr /= 0) call payload_error("BlockOffset truncado")

    ! Valida cada bloco antes de alocar/ler o payload e de escrever em Data
    if (BlockOffset(1) /= 0_int64) call payload_error("BlockOffset(1) /= 0")
    do k = 1, nBlocks
        a1 = BlockIndex(1, k); a2 = BlockIndex(2, k); Rn = BlockIndex(3, k)
        if (a1 < 1 .or. a1 > this%nAtoms .or. a2 < 1 .or. a2 > this%nAtoms &
            .or. Rn < 1 .or. Rn > this%Neighbors) &
            call payload_error("indice de bloco fora do intervalo")
        if (BlockOffset(k + 1) - BlockOffset(k) /= int(this%Orbs(a1), int64) * int(this%Orbs(a2), int64)) &
            call payload_error("tamanho de bloco inconsistente com Orbs")
    end do

    allocate(BlockData(BlockOffset(nBlocks + 1)))
    read(File, iostat=ierr) BlockData
    if (ierr /= 0) call payload_error("BlockData truncado")

    allocate(AtomOffset(this%nAtoms + 1))
    AtomOffset(1) = 0
    do a1 = 1, this%nAtoms
        AtomOffset(a1 + 1) = AtomOffset(a1) + this%Orbs(a1)
    end do

    allocate(this%Data(this%NOrbs, this%NOrbs, this%Neighbors))
    this%Data = cmplx(0.0_dp, 0.0_dp, kind=dp)

    do k = 1, nBlocks
        a1 = BlockIndex(1, k); a2 = BlockIndex(2, k); Rn = BlockIndex(3, k)
        off1 = AtomOffset(a1); n1 = this%Orbs(a1)
        off2 = AtomOffset(a2); n2 = this%Orbs(a2)
        bpos = BlockOffset(k)
        do j = 1, n2
            do i = 1, n1
                this%Data(off1 + i, off2 + j, Rn) = BlockData(bpos + int(j - 1, int64) * int(n1, int64) + i)
            end do
        end do
    end do

contains

    subroutine payload_error(msg)
        character(len=*), intent(in) :: msg
        print '(A,A)', "[HS_io] Erro: arquivo .hsio invalido: ", msg
        stop 1
    end subroutine payload_error

end subroutine read_blocksparse_hs_payload

! Heuristica so' para dar um nome amigavel ao dataset no TOC (nao afeta
! a leitura, que busca por "kind"): S costuma ter Data(i,i,cel.central)
! proximo de 1 na diagonal; H nao. Best-effort, cosmetico.
function is_overlap_like(this) result(is_S)
    class(HSDataFormat), intent(in) :: this
    logical :: is_S
    integer :: center
    is_S = .false.
    if (.not. allocated(this%Data) .or. .not. allocated(this%IndexList)) return
    if (this%NOrbs < 1) return
    do center = 1, this%Neighbors
        if (all(this%IndexList(:, center) == 0)) then
            is_S = abs(real(this%Data(1, 1, center), dp) - 1.0_dp) < 0.3_dp
            return
        end if
    end do
end function is_overlap_like

end module HS_io
