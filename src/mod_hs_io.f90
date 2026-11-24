! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! mod_hs_io.f90 -- Formato .hsio v2 ("LFIOFMT")
!
! DESENHO GERAL
! -------------
! O arquivo em disco NÃO guarda mais H(r)/S(r) como array denso
! (NOrbs,NOrbs,Neighbors). Em vez disso guarda apenas os blocos
! átomo-átomo (TNOa1 x TNOa2) que realmente existem para cada vizinho R
! -- "bloco-esparso", no espírito de BSR (Block Sparse Row) / do próprio
! formato nativo do OpenMX .scfout (natn/ncn/FNAN), do SIESTA (listh) e
! do DBCSR do CP2K. Isso é o denominador comum de qualquer código DFT
! de base local, então serve de base para conversores futuros além do
! OpenMX sem reprojetar o formato.
!
! O arquivo tem uma estrutura de CONTAINER generico e versionado:
!
!   [ Header fixo ] [ payload dataset 1 ] [ payload dataset 2 ] ... [ TOC ]
!
! O Header aponta para o offset do TOC (tabela de conteúdo), que lista,
! para cada dataset gravado, nome/tipo/versão-de-schema/offset/tamanho.
! Hoje só gravamos um dataset por arquivo (kind=DATASET_BLOCKSPARSE_HS,
! nome "H" ou "S"), mas o container já suporta múltiplos datasets no
! mesmo arquivo -- é a base para, no futuro, um .hsio guardar H+S+
! geometria+resultados de unfolding juntos, sem quebrar leitores antigos
! (um leitor que não reconhece um "kind" simplesmente pula via
! offset+nbytes).
!
! COMPATIBILIDADE DE API EM MEMÓRIA
! ----------------------------------
! O tipo público continua se chamando HSDataFormat, com os MESMOS campos
! que o formato v1 denso expunha (Title, Neighbors, nAtoms, NOrbs, Orbs,
! IndexList, Data, Degenerecency) -- ReadFile/WriteFile tem a mesma
! assinatura de antes. Isso significa que mod_core.f90, mod_unfold.f90,
! tools/unfold.f90 e tools/openmx2hsio.f90 NÃO precisam mudar uma linha:
! WriteFile converte o Data denso passado pelo chamador para bloco-
! esparso internamente (descartando blocos átomo-átomo inteiramente
! nulos) antes de gravar; ReadFile lê o bloco-esparso do disco e
! reconstroi Data denso em memória (zero-preenchido) antes de devolver.
!
! O ganho desta etapa é de ESPAÇO EM DISCO (arquivos MUITO menores para
! sistemas com raio de corte curto e muitos átomos) e de base para um
! formato unificado. O ganho de CPU/RAM em tempo de execução (evitar
! materializar o array denso e somar só blocos não-nulos no somatório
! de Fourier de compute_Hk_at_kpoint) fica para uma próxima etapa, feita
! separadamente dos módulos numéricos (mod_core/mod_unfold) para não
! arriscar a física já validada sem testes dedicados.
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
    public :: HSBlockSparse, read_hs_blocksparse, write_hs_blocksparse

    ! ------------------------------------------------------------------
    ! Constantes do container
    ! ------------------------------------------------------------------
    character(len=8), parameter :: HSIO_MAGIC = "LFIOFMT" // achar(0)
    integer,          parameter :: HSIO_VER_MAJOR = 2                              ! incremento: v1 era o formato denso antigo
    integer,          parameter :: HSIO_VER_MINOR = 0
    integer,          parameter :: HSIO_ENDIAN_CANARY = 168496141                  ! = 0x0A0B0C0D, detecção de endianness
    integer,          parameter :: HSIO_HEADER_BYTES  = 8 + 4 + 4 + 4 + 4 + 8 + 4  ! ver header abaixo

    integer, parameter :: TOC_NAME_LEN = 64
    integer, parameter :: HSIO_TOC_ENTRY_BYTES = TOC_NAME_LEN + 4 + 4 + 8 + 8

    integer, parameter, public :: DATASET_BLOCKSPARSE_HS = 1
    integer, parameter         :: SCHEMA_BLOCKSPARSE_HS_V1 = 1
    !> Dataset opcional de metadados de energia (E_F e referência). Leitores
    !> antigos procuram apenas DATASET_BLOCKSPARSE_HS e ignoram este.
    integer, parameter, public :: DATASET_META = 2
    integer, parameter         :: SCHEMA_META_V1 = 1
    integer, parameter, public :: HSIO_EREF_LEN = 32

    ! ------------------------------------------------------------------
    ! Entrada da tabela de conteúdo (TOC). Tamanho fixo em disco =
    ! HSIO_TOC_ENTRY_BYTES, para permitir acesso aleatório (seek) sem
    ! precisar ler o arquivo inteiro só para descobrir onde estão os
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
    ! Tipo público. Campos "de fora" (Title..Data) são a MESMA interface
    ! do formato v1 -- ver nota de compatibilidade no cabeçalho do
    ! módulo. Os campos de baixo de "-- representação esparsa interna --"
    ! são populados/consumidos só internamente por ReadFile/WriteFile.
    ! ------------------------------------------------------------------
    type :: HSDataFormat
        character(len=256)              :: Title                 = ""
        integer                         :: Neighbors             = 0
        integer                         :: nAtoms                = 0
        integer                         :: NOrbs                 = 0
        integer,            allocatable :: Degenerecency(:)
        integer,            allocatable :: Orbs(:)
        integer,            allocatable :: IndexList(:, :)
        complex(dp),        allocatable :: Data(:, :, :)
        ! -- metadados opcionais (dataset "meta"); has_meta = .false. em
        !    arquivos antigos ou gravados sem essa informação
        logical                         :: has_meta              = .false.
        real(dp)                        :: fermi_ev              = 0.0_dp   ! E_F do cálculo de origem, eV, na escala absoluta
        logical                         :: fermi_shift_applied   = .false.  ! .true. se H já contém -E_F*S
        character(len=HSIO_EREF_LEN)    :: energy_reference      = ""       ! 'E-E_F' ou 'absolute'

        contains
            procedure :: ReadFile
            procedure :: WriteFile
            procedure :: DisplayInfo
            procedure :: Clear
    end type HSDataFormat

    ! ------------------------------------------------------------------
    ! Representação bloco-esparsa em memória (sem Data denso): os mesmos
    ! campos do payload em disco. Usada por quem não pode materializar
    ! H(R) denso (supercélulas com milhares de orbitais e centenas de
    ! células) -- ver read_hs_blocksparse/write_hs_blocksparse. Bloco k:
    ! átomos (BlockIndex(1,k), BlockIndex(2,k)), célula IndexList(:,
    ! BlockIndex(3,k)); elementos em BlockData(BlockOffset(k)+1 :
    ! BlockOffset(k+1)), ordem de colunas (i mais rápido).
    ! ------------------------------------------------------------------
    type :: HSBlockSparse
        character(len=256)              :: Title                 = ""
        integer                         :: Neighbors             = 0
        integer                         :: nAtoms                = 0
        integer                         :: NOrbs                 = 0
        integer,            allocatable :: Degenerecency(:)
        integer,            allocatable :: Orbs(:)
        integer,            allocatable :: IndexList(:, :)
        integer(int64)                  :: nBlocks               = 0
        integer,            allocatable :: BlockIndex(:, :)
        integer(int64),     allocatable :: BlockOffset(:)
        complex(dp),        allocatable :: BlockData(:)
        logical                         :: has_meta              = .false.
        real(dp)                        :: fermi_ev              = 0.0_dp
        logical                         :: fermi_shift_applied   = .false.
        character(len=HSIO_EREF_LEN)    :: energy_reference      = ""
    end type HSBlockSparse

contains

    ! ==================================================================
    ! API pública (compatível com o formato v1)
    ! ==================================================================

    subroutine ReadFile(this, filename)
        class(HSDataFormat), intent(inout) :: this
        character(len=*),    intent(in)    :: Filename
        integer                            :: File, ierr
        character(len=8)                   :: magic_read
        integer                            :: header_bytes, ver_major, ver_minor, canary
        integer(int64)                     :: toc_offset
        integer                            :: toc_count
        type(TOCEntry),      allocatable   :: toc(:)
        integer                            :: i, ds_idx

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
        ! int64 não formam um bloco homogêneo, então lê-se entrada por
        ! entrada em vez de um único READ do array de derived type)
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

        ! Dataset opcional de metadados (schemas desconhecidos são ignorados)
        do i = 1, toc_count
            if (toc(i)%kind == DATASET_META .and. toc(i)%schema_ver == SCHEMA_META_V1) then
                call read_meta_payload(File, toc(i)%offset, this)
                exit
            end if
        end do

        close(File)
    end subroutine ReadFile

    subroutine WriteFile(this, filename)
        class(HSDataFormat), intent(in) :: this
        character(len=*),    intent(in) :: Filename
        integer                         :: File, ierr
        integer(int64)                  :: payload_offset, toc_offset, cur_pos
        integer(int64)                  :: meta_offset, meta_nbytes
        type(TOCEntry)                  :: entry

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
        write(File) 0_int64  ! toc_offset (placeholder)
        write(File) 0        ! toc_count  (placeholder)

        ! -- Payload: dataset único, H ou S bloco-esparso -------------------
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
        entry%nbytes     = cur_pos - payload_offset   ! (nbytes ainda não usado por este leitor, mas
                                                       !  fica disponível para leitores futuros pularem
                                                       !  datasets desconhecidos sem parsear o payload)
        ! -- Dataset opcional de metadados (antes do TOC) ----------------
        if (this%has_meta) then
            meta_offset = toc_offset
            call write_meta_payload(File, this)
            inquire(File, pos=cur_pos)
            meta_nbytes = cur_pos - meta_offset
            toc_offset = cur_pos
        end if

        call write_toc_entry(File, entry)
        if (this%has_meta) then
            entry%name       = "meta"
            entry%kind       = DATASET_META
            entry%schema_ver = SCHEMA_META_V1
            entry%offset     = meta_offset
            entry%nbytes     = meta_nbytes
            call write_toc_entry(File, entry)
        end if

        ! -- Corrige o header com o offset/contagem reais do TOC ------------
        write(File, pos=1 + 8 + 4 + 4 + 4 + 4) toc_offset
        write(File) merge(2, 1, this%has_meta)   ! toc_count

        close(File)
    end subroutine WriteFile

    subroutine Clear(this)
        class(HSDataFormat), intent(inout)   :: this
        if (allocated(this%Orbs))           deallocate(this%Orbs)
        if (allocated(this%IndexList))      deallocate(this%IndexList)
        if (allocated(this%Data))           deallocate(this%Data)
        if (allocated(this%Degenerecency))  deallocate(this%Degenerecency)

        this%Title     = ""
        this%has_meta  = .false.
        this%fermi_ev  = 0.0_dp
        this%fermi_shift_applied = .false.
        this%energy_reference = ""
        this%Neighbors = 0
        this%nAtoms    = 0
        this%NOrbs     = 0
    end subroutine Clear

    ! Sumário do conteúdo do arquivo (para depuração e verificação)
    subroutine DisplayInfo(this)
        class(HSDataFormat), intent(in)  :: this
        integer                          :: i, j
        integer(int64)                   :: nElemTotal, nElemNZ
        integer                          :: a1, a2, off1, off2, n1, n2
        integer,             allocatable :: AtomOffset(:)
        logical                          :: any_nz
        integer(int64)                   :: nBlocksTotal, nBlocksNZ

        print '(A)', "================================================================="
        print '(A)', "                      INFORMAÇÕES DO ARQUIVO                     "
        print '(A)', "================================================================="

        print '(A, A)',       "  Título:         ", trim(adjustl(this%Title))
        print '(A, I8)',      "  Nº de Vizinhos: ", this%Neighbors
        print '(A, I8)',      "  Átomos / Célula:", this%nAtoms
        if (this%has_meta) then
            print '(A, F14.6, A)', "  E_F (absoluto): ", this%fermi_ev, " eV"
            print '(A, A, A, L1, A)', "  Referência:     ", trim(this%energy_reference), &
                                      " (shift aplicado: ", this%fermi_shift_applied, ")"
        end if
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

        ! -- Estatística de esparsidade (elemento e bloco átomo-átomo) ------
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
    ! Internals do container: dataset "meta" (schema 1)
    !
    !   fermi_ev             real64   E_F do cálculo de origem (eV, absoluto)
    !   fermi_shift_applied  int32    1 se H já contém -E_F*S, 0 caso contrário
    !   energy_reference     char(32) 'E-E_F' ou 'absolute'
    ! ==================================================================

    subroutine write_meta_payload(File, this)
        integer,             intent(in) :: File
        class(HSDataFormat), intent(in) :: this
        write(File) this%fermi_ev
        write(File) merge(1, 0, this%fermi_shift_applied)
        write(File) this%energy_reference
    end subroutine write_meta_payload

    subroutine read_meta_payload(File, offset, this)
        integer,             intent(in)    :: File
        integer(int64),      intent(in)    :: offset
        class(HSDataFormat), intent(inout) :: this
        integer :: ierr, flag
        real(dp) :: ef
        character(len=HSIO_EREF_LEN) :: eref
        read(File, pos=offset, iostat=ierr) ef
        if (ierr == 0) read(File, iostat=ierr) flag
        if (ierr == 0) read(File, iostat=ierr) eref
        if (ierr /= 0) then
            print '(A)', "[HS_io] Aviso: dataset 'meta' ilegivel; ignorado."
            return
        end if
        this%has_meta = .true.
        this%fermi_ev = ef
        this%fermi_shift_applied = (flag /= 0)
        this%energy_reference = eref
    end subroutine read_meta_payload

    ! ==================================================================
    ! Internals do container: TOC
    ! ==================================================================

    subroutine write_toc_entry(File, entry)
        integer,        intent(in) :: File
        type(TOCEntry), intent(in) :: entry

        write(File) entry%name
        write(File) entry%kind
        write(File) entry%schema_ver
        write(File) entry%offset
        write(File) entry%nbytes
    end subroutine write_toc_entry

    subroutine read_toc_entries(File, toc_offset, toc_count, toc)
        integer,        intent(in)  :: File
        integer(int64), intent(in)  :: toc_offset
        integer,        intent(in)  :: toc_count
        type(TOCEntry), intent(out) :: toc(toc_count)
        integer                     :: i, ierr
        integer(int64)              :: pos_i

        do i = 1, toc_count
            ! toc_offset já vem de INQUIRE(POS=) (1-based, convenção do
            ! próprio Fortran), então NÃO se soma +1 aqui -- ver nota em
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
        integer,             intent(in) :: File
        class(HSDataFormat), intent(in) :: this

        integer,        allocatable :: AtomOffset(:)
        integer,        allocatable :: BlockIndex(:, :)
        integer(int64), allocatable :: BlockOffset(:)
        complex(dp),    allocatable :: BlockData(:)
        integer(int64)              :: nBlocks
        integer                     :: a1, a2, Rn, off1, off2, n1, n2, i, j
        integer(int64)              :: bpos, k

        allocate(AtomOffset(this%nAtoms + 1))
        AtomOffset(1) = 0
        do a1 = 1, this%nAtoms
            AtomOffset(a1 + 1) = AtomOffset(a1) + this%Orbs(a1)
        end do

        ! Passo 1: conta blocos e escalares não-nulos
        call count_nonzero_blocks(this, AtomOffset, nBlocks)

        allocate(BlockIndex(3, max(nBlocks, 1_int64)))
        allocate(BlockOffset(nBlocks + 1))
        BlockOffset(1) = 0

        ! Passo 2: preenche índices/offsets
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
                    ! ordem column-major dentro do bloco (i varia mais rápido) para casar com o
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
        class(HSDataFormat), intent(in)  :: this
        integer,             intent(in)  :: AtomOffset(:)
        integer(int64),      intent(out) :: nBlocks
        integer                          :: a1, a2, Rn, off1, off2, n1, n2

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
        integer,             intent(in)    :: File
        integer(int64),      intent(in)    :: offset
        integer,             intent(in)    :: schema_ver
        class(HSDataFormat), intent(inout) :: this

        integer,        allocatable :: AtomOffset(:)
        integer,        allocatable :: BlockIndex(:, :)
        integer(int64), allocatable :: BlockOffset(:)
        complex(dp),    allocatable :: BlockData(:)
        integer(int64)              :: nBlocks, bpos, k
        integer                     :: a1, a2, Rn, off1, off2, n1, n2, i, j, ierr

        if (schema_ver /= SCHEMA_BLOCKSPARSE_HS_V1) then
            print '(A,I0)', "[HS_io] Erro: schema do dataset H/S bloco-esparso desconhecido: ", schema_ver
            stop 1
        end if

        ! offset já vem de INQUIRE(POS=) (1-based) -- mesma convenção acima
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
                .or. Rn < 1 .or. Rn > this%Neighbors)                          &
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

    ! Heurística só para dar um nome amigavel ao dataset no TOC (não afeta
    ! a leitura, que busca por "kind"): S costuma ter Data(i,i,cel.central)
    ! próximo de 1 na diagonal; H não. Melhor esforço, cosmético.
    function is_overlap_like(this) result(is_S)
        class(HSDataFormat), intent(in) :: this
        logical                         :: is_S
        integer                         :: center
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

    ! ==================================================================
    ! API bloco-esparsa (sem Data denso; erros por status, sem stop)
    ! ==================================================================

    !> Lê um .hsio v2 na forma bloco-esparsa, sem montar Data denso.
    !> ierr: 1 abertura, 2 cabeçalho/TOC inválido, 3 payload inválido.
    subroutine read_hs_blocksparse(filename, bs, ierr, errmsg)
        character(len=*),    intent(in)  :: filename
        type(HSBlockSparse), intent(out) :: bs
        integer,             intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: u, ios, header_bytes, ver_major, ver_minor, canary, toc_count, i, ds_idx, flag
        integer(int64) :: toc_offset, pos_i, k
        character(len=8) :: magic_read
        type(TOCEntry), allocatable :: toc(:)
        real(dp) :: ef
        character(len=HSIO_EREF_LEN) :: eref

        ierr = 0
        errmsg = ''
        open(newunit=u, file=filename, status='old', action='read', iostat=ios, &
             form='unformatted', access='stream')
        if (ios /= 0) then
            ierr = 1; errmsg = "cannot open .hsio file: '" // trim(filename) // "'"; return
        end if
        read(u, pos=1, iostat=ios) magic_read
        if (ios /= 0 .or. magic_read /= HSIO_MAGIC) then
            call fail(2, 'not a .hsio v2 file (LFIOFMT signature missing)'); return
        end if
        read(u, iostat=ios) header_bytes, ver_major, ver_minor, canary, toc_offset, toc_count
        if (ios /= 0) then
            call fail(2, 'truncated header'); return
        end if
        if (canary /= HSIO_ENDIAN_CANARY) then
            call fail(2, 'unexpected endianness'); return
        end if
        if (ver_major /= HSIO_VER_MAJOR) then
            call fail(2, 'incompatible .hsio container version'); return
        end if
        if (toc_count < 1) then
            call fail(2, 'empty TOC'); return
        end if
        allocate(toc(toc_count))
        do i = 1, toc_count
            pos_i = toc_offset + int(i - 1, int64) * int(HSIO_TOC_ENTRY_BYTES, int64)
            read(u, pos=pos_i, iostat=ios) toc(i)%name
            if (ios == 0) read(u, iostat=ios) toc(i)%kind
            if (ios == 0) read(u, iostat=ios) toc(i)%schema_ver
            if (ios == 0) read(u, iostat=ios) toc(i)%offset
            if (ios == 0) read(u, iostat=ios) toc(i)%nbytes
            if (ios /= 0) then
                call fail(2, 'truncated TOC'); return
            end if
        end do
        ds_idx = -1
        do i = 1, toc_count
            if (toc(i)%kind == DATASET_BLOCKSPARSE_HS) then
                ds_idx = i
                exit
            end if
        end do
        if (ds_idx < 0) then
            call fail(2, 'no block-sparse H/S dataset'); return
        end if
        if (toc(ds_idx)%schema_ver /= SCHEMA_BLOCKSPARSE_HS_V1) then
            call fail(3, 'unknown block-sparse schema'); return
        end if

        read(u, pos=toc(ds_idx)%offset, iostat=ios) bs%Title, bs%Neighbors, bs%nAtoms, bs%NOrbs
        if (ios /= 0 .or. bs%Neighbors < 1 .or. bs%nAtoms < 1 .or. bs%NOrbs < 1) then
            call fail(3, 'invalid dataset header'); return
        end if
        allocate(bs%Degenerecency(bs%Neighbors), bs%Orbs(bs%nAtoms), bs%IndexList(3, bs%Neighbors))
        read(u, iostat=ios) bs%Degenerecency
        if (ios == 0) read(u, iostat=ios) bs%Orbs
        if (ios == 0) read(u, iostat=ios) bs%IndexList
        if (ios == 0) read(u, iostat=ios) bs%nBlocks
        if (ios /= 0) then
            call fail(3, 'truncated dataset (Degenerecency/Orbs/IndexList/nBlocks)'); return
        end if
        if (any(bs%Orbs < 1) .or. sum(bs%Orbs) /= bs%NOrbs .or. bs%nBlocks < 0) then
            call fail(3, 'inconsistent Orbs/NOrbs/nBlocks'); return
        end if
        allocate(bs%BlockIndex(3, bs%nBlocks), bs%BlockOffset(bs%nBlocks + 1))
        if (bs%nBlocks > 0) read(u, iostat=ios) bs%BlockIndex
        if (ios == 0) read(u, iostat=ios) bs%BlockOffset
        if (ios /= 0) then
            call fail(3, 'truncated BlockIndex/BlockOffset'); return
        end if
        if (bs%BlockOffset(1) /= 0_int64) then
            call fail(3, 'BlockOffset(1) /= 0'); return
        end if
        do k = 1, bs%nBlocks
            if (any(bs%BlockIndex(1:2, k) < 1) .or. any(bs%BlockIndex(1:2, k) > bs%nAtoms) .or. &
                bs%BlockIndex(3, k) < 1 .or. bs%BlockIndex(3, k) > bs%Neighbors) then
                call fail(3, 'block index out of range'); return
            end if
            if (bs%BlockOffset(k + 1) - bs%BlockOffset(k) /= &
                int(bs%Orbs(bs%BlockIndex(1, k)), int64) * int(bs%Orbs(bs%BlockIndex(2, k)), int64)) then
                call fail(3, 'block size inconsistent with Orbs'); return
            end if
        end do
        allocate(bs%BlockData(bs%BlockOffset(bs%nBlocks + 1)))
        read(u, iostat=ios) bs%BlockData
        if (ios /= 0) then
            call fail(3, 'truncated BlockData'); return
        end if
        do i = 1, toc_count
            if (toc(i)%kind == DATASET_META .and. toc(i)%schema_ver == SCHEMA_META_V1) then
                read(u, pos=toc(i)%offset, iostat=ios) ef
                if (ios == 0) read(u, iostat=ios) flag
                if (ios == 0) read(u, iostat=ios) eref
                if (ios == 0) then
                    bs%has_meta = .true.
                    bs%fermi_ev = ef
                    bs%fermi_shift_applied = (flag /= 0)
                    bs%energy_reference = eref
                end if
                exit
            end if
        end do
        close(u)

    contains

        subroutine fail(code, msg)
            integer,          intent(in) :: code
            character(len=*), intent(in) :: msg
            ierr = code
            errmsg = "invalid .hsio file '" // trim(filename) // "': " // msg
            close(u)
        end subroutine fail

    end subroutine read_hs_blocksparse

    !> Grava a forma bloco-esparsa no mesmo layout de WriteFile (dataset H/S
    !> + meta opcional). dsname é o nome do dataset no TOC ('H' ou 'S').
    subroutine write_hs_blocksparse(filename, bs, dsname, ierr, errmsg)
        character(len=*),    intent(in)  :: filename, dsname
        type(HSBlockSparse), intent(in)  :: bs
        integer,             intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: u, ios, i
        integer(int64) :: payload_offset, toc_offset, cur_pos, meta_offset, meta_nbytes
        type(TOCEntry) :: entry

        ierr = 0
        errmsg = ''
        if (.not. allocated(bs%Orbs) .or. .not. allocated(bs%IndexList) .or. &
            .not. allocated(bs%BlockIndex) .or. .not. allocated(bs%BlockOffset) .or. &
            .not. allocated(bs%BlockData)) then
            ierr = 1; errmsg = 'write_hs_blocksparse: incomplete block-sparse data'; return
        end if
        open(newunit=u, file=filename, status='replace', action='write', iostat=ios, &
             form='unformatted', access='stream')
        if (ios /= 0) then
            ierr = 1; errmsg = "cannot write .hsio file: '" // trim(filename) // "'"; return
        end if
        write(u, pos=1) HSIO_MAGIC
        write(u) HSIO_HEADER_BYTES
        write(u) HSIO_VER_MAJOR
        write(u) HSIO_VER_MINOR
        write(u) HSIO_ENDIAN_CANARY
        write(u) 0_int64
        write(u) 0
        inquire(u, pos=cur_pos)
        payload_offset = cur_pos
        write(u) bs%Title
        write(u) bs%Neighbors
        write(u) bs%nAtoms
        write(u) bs%NOrbs
        if (allocated(bs%Degenerecency)) then
            write(u) bs%Degenerecency
        else
            write(u) [(1, i = 1, bs%Neighbors)]
        end if
        write(u) bs%Orbs
        write(u) bs%IndexList
        write(u) bs%nBlocks
        if (bs%nBlocks > 0) write(u) bs%BlockIndex(:, 1:bs%nBlocks)
        write(u) bs%BlockOffset(1:bs%nBlocks + 1)
        write(u) bs%BlockData(1:bs%BlockOffset(bs%nBlocks + 1))
        inquire(u, pos=cur_pos)
        toc_offset = cur_pos
        entry%name       = dsname
        entry%kind       = DATASET_BLOCKSPARSE_HS
        entry%schema_ver = SCHEMA_BLOCKSPARSE_HS_V1
        entry%offset     = payload_offset
        entry%nbytes     = cur_pos - payload_offset
        if (bs%has_meta) then
            meta_offset = toc_offset
            write(u) bs%fermi_ev
            write(u) merge(1, 0, bs%fermi_shift_applied)
            write(u) bs%energy_reference
            inquire(u, pos=cur_pos)
            meta_nbytes = cur_pos - meta_offset
            toc_offset = cur_pos
        end if
        call write_toc_entry(u, entry)
        if (bs%has_meta) then
            entry%name       = "meta"
            entry%kind       = DATASET_META
            entry%schema_ver = SCHEMA_META_V1
            entry%offset     = meta_offset
            entry%nbytes     = meta_nbytes
            call write_toc_entry(u, entry)
        end if
        write(u, pos=1 + 8 + 4 + 4 + 4 + 4) toc_offset
        write(u, iostat=ios) merge(2, 1, bs%has_meta)
        close(u)
        if (ios /= 0) then
            ierr = 2; errmsg = "error writing .hsio file: '" // trim(filename) // "'"
        end if
    end subroutine write_hs_blocksparse

end module HS_io
