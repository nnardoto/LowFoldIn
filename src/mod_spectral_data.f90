! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Dados espectrais discretos (energias de banda + pesos de projeção por região
!> atômica ao longo de um caminho k) e seus formatos intermediários em disco
!> (docs/design/contracts.md §5).
!>
!> Posição no pipeline: o motor espectral produz um spectral_data_t e
!> write_spectral_file o armazena; as ferramentas posteriores (specmap) o leem de
!> volta com read_spectral_file, portanto não precisam do Hamiltoniano nem da
!> estrutura.
!>
!>   - arquivo binário em stream (write_spectral_file / read_spectral_file)
!>   - exportação em texto para depuração (write_spectral_text, nunca relida)
!>
!> Convenções
!>   - As unidades são fixadas pelo formato: energias em eV, distâncias k em 1/Angstrom.
!>   - Os formatos dos arrays seguem a ordem Fortran, com o índice de banda primeiro:
!>     energy(nband,nk,nspin), weight(nband,nk,nregion,nspin).
!>   - Todos os índices são de base 1. A região 1 sempre se chama 'total'.
!>   - nspin = 2 somente para spin colinear (spin_mode = SPIN_COLLINEAR).
!>   - Os arquivos binários usam a ordem de bytes nativa; um canário de endianness
!>     detecta incompatibilidade. No disco, inteiros são int32 e reais são real64,
!>     independentemente dos kinds padrão usados em memória.
!>
!> As rotinas da biblioteca nunca param: os erros são devolvidos via ierr/errmsg
!> (ierr = 0 significa sucesso).
module spectral_data
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use types, only: dp
    implicit none
    private

    ! ==================================================================
    ! Enumerações
    ! ==================================================================
    integer, parameter, public :: SPIN_NONE    = 0, SPIN_COLLINEAR = 1, SPIN_NONCOLLINEAR = 2
    integer, parameter, public :: METHOD_MULLIKEN = 1, METHOD_LOWDIN = 2                    ! esquema de projeção
    integer, parameter, public :: EREF_UNKNOWN = 0, EREF_FERMI = 1, EREF_ABSOLUTE = 2       ! zero de energia

    ! ==================================================================
    ! Tipo de dados
    ! ==================================================================

    !> Imagem em memória de um arquivo de dados espectrais. Os arrays alocáveis
    !> devem ser alocados com os formatos abaixo antes da escrita; a leitura
    !> os aloca. vertex_k/vertex_label são opcionais quando nvertex = 0.
    !> fermi_energy_ev só tem significado se fermi_known; energy_ref indica
    !> se as energias já estão deslocadas por E_F (EREF_FERMI) ou são absolutas.
    type, public :: spectral_data_t
        integer                        :: version_major = 1, version_minor = 0
        integer                        :: nk = 0, nband = 0, nspin = 1, nregion = 0, nvertex = 0
        integer                        :: spin_mode = SPIN_NONE, method = METHOD_MULLIKEN, energy_ref = EREF_UNKNOWN
        logical                        :: fermi_known = .false.
        real(dp)                       :: fermi_energy_ev = 0
        character(len=16)              :: fingerprint = ''
        character(len=512)             :: source = ''        ! texto livre: nomes dos arquivos de entrada
        real(dp), allocatable          :: kfrac(:,:)         ! (3,nk)
        real(dp), allocatable          :: kdist(:)           ! (nk) 1/Angstrom
        integer,  allocatable          :: vertex_k(:)        ! (nvertex)
        character(len=16), allocatable :: vertex_label(:)    ! (nvertex)
        character(len=64), allocatable :: region_name(:)     ! (nregion); região 1 = 'total'
        integer,  allocatable          :: region_natoms(:)   ! (nregion)
        real(dp), allocatable          :: energy(:,:,:)      ! (nband, nk, nspin) eV
        real(dp), allocatable          :: weight(:,:,:,:)    ! (nband, nk, nregion, nspin)
    end type spectral_data_t

    public :: write_spectral_file, read_spectral_file, write_spectral_text, find_region

    ! ==================================================================
    ! Constantes do formato binário (tabela de layout em contracts.md §5)
    ! ==================================================================
    character(len=8), parameter :: SD_MAGIC   = "LFSPECT" // achar(0)
    character(len=8), parameter :: SD_TRAILER = "LFSPEND" // achar(0)
    integer(int32),   parameter :: SD_VER_MAJOR = 1_int32, SD_VER_MINOR = 0_int32
    integer(int32),   parameter :: SD_CANARY    = 168496141_int32   ! 0x0A0B0C0D, sonda de ordem de bytes

    !> Tamanho fixo do cabeçalho em bytes (itens 1..9 da tabela de layout):
    !> magic + versões + canário + tamanhos de tipo + 5 dimensões + 4 flags + E_F +
    !> fingerprint + source.
    integer(int64), parameter :: SD_HEADER_BYTES = 8 + 2*4 + 4 + 3*4 + 5*4 + 4*4 + 8 + 16 + 512

contains

    ! ==================================================================
    ! API pública
    ! ==================================================================

    !> Índice da região cujo nome aparado é igual a trim(name); 0 se não existir.
    integer function find_region(sd, name) result(idx)
        type(spectral_data_t), intent(in) :: sd
        character(len=*),      intent(in) :: name
        integer                           :: i
        idx = 0
        if (.not. allocated(sd%region_name)) return
        do i = 1, size(sd%region_name)
            if (trim(adjustl(sd%region_name(i))) == trim(adjustl(name))) then
                idx = i
                return
            end if
        end do
    end function find_region

    !> Escreve sd em um arquivo binário em stream (layout de contracts.md §5).
    !> O arquivo é sobrescrito se existir. sd é validado antes (formatos,
    !> enumerações, região 1 = 'total'); as mesmas verificações valem para o
    !> escritor de texto.
    !> ierr: 0 ok; 1 sd inválido ou não abre o arquivo; 2 erro de E/S durante a
    !> escrita; 3 erro no fechamento. errmsg descreve a falha.
    subroutine write_spectral_file(filename, sd, ierr, errmsg)
        character(len=*),      intent(in)  :: filename
        type(spectral_data_t), intent(in)  :: sd
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        integer        :: u, ios, i
        integer(int32) :: fk

        call validate(sd, ierr, errmsg)
        if (ierr /= 0) return

        open(newunit=u, file=filename, access='stream', form='unformatted', &
             status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            errmsg = "write_spectral_file: cannot open '" // trim(filename) // "' for writing"
            return
        end if

        fk = 0_int32
        if (sd%fermi_known) fk = 1_int32

        ! Os registros seguem a ordem da tabela de layout. Qualquer erro de escrita
        ! deixa ios diferente de zero, o que interrompe as escritas restantes.
        ! Cabeçalho (itens 1..9)
        write(u, iostat=ios) SD_MAGIC, SD_VER_MAJOR, SD_VER_MINOR, SD_CANARY, &
            int(storage_size(0_int32)/8, int32), int(storage_size(0_int64)/8, int32), &
            int(storage_size(0.0_real64)/8, int32), &
            int(sd%nk, int32), int(sd%nband, int32), int(sd%nspin, int32), &
            int(sd%nregion, int32), int(sd%nvertex, int32), &
            int(sd%spin_mode, int32), int(sd%method, int32), int(sd%energy_ref, int32), fk, &
            real(sd%fermi_energy_ev, real64), sd%fingerprint, sd%source
        ! Pontos k (item 10), vértices (11-12), regiões (13-14)
        if (ios == 0) write(u, iostat=ios) real(sd%kfrac, real64), real(sd%kdist, real64)
        if (ios == 0 .and. sd%nvertex > 0) then
            write(u, iostat=ios) int(sd%vertex_k, int32)
            do i = 1, sd%nvertex
                if (ios == 0) write(u, iostat=ios) sd%vertex_label(i)
            end do
        end if
        do i = 1, sd%nregion
            if (ios == 0) write(u, iostat=ios) sd%region_name(i)
        end do
        if (ios == 0) write(u, iostat=ios) int(sd%region_natoms, int32)
        ! Arrays volumosos (15-16), depois n_ext = 0 (sem blocos de extensão na v1.0) e trailer
        if (ios == 0) write(u, iostat=ios) real(sd%energy, real64)
        if (ios == 0) write(u, iostat=ios) real(sd%weight, real64)
        if (ios == 0) write(u, iostat=ios) 0_int32, SD_TRAILER
        if (ios /= 0) then
            close(u, iostat=i)
            ierr = 2
            errmsg = "write_spectral_file: I/O error while writing '" // trim(filename) // "'"
            return
        end if
        close(u, iostat=ios)
        if (ios /= 0) then
            ierr = 3
            errmsg = "write_spectral_file: error closing '" // trim(filename) // "'"
            return
        end if
        ierr = 0
        errmsg = ''
    end subroutine write_spectral_file

    !> Lê um arquivo espectral binário, validando o layout completo.
    !> Todos os arrays são alocados por esta rotina. Antes de cada leitura, o número
    !> de bytes ainda disponíveis é conferido com o necessário, de modo que um
    !> arquivo truncado ou corrompido gera um erro em vez de uma falha.
    !> Em qualquer erro, sd é redefinido para o estado padrão (vazio).
    !> ierr: 0 ok; 1 arquivo ausente/não abre; 2 magic inválido; 3 versão,
    !> endianness ou tamanhos de tipo não suportados; 4 dimensões/flags inválidas;
    !> 5 conteúdo de vértice/região inválido; 6 bloco de extensão inválido;
    !> 7 trailer ausente ou dados extras; 8 arquivo truncado; 9 erro de E/S na leitura.
    !> Só version_major = 1 é aceita; qualquer versão minor é lida, e blocos de
    !> extensão desconhecidos são ignorados usando seu campo de tamanho.
    subroutine read_spectral_file(filename, sd, ierr, errmsg)
        character(len=*),      intent(in)  :: filename
        type(spectral_data_t), intent(out) :: sd
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        integer                     :: u, ios, i
        logical                     :: exists
        integer(int64)              :: fsize, pos, nbytes_ext
        character(len=8)            :: magic, trailer
        character(len=16)           :: tag
        integer(int32)              :: vmaj, vmin, canary, sz(3), dims(5), flags(4), n_ext
        integer(int32), allocatable :: itmp(:)
        real(real64)                :: ef

        ierr = 0
        errmsg = ''

        inquire(file=filename, exist=exists, size=fsize)
        if (.not. exists) then
            call fail(1, "file '" // trim(filename) // "' does not exist")
            return
        end if
        if (fsize < 0) then
            call fail(1, "cannot determine size of '" // trim(filename) // "'")
            return
        end if

        open(newunit=u, file=filename, access='stream', form='unformatted', &
             status='old', action='read', iostat=ios)
        if (ios /= 0) then
            call fail(1, "cannot open '" // trim(filename) // "' for reading")
            return
        end if

        ! `pos` é o deslocamento em bytes (base 1) do próximo byte não lido. É
        ! avançado manualmente após cada leitura para que have() possa conferir
        ! o tamanho restante; as leituras sem pos= continuam sequencialmente a
        ! partir da anterior.
        pos = 1
        ! ---- cabeçalho fixo ------------------------------------------
        if (.not. have(8_int64, "magic")) return
        read(u, pos=pos, iostat=ios) magic
        if (.not. ok_read("magic")) return
        if (magic /= SD_MAGIC) then
            call fail_close(2, "bad magic: not a LowFoldIn spectral file")
            return
        end if
        pos = pos + 8

        if (.not. have(SD_HEADER_BYTES - 8, "header")) return
        read(u, iostat=ios) vmaj, vmin, canary
        if (.not. ok_read("version/canary")) return
        if (canary /= SD_CANARY) then
            call fail_close(3, "endian canary mismatch (file written with different byte order?)")
            return
        end if
        if (vmaj /= SD_VER_MAJOR) then
            call fail_close(3, "unsupported version_major (expected 1)")
            return
        end if
        read(u, iostat=ios) sz
        if (.not. ok_read("type sizes")) return
        if (sz(1) /= 4 .or. sz(2) /= 8 .or. sz(3) /= 8) then
            call fail_close(3, "unexpected on-disk type sizes (expected int32=4, int64=8, real64=8)")
            return
        end if
        read(u, iostat=ios) dims, flags, ef
        if (.not. ok_read("dimensions")) return
        read(u, iostat=ios) sd%fingerprint, sd%source
        if (.not. ok_read("fingerprint/source")) return
        pos = pos + (SD_HEADER_BYTES - 8)

        sd%version_major = int(vmaj)
        sd%version_minor = int(vmin)
        sd%nk = int(dims(1)); sd%nband = int(dims(2)); sd%nspin = int(dims(3))
        sd%nregion = int(dims(4)); sd%nvertex = int(dims(5))
        sd%spin_mode = int(flags(1)); sd%method = int(flags(2)); sd%energy_ref = int(flags(3))
        sd%fermi_energy_ev = real(ef, dp)

        if (sd%nk < 1 .or. sd%nband < 1 .or. sd%nregion < 1 .or. sd%nvertex < 0) then
            call fail_close(4, "invalid dimensions (need nk>=1, nband>=1, nregion>=1, nvertex>=0)")
            return
        end if
        if (sd%nspin /= 1 .and. sd%nspin /= 2) then
            call fail_close(4, "invalid nspin (must be 1 or 2)")
            return
        end if
        if (sd%spin_mode < SPIN_NONE .or. sd%spin_mode > SPIN_NONCOLLINEAR) then
            call fail_close(4, "invalid spin_mode")
            return
        end if
        if ((sd%nspin == 2) .neqv. (sd%spin_mode == SPIN_COLLINEAR)) then
            call fail_close(4, "inconsistent nspin and spin_mode (nspin=2 iff collinear)")
            return
        end if
        if (sd%method /= METHOD_MULLIKEN .and. sd%method /= METHOD_LOWDIN) then
            call fail_close(4, "invalid method")
            return
        end if
        if (sd%energy_ref < EREF_UNKNOWN .or. sd%energy_ref > EREF_ABSOLUTE) then
            call fail_close(4, "invalid energy_ref")
            return
        end if
        if (flags(4) /= 0 .and. flags(4) /= 1) then
            call fail_close(4, "invalid fermi_known flag (must be 0 or 1)")
            return
        end if
        sd%fermi_known = (flags(4) == 1)

        ! ---- pontos k ------------------------------------------------
        ! 8 bytes * (3 kfrac + 1 kdist) por ponto k
        if (.not. have(8_int64 * 4_int64 * int(sd%nk, int64), "kfrac/kdist")) return
        allocate(sd%kfrac(3, sd%nk), sd%kdist(sd%nk))
        read(u, iostat=ios) sd%kfrac, sd%kdist
        if (.not. ok_read("kfrac/kdist")) return
        pos = pos + 32_int64 * int(sd%nk, int64)

        ! ---- vértices ------------------------------------------------
        ! 20 bytes por vértice: índice int32 + rótulo char(16)
        if (.not. have(20_int64 * int(sd%nvertex, int64), "vertices")) return
        allocate(itmp(sd%nvertex), sd%vertex_k(sd%nvertex), sd%vertex_label(sd%nvertex))
        if (sd%nvertex > 0) then
            read(u, iostat=ios) itmp
            if (.not. ok_read("vertex_k")) return
            sd%vertex_k = int(itmp)
            do i = 1, sd%nvertex
                read(u, iostat=ios) sd%vertex_label(i)
                if (.not. ok_read("vertex_label")) return
            end do
        end if
        deallocate(itmp)
        pos = pos + 20_int64 * int(sd%nvertex, int64)
        if (sd%nvertex > 0) then
            if (any(sd%vertex_k < 1) .or. any(sd%vertex_k > sd%nk)) then
                call fail_close(5, "vertex_k out of range 1..nk")
                return
            end if
        end if

        ! ---- regiões -------------------------------------------------
        ! 68 bytes por região: nome char(64) + natoms int32
        if (.not. have(68_int64 * int(sd%nregion, int64), "regions")) return
        allocate(sd%region_name(sd%nregion), sd%region_natoms(sd%nregion), itmp(sd%nregion))
        do i = 1, sd%nregion
            read(u, iostat=ios) sd%region_name(i)
            if (.not. ok_read("region_name")) return
        end do
        read(u, iostat=ios) itmp
        if (.not. ok_read("region_natoms")) return
        sd%region_natoms = int(itmp)
        deallocate(itmp)
        pos = pos + 68_int64 * int(sd%nregion, int64)
        if (trim(adjustl(sd%region_name(1))) /= 'total') then
            call fail_close(5, "region 1 must be named 'total'")
            return
        end if

        ! ---- energias e pesos ------------------------------------------
        nbytes_ext = 8_int64 * int(sd%nband, int64) * int(sd%nk, int64) * int(sd%nspin, int64)
        if (.not. have(nbytes_ext, "energy")) return
        allocate(sd%energy(sd%nband, sd%nk, sd%nspin))
        read(u, iostat=ios) sd%energy
        if (.not. ok_read("energy")) return
        pos = pos + nbytes_ext

        nbytes_ext = nbytes_ext * int(sd%nregion, int64)
        if (.not. have(nbytes_ext, "weight")) return
        allocate(sd%weight(sd%nband, sd%nk, sd%nregion, sd%nspin))
        read(u, iostat=ios) sd%weight
        if (.not. ok_read("weight")) return
        pos = pos + nbytes_ext

        ! ---- blocos de extensão ----------------------------------------
        ! Gancho de compatibilidade futura: cada bloco é {tag char(16), nbytes
        ! int64, payload}. Aqui os payloads são ignorados, nunca interpretados.
        if (.not. have(4_int64, "n_ext")) return
        read(u, iostat=ios) n_ext
        if (.not. ok_read("n_ext")) return
        pos = pos + 4
        if (n_ext < 0) then
            call fail_close(6, "negative n_ext")
            return
        end if
        do i = 1, int(n_ext)
            if (.not. have(24_int64, "extension block header")) return
            read(u, pos=pos, iostat=ios) tag, nbytes_ext
            if (.not. ok_read("extension block header")) return
            pos = pos + 24
            if (nbytes_ext < 0) then
                call fail_close(6, "extension block '" // trim(tag) // "' has negative size")
                return
            end if
            if (.not. have(nbytes_ext, "extension block payload")) return
            pos = pos + nbytes_ext
        end do

        ! ---- trailer ---------------------------------------------------
        if (.not. have(8_int64, "trailer")) return
        read(u, pos=pos, iostat=ios) trailer
        if (.not. ok_read("trailer")) return
        if (trailer /= SD_TRAILER) then
            call fail_close(7, "trailer not found at expected position (corrupted file)")
            return
        end if
        pos = pos + 8
        if (pos /= fsize + 1) then
            call fail_close(7, "unexpected trailing data after trailer")
            return
        end if

        close(u, iostat=ios)

    contains

        !> Verdadeiro se restam pelo menos n bytes a partir de pos; caso contrário
        !> fecha o arquivo e define ierr = 8.
        logical function have(n, what)
            integer(int64),   intent(in) :: n
            character(len=*), intent(in) :: what
            have = (fsize - pos + 1 >= n)
            if (.not. have) call fail_close(8, "file truncated: not enough bytes for " // what)
        end function have

        !> Verdadeiro se a última leitura teve sucesso (ios == 0); caso contrário
        !> fecha o arquivo e define ierr = 9.
        logical function ok_read(what)
            character(len=*), intent(in) :: what
            ok_read = (ios == 0)
            if (.not. ok_read) call fail_close(9, "I/O error reading " // what)
        end function ok_read

        subroutine fail_close(code, msg)
            integer,          intent(in) :: code
            character(len=*), intent(in) :: msg
            integer                      :: ios2
            close(u, iostat=ios2)
            call fail(code, msg)
        end subroutine fail_close

        subroutine fail(code, msg)
            integer,          intent(in) :: code
            character(len=*), intent(in) :: msg
            ierr = code
            errmsg = "read_spectral_file: " // msg
            call reset_sd(sd)
        end subroutine fail

    end subroutine read_spectral_file

    !> Exportação legível por humanos (somente depuração, nunca relida).
    !> Layout: cabeçalho de metadados com '#', uma linha '#' nomeando as colunas e
    !> depois uma linha por (spin, k, banda):
    !>   ispin ik kdist k1 k2 k3 iband energy w(1..nregion)
    !> kdist em 1/Angstrom, k1..k3 fracionários, energy em eV.
    !> Mesma validação e convenção de ierr de write_spectral_file
    !> (0 ok; 1 sd inválido ou não abre; 2 erro de E/S; 3 erro no fechamento).
    subroutine write_spectral_text(filename, sd, ierr, errmsg)
        character(len=*),      intent(in)  :: filename
        type(spectral_data_t), intent(in)  :: sd
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        integer                       :: u, ios, i, is, ik, ib
        character(len=:), allocatable :: line

        call validate(sd, ierr, errmsg)
        if (ierr /= 0) return

        open(newunit=u, file=filename, access='sequential', form='formatted', &
             status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            errmsg = "write_spectral_text: cannot open '" // trim(filename) // "' for writing"
            return
        end if

        write(u, '(A)', iostat=ios) "# LowFoldIn spectral data (text export, not read back)"
        if (ios == 0) write(u, '(A,I0,A,I0)', iostat=ios) "# version ", &
            sd%version_major, ".", sd%version_minor
        if (ios == 0) write(u, '(A,I0,A,I0,A,I0,A,I0,A,I0)', iostat=ios) &
            "# nk ", sd%nk, "  nband ", sd%nband, "  nspin ", sd%nspin, &
            "  nregion ", sd%nregion, "  nvertex ", sd%nvertex
        if (ios == 0) write(u, '(A,A)', iostat=ios) "# spin_mode ", spin_mode_name(sd%spin_mode)
        if (ios == 0) write(u, '(A,A)', iostat=ios) "# method ", method_name(sd%method)
        if (ios == 0) write(u, '(A,A)', iostat=ios) "# energy_reference ", eref_name(sd%energy_ref)
        if (ios == 0) write(u, '(A,L1,A,ES23.15)', iostat=ios) "# fermi_known ", sd%fermi_known, &
            "  fermi_energy_ev ", sd%fermi_energy_ev
        if (ios == 0) write(u, '(A,A)', iostat=ios) "# fingerprint ", trim(sd%fingerprint)
        if (ios == 0) write(u, '(A,A)', iostat=ios) "# source ", trim(sd%source)
        if (ios == 0) write(u, '(A)', iostat=ios) "# units energy eV, k_distance 1/Angstrom"
        do i = 1, sd%nvertex
            if (ios == 0) write(u, '(A,I0,A,A,A,I0,A,ES23.15)', iostat=ios) "# vertex ", i, &
                " label ", trim(sd%vertex_label(i)), " ik ", sd%vertex_k(i), &
                " kdist ", sd%kdist(sd%vertex_k(i))
        end do
        do i = 1, sd%nregion
            if (ios == 0) write(u, '(A,I0,A,A,A,I0)', iostat=ios) "# region ", i, &
                " name ", trim(sd%region_name(i)), " natoms ", sd%region_natoms(i)
        end do
        line = "# ispin ik kdist k1 k2 k3 iband energy"
        do i = 1, sd%nregion
            line = line // " " // trim(adjustl(sd%region_name(i)))
        end do
        if (ios == 0) write(u, '(A)', iostat=ios) line

        do is = 1, sd%nspin
            do ik = 1, sd%nk
                do ib = 1, sd%nband
                    if (ios /= 0) exit
                    write(u, '(I0,1X,I0,4(1X,ES23.15),1X,I0,1X,ES23.15,*(1X,ES23.15))', &
                        iostat=ios) is, ik, sd%kdist(ik), sd%kfrac(:, ik), ib, &
                        sd%energy(ib, ik, is), sd%weight(ib, ik, :, is)
                end do
            end do
        end do

        if (ios /= 0) then
            close(u, iostat=i)
            ierr = 2
            errmsg = "write_spectral_text: I/O error while writing '" // trim(filename) // "'"
            return
        end if
        close(u, iostat=ios)
        if (ios /= 0) then
            ierr = 3
            errmsg = "write_spectral_text: error closing '" // trim(filename) // "'"
            return
        end if
        ierr = 0
        errmsg = ''
    end subroutine write_spectral_text

    ! ==================================================================
    ! Internos
    ! ==================================================================

    !> Desaloca todos os arrays de sd e restaura os escalares padrão, de modo que
    !> uma leitura que falhou não deixe nada alocado.
    subroutine reset_sd(sd)
        type(spectral_data_t), intent(inout) :: sd
        sd = spectral_data_t()
    end subroutine reset_sd

    !> Verificações de consistência compartilhadas pelos escritores: dimensões,
    !> enumerações, concordância entre nspin/spin_mode, alocação e formatos dos
    !> arrays, faixa dos vértices e região 1 = 'total'. Espelha as verificações que
    !> read_spectral_file aplica ao conteúdo do arquivo. Define ierr = 1 e errmsg
    !> na primeira falha.
    subroutine validate(sd, ierr, errmsg)
        type(spectral_data_t), intent(in)  :: sd
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        ierr = 1
        errmsg = ''
        if (sd%nk < 1 .or. sd%nband < 1 .or. sd%nregion < 1 .or. sd%nvertex < 0) then
            errmsg = "spectral_data: invalid dimensions (need nk>=1, nband>=1, nregion>=1, nvertex>=0)"
            return
        end if
        if (sd%nspin /= 1 .and. sd%nspin /= 2) then
            errmsg = "spectral_data: nspin must be 1 or 2"
            return
        end if
        if (sd%spin_mode < SPIN_NONE .or. sd%spin_mode > SPIN_NONCOLLINEAR) then
            errmsg = "spectral_data: invalid spin_mode"
            return
        end if
        if ((sd%nspin == 2) .neqv. (sd%spin_mode == SPIN_COLLINEAR)) then
            errmsg = "spectral_data: nspin=2 iff spin_mode=SPIN_COLLINEAR"
            return
        end if
        if (sd%method /= METHOD_MULLIKEN .and. sd%method /= METHOD_LOWDIN) then
            errmsg = "spectral_data: invalid method"
            return
        end if
        if (sd%energy_ref < EREF_UNKNOWN .or. sd%energy_ref > EREF_ABSOLUTE) then
            errmsg = "spectral_data: invalid energy_ref"
            return
        end if
        if (.not. (allocated(sd%kfrac) .and. allocated(sd%kdist) .and. &
                   allocated(sd%region_name) .and. allocated(sd%region_natoms) .and. &
                   allocated(sd%energy) .and. allocated(sd%weight))) then
            errmsg = "spectral_data: required arrays not allocated"
            return
        end if
        if (size(sd%kfrac, 1) /= 3 .or. size(sd%kfrac, 2) /= sd%nk .or. size(sd%kdist) /= sd%nk) then
            errmsg = "spectral_data: kfrac/kdist shape mismatch"
            return
        end if
        if (size(sd%region_name) /= sd%nregion .or. size(sd%region_natoms) /= sd%nregion) then
            errmsg = "spectral_data: region arrays shape mismatch"
            return
        end if
        if (any(shape(sd%energy) /= [sd%nband, sd%nk, sd%nspin])) then
            errmsg = "spectral_data: energy shape mismatch"
            return
        end if
        if (any(shape(sd%weight) /= [sd%nband, sd%nk, sd%nregion, sd%nspin])) then
            errmsg = "spectral_data: weight shape mismatch"
            return
        end if
        if (sd%nvertex > 0) then
            if (.not. (allocated(sd%vertex_k) .and. allocated(sd%vertex_label))) then
                errmsg = "spectral_data: vertex arrays not allocated"
                return
            end if
        end if
        if (allocated(sd%vertex_k)) then
            if (size(sd%vertex_k) /= sd%nvertex) then
                errmsg = "spectral_data: vertex_k shape mismatch"
                return
            end if
            if (any(sd%vertex_k < 1) .or. any(sd%vertex_k > sd%nk)) then
                errmsg = "spectral_data: vertex_k out of range 1..nk"
                return
            end if
        end if
        if (allocated(sd%vertex_label)) then
            if (size(sd%vertex_label) /= sd%nvertex) then
                errmsg = "spectral_data: vertex_label shape mismatch"
                return
            end if
        end if
        if (trim(adjustl(sd%region_name(1))) /= 'total') then
            errmsg = "spectral_data: region 1 must be named 'total'"
            return
        end if
        ierr = 0
    end subroutine validate

    ! Enumeração -> nome em texto, usado pela exportação em texto ("invalid" se desconhecida).
    pure function spin_mode_name(m) result(s)
        integer, intent(in)           :: m
        character(len=:), allocatable :: s
        select case (m)
        case (SPIN_NONE);         s = "none"
        case (SPIN_COLLINEAR);    s = "collinear"
        case (SPIN_NONCOLLINEAR); s = "noncollinear"
        case default;             s = "invalid"
        end select
    end function spin_mode_name

    pure function method_name(m) result(s)
        integer, intent(in)           :: m
        character(len=:), allocatable :: s
        select case (m)
        case (METHOD_MULLIKEN); s = "mulliken"
        case (METHOD_LOWDIN);   s = "lowdin"
        case default;           s = "invalid"
        end select
    end function method_name

    pure function eref_name(m) result(s)
        integer, intent(in)           :: m
        character(len=:), allocatable :: s
        select case (m)
        case (EREF_UNKNOWN);  s = "unknown"
        case (EREF_FERMI);    s = "fermi"
        case (EREF_ABSOLUTE); s = "absolute"
        case default;         s = "invalid"
        end select
    end function eref_name

end module spectral_data
