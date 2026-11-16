! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! mod_scfout_io.f90
!
! Leitor do formato binário .scfout do OpenMX (versão SCFOUT_VERSION = 3,
! OpenMX 3.9/4.0), porte direto de read_scfout.c (T. Ozaki / H. Kotaka /
! N. Yamaguchi) para Fortran, usado como front-end da ferramenta
! tools/openmx2hsio.f90.
!
! A ordem de leitura segue EXATAMENTE a funcao Input() de read_scfout.c:
!
!   cabeçalho(6 int) -> order_max -> atv -> atv_ijk -> Total_NumOrbs
!   -> FNAN -> natn -> ncn -> tv -> rtv -> Gxyz -> Hks -> iHks (se SO)
!   -> OLP -> OLPpo -> OLPmo -> DM -> iDM -> Solver -> ChemP/E_Temp/...
!   -> (arquivo de input embutido: só usado para os rótulos das especies)
!
! Geometria: tv (rede direta, linhas = a1,a2,a3) e Gxyz (posicoes
! cartesianas) são guardados em Bohr (lattice_bohr, cart_bohr); a conversão
! para Ångström fica a cargo de quem chama. Os rótulos de espécie vêm do
! bloco <Atoms.SpeciesAndCoordinates do input embutido (atom_label).
!
! OLPpo, OLPmo, DM e iDM (overlap com operador posição/momento, matriz de
! densidade) não são usados por este conversor; seus bytes são pulados
! via leitura-e-descarte em blocos, sem alocar nada, após calcular o
! número exato de elementos pela mesma estrutura átomo-vizinho usada em
! Hks/OLP.
!
! Nathanael N. Batista / porte para Fortran, Agosto de 2026
! ----------------------------------------------------------------------------------------------------------
module scfout_io
    use types, only: dp
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private

    public :: ScfoutData, ReadScfoutFile, ParseEmbeddedAtomLabels

    integer, parameter :: SCFOUT_VERSION = 3
    integer, parameter :: SKIP_CHUNK     = 65536   ! doubles por bloco ao descartar OLPpo/OLPmo/DM/iDM

    type :: ScfoutData
        integer  :: atomnum      = 0
        integer  :: SpinP_switch = 0
        integer  :: version      = 0
        integer  :: Catomnum     = 0
        integer  :: Latomnum     = 0
        integer  :: Ratomnum     = 0
        integer  :: TCpyCell     = 0
        integer  :: order_max    = 0
        integer  :: Solver       = 0
        real(dp) :: ChemP        = 0.0_dp
        real(dp) :: E_Temp       = 0.0_dp

        ! Geometria do .scfout (NÃO convertida: o OpenMX guarda em Bohr).
        ! As linhas de lattice_bohr são os vetores da rede direta a1,a2,a3;
        ! as colunas de cart_bohr são átomos (mesma ordem dos blocos do Hamiltoniano).
        real(dp)              :: lattice_bohr(3, 3) = 0.0_dp
        real(dp), allocatable :: cart_bohr(:, :)              ! (3, atomnum), Bohr
        ! Rótulo de espécie por átomo, extraído do input embutido do OpenMX.
        ! atom_label só é alocado se esse input pôde ser lido, e seu conteúdo
        ! só tem significado quando has_atom_labels é .true.
        character(len=8), allocatable :: atom_label(:)
        logical :: has_atom_labels = .false.

        integer, allocatable :: TotalNumOrbs(:)     ! 1:atomnum
        integer, allocatable :: FNAN(:)             ! 1:atomnum

        ! Listas de vizinhos em formato "flat" (estilo CSR): o átomo ct_AN
        ! tem FNAN(ct_AN)+1 vizinhos (h_AN=0..FNAN(ct_AN)), armazenados em
        ! natn_flat(neighbor_offset(ct_AN)+1 : neighbor_offset(ct_AN)+FNAN(ct_AN)+1)
        integer, allocatable :: neighbor_offset(:)  ! 1:atomnum (base 0)
        integer, allocatable :: natn_flat(:)
        integer, allocatable :: ncn_flat(:)
        integer              :: n_pairs = 0

        integer, allocatable :: atv_ijk(:, :)       ! (3, 0:TCpyCell) -> IndexList (h,k,l)

        ! Blocos densos H(R)/S(R) por par (ct_AN,h_AN), em formato "flat"
        ! (CSR de novo): block_offset(p) dá o offset (0-based, em
        ! elementos escalares) do par p=1..n_pairs dentro de Hks_flat(:,ispin)
        ! ou OLP_flat(:). block_offset(n_pairs+1) = número total de
        ! elementos escalares (n_scalar).
        integer, allocatable :: block_offset(:)     ! 1:n_pairs+1
        integer              :: n_scalar = 0

        real(dp), allocatable :: Hks_flat(:, :)     ! (n_scalar, 0:SpinP_switch)
        real(dp), allocatable :: iHks_flat(:, :)    ! (n_scalar, 0:2), só se SpinP_switch==3
        real(dp), allocatable :: OLP_flat(:)        ! (n_scalar)
    end type ScfoutData

contains

    ! Lê o .scfout inteiro em `this` (stream binário, ordem de Input() de
    ! read_scfout.c). Não retorna erro: arquivo ausente, versão /= 3 ou EOF
    ! inesperado fazem `stop 1` com mensagem. Única exceção: o input embutido
    ! no final é opcional (ver has_atom_labels).
    subroutine ReadScfoutFile(this, filename)
        type(ScfoutData), intent(inout) :: this
        character(len=*), intent(in)    :: filename

        integer  :: iunit, ios
        integer  :: i_vec(6)
        integer  :: packed
        integer  :: ct_AN, h_AN, p, off
        integer  :: Rn, ncell
        integer  :: TNO1, TNO2, Gh_AN
        integer  :: n_block
        real(dp) :: d_vec(10)

        open(newunit=iunit, file=filename, status='old', action='read', &
             iostat=ios, form='unformatted', access='stream')
        if (ios /= 0) then
            print *, "[scfout_io] Erro: nao foi possivel abrir ", trim(filename)
            stop 1
        end if

        ! -------------------------------------------------------------
        ! Cabecalho: atomnum, SpinP_switch+versão*4, Catomnum, Latomnum,
        ! Ratomnum, TCpyCell (6 inteiros)
        ! -------------------------------------------------------------
        read(iunit, iostat=ios) i_vec
        call check_io(ios, 'cabecalho')
        this%atomnum      = i_vec(1)
        packed            = i_vec(2)
        this%SpinP_switch = mod(packed, 4)
        this%version      = packed / 4
        this%Catomnum     = i_vec(3)
        this%Latomnum     = i_vec(4)
        this%Ratomnum     = i_vec(5)
        this%TCpyCell     = i_vec(6)

        if (this%version /= SCFOUT_VERSION) then
            print '(A,I0)', "[scfout_io] Erro: versao do .scfout = ", this%version
            print '(A,I0)', "[scfout_io] Este leitor so suporta a versao ", SCFOUT_VERSION
            print '(A)',    "[scfout_io] (OpenMX 3.9/4.0). Avise se precisar de outra versao."
            stop 1
        end if

        ! Valida o cabeçalho antes de usá-lo como tamanho/índice
        if (this%atomnum < 1) call fail_header('atomnum', this%atomnum)
        if (this%TCpyCell < 0) call fail_header('TCpyCell', this%TCpyCell)
        if (this%SpinP_switch /= 0 .and. this%SpinP_switch /= 1 .and. this%SpinP_switch /= 3) &
            call fail_header('SpinP_switch', this%SpinP_switch)

        ! -------------------------------------------------------------
        ! order_max
        ! -------------------------------------------------------------
        block
            integer :: tmp(1)
            read(iunit, iostat=ios) tmp
            call check_io(ios, 'order_max')
            this%order_max = tmp(1)
        end block

        ncell = this%TCpyCell + 1

        ! -------------------------------------------------------------
        ! atv[0..TCpyCell][4] (double) -- não usado por este conversor,
        ! mas precisa ser consumido do stream na ordem certa
        ! -------------------------------------------------------------
        block
            real(dp) :: atv_junk(4, ncell)
            read(iunit, iostat=ios) atv_junk
            call check_io(ios, 'atv')
        end block

        ! -------------------------------------------------------------
        ! atv_ijk[0..TCpyCell][4] (int) -- guardamos só (i,j,k), o
        ! primeiro elemento de cada linha (índice da propria celula) é
        ! descartado
        ! -------------------------------------------------------------
        block
            integer :: atv_ijk_raw(4, ncell)
            read(iunit, iostat=ios) atv_ijk_raw
            call check_io(ios, 'atv_ijk')
            allocate(this%atv_ijk(3, 0:this%TCpyCell))
            do Rn = 0, this%TCpyCell
                this%atv_ijk(:, Rn) = atv_ijk_raw(2:4, Rn + 1)
            end do
        end block

        ! -------------------------------------------------------------
        ! Total_NumOrbs[1..atomnum], FNAN[1..atomnum]
        ! -------------------------------------------------------------
        allocate(this%TotalNumOrbs(this%atomnum))
        read(iunit, iostat=ios) this%TotalNumOrbs
        call check_io(ios, 'Total_NumOrbs')
        if (any(this%TotalNumOrbs < 1)) call fail_header('Total_NumOrbs', minval(this%TotalNumOrbs))

        allocate(this%FNAN(this%atomnum))
        read(iunit, iostat=ios) this%FNAN
        call check_io(ios, 'FNAN')
        if (any(this%FNAN < 0)) call fail_header('FNAN', minval(this%FNAN))

        ! -------------------------------------------------------------
        ! natn[][], ncn[][] -- ragged, guardados em formato flat (CSR)
        ! -------------------------------------------------------------
        allocate(this%neighbor_offset(this%atomnum))
        off = 0
        do ct_AN = 1, this%atomnum
            this%neighbor_offset(ct_AN) = off
            off = off + this%FNAN(ct_AN) + 1
        end do
        this%n_pairs = off

        allocate(this%natn_flat(this%n_pairs))
        do ct_AN = 1, this%atomnum
            off = this%neighbor_offset(ct_AN)
            read(iunit, iostat=ios) this%natn_flat(off + 1 : off + this%FNAN(ct_AN) + 1)
            call check_io(ios, 'natn')
            if (any(this%natn_flat(off + 1 : off + this%FNAN(ct_AN) + 1) < 1 .or. &
                    this%natn_flat(off + 1 : off + this%FNAN(ct_AN) + 1) > this%atomnum)) &
                call fail_header('natn (indice de vizinho fora de 1..atomnum), atomo', ct_AN)
        end do

        allocate(this%ncn_flat(this%n_pairs))
        do ct_AN = 1, this%atomnum
            off = this%neighbor_offset(ct_AN)
            read(iunit, iostat=ios) this%ncn_flat(off + 1 : off + this%FNAN(ct_AN) + 1)
            call check_io(ios, 'ncn')
            if (any(this%ncn_flat(off + 1 : off + this%FNAN(ct_AN) + 1) < 0 .or. &
                    this%ncn_flat(off + 1 : off + this%FNAN(ct_AN) + 1) > this%TCpyCell)) &
                call fail_header('ncn (indice de celula fora de 0..TCpyCell), atomo', ct_AN)
        end do

        ! -------------------------------------------------------------
        ! tv[1..3][4], rtv[1..3][4]. O OpenMX guarda os vetores da rede direta
        ! como linhas, em Bohr; o primeiro valor de cada linha não é usado.
        ! -------------------------------------------------------------
        block
            real(dp) :: tv_raw(4, 3), rtv_junk(4, 3)
            read(iunit, iostat=ios) tv_raw
            call check_io(ios, 'tv')
            read(iunit, iostat=ios) rtv_junk
            call check_io(ios, 'rtv')
            do Rn = 1, 3
                this%lattice_bohr(Rn, :) = tv_raw(2:4, Rn)
            end do
        end block

        ! -------------------------------------------------------------
        ! Gxyz[1..atomnum][0..3]: coordenadas cartesianas em Bohr. O primeiro
        ! valor de cada registro de átomo não é usado.
        ! -------------------------------------------------------------
        block
            real(dp) :: gxyz_raw(4, this%atomnum)
            read(iunit, iostat=ios) gxyz_raw
            call check_io(ios, 'Gxyz')
            allocate(this%cart_bohr(3, this%atomnum))
            this%cart_bohr = gxyz_raw(2:4, :)
        end block

        ! -------------------------------------------------------------
        ! block_offset: offset (0-based, em escalares) de cada par
        ! (ct_AN,h_AN) dentro dos blocos Hks/iHks/OLP -- calculado ANTES
        ! de ler esses blocos, pois só depende de natn/TotalNumOrbs
        ! -------------------------------------------------------------
        allocate(this%block_offset(this%n_pairs + 1))
        p = 0
        n_block = 0
        do ct_AN = 1, this%atomnum
            TNO1 = this%TotalNumOrbs(ct_AN)
            do h_AN = 0, this%FNAN(ct_AN)
                p = p + 1
                Gh_AN = this%natn_flat(this%neighbor_offset(ct_AN) + h_AN + 1)
                TNO2 = this%TotalNumOrbs(Gh_AN)
                this%block_offset(p) = n_block
                n_block = n_block + TNO1 * TNO2
            end do
        end do
        this%block_offset(this%n_pairs + 1) = n_block
        this%n_scalar = n_block

        ! -------------------------------------------------------------
        ! Hks[spin][ct_AN][h_AN][TNO1][TNO2], spin=0..SpinP_switch
        ! Cada par (ct_AN,h_AN) é TNO1*TNO2 doubles CONTÍGUOS no stream
        ! (loop "for i in TNO1: read TNO2 doubles" = flatten row-major),
        ! e os pares são lidos em sequência -> uma única leitura de
        ! n_scalar doubles cobre o canal de spin inteiro.
        ! -------------------------------------------------------------
        allocate(this%Hks_flat(this%n_scalar, 0:this%SpinP_switch))
        block
            integer :: spin
            do spin = 0, this%SpinP_switch
                read(iunit, iostat=ios) this%Hks_flat(:, spin)
                call check_io(ios, 'Hks')
            end do
        end block

        ! -------------------------------------------------------------
        ! iHks[spin][...], spin=0..2 -- só existe se SpinP_switch==3
        ! -------------------------------------------------------------
        if (this%SpinP_switch == 3) then
            allocate(this%iHks_flat(this%n_scalar, 0:2))
            block
                integer :: spin
                do spin = 0, 2
                    read(iunit, iostat=ios) this%iHks_flat(:, spin)
                    call check_io(ios, 'iHks')
                end do
            end block
        end if

        ! -------------------------------------------------------------
        ! OLP[ct_AN][h_AN][TNO1][TNO2] -- um único canal, mesma logica
        ! -------------------------------------------------------------
        allocate(this%OLP_flat(this%n_scalar))
        read(iunit, iostat=ios) this%OLP_flat
        call check_io(ios, 'OLP')

        ! -------------------------------------------------------------
        ! Blocos que NÃO precisamos materializar, mas que ficam no meio
        ! do caminho até ChemP: OLPpo, OLPmo, DM, iDM. Todos seguem a
        ! MESMA estrutura átomo-vizinho (n_scalar elementos por "canal"),
        ! só que repetida por direção (x,y,z), por ordem (OLPpo), ou por
        ! canal de spin (DM/iDM). Pulamos lendo-e-descartando em blocos.
        ! -------------------------------------------------------------

        ! OLPpo: 3 direções * order_max ordens * n_scalar (vazio se order_max==0)
        call skip_doubles(iunit, 3_int64 * int(this%order_max, int64) * int(this%n_scalar, int64))

        ! OLPmo: 3 direções * n_scalar -- SEMPRE presente, independente de order_max
        call skip_doubles(iunit, 3_int64 * int(this%n_scalar, int64))

        ! DM: (SpinP_switch+1) * n_scalar
        call skip_doubles(iunit, int(this%SpinP_switch + 1, int64) * int(this%n_scalar, int64))

        ! iDM: 2 * n_scalar
        call skip_doubles(iunit, 2_int64 * int(this%n_scalar, int64))

        ! -------------------------------------------------------------
        ! Solver
        ! -------------------------------------------------------------
        block
            integer :: tmp(1)
            read(iunit, iostat=ios) tmp
            call check_io(ios, 'Solver')
            this%Solver = tmp(1)
        end block

        ! -------------------------------------------------------------
        ! ChemP, E_Temp, dipolos, Valence_Electrons, Total_SpinS
        ! (10 doubles, ordem exata de read_scfout.c: d_vec[0..9])
        ! -------------------------------------------------------------
        read(iunit, iostat=ios) d_vec
        call check_io(ios, 'ChemP/E_Temp')
        this%ChemP  = d_vec(1)
        this%E_Temp = d_vec(2)
        ! d_vec(3:5) = dipole_moment_core, d_vec(6:8) = dipole_moment_background,
        ! d_vec(9) = Valence_Electrons, d_vec(10) = Total_SpinS -- não usados aqui.

        ! -------------------------------------------------------------
        ! Input embutido do OpenMX: num_lines (int32) seguido de registros fixos
        ! de 256 bytes, anexados após os dados numéricos. Fornece os rótulos de
        ! espécie, que não estão nos próprios registros binários de geometria.
        ! Input ausente/malformado (EOF, num_lines fora de 1..100000, falha de
        ! alocação) não é erro: a conversão de H/S continua funcionando e
        ! has_atom_labels permanece falso.
        ! -------------------------------------------------------------
        block
            integer :: num_lines, iline, alloc_stat
            character(len=256), allocatable :: input_lines(:)

            num_lines = 0
            read(iunit, iostat=ios) num_lines
            if (ios == 0) then
                if (num_lines > 0 .and. num_lines <= 100000) then
                    allocate(input_lines(num_lines), stat=alloc_stat)
                    if (alloc_stat == 0) then
                        do iline = 1, num_lines
                            read(iunit, iostat=ios) input_lines(iline)
                            if (ios /= 0) exit
                        end do
                        if (ios == 0) then
                            allocate(this%atom_label(this%atomnum))
                            call ParseEmbeddedAtomLabels(input_lines, this%atomnum, this%atom_label, &
                                                         this%has_atom_labels)
                        end if
                    end if
                end if
            end if
        end block

        close(iunit, iostat=ios)
        if (ios /= 0) then
            print *, "[scfout_io] Erro: falha ao fechar ", trim(filename)
            stop 1
        end if
    end subroutine ReadScfoutFile

    !> Extrai os rótulos de espécie por átomo do texto do input embutido do OpenMX.
    !>
    !> Procura o bloco delimitado por `<Atoms.SpeciesAndCoordinates` e
    !> `Atoms.SpeciesAndCoordinates>` (sem distinguir maiúsculas) e lê, de cada
    !> linha interna, o índice do átomo (base 1) e o rótulo da espécie (os dois
    !> primeiros campos em formato livre; coordenadas etc. são ignoradas). Linhas
    !> em branco, linhas de comentário '#'/'!' e o preenchimento com NUL (os
    !> registros têm 256 bytes fixos) são ignorados.
    !>
    !> `found` é .true. somente se o bloco foi fechado e cada átomo 1..nat
    !> recebeu exatamente um rótulo de 1..8 caracteres. Em qualquer linha
    !> malformada (índice inválido, duplicado, rótulo longo demais, ilegível)
    !> retorna cedo com found = .false. e `labels` parcialmente preenchido: não
    !> use os rótulos nesse caso. Nunca para o programa.
    subroutine ParseEmbeddedAtomLabels(lines, nat, labels, found)
        character(len=*), intent(in)  :: lines(:)
        integer,          intent(in)  :: nat
        character(len=8), intent(out) :: labels(nat)
        logical,          intent(out) :: found

        character(len=256)   :: line
        character(len=64)    :: label
        logical, allocatable :: seen(:)
        logical              :: in_atoms, closed
        integer              :: i, idx, ios, nul

        labels = ''
        found = .false.
        in_atoms = .false.
        closed = .false.
        allocate(seen(nat))
        seen = .false.

        do i = 1, size(lines)
            line = lines(i)
            nul = index(line, achar(0))
            if (nul > 0) line(nul:) = ' '
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle

            if (.not. in_atoms) then
                if (index(to_lower(line), '<atoms.speciesandcoordinates') > 0) in_atoms = .true.
                cycle
            end if

            if (index(to_lower(line), 'atoms.speciesandcoordinates>') > 0) then
                closed = .true.
                exit
            end if

            idx = 0
            label = ''
            read(line, *, iostat=ios) idx, label
            if (ios /= 0) return
            if (idx < 1 .or. idx > nat) return
            if (seen(idx)) return
            if (len_trim(label) < 1 .or. len_trim(label) > len(labels(idx))) return
            labels(idx) = trim(label)
            seen(idx) = .true.
        end do

        found = closed .and. all(seen)
    end subroutine ParseEmbeddedAtomLabels

    !> Cópia de text em minúsculas ASCII (usada na busca de palavras-chave sem distinguir maiúsculas).
    pure function to_lower(text) result(lower)
        character(len=*), intent(in) :: text
        character(len=len(text))     :: lower
        integer                      :: i, code

        lower = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
        end do
    end function to_lower

    ! Aborta com mensagem se uma leitura do stream falhou (EOF/arquivo truncado).
    subroutine check_io(ios, what)
        integer,          intent(in) :: ios
        character(len=*), intent(in) :: what

        if (ios /= 0) then
            print '(A)', "[scfout_io] Erro: falha lendo " // what // " (EOF/arquivo truncado?)"
            stop 1
        end if
    end subroutine check_io

    ! Aborta com mensagem se um valor do cabeçalho/vizinhança é inválido.
    subroutine fail_header(what, val)
        character(len=*), intent(in) :: what
        integer,          intent(in) :: val

        print '(A,I0)', "[scfout_io] Erro: valor invalido em " // what // " = ", val
        print '(A)',    "[scfout_io] (arquivo corrompido ou nao e' um .scfout valido?)"
        stop 1
    end subroutine fail_header

    ! Lê e descarta n doubles do stream, em blocos, sem alocar um buffer
    ! do tamanho total (que pode ser grande: OLPmo/DM/iDM são O(NOrbs^2)).
    ! Usa iostat para detectar EOF inesperado (cálculo de offset errado)
    ! em vez de falhar silenciosamente lendo ChemP do lugar errado.
    subroutine skip_doubles(iunit, n)
        integer,        intent(in) :: iunit
        integer(int64), intent(in) :: n

        real(dp)       :: buf(SKIP_CHUNK)
        integer(int64) :: remaining
        integer        :: this_chunk, ios

        if (n <= 0_int64) return

        remaining = n
        do while (remaining > 0_int64)
            this_chunk = int(min(remaining, int(SKIP_CHUNK, int64)))
            read(iunit, iostat=ios) buf(1:this_chunk)
            if (ios /= 0) then
                print *, "[scfout_io] Erro: EOF inesperado pulando blocos " // &
                         "OLPpo/OLPmo/DM/iDM (calculo de offset pode estar errado)."
                stop 1
            end if
            remaining = remaining - this_chunk
        end do
    end subroutine skip_doubles

end module scfout_io
