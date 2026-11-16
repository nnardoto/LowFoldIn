! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! mod_scfout_io.f90
!
! Leitor do formato binario .scfout do OpenMX (versao SCFOUT_VERSION = 3,
! OpenMX 3.9/4.0), porte direto de read_scfout.c (T. Ozaki / H. Kotaka /
! N. Yamaguchi) para Fortran, usado como front-end da ferramenta
! tools/openmx2hsio.f90.
!
! A ordem de leitura segue EXATAMENTE a funcao Input() de read_scfout.c:
!
!   cabecalho(6 int) -> order_max -> atv -> atv_ijk -> Total_NumOrbs
!   -> FNAN -> natn -> ncn -> tv -> rtv -> Gxyz -> Hks -> iHks (se SO)
!   -> OLP -> OLPpo -> OLPmo -> DM -> iDM -> Solver -> ChemP/E_Temp/...
!   -> (arquivo de input embutido, ignorado)
!
! OLPpo, OLPmo, DM e iDM (overlap com operador posicao/momento, matriz de
! densidade) nao sao usados por este conversor; seus bytes sao pulados
! via leitura-e-descarte em blocos, sem alocar nada, apos calcular o
! numero exato de elementos pela mesma estrutura atomo-vizinho usada em
! Hks/OLP.
!
! Nathanael N. Batista / porte para Fortran, Agosto de 2026
! ----------------------------------------------------------------------------------------------------------
module scfout_io
    use types, only: dp
    implicit none
    private

    public :: ScfoutData, ReadScfoutFile, ParseEmbeddedAtomLabels

    integer, parameter :: SCFOUT_VERSION = 3
    integer, parameter :: SKIP_CHUNK = 65536   ! doubles por bloco ao descartar OLPpo/OLPmo/DM/iDM

    type :: ScfoutData
        integer :: atomnum       = 0
        integer :: SpinP_switch  = 0
        integer :: version       = 0
        integer :: Catomnum      = 0
        integer :: Latomnum      = 0
        integer :: Ratomnum      = 0
        integer :: TCpyCell      = 0
        integer :: order_max     = 0
        integer :: Solver        = 0
        real(dp) :: ChemP        = 0.0_dp
        real(dp) :: E_Temp       = 0.0_dp

        ! OpenMX stores both quantities in Bohr.  Rows of lattice_bohr are
        ! the direct lattice vectors; columns of cart_bohr are atoms.
        real(dp) :: lattice_bohr(3, 3) = 0.0_dp
        real(dp), allocatable :: cart_bohr(:, :)       ! (3, atomnum)
        character(len=8), allocatable :: atom_label(:)
        logical :: has_atom_labels = .false.

        integer, allocatable :: TotalNumOrbs(:)     ! 1:atomnum
        integer, allocatable :: FNAN(:)             ! 1:atomnum

        ! Listas de vizinhos em formato "flat" (estilo CSR): o atomo ct_AN
        ! tem FNAN(ct_AN)+1 vizinhos (h_AN=0..FNAN(ct_AN)), armazenados em
        ! natn_flat(neighbor_offset(ct_AN)+1 : neighbor_offset(ct_AN)+FNAN(ct_AN)+1)
        integer, allocatable :: neighbor_offset(:)  ! 1:atomnum (0-based)
        integer, allocatable :: natn_flat(:)
        integer, allocatable :: ncn_flat(:)
        integer :: n_pairs = 0

        integer, allocatable :: atv_ijk(:, :)       ! (3, 0:TCpyCell) -> IndexList (h,k,l)

        ! Blocos densos H(R)/S(R) por par (ct_AN,h_AN), em formato "flat"
        ! (CSR de novo): block_offset(p) da' o offset (0-based, em
        ! elementos escalares) do par p=1..n_pairs dentro de Hks_flat(:,ispin)
        ! ou OLP_flat(:). block_offset(n_pairs+1) = numero total de
        ! elementos escalares (n_scalar).
        integer, allocatable :: block_offset(:)     ! 1:n_pairs+1
        integer :: n_scalar = 0

        real(dp), allocatable :: Hks_flat(:, :)     ! (n_scalar, 0:SpinP_switch)
        real(dp), allocatable :: iHks_flat(:, :)    ! (n_scalar, 0:2), so' se SpinP_switch==3
        real(dp), allocatable :: OLP_flat(:)        ! (n_scalar)
    end type ScfoutData

contains

    subroutine ReadScfoutFile(this, filename)
        type(ScfoutData), intent(inout) :: this
        character(len=*),  intent(in) :: filename

        integer :: iunit, ios
        integer :: i_vec(6)
        integer :: packed
        integer :: ct_AN, h_AN, p, off
        integer :: Rn, ncell
        integer :: TNO1, TNO2, Gh_AN
        integer :: n_block
        real(dp) :: d_vec(10)

        open(newunit=iunit, file=filename, status='old', action='read', &
             iostat=ios, form='unformatted', access='stream')
        if (ios /= 0) then
            print *, "[scfout_io] Erro: nao foi possivel abrir ", trim(filename)
            stop 1
        end if

        ! -------------------------------------------------------------
        ! Cabecalho: atomnum, SpinP_switch+versao*4, Catomnum, Latomnum,
        ! Ratomnum, TCpyCell (6 inteiros)
        ! -------------------------------------------------------------
        read(iunit) i_vec
        this%atomnum  = i_vec(1)
        packed        = i_vec(2)
        this%SpinP_switch = mod(packed, 4)
        this%version       = packed / 4
        this%Catomnum = i_vec(3)
        this%Latomnum = i_vec(4)
        this%Ratomnum = i_vec(5)
        this%TCpyCell = i_vec(6)

        if (this%version /= SCFOUT_VERSION) then
            print '(A,I0)', "[scfout_io] Erro: versao do .scfout = ", this%version
            print '(A,I0)', "[scfout_io] Este leitor so suporta a versao ", SCFOUT_VERSION
            print '(A)',    "[scfout_io] (OpenMX 3.9/4.0). Avise se precisar de outra versao."
            stop 1
        end if

        ! -------------------------------------------------------------
        ! order_max
        ! -------------------------------------------------------------
        block
            integer :: tmp(1)
            read(iunit) tmp
            this%order_max = tmp(1)
        end block

        ncell = this%TCpyCell + 1

        ! -------------------------------------------------------------
        ! atv[0..TCpyCell][4] (double) -- nao usado por este conversor,
        ! mas precisa ser consumido do stream na ordem certa
        ! -------------------------------------------------------------
        block
            real(dp) :: atv_junk(4, ncell)
            read(iunit) atv_junk
        end block

        ! -------------------------------------------------------------
        ! atv_ijk[0..TCpyCell][4] (int) -- guardamos so' (i,j,k), o
        ! primeiro elemento de cada linha (indice da propria celula) e'
        ! descartado
        ! -------------------------------------------------------------
        block
            integer :: atv_ijk_raw(4, ncell)
            read(iunit) atv_ijk_raw
            allocate(this%atv_ijk(3, 0:this%TCpyCell))
            do Rn = 0, this%TCpyCell
                this%atv_ijk(:, Rn) = atv_ijk_raw(2:4, Rn + 1)
            end do
        end block

        ! -------------------------------------------------------------
        ! Total_NumOrbs[1..atomnum], FNAN[1..atomnum]
        ! -------------------------------------------------------------
        allocate(this%TotalNumOrbs(this%atomnum))
        read(iunit) this%TotalNumOrbs

        allocate(this%FNAN(this%atomnum))
        read(iunit) this%FNAN

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
            read(iunit) this%natn_flat(off + 1 : off + this%FNAN(ct_AN) + 1)
        end do

        allocate(this%ncn_flat(this%n_pairs))
        do ct_AN = 1, this%atomnum
            off = this%neighbor_offset(ct_AN)
            read(iunit) this%ncn_flat(off + 1 : off + this%FNAN(ct_AN) + 1)
        end do

        ! -------------------------------------------------------------
        ! tv[1..3][4], rtv[1..3][4]. OpenMX stores the direct lattice
        ! vectors as rows in Bohr; the first value of each row is unused.
        ! -------------------------------------------------------------
        block
            real(dp) :: tv_raw(4, 3), rtv_junk(4, 3)
            read(iunit) tv_raw
            read(iunit) rtv_junk
            do Rn = 1, 3
                this%lattice_bohr(Rn, :) = tv_raw(2:4, Rn)
            end do
        end block

        ! -------------------------------------------------------------
        ! Gxyz[1..atomnum][0..3]: Cartesian coordinates in Bohr. The first
        ! value of each atom record is unused.
        ! -------------------------------------------------------------
        block
            real(dp) :: gxyz_raw(4, this%atomnum)
            read(iunit) gxyz_raw
            allocate(this%cart_bohr(3, this%atomnum))
            this%cart_bohr = gxyz_raw(2:4, :)
        end block

        ! -------------------------------------------------------------
        ! block_offset: offset (0-based, em escalares) de cada par
        ! (ct_AN,h_AN) dentro dos blocos Hks/iHks/OLP -- calculado ANTES
        ! de ler esses blocos, pois so' depende de natn/TotalNumOrbs
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
        ! Cada par (ct_AN,h_AN) e' TNO1*TNO2 doubles CONTIGUOS no stream
        ! (loop "for i in TNO1: read TNO2 doubles" = flatten row-major),
        ! e os pares sao lidos em sequencia -> uma unica leitura de
        ! n_scalar doubles cobre o canal de spin inteiro.
        ! -------------------------------------------------------------
        allocate(this%Hks_flat(this%n_scalar, 0:this%SpinP_switch))
        block
            integer :: spin
            do spin = 0, this%SpinP_switch
                read(iunit) this%Hks_flat(:, spin)
            end do
        end block

        ! -------------------------------------------------------------
        ! iHks[spin][...], spin=0..2 -- so' existe se SpinP_switch==3
        ! -------------------------------------------------------------
        if (this%SpinP_switch == 3) then
            allocate(this%iHks_flat(this%n_scalar, 0:2))
            block
                integer :: spin
                do spin = 0, 2
                    read(iunit) this%iHks_flat(:, spin)
                end do
            end block
        end if

        ! -------------------------------------------------------------
        ! OLP[ct_AN][h_AN][TNO1][TNO2] -- um unico canal, mesma logica
        ! -------------------------------------------------------------
        allocate(this%OLP_flat(this%n_scalar))
        read(iunit) this%OLP_flat

        ! -------------------------------------------------------------
        ! Blocos que NAO precisamos materializar, mas que ficam no meio
        ! do caminho ate ChemP: OLPpo, OLPmo, DM, iDM. Todos seguem a
        ! MESMA estrutura atomo-vizinho (n_scalar elementos por "canal"),
        ! soh que repetida por direcao (x,y,z), por ordem (OLPpo), ou por
        ! canal de spin (DM/iDM). Pulamos lendo-e-descartando em blocos.
        ! -------------------------------------------------------------

        ! OLPpo: 3 direcoes * order_max ordens * n_scalar (vazio se order_max==0)
        call skip_doubles(iunit, 3_8 * int(this%order_max, 8) * int(this%n_scalar, 8))

        ! OLPmo: 3 direcoes * n_scalar -- SEMPRE presente, independente de order_max
        call skip_doubles(iunit, 3_8 * int(this%n_scalar, 8))

        ! DM: (SpinP_switch+1) * n_scalar
        call skip_doubles(iunit, int(this%SpinP_switch + 1, 8) * int(this%n_scalar, 8))

        ! iDM: 2 * n_scalar
        call skip_doubles(iunit, 2_8 * int(this%n_scalar, 8))

        ! -------------------------------------------------------------
        ! Solver
        ! -------------------------------------------------------------
        block
            integer :: tmp(1)
            read(iunit) tmp
            this%Solver = tmp(1)
        end block

        ! -------------------------------------------------------------
        ! ChemP, E_Temp, dipolos, Valence_Electrons, Total_SpinS
        ! (10 doubles, ordem exata de read_scfout.c: d_vec[0..9])
        ! -------------------------------------------------------------
        read(iunit) d_vec
        this%ChemP  = d_vec(1)
        this%E_Temp = d_vec(2)
        ! d_vec(3:5) = dipole_moment_core, d_vec(6:8) = dipole_moment_background,
        ! d_vec(9) = Valence_Electrons, d_vec(10) = Total_SpinS -- nao usados aqui.

        ! -------------------------------------------------------------
        ! Embedded OpenMX input: num_lines followed by fixed 256-byte
        ! records. It supplies species labels, which are not present in the
        ! binary geometry records themselves. Missing/malformed input does
        ! not prevent H/S conversion; has_atom_labels remains false.
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

        close(iunit)
    end subroutine ReadScfoutFile

    subroutine ParseEmbeddedAtomLabels(lines, nat, labels, found)
        character(len=*), intent(in) :: lines(:)
        integer, intent(in) :: nat
        character(len=8), intent(out) :: labels(nat)
        logical, intent(out) :: found

        character(len=256) :: line
        character(len=64) :: label
        logical, allocatable :: seen(:)
        logical :: in_atoms, closed
        integer :: i, idx, ios, nul

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

    pure function to_lower(text) result(lower)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: lower
        integer :: i, code

        lower = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
        end do
    end function to_lower

    ! Le-e-descarta n doubles do stream, em blocos, sem alocar um buffer
    ! do tamanho total (que pode ser grande: OLPmo/DM/iDM sao O(NOrbs^2)).
    ! Usa iostat para detectar EOF inesperado (calculo de offset errado)
    ! em vez de falhar silenciosamente lendo ChemP do lugar errado.
    subroutine skip_doubles(iunit, n)
        integer,            intent(in) :: iunit
        integer(kind=8),    intent(in) :: n

        real(dp) :: buf(SKIP_CHUNK)
        integer(kind=8) :: remaining
        integer :: this_chunk, ios

        if (n <= 0_8) return

        remaining = n
        do while (remaining > 0_8)
            this_chunk = int(min(remaining, int(SKIP_CHUNK, 8)))
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
