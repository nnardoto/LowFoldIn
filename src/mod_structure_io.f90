! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Geometria da estrutura: leitor de arquivo, fingerprint e auxiliares de rede.
!>
!> Papel no pipeline: o arquivo HSIO não guarda geometria, então toda
!> ferramenta espectral (atomselect, spectral, specmap) obtém rede, posições
!> atômicas e rótulos deste módulo. O formato do arquivo é o de
!> unfold::ReadStructureFile (&STRUCTURE, &SPECIES opcional, &ATOMS). O
!> fingerprint amarra um arquivo de regiões à estrutura da qual ele foi
!> construído (ver atomic_regions).
!>
!> Convenções: átomos começam em 1 e mantêm a ordem do arquivo (que deve ser a
!> ordem dos átomos no HSIO); comprimentos são cartesianos em Å; lattice(i,:) é
!> o i-ésimo vetor de rede (linhas, não colunas); as coordenadas fracionárias f
!> satisfazem r = sum_i f_i a_i. As rotinas da biblioteca nunca param o
!> programa: erros voltam por ierr/errmsg.
!>
!> Ver docs/design/contracts.md §1 e §1.1.
module structure_io
    use types, only: dp
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private

    ! ==================================================================
    ! TIPOS
    ! ==================================================================
    !> Comprimento máximo de um rótulo de átomo/espécie (rótulos maiores são truncados).
    integer, parameter, public :: LABEL_LEN = 8

    !> Estrutura cristalina lida do arquivo de estrutura. nat = 0 significa "vazia
    !> ou leitura falhou". species_* são alocados (tamanho 0 se não houver bloco
    !> &SPECIES) assim que a leitura passa da busca por &SPECIES/&ATOMS.
    type, public :: structure_t
        integer                               :: nat = 0
        real(dp)                              :: lattice(3,3) = 0.0_dp  ! lattice(i,:) = a_i, Å
        real(dp),                 allocatable :: cart(:,:)              ! (3,nat) Å
        character(len=LABEL_LEN), allocatable :: label(:)               ! (nat)
        logical                               :: has_species = .false.
        integer                               :: nspecies = 0
        character(len=LABEL_LEN), allocatable :: species_name(:)
        integer,                  allocatable :: species_norb(:)        ! por espécie, contado na base HSIO (incl. dobro de spin)
    end type structure_t

    ! ==================================================================
    ! API PÚBLICA / CONSTANTES
    ! ==================================================================

    public :: read_structure, structure_fingerprint, fnv1a64_hex
    public :: reciprocal_lattice, cart_to_frac, species_norb_of, write_structure

    real(dp),       parameter :: TWO_PI  = 2.0_dp * acos(-1.0_dp)
    real(dp),       parameter :: DET_TOL = 1.0e-12_dp   ! |det| (Å^3) abaixo do qual a rede é considerada singular
    integer(int64), parameter :: MASK32  = int(z'FFFFFFFF', int64)

contains

    ! ==================================================================
    ! LEITOR
    ! ==================================================================

    !> Lê um arquivo de estrutura em `st`.
    !>
    !> Layout: namelist &STRUCTURE (LatticeVectors, linhas = a1,a2,a3 em Å), depois
    !> um bloco &SPECIES opcional (contagem, depois linhas `label nOrb`), depois
    !> &ATOMS (contagem, depois linhas `label x y z`, Å cartesianos). Linhas em
    !> branco e linhas começando com '!' ou '#' são ignoradas dentro dos blocos.
    !>
    !> Nunca para o programa. ierr = 0 em sucesso; senão ierr /= 0, errmsg nomeia o
    !> arquivo e st%nat = 0. Códigos: 1 não abre, 2 namelist &STRUCTURE inválido,
    !> 3 rede ausente/singular, 4 &ATOMS não encontrado, 5 &SPECIES inválido,
    !> 6 contagem de átomos inválida, 7 linha de átomo inválida/ausente,
    !> 8 rótulo de átomo não está em &SPECIES.
    subroutine read_structure(filename, st, ierr, errmsg)
        character(len=*),  intent(in)  :: filename
        type(structure_t), intent(out) :: st
        integer,           intent(out) :: ierr
        character(len=*),  intent(out) :: errmsg

        integer            :: iu, ios, i, nsp, nat
        character(len=256) :: line, iom
        real(dp)           :: LatticeVectors(3,3), rcart(3)
        logical            :: opened
        namelist /STRUCTURE/ LatticeVectors

        ierr = 0
        errmsg = ''
        opened = .false.

        open(newunit=iu, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            call fail(1, "cannot open structure file: " // trim(filename))
            return
        end if
        opened = .true.

        ! 1. &STRUCTURE namelist
        LatticeVectors = 0.0_dp
        iom = ''
        read(iu, nml=STRUCTURE, iostat=ios, iomsg=iom)
        if (ios /= 0) then
            call fail(2, "failed to read &STRUCTURE in " // trim(filename) // ": " // trim(iom))
            return
        end if
        if (all(abs(LatticeVectors) < DET_TOL)) then
            call fail(3, "LatticeVectors missing or null in " // trim(filename))
            return
        end if
        if (abs(det3(LatticeVectors)) < DET_TOL) then
            call fail(3, "LatticeVectors singular (|det| < 1e-12) in " // trim(filename))
            return
        end if
        st%lattice = LatticeVectors

        ! 2. Avança até a primeira linha marcadora '&SPECIES' ou '&ATOMS'. Esses
        !    blocos são lidos à mão (não como namelist) porque suas linhas são
        !    texto livre.
        do
            read(iu, '(A)', iostat=ios) line
            if (ios /= 0) then
                call fail(4, "block &ATOMS not found in " // trim(filename))
                return
            end if
            if (trim(adjustl(line)) == "&SPECIES" .or. trim(adjustl(line)) == "&ATOMS") exit
        end do

        ! 2b. Bloco &SPECIES opcional; depois avança até &ATOMS
        if (trim(adjustl(line)) == "&SPECIES") then
            call next_data_line(iu, line, ios)
            nsp = 0
            if (ios == 0) read(line, *, iostat=ios) nsp
            if (ios /= 0 .or. nsp < 1) then
                call fail(5, "failed to read number of species in " // trim(filename))
                return
            end if
            st%has_species = .true.
            st%nspecies = nsp
            allocate(st%species_name(nsp), st%species_norb(nsp))
            do i = 1, nsp
                call next_data_line(iu, line, ios)
                if (ios == 0) read(line, *, iostat=ios) st%species_name(i), st%species_norb(i)
                if (ios /= 0) then
                    call fail(5, "failed to read species " // itoa(i) // " in " // trim(filename))
                    return
                end if
                if (st%species_norb(i) < 1) then
                    call fail(5, "invalid nOrb for species " // itoa(i) // " in " // trim(filename))
                    return
                end if
            end do
            do
                read(iu, '(A)', iostat=ios) line
                if (ios /= 0) then
                    call fail(4, "block &ATOMS not found in " // trim(filename))
                    return
                end if
                if (trim(adjustl(line)) == "&ATOMS") exit
            end do
        else
            allocate(st%species_name(0), st%species_norb(0))
        end if

        ! 3. &ATOMS: contagem, depois uma linha `label x y z` por átomo. Com uma
        !    tabela &SPECIES, todo rótulo deve constar nela.
        call next_data_line(iu, line, ios)
        nat = 0
        if (ios == 0) read(line, *, iostat=ios) nat
        if (ios /= 0 .or. nat < 1) then
            call fail(6, "failed to read number of atoms in " // trim(filename))
            return
        end if
        allocate(st%cart(3, nat), st%label(nat))
        do i = 1, nat
            call next_data_line(iu, line, ios)
            if (ios /= 0) then
                call fail(7, "fewer atom lines than nAtoms (" // itoa(nat) // ") in " // trim(filename))
                return
            end if
            read(line, *, iostat=ios) st%label(i), rcart
            if (ios /= 0) then
                call fail(7, "failed to read atom " // itoa(i) // " in " // trim(filename))
                return
            end if
            if (st%has_species) then
                if (.not. any(st%label(i) == st%species_name)) then
                    call fail(8, "atom " // itoa(i) // " (species '" // trim(adjustl(st%label(i))) // &
                        "') not listed in &SPECIES of " // trim(filename))
                    return
                end if
            end if
            st%cart(:, i) = rcart
        end do
        st%nat = nat
        close(iu)

    contains

        !> Registra o erro, invalida st e fecha o arquivo se estiver aberto.
        subroutine fail(code, msg)
            integer,          intent(in) :: code
            character(len=*), intent(in) :: msg
            ierr = code
            errmsg = msg
            st%nat = 0
            if (opened) close(iu)
        end subroutine fail

    end subroutine read_structure

    ! ==================================================================
    ! AUXILIARES DE REDE
    ! ==================================================================

    !> Rede recíproca (linhas): b(i,:) = b_i, com b_i . a_j = 2π δ_ij, em 1/Å e
    !> incluindo o fator 2π. Exige st%lattice não singular.
    !> b = 2π (A^-1)^T, onde A tem os a_i como linhas.
    function reciprocal_lattice(st) result(b)
        type(structure_t), intent(in) :: st
        real(dp)                      :: b(3,3)
        b = TWO_PI * transpose(inv3(st%lattice))
    end function reciprocal_lattice

    !> Coordenadas fracionárias f tais que r = Σ_i f_i a_i, para uma posição
    !> cartesiana r em Å. Sem reembrulhar em [0,1): imagens periódicas são mantidas.
    function cart_to_frac(st, r) result(f)
        type(structure_t), intent(in) :: st
        real(dp),          intent(in) :: r(3)
        real(dp)                      :: f(3)
        ! r = A^T f  =>  f = (A^T)^-1 r = (A^-1)^T r
        f = matmul(r, inv3(st%lattice))
    end function cart_to_frac

    !> Número de orbitais do átomo iat (base 1) a partir da tabela de espécies,
    !> casada pelo rótulo. Devolve -1 se não há tabela, iat está fora da faixa ou
    !> o rótulo não é encontrado (os chamadores devem tratar -1 como
    !> "desconhecido", não como erro).
    function species_norb_of(st, iat) result(n)
        type(structure_t), intent(in) :: st
        integer,           intent(in) :: iat
        integer                       :: n, k
        n = -1
        if (.not. st%has_species) return
        if (iat < 1 .or. iat > st%nat) return
        do k = 1, st%nspecies
            if (st%species_name(k) == st%label(iat)) then
                n = st%species_norb(k)
                return
            end if
        end do
    end function species_norb_of

    ! ==================================================================
    ! FINGERPRINT (IMPRESSÃO DIGITAL)
    ! ==================================================================

    !> Fingerprint da estrutura (contracts §1.1): FNV-1a 64 de um texto canônico
    !> "LFSTRUCT1\n<nat>\n<label x y z>\n..." com coordenadas arredondadas a
    !> 1e-4 Å (inteiros, logo sem "-0"). Entram só rótulos e posições; a rede
    !> não entra. Devolve 16 caracteres hexadecimais minúsculos.
    function structure_fingerprint(st) result(fp)
        type(structure_t),           intent(in) :: st
        character(len=16)                       :: fp
        character(len=:),            allocatable :: text
        character(len=24)                       :: buf
        character(len=1),            parameter  :: NL = achar(10)
        integer                                 :: i, j

        write(buf, '(I0)') st%nat
        text = "LFSTRUCT1" // NL // trim(buf) // NL
        do i = 1, st%nat
            text = text // trim(adjustl(st%label(i)))
            do j = 1, 3
                write(buf, '(I0)') nint(st%cart(j, i) * 1.0e4_dp, int64)
                text = text // " " // trim(buf)
            end do
            text = text // NL
        end do
        fp = fnv1a64_hex(text)
    end function structure_fingerprint

    !> Hash FNV-1a de 64 bits de `text`, como 16 caracteres hexadecimais
    !> minúsculos. O estado de 64 bits é guardado em dois limbs de 32 bits (hi, lo)
    !> em int64 para evitar overflow com sinal (Fortran não tem inteiros sem
    !> sinal). A multiplicação pelo primo FNV é feita limb a limb, mantendo só os
    !> 64 bits baixos do produto.
    function fnv1a64_hex(text) result(hex)
        character(len=*), intent(in) :: text
        character(len=16)            :: hex
        ! offset basis 0xcbf29ce484222325 ; primo 0x100000001b3 = 256*2^32 + 435
        integer(int64), parameter    :: P_HI = 256_int64, P_LO = 435_int64
        integer(int64)               :: hi, lo, t, carry
        integer                      :: i

        hi = int(z'CBF29CE4', int64)
        lo = int(z'84222325', int64)
        do i = 1, len(text)
            lo = ieor(lo, int(iachar(text(i:i)), int64))
            ! (hi*2^32 + lo) * (P_HI*2^32 + P_LO) mod 2^64
            t = lo * P_LO                        ! < 2^41
            carry = shiftr(t, 32)
            hi = iand(hi * P_LO + lo * P_HI + carry, MASK32)
            lo = iand(t, MASK32)
        end do
        hex = hex8(hi) // hex8(lo)
    end function fnv1a64_hex

    ! ==================================================================
    ! AUXILIARES PRIVADAS
    ! ==================================================================

    !> Hexadecimal minúsculo dos 32 bits baixos de v (8 caracteres, mais significativo primeiro).
    function hex8(v) result(s)
        integer(int64),    intent(in) :: v
        character(len=8)              :: s
        character(len=16), parameter  :: DIGITS = '0123456789abcdef'
        integer                       :: k, d
        do k = 1, 8
            d = int(iand(shiftr(v, 4 * (8 - k)), 15_int64))
            s(k:k) = DIGITS(d + 1:d + 1)
        end do
    end function hex8

    !> Próxima linha não vazia e que não seja comentário ('!' ou '#').
    subroutine next_data_line(iu, line, ios)
        integer,          intent(in)  :: iu
        character(len=*), intent(out) :: line
        integer,          intent(out) :: ios
        character(len=len(line))      :: raw
        do
            read(iu, '(A)', iostat=ios) raw
            if (ios /= 0) return
            line = adjustl(raw)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '!' .or. line(1:1) == '#') cycle
            exit
        end do
    end subroutine next_data_line

    !> Determinante de uma matriz 3x3 (expansão em cofatores pela primeira linha).
    pure function det3(A) result(det)
        real(dp), intent(in) :: A(3,3)
        real(dp)             :: det
        det = A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
            - A(1,2) * (A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
            + A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    end function det3

    !> Inversa por cofatores; o chamador garante que não é singular (checado na leitura).
    pure function inv3(A) result(Ainv)
        real(dp), intent(in) :: A(3,3)
        real(dp)             :: Ainv(3,3), det
        det = det3(A)
        Ainv(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2)) / det
        Ainv(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2)) / det
        Ainv(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2)) / det
        Ainv(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1)) / det
        Ainv(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1)) / det
        Ainv(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1)) / det
        Ainv(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1)) / det
        Ainv(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1)) / det
        Ainv(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1)) / det
    end function inv3

    !> Grava a estrutura no formato lido por read_structure (Å, ES24.16, o
    !> mesmo do openmx2hsio --write-structure). Sem stop.
    subroutine write_structure(filename, st, ierr, errmsg)
        character(len=*),  intent(in)  :: filename
        type(structure_t), intent(in)  :: st
        integer,           intent(out) :: ierr
        character(len=*),  intent(out) :: errmsg
        integer :: u, ios, i
        ierr = 0
        errmsg = ''
        open(newunit=u, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = 'cannot write structure file: ' // trim(filename); return
        end if
        write(u, '(A)') '&STRUCTURE'
        do i = 1, 3
            write(u, '(A,I0,A,3(1X,ES24.16))') '    LatticeVectors(', i, ',:) =', st%lattice(i, :)
        end do
        write(u, '(A)') '/'
        write(u, '(A)') ''
        if (st%has_species) then
            write(u, '(A)') '&SPECIES'
            write(u, '(I0)') st%nspecies
            do i = 1, st%nspecies
                write(u, '(A,1X,I0)') trim(st%species_name(i)), st%species_norb(i)
            end do
            write(u, '(A)') '/'
            write(u, '(A)') ''
        end if
        write(u, '(A)') '&ATOMS'
        write(u, '(I0)') st%nat
        do i = 1, st%nat
            write(u, '(A,3(1X,ES24.16))') trim(st%label(i)), st%cart(:, i)
        end do
        write(u, '(A)') '/'
        close(u)
    end subroutine write_structure

    !> Inteiro para string sem brancos (para mensagens de erro).
    pure function itoa(i) result(s)
        integer,          intent(in)  :: i
        character(len=:), allocatable :: s
        character(len=16)             :: buf
        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

end module structure_io
