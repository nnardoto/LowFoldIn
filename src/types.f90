! ----------------------------------------------------------------------------------------------------------
! Descrição Completa do Formato Esperado do Arquivo:
! O arquivo deve conter os seguintes elementos na ordem:
! 1. Título                       (uma linha de texto - até 256 caracteres)
! 2. Número de vizinhos           (um inteiro)
! 3. Número de átomos por célula  (um inteiro)
! 4. Número de total de Orbitais  (um inteiro - Para teste de Sanidade)
! 5. Lista de Degenerescência     (uma matriz de inteiros com dimensão [Neighbors] - Wannier Orbital Degeneracy)
! 6. Lista de Orbitais por Atomo  (uma matriz de inteiros com dimensão [nAtoms])
! 7. Lista de índices             (uma matriz de inteiros com dimensões [3, Neighbors]): hkl de cada vizinho
! 8. Dados                        (uma matriz de números complexos com dimensões [NOrbs, NOrbs, Neighbors])
! ----------------------------------------------------------------------------------------------------------
! 
! Nathanael N. Batista
! Agosto de 2026

module types
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    type :: FileFormat
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
    end type FileFormat

contains

! Subrotina para Leitura de Arquivo no Formato Especificado 
subroutine ReadFile(this, filename)
    class(FileFormat), intent(inout)        :: this
    character(len=*),  intent(in)           :: Filename
    integer                                 :: ierr, File

    open(newunit=File, file=Filename, status='old', action='read', iostat=ierr, form='unformatted', access='stream')
        read(File, iostat=ierr)     this%Title
        read(File, iostat=ierr)     this%Neighbors
        read(File, iostat=ierr)     this%nAtoms
        read(File, iostat=ierr)     this%NOrbs
        
        allocate(this%Degenerecency(this%Neighbors))
        allocate(this%Orbs(this%nAtoms))
        allocate(this%IndexList(3, this%Neighbors))
        allocate(this%Data(this%NOrbs, this%NOrbs, this%Neighbors))

        read(File, iostat=ierr)     this%Degenerecency
        read(File, iostat=ierr)     this%Orbs
        read(File, iostat=ierr)     this%IndexList
        read(File, iostat=ierr)     this%Data
    close(File) 
end subroutine ReadFile

! Subrotina para Escrita de Arquivo no Formato Especificado
subroutine WriteFile(this, filename)
    class(FileFormat), intent(in)       :: this
    character(len=*), intent(in)        :: Filename
    integer                             :: ierr, File

    open(newunit=File, file=Filename, status='replace', action='write', iostat=ierr, form='unformatted', access='stream')
        write(File, iostat=ierr)     this%Title
        write(File, iostat=ierr)     this%Neighbors
        write(File, iostat=ierr)     this%nAtoms
        write(File, iostat=ierr)     this%NOrbs
        write(File, iostat=ierr)     this%Degenerecency
        write(File, iostat=ierr)     this%Orbs
        write(File, iostat=ierr)     this%IndexList
        write(File, iostat=ierr)     this%Data
    close(File) 
end subroutine WriteFile

subroutine Clear(this)
    class(FileFormat), intent(inout)   :: this
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
    class(FileFormat), intent(in) :: this
    integer :: i, j

    ! Cabeçalho Principal
    print '(A)', "================================================================="
    print '(A)', "                      INFORMAÇÕES DO ARQUIVO                     "
    print '(A)', "================================================================="
    
    ! Informações Básicas (Alinhadas em colunas)
    print '(A, A)',       "  Título:         ", trim(adjustl(this%Title))
    print '(A, I8)',      "  Nº de Vizinhos: ", this%Neighbors
    print '(A, I8)',      "  Átomos / Célula:", this%nAtoms
    print '(A, I8)',      "  Total Orbitais: ", this%NOrbs
    print '(A)', "-----------------------------------------------------------------"

    ! Lista de Orbitais por Átomo
    if (allocated(this%Orbs)) then
        print '(A)', "  Orbitais por Átomo:"
        write(*, '("    ", 10(I0, 1X))') this%Orbs
    end if
    
    print '(A)', "-----------------------------------------------------------------"

    ! Lista de Índices hkl (Formatado como Matriz N x 3)
    if (allocated(this%IndexList)) then
        print '(A)', "  Lista de Índices (h, k, l):"
        ! Mostra os primeiros 5 para não lotar a tela se houver muitos vizinhos
        do i = 1, min(this%Neighbors, 5)
            write(*, '("    Vizinho ", I3, ": [", I4, ",", I4, ",", I4, " ]")') &
                i, this%IndexList(1, i), this%IndexList(2, i), this%IndexList(3, i)
        end do
        if (this%Neighbors > 5) then
            print '(A, I0, A)', "    ... e mais ", this%Neighbors - 5, " vizinhos."
        end if
    end if

    print '(A)', "-----------------------------------------------------------------"

    ! Amostra da Matriz de Dados (Complexos)
    if (allocated(this%Data)) then
        print '(A)', "  Dados Complexos Data(1, 1:min(NOrbs,3), 1:min(NOrbs,3)):"
        print '(A)', "  (Amostra do 1º Vizinho):"
        
        do i = 1, min(this%NOrbs, 3)
            write(*, '("    ")', advance='no')
            do j = 1, min(this%NOrbs, 3)
                ! Formata números complexos como: ( 1.23 + 4.56i )
                if (aimag(this%Data(i, j, 1)) >= 0.0_8) then
                    write(*, '("(", F7.3, "+", F7.3, "i) ")', advance='no') &
                        real(this%Data(i, j, 1)), aimag(this%Data(i, j, 1))
                else
                    write(*, '("(", F7.3, "-", F7.3, "i) ")', advance='no') &
                        real(this%Data(i, j, 1)), abs(aimag(this%Data(i, j, 1)))
                end if
            end do
            print * ! Nova linha
        end do
        
        if (this%NOrbs > 3) then
            print '(A)', "    ..."
        end if
    end if

    print '(A)', "================================================================="
end subroutine DisplayInfo

end module types
