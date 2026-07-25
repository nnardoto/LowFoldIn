module toyModels
    use types
    implicit none
    private
    public :: ToyModel

    type :: ToyModel
        character(len=256)         :: Name          = ""
        integer                    :: NumParameters = 0
        real(8), allocatable       :: Parameters(:)
        type(FileFormat)           :: FileData

        contains
            procedure :: InitializeGraphene
            procedure :: DisplayInfo => DisplayToyModelInfo
    end type ToyModel

contains

    ! Subrotina para inicializar o modelo de Tight-Binding do Grafeno
    subroutine InitializeGraphene(this, e0, t)
        class(ToyModel), intent(inout) :: this
        real(8), intent(in), optional  :: e0, t

        this%Name = "Graphene_TB"
        this%NumParameters = 2
        
        if (allocated(this%Parameters)) deallocate(this%Parameters)
        allocate(this%Parameters(this%NumParameters))

        ! Valores padrão do Grafeno: e0 = 0.0 eV, t = -2.7 eV
        this%Parameters(1) = 0.0d0;  if (present(e0)) this%Parameters(1) = e0
        this%Parameters(2) = -2.7d0; if (present(t))  this%Parameters(2) = t

        ! Configurar o Cabeçalho do FileData
        this%FileData%Title     = "Graphene First-Neighbors Tight-Binding Model"
        this%FileData%Neighbors = 3  ! 3 Vizinhos mais próximos
        this%FileData%nAtoms    = 2  ! Átomos A e B por célula unitária
        this%FileData%NOrbs     = 2  ! Total: 1 orbital pz por átomo (2 no total)

        if (allocated(this%FileData%Orbs)) deallocate(this%FileData%Orbs)
        allocate(this%FileData%Orbs(this%FileData%nAtoms))
        this%FileData%Orbs = [1, 1]

        ! Configurar a Lista de Degenerescência
        if (allocated(this%FileData%Degenerecency)) deallocate(this%FileData%Degenerecency)
        allocate(this%FileData%Degenerecency(this%FileData%Neighbors))
        this%FileData%Degenerecency = [1, 1, 1]

        ! Configurar a Lista de Índices dos Vizinhos [3, Neighbors] (h, k, l)
        if (allocated(this%FileData%IndexList)) deallocate(this%FileData%IndexList)
        allocate(this%FileData%IndexList(3, this%FileData%Neighbors))

        this%FileData%IndexList(:, 1) = [ 0,  0, 0]
        this%FileData%IndexList(:, 2) = [-1,  0, 0]
        this%FileData%IndexList(:, 3) = [ 0, -1, 0]

        ! Matriz Hamiltoniana de Hopping Data(Neighbors, NOrbs, NOrbs)
        if (allocated(this%FileData%Data)) deallocate(this%FileData%Data)
        allocate(this%FileData%Data(this%FileData%Neighbors, &
                                    this%FileData%NOrbs,     &
                                    this%FileData%NOrbs))

        this%FileData%Data = cmplx(0.0d0, 0.0d0, kind=8)

        ! Vizinho 1: Célula [0,0,0]
        this%FileData%Data(1, 1, 1) = cmplx(this%Parameters(1), 0.0d0, kind=8) ! <A|H|A>
        this%FileData%Data(1, 2, 2) = cmplx(this%Parameters(1), 0.0d0, kind=8) ! <B|H|B>
        this%FileData%Data(1, 1, 2) = cmplx(this%Parameters(2), 0.0d0, kind=8) ! <A|H|B>
        this%FileData%Data(1, 2, 1) = cmplx(this%Parameters(2), 0.0d0, kind=8) ! <B|H|A>

        ! Vizinho 2: Célula [-1,0,0]
        this%FileData%Data(2, 1, 2) = cmplx(this%Parameters(2), 0.0d0, kind=8)

        ! Vizinho 3: Célula [0,-1,0]
        this%FileData%Data(3, 1, 2) = cmplx(this%Parameters(2), 0.0d0, kind=8)

    end subroutine InitializeGraphene

    ! Exibe informações do Modelo especificando 'class(ToyModel)'
    subroutine DisplayToyModelInfo(this)
        class(ToyModel), intent(in) :: this

        print '(A)', "================================================="
        print '(A, A)', "Modelo: ", trim(this%Name)
        if (allocated(this%Parameters) .and. size(this%Parameters) >= 2) then
            print '(A, F8.3, A, F8.3, A)', "Parâmetros -> On-site (e0): ", &
                  this%Parameters(1), " eV, Hopping (t): ", this%Parameters(2), " eV"
        end if
        print '(A)', "================================================="
        
        ! Chama a rotina de exibição contida no FileData
        call this%FileData%DisplayInfo()
    end subroutine DisplayToyModelInfo

end module toyModels