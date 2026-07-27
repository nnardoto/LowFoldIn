module config
    use types
    implicit none
    private

    public :: AppConfig, ReadConfigFromStdin, DisplayConfig

    ! npaths: Número de caminhos de k-points
    ! total_kpoints: Número total de k-points ao longo de todos os caminhos
    ! nkpts: Número de k-points em cada caminho
    ! kpoints: Lista de k-points para cada caminho (dimensão [3, npaths])
    ! FullPath: Lista completa de k-points ao longo de todos os caminhos
    ! xpoints: Posição acumulada ao longo do caminho de k-points (dimensão [total_kpoints])
    type :: BandsConfig
        integer                 :: npaths = 0
        integer, allocatable    :: total_kpoints
        integer,  allocatable   :: nkpts(:)
        real(dp), allocatable   :: kpoints(:, :)
        real(dp), allocatable   :: FullPath(:, :)
        real(dp), allocatable   :: xpoints(:)  ! Posição acumulada ao longo do caminho de k-points
    end type BandsConfig

    type :: AppConfig
        character(len = 256)              :: ProjectName = ""
        character(len = 256)              :: DataFile    = ""
        character(len = 256)              :: OutputDir   = "./"
        character(len = 256)              :: Calculation = ""

        ! Configuração específica para cálculo de bandas
        type(BandsConfig)                  :: Bands
    end type AppConfig

    contains

        subroutine ReadConfigFromStdin(config)
            use terminal_ui
            implicit none

            type(AppConfig), intent(inout)  :: config
            integer                         :: ierr

            character(len = 256)          :: ProjectName
            character(len = 256)          :: DataFile
            character(len = 256)          :: OutputDir
            character(len = 256)          :: Calculation
            character(len = 256)          :: ErrorMessage

            namelist /SYSTEM/ ProjectName, Calculation, DataFile, OutputDir

            ! 1. Preenche Valores Padrão do Módulo
            ProjectName = config%ProjectName
            Calculation = config%Calculation
            DataFile    = config%DataFile
            OutputDir   = config%OutputDir

            ! 2. Lê o Input Fornecido
            read(*, nml = system, iostat = ierr, iomsg = ErrorMessage)
            if (ierr /= 0) then
                call print_box("Input Error")
                call print_line("!")
                print*, "Erro ao ler o bloco SYSTEM"
                print*, "Mensagem de Erro: ", trim(ErrorMessage)
                call print_line("!")
                stop
            endif

            ! 3. Atualiza Valores a partir do Input
            config%ProjectName     = ProjectName   
            config%Calculation     = Calculation 
            config%DataFile        = DataFile    
            config%OutputDir       = OutputDir
            
            ! 4. Aqui deve vir um resumo do que foi lido, para o usuário ter certeza de que está correto
            call DisplayConfig(config)

            ! 5. Segue a Leitura a depender do Calculation
            select case (trim(Calculation))
            case ("bands")
                call ReadBandsConfig(config%Bands)
                call DisplayBandsConfig(config%Bands)
            case("unfold")
                stop "[Calculation Error] Unfolding não implementado"
            case default
                stop "[Calculation Error] Calculation Inexistente"
            end select
        end subroutine ReadConfigFromStdin

        subroutine DisplayConfig(config)
            use terminal_ui, only: print_box, print_row_split, print_box_close
            type(AppConfig), intent(in) :: config
            
            call print_box("System Configuration")

            ! Agora usando o alinhamento dividido nas pontas
            call print_row_split("Project Name:", config%ProjectName)
            call print_row_split("Calculation:",  config%Calculation)
            call print_row_split("Data File:",    config%DataFile)
            call print_row_split("Output Dir:",   config%OutputDir)

            call print_box_close()
        end subroutine DisplayConfig


        subroutine ReadBandsConfig(bands_cfg)
            type(BandsConfig), intent(inout) :: bands_cfg
            integer                          :: ierr, i
            character(len=256)               :: line
            integer                          :: idx, j
            real(dp)                         :: dk(3), xbuffer, length

            ! 1. Avança até  &BANDS
            do 
                read(*, '(A)', iostat=ierr) line
                if (ierr /= 0) then
                    stop "[BandsConfig Error] Erro ao ler o bloco &BANDS"
                endif
                if (trim(adjustl(line)) == "&BANDS") exit
            end do

            ! 2. Lê o número de caminhos
            read(*, *, iostat=ierr) bands_cfg%npaths
            if (ierr /= 0) then
                stop "[BandsConfig Error] Erro ao ler Número de Caminhos (npaths)"
            endif   

            ! 3. Aloca os arrays de k-points
            allocate(bands_cfg%nkpts(bands_cfg%npaths))
            allocate(bands_cfg%kpoints(3, bands_cfg%npaths))

            ! 4. Lê os k-points para cada caminho
            do i = 1, bands_cfg%npaths
                read(*, *, iostat=ierr) bands_cfg%kpoints(:, i), bands_cfg%nkpts(i)
                
                if (ierr /= 0) then
                    stop "[BandsConfig Error] Erro ao ler o caminho de k-points"
                endif
            end do

        ! 5. Gera a lista completa de k-points ao longo de todos os caminhos
            bands_cfg%total_kpoints = sum(bands_cfg%nkpts(1:bands_cfg%npaths - 1))  ! Exclui o último ponto do último caminho
            allocate(bands_cfg%FullPath(3, bands_cfg%total_kpoints))
            allocate(bands_cfg%xpoints(bands_cfg%total_kpoints))

            idx = 1
            xbuffer = 0.0_dp
            do i = 1 , bands_cfg%npaths - 1
                dk = (bands_cfg%kpoints(:, i + 1) - bands_cfg%kpoints(:, i)) / real(bands_cfg%nkpts(i), kind=dp)
                length = norm2(dk)
                do j = 1, bands_cfg%nkpts(i)
                    bands_cfg%FullPath(:, idx) = bands_cfg%kpoints(:, i) + real(j - 1, kind=dp) * dk
                    bands_cfg%xpoints(idx) = xbuffer  ! Posição acumulada ao longo do caminho
                    xbuffer = xbuffer + length
                    idx = idx + 1
                end do
            end do

            ! Últimos Pontos do caminho (não incluídos no loop acima)
            bands_cfg%xpoints(idx - 1)     = xbuffer                                
            bands_cfg%FullPath(:, idx - 1) = bands_cfg%kpoints(:, bands_cfg%npaths)  
        end subroutine ReadBandsConfig

        subroutine DisplayBandsConfig(bands_cfg)
            use terminal_ui, only: print_box, print_row_split, print_box_close, &
                                   UI_WIDTH, print_line
            type(BandsConfig), intent(in) :: bands_cfg
            integer                       :: i
            character(len=12)             :: npaths_str, ntotal_str

            ! 1. PAINEL DE METADADOS (Caixa Superior)
            call print_box("Bands Configuration")

            write(npaths_str, '(I0)') bands_cfg%npaths
            call print_row_split("Number of Paths:", npaths_str)
            
            write(ntotal_str, '(I0)') sum(bands_cfg%nkpts(1:bands_cfg%npaths-1))
            call print_row_split("Total K-Points:", ntotal_str)
            
            call print_box_close()
            
            print * ! Linha em branco de respiro visual entre as caixas

            ! 2. PAINEL DA TABELA (Caixa Inferior Separada)
            call print_box("Caminho de pontos K")

            ! Linha divisória interna dos títulos das colunas
            print '(A, A, A)', "|", repeat("-", UI_WIDTH - 2), "|"
            
            ! Títulos das colunas alinhados puramente por máscaras matemáticas (Sem espaços nas strings)
            print '(A, 2X, A4, 7X, A2, 12X, A2, 12X, A2, 14X, A5, 1X, A)', &
                "|", "Path", "Kx", "Ky", "Kz", "NKPTS", "|"
                
            print '(A, A, A)', "|", repeat("-", UI_WIDTH - 2), "|"

            ! Linhas de dados da tabela (Largura total rígida de 65 caracteres)
            do i = 1, bands_cfg%npaths
                print '(A, I6, " : ", 3F14.8, I11, A)', "|", &
                    i, bands_cfg%kpoints(:, i), bands_cfg%nkpts(i), " |"
            end do

            ! Fecha a base do painel da tabela
            call print_box_close()
            call print_box("List of all K-Points")
            
            do i = 1, sum(bands_cfg%nkpts(1:bands_cfg%npaths-1))
                print '(3F14.8, F14.8)', bands_cfg%FullPath(:, i), bands_cfg%xpoints(i)
            end do  

        end subroutine DisplayBandsConfig

end module config