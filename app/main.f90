program LowFoldIn
    use core
    use types
    use config
    use bands_io, only: BandsDataFormat
    use terminal_ui, only: print_box, print_box_close
    implicit none

    type(AppConfig)        :: cfg
    type(BandsDataFormat)  :: BandsData
    character(len=:), allocatable :: BaseName, TextFile, BinFile, BlockFile

    ! 1. Lê a configuração de entrada (via stdin) e o bloco &BANDS, se aplicável
    call ReadConfigFromStdin(cfg)

    ! 2. Segue o fluxo de acordo com o tipo de cálculo solicitado
    select case (trim(cfg%Calculation))
    case ("bands")

        ! 2.1 Calcula a estrutura de bandas (rota ortogonal ou não-ortogonal)
        if (cfg%nonOrtho) then
            BandsData%EigenValues = calc_band_structure_non_orth(cfg)
        else
            BandsData%EigenValues = calc_band_structure(cfg)
        end if

        ! 2.2 Monta os nomes dos arquivos de saída a partir do OutputDir/ProjectName
        BaseName  = BuildOutputPath(cfg%OutputDir, cfg%ProjectName)
        TextFile  = BaseName // ".bands.dat"
        BinFile   = BaseName // ".bands.bin"
        BlockFile = BaseName // ".bands.gnuplot.dat"

        ! 2.3 Exporta os resultados em texto (colunas), binário (NumPy) e em blocos (Gnuplot)
        call print_box("Exportando Resultados")
        call BandsData%SaveToText(TextFile, cfg%Bands)
        call BandsData%SaveToBinary(BinFile, cfg%Bands)
        call BandsData%SaveToGnuplotBlocks(BlockFile, cfg%Bands)
        call print_box_close()

    case default
        stop "[Calculation Error] Calculation Inexistente"
    end select

contains

    ! Monta o caminho base "OutputDir/ProjectName", cuidando da barra "/" faltante
    function BuildOutputPath(OutputDir, ProjectName) result(path)
        character(len=*), intent(in)  :: OutputDir
        character(len=*), intent(in)  :: ProjectName
        character(len=:), allocatable :: path
        character(len=:), allocatable :: dir

        dir = trim(adjustl(OutputDir))
        if (len(dir) == 0) dir = "."

        if (dir(len(dir):len(dir)) == "/") then
            path = dir // trim(adjustl(ProjectName))
        else
            path = dir // "/" // trim(adjustl(ProjectName))
        end if
    end function BuildOutputPath

end program LowFoldIn