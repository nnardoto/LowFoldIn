module bands_io
    use types
    implicit none

    private
    public :: BandsDataFormat

    type :: BandsDataFormat
        real(dp), allocatable       :: EigenValues(:, :)     ! Dimensão: [total_kpoints, nBands]
        complex(dp), allocatable    :: Eigenvectors(:, :, :)  ! Dimensão: [total_kpoints, nBands, NOrbs]
        logical                     :: HasEigenvectors = .false.

        contains
            procedure :: SaveToBinary
            procedure :: SaveToText
    end type BandsDataFormat

    contains

    !> Salva as Bandas em Formato de Texto Plano (ASCII)
    subroutine SaveToText(this, filename)
        class(BandsDataFormat), intent(in)  :: this
        character(len=*),       intent(in)  :: filename

        integer :: iunit, ios, nOrbs, nkPoints, ik
        character(len=256) :: FormatString

        if (.not. allocated(this%EigenValues)) then
            print *, "Error: Tentativa de salvar dados sem que os Eigenvalues estejam calculados."
            return
        end if

        nOrbs    = size(this%EigenValues, 1)
        nkPoints = size(this%EigenValues, 2)

        open(newunit=iunit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Não foi possível abrir o arquivo para escrita: ", filename
            return
        endif

        ! Cabeçalho legível por softwares de plot (linhas com # são ignoradas por padrão)
        write(iunit, '(A)') "# Dados de Estrutura de Bandas Gerados no LowFoldIn"
        write(iunit, '(A, I6, A, I6)') "# Total de Pontos K: ", nkPoints, " | Total de Orbitais: ", nOrbs
        write(iunit, '(A)') "# Coluna 1: Indice do Ponto K | Colunas 2+: Energias (eV)"

        ! Monta o formato dinamicamente para acomodar todas as colunas de energia de uma vez
        ! Formato: Um inteiro (I6) seguido por N colunas de ponto flutuante (F14.6)
        write(FormatString, '(A, I0, A)') '(I6, ', nOrbs, 'F14.6)'

        ! Escrita em formato de matriz/colunas
        do ik = 1, nkPoints
            ! Escreve: [Índice do ponto K] [Energia_Orb1] [Energia_Orb2] ... [Energia_OrbN]
            write(iunit, fmt=FormatString) ik, this%EigenValues(:, ik)
            
            ! LINHA EM BRANCO OPCIONAL: Se você usa o Gnuplot tradicional, descomente a linha abaixo 
            ! para separar os blocos de pontos k (útil para alguns tipos de loops de plot).
            ! write(iunit, *)
        end do

        close(iunit)
        print *, "Bandas salvos com sucesso em: ", filename
        print *, "Use o comando abaixo para plotar no Gnuplot:"
        print *, 'plot for [i=2:*] ', filename, ' using 1:i with lines notitle lw 1.5 lc rgb "blue"'
    end subroutine SaveToText

    !> Salva em formato binário bruto (Stream). Perfeito para ler em scripts Python com NumPy (Ultra rápido)
    subroutine SaveToBinary(this, filename)
        class(BandsDataFormat), intent(in) :: this
        character(len=*),       intent(in) :: filename
        integer :: iunit, ios

        open(newunit=iunit, file=filename, action='write', status='replace', &
             form='unformatted', access='stream', iostat=ios)
        if (ios /= 0) then
            print *, "Erro ao criar arquivo binário:", filename
            return
        end if

        ! Escreve os metadados brutos primeiro para o leitor saber o formato da matriz
        write(iunit) size(this%EigenValues, 1) ! nOrbs
        write(iunit) size(this%EigenValues, 2) ! nkPoints
        
        ! Salva o bloco inteiro de memória da matriz de uma só vez
        write(iunit) this%EigenValues

        close(iunit)
        print *, "Bandas salvas com sucesso em binário: ", filename
    end subroutine SaveToBinary

end module bands_io