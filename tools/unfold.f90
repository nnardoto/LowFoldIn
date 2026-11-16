! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! tools/unfold.f90
!
! Executável independente para unfolding de estrutura de bandas de bases LCAO
! não-ortogonais, conforme:
!
!   J. Quan, N. Rybin, M. Scheffler, C. Carbogno,
!   "Efficient band structure unfolding with atom-centered orbitals: General
!    theory and application", Phys. Rev. B 113, 085112 (2026). DOI: 10.1103/7xym-7388
!
! Uso:
!   unfold <arquivo_de_entrada>
!
! O arquivo de entrada é formado por blocos, nesta ordem: &UNFOLD, depois
! (&ATOMMAP  OU  nada, se o modo automático via arquivos de estrutura for
! usado), e por fim &PCPATH.
!
!   &UNFOLD
!       ProjectName = "MeuSistema"
!       HDataFile   = "SC_H.hsio"      ! H(R) da SUPERCÉLULA (mesmo formato usado em "bands")
!       SDataFile   = "SC_S.hsio"      ! S(R) da SUPERCÉLULA
!       OutputDir   = "./"
!
!       ! Matriz de transformação PC->SC completa (Eq. 1), 3x3: escreva LINHA
!       ! por LINHA (forma visual, uma atribuição de seção por linha da
!       ! matriz). M pode ser diagonal (rotina rápida especializada) ou
!       ! NÃO-diagonal (motor geral, mais caro; det(M) /= 0).
!       M(1,:)      = 2, 0, 0
!       M(2,:)      = 0, 2, 0
!       M(3,:)      = 0, 0, 2
!
!       Broadening  = 0.010            ! eV, FWHM Lorentziana da função espectral A(k,E)
!       Emin        = -10.0
!       Emax        =  10.0
!       nE          = 400
!       WriteSpectralFunction = .true.
!
!       ! --- modo automático de mapeamento átomo->PC (RECOMENDADO) ---------
!       ! Fornecendo os dois arquivos de ESTRUTURA abaixo, o bloco &ATOMMAP
!       ! NÃO é lido: o mapeamento é construído por casamento geométrico
!       ! entre as posições da PC e da SC (via M) e o símbolo químico.
!       PCStructureFile = "pc.nml"
!       SCStructureFile = "sc.nml"
!       MatchTolerance  = 1.0e-3       ! coord. fracionária da PC (opcional)
!   /
!
!   Formato de PCStructureFile / SCStructureFile -- MESMO ESTILO .nml do
!   resto do pacote, um arquivo para a PC e outro para a SC (SCStructureFile
!   com átomos na MESMA ORDEM do arquivo .hsio). Até três blocos, nesta
!   ordem: vetores de rede, tabela espécie->número de orbitais (OPCIONAL --
!   se presente, CONFERE contra o .hsio e pega cedo erros de base/ordenacao;
!   se ausente, o casamento geométrico continua funcionando via espécie
!   química, só sem essa checagem extra), e a lista de átomos da célula:
!
!       &STRUCTURE
!           LatticeVectors(1,:) = a1x, a1y, a1z   ! Angstrom -- MESMA convenção
!           LatticeVectors(2,:) = a2x, a2y, a2z   ! de POSCAR/ASE/CIF: cada
!           LatticeVectors(3,:) = a3x, a3y, a3z   ! LINHA é um vetor de rede.
!       /
!
!       &SPECIES              ! OPCIONAL
!       <nEspecies>
!       Espécie  nOrbitais
!       ...
!       /
!
!       &ATOMS
!       <nAtomos>
!       Espécie  x  y  z    ! Cartesiano, Angstrom
!       ...
!       /
!
!   IMPORTANTE:
!     - Namelist do Fortran usa PARÊNTESES para indexar array, "M(1,:)" ou
!       "LatticeVectors(1,:)" -- NÃO colchetes ("M[1,:]" dá erro de leitura).
!     - Para células ORTOGONAIS (todos os ângulos 90 graus) é tentador
!       escrever tudo numa linha só (ex. "LatticeVectors = 10,0,0, 0,10,0,
!       0,0,10)"; isso ainda funciona porque nesse caso especial a matriz é
!       igual a sua transposta. Para uma célula QUALQUER (hexagonal,
!       monoclínica, etc.), use SEMPRE a forma explícita por linha acima --
!       o preenchimento column-major do Fortran para uma lista "achatada"
!       transporia a matriz silenciosamente.
!
!   Exemplo (célula cúbica de 10 A com 2 átomos, Ga na origem com 9 orbitais
!   de base e N em (0.25,0.25,0.25) fracionário com 4 orbitais):
!
!       &STRUCTURE
!           LatticeVectors(1,:) = 10.0, 0.0, 0.0
!           LatticeVectors(2,:) =  0.0,10.0, 0.0
!           LatticeVectors(3,:) =  0.0, 0.0,10.0
!       /
!
!       &SPECIES
!       2
!       Ga   9
!       N    4
!       /
!
!       &ATOMS
!       2
!       Ga   0.00  0.00  0.00
!       N    2.50  2.50  2.50
!       /
!
!   --- OU, se PCStructureFile/SCStructureFile NÃO forem fornecidos em
!       &UNFOLD, o mapeamento deve vir explicitamente aqui: -----------------
!
!   &ATOMMAP
!   ! uma linha por átomo da SC, NA MESMA ORDEM do arquivo .hsio (Orbs(:)):
!   !   PCAtom  n1 n2 n3
!   ! PCAtom  : índice do átomo da célula primitiva (1..nAtomsPC) do qual este
!   !           átomo é réplica (todas as m1*m2*m3 réplicas de um mesmo átomo-PC
!   !           devem ter o MESMO número de orbitais)
!   ! n1 n2 n3: índice inteiro 0-based da célula-PC que este átomo ocupa dentro
!   !           da SC ao longo de a1,a2,a3 (0 <= ni < Mi)
!       1   0 0 0
!       2   0 0 0
!       1   1 0 0
!       2   1 0 0
!       ...
!   /
!
!   &PCPATH
!   ! caminho de pontos k da CÉLULA PRIMITIVA (coordenadas fracionárias na
!   ! base recíproca da PC), EXATAMENTE no mesmo formato do bloco &BANDS
!       <npaths>
!       fx fy fz  nkpts
!       ...
!   /
!
! Saídas (em OutputDir/ProjectName.*):
!   .unfold.dat           -> (x, E, W) formato "scatter" para plot com peso como cor
!   .unfold_spectral.dat  -> A(k,E) numa grade de energias (se WriteSpectralFunction)
!
! Paralelização: o laço sobre os pontos k é paralelizado com OpenMP (mesma
! técnica de calc_band_structure_non_orth em mod_core.f90, usada por
! "bands"), respeitando OMP_NUM_THREADS. Progresso é impresso após cada ponto concluído, com o tempo
! individual.
!
! Nathanael N. Batista, Agosto de 2026
! ----------------------------------------------------------------------------------------------------------
program unfold_tool
    use iso_fortran_env, only: int64
    use types,     only: dp
    use HS_io,     only: HSDataFormat
    use unfold,    only: UnfoldConfig, UnfoldResultFormat, &
                          ReadUnfoldNamelist, ReadAtomMap, ReadPCPath, DisplayUnfoldConfig, &
                          BuildAtomMapFromPositions, UsesAutomaticAtomMap, &
                          RunUnfolding, UnfoldPerformance
    use terminal_ui, only: print_box, print_box_close
    implicit none

    character(len=1024) :: input_path
    integer             :: iu, ios, nargs

    type(UnfoldConfig)       :: cfg
    type(HSDataFormat)       :: HData, SData
    type(UnfoldResultFormat) :: Result
    type(UnfoldPerformance)  :: Performance

    character(len=:), allocatable :: BaseName, ScatterFile, SpectralFile
    integer                       :: ik
    integer                       :: nOrbitsPC_expected
    real(dp)                      :: worst_dev
    integer(int64)                :: total_start, preparation_end, writing_start, writing_end, total_end, clock_rate
    real(dp)                      :: preparation_time, writing_time, total_time

    call system_clock(count_rate=clock_rate)
    call system_clock(count=total_start)

    nargs = command_argument_count()
    if (nargs < 1) then
        print '(A)', "Uso: unfold <arquivo_de_entrada>"
        stop 1
    end if
    call get_command_argument(1, input_path)

    ! ------------------------------------------------------------------
    ! 1. Lê o bloco &UNFOLD (nomes de arquivo, tamanho da SC, saída)
    ! ------------------------------------------------------------------
    open(newunit=iu, file=trim(input_path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print '(A, A)', "Erro: nao foi possivel abrir o arquivo de entrada: ", trim(input_path)
        stop 1
    end if

    call ReadUnfoldNamelist(iu, cfg)

    ! ------------------------------------------------------------------
    ! 2. Carrega H(R) e S(R) da supercélula
    ! ------------------------------------------------------------------
    call print_box("Carregando dados da Supercelula (SC)")
    call HData%ReadFile(trim(cfg%HDataFile))
    call SData%ReadFile(trim(cfg%SDataFile))
    call print_box_close()

    ! ------------------------------------------------------------------
    ! 3. Determina o mapeamento átomo(SC) -> (átomo-PC, célula):
    !    - se &UNFOLD forneceu PCStructureFile/SCStructureFile, constrói
    !      automaticamente por casamento geométrico (RECOMENDADO);
    !    - caso contrário, lê o bloco manual &ATOMMAP do próprio arquivo
    !      de entrada (modo antigo, útil para casos pequenos/depuração).
    !    Em qualquer um dos dois casos, o caminho k da PC vem de &PCPATH.
    ! ------------------------------------------------------------------
    if (UsesAutomaticAtomMap(cfg)) then
        call print_box("Construindo Mapa Atomo-PC (automatico, via posicoes)")
        call BuildAtomMapFromPositions(cfg, HData)
        call print_box_close()
    else
        call ReadAtomMap(iu, cfg, HData%nAtoms)
    end if
    call ReadPCPath(iu, cfg%PCPath)
    close(iu)

    call DisplayUnfoldConfig(cfg, HData%nAtoms, HData%NOrbs)

    ! ------------------------------------------------------------------
    ! 4. Executa o unfolding
    ! ------------------------------------------------------------------
    call print_box("Executando Unfolding")
    call system_clock(count=preparation_end)
    call RunUnfolding(cfg, HData, SData, Result, Performance)

    ! Diagnóstico: soma_N W^k_KN deveria ser p, o número de orbitais-PC por
    ! célula (= NOrbsSC / (m1*m2*m3)) -- NÃO 1, exceto no caso especial p=1
    ! [Eq. (28)+(51) do artigo: a soma em `orbit` na definição do peso já
    ! projeta no subespaço degenerado inteiro de dimensão p, então somar
    ! ainda sobre N acumula um fator p adicional por completude].
    nOrbitsPC_expected = HData%NOrbs / cfg%nCells
    worst_dev = 0.0_dp
    do ik = 1, cfg%PCPath%total_kpoints
        worst_dev = max(worst_dev, abs(Result%SumCheck(ik) - real(nOrbitsPC_expected, dp)))
    end do
    print '(A, I0, A, ES12.4)', "  Regra de soma esperada: sum_N W = p = ", nOrbitsPC_expected, &
                                 " ; maior desvio observado: ", worst_dev
    if (worst_dev > 1.0e-6_dp) then
        print '(A)', "  Aviso: desvio da regra de soma acima de 1e-6. Verifique &ATOMMAP e M."
    end if
    call print_box_close()

    ! ------------------------------------------------------------------
    ! 5. Exporta os resultados
    ! ------------------------------------------------------------------
    BaseName     = BuildOutputPath(cfg%OutputDir, cfg%ProjectName)
    ScatterFile  = BaseName // ".unfold.dat"
    SpectralFile = BaseName // ".unfold_spectral.dat"

    call print_box("Exportando Resultados")
    call system_clock(count=writing_start)
    call Result%SaveScatter(ScatterFile, cfg%PCPath)
    if (cfg%WriteSpectralFunction) then
        call Result%SaveSpectralFunction(SpectralFile, cfg%PCPath, cfg%Emin, cfg%Emax, cfg%nE, cfg%Broadening)
    end if
    call system_clock(count=writing_end)
    call print_box_close()

    call system_clock(count=total_end)
    preparation_time = real(preparation_end - total_start, dp) / real(clock_rate, dp)
    writing_time = real(writing_end - writing_start, dp) / real(clock_rate, dp)
    total_time = real(total_end - total_start, dp) / real(clock_rate, dp)
    call PrintPerformanceSummary(Performance, preparation_time, writing_time, total_time)

contains

    subroutine PrintPerformanceSummary(perf, prep_time, write_time, elapsed_total)
        type(UnfoldPerformance), intent(in) :: perf
        real(dp),                intent(in) :: prep_time, write_time, elapsed_total

        call print_box("Resumo de Desempenho")
        print '(A, I12)',      "  Pontos k calculados             : ", perf%PointsCalculated
        print '(A, I12)',      "  Threads OpenMP                  : ", perf%ThreadsUsed
        print '(A, F12.3, A)', "  Tempo de preparacao             : ", prep_time, " s"
        print '(A, F12.3, A)', "  Tempo de calculo                : ", perf%CalculationTime, " s"
        print '(A, F12.3, A)', "  Tempo de escrita                : ", write_time, " s"
        print '(A, F12.3, A)', "  Tempo total                     : ", elapsed_total, " s"
        print '(A)',           "  ----------------------------------------------------------------"
        print '(A, F12.3, A)', "  Tempo minimo por ponto          : ", perf%PointTimeMin, " s"
        print '(A, F12.3, A)', "  Tempo medio por ponto           : ", perf%PointTimeMean, " s"
        print '(A, F12.3, A)', "  Tempo maximo por ponto          : ", perf%PointTimeMax, " s"
        print '(A, F12.3, A)', "  Desvio-padrao                   : ", perf%PointTimeStdDev, " s"
        print '(A, F12.3, A)', "  Tempo efetivo por ponto         : ", perf%EffectiveTimePerPoint, " s"
        print '(A, F12.4, A)', "  Vazao                           : ", perf%ThroughputPointsPerSecond, " pontos/s"
        print '(A, F12.2, A)', "  Vazao                           : ", 60.0_dp * perf%ThroughputPointsPerSecond, &
                                  " pontos/min"
        print '(A, F12.1, A)', "  Eficiencia paralela aproximada : ", &
                                  100.0_dp * perf%ParallelEfficiency, " %"
        call print_box_close()
    end subroutine PrintPerformanceSummary

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

end program unfold_tool
