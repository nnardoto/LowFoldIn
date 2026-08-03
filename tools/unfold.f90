! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! tools/unfold.f90
!
! Executavel independente para unfolding de estrutura de bandas de bases LCAO
! nao-ortogonais, conforme:
!
!   J. Quan, N. Rybin, M. Scheffler, C. Carbogno,
!   "Efficient band structure unfolding with atom-centered orbitals: General
!    theory and application", Phys. Rev. B 113, 085112 (2026). DOI: 10.1103/7xym-7388
!
! Uso:
!   unfold <arquivo_de_entrada>
!
! O arquivo de entrada e' formado por blocos, nesta ordem: &UNFOLD, depois
! (&ATOMMAP  OU  nada, se o modo automatico via arquivos de estrutura for
! usado), e por fim &PCPATH.
!
!   &UNFOLD
!       ProjectName = "MeuSistema"
!       HDataFile   = "SC_H.hsio"      ! H(R) da SUPERCELULA (mesmo formato usado em "bands")
!       SDataFile   = "SC_S.hsio"      ! S(R) da SUPERCELULA
!       OutputDir   = "./"
!
!       ! Matriz de transformacao PC->SC completa (Eq. 1), 3x3: escreva LINHA
!       ! por LINHA (forma visual, uma atribuicao de secao por linha da
!       ! matriz). O MOTOR atual so' aceita M DIAGONAL (SC "ortogonal"); M
!       ! nao-diagonal e' rejeitado com um erro claro (suporte planejado
!       ! para o futuro).
!       M(1,:)      = 2, 0, 0
!       M(2,:)      = 0, 2, 0
!       M(3,:)      = 0, 0, 2
!
!       Broadening  = 0.010            ! eV, largura Lorentziana da funcao espectral A(k,E)
!       Emin        = -10.0
!       Emax        =  10.0
!       nE          = 400
!       WriteSpectralFunction = .true.
!
!       ! --- modo automatico de mapeamento atomo->PC (RECOMENDADO) ---------
!       ! Fornecendo os dois arquivos de ESTRUTURA abaixo, o bloco &ATOMMAP
!       ! NAO e' lido: o mapeamento e' construido por casamento geometrico
!       ! entre as posicoes da PC e da SC (via M) e o simbolo quimico.
!       PCStructureFile = "pc.nml"
!       SCStructureFile = "sc.nml"
!       MatchTolerance  = 1.0e-3       ! coord. fracionaria da PC (opcional)
!   /
!
!   Formato de PCStructureFile / SCStructureFile -- MESMO ESTILO .nml do
!   resto do pacote, um arquivo para a PC e outro para a SC (SCStructureFile
!   com atomos na MESMA ORDEM do arquivo .hsio). Ate' tres blocos, nesta
!   ordem: vetores de rede, tabela especie->numero de orbitais (OPCIONAL --
!   se presente, CONFERE contra o .hsio e pega cedo erros de base/ordenacao;
!   se ausente, o casamento geometrico continua funcionando via especie
!   quimica, so' sem essa checagem extra), e a lista de atomos da celula:
!
!       &STRUCTURE
!           LatticeVectors(1,:) = a1x, a1y, a1z   ! Angstrom -- MESMA convencao
!           LatticeVectors(2,:) = a2x, a2y, a2z   ! de POSCAR/ASE/CIF: cada
!           LatticeVectors(3,:) = a3x, a3y, a3z   ! LINHA e' um vetor de rede.
!       /
!
!       &SPECIES              ! OPCIONAL
!       <nEspecies>
!       Especie  nOrbitais
!       ...
!       /
!
!       &ATOMS
!       <nAtomos>
!       Especie  x  y  z    ! Cartesiano, Angstrom
!       ...
!       /
!
!   IMPORTANTE:
!     - Namelist do Fortran usa PARENTESES para indexar array, "M(1,:)" ou
!       "LatticeVectors(1,:)" -- NAO colchetes ("M[1,:]" da' erro de leitura).
!     - Para celulas ORTOGONAIS (todos os angulos 90 graus) e' tentador
!       escrever tudo numa linha so' (ex. "LatticeVectors = 10,0,0, 0,10,0,
!       0,0,10)"; isso ainda funciona porque nesse caso especial a matriz e'
!       igual a sua transposta. Para uma celula QUALQUER (hexagonal,
!       monoclinica, etc.), use SEMPRE a forma explicita por linha acima --
!       o preenchimento column-major do Fortran para uma lista "achatada"
!       transporia a matriz silenciosamente.
!
!   Exemplo (celula cubica de 10 A com 2 atomos, Ga na origem com 9 orbitais
!   de base e N em (0.25,0.25,0.25) fracionario com 4 orbitais):
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
!   --- OU, se PCStructureFile/SCStructureFile NAO forem fornecidos em
!       &UNFOLD, o mapeamento deve vir explicitamente aqui: -----------------
!
!   &ATOMMAP
!   ! uma linha por atomo da SC, NA MESMA ORDEM do arquivo .hsio (Orbs(:)):
!   !   PCAtom  n1 n2 n3
!   ! PCAtom  : indice do atomo da celula primitiva (1..nAtomsPC) do qual este
!   !           atomo e' replica (todas as m1*m2*m3 replicas de um mesmo atomo-PC
!   !           devem ter o MESMO numero de orbitais)
!   ! n1 n2 n3: indice inteiro 0-based da celula-PC que este atomo ocupa dentro
!   !           da SC ao longo de a1,a2,a3 (0 <= ni < Mi)
!       1   0 0 0
!       2   0 0 0
!       1   1 0 0
!       2   1 0 0
!       ...
!   /
!
!   &PCPATH
!   ! caminho de pontos k da CELULA PRIMITIVA (coordenadas fracionarias na
!   ! base reciproca da PC), EXATAMENTE no mesmo formato do bloco &BANDS
!       <npaths>
!       fx fy fz  nkpts
!       ...
!   /
!
! Saidas (em OutputDir/ProjectName.*):
!   .unfold.dat           -> (x, E, W) formato "scatter" para plot com peso como cor
!   .unfold_spectral.dat  -> A(k,E) numa grade de energias (se WriteSpectralFunction)
!
! Paralelizacao: o laco sobre os pontos k e' paralelizado com OpenMP (mesma
! tecnica de calc_band_structure_non_orth em mod_core.f90, usada por
! "bands"), respeitando OMP_NUM_THREADS. Progresso e' impresso apos cada ponto concluido, com o tempo individual.
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
    integer :: iu, ios, nargs

    type(UnfoldConfig)      :: cfg
    type(HSDataFormat)      :: HData, SData
    type(UnfoldResultFormat) :: Result
    type(UnfoldPerformance)  :: Performance

    character(len=:), allocatable :: BaseName, ScatterFile, SpectralFile
    integer :: ik
    integer :: nOrbitsPC_expected
    real(dp) :: worst_dev
    integer(int64) :: total_start, preparation_end, writing_start, writing_end, total_end, clock_rate
    real(dp) :: preparation_time, writing_time, total_time

    call system_clock(count_rate=clock_rate)
    call system_clock(count=total_start)

    nargs = command_argument_count()
    if (nargs < 1) then
        print '(A)', "Uso: unfold <arquivo_de_entrada>"
        stop 1
    end if
    call get_command_argument(1, input_path)

    ! ------------------------------------------------------------------
    ! 1. Le o bloco &UNFOLD (nomes de arquivo, tamanho da SC, saida)
    ! ------------------------------------------------------------------
    open(newunit=iu, file=trim(input_path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print '(A, A)', "Erro: nao foi possivel abrir o arquivo de entrada: ", trim(input_path)
        stop 1
    end if

    call ReadUnfoldNamelist(iu, cfg)

    ! ------------------------------------------------------------------
    ! 2. Carrega H(R) e S(R) da supercelula
    ! ------------------------------------------------------------------
    call print_box("Carregando dados da Supercelula (SC)")
    call HData%ReadFile(trim(cfg%HDataFile))
    call SData%ReadFile(trim(cfg%SDataFile))
    call print_box_close()

    ! ------------------------------------------------------------------
    ! 3. Le o mapeamento atomo(SC) -> (atomo-PC, celula) e o caminho k da PC
    ! ------------------------------------------------------------------
    ! ------------------------------------------------------------------
    ! 3. Determina o mapeamento atomo(SC) -> (atomo-PC, celula):
    !    - se &UNFOLD forneceu PCStructureFile/SCStructureFile, constroi
    !      automaticamente por casamento geometrico (RECOMENDADO);
    !    - caso contrario, le o bloco manual &ATOMMAP do proprio arquivo
    !      de entrada (modo antigo, util para casos pequenos/depuracao).
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

    ! Diagnostico: soma_N W^k_KN deveria ser p, o numero de orbitais-PC por
    ! celula (= NOrbsSC / (m1*m2*m3)) -- NAO 1, exceto no caso especial p=1
    ! [Eq. (28)+(51) do artigo: a soma em `orbit` na definicao do peso ja'
    ! projeta no subespaco degenerado inteiro de dimensao p, entao somar
    ! ainda sobre N acumula um fator p adicional por completude].
    nOrbitsPC_expected = HData%NOrbs / (cfg%MDiag(1) * cfg%MDiag(2) * cfg%MDiag(3))
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
        real(dp), intent(in) :: prep_time, write_time, elapsed_total

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
