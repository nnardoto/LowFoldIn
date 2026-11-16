! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! mod_unfold.f90
!
! Implementa o unfolding de estrutura de bandas para bases LCAO não-ortogonais,
! conforme:
!
!   J. Quan, N. Rybin, M. Scheffler, C. Carbogno,
!   "Efficient band structure unfolding with atom-centered orbitals: General
!    theory and application", Phys. Rev. B 113, 085112 (2026). DOI: 10.1103/7xym-7388
!
! Formalismo utilizado (ver Sec. II B-D do artigo):
!
!   1. Transformação de Lowdin: C'_KN = S(K)^(1/2) C_KN  [Eq. (24)]
!   2. Os autovetores F_kn do operador de translação da célula primitiva (PC),
!      representado na base da supercelula (SC), tem forma analítica fechada
!      quando a matriz de transformação M é diagonal (SC = M1 x M2 x M3
!      réplicas "ortogonais" da PC, o caso mais comum na prática). Combinando
!      a recursão das Eqs. (33)-(37) nas 3 direções (Eq. 48), obtém-se, para
!      cada "orbital primitivo" (uma identidade de orbital-PC que se repete
!      m1*m2*m3 vezes dentro da SC) o autovetor
!
!         F_{k}(célula L=(n1,n2,n3)) = (1/sqrt(m1*m2*m3)) * exp(i * 2*pi *
!                                        (fk1*n1 + fk2*n2 + fk3*n3))
!
!      onde fk = (fk1,fk2,fk3) são as coordenadas fracionárias do ponto k da
!      PC na base reciproca da PC, e (n1,n2,n3) é o índice inteiro (0-based)
!      da célula-PC replicada dentro da SC ao longo de a1,a2,a3.
!
!   3. O peso de "unfolding" de cada autoestado |KN> da SC no ponto k da PC
!      é então [Eq. (28)]:
!
!         W^k_KN = sum_{n=1}^{p} | F_{k,n}^dagger . C'_KN |^2
!
!      onde a soma em n percorre os p "orbitais primitivos" distintos da PC.
!
!   4. O ponto K da SC ao qual um dado k da PC é mapeado segue da Eq. (7).
!      Com a convenção deste código (vetores de rede como LINHAS,
!      A_i = sum_j M_ij a_j), as coordenadas fracionárias se transformam de
!      formas DIFERENTES no espaço real e no recíproco:
!         posições  (covariantes):     xi_PC = xi_SC . M   (vetor-linha)
!         pontos k  (contravariantes): F_K   = M . f_k     (vetor-coluna)
!      pois F_K,i = K.A_i/2pi = sum_j M_ij f_j. Para M diagonal ambas se
!      reduzem a FK_i = mod(fk_i * m_i, 1).
!
!   5. M diagonal usa a fórmula fechada das Eqs. (33)-(38); M não-diagonal usa
!      o projetor de grupo sobre Z^3/M.Z^3 (Eq. 14). Nesse grupo, os vetores
!      da SC em coordenadas inteiras da PC são as LINHAS de M, isto é, as
!      colunas de M^T: rotulo ~ rotulo + M^T.s.
!
! Nathanael N. Batista, Agosto de 2026
! ----------------------------------------------------------------------------------------------------------
module unfold
    use iso_fortran_env, only: int64, output_unit
    use types,  only: dp
    use HS_io,  only: HSDataFormat
    use core,   only: compute_Hk_at_kpoint
    implicit none
    private

    public :: PCPathConfig, UnfoldConfig, UnfoldResultFormat, UnfoldPerformance
    public :: ReadUnfoldNamelist, ReadAtomMap, ReadPCPath, DisplayUnfoldConfig
    public :: BuildAtomMapFromPositions, UsesAutomaticAtomMap, ReadStructureFile
    public :: RunUnfolding

    real(dp), parameter :: TWO_PI = 6.283185307179586476925286766559_dp

    ! ----------------------------------------------------------------------
    ! Caminho de pontos k da CÉLULA PRIMITIVA (coordenadas fracionárias na
    ! base reciproca da PC). Estrutura idêntica, em espírito, a BandsConfig
    ! de mod_config.f90, mas mantida separada para não acoplar os dois
    ! módulos (e para deixar claro que aqui os pontos são da PC, não da SC).
    ! ----------------------------------------------------------------------
    type :: UnfoldPerformance
        integer  :: PointsCalculated          = 0
        integer  :: ThreadsUsed               = 1
        real(dp) :: CalculationTime           = 0.0_dp
        real(dp) :: PointTimeMin              = 0.0_dp
        real(dp) :: PointTimeMean             = 0.0_dp
        real(dp) :: PointTimeMax              = 0.0_dp
        real(dp) :: PointTimeStdDev           = 0.0_dp
        real(dp) :: EffectiveTimePerPoint     = 0.0_dp
        real(dp) :: ThroughputPointsPerSecond = 0.0_dp
        real(dp) :: ParallelEfficiency        = 0.0_dp
    end type UnfoldPerformance

    type :: PCPathConfig
        integer               :: npaths = 0
        integer,  allocatable :: nkpts(:)
        real(dp), allocatable :: kpoints(:, :)   ! (3, npaths) vértices do caminho
        real(dp), allocatable :: FullPath(:, :)  ! (3, total_kpoints) caminho completo
        real(dp), allocatable :: xpoints(:)      ! (total_kpoints) posição acumulada
        integer               :: total_kpoints = 0
    end type PCPathConfig

    ! ----------------------------------------------------------------------
    ! Configuração completa do cálculo de unfolding
    ! ----------------------------------------------------------------------
    type :: UnfoldConfig
        character(len=256) :: ProjectName = ""
        character(len=256) :: HDataFile   = ""
        character(len=256) :: SDataFile   = ""
        character(len=256) :: OutputDir   = "./"

        ! Matriz de transformação PC->SC completa [Eq. (1) do artigo, A=a.M].
        ! M pode ser QUALQUER matriz 3x3 de inteiros não-singular (SC
        ! "rodada"/não-ortogonal incluída, Eq. 9 e Sec. II D do artigo).
        ! DiagonalM é calculado uma vez em ReadUnfoldNamelist e usado para
        ! escolher, em RunUnfolding, entre a rotina RÁPIDA especializada
        ! (M diagonal, O(nCells) por órbita-PC por ponto k -- fórmula
        ! fechada Eqs. 33-38) e a rotina GERAL (M qualquer, O(nCells^2) por
        ! órbita-PC por ponto k -- fórmula direta do operador de projeção,
        ! Eq. 14 + Eq. 10, provada equivalente a Eq. 28 no Apêndice A do
        ! artigo). MDiag(3) só é válido/usado quando DiagonalM = .true.
        integer  :: M(3, 3)   = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
        integer  :: MDiag(3)  = (/1, 1, 1/)
        logical  :: DiagonalM = .true.
        integer  :: nCells    = 1  ! = |det(M)|, número de células-PC na SC

        real(dp) :: Broadening            = 0.010_dp  ! eV (FWHM Lorentziana p/ A(k,E))
        real(dp) :: Emin                  = -10.0_dp
        real(dp) :: Emax                  =  10.0_dp
        integer  :: nE                    = 400
        logical  :: WriteSpectralFunction = .true.

        ! --------------------------------------------------------------
        ! Modo automático (RECOMENDADO): fornecendo os dois arquivos de
        ! ESTRUTURA abaixo (PC e SC, no formato .nml do próprio pacote --
        ! ver ReadStructureFile), o mapeamento átomo(SC)->(átomo-PC, célula)
        ! é construído automaticamente por casamento geométrico (ver
        ! BuildAtomMapFromPositions), dispensando o bloco &ATOMMAP manual.
        ! Cada arquivo de estrutura contém os vetores de rede (&STRUCTURE)
        ! e a lista de átomos da célula, Espécie X Y Z em Angstrom (&ATOMS).
        ! Os símbolos químicos são usados como checagem extra no casamento
        ! geométrico (um átomo de Ga nunca casa com um N).
        ! --------------------------------------------------------------
        character(len=256) :: PCStructureFile = ""
        character(len=256) :: SCStructureFile = ""
        real(dp)           :: MatchTolerance  = 1.0e-3_dp  ! coord. fracionária da PC

        ! Mapeamento átomo(SC) -> (átomo-PC, índice de célula-PC dentro da SC).
        ! Preenchido OU pelo bloco manual &ATOMMAP (ReadAtomMap), OU
        ! automaticamente a partir das posições (BuildAtomMapFromPositions).
        integer, allocatable :: PCOrbit(:)       ! (nAtomsSC) rótulo do átomo-PC (1..nAtomsPC)
        integer, allocatable :: CellIndex(:, :)  ! (3, nAtomsSC) índice (n1,n2,n3), 0-based

        type(PCPathConfig)   :: PCPath
    end type UnfoldConfig

    ! ----------------------------------------------------------------------
    ! Resultado do unfolding: para cada ponto k da PC, todas as energias e
    ! pesos W^k_KN dos estados da SC no K correspondente.
    ! ----------------------------------------------------------------------
    type :: UnfoldResultFormat
        real(dp), allocatable :: Energies(:, :)  ! (NOrbsSC, total_kpoints)
        real(dp), allocatable :: Weights(:, :)   ! (NOrbsSC, total_kpoints)
        real(dp), allocatable :: SumCheck(:)     ! (total_kpoints) soma_N W^k_KN (deveria ser p, o
                                                 ! número de orbitais-PC por célula; ver RunUnfolding)

        contains
            procedure :: SaveScatter
            procedure :: SaveSpectralFunction
    end type UnfoldResultFormat

contains

    ! ======================================================================
    ! LEITURA DE ENTRADA
    ! ======================================================================

    !> Lê a próxima linha "útil" da unidade `iu`, pulando linhas em branco e
    !> linhas de comentário puro (primeiro caractere não-branco "!" ou "#").
    !> Comentários ao FINAL de uma linha de dados (ex.: "0.5 0.0 0.0 10 !X")
    !> não precisam disso: list-directed READ já os ignora nativamente.
    subroutine NextDataLine(iu, line, ierr)
        integer,                  intent(in)  :: iu
        character(len=*),         intent(out) :: line
        integer,                  intent(out) :: ierr
        character(len=len(line))              :: raw

        do
            read(iu, '(A)', iostat=ierr) raw
            if (ierr /= 0) return
            line = adjustl(raw)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '!' .or. line(1:1) == '#') cycle
            exit
        end do
    end subroutine NextDataLine

    !> Lê o bloco de namelist &UNFOLD a partir da unidade já aberta `iu`.
    subroutine ReadUnfoldNamelist(iu, cfg)
        integer,            intent(in)    :: iu
        type(UnfoldConfig), intent(inout) :: cfg
        integer                           :: ierr
        character(len=256)                :: ErrorMessage

        character(len=256) :: ProjectName, HDataFile, SDataFile, OutputDir
        character(len=256) :: PCStructureFile, SCStructureFile
        integer            :: M(3, 3), nE
        real(dp)           :: Broadening, Emin, Emax, MatchTolerance
        logical            :: WriteSpectralFunction

        namelist /UNFOLD/ ProjectName, HDataFile, SDataFile, OutputDir, M,    &
                           Broadening, Emin, Emax, nE, WriteSpectralFunction, &
                           PCStructureFile, SCStructureFile, MatchTolerance

        ! Valores padrão
        ProjectName      = cfg%ProjectName
        HDataFile        = cfg%HDataFile
        SDataFile        = cfg%SDataFile
        OutputDir        = cfg%OutputDir
        M                = cfg%M
        Broadening       = cfg%Broadening
        Emin             = cfg%Emin
        Emax             = cfg%Emax
        nE               = cfg%nE
        WriteSpectralFunction = cfg%WriteSpectralFunction
        PCStructureFile  = cfg%PCStructureFile
        SCStructureFile  = cfg%SCStructureFile
        MatchTolerance   = cfg%MatchTolerance

        read(iu, nml=UNFOLD, iostat=ierr, iomsg=ErrorMessage)
        if (ierr /= 0) then
            print *, "[Unfold Input Error] Falha ao ler o bloco &UNFOLD"
            print *, "Mensagem: ", trim(ErrorMessage)
            stop 1
        end if

        cfg%ProjectName = ProjectName
        cfg%HDataFile   = HDataFile
        cfg%SDataFile   = SDataFile
        cfg%OutputDir   = OutputDir
        cfg%M           = M
        cfg%Broadening  = Broadening
        cfg%Emin        = Emin
        cfg%Emax        = Emax
        cfg%nE          = nE
        cfg%WriteSpectralFunction = WriteSpectralFunction
        cfg%PCStructureFile = PCStructureFile
        cfg%SCStructureFile = SCStructureFile
        cfg%MatchTolerance  = MatchTolerance

        ! Parâmetros numéricos de &UNFOLD (grade de energia, largura, tolerância)
        if (cfg%Broadening <= 0.0_dp) then
            print '(A)', "[Unfold Input Error] Broadening deve ser > 0"
            stop 1
        end if
        if (cfg%nE < 2) then
            print '(A)', "[Unfold Input Error] nE deve ser >= 2"
            stop 1
        end if
        if (cfg%Emax <= cfg%Emin) then
            print '(A)', "[Unfold Input Error] Emax deve ser > Emin"
            stop 1
        end if
        if (cfg%MatchTolerance <= 0.0_dp .or. cfg%MatchTolerance >= 0.5_dp) then
            print '(A)', "[Unfold Input Error] MatchTolerance deve estar em (0, 0.5)"
            stop 1
        end if

        call ClassifyM(cfg)

        if (cfg%DiagonalM .and. any(cfg%MDiag < 1)) then
            print *, "[Unfold Input Error] Os elementos da diagonal de M devem ser >= 1"
            stop 1
        end if
    end subroutine ReadUnfoldNamelist

    !> Classifica M (diagonal ou não) e calcula nCells = |det(M)| [Eq. (1):
    !> "the volume of the supercell is m = |det(M)| times larger than that
    !> of the PC"]. Quando M é diagonal, também preenche MDiag(3) para a
    !> rotina rápida especializada de RunUnfolding continuar funcionando
    !> exatamente como antes (mesmo código, mesmo desempenho).
    subroutine ClassifyM(cfg)
        type(UnfoldConfig), intent(inout) :: cfg
        integer                           :: i, j
        integer                           :: det

        cfg%DiagonalM = .true.
        do i = 1, 3
            do j = 1, 3
                if (i /= j .and. cfg%M(i, j) /= 0) cfg%DiagonalM = .false.
            end do
        end do

        det = cfg%M(1,1) * (cfg%M(2,2)*cfg%M(3,3) - cfg%M(2,3)*cfg%M(3,2)) &
            - cfg%M(1,2) * (cfg%M(2,1)*cfg%M(3,3) - cfg%M(2,3)*cfg%M(3,1)) &
            + cfg%M(1,3) * (cfg%M(2,1)*cfg%M(3,2) - cfg%M(2,2)*cfg%M(3,1))

        if (det == 0) then
            print '(A)', "[Unfold Input Error] M e' singular (det(M) = 0) -- nao define uma supercelula valida."
            stop 1
        end if
        cfg%nCells = abs(det)

        if (cfg%DiagonalM) then
            cfg%MDiag = (/ cfg%M(1,1), cfg%M(2,2), cfg%M(3,3) /)
        else
            cfg%MDiag = 0   ! não usado fora do caso diagonal; deixado obviamente inválido
            print '(A)', "  M nao e' diagonal: usando o motor de unfolding GERAL (Sec. II D do artigo)," // &
                          " mais caro (O(nCells) extra por ponto k) que a rotina especializada" //          &
                          " para M diagonal."
        end if
    end subroutine ClassifyM

    !> .true. se o usuário forneceu os dois arquivos de estrutura (PC e SC),
    !> modo automático de construção do &ATOMMAP; .false. se deve usar o
    !> bloco manual &ATOMMAP no próprio arquivo de entrada.
    function UsesAutomaticAtomMap(cfg) result(auto)
        type(UnfoldConfig), intent(in) :: cfg
        logical                        :: auto
        auto = (len_trim(cfg%PCStructureFile) > 0) .and. (len_trim(cfg%SCStructureFile) > 0)
    end function UsesAutomaticAtomMap

    !> Lê o bloco &ATOMMAP: exatamente nAtoms linhas "PCAtom  n1 n2 n3",
    !> na MESMA ordem dos átomos do arquivo .hsio (ordem de Orbs(:)).
    subroutine ReadAtomMap(iu, cfg, nAtoms)
        integer,            intent(in)    :: iu
        type(UnfoldConfig), intent(inout) :: cfg
        integer,            intent(in)    :: nAtoms

        character(len=256) :: line
        integer            :: ierr, i

        ! 1. Avança até a linha "&ATOMMAP"
        do
            read(iu, '(A)', iostat=ierr) line
            if (ierr /= 0) then
                print '(A)', "[Unfold Input Error] Bloco &ATOMMAP nao encontrado"
                stop 1
            end if
            if (trim(adjustl(line)) == "&ATOMMAP") exit
        end do

        allocate(cfg%PCOrbit(nAtoms))
        allocate(cfg%CellIndex(3, nAtoms))

        do i = 1, nAtoms
            call NextDataLine(iu, line, ierr)
            if (ierr == 0) read(line, *, iostat=ierr) cfg%PCOrbit(i), cfg%CellIndex(:, i)
            if (ierr /= 0) then
                print '(A, I0)', "[Unfold Input Error] Falha ao ler &ATOMMAP na linha do atomo ", i
                stop 1
            end if
        end do

        if (any(cfg%CellIndex < 0)) then
            print '(A)', "[Unfold Input Error] Indices de celula em &ATOMMAP devem ser >= 0 (0-based)"
            stop 1
        end if
        if (any(cfg%PCOrbit < 1)) then
            print '(A)', "[Unfold Input Error] Rotulos de atomo-PC em &ATOMMAP devem ser >= 1"
            stop 1
        end if
        if (cfg%DiagonalM) then
            do i = 1, 3
                if (any(cfg%CellIndex(i, :) >= cfg%MDiag(i))) then
                    print '(A, I0, A)', "[Unfold Input Error] Indice de celula na direcao ", i, &
                                         " excede M(i,i)-1"
                    stop 1
                end if
            end do
        end if
        ! Para M não-diagonal não há um limite por-eixo simples (os
        ! rótulos de célula válidos formam um paralelepípedo torto, não
        ! uma caixa); a consistência real é verificada de forma GERAL em
        ! RunUnfolding (cada orbital-PC precisa de exatamente nCells
        ! réplicas com rótulos de célula distintos).
    end subroutine ReadAtomMap

    !> Constroi automaticamente o mapeamento átomo(SC) -> (átomo-PC, célula)
    !> por CASAMENTO GEOMÉTRICO, lendo PCStructureFile e SCStructureFile
    !> (formato .nml próprio do pacote, ver ReadStructureFile).
    !>
    !> Para M diagonal, a coordenada fracionária de um átomo da SC na base da
    !> PC é obtida componente-a-componente por xi_PC(i) = M(i,i) * xi_SC(i)
    !> [caso diagonal da relação geral xi_PC = xi_SC . M, vetor-LINHA --
    !> consequência de A = M.a com LatticeVectors em linhas, convenção
    !> POSCAR/ASE/CIF desta base de código (ver ReadStructureFile); NÃO a
    !> notação vetor-coluna A=a.M do artigo, Eq. (1), que só coincide com
    !> esta quando M é diagonal (M=M^T)]. A parte inteira de xi_PC(i) é o
    !> índice de célula n_i (0-based), e a parte fracionária é comparada (com
    !> imagem mínima periódica, e mesmo símbolo químico) contra a lista de
    !> posições da PC para identificar de qual átomo-PC aquele átomo da SC
    !> é réplica.
    !>
    !> Além do casamento geométrico, esta rotina também CONFERE, para cada
    !> átomo, que o número de orbitais realmente presente no .hsio
    !> (HData%Orbs) bate com o número declarado no bloco &SPECIES do arquivo
    !> de estrutura correspondente (PC e SC) -- um jeito barato de pegar
    !> cedo erros de ordenação de átomos ou de conjunto de base usado.
    !>
    !> Esta rotina evita ter que digitar manualmente uma linha por átomo em
    !> &ATOMMAP -- basta fornecer a estrutura da PC e da SC (que o usuário
    !> já tem, vindas de qualquer ferramenta usada para construir a SC) e M.
    subroutine BuildAtomMapFromPositions(cfg, HData)
        type(UnfoldConfig), intent(inout) :: cfg
        type(HSDataFormat), intent(in)    :: HData

        real(dp),         allocatable :: PCPos(:, :)  ! (3, nAtomsPC)
        real(dp),         allocatable :: SCPos(:, :)  ! (3, nAtoms)
        character(len=8), allocatable :: PCEl(:), SCEl(:)
        character(len=8), allocatable :: PCSpecies(:), SCSpecies(:)
        integer,          allocatable :: PCSpeciesNOrbs(:), SCSpeciesNOrbs(:)
        integer                       :: nAtoms, nAtomsPC, i, j, jbest, expectedOrbs
        real(dp)                      :: xiSC(3), xiPC(3), remainder(3), diff(3), dist, distbest
        integer                       :: cellIndex(3)

        nAtoms = HData%nAtoms

        call ReadStructureFile(cfg%PCStructureFile, PCPos, PCEl, nAtomsPC, PCSpecies, PCSpeciesNOrbs)
        call ReadStructureFile(cfg%SCStructureFile, SCPos, SCEl, i, SCSpecies, SCSpeciesNOrbs)

        if (i /= nAtoms) then
            print '(A, I0, A, I0)', "[Unfold Input Error] " // trim(cfg%SCStructureFile) // &
                " tem numero de atomos diferente do .hsio: ", i, " vs ", nAtoms
            stop 1
        end if

        allocate(cfg%PCOrbit(nAtoms))
        allocate(cfg%CellIndex(3, nAtoms))

        do i = 1, nAtoms
            xiSC = modulo(SCPos(:, i), 1.0_dp)                 ! defensivo: garante [0,1)
            ! xi_PC = xi_SC . M (vetor-LINHA) -- NÃO matmul(M,xiSC) (vetor-coluna).
            ! Esta base de código guarda LatticeVectors como LINHAS (convenção
            ! POSCAR/ASE/CIF, ver ReadStructureFile): r = sum_i xiSC_i A_i =
            ! sum_ij xiSC_i M_ij a_j. ATENÇÃO: pontos k se transformam com a
            ! TRANSPOSTA (F_K = M . f_k, ver seção 5.1 de RunUnfolding), porque
            ! coordenadas recíprocas são contravariantes. Para M simétrica as
            ! duas formas coincidem, por isso a diferença só aparece com M
            ! não-diagonal E não-simétrica (ex.: grafeno sqrt3 x sqrt3 R30).
            xiPC = matmul(xiSC, real(cfg%M, kind=dp))

            ! "Arredonda para a fronteira" quando xiPC(j) está a menos de
            ! MatchTolerance de um inteiro, em vez de só floor() puro com um
            ! epsilon fixo minúsculo. Motivo: coordenadas de estruturas DFT
            ! reais tem ruído numérico (relaxação, precisão de escrita do
            ! arquivo) da ORDEM de MatchTolerance -- maior que qualquer
            ! epsilon fixo pequeno. Sem isso, dois átomos que fisicamente
            ! pertencem à MESMA cópia da PC (mesmo rótulo de célula) podem
            ! computar xiPC(j) em lados OPOSTOS de uma fronteira inteira por
            ! ruído independente em cada átomo (ex.: 1.999950 vs 2.000010),
            ! dando cellIndex(j) diferente para os dois -- o que quebra a
            ! consistência exigida pela tabela de grupo do motor GERAL (M
            ! não-diagonal) em RunUnfolding ("Rótulo de célula inconsistente
            ! entre orbitais-PC"). Reusa MatchTolerance (já configurável via
            ! &UNFOLD) em vez de introduzir mais um parâmetro.
            do j = 1, 3
                if (xiPC(j) - floor(xiPC(j)) > 1.0_dp - cfg%MatchTolerance) then
                    cellIndex(j) = floor(xiPC(j)) + 1   ! bem perto do próximo inteiro: arredonda p/ cima
                else
                    cellIndex(j) = floor(xiPC(j))
                end if
            end do
            remainder = xiPC - real(cellIndex, kind=dp)

            ! Casa `remainder` (posição fracionária dentro da célula-PC) contra
            ! a lista de átomos-template da PC, usando imagem mínima periódica
            ! e exigindo o MESMO símbolo químico.
            jbest = -1
            distbest = huge(1.0_dp)
            do j = 1, nAtomsPC
                if (trim(adjustl(SCEl(i))) /= trim(adjustl(PCEl(j)))) cycle
                diff = remainder - PCPos(:, j)
                diff = diff - nint(diff)       ! imagem mínima: dobra para [-0.5, 0.5)
                dist = norm2(diff)
                if (dist < distbest) then
                    distbest = dist
                    jbest = j
                end if
            end do

            if (jbest < 0 .or. distbest > cfg%MatchTolerance) then
                print '(A, I0)', "[Unfold Input Error] Nao foi possivel casar geometricamente o atomo da SC #", i
                print '(A, A)', "  elemento: ", trim(adjustl(SCEl(i)))
                print '(A, 3F10.6)', "  posicao fracionaria (SC): ", SCPos(:, i)
                print '(A, 3F10.6)', "  posicao fracionaria (PC) reduzida: ", remainder
                if (jbest > 0) then
                    print '(A, I0, A, F10.6, A, F10.6, A)', "  atomo-PC mais proximo (mesmo elemento): ", jbest, &
                        " (distancia = ", distbest, ", tolerancia = ", cfg%MatchTolerance, ")"
                else
                    print '(A)', "  Nenhum atomo-PC com o mesmo elemento quimico foi encontrado."
                end if
                print '(A)', "  Verifique se M, os vetores de rede/posicoes/simbolos e a ordem" // &
                             " dos atomos no .hsio estao consistentes, ou aumente MatchTolerance em &UNFOLD."
                stop 1
            end if

            ! Confere o número de orbitais declarado no &SPECIES (tanto do
            ! lado PC quanto do lado SC, que devem concordar entre si) contra
            ! o que o .hsio realmente tem para este átomo -- SOMENTE se
            ! ambos os arquivos de estrutura tiverem o bloco &SPECIES
            ! (opcional); caso contrário, esta checagem extra é pulada.
            if (size(PCSpeciesNOrbs) > 0 .and. size(SCSpeciesNOrbs) > 0) then
                expectedOrbs = SpeciesNOrbsLookup(PCSpecies, PCSpeciesNOrbs, PCEl(jbest))
                if (expectedOrbs /= SpeciesNOrbsLookup(SCSpecies, SCSpeciesNOrbs, SCEl(i))) then
                    print '(A, A)', "[Unfold Input Error] &SPECIES da PC e da SC discordam para a especie ", &
                        trim(adjustl(SCEl(i)))
                    stop 1
                end if
                if (expectedOrbs /= HData%Orbs(i)) then
                    print '(A, I0)', "[Unfold Input Error] Numero de orbitais inconsistente no atomo da SC #", i
                    print '(A, A)', "  especie: ", trim(adjustl(SCEl(i)))
                    print '(A, I0)', "  declarado em &SPECIES (PC/SC): ", expectedOrbs
                    print '(A, I0)', "  numero real no .hsio (Orbs): ", HData%Orbs(i)
                    print '(A)', "  Verifique se o conjunto de base usado para gerar o .hsio bate com" // &
                                 " o declarado no bloco &SPECIES, e se a ordem dos atomos esta correta."
                    stop 1
                end if
            end if

            cfg%PCOrbit(i)      = jbest
            cfg%CellIndex(:, i) = cellIndex
        end do

        if (size(PCSpeciesNOrbs) > 0 .and. size(SCSpeciesNOrbs) > 0) then
            print '(A, I0, A)', "  Mapeamento atomo->PC construido automaticamente (", nAtoms, " atomos, " // &
                                 "casamento geometrico via estrutura, orbitais/especie conferidos)."
        else
            print '(A, I0, A)', "  Mapeamento atomo->PC construido automaticamente (", nAtoms, " atomos, " // &
                                 "casamento geometrico via estrutura; &SPECIES ausente, sem checagem de orbitais)."
        end if
    end subroutine BuildAtomMapFromPositions

    !> Procura `species` na tabela `names(:)` e devolve o número de orbitais
    !> correspondente em `norbs(:)`; para com erro claro se não encontrar.
    function SpeciesNOrbsLookup(names, norbs, species) result(n)
        character(len=8), intent(in) :: names(:)
        integer,          intent(in) :: norbs(:)
        character(len=8), intent(in) :: species
        integer                      :: n, k

        n = -1
        do k = 1, size(names)
            if (trim(adjustl(names(k))) == trim(adjustl(species))) then
                n = norbs(k)
                return
            end if
        end do
        print '(A, A)', "[Unfold Input Error] Especie nao encontrada em &SPECIES: ", trim(adjustl(species))
        stop 1
    end function SpeciesNOrbsLookup

    !> Lê um arquivo de ESTRUTURA no formato .nml do próprio pacote: vetores
    !> de rede (bloco &STRUCTURE, namelist), a associação espécie->número de
    !> orbitais (bloco &SPECIES) e a lista de átomos da célula (bloco &ATOMS:
    !> Espécie X Y Z, Cartesiano, em Angstrom):
    !>
    !>   &STRUCTURE
    !>       LatticeVectors(1,:) = a1x, a1y, a1z   ! Angstrom
    !>       LatticeVectors(2,:) = a2x, a2y, a2z
    !>       LatticeVectors(3,:) = a3x, a3y, a3z
    !>   /
    !>
    !>   &SPECIES
    !>   <nEspecies>
    !>   Espécie  nOrbitais
    !>   ...
    !>   /
    !>
    !>   &ATOMS
    !>   <nAtomos>
    !>   Espécie  x  y  z    ! Cartesiano, Angstrom
    !>   ...
    !>   /
    !>
    !> Devolve as posições já convertidas para coordenadas FRACIONÁRIAS.
    !> LatticeVectors tem uma LINHA por vetor de rede (LatticeVectors(1,:)=a1,
    !> (2,:)=a2, (3,:)=a3 -- a mesma convenção de POSCAR/ASE/CIF). Como
    !> r_cart = f_frac . LatticeVectors (f_frac linha, vetores em linha),
    !> tem-se f_frac = r_cart . LatticeVectors^-1, i.e., em notação
    !> matriz-vetor coluna, f_frac = (LatticeVectors^T)^-1 . r_cart.
    !> Devolve também os símbolos químicos de cada átomo, e a tabela
    !> espécie->número de orbitais (usada por BuildAtomMapFromPositions para
    !> conferir, contra o .hsio, que a base realmente usada tem o número de
    !> orbitais esperado para cada espécie -- um jeito barato de pegar cedo
    !> erros de ordenação ou de conjunto de base).
    subroutine ReadStructureFile(filename, Pos, Elements, nAtoms, SpeciesNames, SpeciesNOrbs)
        character(len=*), intent(in)                        :: filename
        real(dp),         allocatable,          intent(out) :: Pos(:, :)
        character(len=8), allocatable,  intent(out)         :: Elements(:)
        integer,          intent(out)                       :: nAtoms
        character(len=8), allocatable,  intent(out)         :: SpeciesNames(:)
        integer,          allocatable,  intent(out)         :: SpeciesNOrbs(:)

        integer            :: iu, ierr, i, nSpecies
        character(len=256) :: line, ErrorMessage
        real(dp)           :: LatticeVectors(3, 3), InvLattice(3, 3), rcart(3)
        namelist /STRUCTURE/ LatticeVectors

        open(newunit=iu, file=trim(filename), status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
            print '(A, A)', "[Unfold Input Error] Nao foi possivel abrir arquivo de estrutura: ", trim(filename)
            stop 1
        end if

        ! 1. Bloco &STRUCTURE (namelist): LatticeVectors(3,3), LINHAS = a1,a2,a3
        !    (convenção POSCAR/ASE/CIF). f_frac = (LatticeVectors^T)^-1 . r_cart
        LatticeVectors = 0.0_dp
        read(iu, nml=STRUCTURE, iostat=ierr, iomsg=ErrorMessage)
        if (ierr /= 0) then
            print '(A, A)', "[Unfold Input Error] Falha ao ler o bloco &STRUCTURE em ", trim(filename)
            print '(A, A)', "  Mensagem: ", trim(ErrorMessage)
            stop 1
        end if
        if (all(abs(LatticeVectors) < 1.0e-12_dp)) then
            print '(A, A)', "[Unfold Input Error] LatticeVectors ausente ou nulo em ", trim(filename)
            stop 1
        end if
        InvLattice = Invert3x3(transpose(LatticeVectors))

        ! 2. Bloco &SPECIES (OPCIONAL): <nEspecies> seguido de "Espécie
        !    nOrbitais" por linha. Se ausente, pula direto para &ATOMS e a
        !    checagem espécie->número-de-orbitais é simplesmente ignorada
        !    (BuildAtomMapFromPositions detecta a tabela vazia e não confere).
        do
            read(iu, '(A)', iostat=ierr) line
            if (ierr /= 0) then
                print '(A, A)', "[Unfold Input Error] Bloco &ATOMS nao encontrado em ", trim(filename)
                stop 1
            end if
            if (trim(adjustl(line)) == "&SPECIES" .or. trim(adjustl(line)) == "&ATOMS") exit
        end do

        if (trim(adjustl(line)) == "&SPECIES") then
            call NextDataLine(iu, line, ierr)
            if (ierr == 0) read(line, *, iostat=ierr) nSpecies
            if (ierr /= 0 .or. nSpecies < 1) then
                print '(A, A)', "[Unfold Input Error] Falha ao ler o numero de especies em ", trim(filename)
                stop 1
            end if

            allocate(SpeciesNames(nSpecies))
            allocate(SpeciesNOrbs(nSpecies))
            do i = 1, nSpecies
                call NextDataLine(iu, line, ierr)
                if (ierr == 0) read(line, *, iostat=ierr) SpeciesNames(i), SpeciesNOrbs(i)
                if (ierr /= 0 .or. SpeciesNOrbs(i) < 1) then
                    print '(A, I0, A, A)', "[Unfold Input Error] Falha ao ler a especie ", i, " em ", trim(filename)
                    stop 1
                end if
            end do

            ! 3. Bloco &ATOMS: <nAtomos> seguido de "Espécie x y z" por linha
            do
                read(iu, '(A)', iostat=ierr) line
                if (ierr /= 0) then
                    print '(A, A)', "[Unfold Input Error] Bloco &ATOMS nao encontrado em ", trim(filename)
                    stop 1
                end if
                if (trim(adjustl(line)) == "&ATOMS") exit
            end do
        else
            ! &SPECIES ausente: a linha atual JÁ é "&ATOMS" (não há o que
            ! avançar); tabela de espécies fica vazia.
            allocate(SpeciesNames(0))
            allocate(SpeciesNOrbs(0))
        end if

        call NextDataLine(iu, line, ierr)
        if (ierr == 0) read(line, *, iostat=ierr) nAtoms
        if (ierr /= 0 .or. nAtoms < 1) then
            print '(A, A)', "[Unfold Input Error] Falha ao ler o numero de atomos em ", trim(filename)
            stop 1
        end if

        allocate(Pos(3, nAtoms))
        allocate(Elements(nAtoms))
        do i = 1, nAtoms
            call NextDataLine(iu, line, ierr)
            if (ierr == 0) read(line, *, iostat=ierr) Elements(i), rcart
            if (ierr /= 0) then
                print '(A, I0, A, A)', "[Unfold Input Error] Falha ao ler o atomo ", i, " em ", trim(filename)
                stop 1
            end if
            if (size(SpeciesNames) > 0 .and. .not. any(Elements(i) == SpeciesNames)) then
                print '(A, I0, A, A, A, A)', "[Unfold Input Error] Atomo ", i, " (especie '", &
                    trim(adjustl(Elements(i))), "') nao consta no bloco &SPECIES de ", trim(filename)
                stop 1
            end if
            Pos(:, i) = matmul(InvLattice, rcart)
        end do
        close(iu)
    end subroutine ReadStructureFile


    !> Inversa de uma matriz 3x3 (fórmula de cofatores). Usada para converter
    !> posições Cartesianas de arquivos .xyz em coordenadas fracionárias.
    function Invert3x3(A) result(Ainv)
        real(dp), intent(in) :: A(3, 3)
        real(dp)             :: Ainv(3, 3), det

        det = A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
            - A(1,2) * (A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
            + A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))

        if (abs(det) < 1.0e-12_dp) then
            print '(A)', "[Unfold Error] Matriz de rede singular (determinante ~ 0)."
            stop 1
        end if

        Ainv(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2)) / det
        Ainv(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2)) / det
        Ainv(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2)) / det
        Ainv(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1)) / det
        Ainv(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1)) / det
        Ainv(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1)) / det
        Ainv(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1)) / det
        Ainv(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1)) / det
        Ainv(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1)) / det
    end function Invert3x3

    !> Lê o bloco &PCPATH, no MESMO formato do bloco &BANDS de mod_config.f90,
    !> mas interpretado como coordenadas fracionárias da PC.
    subroutine ReadPCPath(iu, path)
        integer,            intent(in)    :: iu
        type(PCPathConfig), intent(inout) :: path

        character(len=256) :: line
        integer            :: ierr, i, j, idx
        real(dp)           :: dk(3), xbuffer, length

        do
            read(iu, '(A)', iostat=ierr) line
            if (ierr /= 0) then
                print '(A)', "[Unfold Input Error] Bloco &PCPATH nao encontrado"
                stop 1
            end if
            if (trim(adjustl(line)) == "&PCPATH") exit
        end do

        call NextDataLine(iu, line, ierr)
        if (ierr == 0) read(line, *, iostat=ierr) path%npaths
        if (ierr /= 0) then
            print '(A)', "[Unfold Input Error] Falha ao ler o numero de caminhos em &PCPATH"
            stop 1
        end if
        if (path%npaths < 2) then
            print '(A)', "[Unfold Input Error] &PCPATH precisa de ao menos 2 vertices"
            stop 1
        end if

        allocate(path%nkpts(path%npaths))
        allocate(path%kpoints(3, path%npaths))

        do i = 1, path%npaths
            call NextDataLine(iu, line, ierr)
            if (ierr == 0) read(line, *, iostat=ierr) path%kpoints(:, i), path%nkpts(i)
            if (ierr /= 0) then
                print '(A, I0)', "[Unfold Input Error] Falha ao ler o vertice ", i, " de &PCPATH"
                stop 1
            end if
        end do

        if (any(path%nkpts(1:path%npaths - 1) < 1)) then
            print '(A)', "[Unfold Input Error] &PCPATH: o numero de pontos de cada segmento deve ser >= 1"
            stop 1
        end if

        path%total_kpoints = sum(path%nkpts(1:path%npaths - 1))
        allocate(path%FullPath(3, path%total_kpoints))
        allocate(path%xpoints(path%total_kpoints))

        idx = 1
        xbuffer = 0.0_dp
        do i = 1, path%npaths - 1
            dk = (path%kpoints(:, i + 1) - path%kpoints(:, i)) / real(path%nkpts(i), kind=dp)
            length = norm2(dk)
            do j = 1, path%nkpts(i)
                path%FullPath(:, idx) = path%kpoints(:, i) + real(j - 1, kind=dp) * dk
                path%xpoints(idx) = xbuffer
                xbuffer = xbuffer + length
                idx = idx + 1
            end do
        end do
        path%xpoints(idx - 1)     = xbuffer
        path%FullPath(:, idx - 1) = path%kpoints(:, path%npaths)
    end subroutine ReadPCPath

    subroutine DisplayUnfoldConfig(cfg, nAtoms, nOrbs)
        use terminal_ui, only: print_box, print_row_split, print_box_close
        type(UnfoldConfig), intent(in) :: cfg
        integer,            intent(in) :: nAtoms, nOrbs
        character(len=64)              :: buf

        call print_box("Unfold Configuration")
        call print_row_split("Project Name:", cfg%ProjectName)
        call print_row_split("HData File:",   cfg%HDataFile)
        call print_row_split("SData File:",   cfg%SDataFile)
        if (cfg%DiagonalM) then
            write(buf, '(I0, 1X, I0, 1X, I0)') cfg%MDiag
            call print_row_split("SC multiplicity (m1,m2,m3):", buf)
        else
            write(buf, '(A, 3(I0,1X), A, 3(I0,1X), A, 3(I0,1X), A)') "[", cfg%M(1,:), "/ ", cfg%M(2,:), &
                "/ ", cfg%M(3,:), "]"
            call print_row_split("SC transform M (nao-diagonal):", buf)
            write(buf, '(I0)') cfg%nCells
            call print_row_split("nCells = |det(M)|:", buf)
        end if
        write(buf, '(I0)') nAtoms
        call print_row_split("Atoms in SC:", buf)
        write(buf, '(I0)') nOrbs
        call print_row_split("Orbitals in SC:", buf)
        write(buf, '(I0)') cfg%PCPath%total_kpoints
        call print_row_split("PC k-points to unfold:", buf)
        write(buf, '(F0.4, A)') cfg%Broadening, " eV"
        call print_row_split("Broadening:", buf)
        call print_box_close()
    end subroutine DisplayUnfoldConfig

    ! ======================================================================
    ! MOTOR DE UNFOLDING
    ! ======================================================================

    !> Executa o unfolding completo: monta as tabelas de órbita, faz a
    !> diagonalização generalizada em cada K da SC e projeta o resultado
    !> nos "orbitais primitivos" para obter os pesos W^k_KN.
    !>
    !> Paralelizado sobre os pontos k da PC com OpenMP, na MESMA técnica
    !> usada em calc_band_structure_non_orth (mod_core.f90): cada thread
    !> recebe seu próprio conjunto de matrizes de trabalho (pool indexado
    !> pelo id da thread), evitando qualquer condição de corrida -- cada
    !> iteração ik só escreve nas colunas ik de Result%Energies/Weights/
    !> SumCheck, que são posições disjuntas entre threads.
    subroutine RunUnfolding(cfg, HData, SData, Result, Performance)
        use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
        type(UnfoldConfig),       intent(in)  :: cfg
        type(HSDataFormat),       intent(in)  :: HData, SData
        type(UnfoldResultFormat), intent(out) :: Result
        type(UnfoldPerformance),  intent(out) :: Performance

        integer              :: nAtoms, nOrbs, nCells, nAtomsPC, nOrbitsPC
        integer, allocatable :: AtomOrbOffset(:), nOrbsPerPCAtom(:), PCAtomOrbOffset(:)
        integer, allocatable :: OrbitID(:)                ! (nOrbs)
        integer, allocatable :: OrbitCell(:, :)           ! (3, nOrbs)
        integer, allocatable :: OrbitMemberOrb(:, :)      ! (nCells, nOrbitsPC)
        integer, allocatable :: OrbitMemberCell(:, :, :)  ! (3, nCells, nOrbitsPC)
        integer, allocatable :: fillcount(:)

        ! -- SÓ usados no caso M NÃO-diagonal (ver seção 4b/5.5 abaixo) --
        integer,     allocatable :: Reps(:, :)          ! (3, nCells) rótulos-canônicos de célula
        integer,     allocatable :: OrbCellCanon(:, :)  ! (nCells, nOrbitsPC) reindexado p/ Reps
        integer,     allocatable :: OrbShift(:, :, :)   ! (3, nCells, nOrbitsPC) s tal que rótulo = Reps + M.s
        integer,     allocatable :: DestIndex(:, :)     ! (nCells, nCells) aritmética do grupo Z^3/M.Z^3
        integer,     allocatable :: WrapT(:, :, :)      ! (3, nCells, nCells) wraparound (Eq. 30)
        integer                  :: rep_idx, ga, gb, gc
        complex(dp)              :: g_inner, g_W, g_wrap_phase, g_phase_a, g_phase_b
        integer                  :: shift_s(3)
        logical,     allocatable :: filled(:, :), seen(:, :, :, :)
        integer                  :: cdiag(3)

        integer                         :: iAtom, iOrb, gOrb, orbit, L, ik, ic, iN
        integer                         :: total_kpoints, max_threads, thread_id, completed_count, progress_width
        integer                         :: threads_used
        integer(int64)                  :: clock_start, clock_end, clock_rate, calculation_start, calculation_end
        real(dp)                        :: fk(3), ScK(3), point_elapsed, point_time_mean
        real(dp),           allocatable :: point_times(:)
        character(len=32)               :: total_kpoints_text
        character(len=128)              :: progress_format
        complex(dp)                     :: proj

        ! Pools de matrizes de trabalho: uma "fatia" por thread (última
        ! dimensão = índice da thread), exatamente como em
        ! calc_band_structure_non_orth (mod_core.f90)
        complex(dp), allocatable :: Hk_pool(:, :, :), Sk_pool(:, :, :)
        complex(dp), allocatable :: Cvec_pool(:, :, :), Cprime_pool(:, :, :)
        complex(dp), allocatable :: Ssqrt_pool(:, :, :), Usq_pool(:, :, :)
        real(dp),    allocatable :: Eval_pool(:, :), wS_pool(:, :)
        complex(dp), allocatable :: phase_pool(:, :)

        nAtoms = HData%nAtoms
        nOrbs  = HData%NOrbs
        nCells = cfg%nCells

        if (nAtoms /= size(cfg%PCOrbit)) then
            print '(A)', "[Unfold Error] &ATOMMAP incompatível com o numero de atomos do .hsio"
            stop 1
        end if
        if (SData%NOrbs /= nOrbs .or. SData%nAtoms /= nAtoms) then
            print '(A)', "[Unfold Error] Arquivos H e S sao incompatíveis (NOrbs/nAtoms diferentes)"
            stop 1
        end if

        ! --------------------------------------------------------------
        ! 1. Deslocamentos de orbital por átomo (blocos contíguos, na
        !    mesma ordem de HData%Orbs, exatamente como AssembleDense usa)
        ! --------------------------------------------------------------
        allocate(AtomOrbOffset(nAtoms))
        AtomOrbOffset(1) = 0
        do iAtom = 2, nAtoms
            AtomOrbOffset(iAtom) = AtomOrbOffset(iAtom - 1) + HData%Orbs(iAtom - 1)
        end do
        if (AtomOrbOffset(nAtoms) + HData%Orbs(nAtoms) /= nOrbs) then
            print '(A)', "[Unfold Error] Soma de Orbs(:) nao bate com NOrbs do .hsio"
            stop 1
        end if

        ! --------------------------------------------------------------
        ! 2. Identifica quantos "átomos-PC" distintos existem e quantos
        !    orbitais cada um tem, e monta os deslocamentos no espaço de
        !    "identidades de orbital-PC" (1..nOrbitsPC)
        ! --------------------------------------------------------------
        nAtomsPC = maxval(cfg%PCOrbit)
        allocate(nOrbsPerPCAtom(nAtomsPC), source=0)
        do iAtom = 1, nAtoms
            if (nOrbsPerPCAtom(cfg%PCOrbit(iAtom)) == 0) then
                nOrbsPerPCAtom(cfg%PCOrbit(iAtom)) = HData%Orbs(iAtom)
            else if (nOrbsPerPCAtom(cfg%PCOrbit(iAtom)) /= HData%Orbs(iAtom)) then
                print '(A)', "[Unfold Error] Atomos com o mesmo rotulo PCOrbit tem numeros diferentes de orbitais"
                stop 1
            end if
        end do
        if (any(nOrbsPerPCAtom == 0)) then
            print '(A)', "[Unfold Error] Rotulo de atomo-PC sem nenhum atomo associado em &ATOMMAP (1..max devem ser todos usados)"
            stop 1
        end if

        allocate(PCAtomOrbOffset(nAtomsPC))
        PCAtomOrbOffset(1) = 0
        do iAtom = 2, nAtomsPC
            PCAtomOrbOffset(iAtom) = PCAtomOrbOffset(iAtom - 1) + nOrbsPerPCAtom(iAtom - 1)
        end do
        nOrbitsPC = sum(nOrbsPerPCAtom)

        if (nOrbitsPC * nCells /= nOrbs) then
            print '(A, I0, A, I0, A, I0)',                                             &
                "[Unfold Error] Incompatibilidade de contagem: nOrbitsPC(", nOrbitsPC, &
                ") * nCells(", nCells, ") != NOrbs(", nOrbs, ")"
            stop 1
        end if

        ! --------------------------------------------------------------
        ! 3. Expande o mapa átomo->(PCOrbit, célula) para o nível de
        !    orbital, atribuindo a cada orbital da SC sua "identidade de
        !    orbital-PC" (1..nOrbitsPC) e sua célula (n1,n2,n3)
        ! --------------------------------------------------------------
        allocate(OrbitID(nOrbs))
        allocate(OrbitCell(3, nOrbs))
        do iAtom = 1, nAtoms
            do iOrb = 1, HData%Orbs(iAtom)
                gOrb = AtomOrbOffset(iAtom) + iOrb
                OrbitID(gOrb)    = PCAtomOrbOffset(cfg%PCOrbit(iAtom)) + iOrb
                OrbitCell(:, gOrb) = cfg%CellIndex(:, iAtom)
            end do
        end do

        ! --------------------------------------------------------------
        ! 4. Tabela inversa: para cada identidade de orbital-PC, a lista
        !    dos nCells orbitais da SC que são suas réplicas, e a célula
        !    de cada um
        ! --------------------------------------------------------------
        allocate(OrbitMemberOrb(nCells, nOrbitsPC))
        allocate(OrbitMemberCell(3, nCells, nOrbitsPC))
        allocate(fillcount(nOrbitsPC), source=0)
        if (cfg%DiagonalM) &
            allocate(seen(0:cfg%MDiag(1)-1, 0:cfg%MDiag(2)-1, 0:cfg%MDiag(3)-1, nOrbitsPC), source=.false.)

        do gOrb = 1, nOrbs
            orbit = OrbitID(gOrb)
            fillcount(orbit) = fillcount(orbit) + 1
            if (fillcount(orbit) > nCells) then
                print '(A, I0)', "[Unfold Error] Orbital-PC com mais de nCells replicas: orbit ", orbit
                stop 1
            end if
            if (cfg%DiagonalM) then
                ! Rótulos podem estar fora de 0..MDiag-1 (representantes deslocados
                ! por um vetor da SC): reduz módulo MDiag antes de indexar seen.
                cdiag = modulo(OrbitCell(:, gOrb), cfg%MDiag)
                if (seen(cdiag(1), cdiag(2), cdiag(3), orbit)) then
                    print '(A, I0)', "[Unfold Error] Rotulo de celula duplicado dentro do orbital-PC ", orbit
                    stop 1
                end if
                seen(cdiag(1), cdiag(2), cdiag(3), orbit) = .true.
            end if
            OrbitMemberOrb(fillcount(orbit), orbit)    = gOrb
            OrbitMemberCell(:, fillcount(orbit), orbit) = OrbitCell(:, gOrb)
        end do
        if (any(fillcount /= nCells)) then
            print '(A)', "[Unfold Error] Alguma identidade de orbital-PC nao tem exatamente nCells replicas"
            stop 1
        end if

        ! --------------------------------------------------------------
        ! 4b. [SÓ PARA M NÃO-DIAGONAL] Pré-computa a "aritmética de
        !     grupo" dos rótulos de célula-PC, elementos do grupo
        !     Z^3/M.Z^3 de ordem nCells [Sec. II A/D do artigo].
        !
        !     Fisicamente, cada réplica de célula-PC dentro da SC contém
        !     exatamente um átomo de cada identidade-PC, então o CONJUNTO
        !     de rótulos de célula (Reps) é o MESMO para qualquer
        !     orbital-PC -- usamos o do orbital 1 como referência
        !     canônica e reindexamos os demais orbitais contra ele
        !     (OrbCellCanon), já que a ORDEM de preenchimento em
        !     OrbitMemberCell/OrbitMemberOrb (feita em 3-4 acima) não é
        !     necessariamente a mesma entre orbitais distintos.
        !
        !     DestIndex(a,c)/WrapT(:,a,c) respondem: "partindo do rótulo
        !     Reps(:,a), transladando pelo vetor-PC Reps(:,c) (tratado
        !     como elemento do grupo), a qual rótulo Reps(:,b) eu chego,
        !     e quantas vezes cruzei a fronteira da SC em cada direção
        !     (WrapT, usado na fase de wraparound e^{iK.T} da Eq. 30)?"
        !     Construído UMA VEZ (não depende do ponto k), reusado em
        !     todos os pontos k do loop paralelo abaixo (seção 5.5).
        ! --------------------------------------------------------------
        if (.not. cfg%DiagonalM) then
            allocate(Reps(3, nCells))
            Reps = OrbitMemberCell(:, :, 1)

            allocate(OrbCellCanon(nCells, nOrbitsPC))
            allocate(OrbShift(3, nCells, nOrbitsPC), source=0)
            allocate(filled(nCells, nOrbitsPC), source=.false.)
            do orbit = 1, nOrbitsPC
                do ic = 1, nCells
                    rep_idx = FindCellIndex(Reps, OrbitMemberCell(:, ic, orbit), nCells, transpose(cfg%M), shift_s)
                    if (rep_idx < 0) then
                        print '(A)', "[Unfold Error] Rotulo de celula inconsistente entre orbitais-PC " // &
                                       "(M nao-diagonal) -- verifique o &ATOMMAP"
                        stop 1
                    end if
                    if (filled(rep_idx, orbit)) then
                        print '(A, I0)', "[Unfold Error] Rotulo de celula duplicado dentro do orbital-PC ", orbit
                        stop 1
                    end if
                    filled(rep_idx, orbit) = .true.
                    OrbCellCanon(rep_idx, orbit) = OrbitMemberOrb(ic, orbit)
                    ! o rótulo real é Reps(rep_idx) + M.s: guarda s para a fase de
                    ! gauge exp(-i.2pi.ScK.s) aplicada a C' no loop k (seção 5.5)
                    OrbShift(:, rep_idx, orbit) = shift_s
                end do
            end do

            deallocate(filled)

            call BuildGroupTable(transpose(cfg%M), Reps, nCells, DestIndex, WrapT)
        end if

        ! --------------------------------------------------------------
        ! 5. Loop principal sobre os pontos k da PC (paralelo, OpenMP)
        ! --------------------------------------------------------------
        total_kpoints = cfg%PCPath%total_kpoints
        allocate(Result%Energies(nOrbs, total_kpoints))
        allocate(Result%Weights(nOrbs, total_kpoints))
        allocate(Result%SumCheck(total_kpoints))
        allocate(point_times(total_kpoints), source=0.0_dp)

        max_threads = 1
        !$ max_threads = omp_get_max_threads()

        allocate(Hk_pool(nOrbs, nOrbs, max_threads))
        allocate(Sk_pool(nOrbs, nOrbs, max_threads))
        allocate(Cvec_pool(nOrbs, nOrbs, max_threads))
        allocate(Cprime_pool(nOrbs, nOrbs, max_threads))
        allocate(Ssqrt_pool(nOrbs, nOrbs, max_threads))
        allocate(Usq_pool(nOrbs, nOrbs, max_threads))
        allocate(Eval_pool(nOrbs, max_threads))
        allocate(wS_pool(nOrbs, max_threads))
        allocate(phase_pool(nCells, max_threads))

        ! Relata o progresso imediatamente após a conclusão de cada ponto k.
        ! Como o loop é paralelo, a ordem de conclusão pode diferir de ik.
        completed_count = 0
        threads_used = 1
        write(total_kpoints_text, '(I0)') total_kpoints
        progress_width = len_trim(total_kpoints_text)
        write(progress_format, '("(A,I",I0,",A,I",I0,",A,F6.1,A,F10.3,A)")') &
            progress_width, progress_width
        call system_clock(count_rate=clock_rate)
        print '(A, I0, A, I0, A)', "  Diagonalizando ", total_kpoints, " pontos k (ate ", &
                                    max_threads, " threads OpenMP)..."
        call system_clock(count=calculation_start)

        !$omp parallel default(shared) &
        !$omp private(ik, thread_id, fk, ScK, L, iN, orbit, ic, proj, clock_start, clock_end, point_elapsed) &
        !$omp private(ga, gb, gc, g_inner, g_W, g_wrap_phase)
        !$omp single
        !$ threads_used = omp_get_num_threads()
        !$omp end single
        !$omp do
        do ik = 1, total_kpoints
            call system_clock(count=clock_start)

            thread_id = 0
            !$ thread_id = omp_get_thread_num()
            thread_id = thread_id + 1

            fk = cfg%PCPath%FullPath(:, ik)

            ! 5.1 Ponto K da SC ao qual este k da PC é mapeado [Eq. (7)]. Com
            !     A_i = sum_j M_ij a_j (linhas), F_K = M . f_k (vetor-COLUNA);
            !     NÃO f_k . M, que só coincide quando M é simétrica.
            if (cfg%DiagonalM) then
                do L = 1, 3
                    ScK(L) = modulo(fk(L) * real(cfg%MDiag(L), kind=dp), 1.0_dp)
                end do
            else
                ScK = modulo(matmul(real(cfg%M, kind=dp), fk), 1.0_dp)
            end if

            ! 5.2 Monta H(K) e S(K) por transformada de Fourier do espaço real
            call compute_Hk_at_kpoint(HData%Data, HData%IndexList, Hk_pool(:, :, thread_id), ScK)
            call compute_Hk_at_kpoint(SData%Data, SData%IndexList, Sk_pool(:, :, thread_id), ScK)

            ! 5.3 Resolve o problema de autovalor generalizado H(K) C = E S(K) C
            call SolveGeneralizedEigenproblem(Hk_pool(:, :, thread_id), Sk_pool(:, :, thread_id), &
                nOrbs, Eval_pool(:, thread_id), Cvec_pool(:, :, thread_id))

            ! 5.4 Transformação de Lowdin: C' = S(K)^(1/2) C  [Eq. (24)]
            call MatrixSqrtHermitian(Sk_pool(:, :, thread_id), nOrbs, &
                Ssqrt_pool(:, :, thread_id), Usq_pool(:, :, thread_id), wS_pool(:, thread_id))
            Cprime_pool(:, :, thread_id) = matmul(Ssqrt_pool(:, :, thread_id), Cvec_pool(:, :, thread_id))

            ! 5.5 Projeta cada autoestado da SC nos p orbitais-PC.
            !
            !     CASO M DIAGONAL [Eq. (28)]: W^k_KN = sum_orbit |F_{k,orbit}^dagger . C'_KN|^2,
            !     com F_k na forma fechada das Eqs. (33)-(38) (produto de exponenciais
            !     independente por direção) -- O(nCells) por orbital-PC.
            !
            !     CASO M GERAL: essa forma fechada NÃO vale (o grupo Z^3/M.Z^3 não se
            !     fatora em 3 ciclos independentes por eixo quando M não é diagonal --
            !     ver Sec. II D do artigo, Eqs. 39-41). Em vez de diagonalizar as
            !     matrizes de translação passo-a-passo como o artigo faz na sua
            !     implementação (FHI-aims), usamos a fórmula do OPERADOR DE PROJEÇÃO
            !     diretamente [Eq. (14) substituida em Eq. (10), equivalência provada
            !     no Apêndice A do artigo]:
            !
            !         W^k_KN = (1/nCells) * sum_c exp(-i.2pi.fk.Reps(c)) * <KN|t_c|KN>
            !
            !     onde t_c é a translação (na base de Bloch de Lowdin, já ortonormal
            !     -- Eq. 27) pelo elemento de grupo Reps(c), cuja AÇÃO sobre C' é
            !     dada pela tabela pre-computada DestIndex/WrapT (seção 4b): manda o
            !     orbital na célula Reps(a) para a célula Reps(DestIndex(a,c)), com
            !     fase de wraparound exp(i.2pi.ScK.WrapT(:,a,c)) [Eq. 30] sempre que a
            !     translação cruza a fronteira da SC. Isso é O(nCells) POR PAR (a,c)
            !     -> O(nCells^2) por orbital-PC (mais caro que o caso diagonal, mas
            !     matematicamente exato para qualquer M -- ver validação no exemplo
            !     sintético UnfoldSquare2AtomNonDiag).
            do iN = 1, nOrbs
                Result%Energies(iN, ik) = Eval_pool(iN, thread_id)
                Result%Weights(iN, ik)  = 0.0_dp

                if (cfg%DiagonalM) then
                    do orbit = 1, nOrbitsPC
                        do ic = 1, nCells
                            associate(nvec => OrbitMemberCell(:, ic, orbit))
                                phase_pool(ic, thread_id) =                                                     &
                                    exp(-cmplx(0.0_dp, TWO_PI * dot_product(fk, real(nvec, kind=dp)), kind=dp)) &
                                    / sqrt(real(nCells, kind=dp))
                            end associate
                        end do
                        proj = cmplx(0.0_dp, 0.0_dp, kind=dp)
                        do ic = 1, nCells
                            proj = proj + phase_pool(ic, thread_id) &
                                   * Cprime_pool(OrbitMemberOrb(ic, orbit), iN, thread_id)
                        end do
                        Result%Weights(iN, ik) = Result%Weights(iN, ik) + abs(proj)**2
                    end do
                else
                    g_W = cmplx(0.0_dp, 0.0_dp, kind=dp)
                    do gc = 1, nCells
                        g_inner = cmplx(0.0_dp, 0.0_dp, kind=dp)
                        do orbit = 1, nOrbitsPC
                            do ga = 1, nCells
                                gb = DestIndex(ga, gc)
                                ! Leva C' do orbital no rótulo real Reps+M.s ao representante
                                ! Reps (fase de gauge, como exp(-i.2pi.fk.n) no caso diagonal)
                                g_phase_a = exp(-cmplx(0.0_dp, &
                                    TWO_PI * dot_product(ScK, real(OrbShift(:, ga, orbit), kind=dp)), kind=dp))
                                g_phase_b = exp(-cmplx(0.0_dp, &
                                    TWO_PI * dot_product(ScK, real(OrbShift(:, gb, orbit), kind=dp)), kind=dp))
                                g_wrap_phase = exp(cmplx(0.0_dp, &
                                    TWO_PI * dot_product(ScK, real(WrapT(:, ga, gc), kind=dp)), kind=dp))
                                g_inner = g_inner                                                            &
                                    + conjg(g_phase_a * Cprime_pool(OrbCellCanon(ga, orbit), iN, thread_id)) &
                                    * g_wrap_phase * g_phase_b * Cprime_pool(OrbCellCanon(gb, orbit), iN, thread_id)
                            end do
                        end do
                        g_W = g_W + exp(-cmplx(0.0_dp, &
                            TWO_PI * dot_product(fk, real(Reps(:, gc), kind=dp)), kind=dp)) * g_inner
                    end do
                    Result%Weights(iN, ik) = real(g_W / real(nCells, kind=dp), kind=dp)
                end if
            end do

            Result%SumCheck(ik) = sum(Result%Weights(:, ik))

            call system_clock(count=clock_end)
            point_elapsed = real(clock_end - clock_start, dp) / real(clock_rate, dp)
            point_times(ik) = point_elapsed

            !$omp critical (unfold_progress_report)
            completed_count = completed_count + 1
            write(*, progress_format)                                            &
                "    Progresso [", completed_count, " / ", total_kpoints, " | ", &
                100.0_dp * real(completed_count, dp) / real(total_kpoints, dp),  &
                "% ]  tempo do ponto: ", point_elapsed, " s"
            flush(output_unit)
            !$omp end critical (unfold_progress_report)
        end do
        !$omp end do
        !$omp end parallel
        call system_clock(count=calculation_end)

        Performance%PointsCalculated = total_kpoints
        Performance%ThreadsUsed = threads_used
        Performance%CalculationTime = real(calculation_end - calculation_start, dp) / real(clock_rate, dp)

        if (total_kpoints > 0) then
            point_time_mean = sum(point_times) / real(total_kpoints, dp)
            Performance%PointTimeMin = minval(point_times)
            Performance%PointTimeMean = point_time_mean
            Performance%PointTimeMax = maxval(point_times)
            Performance%PointTimeStdDev = sqrt(sum((point_times - point_time_mean)**2) / real(total_kpoints, dp))

            if (Performance%CalculationTime > 0.0_dp) then
                Performance%EffectiveTimePerPoint = Performance%CalculationTime / real(total_kpoints, dp)
                Performance%ThroughputPointsPerSecond = real(total_kpoints, dp) / Performance%CalculationTime
                Performance%ParallelEfficiency = sum(point_times) / &
                    (real(max(1, threads_used), dp) * Performance%CalculationTime)
            end if
        end if

        deallocate(Hk_pool, Sk_pool, Cvec_pool, Cprime_pool, Ssqrt_pool, Usq_pool, &
                   Eval_pool, wS_pool, phase_pool, point_times)
        if (allocated(Reps))          deallocate(Reps)
        if (allocated(OrbCellCanon))  deallocate(OrbCellCanon)
        if (allocated(OrbShift))      deallocate(OrbShift)
        if (allocated(DestIndex))     deallocate(DestIndex)
        if (allocated(WrapT))         deallocate(WrapT)
    end subroutine RunUnfolding

    !> Procura `target` entre as colunas de `Reps(3,nCells)`; devolve o
    !> índice (1..nCells) ou -1 se não encontrado. Usado só no caso M
    !> não-diagonal para reindexar OrbitMemberCell/Orb contra o conjunto
    !> canônico de rótulos de célula (ver seção 4b de RunUnfolding).
    !> Procura, entre as colunas de `Reps(3,nCells)`, aquela que representa a
    !> MESMA classe lateral (coset) de Z^3/M.Z^3 que `target` -- ou seja,
    !> `target - Reps(:,k)` é um vetor de rede da SC (em unidades da PC).
    !> Devolve o índice (1..nCells) ou -1 se não encontrado.
    !>
    !> Usa equivalência de GRUPO (tolerante, via inversa de M), não
    !> igualdade exata: rótulos de célula calculados independentemente por
    !> átomo (BuildAtomMapFromPositions, um floor() por átomo) podem cair em
    !> QUALQUER representante válido da mesma classe lateral -- não
    !> necessariamente o mesmo escolhido para outro átomo-PC -- em
    !> estruturas reais onde os átomos foram "embrulhados" de volta pra
    !> célula unitária da SC (comportamento padrão de praticamente
    !> qualquer gerador de supercelula). Comparar só por igualdade exata
    !> rejeita rótulos legítimos com o erro espúrio "Rótulo de célula
    !> inconsistente entre orbitais-PC".
    function FindCellIndex(Reps, target, nCells, M, shift) result(idx)
        integer, intent(in) :: nCells
        !> M aqui é a matriz cujas COLUNAS são os vetores da SC em unidades
        !> inteiras da PC, i.e. transpose(cfg%M) (os vetores da SC são as
        !> LINHAS de cfg%M). O chamador passa transpose(cfg%M).
        integer, intent(in) :: Reps(3, nCells), target(3), M(3, 3)
        !> Opcional: s inteiro tal que target = Reps(:,idx) + M.s (0 se igualdade exata)
        integer,  intent(out), optional :: shift(3)
        integer                         :: idx, k
        real(dp)                        :: Minv(3, 3), tvec(3)

        Minv = Invert3x3(real(M, kind=dp))
        idx = -1
        if (present(shift)) shift = 0
        do k = 1, nCells
            ! mesma convenção de BuildGroupTable (Minv@diff, NÃO diff@Minv --
            ! ver nota lá)
            tvec = matmul(Minv, real(target - Reps(:, k), kind=dp))
            if (all(abs(tvec - real(nint(tvec), kind=dp)) < 1.0e-6_dp)) then
                idx = k
                if (present(shift)) shift = nint(tvec)
                return
            end if
        end do
    end function FindCellIndex

    !> Constroi a tabela de "aritmética de grupo" de Z^3/M.Z^3 (ordem
    !> nCells = |det M|) sobre o conjunto de rótulos canônicos Reps:
    !> para cada par (a,c), encontra b e o inteiro t tais que
    !>   Reps(:,a) + Reps(:,c) - matmul(M, t) = Reps(:,b)
    !> (i.e., soma os dois rótulos como elementos do grupo e reduz de
    !> volta ao conjunto de representantes). t é o número de vezes que a
    !> translação cruzou a fronteira da SC em cada direção, usado para a
    !> fase de wraparound e^{i.K.T} [Eq. (30) do artigo] em RunUnfolding.
    !>
    !> Busca por força bruta (O(nCells) candidatos por par, O(nCells^2)
    !> pares -> O(nCells^3) no total) via inversa real de M: como M é
    !> inteira e não-singular, para o candidato correto o vetor
    !> M^{-1}.(soma - Reps(:,b)) fica bem próximo de um inteiro; para os
    !> demais candidatos, longe. Só roda UMA vez por cálculo (não por
    !> ponto k), então o custo O(nCells^3) é irrelevante na prática.
    subroutine BuildGroupTable(M, Reps, nCells, DestIndex, WrapT)
        integer, intent(in)               :: M(3, 3), nCells
        integer, intent(in)               :: Reps(3, nCells)
        integer, allocatable, intent(out) :: DestIndex(:, :)
        integer, allocatable, intent(out) :: WrapT(:, :, :)

        real(dp) :: Minv(3, 3), tvec(3)
        integer  :: a, c, b, raw(3)
        logical  :: found

        Minv = Invert3x3(real(M, kind=dp))
        allocate(DestIndex(nCells, nCells))
        allocate(WrapT(3, nCells, nCells))

        do a = 1, nCells
            do c = 1, nCells
                raw = Reps(:, a) + Reps(:, c)
                found = .false.
                do b = 1, nCells
                    ! diff = raw - Reps(b) é o vetor SC-latice (candidato) EXPRESSO
                    ! em coeficientes INTEIROS da base da PC (mesmas unidades de
                    ! Reps). tvec = Minv . diff dá os coeficientes inteiros s tais
                    ! que diff = M . s, onde o argumento M é transpose(cfg%M): suas
                    ! COLUNAS são os vetores da SC em unidades da PC. (Com cfg%M
                    ! sem transpor, o resultado só é correto para M simétrica --
                    ! erro que a cadeia 1D cisalhada não detecta, pois y é
                    ! desacoplado; o teste test_unfold_core com grafeno sqrt3
                    ! detecta.)
                    tvec = matmul(Minv, real(raw - Reps(:, b), kind=dp))
                    if (all(abs(tvec - real(nint(tvec), kind=dp)) < 1.0e-6_dp)) then
                        DestIndex(a, c) = b
                        WrapT(:, a, c) = nint(tvec)
                        found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    print '(A, I0, A, I0)',                                                                   &
                        "[Unfold Error] Falha ao construir a tabela de grupo (M nao-diagonal) para (a,c) = ", &
                        a, ", ", c
                    stop 1
                end if
            end do
        end do
    end subroutine BuildGroupTable

    ! ======================================================================
    ! ÁLGEBRA LINEAR AUXILIAR (LAPACK)
    ! ======================================================================

    !> Resolve H C = E S C (Hermitiano generalizado) retornando TODOS os
    !> autovalores e autovetores (necessário para o unfolding, ao contrário
    !> de calc_band_structure_non_orth em mod_core.f90, que só guarda E).
    subroutine SolveGeneralizedEigenproblem(H, S, n, Eval, Cvec)
        integer,     intent(in)  :: n
        complex(dp), intent(in)  :: H(n, n), S(n, n)
        real(dp),    intent(out) :: Eval(n)
        complex(dp), intent(out) :: Cvec(n, n)

        complex(dp), allocatable :: Awork(:, :), Bwork(:, :), work(:)
        real(dp),    allocatable :: rwork(:)
        integer                  :: lwork, info

        interface
            subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
                use types, only: dp
                character   :: jobz, uplo
                integer     :: itype, n, lda, ldb, lwork, info
                complex(dp) :: a(lda, *), b(ldb, *), work(*)
                real(dp)    :: w(*), rwork(*)
            end subroutine zhegv
        end interface

        allocate(Awork(n, n), Bwork(n, n))
        Awork = H
        Bwork = S

        lwork = max(1, 2 * n - 1)
        allocate(work(lwork), rwork(max(1, 3 * n - 2)))

        call zhegv(1, 'V', 'U', n, Awork, n, Bwork, n, Eval, work, lwork, rwork, info)
        if (info /= 0) then
            print '(A, I0)', "[Unfold Error] zhegv falhou com info = ", info
            stop 1
        end if

        Cvec = Awork
        deallocate(Awork, Bwork, work, rwork)
    end subroutine SolveGeneralizedEigenproblem

    !> Diagonaliza a matriz Hermitiana S e devolve S^(1/2) = U * sqrt(w) * U^dagger,
    !> junto com os autovetores U e autovalores w (usados só para diagnóstico).
    subroutine MatrixSqrtHermitian(S, n, Ssqrt, U, w)
        integer,     intent(in)  :: n
        complex(dp), intent(in)  :: S(n, n)
        complex(dp), intent(out) :: Ssqrt(n, n), U(n, n)
        real(dp),    intent(out) :: w(n)

        complex(dp), allocatable :: Awork(:, :), work(:), scaled(:, :)
        real(dp),    allocatable :: rwork(:)
        integer                  :: lwork, info, i
        real(dp),    parameter   :: EPS_FLOOR = 1.0e-10_dp

        interface
            subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
                use types, only: dp
                character   :: jobz, uplo
                integer     :: n, lda, lwork, info
                complex(dp) :: a(lda, *), work(*)
                real(dp)    :: w(*), rwork(*)
            end subroutine zheev
        end interface

        allocate(Awork(n, n))
        Awork = S

        lwork = max(1, 2 * n - 1)
        allocate(work(lwork), rwork(max(1, 3 * n - 2)))

        call zheev('V', 'U', n, Awork, n, w, work, lwork, rwork, info)
        if (info /= 0) then
            print '(A, I0)', "[Unfold Error] zheev (S^1/2) falhou com info = ", info
            stop 1
        end if
        U = Awork

        ! Autovalores do overlap devem ser positivos; ruído numérico próximo
        ! de zero é truncado para evitar NaN/instabilidade em sqrt().
        do i = 1, n
            if (w(i) < EPS_FLOOR) then
                w(i) = EPS_FLOOR
            end if
        end do

        allocate(scaled(n, n))
        do i = 1, n
            scaled(:, i) = U(:, i) * sqrt(w(i))
        end do
        Ssqrt = matmul(scaled, conjg(transpose(U)))

        deallocate(Awork, work, rwork, scaled)
    end subroutine MatrixSqrtHermitian

    ! ======================================================================
    ! SAÍDA DE RESULTADOS
    ! ======================================================================

    !> Salva (x, E, W) em formato "scatter": um bloco por ponto k, separado
    !> por linha em branco (compatível com `splot ... with points` no gnuplot,
    !> usando a coluna 3 (peso) como cor via `lc palette`).
    subroutine SaveScatter(this, filename, PCPath)
        class(UnfoldResultFormat), intent(in) :: this
        character(len=*),          intent(in) :: filename
        type(PCPathConfig),        intent(in) :: PCPath

        integer :: iunit, ios, nOrbs, nk, ik, iN

        nOrbs = size(this%Energies, 1)
        nk    = size(this%Energies, 2)

        open(newunit=iunit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            print *, "Error: nao foi possivel abrir para escrita: ", filename
            stop 1
        end if

        write(iunit, '(A)') "# Unfolding LowFoldIn -- PHYSICAL REVIEW B 113, 085112 (2026)"
        write(iunit, '(A)') "# Colunas: [1] posicao x no caminho PC | [2] Energia (eV) | [3] Peso W^k_KN"
        do ik = 1, nk
            do iN = 1, nOrbs
                write(iunit, '(F14.8, F14.6, ES16.6)') PCPath%xpoints(ik), &
                    this%Energies(iN, ik), this%Weights(iN, ik)
            end do
            if (ik < nk) write(iunit, *)
        end do
        close(iunit)

        print *, "Unfolding (scatter) salvo em: ", filename
        print *, 'Use no gnuplot: '
        print *, '  set palette defined (0 "white", 1 "blue")'
        print *, '  plot "', trim(filename), '" using 1:2:3 with points pt 7 ps 0.6 lc palette notitle'
    end subroutine SaveScatter

    !> Calcula e salva a função espectral A(k,E) [Eq. (12)] em uma grade de
    !> energias, aproximando a delta de Dirac por uma Lorentziana de largura
    !> `cfg%Broadening`, no formato de triplas (x, E, A) pronto para
    !> `splot ... with pm3d` / `plot ... with image`.
    subroutine SaveSpectralFunction(this, filename, PCPath, Emin, Emax, nE, Broadening)
        class(UnfoldResultFormat), intent(in) :: this
        character(len=*),          intent(in) :: filename
        type(PCPathConfig),        intent(in) :: PCPath
        real(dp),                  intent(in) :: Emin, Emax, Broadening
        integer,                   intent(in) :: nE

        integer  :: iunit, ios, nOrbs, nk, ik, iN, iE
        real(dp) :: E, dE, A, gamma_half

        nOrbs = size(this%Energies, 1)
        nk    = size(this%Energies, 2)
        dE = (Emax - Emin) / real(max(nE - 1, 1), kind=dp)
        gamma_half = 0.5_dp * Broadening

        open(newunit=iunit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            print *, "Error: nao foi possivel abrir para escrita: ", filename
            stop 1
        end if

        write(iunit, '(A)') "# Funcao espectral A(k,E) [Eq. (12), PRB 113, 085112 (2026)]"
        write(iunit, '(A, F0.4, A)') "# Broadening Lorentziano: ", Broadening, " eV"
        write(iunit, '(A)') "# Colunas: [1] posicao x no caminho PC | [2] Energia (eV) | [3] A(k,E)"

        do ik = 1, nk
            do iE = 1, nE
                E = Emin + real(iE - 1, kind=dp) * dE
                A = 0.0_dp
                do iN = 1, nOrbs
                    A = A + this%Weights(iN, ik) * (gamma_half / 3.14159265358979323846_dp) &
                            / ((E - this%Energies(iN, ik))**2 + gamma_half**2)
                end do
                write(iunit, '(F14.8, F14.6, ES16.6)') PCPath%xpoints(ik), E, A
            end do
            write(iunit, *)
        end do
        close(iunit)

        print *, "Funcao espectral salva em: ", filename
        print *, 'Use no gnuplot: '
        print *, '  set pm3d map'
        print *, '  splot "', trim(filename), '" using 1:2:3 notitle'
    end subroutine SaveSpectralFunction

end module unfold
