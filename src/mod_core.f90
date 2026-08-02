module core
    use types, only: sp, dp, qp
    use hs_io
    use config
    use terminal_ui
    
    private
    public :: calc_band_structure, calc_band_structure_non_orth, compute_Hk_at_kpoint

    ! ==================================================================
    ! INTERFACES EXPLÍCITAS PARA AS ROTINAS DO LAPACK (Double Precision)
    ! ==================================================================
    interface
        ! 1. Ortogonal Clássica (QR/QL)
        subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
            use types, only: dp
            character   :: jobz, uplo
            integer     :: n, lda, lwork, info
            complex(dp) :: a(lda, *), work(*)
            real(dp)    :: w(*), rwork(*)
        end subroutine zheev

        ! 2. Ortogonal Ultra-Rápida (Divide & Conquer) - Recomendada para NOrbs ~ 2000
        subroutine zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
            use types, only: dp
            character   :: jobz, uplo
            integer     :: n, lda, lwork, lrwork, liwork, info
            complex(dp) :: a(lda, *), work(*)
            real(dp)    :: w(*), rwork(*)
            integer     :: iwork(*)
        end subroutine zheevd

        ! 3. Não Ortogonal Clássica (Problema de Autovalores Generalizado)
        subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
            use types, only: dp
            character   :: jobz, uplo
            integer     :: itype, n, lda, ldb, lwork, info
            complex(dp) :: a(lda, *), b(ldb, *), work(*)
            real(dp)    :: w(*), rwork(*)
        end subroutine zhegv

        ! 4. Não Ortogonal Ultra-Rápida (Divide & Conquer Generalizado)
        subroutine zhegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info)
            use types, only: dp
            character   :: jobz, uplo
            integer     :: itype, n, lda, ldb, lwork, lrwork, liwork, info
            complex(dp) :: a(lda, *), b(ldb, *), work(*)
            real(dp)    :: w(*), rwork(*)
            integer     :: iwork(*)
        end subroutine zhegvd

    end interface



contains
    
    function calc_band_structure(cfg) result(Energies)
        use types
        use omp_lib, only: omp_get_max_threads, omp_get_thread_num

        type(AppConfig),  intent(in)         :: cfg
        real(kind=dp),    allocatable        :: Energies(:, :)
        
        type(HSDataFormat) :: FileData
        integer            :: ik

        ! 1. Variaveis para Gerenciar a Pré-alocação alinhado com openmp
        integer :: max_threads, thread_id

        ! 2. Matrizes de Trabalho para LAPACK com Dimensão extra para openmp
        complex(dp), allocatable :: Hk_pool(:, :, :)
        real(dp),    allocatable :: w_pool(:, :)
        complex(dp), allocatable :: work_pool(:, :)
        real(dp),    allocatable :: rwork_pool(:, :)
        integer,     allocatable :: info_pool(:)
        integer                  :: lwork

        ! 3. Determina o número máximo de threads disponíveis para paralelização
        !    Se o OpenMP não estiver habilitado, max_threads será definido como 1
        !    Caso contrário, será definido como o número máximo de threads disponíveis
        !    Isso garante que cada thread tenha suas próprias matrizes de trabalho para evitar condições de corrida
        !    e problemas de concorrência durante a diagonalização paralela
        !    A variável max_threads é usada para alocar matrizes de trabalho separadas para cada thread, 
        !garantindo que cada thread tenha seu próprio espaço de memória para cálculos independentes
        max_threads = 1
        !$ max_threads = omp_get_max_threads()
        
        ! 4. Carrega os dados do arquivo
        call FileData%ReadFile(cfg%HDataFile)

        ! Tamanho de lwork = 3 * NOrbs - 1, conforme exigido pelo LAPACK
        lwork = 3 * FileData%NOrbs - 1

        associate( Hr        => FileData%Data,           &
                   nOrbs     => FileData%NOrbs,          & 
                   nkPoints  => cfg%Bands%total_kpoints, &
                   FullPath  => cfg%Bands%FullPath,      &
                   Rpt       => FileData%IndexList       & ! Rpt é a lista de vetores de vizinhos (h, k, l)
        )
            ! Aloca a matriz completa de uma só vez
            allocate(Energies(nOrbs, nkPoints))
            Energies = 0.0d0

            ! 5. Aloca as matrizes de trabalho para cada thread uma unica vez
            allocate(Hk_pool(nOrbs, nOrbs, max_threads))
            allocate(w_pool(nOrbs, max_threads))
            allocate(work_pool(3*nOrbs - 1, max_threads))
            allocate(rwork_pool(3*nOrbs - 2, max_threads))
            allocate(info_pool(max_threads))

            ! 6. Loop paralelo sobre os pontos k
            !$omp parallel do default(shared) private(ik, thread_id)
            do ik = 1, nkPoints

                ! 1. Obtém o ID da thread atual (0-based)
                thread_id = 0
                !$ thread_id = omp_get_thread_num()

                ! 2. Ajusta o indice para fortran (1-based) para acessar a matriz de trabalho da thread correta
                thread_id = thread_id + 1

                ! 3. Calcula a matriz H(k) para o ponto k atual
                call compute_Hk_at_kpoint(Hr, Rpt, Hk_pool(:, :, thread_id), FullPath(:, ik))

                ! 4. Diagonaliza H(k) usando LAPACK
                call zheev('N', 'U', nOrbs, Hk_pool(:, :, thread_id), nOrbs, & 
                            w_pool(:, thread_id), work_pool(:, thread_id),   &
                            lwork, rwork_pool(:, thread_id), info_pool(thread_id))

                ! 5. Verifica se a diagonalização foi bem-sucedida
                if (info_pool(thread_id) /= 0) then
                    print *, "Erro: LAPACK falhou no ponto k:", ik, " Thread:", thread_id
                endif

                ! 6. Armazena os autovalores (energias) no array de saída
                Energies(:, ik) = w_pool(:, thread_id)

            end do
            !$omp end parallel do

            ! 7. Desaloca as matrizes de trabalho após o uso
            deallocate(Hk_pool, w_pool, work_pool, rwork_pool, info_pool)
        end associate
    end function calc_band_structure

    function calc_band_structure_non_orth(cfg) result(Energies)
        use types
        use omp_lib, only: omp_get_max_threads, omp_get_thread_num

        type(AppConfig),  intent(in)         :: cfg
        real(kind=dp),    allocatable        :: Energies(:, :)
        
        type(HSDataFormat) :: HData, SData
        integer          :: ik, max_threads, thread_id

        ! Pools de matrizes temporárias (Hk e Sk) por thread
        complex(dp), allocatable :: Hk_pool(:, :, :)
        complex(dp), allocatable :: Sk_pool(:, :, :) ! <- Nova matriz de Overlap
        real(dp),    allocatable :: w_pool(:, :)
        
        ! Espaços de trabalho expandidos para a zhegv
        complex(dp), allocatable :: work_pool(:, :)
        real(dp),    allocatable :: rwork_pool(:, :)
        integer,     allocatable :: info_pool(:)
        integer                  :: lwork, lrwork

        ! 1. Carrega os dados (o arquivo agora deve conter Hr e Sr)
        call HData%ReadFile(cfg%HDataFile)
        call SData%ReadFile(cfg%SDataFile)

        max_threads = 1
        !$ max_threads = omp_get_max_threads()
        
        ! Regra do LAPACK zhegv para o tamanho mínimo de trabalho (jobz='N')
        lwork  = max(1, 2 * HData%NOrbs - 1)
        lrwork = max(1, 3 * HData%NOrbs - 2)

        associate( Hr        => HData%Data,              & ! H(r)
                   Sr        => SData%Data,              & ! S(r) - Novo dado do arquivo
                   nOrbs     => HData%NOrbs,             & 
                   nkPoints  => cfg%Bands%total_kpoints, &
                   FullPath  => cfg%Bands%FullPath,      &
                   Rpt       => HData%IndexList          & 
        )
            allocate(Energies(nOrbs, nkPoints))
            Energies = 0.0_dp

            ! Aloca os pools incluindo a matriz Sk e os novos tamanhos de work
            allocate(Hk_pool(nOrbs, nOrbs, max_threads))
            allocate(Sk_pool(nOrbs, nOrbs, max_threads)) ! Pool para S(k)
            allocate(w_pool(nOrbs, max_threads))
            allocate(work_pool(lwork, max_threads))
            allocate(rwork_pool(lrwork, max_threads))
            allocate(info_pool(max_threads))

            ! Loop sobre os pontos k (OpenMP)
            !$omp parallel do default(shared) private(ik, thread_id)
            do ik = 1, nkPoints
                thread_id = 0
                !$ thread_id = omp_get_thread_num()
                thread_id = thread_id + 1

                ! 2. Calcula as matrizes H(k) e S(k) para o ponto k atual
                call compute_Hk_at_kpoint(Hr, Rpt, Hk_pool(:, :, thread_id), FullPath(:, ik))
                call compute_Hk_at_kpoint(Sr, Rpt, Sk_pool(:, :, thread_id), FullPath(:, ik)) ! Reutiliza a função de Fourier para S(k)

                ! 3. Diagonalização Generalizada usando LAPACK zhegv
                ! '1' significa resolver o problema do tipo H*c = E*S*c
                ! 'N' calcula apenas autovalores (Energias)
                ! 'U' usa a parte triangular superior das matrizes
                call zhegv(1, 'N', 'U', nOrbs, Hk_pool(:, :, thread_id), nOrbs, & 
                            Sk_pool(:, :, thread_id), nOrbs, w_pool(:, thread_id), &
                            work_pool(:, thread_id), lwork, rwork_pool(:, thread_id), &
                            info_pool(thread_id))

                if (info_pool(thread_id) /= 0) then
                    print *, "Erro: zhegv falhou no ponto k:", ik, " Erro:", info_pool(thread_id)
                endif

                ! Armazena os autovalores generalizados (energias)
                Energies(:, ik) = w_pool(:, thread_id)
            end do
            !$omp end parallel do
            
            deallocate(Hk_pool, Sk_pool, w_pool, work_pool, rwork_pool, info_pool)
        end associate
    end function calc_band_structure_non_orth


    pure subroutine compute_Hk_at_kpoint(Hr, Rpt, Hk_k, kpt)
        use types

        ! Garante contiguidade de memória para o compilador vetorizar com segurança
        complex(dp), intent(in),  contiguous :: Hr(:, :, :)
        complex(dp), intent(out), contiguous :: Hk_k(:, :)
        real(dp),    intent(in)              :: kpt(3)
        integer,     intent(in),  contiguous :: Rpt(:, :)

        ! Constante 2*pi: converte k (coordenadas reduzidas/fracionárias da
        ! rede recíproca) e R (múltiplos inteiros dos vetores de rede direta)
        ! em uma fase física e-i.k.R, já que b_i . a_j = 2*pi * delta_ij.
        real(dp), parameter :: TWO_PI = 6.283185307179586476925286766559_dp

        integer     :: nNeighbors, ir
        real(dp)    :: k_dot_R
        complex(dp) :: phase_factor

        ! Inicializa a matriz deste ponto k
        Hk_k = cmplx(0.0d0, 0.0d0, kind=dp)
        
        ! Obtém o número total de vizinhos / vetores R a partir das dimensões de Hr
        nNeighbors = size(Hr, 3)

        ! Exemplo do loop de Fourier espacial ultra-eficiente
        do ir = 1, nNeighbors
            ! k . R (kpt em coordenadas reduzidas, Rpt em unidades inteiras de
            ! vetores de rede direta -> precisa do fator 2*pi para virar fase)
            k_dot_R = TWO_PI * dot_product(kpt, real(Rpt(:, ir), kind=dp))
            phase_factor = exp(cmplx(0.0d0, k_dot_R,    kind=dp))

            ! Auxilia no uso de operações vetorizadas e contíguas de memória
            Hk_k = Hk_k + Hr(:, :, ir) * phase_factor
        end do
    end subroutine compute_Hk_at_kpoint
end module core