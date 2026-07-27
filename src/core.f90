module core
    use types
    use config
    use terminal_ui
    
    private
    public :: calc_band_structure, compute_Hk_at_kpoint

    ! Interface para a rotina de diagonalização do LAPACK (Double Precision Complex)
    interface
        subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
            use iso_fortran_env, only: dp => real64
            character          :: jobz, uplo
            integer            :: n, lda, lwork, info
            complex(dp)        :: a(lda, *), work(*)
            real(dp)           :: w(*), rwork(*)
        end subroutine zheev
    end interface


contains
    
    function calc_band_structure(cfg) result(Energies)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        use omp_lib, only: omp_get_max_threads, omp_get_thread_num

        type(AppConfig),  intent(in)         :: cfg
        real(kind=dp),    allocatable        :: Energies(:, :)
        
        type(FileFormat) :: FileData
        integer          :: ik

        ! Variaveis para Gerenciar a Pré-alocação alinhado com openmp
        integer :: max_threads, thread_id

        ! Matrizes de Trabalho para LAPACK com Dimensão extra para openmp
        complex(dp), allocatable :: Hk_pool(:, :, :)
        real(dp),    allocatable :: w_pool(:, :)
        complex(dp), allocatable :: work_pool(:, :)
        real(dp),    allocatable :: rwork_pool(:, :)
        integer,     allocatable :: info_pool(:)
        integer                  :: lwork

        max_threads = 1
        !$ max_threads = omp_get_max_threads()
        
        ! Tamanho de lwork = 3 * NOrbs - 1, conforme exigido pelo LAPACK
        lwork = 3 * FileData%NOrbs - 1

        ! 1. Carrega os dados do arquivo
        call FileData%ReadFile(cfg%DataFile)

        associate( Hr        => FileData%Data,           &
                   nOrbs     => FileData%NOrbs,          & 
                   nkPoints  => cfg%Bands%total_kpoints, &
                   FullPath  => cfg%Bands%FullPath,      &
                   Rpt       => FileData%IndexList       & ! Rpt é a lista de vetores de vizinhos (h, k, l)
        )
            ! Aloca a matriz completa de uma só vez
            allocate(Energies(nOrbs, nkPoints))
            Energies = 0.0d0

            ! Aloca as matrizes de trabalho para cada thread uma unica vez
            allocate(Hk_pool(nOrbs, nOrbs, max_threads))
            allocate(w_pool(nOrbs, max_threads))
            allocate(work_pool(3*nOrbs - 1, max_threads))
            allocate(rwork_pool(3*nOrbs - 2, max_threads))
            allocate(info_pool(max_threads))

            ! Loop sobre os pontos k
            !$omp parallel do default(shared) private(ik, thread_id)
            do ik = 1, nkPoints
                thread_id = 0
                !$ thread_id = omp_get_thread_num()

                ! Ajusta o indice para fortran (1-based) para acessar a matriz de trabalho da thread correta
                thread_id = thread_id + 1

                ! 2. Calcula a matriz H(k) para o ponto k atual
                call compute_Hk_at_kpoint(Hr, Rpt, Hk_pool(:, :, thread_id), FullPath(:, ik))

                ! 3. Diagonaliza H(k) usando LAPACK
                call zheev('N', 'U', nOrbs, Hk_pool(:, :, thread_id), nOrbs, & 
                            w_pool(:, thread_id), work_pool(:, thread_id),   &
                            lwork, rwork_pool(:, thread_id), info_pool(thread_id))

                if (info_pool(thread_id) /= 0) then
                    print *, "Erro: LAPACK falhou no ponto k:", ik, " Thread:", thread_id
                endif

            end do
            !$omp end parallel do
            deallocate(Hk_pool, w_pool, work_pool, rwork_pool, info_pool)
        end associate
    end function calc_band_structure

    pure subroutine compute_Hk_at_kpoint(Hr, Rpt, Hk_k, kpt)
        use, intrinsic :: iso_fortran_env, only: dp => real64

        ! Garante contiguidade de memória para o compilador vetorizar com segurança
        complex(dp), intent(in),  contiguous :: Hr(:, :, :)
        complex(dp), intent(out), contiguous :: Hk_k(:, :)
        real(dp),    intent(in)              :: kpt(3)
        integer,     intent(in),  contiguous :: Rpt(:, :)
        
        integer     :: nNeighbors, ir
        real(dp)    :: k_dot_R
        complex(dp) :: phase_factor

        ! Inicializa a matriz deste ponto k
        Hk_k = cmplx(0.0d0, 0.0d0, kind=dp)
        
        ! Obtém o número total de vizinhos / vetores R a partir das dimensões de Hr
        nNeighbors = size(Hr, 3)

        ! Exemplo do loop de Fourier espacial ultra-eficiente
        do ir = 1, nNeighbors
            ! k . R (Assumindo cálculo do produto escalar com vetor R)
            k_dot_R = dot_product(kpt, real(Rpt(:, ir), kind=dp))
            phase_factor = exp(cmplx(0.0d0, k_dot_R,    kind=dp))

            ! Auxilia no uso de operações vetorizadas e contíguas de memória
            Hk_k = Hk_k + Hr(:, :, ir) * phase_factor
        end do
    end subroutine compute_Hk_at_kpoint

end module core