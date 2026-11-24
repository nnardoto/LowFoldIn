! paoeval_bin <in.bin> <out.bin> [pao]: in = int32 n + n*3 real64 (bohr); out = n*29 real64 (pao_eval, s3p3d2f1, rcut 7)
! pao: arquivo .pao do Te (padrão test/data/pao/Te7.0.pao, relativo à raiz do repositório)
program paoeval_bin
    use types, only: dp
    use pao_basis
    implicit none
    type(pao_species_t) :: sp
    integer :: ierr, n, i, u
    character(len=:), allocatable :: msg
    character(len=512) :: fi, fo, pao
    real(dp), allocatable :: d(:, :), phi(:, :)
    call get_command_argument(1, fi); call get_command_argument(2, fo)
    pao = 'test/data/pao/Te7.0.pao'
    if (command_argument_count() >= 3) call get_command_argument(3, pao)
    call read_pao_species(trim(pao), &
                          's3p3d2f1', 7.0_dp, sp, ierr, msg)
    if (ierr /= 0) error stop 1
    open(newunit=u, file=trim(fi), access='stream', form='unformatted', status='old')
    read(u) n; allocate(d(3, n), phi(29, n)); read(u) d; close(u)
    !$omp parallel do
    do i = 1, n
        call pao_eval(sp, d(:, i), phi(:, i))
    end do
    open(newunit=u, file=trim(fo), access='stream', form='unformatted', status='replace')
    write(u) phi; close(u)
end program
