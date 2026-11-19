! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> DOS de superficie semi-infinita a partir de uma Hamiltoniana `.hsio`
!> ortonormal. A saida e uma tabela simples (k1, k2, E, DOS).
program surface_dos
    use types, only: dp
    use hs_io, only: HSDataFormat
    use surface_green, only: build_principal_layer, surface_greens_function, trace_surface_dos
    implicit none

    integer, parameter :: PATH_LEN = 1024
    character(len=PATH_LEN) :: h_file, output_file, input_file
    character(len=8) :: side
    integer :: normal_axis, nk1, nk2, energy_num, max_iter
    real(dp) :: k_origin(3), k_span1(3), k_span2(3)
    real(dp) :: energy_min, energy_max, eta, tolerance
    namelist /SURFACE_INPUT/ h_file, output_file, side, normal_axis, nk1, nk2, &
        k_origin, k_span1, k_span2, energy_min, energy_max, energy_num, eta, tolerance, max_iter

    type(HSDataFormat) :: hr
    complex(dp), allocatable :: h00(:, :), h01(:, :), gs(:, :)
    character(len=:), allocatable :: errmsg
    character(len=512) :: parallel_errmsg
    real(dp), allocatable :: dos_grid(:, :, :)
    integer, allocatable :: iter_grid(:, :, :)
    real(dp) :: u1, u2, energy, kfrac(3), dos
    integer :: iu, ios, ik, ik1, ik2, ie, ierr, thickness, thickness_k, niter, out
    integer :: parallel_error

    if (command_argument_count() /= 1) then
        write(*, '(A)') 'Usage: surface_dos INPUT.nml'
        stop 1
    end if
    call get_command_argument(1, input_file)

    h_file = ''
    output_file = 'surface_dos.dat'
    side = 'left'
    normal_axis = 3
    nk1 = 51
    nk2 = 51
    k_origin = [-0.5_dp, -0.5_dp, 0.0_dp]
    k_span1 = [1.0_dp, 0.0_dp, 0.0_dp]
    k_span2 = [0.0_dp, 1.0_dp, 0.0_dp]
    energy_min = -1.0_dp
    energy_max = 1.0_dp
    energy_num = 201
    eta = 0.01_dp
    tolerance = 1.0e-10_dp
    max_iter = 100

    open(newunit=iu, file=trim(input_file), status='old', action='read', iostat=ios)
    if (ios /= 0) call fatal("cannot open input file '" // trim(input_file) // "'")
    read(iu, nml=SURFACE_INPUT, iostat=ios)
    close(iu)
    if (ios /= 0) call fatal('cannot parse &SURFACE_INPUT namelist')
    if (len_trim(h_file) == 0) call fatal('h_file is required')
    if (nk1 < 1 .or. nk2 < 1 .or. energy_num < 1) call fatal('grid dimensions must be positive')
    if (eta <= 0.0_dp .or. tolerance <= 0.0_dp .or. max_iter < 1) &
        call fatal('eta, tolerance and max_iter must be positive')
    if (normal_axis < 1 .or. normal_axis > 3) call fatal('normal_axis must be 1, 2 or 3')

    call hr%ReadFile(trim(h_file))

    ! Validate the geometry once and determine the principal-layer thickness.
    kfrac = k_origin
    call build_principal_layer(hr%Data, hr%IndexList, kfrac, normal_axis, &
                               h00, h01, thickness, ierr, errmsg)
    if (ierr /= 0) call fatal(errmsg)
    deallocate(h00, h01)

    allocate(dos_grid(energy_num, nk1, nk2), iter_grid(energy_num, nk1, nk2))
    dos_grid = 0.0_dp
    iter_grid = 0
    parallel_error = 0
    parallel_errmsg = ''

    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp& private(ik, ik1, ik2, u1, u2, kfrac, h00, h01, gs, thickness_k, &
    !$omp&         ierr, errmsg, ie, energy, niter, dos)
    do ik = 1, nk1 * nk2
        ik1 = modulo(ik - 1, nk1) + 1
        ik2 = (ik - 1) / nk1 + 1
        u1 = grid_fraction(ik1, nk1)
        u2 = grid_fraction(ik2, nk2)
        kfrac = k_origin + u1 * k_span1 + u2 * k_span2
        call build_principal_layer(hr%Data, hr%IndexList, kfrac, normal_axis, &
                                   h00, h01, thickness_k, ierr, errmsg)
        if (ierr == 0) then
            allocate(gs(size(h00, 1), size(h00, 2)))
            do ie = 1, energy_num
                energy = energy_min + grid_fraction(ie, energy_num) * (energy_max - energy_min)
                call surface_greens_function(h00, h01, energy, eta, trim(side), tolerance, &
                                             max_iter, gs, niter, ierr, errmsg)
                if (ierr /= 0) exit
                dos = trace_surface_dos(gs)
                dos_grid(ie, ik1, ik2) = dos
                iter_grid(ie, ik1, ik2) = niter
            end do
            deallocate(gs)
        end if
        if (allocated(h00)) deallocate(h00)
        if (allocated(h01)) deallocate(h01)
        if (ierr /= 0) then
            !$omp critical (surface_dos_error)
            if (parallel_error == 0) then
                parallel_error = ierr
                parallel_errmsg = errmsg
            end if
            !$omp end critical (surface_dos_error)
        end if
    end do
    !$omp end parallel do

    if (parallel_error /= 0) call fatal(trim(parallel_errmsg))

    open(newunit=out, file=trim(output_file), status='replace', action='write', iostat=ios)
    if (ios /= 0) call fatal("cannot open output file '" // trim(output_file) // "'")
    write(out, '(A)') '# LowFoldIn semi-infinite surface DOS'
    write(out, '(A,A)') '# h_file = ', trim(h_file)
    write(out, '(A,A)') '# side = ', trim(side)
    write(out, '(A,I0)') '# normal_axis = ', normal_axis
    write(out, '(A,ES16.8)') '# eta_eV = ', eta
    write(out, '(A)') '# k1_fraction  k2_fraction  energy_eV  dos_states_per_eV  iterations'

    do ik2 = 1, nk2
        u2 = grid_fraction(ik2, nk2)
        do ik1 = 1, nk1
            u1 = grid_fraction(ik1, nk1)
            do ie = 1, energy_num
                energy = energy_min + grid_fraction(ie, energy_num) * (energy_max - energy_min)
                write(out, '(4(ES22.13,1X),I0)') u1, u2, energy, &
                    dos_grid(ie, ik1, ik2), iter_grid(ie, ik1, ik2)
            end do
            write(out, '()')
        end do
    end do
    close(out)

    write(*, '(A,I0)') 'Principal-layer thickness (unit cells): ', thickness
    write(*, '(A)') 'Surface DOS written to: ' // trim(output_file)

contains

    pure real(dp) function grid_fraction(index, count) result(value)
        integer, intent(in) :: index, count
        if (count == 1) then
            value = 0.0_dp
        else
            value = real(index - 1, dp) / real(count - 1, dp)
        end if
    end function grid_fraction

    subroutine fatal(message)
        character(len=*), intent(in) :: message
        write(*, '(A)') 'Erro: ' // trim(message)
        stop 1
    end subroutine fatal
end program surface_dos
