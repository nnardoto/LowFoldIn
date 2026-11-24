! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> sktrench: trincheira (defeito de linha) num slab de base ortogonal,
!> embutida entre dois meios semi-infinitos do próprio slab (módulo
!> trench_green). Entrada: namelists &TRENCH, &KZ e &WINDOWS (ver
!> docs/sktrench.md). Saídas: PREFIX.green.bin (A por átomo), PREFIX.edges.dat,
!> PREFIX.atoms.dat (pesos por janela), PREFIX.rho_onsite.dat,
!> PREFIX.rho_pairs.dat, PREFIX.central.structure.nml e, com complex_band,
!> PREFIX.lambda.dat.
!>
!> Uso: sktrench entrada.nml
program sktrench_main
    use, intrinsic :: iso_fortran_env, only: error_unit
    use trench_green, only: trench_input_t, trench_result_t, trench_read_input, trench_run, trench_write_outputs
    !$ use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none
    type(trench_input_t) :: inp
    type(trench_result_t) :: res
    character(len=512) :: nml
    character(len=:), allocatable :: errmsg
    integer :: ierr, i, nthr
    real :: t0, t1
    double precision :: w0, w1

    if (command_argument_count() < 1) then
        write(error_unit, '(a)') 'uso: sktrench entrada.nml'
        stop 1
    end if
    call get_command_argument(1, nml)
    call trench_read_input(trim(nml), inp, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    nthr = 1
    w0 = 0.0d0
    !$ nthr = omp_get_max_threads()
    !$ w0 = omp_get_wtime()
    call cpu_time(t0)
    call trench_run(inp, res, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    call trench_write_outputs(inp, res, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    call cpu_time(t1)
    w1 = w0
    !$ w1 = omp_get_wtime()
    print '(a)', 'sktrench: ' // trim(inp%output_prefix)
    do i = 1, size(res%info)
        print '(2x,a)', trim(res%info(i))
    end do
    print '(a,i0,a,f10.2,a,f10.2,a)', '  threads = ', nthr, ', wall = ', w1 - w0, ' s, cpu = ', t1 - t0, ' s'

contains

    subroutine die(msg)
        character(len=*), intent(in) :: msg
        write(error_unit, '(a)') 'sktrench: ' // trim(msg)
        stop 1
    end subroutine die

end program sktrench_main
