! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> spectral: energias de banda e pesos de projeção por região ao longo de um caminho k
!> (contracts §3, §5). Driver fino em torno do módulo spectral_engine: este arquivo
!> trata apenas a CLI, a saída de progresso/resumo e a escrita dos arquivos de resultado.
!>
!> Usage:
!>   spectral <spectral.nml>
!>
!> Entrada: um arquivo namelist (veja docs/spectral.md para todas as opções), que por sua
!>          vez aponta para os dados HSIO, a estrutura, o caminho k e as definições de região.
!> Saída  : <output_prefix>.spectral.bin (sempre) e <output_prefix>.spectral.txt
!>          (somente se write_text = .true.). Estes são a entrada do specmap.
!> Código : 0 em caso de sucesso; 1 em uso incorreto ou qualquer erro (mensagem com
!>          prefixo "[spectral error]" no stderr).
!>
!> Documentação de usuário: docs/spectral.md
program spectral
    use, intrinsic :: iso_fortran_env, only: int64, error_unit
    use types, only: dp
    use spectral_engine, only: spectral_input_t, spectral_timing_t, read_spectral_input, run_spectral
    use spectral_data, only: spectral_data_t, write_spectral_file, write_spectral_text, &
                             SPIN_COLLINEAR
    use terminal_ui, only: print_box, print_box_close
    !$ use omp_lib, only: omp_get_max_threads
    implicit none

    type(spectral_input_t)        :: inp
    type(spectral_data_t)         :: sd
    type(spectral_timing_t)       :: timing
    character(len=1024)           :: nml_file, errmsg
    character(len=:), allocatable :: bin_file, txt_file
    integer                       :: ierr, ir, max_thr
    integer(int64)                :: c0, c1, crate

    ! 1. Leitura dos argumentos
    if (command_argument_count() /= 1) then
        write(error_unit, '(A)') 'Usage: spectral <spectral.nml>'
        write(error_unit, '(A)') '  Computes band energies and Mulliken or Lowdin weights per region along a k-path'
        write(error_unit, '(A)') '  and writes <output_prefix>.spectral.bin (+ .spectral.txt if write_text).'
        stop 1
    end if
    call get_command_argument(1, nml_file)

    ! 2. Banner de início: informação de threads (o código nunca define variáveis de ambiente).
    max_thr = 1
    !$ max_thr = omp_get_max_threads()
    call print_box('spectral | projected band calculation')
    print '(A,A)',  '  Input file       : ', trim(nml_file)
    print '(A,I0)', '  OpenMP threads   : ', max_thr
    if (max_thr > 1) then
        print '(A)', '  Note: keep BLAS/LAPACK single-threaded to avoid oversubscription.'
    end if
    call print_box_close()
    print '(A)', ''

    ! 3. Leitura e validação da entrada namelist
    call read_spectral_input(trim(nml_file), inp, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)

    ! 4. Cálculo das energias e dos pesos de projeção (laço em k, paralelo em k)
    call run_spectral(inp, sd, timing, ierr, errmsg, show_progress=.true.)
    if (ierr /= 0) call die(errmsg)

    ! 5. Escrita dos arquivos de resultado (o tempo de escrita é somado ao total do motor)
    bin_file = trim(inp%output_prefix) // '.spectral.bin'
    txt_file = trim(inp%output_prefix) // '.spectral.txt'
    print '(A)', '[output] Writing result files'
    call system_clock(c0, crate)
    call write_spectral_file(bin_file, sd, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    if (inp%write_text) then
        call write_spectral_text(txt_file, sd, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
    end if
    call system_clock(c1)
    timing%write_s = real(c1 - c0, dp) / real(max(crate, 1_int64), dp)
    timing%total_s = timing%total_s + timing%write_s

    ! 6. Resumo
    print '(A)', ''
    call print_box('spectral | summary')
    print '(A,I0,A,I0,A,I0,A,I0)', '  nk = ', sd%nk, '   nband = ', sd%nband, &
        '   nspin = ', sd%nspin, '   nregion = ', sd%nregion
    print '(A,A)', '  spin mode        : ', trim(inp%spin_mode)
    print '(A,A)', '  method           : ', trim(inp%method)
    print '(A,A)', '  energy reference : ', eref_text()
    print '(A,A)', '  fingerprint      : ', sd%fingerprint
    print '(A)',   '  regions          :'
    do ir = 1, sd%nregion
        print '(A,I3,2X,A,A,I0,A)', '    ', ir, trim(sd%region_name(ir)), ' (', &
            sd%region_natoms(ir), ' atoms)'
    end do
    print '(A)', '  files written    :'
    print '(A,A)', '    ', bin_file
    if (inp%write_text) print '(A,A)', '    ', txt_file
    print '(A)',   '  timing (measured with system_clock, s):'
    print '(A,I0)',    '    threads (k-loop)                ', timing%nthreads
    print '(A,F12.4)', '    read                    [wall]  ', timing%read_s
    print '(A,F12.4)', '    k-loop                  [wall]  ', timing%loop_wall_s
    print '(A,F12.4)', '      build H(k),S(k) [accumulated] ', timing%build_s
    print '(A,F12.4)', '      diagonalization [accumulated] ', timing%diag_s
    print '(A,F12.4)', '      projection      [accumulated] ', timing%proj_s
    print '(A,F12.4)', '    write                   [wall]  ', timing%write_s
    print '(A,F12.4)', '    total                   [wall]  ', timing%total_s
    print '(A)',   '  [accumulated] = sum over threads of per-thread wall time (thread-seconds);'
    print '(A)',   '  it can exceed the k-loop wall time when more than one thread is used.'
    if (sd%spin_mode == SPIN_COLLINEAR) print '(A)', '  (spin channel 1 = up, 2 = down)'
    call print_box_close()

contains

    !> Descrição legível da referência de energia para a caixa de resumo.
    !> 'fermi' significa que o motor já entrega E - E_F; 'absolute' mantém as
    !> energias brutas e E_F (se conhecido) é apenas guardado como metadado.
    function eref_text() result(s)
        character(len=:), allocatable :: s
        character(len=32) :: buf
        select case (trim(inp%energy_reference))
        case ('fermi')
            s = 'fermi: E - E_F (eV), E_F = 0 by definition (energies not shifted by spectral)'
        case ('absolute')
            if (sd%fermi_known) then
                write(buf, '(F12.6)') sd%fermi_energy_ev
                s = 'absolute: E (eV), E_F = ' // trim(adjustl(buf)) // ' eV (metadata only)'
            else
                s = 'absolute: E (eV), E_F unknown'
            end if
        case default
            s = 'unknown: E (eV, reference unknown)'
        end select
    end function eref_text

    !> Imprime a mensagem no stderr e encerra com código 1.
    subroutine die(msg)
        character(len=*), intent(in) :: msg
        write(error_unit, '(A,A)') '[spectral error] ', trim(msg)
        stop 1
    end subroutine die

end program spectral
