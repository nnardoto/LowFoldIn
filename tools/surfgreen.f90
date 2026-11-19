! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> surfgreen: funções de Green de superfícies, filmes e junções em base
!> não ortogonal (H(R), S(R) do OpenMX ou de Wannier90).
!>
!> Usage:
!>   surfgreen INPUT.nml                    cálculo; grava <prefix>.green.bin e <prefix>_central.nml
!>   surfgreen --write-structure INPUT.nml  grava só <prefix>_central.nml (para o atomselect)
!>   surfgreen --info FILE.green.bin        resumo do arquivo
!>   surfgreen --slice FILE.green.bin [--map total|REGION|coh|incoh] [--iq N]
!>                     [--energy E] [--out PREFIX]
!>       exporta um corte em texto compatível com tools/gnuplot:
!>       caminho k -> PREFIX.map.dat + PREFIX.kticks.gp (map_positive.gp);
!>       grade k   -> PREFIX.cemap.dat + PREFIX.cegrid.gp (cemap.gp), na energia
!>                    mais próxima de E.
!>
!> Documentação: docs/surfgreen.md
program surfgreen
    use, intrinsic :: iso_fortran_env, only: error_unit, int64
    use types, only: dp
    use green_engine, only: green_input_t, read_green_input, run_green
    use green_data, only: green_data_t, write_green_file, read_green_file, meta_get, GREEN_KPATH
    use terminal_ui, only: print_box, print_box_close
    implicit none

    character(len=1024) :: arg1, arg2
    integer :: nargs

    nargs = command_argument_count()
    if (nargs < 1) call usage()
    call get_command_argument(1, arg1)
    select case (trim(arg1))
    case ('--write-structure')
        if (nargs /= 2) call usage()
        call get_command_argument(2, arg2)
        call run(trim(arg2), .true.)
    case ('--info')
        if (nargs /= 2) call usage()
        call get_command_argument(2, arg2)
        call info(trim(arg2))
    case ('--slice')
        if (nargs < 2) call usage()
        call get_command_argument(2, arg2)
        call slice(trim(arg2))
    case ('-h', '--help')
        call usage()
    case default
        if (nargs /= 1) call usage()
        call run(trim(arg1), .false.)
    end select

contains

    subroutine usage()
        write(error_unit, '(A)') 'Usage: surfgreen INPUT.nml'
        write(error_unit, '(A)') '       surfgreen --write-structure INPUT.nml'
        write(error_unit, '(A)') '       surfgreen --info FILE.green.bin'
        write(error_unit, '(A)') '       surfgreen --slice FILE.green.bin [--map total|REGION|coh|incoh] [--iq N]'
        write(error_unit, '(A)') '                 [--energy E] [--out PREFIX]'
        stop 1
    end subroutine usage

    subroutine die(msg)
        character(len=*), intent(in) :: msg
        write(error_unit, '(A,A)') '[surfgreen error] ', trim(msg)
        stop 1
    end subroutine die

    subroutine run(nml, structure_only)
        character(len=*), intent(in) :: nml
        logical,          intent(in) :: structure_only
        type(green_input_t) :: inp
        type(green_data_t) :: gd
        character(len=:), allocatable :: errmsg
        integer :: ierr, nbad, nthr
        integer(int64) :: c0, c1, crate
        character(len=:), allocatable :: out

        call read_green_input(nml, inp, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        call system_clock(c0, crate)
        if (.not. structure_only) then
            call print_box('surfgreen | Green functions of surfaces and interfaces')
            print '(A,A)', '  Input file       : ', nml
            print '(A,A)', '  Mode             : ', trim(inp%mode)
            print '(A,A)', '  Projection       : ', trim(inp%projection)
            call print_box_close()
        end if
        call run_green(inp, gd, ierr, errmsg, structure_only=structure_only, nonconverged=nbad, nthreads=nthr)
        if (ierr /= 0) call die(errmsg)
        if (structure_only) then
            print '(A,A)', 'Central structure written: ', trim(inp%output_prefix) // '_central.nml'
            return
        end if
        out = trim(inp%output_prefix) // '.green.bin'
        call write_green_file(out, gd, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        call system_clock(c1)
        print '(A,I0,A,I0,A,I0,A)', '  Points           : ', gd%nk, ' k x ', gd%ne, ' E (', gd%nreg, ' regions)'
        print '(A,I0)',  '  OpenMP threads   : ', nthr
        print '(A,I0,A)', '  Decimation iters : max ', maxval(abs(gd%niter)), ' per point'
        print '(A,ES10.2)', '  Max residual     : ', maxval(gd%residual)
        print '(A,ES10.2)', '  min(-Im G_ii)    : ', minval(gd%min_neg_im_gii)
        print '(A,ES10.2)', '  min eig(S_C)     : ', minval(gd%smin)
        if (nbad > 0) print '(A,I0,A)', '  WARNING: decimation did not converge at ', nbad, &
            ' point(s) (negative niter in diagnostics); raise decim_maxit or eta_ev'
        if (minval(gd%min_neg_im_gii) < -1.0e-8_dp) &
            print '(A)', '  WARNING: negative -Im G_ii found (check eta, decimation and contact blocks)'
        if (len(meta_get(gd, 'warning')) > 0) print '(A,A)', '  WARNING: ', meta_get(gd, 'warning')
        print '(A,A)', '  Output           : ', out
        print '(A,A)', '  Central structure: ', trim(inp%output_prefix) // '_central.nml'
        print '(A,F10.2,A)', '  Wall time        : ', real(c1 - c0, dp) / real(crate, dp), ' s'
    end subroutine run

    subroutine info(fname)
        character(len=*), intent(in) :: fname
        type(green_data_t) :: gd
        character(len=:), allocatable :: errmsg
        integer :: ierr, i
        call read_green_file(fname, gd, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        print '(A,A)', 'File: ', fname
        do i = 1, gd%nmeta
            print '(2X,A)', trim(gd%meta(i))
        end do
        print '(A,A,A,I0,A,I0,A,I0)', '  k mode = ', merge('path', 'grid', gd%kmode == GREEN_KPATH), &
            ', nk = ', gd%nk, ', ne = ', gd%ne, ', nqz = ', gd%nq
        print '(A,F10.4,A,F10.4)', '  energy range = ', gd%energy(1), ' .. ', gd%energy(gd%ne)
        print '(A)', '  regions:'
        do i = 1, gd%nreg
            print '(4X,A)', trim(gd%region(i))
        end do
        print '(A,I0,A,ES10.2)', '  max decimation iterations = ', maxval(abs(gd%niter)), &
            ', max residual = ', maxval(gd%residual)
        print '(A,ES10.2,A,ES10.2)', '  min(-Im G_ii) = ', minval(gd%min_neg_im_gii), &
            ', min eig(S_C) = ', minval(gd%smin)
    end subroutine info

    subroutine slice(fname)
        character(len=*), intent(in) :: fname
        type(green_data_t) :: gd
        character(len=:), allocatable :: errmsg, eref
        character(len=1024) :: a, mapname, outp
        real(dp) :: energy
        logical :: have_energy
        integer :: ierr, i, iq, ie, ik, ir, u, j, i1
        real(dp), allocatable :: data(:, :)   ! (ne, nk)
        real(dp) :: eu(3), ev(3)

        mapname = 'total'
        outp = ''
        iq = 1
        have_energy = .false.
        energy = 0.0_dp
        i = 3
        do while (i <= command_argument_count())
            call get_command_argument(i, a)
            if (i + 1 > command_argument_count()) call usage()
            select case (trim(a))
            case ('--map');    call get_command_argument(i + 1, mapname)
            case ('--out');    call get_command_argument(i + 1, outp)
            case ('--iq')
                call get_command_argument(i + 1, a)
                read(a, *, iostat=ierr) iq
                if (ierr /= 0) call die('invalid --iq value')
            case ('--energy')
                call get_command_argument(i + 1, a)
                read(a, *, iostat=ierr) energy
                if (ierr /= 0) call die('invalid --energy value')
                have_energy = .true.
            case default
                call usage()
            end select
            i = i + 2
        end do
        call read_green_file(fname, gd, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        if (len_trim(outp) == 0) then
            j = index(fname, '.green.bin', back=.true.)
            if (j > 1) then
                outp = fname(1:j - 1) // '.' // trim(mapname)
            else
                outp = fname // '.' // trim(mapname)
            end if
        end if

        allocate(data(gd%ne, gd%nk))
        select case (trim(mapname))
        case ('total')
            data = gd%total
        case ('coh', 'incoh')
            if (.not. gd%has_coherent) call die('file has no coherent maps')
            if (iq < 1 .or. iq > gd%nq) call die('--iq out of range')
            if (trim(mapname) == 'coh') then
                data = gd%icoh(:, :, iq)
            else
                data = gd%iinc(:, :, iq)
            end if
        case default
            ir = 0
            do j = 1, gd%nreg
                if (trim(gd%region(j)) == trim(mapname)) ir = j
            end do
            if (ir == 0) call die("region '" // trim(mapname) // "' not found (see surfgreen --info)")
            data = gd%amap(:, :, ir)
        end select

        eref = meta_get(gd, 'energy_reference')
        if (gd%kmode == GREEN_KPATH) then
            open(newunit=u, file=trim(outp) // '.map.dat', status='replace', action='write')
            write(u, '(A)') '# surfgreen slice of ' // fname // ', map = ' // trim(mapname)
            write(u, '(A)') '# columns: k_distance energy_eV intensity'
            do ik = 1, gd%nk
                do ie = 1, gd%ne
                    write(u, '(3ES18.9)') gd%kdist(ik), gd%energy(ie), data(ie, ik)
                end do
                write(u, '(A)') ''
            end do
            close(u)
            open(newunit=u, file=trim(outp) // '.kticks.gp', status='replace', action='write')
            write(u, '(A)') '# gnuplot include written by surfgreen --slice'
            write(u, '(A,ES16.8,A)') 'set xrange [0:', gd%kdist(gd%nk), ']'
            if (gd%nvertex > 0) then
                write(u, '(A)', advance='no') 'set xtics ('
                do j = 1, gd%nvertex
                    if (j > 1) write(u, '(A)', advance='no') ', '
                    write(u, '(A,A,A,ES16.8)', advance='no') '"', trim(gd%label(j)), '" ', &
                        gd%kdist(gd%vertex_k(j))
                end do
                write(u, '(A)') ')'
                write(u, '(A)') 'if (!exists("VLINE_COLOR")) VLINE_COLOR = "white"'
                do j = 2, gd%nvertex - 1
                    write(u, '(A,ES16.8,A,ES16.8,A)') 'set arrow from ', gd%kdist(gd%vertex_k(j)), &
                        ', graph 0 to ', gd%kdist(gd%vertex_k(j)), ', graph 1 nohead dt 2 lw 1 lc rgb VLINE_COLOR front'
                end do
            end if
            call write_eref(u, eref)
            close(u)
            print '(A)', 'Wrote ' // trim(outp) // '.map.dat and ' // trim(outp) // '.kticks.gp'
        else
            ie = 1
            if (have_energy) ie = minloc(abs(gd%energy - energy), 1)
            if (.not. have_energy .and. gd%ne > 1) call die('k-grid file with several energies: give --energy')
            ! base no plano (e_u || primeiro passo da grade)
            eu = gd%kcart(:, min(2, gd%nk)) - gd%kcart(:, 1)
            if (norm2(eu) > 0.0_dp) eu = eu / norm2(eu)
            i1 = min(1 + gd%n1, gd%nk)
            ev = gd%kcart(:, i1) - gd%kcart(:, 1)
            ev = ev - dot_product(ev, eu) * eu
            if (norm2(ev) > 0.0_dp) ev = ev / norm2(ev)
            open(newunit=u, file=trim(outp) // '.cemap.dat', status='replace', action='write')
            write(u, '(A)') '# surfgreen slice of ' // fname // ', map = ' // trim(mapname)
            write(u, '(A,F12.6,A)') '# energy = ', gd%energy(ie), ' eV'
            write(u, '(A)') '# columns: kx_invA ky_invA kz_invA k_u_invA k_v_invA intensity'
            do j = 1, gd%n2
                do i = 1, gd%n1
                    ik = i + (j - 1) * gd%n1
                    write(u, '(6ES18.9)') gd%kcart(:, ik), dot_product(gd%kcart(:, ik), eu), &
                        dot_product(gd%kcart(:, ik), ev), data(ie, ik)
                end do
                write(u, '(A)') ''
            end do
            close(u)
            open(newunit=u, file=trim(outp) // '.cegrid.gp', status='replace', action='write')
            write(u, '(A)') '# gnuplot include for constant-energy maps (written by surfgreen --slice)'
            if (maxval(abs(gd%kcart(3, :))) < 1.0e-10_dp) then
                write(u, '(A)') 'COLS = "1:2:6"'
                write(u, '(A,ES16.8,A,ES16.8,A)') 'set xrange [', minval(gd%kcart(1, :)), ':', maxval(gd%kcart(1, :)), ']'
                write(u, '(A,ES16.8,A,ES16.8,A)') 'set yrange [', minval(gd%kcart(2, :)), ':', maxval(gd%kcart(2, :)), ']'
                write(u, '(A)') 'set xlabel "k_x (1/Angstrom)"'
                write(u, '(A)') 'set ylabel "k_y (1/Angstrom)"'
            else
                write(u, '(A)') 'COLS = "4:5:6"'
                write(u, '(A)') 'set xlabel "k_u (1/Angstrom)"'
                write(u, '(A)') 'set ylabel "k_v (1/Angstrom)"'
            end if
            write(u, '(A)') 'set size ratio -1'
            call write_eref(u, eref)
            close(u)
            print '(A,F10.4,A)', 'Wrote ' // trim(outp) // '.cemap.dat (E = ', gd%energy(ie), &
                ' eV) and ' // trim(outp) // '.cegrid.gp'
        end if
    end subroutine slice

    subroutine write_eref(u, eref)
        integer, intent(in) :: u
        character(len=*), intent(in) :: eref
        select case (eref)
        case ('E-E_F')
            write(u, '(A)') 'EF = 0.0'
            write(u, '(A)') 'ELABEL = "E - E_F (eV)"'
        case ('absolute')
            write(u, '(A)') 'EF = NaN'
            write(u, '(A)') 'ELABEL = "E (eV)"'
        case default
            write(u, '(A)') 'EF = NaN'
            write(u, '(A)') 'ELABEL = "E (eV, reference unknown)"'
        end select
    end subroutine write_eref

end program surfgreen
