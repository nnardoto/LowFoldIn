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
!>   surfgreen --kperp-slice FILE.green.bin [--map coh|incoh] [--ik N] [--out PREFIX]
!>   surfgreen --photon-slice FILE.green.bin --photon-energy HV --work-function PHI
!>                     --inner-potential V0 [--energy-reference auto|absolute|E-E_F]
!>                     [--fermi-energy EF] [--map coh|incoh] [--ik N] [--out PREFIX]
!>       exporta um corte em texto compatível com tools/gnuplot:
!>       caminho k -> PREFIX.map.dat + PREFIX.kticks.gp (map_positive.gp);
!>       grade k   -> PREFIX.cemap.dat + PREFIX.cegrid.gp (cemap.gp), na energia
!>                    mais próxima de E.
!>
!> Documentação: docs/surfgreen.md
program surfgreen
    use, intrinsic :: iso_fortran_env, only: error_unit, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use types, only: dp
    use arpes_conversion, only: photon_to_kperp
    use green_engine, only: green_input_t, read_green_input, run_green
    use extracted_green, only: extracted_input_t, input_is_extracted, read_extracted_input, run_extracted
    use scc, only: scc_input_t, input_has_scc, read_scc_input, run_scc
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
    case ('--kperp-slice')
        if (nargs < 2) call usage()
        call get_command_argument(2, arg2)
        call kperp_slice(trim(arg2))
    case ('--photon-slice')
        if (nargs < 2) call usage()
        call get_command_argument(2, arg2)
        call photon_slice(trim(arg2))
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
        write(error_unit, '(A)') '       surfgreen --kperp-slice FILE.green.bin [--map coh|incoh] [--ik N] [--out PREFIX]'
        write(error_unit, '(A)') '       surfgreen --photon-slice FILE.green.bin --photon-energy HV --work-function PHI'
        write(error_unit, '(A)') '                 --inner-potential V0 [--energy-reference auto|absolute|E-E_F]'
        write(error_unit, '(A)') '                 [--fermi-energy EF] [--map coh|incoh] [--ik N] [--out PREFIX]'
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
        type(green_data_t), allocatable :: gs(:)
        character(len=:), allocatable :: errmsg
        integer :: ierr, nbad, nthr, is, u
        character(len=16) :: tag
        integer(int64) :: c0, c1, crate
        character(len=:), allocatable :: out

        if (input_is_extracted(nml)) then
            call run_x(nml, structure_only)
            return
        end if
        call read_green_input(nml, inp, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        call system_clock(c0, crate)
        if (.not. structure_only) then
            call print_box('surfgreen | Green functions of surfaces and interfaces')
            print '(A,A)', '  Input file       : ', nml
            print '(A,A)', '  Mode             : ', trim(inp%mode)
            print '(A,A)', '  Projection       : ', trim(inp%projection)
            if (inp%have_film) print '(A,A)', '  Contact          : ', trim(inp%contact_model)
            if (inp%nscan > 0) print '(A,A,A,I0,A)', '  Scan             : ', trim(inp%scan_parameter), ' (', &
                inp%nscan, ' values)'
            call print_box_close()
        end if
        if (inp%nscan > 0 .and. .not. structure_only) then
            call run_green(inp, gd, ierr, errmsg, nonconverged=nbad, nthreads=nthr, gscan=gs)
        else
            call run_green(inp, gd, ierr, errmsg, structure_only=structure_only, nonconverged=nbad, nthreads=nthr)
        end if
        if (ierr /= 0) call die(errmsg)
        if (structure_only) then
            print '(A,A)', 'Central structure written: ', trim(inp%output_prefix) // '_central.nml'
            return
        end if
        if (inp%nscan > 0) then
            ! um arquivo por valor: PREFIX_sNNN.green.bin, e o índice PREFIX.scan.txt
            open(newunit=u, file=trim(inp%output_prefix) // '.scan.txt', status='replace', action='write')
            write(u, '(A)') '# surfgreen scan of ' // trim(inp%scan_parameter)
            write(u, '(A)') '# index  value  file'
            do is = 1, size(gs)
                write(tag, '(A,I3.3)') '_s', is
                out = trim(inp%output_prefix) // trim(tag) // '.green.bin'
                call write_green_file(out, gs(is), ierr, errmsg)
                if (ierr /= 0) call die(errmsg)
                write(u, '(I5,ES20.10,2X,A)') is, inp%scan_values(is), out
            end do
            close(u)
            out = trim(inp%output_prefix) // '_sNNN.green.bin (' // trim(inp%output_prefix) // '.scan.txt)'
        else
            out = trim(inp%output_prefix) // '.green.bin'
            call write_green_file(out, gd, ierr, errmsg)
            if (ierr /= 0) call die(errmsg)
        end if
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

    !> Modo 'extracted': região central tirada de uma supercélula DFT.
    subroutine run_x(nml, structure_only)
        character(len=*), intent(in) :: nml
        logical,          intent(in) :: structure_only
        type(extracted_input_t) :: inp
        type(green_data_t) :: gd
        character(len=:), allocatable :: errmsg, out
        character(len=256), allocatable :: notes(:)
        integer :: ierr, nbad, nthr, i
        integer(int64) :: c0, c1, crate
        type(scc_input_t) :: si
        logical :: do_scc, conv

        call read_extracted_input(nml, inp, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        call system_clock(c0, crate)
        do_scc = input_has_scc(nml) .and. .not. structure_only
        if (do_scc) then
            ! &SCC: ciclo eletrostático autoconsistente (mod_scc) e mapa com o V convergido
            call read_scc_input(nml, si, ierr, errmsg)
            if (ierr /= 0) call die(errmsg)
            call print_box('surfgreen | extracted + self-consistent electrostatics (&SCC)')
            print '(A,A)', '  Input file       : ', nml
            call print_box_close()
            call run_scc(inp, si, gd, ierr, errmsg, notes, conv)
            if (allocated(notes)) then
                do i = 1, size(notes)
                    print '(2X,A)', trim(notes(i))
                end do
            end if
            if (ierr /= 0) call die(errmsg)
            print '(A,L1)', '  SCC converged    : ', conv
            print '(A,A)', '  SCC log/profile  : ', trim(inp%output_prefix) // '.scc_log.dat / .scc_profile.dat'
            if (.not. si%final_map) return
            out = trim(inp%output_prefix) // '.green.bin'
            call write_green_file(out, gd, ierr, errmsg)
            if (ierr /= 0) call die(errmsg)
            call system_clock(c1)
            print '(A,A)', '  Output           : ', out
            print '(A,F10.2,A)', '  Wall time        : ', real(c1 - c0, dp) / real(crate, dp), ' s'
            return
        end if
        if (.not. structure_only) then
            call print_box('surfgreen | central region extracted from a DFT supercell')
            print '(A,A)', '  Input file       : ', nml
            print '(A,A)', '  Mode             : extracted, leads ', trim(inp%left) // ' / ' // trim(inp%right)
            print '(A,A)', '  Projection       : ', trim(inp%projection)
            call print_box_close()
        end if
        call run_extracted(inp, gd, ierr, errmsg, structure_only=structure_only, nonconverged=nbad, &
                           nthreads=nthr, notes=notes)
        if (allocated(notes)) then
            do i = 1, size(notes)
                print '(2X,A)', trim(notes(i))
            end do
        end if
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
        if (inp%overlap_diagnostic) print '(A,ES10.2)', '  min eig(S_C)     : ', minval(gd%smin)
        if (nbad > 0) print '(A,I0,A)', '  WARNING: decimation did not converge at ', nbad, &
            ' point(s) (negative niter in diagnostics); raise decim_maxit or eta_ev'
        if (minval(gd%min_neg_im_gii) < -1.0e-8_dp) &
            print '(A)', '  WARNING: negative -Im G_ii found (check eta, decimation and seam blocks)'
        print '(A,A)', '  Output           : ', out
        print '(A,A)', '  Central structure: ', trim(inp%output_prefix) // '_central.nml'
        print '(A,F10.2,A)', '  Wall time        : ', real(c1 - c0, dp) / real(crate, dp), ' s'
    end subroutine run_x

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
        if (allocated(gd%kperp)) print '(A,F10.4,A,F10.4,A)', '  k_perp range = ', &
            gd%kperp(1), ' .. ', gd%kperp(gd%nq), ' 1/angstrom'
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
            if (allocated(gd%kperp) .and. (trim(mapname) == 'coh' .or. trim(mapname) == 'incoh')) &
                write(u, '(A,ES16.8,A)') '# k_perp = ', gd%kperp(iq), ' 1/angstrom'
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
            if (allocated(gd%kperp) .and. (trim(mapname) == 'coh' .or. trim(mapname) == 'incoh')) &
                write(u, '(A,ES16.8,A)') '# k_perp = ', gd%kperp(iq), ' 1/angstrom'
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

    subroutine kperp_slice(fname)
        character(len=*), intent(in) :: fname
        type(green_data_t) :: gd
        character(len=:), allocatable :: errmsg
        character(len=1024) :: a, mapname, outp
        real(dp), allocatable :: data(:, :)
        integer :: ierr, i, ik, iq, ie, j, u

        mapname = 'coh'
        outp = ''
        ik = 1
        i = 3
        do while (i <= command_argument_count())
            call get_command_argument(i, a)
            if (i + 1 > command_argument_count()) call usage()
            select case (trim(a))
            case ('--map');    call get_command_argument(i + 1, mapname)
            case ('--out');    call get_command_argument(i + 1, outp)
            case ('--ik')
                call get_command_argument(i + 1, a)
                read(a, *, iostat=ierr) ik
                if (ierr /= 0) call die('invalid --ik value')
            case default
                call usage()
            end select
            i = i + 2
        end do
        call read_green_file(fname, gd, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        if (.not. allocated(gd%kperp)) call die('file has no ARPES k_perp axis')
        if (ik < 1 .or. ik > gd%nk) call die('--ik out of range')
        if (trim(mapname) /= 'coh' .and. trim(mapname) /= 'incoh') call die("--map must be 'coh' or 'incoh'")
        allocate(data(gd%ne, gd%nq))
        if (trim(mapname) == 'coh') then
            data = gd%icoh(:, ik, :)
        else
            data = gd%iinc(:, ik, :)
        end if
        if (len_trim(outp) == 0) then
            j = index(fname, '.green.bin', back=.true.)
            if (j > 1) then
                outp = fname(1:j - 1) // '.' // trim(mapname)
            else
                outp = fname // '.' // trim(mapname)
            end if
        end if
        open(newunit=u, file=trim(outp) // '.kperp.dat', status='replace', action='write')
        write(u, '(A)') '# surfgreen normal-emission ARPES proxy of ' // fname // ', map = ' // trim(mapname)
        write(u, '(A,I0)') '# k_point_index = ', ik
        write(u, '(A,3ES18.9)') '# k_cart_invA = ', gd%kcart(:, ik)
        write(u, '(A)') '# columns: k_perp_invA energy_eV intensity'
        do iq = 1, gd%nq
            do ie = 1, gd%ne
                write(u, '(3ES18.9)') gd%kperp(iq), gd%energy(ie), data(ie, iq)
            end do
            write(u, '(A)') ''
        end do
        close(u)
        print '(A)', 'Wrote ' // trim(outp) // '.kperp.dat'
    end subroutine kperp_slice

    subroutine photon_slice(fname)
        character(len=*), intent(in) :: fname
        type(green_data_t) :: gd
        character(len=:), allocatable :: errmsg, eref, file_eref
        character(len=1024) :: a, mapname, outp
        character(len=16) :: reference
        real(dp), allocatable :: data(:, :)
        real(dp), allocatable :: out_energy(:), out_kinetic(:), out_kperp(:), out_intensity(:)
        real(dp) :: photon_ev, work_function_ev, inner_potential_ev, fermi_energy_ev
        real(dp) :: kperp, kinetic, fraction
        logical :: have_photon, have_work_function, have_inner_potential, have_fermi
        integer :: ierr, i, ik, iq, ie, j, u, nvalid, nthreshold, noutside

        mapname = 'coh'
        outp = ''
        reference = 'auto'
        ik = 1
        photon_ev = 0.0_dp
        work_function_ev = 0.0_dp
        inner_potential_ev = 0.0_dp
        fermi_energy_ev = 0.0_dp
        have_photon = .false.
        have_work_function = .false.
        have_inner_potential = .false.
        have_fermi = .false.
        i = 3
        do while (i <= command_argument_count())
            call get_command_argument(i, a)
            if (i + 1 > command_argument_count()) call usage()
            select case (trim(a))
            case ('--map');    call get_command_argument(i + 1, mapname)
            case ('--out');    call get_command_argument(i + 1, outp)
            case ('--energy-reference')
                call get_command_argument(i + 1, a)
                reference = trim(a)
            case ('--ik')
                call get_command_argument(i + 1, a)
                read(a, *, iostat=ierr) ik
                if (ierr /= 0) call die('invalid --ik value')
            case ('--photon-energy')
                call get_command_argument(i + 1, a)
                read(a, *, iostat=ierr) photon_ev
                if (ierr /= 0) call die('invalid --photon-energy value')
                have_photon = .true.
            case ('--work-function')
                call get_command_argument(i + 1, a)
                read(a, *, iostat=ierr) work_function_ev
                if (ierr /= 0) call die('invalid --work-function value')
                have_work_function = .true.
            case ('--inner-potential')
                call get_command_argument(i + 1, a)
                read(a, *, iostat=ierr) inner_potential_ev
                if (ierr /= 0) call die('invalid --inner-potential value')
                have_inner_potential = .true.
            case ('--fermi-energy')
                call get_command_argument(i + 1, a)
                read(a, *, iostat=ierr) fermi_energy_ev
                if (ierr /= 0) call die('invalid --fermi-energy value')
                have_fermi = .true.
            case default
                call usage()
            end select
            i = i + 2
        end do
        if (.not. have_photon) call die('--photon-slice requires --photon-energy')
        if (.not. have_work_function) call die('--photon-slice requires --work-function')
        if (.not. have_inner_potential) call die('--photon-slice requires --inner-potential')
        if (.not. ieee_is_finite(photon_ev) .or. photon_ev <= 0.0_dp) call die('--photon-energy must be finite and > 0')
        if (.not. ieee_is_finite(work_function_ev) .or. work_function_ev < 0.0_dp) &
            call die('--work-function must be finite and >= 0')
        if (.not. ieee_is_finite(inner_potential_ev) .or. inner_potential_ev < 0.0_dp) &
            call die('--inner-potential must be finite and >= 0 (positive below vacuum)')

        call read_green_file(fname, gd, ierr, errmsg)
        if (ierr /= 0) call die(errmsg)
        if (.not. allocated(gd%kperp)) call die('file has no ARPES k_perp axis')
        if (ik < 1 .or. ik > gd%nk) call die('--ik out of range')
        if (trim(mapname) /= 'coh' .and. trim(mapname) /= 'incoh') call die("--map must be 'coh' or 'incoh'")
        if (maxval(abs(gd%kfrac(1:2, ik) - real(nint(gd%kfrac(1:2, ik)), dp))) > 1.0e-8_dp) &
            call die('photon conversion requires normal emission: select k_par = 0 modulo a surface reciprocal vector')

        file_eref = meta_get(gd, 'energy_reference')
        eref = file_eref
        if (trim(reference) /= 'auto') eref = trim(reference)
        select case (trim(eref))
        case ('E-E_F')
            if (have_fermi) call die('--fermi-energy is incompatible with --energy-reference E-E_F')
            fermi_energy_ev = 0.0_dp
        case ('absolute')
            if (.not. have_fermi) then
                if (trim(meta_get(gd, 'mode')) /= 'surface') &
                    call die('absolute energies in a film require an explicit --fermi-energy')
                a = meta_get(gd, 'substrate_fermi_ev')
                read(a, *, iostat=ierr) fermi_energy_ev
                if (ierr /= 0) call die('absolute energies require --fermi-energy (no valid substrate E_F in file)')
            end if
            if (.not. ieee_is_finite(fermi_energy_ev)) call die('--fermi-energy must be finite')
        case default
            call die('unknown energy reference: use --energy-reference absolute or E-E_F')
        end select

        allocate(data(gd%ne, gd%nq))
        if (trim(mapname) == 'coh') then
            data = gd%icoh(:, ik, :)
        else
            data = gd%iinc(:, ik, :)
        end if
        allocate(out_energy(gd%ne), out_kinetic(gd%ne), out_kperp(gd%ne), out_intensity(gd%ne))
        nvalid = 0
        nthreshold = 0
        noutside = 0
        do ie = 1, gd%ne
            call photon_to_kperp(photon_ev, gd%energy(ie) - fermi_energy_ev, work_function_ev, &
                                 inner_potential_ev, kperp, kinetic, ierr, errmsg)
            if (ierr == 2) then
                nthreshold = nthreshold + 1
                cycle
            else if (ierr /= 0) then
                call die(errmsg)
            end if
            if (gd%nq == 1) then
                if (kperp /= gd%kperp(1)) then
                    noutside = noutside + 1
                    cycle
                end if
                j = 1
                fraction = 0.0_dp
            else
                j = 0
                do iq = 1, gd%nq - 1
                    if (kperp >= gd%kperp(iq) .and. kperp <= gd%kperp(iq + 1)) then
                        j = iq
                        exit
                    end if
                end do
                if (j == 0) then
                    noutside = noutside + 1
                    cycle
                end if
                fraction = (kperp - gd%kperp(j)) / (gd%kperp(j + 1) - gd%kperp(j))
            end if
            nvalid = nvalid + 1
            out_energy(nvalid) = gd%energy(ie)
            out_kinetic(nvalid) = kinetic
            out_kperp(nvalid) = kperp
            if (gd%nq == 1) then
                out_intensity(nvalid) = data(ie, j)
            else
                out_intensity(nvalid) = (1.0_dp - fraction) * data(ie, j) + fraction * data(ie, j + 1)
            end if
        end do
        if (nvalid == 0) call die('no energies overlap the sampled k_perp axis; adjust &ARPES range or photon energy')

        if (len_trim(outp) == 0) then
            j = index(fname, '.green.bin', back=.true.)
            if (j > 1) then
                outp = fname(1:j - 1) // '.' // trim(mapname)
            else
                outp = fname // '.' // trim(mapname)
            end if
        end if
        open(newunit=u, file=trim(outp) // '.photon.dat', status='replace', action='write')
        write(u, '(A)') '# surfgreen normal-emission ARPES proxy of ' // fname // ', map = ' // trim(mapname)
        write(u, '(A)') '# final_state = nearly_free_electron; intensity interpolation = linear in k_perp; no extrapolation'
        write(u, '(A,ES18.9)') '# photon_energy_ev = ', photon_ev
        write(u, '(A,ES18.9)') '# work_function_ev = ', work_function_ev
        write(u, '(A,ES18.9)') '# inner_potential_V0_ev = ', inner_potential_ev
        write(u, '(A,A)') '# energy_reference = ', trim(eref)
        write(u, '(A,A)') '# file_energy_reference = ', trim(file_eref)
        write(u, '(A,ES18.9)') '# fermi_energy_ev = ', fermi_energy_ev
        write(u, '(A,I0)') '# k_point_index = ', ik
        write(u, '(A,3ES18.9)') '# k_fractional = ', gd%kfrac(:, ik)
        write(u, '(A,3ES18.9)') '# k_cart_invA = ', gd%kcart(:, ik)
        write(u, '(A,I0,A,I0)') '# skipped_below_threshold = ', nthreshold, '; skipped_outside_kperp = ', noutside
        write(u, '(A)') '# columns: initial_energy_eV kinetic_energy_eV k_perp_invA intensity'
        do ie = 1, nvalid
            write(u, '(4ES18.9)') out_energy(ie), out_kinetic(ie), out_kperp(ie), out_intensity(ie)
        end do
        close(u)
        print '(A,I0,A,I0,A)', 'Wrote ' // trim(outp) // '.photon.dat (', nvalid, &
            ' energies; ', noutside, ' outside k_perp range)'
    end subroutine photon_slice

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
