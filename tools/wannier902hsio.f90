! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Converte `seedname_hr.dat` do Wannier90 em H/S `.hsio`.
program wannier902hsio
    use types, only: dp
    use hs_io, only: HSDataFormat
    use wannier90_io, only: read_wannier90_hr, build_identity_overlap
    implicit none

    type(HSDataFormat) :: hr, sr
    character(len=1024) :: input_path, out_prefix, arg, orbital_spec
    character(len=:), allocatable :: errmsg
    integer, allocatable :: orbs(:)
    integer :: nargs, iarg, ierr, ios, ir0, i
    real(dp) :: fermi
    logical :: have_orbital_spec

    nargs = command_argument_count()
    if (nargs < 1) then
        call print_usage()
        stop 1
    end if
    call get_command_argument(1, input_path)
    if (trim(input_path) == '-h' .or. trim(input_path) == '--help') then
        call print_usage()
        stop 0
    end if
    out_prefix = default_prefix(trim(input_path))
    orbital_spec = ''
    have_orbital_spec = .false.
    fermi = 0.0_dp

    iarg = 2
    do while (iarg <= nargs)
        call get_command_argument(iarg, arg)
        select case (trim(arg))
        case ('--out-prefix')
            iarg = iarg + 1
            if (iarg > nargs) call argument_error('--out-prefix needs a value')
            call get_command_argument(iarg, out_prefix)
        case ('--orbs-per-atom')
            iarg = iarg + 1
            if (iarg > nargs) call argument_error('--orbs-per-atom needs a comma-separated list')
            call get_command_argument(iarg, orbital_spec)
            have_orbital_spec = .true.
        case ('--fermi')
            iarg = iarg + 1
            if (iarg > nargs) call argument_error('--fermi needs an energy in eV')
            call get_command_argument(iarg, arg)
            read(arg, *, iostat=ios) fermi
            if (ios /= 0) call argument_error('invalid value for --fermi')
        case ('-h', '--help')
            call print_usage()
            stop 0
        case default
            call argument_error("unknown argument '" // trim(arg) // "'")
        end select
        iarg = iarg + 1
    end do

    if (trim(input_path) == trim(out_prefix) // '_H.hsio' .or. &
        trim(input_path) == trim(out_prefix) // '_S.hsio') then
        call argument_error('an output path would overwrite the input')
    end if

    call read_wannier90_hr(trim(input_path), hr, ierr, errmsg)
    if (ierr /= 0) then
        write(*, '(A)') 'Erro: ' // errmsg
        stop 1
    end if

    if (have_orbital_spec) then
        call parse_positive_integer_list(trim(orbital_spec), orbs, ierr)
        if (ierr /= 0) then
            write(*, '(A)') 'Erro: --orbs-per-atom must be a list of positive integers'
            stop 1
        end if
        if (sum(orbs) /= hr%NOrbs) then
            write(*, '(A,I0)') 'Erro: --orbs-per-atom must be positive and sum to ', hr%NOrbs
            stop 1
        end if
        hr%nAtoms = size(orbs)
        deallocate(hr%Orbs)
        allocate(hr%Orbs(size(orbs)))
        hr%Orbs = orbs
    else
        write(*, '(A)') 'Aviso: sem --orbs-per-atom; cada funcao de Wannier sera tratada como um atomo.'
    end if

    if (fermi /= 0.0_dp) then
        ir0 = 0
        do i = 1, hr%Neighbors
            if (all(hr%IndexList(:, i) == 0)) ir0 = i
        end do
        if (ir0 == 0) then
            write(*, '(A)') 'Erro: nao ha bloco R=(0,0,0) para aplicar --fermi.'
            stop 1
        end if
        do i = 1, hr%NOrbs
            hr%Data(i, i, ir0) = hr%Data(i, i, ir0) - cmplx(fermi, 0.0_dp, dp)
        end do
    end if

    call build_identity_overlap(hr, sr, ierr, errmsg)
    if (ierr /= 0) then
        write(*, '(A)') 'Erro: ' // errmsg
        stop 1
    end if
    call hr%WriteFile(trim(out_prefix) // '_H.hsio')
    call sr%WriteFile(trim(out_prefix) // '_S.hsio')

    write(*, '(A,I0)') 'Funcoes de Wannier : ', hr%NOrbs
    write(*, '(A,I0)') 'Vetores R          : ', hr%Neighbors
    write(*, '(A,I0)') 'Atomos logicos     : ', hr%nAtoms
    write(*, '(A,ES14.6,A)') 'Shift de Fermi     : ', fermi, ' eV'
    write(*, '(A)') 'H                  : ' // trim(out_prefix) // '_H.hsio'
    write(*, '(A)') 'S=I                : ' // trim(out_prefix) // '_S.hsio'

contains

    subroutine parse_positive_integer_list(text, values, stat)
        character(len=*), intent(in) :: text
        integer, allocatable, intent(out) :: values(:)
        integer, intent(out) :: stat
        character(len=:), allocatable :: work
        integer :: i, n, ios_local

        work = trim(text)
        n = 1
        do i = 1, len(work)
            if (work(i:i) == ',') then
                work(i:i) = ' '
                n = n + 1
            end if
        end do
        allocate(values(n))
        read(work, *, iostat=ios_local) values
        stat = ios_local
        if (stat == 0) then
            if (any(values <= 0)) stat = 1
        end if
    end subroutine parse_positive_integer_list

    function default_prefix(path) result(prefix)
        character(len=*), intent(in) :: path
        character(len=1024) :: prefix
        integer :: last_slash, suffix
        prefix = trim(path)
        suffix = index(prefix, '_hr.dat', back=.true.)
        last_slash = index(prefix, '/', back=.true.)
        if (suffix > last_slash .and. suffix + 6 == len_trim(prefix)) then
            prefix = prefix(:suffix - 1)
        else if (len_trim(prefix) >= 4 .and. prefix(len_trim(prefix)-3:len_trim(prefix)) == '.dat') then
            prefix = prefix(:len_trim(prefix)-4)
        end if
    end function default_prefix

    subroutine argument_error(message)
        character(len=*), intent(in) :: message
        write(*, '(A)') 'Erro: ' // message
        call print_usage()
        stop 1
    end subroutine argument_error

    subroutine print_usage()
        write(*, '(A)') 'Usage: wannier902hsio seedname_hr.dat [options]'
        write(*, '(A)') '  --out-prefix PREFIX'
        write(*, '(A)') '  --orbs-per-atom "N1,N2,..."   (must sum to num_wann)'
        write(*, '(A)') '  --fermi E_EV                   subtract E_F from H(R=0)'
        write(*, '(A)') ''
        write(*, '(A)') 'The converter writes PREFIX_H.hsio and PREFIX_S.hsio (S=I).'
    end subroutine print_usage
end program wannier902hsio
