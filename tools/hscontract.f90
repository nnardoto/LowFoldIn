! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> hscontract — contração da base de H(R)/S(R) em arquivos .hsio.
!>
!>   hscontract contract.nml
!>
!> Lê H(R) e S(R) no formato .hsio, aplica uma transformação T bloco-diagonal
!> por átomo definida por regras radiais (&RADIAL_RULE, ver docs/hscontract.md)
!> e grava SOMENTE os arquivos reduzidos:
!>   spin none / noncollinear : <prefix>_H.hsio, <prefix>_S.hsio
!>   spin collinear          : <prefix>_H_up.hsio, <prefix>_H_dn.hsio, <prefix>_S.hsio
!> Neighbors/IndexList/Degenerecency, metadados de energia (E_F/referência),
!> unidades e a ordem atômica são preservados; Orbs/NOrbs e Title são
!> atualizados com proveniência. Não gera nem modifica structure.nml nem
!> arquivos laterais. Ferramentas que exigem structure.nml (spectral,
!> surfgreen, ...) precisam de uma contagem de orbitais compatível com a base
!> contratada — o usuário deve ajustar o bloco &SPECIES delas.
program hscontract
    use, intrinsic :: iso_fortran_env, only: error_unit
    use types, only: dp
    use HS_io, only: HSDataFormat
    use structure_io, only: structure_t, read_structure
    use basis_contract, only: contract_input_t, contraction_t, read_contract_input, build_contraction, &
        contract_hs, check_hs_compatible, check_pair_symmetry, lambda_min_grid, BC_PATH_LEN, BC_TOL_PAIR
    implicit none

    character(len=1024) :: arg
    character(len=512) :: errmsg
    integer :: ierr, nfile, i, j, na
    logical :: ex

    type(contract_input_t) :: inp
    type(contraction_t) :: ct
    type(structure_t) :: st
    type(HSDataFormat) :: hs(3), hs_red(3)
    character(len=BC_PATH_LEN) :: in_name(3), out_name(3)
    real(dp) :: smin, kmin(3)

    if (command_argument_count() /= 1) call usage()
    call get_command_argument(1, arg)
    if (trim(arg) == '-h' .or. trim(arg) == '--help') call usage()

    call read_contract_input(trim(arg), inp, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)

    call read_structure(trim(inp%structure_file), st, ierr, errmsg)
    if (ierr /= 0) call die('structure: ' // trim(errmsg))
    if (st%nat < 1) call die('structure has no atoms')

    call build_contraction(inp, st%label, ct, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)

    ! -- tabela de arquivos por modo de spin --------------------------------
    nfile = 0
    in_name = ''
    out_name = ''
    if (trim(inp%spin_mode) == 'collinear') then
        nfile = 3
        in_name(1) = inp%h_up_file; out_name(1) = trim(inp%output_prefix) // '_H_up.hsio'
        in_name(2) = inp%h_dn_file; out_name(2) = trim(inp%output_prefix) // '_H_dn.hsio'
        in_name(3) = inp%s_file;   out_name(3) = trim(inp%output_prefix) // '_S.hsio'
    else
        nfile = 2
        in_name(1) = inp%h_file; out_name(1) = trim(inp%output_prefix) // '_H.hsio'
        in_name(2) = inp%s_file; out_name(2) = trim(inp%output_prefix) // '_S.hsio'
    end if

    ! -- colisão de saídas com entradas -------------------------------------
    do i = 1, nfile
        if (same_path(out_name(i), inp%h_file) .or. same_path(out_name(i), inp%h_up_file) &
            .or. same_path(out_name(i), inp%h_dn_file) .or. same_path(out_name(i), inp%s_file) &
            .or. same_path(out_name(i), inp%structure_file) .or. same_path(out_name(i), trim(arg))) then
            call die("output '" // trim(out_name(i)) // "' would overwrite an input file" // &
                     ' (choose another output_prefix)')
        end if
    end do
    do i = 1, nfile
        do j = 1, i - 1
            if (same_path(out_name(i), out_name(j))) &
                call die("duplicate output name '" // trim(out_name(i)) // "'")
        end do
    end do

    ! -- leitura e validação -------------------------------------------------
    do i = 1, nfile
        inquire(file=trim(in_name(i)), exist=ex)
        if (.not. ex) call die("input file not found: '" // trim(in_name(i)) // "'")
        call hs(i)%ReadFile(trim(in_name(i)))
        if (hs(i)%nAtoms /= st%nat) then
            call die("'" // trim(in_name(i)) // "': nAtoms = " // itoa(hs(i)%nAtoms) // &
                     ' differs from the structure nat = ' // itoa(st%nat))
        end if
        if (i > 1) then
            call check_hs_compatible(hs(1), trim(in_name(1)), hs(i), trim(in_name(i)), ierr, errmsg)
            if (ierr /= 0) call die(errmsg)
        end if
        call check_pair_symmetry(hs(i), BC_TOL_PAIR, ierr, errmsg)
        if (ierr /= 0) call die("'" // trim(in_name(i)) // "': " // trim(errmsg))
    end do

    ! -- contração ------------------------------------------------------------
    do i = 1, nfile
        call contract_hs(hs(i), ct, hs_red(i), ierr, errmsg)
        if (ierr /= 0) call die("'" // trim(in_name(i)) // "': " // trim(errmsg))
    end do

    ! -- lambda_min(S_r(k)) na grade configurável, sem piso artificial ---------
    i = nfile   ! S é o último arquivo da tabela
    call lambda_min_grid(hs_red(i), inp%check_kgrid, inp%tol_s_min, smin, kmin, ierr, errmsg)
    if (ierr /= 0) call die('reduced S: ' // trim(errmsg))

    ! -- saídas -----------------------------------------------------------------
    do i = 1, nfile
        call hs_red(i)%WriteFile(trim(out_name(i)))
    end do

    ! -- resumo ------------------------------------------------------------------
    print '(A)', 'hscontract | basis contraction of H(R)/S(R)'
    print '(A,A)',   '  spin mode      : ', trim(inp%spin_mode)
    do i = 1, nfile
        print '(A,A,A,A)', '  ', trim(in_name(i)), ' -> ', trim(out_name(i))
    end do
    print '(A,I0,A,I0)', '  orbitals       : ', ct%norb_full, ' -> ', ct%norb_red
    print '(A,I0)', '  atoms          : ', ct%nat
    print '(A)', '  orbitals/atom  :'
    na = min(ct%nat, 20)
    do i = 1, na
        print '(A,I0,A,I0)', '    atom ', i, ' : ', ct%orbs_red(i)
    end do
    if (ct%nat > 20) print '(A)', '    ... (first 20 shown)'
    print '(A,I0,A,I0,A,I0)', '  check grid     : ', inp%check_kgrid(1), ' x ', inp%check_kgrid(2), &
        ' x ', inp%check_kgrid(3)
    print '(A,ES12.4,A,3F10.6)', '  lambda_min(S_r): ', smin, '  at k = ', kmin

contains

    subroutine usage()
        write(error_unit, '(A)') 'Usage: hscontract CONTRACT.nml   (&HSCONTRACT + &RADIAL_RULE blocks; see docs/hscontract.md)'
        stop 1
    end subroutine usage

    subroutine die(msg)
        character(len=*), intent(in) :: msg
        write(error_unit, '(A,A)') '[hscontract error] ', trim(msg)
        stop 1
    end subroutine die

    !> Igualdade exata de caminhos trim(adjustl(.)); a comparação não
    !> canonicaliza caminhos (documentado em docs/hscontract.md).
    logical function same_path(a, b)
        character(len=*), intent(in) :: a, b
        same_path = (len_trim(a) > 0 .and. trim(adjustl(a)) == trim(adjustl(b)))
    end function same_path

    pure function itoa(i)
        integer, intent(in) :: i
        character(len=16) :: itoa
        write(itoa, '(I0)') i
    end function itoa

end program hscontract
