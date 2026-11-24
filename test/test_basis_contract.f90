! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Testes do módulo basis_contract (hscontract):
!>  - identidade (sem regras) e bit-exatidão;
!>  - poda de shell completa (n_out = 0);
!>  - contração analítica com base NÃO ortogonal (valores calculados à mão);
!>  - metadados (Title/meta/Neighbors/IndexList/Degenerecency) e roundtrip HSIO;
!>  - índices/Orbs com múltiplas espécies e regras por l;
!>  - collinear e noncollinear (mesma T espacial nos setores/variáveis);
!>  - erros de input (namelist) e de posto de coeficientes;
!>  - lambda_min(S_r(k)): sucesso e falha com k identificado;
!>  - com `gen <dir>`: escreve um caso sintético para o teste CLI (bash).
program test_basis_contract
    use, intrinsic :: iso_fortran_env, only: int64
    use types, only: dp
    use testing, only: check, check_close, test_summary
    use HS_io, only: HSDataFormat
    use basis_contract, only: contract_input_t, contraction_t, radial_rule_t, read_contract_input, &
        build_contraction, contract_hs, check_hs_compatible, check_pair_symmetry, lambda_min_grid
    implicit none

    character(len=256) :: arg1, arg2
    integer :: nargs

    nargs = command_argument_count()
    arg1 = ''
    arg2 = ''
    if (nargs >= 1) call get_command_argument(1, arg1)
    if (nargs >= 2) call get_command_argument(2, arg2)
    if (trim(arg1) == 'gen') then
        call write_case(trim(arg2))
        print '(A)', 'test_basis_contract: synthetic case written to ' // trim(arg2)
        stop
    end if
    if (trim(arg1) == 'genbands') then
        call get_command_argument(2, arg2)
        call write_bands_case(trim(arg2))
        print '(A)', 'test_basis_contract: bands case written to ' // trim(arg2)
        stop
    end if
    if (trim(arg1) == 'bandcmp') then
        call run_bandcmp()
        stop
    end if

    call test_reader_ok()
    call test_reader_errors()
    call test_identity()
    call test_prune()
    call test_analytic()
    call test_metadata_roundtrip()
    call test_indices()
    call test_noncollinear()
    call test_collinear()
    call test_d_shell()
    call test_noncollinear_2atoms()
    call test_smin()
    call test_smin_hermiticity()
    call test_build_errors()
    call test_apply_errors()
    call test_title_long()
    call test_summary('test_basis_contract')

contains

    ! ------------------------------------------------------------------
    ! Auxiliares
    ! ------------------------------------------------------------------

    subroutine init_hs(h, orbs, rl, title)
        type(HSDataFormat), intent(out) :: h
        integer,            intent(in)  :: orbs(:)
        integer,            intent(in)  :: rl(:, :)
        character(len=*),   intent(in)  :: title
        integer :: nat, nr
        nat = size(orbs)
        nr = size(rl, 2)
        h%Title = title
        h%Neighbors = nr
        h%nAtoms = nat
        h%NOrbs = sum(orbs)
        allocate(h%Orbs(nat), h%IndexList(3, nr), h%Degenerecency(nr), h%Data(h%NOrbs, h%NOrbs, nr))
        h%Orbs = orbs
        h%IndexList = rl
        h%Degenerecency = 1
        h%Data = (0.0_dp, 0.0_dp)
    end subroutine init_hs

    subroutine set_rule(inp, ir, species, l, n_out, coeff)
        type(contract_input_t), intent(inout) :: inp
        integer,                intent(in)    :: ir, l, n_out
        character(len=*),       intent(in)    :: species
        real(dp),               intent(in)    :: coeff(:)
        integer :: i2, need
        if (.not. allocated(inp%rule)) allocate(inp%rule(0))
        need = ir - size(inp%rule)
        if (need > 0) inp%rule = [inp%rule, (radial_rule_t(), i2 = 1, need)]
        inp%rule(ir) = radial_rule_t(species, l, n_out)
        inp%rule(ir)%ncoeff = size(coeff)
        inp%rule(ir)%coeff(1:size(coeff)) = coeff
        inp%nrule = max(inp%nrule, ir)
    end subroutine set_rule

    subroutine expect_reader_error(text, msg)
        character(len=*), intent(in) :: text, msg
        type(contract_input_t) :: inp
        character(len=512) :: errmsg
        integer :: u, ierr
        open(newunit=u, file='test_basis_contract_tmp.nml', status='replace', action='write')
        write(u, '(A)') text
        close(u)
        call read_contract_input('test_basis_contract_tmp.nml', inp, ierr, errmsg)
        call check(ierr /= 0, msg)
        if (ierr == 0) print '(A)', '  note: reader accepted invalid input'
        open(newunit=u, file='test_basis_contract_tmp.nml', status='old')
        close(u, status='delete')
    end subroutine expect_reader_error

    ! ------------------------------------------------------------------
    ! Leitura de contract.nml
    ! ------------------------------------------------------------------

    subroutine test_reader_ok()
        type(contract_input_t) :: inp
        character(len=512) :: errmsg
        character(len=512) :: text
        integer :: u, ierr

        text = ''
        text = trim(text) // new_line('a') // '&HSCONTRACT'
        text = trim(text) // new_line('a') // "    spin_mode = 'collinear'"
        text = trim(text) // new_line('a') // "    h_up_file = 'a_H_up.hsio'"
        text = trim(text) // new_line('a') // "    h_dn_file = 'a_H_dn.hsio'"
        text = trim(text) // new_line('a') // "    s_file = 'a_S.hsio'"
        text = trim(text) // new_line('a') // "    structure_file = 'a.nml'"
        text = trim(text) // new_line('a') // "    output_prefix = 'red'"
        text = trim(text) // new_line('a') // "    basis = 'Bi s3p2d2', 'Sb s3p2d2'"
        text = trim(text) // new_line('a') // '/'
        text = trim(text) // new_line('a') // '&RADIAL_RULE'
        text = trim(text) // new_line('a') // "    species = 'Bi'"
        text = trim(text) // new_line('a') // '    l = 2'
        text = trim(text) // new_line('a') // '    n_out = 1'
        text = trim(text) // new_line('a') // '    coefficients = 0.5, 0.4'
        text = trim(text) // new_line('a') // '/'
        open(newunit=u, file='test_basis_contract_tmp.nml', status='replace', action='write')
        write(u, '(A)') trim(text)
        close(u)

        call read_contract_input('test_basis_contract_tmp.nml', inp, ierr, errmsg)
        call check(ierr == 0, 'reader: valid contract.nml')
        if (ierr == 0) then
            call check(trim(inp%spin_mode) == 'collinear', 'reader: spin_mode collinear')
            call check(trim(inp%output_prefix) == 'red', 'reader: output_prefix')
            call check(inp%nrule == 1 .and. inp%rule(1)%species == 'Bi', 'reader: rule species')
            call check(inp%rule(1)%l == 2 .and. inp%rule(1)%n_out == 1, 'reader: rule l/n_out')
            call check(inp%rule(1)%ncoeff == 2, 'reader: ncoeff = n_out*nzeta (d2)')
            call check_close(inp%rule(1)%coeff(2), 0.4_dp, 0.0_dp, 'reader: coeff order')
            call check(all(inp%check_kgrid == 2), 'reader: default check_kgrid = 2,2,2')
            call check_close(inp%tol_s_min, 1.0e-10_dp, 0.0_dp, 'reader: default tol_s_min')
            call check(trim(inp%basis(2)) == 'Sb s3p2d2', 'reader: basis entries')
        end if
        open(newunit=u, file='test_basis_contract_tmp.nml', status='old')
        close(u, status='delete')
    end subroutine test_reader_ok

    subroutine test_reader_errors()
        character(len=512) :: base
        base = ''
        base = trim(base) // new_line('a') // '&HSCONTRACT'
        base = trim(base) // new_line('a') // "    spin_mode = 'none'"
        base = trim(base) // new_line('a') // "    h_file = 'in_H.hsio'"
        base = trim(base) // new_line('a') // "    s_file = 'in_S.hsio'"
        base = trim(base) // new_line('a') // "    structure_file = 'st.nml'"
        base = trim(base) // new_line('a') // "    output_prefix = 'out'"
        base = trim(base) // new_line('a') // "    basis = 'Aa s2', 'Bb s2'"
        base = trim(base) // new_line('a') // '/'

        call expect_reader_error(replace(base, "spin_mode = 'none'", "spin_mode = 'diagonal'"), &
            'reader: invalid spin_mode')
        call expect_reader_error(replace(base, "    h_file = 'in_H.hsio'", ''), &
            'reader: h_file required in spin none')
        call expect_reader_error(replace(replace(base, "spin_mode = 'none'", "spin_mode = 'collinear'"), &
            "    h_file = 'in_H.hsio'", ''), 'reader: h_up/h_dn required in collinear')
        call expect_reader_error(replace(base, "    s_file = 'in_S.hsio'", ''), &
            'reader: s_file required')
        call expect_reader_error(replace(base, "    structure_file = 'st.nml'", ''), &
            'reader: structure_file required')
        call expect_reader_error(replace(base, "    output_prefix = 'out'", ''), &
            'reader: output_prefix required')
        call expect_reader_error(replace(base, "    basis = 'Aa s2', 'Bb s2'", ''), &
            'reader: basis required')
        call expect_reader_error(replace(base, "'Aa s2'", "'Aa s2x1'"), 'reader: invalid basis string')
        call expect_reader_error(replace(base, "'Aa s2', 'Bb s2'", "'Aa s2', 'Aa s2'"), &
            'reader: duplicate basis entry')
        call expect_reader_error(replace(base, "output_prefix = 'out'", "output_prefix = 'out'" // &
            new_line('a') // '    check_kgrid = 2, 0, 2'), 'reader: check_kgrid >= 1')
        call expect_reader_error(replace(base, "output_prefix = 'out'", "output_prefix = 'out'" // &
            new_line('a') // '    tol_s_min = -0.5'), 'reader: tol_s_min > 0')

        ! regras
        call expect_reader_error(trim(base) // new_line('a') // '&RADIAL_RULE' // new_line('a') // &
            "    species = 'Cc'" // new_line('a') // '    l = 0' // new_line('a') // '    n_out = 1' // &
            new_line('a') // '    coefficients = 1.0' // new_line('a') // '/', &
            'reader: rule species without basis entry')
        call expect_reader_error(trim(base) // new_line('a') // '&RADIAL_RULE' // new_line('a') // &
            "    species = 'Bb'" // new_line('a') // '    l = 1' // new_line('a') // '    n_out = 1' // &
            new_line('a') // '    coefficients = 1.0' // new_line('a') // '/', &
            'reader: rule for absent shell')
        call expect_reader_error(trim(base) // new_line('a') // '&RADIAL_RULE' // new_line('a') // &
            "    species = 'Aa'" // new_line('a') // '    l = 0' // new_line('a') // '    n_out = 3' // &
            new_line('a') // '    coefficients = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0' // new_line('a') // '/', &
            'reader: n_out > nzeta')
        call expect_reader_error(trim(base) // new_line('a') // '&RADIAL_RULE' // new_line('a') // &
            "    species = 'Aa'" // new_line('a') // '    l = 0' // new_line('a') // '    n_out = 2' // &
            new_line('a') // '    coefficients = 1.0, 1.0, 1.0' // new_line('a') // '/', &
            'reader: ncoeff /= n_out*nzeta')
        call expect_reader_error(trim(base) // new_line('a') // '&RADIAL_RULE' // new_line('a') // &
            "    species = 'Aa'" // new_line('a') // '    l = 0' // new_line('a') // '    n_out = 1' // &
            new_line('a') // '    coefficients = 1.0, 1.0' // new_line('a') // '/' // new_line('a') // &
            '&RADIAL_RULE' // new_line('a') // "    species = 'Aa'" // new_line('a') // '    l = 0' // &
            new_line('a') // '    n_out = 0' // new_line('a') // '/', &
            'reader: duplicate rule (species, l)')
        call expect_reader_error(trim(base) // new_line('a') // '&RADIAL_RULE' // new_line('a') // &
            "    species = 'Aa'" // new_line('a') // '    l = 5' // new_line('a') // '    n_out = 0' // &
            new_line('a') // '/', 'reader: l out of range')
        call expect_reader_error('&RADIAL_RULE' // new_line('a') // "    species = 'Aa'" // new_line('a') // &
            '/' // trim(base), 'reader: &HSCONTRACT must be first')
        ! typo de nome de grupo: sem a pré-varredura viraria identidade silenciosa
        call expect_reader_error(trim(base) // new_line('a') // '&RADIAL_RULES' // new_line('a') // &
            "    species = 'Aa'" // new_line('a') // '    l = 0' // new_line('a') // '    n_out = 0' // &
            new_line('a') // '/', 'reader: typo &RADIAL_RULES is rejected')
        call expect_reader_error(trim(base) // new_line('a') // '&EXTRA' // new_line('a') // '/', &
            'reader: unexpected namelist group is rejected')
        ! tol_s_min <= 0 explicitamente digitado (antes -1 era confundido com ausente)
        call expect_reader_error(replace(base, "output_prefix = 'out'", "output_prefix = 'out'" // &
            new_line('a') // '    tol_s_min = -1.0'), 'reader: tol_s_min = -1 is rejected')
        call expect_reader_error(replace(base, "output_prefix = 'out'", "output_prefix = 'out'" // &
            new_line('a') // '    tol_s_min = 0.0'), 'reader: tol_s_min = 0 is rejected')
    end subroutine test_reader_errors

    pure function replace(text, old, new) result(out)
        character(len=*), intent(in) :: text, old, new
        character(len=1024) :: out
        integer :: p
        out = text
        p = index(out, old)
        if (p > 0) out = out(1:p - 1) // new // out(p + len(old):)
    end function replace

    ! ------------------------------------------------------------------
    ! Contração
    ! ------------------------------------------------------------------

    subroutine test_identity()
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: h, hr
        character(len=8) :: labels(2)
        character(len=512) :: errmsg
        integer :: ierr, i, j, r
        real(dp) :: val

        labels = ['Aa', 'Aa']
        inp%spin_mode = 'none'
        inp%basis(1) = 'Aa s1p1'
        inp%nrule = 0
        call init_hs(h, [4, 4], reshape([0, 0, 0, 1, 0, 0, -1, 0, 0], [3, 3]), 'ident')
        do r = 1, 3
            do j = 1, 8
                do i = 1, 8
                    if (r == 1) then
                        ! R = 0: real e simétrica (hermitiana)
                        val = 0.05_dp * (i + j)
                        h%Data(i, j, r) = cmplx(val, 0.0_dp, kind=dp)
                    else
                        val = 0.1_dp * (10 * r + i) + 0.01_dp * j
                        h%Data(i, j, r) = cmplx(val, 0.5_dp * val, kind=dp)
                    end if
                end do
            end do
        end do
        ! estrutura algébrica: H(-R) = H(R)^H
        h%Data(:, :, 3) = transpose(conjg(h%Data(:, :, 2)))

        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'identity: build')
        if (ierr /= 0) return
        call check(all(ct%orbs_red == [4, 4]) .and. ct%norb_red == 8, 'identity: orbs preserved')
        call check(ct%norb_full == 8, 'identity: norb_full')
        call contract_hs(h, ct, hr, ierr, errmsg)
        call check(ierr == 0, 'identity: contract')
        if (ierr /= 0) return
        call check(all(hr%Data == h%Data), 'identity: Data bit-exact')
        call check(all(hr%Orbs == h%Orbs) .and. hr%NOrbs == h%NOrbs, 'identity: Orbs/NOrbs')
        call check(index(hr%Title, 'hscontract') > 0, 'identity: Title provenance')
        call check(all(hr%IndexList == h%IndexList) .and. hr%Neighbors == h%Neighbors, &
            'identity: IndexList/Neighbors')
        call check_pair_symmetry(h, 1.0e-9_dp, ierr, errmsg)
        call check(ierr == 0, 'identity: pair symmetry ok')
    end subroutine test_identity

    subroutine test_prune()
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: h, hr
        character(len=8) :: labels(1)
        character(len=512) :: errmsg
        integer :: ierr, i

        labels = ['Aa']
        inp%spin_mode = 'none'
        inp%basis(1) = 'Aa s1p1'
        call set_rule(inp, 1, 'Aa', 1, 0, [real(dp) ::])
        call init_hs(h, [4], reshape([0, 0, 0], [3, 1]), 'prune')
        do i = 1, 4
            h%Data(i, i, 1) = cmplx(real(i, dp), 0.0_dp, kind=dp)
        end do

        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'prune: build')
        if (ierr /= 0) return
        call check(ct%orbs_red(1) == 1, 'prune: p shell removed (4 -> 1)')
        call contract_hs(h, ct, hr, ierr, errmsg)
        call check(ierr == 0, 'prune: contract')
        if (ierr /= 0) return
        call check(hr%NOrbs == 1 .and. hr%Orbs(1) == 1, 'prune: NOrbs/Orbs updated')
        call check_close(real(hr%Data(1, 1, 1)), 1.0_dp, 1.0e-14_dp, 'prune: s orbital kept')
    end subroutine test_prune

    subroutine test_analytic()
        ! 2 átomos, base NÃO ortogonal; valores calculados à mão:
        !   Aa: c=[0.8,0.5], Bb: c=[0.6,0.3]
        !   S_A=[[1,.3],[.3,.9]] -> Sr_AA=1.105; S_B=[[1,.2],[.2,.8]] -> Sr_BB=.504
        !   S_AB=[[.1,.2],[.3,.4]] -> Sr_AB=.246
        !   H_A=[[2,-1],[-1,3]] -> 1.23; H_B=[[.5,.2],[.2,-.7]] -> .189
        !   H_AB=[[.05,0],[.1,-.05]] -> 0.0465
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: s, sr
        character(len=8) :: labels(2)
        character(len=512) :: errmsg
        real(dp) :: sa(2, 2), sb(2, 2), sab(2, 2), ha(2, 2), hb(2, 2), hab(2, 2)
        integer :: ierr

        labels = ['Aa', 'Bb']
        inp%spin_mode = 'none'
        inp%basis(1) = 'Aa s2'
        inp%basis(2) = 'Bb s2'
        call set_rule(inp, 1, 'Aa', 0, 1, [0.8_dp, 0.5_dp])
        call set_rule(inp, 2, 'Bb', 0, 1, [0.6_dp, 0.3_dp])

        sa = reshape([1.0_dp, 0.3_dp, 0.3_dp, 0.9_dp], [2, 2])
        sb = reshape([1.0_dp, 0.2_dp, 0.2_dp, 0.8_dp], [2, 2])
        sab = reshape([0.1_dp, 0.3_dp, 0.2_dp, 0.4_dp], [2, 2])
        ha = reshape([2.0_dp, -1.0_dp, -1.0_dp, 3.0_dp], [2, 2])
        hb = reshape([0.5_dp, 0.2_dp, 0.2_dp, -0.7_dp], [2, 2])
        hab = reshape([0.05_dp, 0.1_dp, 0.0_dp, -0.05_dp], [2, 2])

        call init_hs(s, [2, 2], reshape([0, 0, 0], [3, 1]), 'S')
        s%Data(1:2, 1:2, 1) = cmplx(sa, 0.0_dp, kind=dp)
        s%Data(3:4, 3:4, 1) = cmplx(sb, 0.0_dp, kind=dp)
        s%Data(1:2, 3:4, 1) = cmplx(sab, 0.0_dp, kind=dp)
        s%Data(3:4, 1:2, 1) = cmplx(transpose(sab), 0.0_dp, kind=dp)

        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'analytic: build')
        if (ierr /= 0) return
        call check(all(ct%orbs_red == [1, 1]) .and. ct%norb_red == 2, 'analytic: orbs 2,2 -> 1,1')
        call contract_hs(s, ct, sr, ierr, errmsg)
        call check(ierr == 0, 'analytic: contract S')
        if (ierr /= 0) return
        call check_close(real(sr%Data(1, 1, 1)), 1.105_dp, 1.0e-13_dp, 'analytic: Sr_AA = 1.105')
        call check_close(real(sr%Data(2, 2, 1)), 0.504_dp, 1.0e-13_dp, 'analytic: Sr_BB = 0.504')
        call check_close(real(sr%Data(1, 2, 1)), 0.246_dp, 1.0e-13_dp, 'analytic: Sr_AB = 0.246')
        call check_close(real(sr%Data(2, 1, 1)), 0.246_dp, 1.0e-13_dp, 'analytic: Sr_BA = 0.246')

        ! H com a mesma contração
        block
            type(HSDataFormat) :: h, hred
            call init_hs(h, [2, 2], reshape([0, 0, 0], [3, 1]), 'H')
            h%Data(1:2, 1:2, 1) = cmplx(ha, 0.0_dp, kind=dp)
            h%Data(3:4, 3:4, 1) = cmplx(hb, 0.0_dp, kind=dp)
            h%Data(1:2, 3:4, 1) = cmplx(hab, 0.0_dp, kind=dp)
            h%Data(3:4, 1:2, 1) = cmplx(transpose(hab), 0.0_dp, kind=dp)
            call contract_hs(h, ct, hred, ierr, errmsg)
            call check(ierr == 0, 'analytic: contract H')
            if (ierr /= 0) return
            call check_close(real(hred%Data(1, 1, 1)), 1.23_dp, 1.0e-13_dp, 'analytic: Hr_AA = 1.23')
            call check_close(real(hred%Data(2, 2, 1)), 0.189_dp, 1.0e-13_dp, 'analytic: Hr_BB = 0.189')
            call check_close(real(hred%Data(1, 2, 1)), 0.0465_dp, 1.0e-13_dp, 'analytic: Hr_AB = 0.0465')
        end block
    end subroutine test_analytic

    subroutine test_metadata_roundtrip()
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: s, sr, back
        character(len=8) :: labels(2)
        character(len=512) :: errmsg
        character(len=*), parameter :: fname = 'test_basis_contract_rt.hsio'
        integer :: ierr, u

        labels = ['Aa', 'Bb']
        inp%spin_mode = 'none'
        inp%basis(1) = 'Aa s2'
        inp%basis(2) = 'Bb s2'
        call set_rule(inp, 1, 'Aa', 0, 1, [0.8_dp, 0.5_dp])
        call set_rule(inp, 2, 'Bb', 0, 1, [0.6_dp, 0.3_dp])

        call init_hs(s, [2, 2], reshape([0, 0, 0, 1, 0, 0, -1, 0, 0], [3, 3]), 'orig title')
        s%Degenerecency = [1, 3, 2]
        s%Data(1, 1, 1) = (1.0_dp, 0.0_dp)
        s%Data(2, 2, 1) = (1.0_dp, 0.0_dp)
        s%Data(3, 3, 1) = (1.0_dp, 0.0_dp)
        s%Data(4, 4, 1) = (1.0_dp, 0.0_dp)
        s%Data(1, 2, 2) = (0.2_dp, 0.1_dp)
        s%Data(2, 1, 3) = (0.2_dp, -0.1_dp)
        s%has_meta = .true.
        s%fermi_ev = -3.456789_dp
        s%fermi_shift_applied = .true.
        s%energy_reference = 'E-E_F'

        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'meta: build')
        if (ierr /= 0) return
        call contract_hs(s, ct, sr, ierr, errmsg)
        call check(ierr == 0, 'meta: contract')
        if (ierr /= 0) return
        call check(sr%has_meta .and. sr%fermi_shift_applied, 'meta: has_meta/shift preserved')
        call check_close(sr%fermi_ev, -3.456789_dp, 0.0_dp, 'meta: fermi_ev preserved')
        call check(trim(sr%energy_reference) == 'E-E_F', 'meta: energy_reference preserved')
        call check(all(sr%Degenerecency == [1, 3, 2]), 'meta: Degenerecency preserved')
        call check(all(sr%IndexList == s%IndexList) .and. sr%Neighbors == 3, 'meta: IndexList preserved')
        call check(sr%nAtoms == 2 .and. all(sr%Orbs == [1, 1]), 'meta: atom order/Orbs updated')
        call check(index(sr%Title, 'orig title') == 1 .and. index(sr%Title, 'hscontract') > 0, &
            'meta: Title keeps original + provenance')
        call check_close(real(sr%Data(1, 1, 1)), 0.89_dp, 1.0e-13_dp, 'meta: contracted value')

        ! roundtrip HSIO: escreve o reduzido e lê de volta
        call sr%WriteFile(fname)
        call back%ReadFile(fname)
        call check(all(back%Data == sr%Data), 'roundtrip: Data')
        call check(back%NOrbs == 2 .and. all(back%Orbs == [1, 1]), 'roundtrip: Orbs')
        call check(all(back%IndexList == sr%IndexList), 'roundtrip: IndexList')
        call check(back%Title == sr%Title, 'roundtrip: Title')
        call check(back%has_meta .and. back%fermi_ev == sr%fermi_ev, 'roundtrip: meta')
        open(newunit=u, file=fname, status='old')
        close(u, status='delete')
    end subroutine test_metadata_roundtrip

    subroutine test_indices()
        ! Aa s2p1 (regra s: [0.7,0.7]) e Bb s1p2 (regra p: [0.6,0.8]):
        ! H(i,i) = i; H_Aa(s1,s2) = 0.2; H_Bb(p_z1, p_z2 mesmo m) = 0.1
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: h, hr
        character(len=8) :: labels(2)
        character(len=512) :: errmsg
        integer :: ierr, i

        labels = ['Aa', 'Bb']
        inp%spin_mode = 'none'
        inp%basis(1) = 'Aa s2p1'
        inp%basis(2) = 'Bb s1p2'
        call set_rule(inp, 1, 'Aa', 0, 1, [0.7_dp, 0.7_dp])
        call set_rule(inp, 2, 'Bb', 1, 1, [0.6_dp, 0.8_dp])

        call init_hs(h, [5, 7], reshape([0, 0, 0], [3, 1]), 'idx')
        do i = 1, 12
            h%Data(i, i, 1) = cmplx(real(i, dp), 0.0_dp, kind=dp)
        end do
        h%Data(1, 2, 1) = (0.2_dp, 0.0_dp)
        h%Data(2, 1, 1) = (0.2_dp, 0.0_dp)
        ! Bb: s = 6; p zeta 1 = 7, 8, 9; p zeta 2 = 10, 11, 12 (mesmo m acoplado)
        h%Data(7, 10, 1) = (0.1_dp, 0.0_dp);  h%Data(10, 7, 1) = (0.1_dp, 0.0_dp)
        h%Data(8, 11, 1) = (0.1_dp, 0.0_dp);  h%Data(11, 8, 1) = (0.1_dp, 0.0_dp)
        h%Data(9, 12, 1) = (0.1_dp, 0.0_dp);  h%Data(12, 9, 1) = (0.1_dp, 0.0_dp)

        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'idx: build')
        if (ierr /= 0) return
        call check(all(ct%orbs_full == [5, 7]) .and. ct%norb_full == 12, 'idx: orbs_full from basis')
        call check(all(ct%orbs_red == [4, 4]) .and. ct%norb_red == 8, 'idx: orbs_red = [4, 4]')
        call check(all(ct%off_red == [0, 4, 8]), 'idx: reduced offsets')
        call contract_hs(h, ct, hr, ierr, errmsg)
        call check(ierr == 0, 'idx: contract')
        if (ierr /= 0) return
        call check(all(hr%Orbs == [4, 4]) .and. hr%NOrbs == 8, 'idx: hs_red Orbs/NOrbs')
        call check_close(real(hr%Data(1, 1, 1)), 0.49_dp * 3.4_dp, 1.0e-13_dp, 'idx: Aa s contraction')
        call check_close(real(hr%Data(2, 2, 1)), 3.0_dp, 1.0e-13_dp, 'idx: Aa px identity')
        call check_close(real(hr%Data(5, 5, 1)), 6.0_dp, 1.0e-13_dp, 'idx: Bb s identity')
        call check_close(real(hr%Data(6, 6, 1)), 9.016_dp, 1.0e-13_dp, 'idx: Bb px contraction')
        call check_close(real(hr%Data(8, 8, 1)), 11.016_dp, 1.0e-13_dp, 'idx: Bb pz contraction')
    end subroutine test_indices

    subroutine test_noncollinear()
        ! átomo único Cc s2, spinor [up|dn]; T = blkdiag(c, c), c = [0.9, 0.2]
        ! H(0) = [[Hu, X], [X^T, Hd]] com Hu, Hd simétricas:
        !   c Hu c = 0.902, c Hd c = 1.616, c X c = 0.278
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: h, s, hred, sred
        character(len=8) :: labels(1)
        character(len=512) :: errmsg
        real(dp) :: hu(2, 2), x(2, 2), hd(2, 2), smin, kmin(3)
        integer :: ierr

        labels = ['Cc']
        inp%spin_mode = 'noncollinear'
        inp%basis(1) = 'Cc s2'
        call set_rule(inp, 1, 'Cc', 0, 1, [0.9_dp, 0.2_dp])

        hu = reshape([1.0_dp, 0.2_dp, 0.2_dp, 0.5_dp], [2, 2])
        x = reshape([0.3_dp, 0.1_dp, 0.05_dp, 0.2_dp], [2, 2])
        hd = reshape([2.0_dp, 0.1_dp, 0.1_dp, -1.0_dp], [2, 2])

        call init_hs(h, [4], reshape([0, 0, 0, 1, 0, 0, -1, 0, 0], [3, 3]), 'nc')
        h%Data(1:2, 1:2, 1) = cmplx(hu, 0.0_dp, kind=dp)
        h%Data(3:4, 3:4, 1) = cmplx(hd, 0.0_dp, kind=dp)
        h%Data(1:2, 3:4, 1) = cmplx(x, 0.0_dp, kind=dp)
        h%Data(3:4, 1:2, 1) = cmplx(transpose(x), 0.0_dp, kind=dp)
        ! hopping em +-a1 nos setores up e down (mesma T espacial)
        h%Data(1, 1, 2) = (0.4_dp, 0.0_dp); h%Data(1, 1, 3) = (0.4_dp, 0.0_dp)
        h%Data(3, 3, 2) = (0.8_dp, 0.0_dp); h%Data(3, 3, 3) = (0.8_dp, 0.0_dp)

        call init_hs(s, [4], reshape([0, 0, 0, 1, 0, 0, -1, 0, 0], [3, 3]), 'nc S')
        s%Data(1, 1, 1) = (1.0_dp, 0.0_dp); s%Data(2, 2, 1) = (1.0_dp, 0.0_dp)
        s%Data(3, 3, 1) = (1.0_dp, 0.0_dp); s%Data(4, 4, 1) = (1.0_dp, 0.0_dp)

        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'nc: build')
        if (ierr /= 0) return
        call check(ct%spinor == 2 .and. ct%orbs_full(1) == 4, 'nc: spinor doubling')
        call check(ct%orbs_red(1) == 2, 'nc: orbs_red = 2')
        call check_pair_symmetry(h, 1.0e-9_dp, ierr, errmsg)
        call check(ierr == 0, 'nc: pair symmetry (spinor)')

        call contract_hs(h, ct, hred, ierr, errmsg)
        call check(ierr == 0, 'nc: contract H')
        if (ierr /= 0) return
        call check_close(real(hred%Data(1, 1, 1)), 0.902_dp, 1.0e-13_dp, 'nc: up-up block')
        call check_close(real(hred%Data(2, 2, 1)), 1.616_dp, 1.0e-13_dp, 'nc: dn-dn block')
        call check_close(real(hred%Data(1, 2, 1)), 0.278_dp, 1.0e-13_dp, 'nc: spin-flip block')
        call check_close(real(hred%Data(1, 1, 2)), 0.81_dp * 0.4_dp, 1.0e-13_dp, 'nc: up hopping')
        call check_close(real(hred%Data(2, 2, 2)), 0.81_dp * 0.8_dp, 1.0e-13_dp, 'nc: dn hopping')

        ! S = I -> S_r = 0.85 I (blkdiag(c, c)); lambda_min na grade 2x2x2
        call contract_hs(s, ct, sred, ierr, errmsg)
        call check(ierr == 0, 'nc: contract S')
        if (ierr /= 0) return
        call lambda_min_grid(sred, [2, 2, 2], 1.0e-10_dp, smin, kmin, ierr, errmsg)
        call check(ierr == 0, 'nc: lambda_min ok')
        call check_close(smin, 0.85_dp, 1.0e-12_dp, 'nc: lambda_min = 0.85')
    end subroutine test_noncollinear

    subroutine test_collinear()
        ! mesma T espacial em H_up, H_dn e S (canais independentes)
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: hu, hd, s, hur, hdr, sr
        character(len=8) :: labels(1)
        character(len=512) :: errmsg
        real(dp) :: m
        integer :: ierr

        labels = ['Dd']
        inp%spin_mode = 'collinear'
        inp%basis(1) = 'Dd s2'
        call set_rule(inp, 1, 'Dd', 0, 1, [0.8_dp, 0.5_dp])

        call init_hs(hu, [2], reshape([0, 0, 0], [3, 1]), 'up')
        hu%Data(1:2, 1:2, 1) = cmplx(reshape([1.0_dp, 0.2_dp, 0.2_dp, 0.5_dp], [2, 2]), 0.0_dp, kind=dp)
        call init_hs(hd, [2], reshape([0, 0, 0], [3, 1]), 'dn')
        hd%Data(1:2, 1:2, 1) = cmplx(reshape([2.0_dp, 0.1_dp, 0.1_dp, -1.0_dp], [2, 2]), 0.0_dp, kind=dp)
        call init_hs(s, [2], reshape([0, 0, 0], [3, 1]), 'S')
        s%Data(1:2, 1:2, 1) = cmplx(reshape([1.0_dp, 0.3_dp, 0.3_dp, 0.9_dp], [2, 2]), 0.0_dp, kind=dp)

        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'coll: build')
        if (ierr /= 0) return
        call check(ct%spinor == 1 .and. all(ct%orbs_full == [2]), 'coll: spatial T (no doubling)')
        call contract_hs(hu, ct, hur, ierr, errmsg)
        if (ierr /= 0) then; call check(.false., 'coll: contract up'); return; end if
        call contract_hs(hd, ct, hdr, ierr, errmsg)
        if (ierr /= 0) then; call check(.false., 'coll: contract dn'); return; end if
        call contract_hs(s, ct, sr, ierr, errmsg)
        if (ierr /= 0) then; call check(.false., 'coll: contract S'); return; end if
        m = 0.8_dp * 1.0_dp + 0.5_dp * 0.2_dp
        call check_close(real(hur%Data(1, 1, 1)), 0.8_dp * m + 0.5_dp * (0.8_dp * 0.2_dp + 0.5_dp * 0.5_dp), &
            1.0e-13_dp, 'coll: H_up contracted')
        m = 0.8_dp * 2.0_dp + 0.5_dp * 0.1_dp
        call check_close(real(hdr%Data(1, 1, 1)), 0.8_dp * m + 0.5_dp * (0.8_dp * 0.1_dp - 0.5_dp * 1.0_dp), &
            1.0e-13_dp, 'coll: H_dn contracted')
        call check_close(real(sr%Data(1, 1, 1)), 1.105_dp, 1.0e-13_dp, 'coll: S contracted')
    end subroutine test_collinear

    ! ------------------------------------------------------------------
    ! lambda_min(S_r(k))
    ! ------------------------------------------------------------------

    subroutine test_smin()
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: s, sr
        character(len=8) :: labels(1)
        character(len=512) :: errmsg
        real(dp) :: eps, smin, kmin(3)
        integer :: ierr

        ! S quase singular na direção [1,-1]: S_r = 2*eps <= tol_s_min
        eps = 1.0e-12_dp
        labels = ['Yy']
        inp%spin_mode = 'none'
        inp%basis(1) = 'Yy s2'
        call set_rule(inp, 1, 'Yy', 0, 1, [1.0_dp, -1.0_dp])
        call init_hs(s, [2], reshape([0, 0, 0], [3, 1]), 'S')
        s%Data(1:2, 1:2, 1) = cmplx(reshape([1.0_dp, 1.0_dp - eps, 1.0_dp - eps, 1.0_dp], [2, 2]), 0.0_dp, kind=dp)
        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'smin: build')
        if (ierr /= 0) return
        call contract_hs(s, ct, sr, ierr, errmsg)
        call check(ierr == 0, 'smin: contract')
        if (ierr /= 0) return
        call lambda_min_grid(sr, [2, 2, 2], 1.0e-10_dp, smin, kmin, ierr, errmsg)
        call check(ierr == 1, 'smin: non-positive S_r fails')
        call check(index(errmsg, 'k =') > 0, 'smin: failing k is identified')
        call check_close(smin, 2.0_dp * eps, 1.0e-16_dp, 'smin: lambda value reported')

        ! caso positivo
        call set_rule(inp, 1, 'Yy', 0, 1, [1.0_dp, 1.0_dp])
        call init_hs(s, [2], reshape([0, 0, 0], [3, 1]), 'S')
        s%Data(1, 1, 1) = (1.0_dp, 0.0_dp)
        s%Data(2, 2, 1) = (1.0_dp, 0.0_dp)
        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'smin: build (positive)')
        if (ierr /= 0) return
        call contract_hs(s, ct, sr, ierr, errmsg)
        call check(ierr == 0, 'smin: contract (positive)')
        if (ierr /= 0) return
        call lambda_min_grid(sr, [2, 2, 2], 1.0e-10_dp, smin, kmin, ierr, errmsg)
        call check(ierr == 0, 'smin: positive S_r passes')
        call check_close(smin, 2.0_dp, 1.0e-12_dp, 'smin: lambda = 2 at k = 0')
    end subroutine test_smin

    subroutine test_smin_hermiticity()
        ! listas R sem o par -R dão S(k) nao hermitiana em k generico (a fase
        ! e^{i2k} nao e cancelada); a checagem deve rejeitar nomeando o k.
        ! Com os pares +-R completos e fases complexas, passa em todo k != 0.
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: s, sr
        character(len=8) :: labels(2)
        character(len=512) :: errmsg
        real(dp) :: smin, kmin(3)
        integer :: ierr, ir, nR

        labels = ['Aa', 'Aa']
        inp%spin_mode = 'none'
        inp%basis(1) = 'Aa s2'
        call set_rule(inp, 1, 'Aa', 0, 1, [0.8_dp, 0.5_dp])

        ! -- caso correto: R = 0, +-a1, +-a2 com S(-R) = S(R)^H (fases complexas)
        nR = 5
        call init_hs(s, [2, 2], reshape([0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0], [3, nR]), &
                     'Sc')
        s%Data(1:2, 1:2, 1) = cmplx(reshape([1.0_dp, 0.3_dp, 0.3_dp, 1.0_dp], [2, 2]), 0.0_dp, kind=dp)
        s%Data(3:4, 3:4, 1) = cmplx(reshape([1.0_dp, 0.2_dp, 0.2_dp, 1.0_dp], [2, 2]), 0.0_dp, kind=dp)
        do ir = 1, 2
            s%Data(1, 1, ir + 1) = cmplx(0.05_dp, 0.02_dp * ir, kind=dp)
            s%Data(2, 2, ir + 1) = cmplx(0.02_dp, -0.01_dp * ir, kind=dp)
        end do
        s%Data(:, :, 3) = transpose(conjg(s%Data(:, :, 2)))
        s%Data(:, :, 5) = transpose(conjg(s%Data(:, :, 4)))
        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'herm: build')
        if (ierr /= 0) return
        call contract_hs(s, ct, sr, ierr, errmsg)
        call check(ierr == 0, 'herm: contract')
        if (ierr /= 0) return
        call lambda_min_grid(sr, [3, 3, 1], 1.0e-10_dp, smin, kmin, ierr, errmsg)
        call check(ierr == 0, 'herm: full +-R list passes at k = 1/3, 2/3')
        call check(smin > 0.0_dp, 'herm: smin > 0 with complex phases')

        ! -- caso sem o par -R: k = 1/3 tem fase e^{i2pi/3} nao cancelada
        nR = 2
        call init_hs(s, [2, 2], reshape([0, 0, 0, 1, 0, 0], [3, nR]), 'Sx')
        s%Data(1:2, 1:2, 1) = cmplx(reshape([1.0_dp, 0.3_dp, 0.3_dp, 1.0_dp], [2, 2]), 0.0_dp, kind=dp)
        s%Data(3:4, 3:4, 1) = cmplx(reshape([1.0_dp, 0.2_dp, 0.2_dp, 1.0_dp], [2, 2]), 0.0_dp, kind=dp)
        s%Data(1, 1, 2) = cmplx(0.05_dp, 0.02_dp, kind=dp)
        s%Data(2, 2, 2) = cmplx(0.02_dp, -0.01_dp, kind=dp)
        call contract_hs(s, ct, sr, ierr, errmsg)
        call check(ierr == 0, 'herm: contract (no -R)')
        if (ierr /= 0) return
        call lambda_min_grid(sr, [3, 3, 1], 1.0e-10_dp, smin, kmin, ierr, errmsg)
        call check(ierr == 1, 'herm: no -R list is rejected at k != 0')
        call check(index(errmsg, 'hermitian') > 0, 'herm: message names hermiticity')
        call check(index(errmsg, 'k =') > 0, 'herm: message names the k point')
    end subroutine test_smin_hermiticity


    subroutine test_d_shell()
        ! shell d completa (todos os m) com contração d2 -> 1:
        ! red_d(m) = 0.36*v1(m) + 2*0.18*0.1 + 0.09*v2(m), mesma fórmula para
        ! cada m; shell s preservada (identidade).
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: h, hr
        character(len=8) :: labels(1)
        character(len=512) :: errmsg
        real(dp) :: expect
        integer :: ierr, i, m

        labels = ['Ee']
        inp%spin_mode = 'none'
        inp%basis(1) = 'Ee s1d2'
        call set_rule(inp, 1, 'Ee', 2, 1, [0.6_dp, 0.3_dp])

        call init_hs(h, [11], reshape([0, 0, 0], [3, 1]), 'd')
        do i = 1, 11
            h%Data(i, i, 1) = cmplx(real(i, dp), 0.0_dp, kind=dp)
        end do
        ! d zeta 1 (orbitais 2..6) acopla com d zeta 2 (7..11) no mesmo m
        do m = 1, 5
            h%Data(1 + m, 6 + m, 1) = (0.1_dp, 0.0_dp)
            h%Data(6 + m, 1 + m, 1) = (0.1_dp, 0.0_dp)
        end do

        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'd: build')
        if (ierr /= 0) return
        call check(ct%orbs_full(1) == 11 .and. ct%orbs_red(1) == 6, 'd: 11 -> 6 orbitals')
        call contract_hs(h, ct, hr, ierr, errmsg)
        call check(ierr == 0, 'd: contract')
        if (ierr /= 0) return
        call check_close(real(hr%Data(1, 1, 1)), 1.0_dp, 1.0e-14_dp, 'd: s identity')
        do m = 1, 5
            expect = 0.36_dp * (1 + m) + 2.0_dp * 0.6_dp * 0.3_dp * 0.1_dp + 0.09_dp * (6 + m)
            call check_close(real(hr%Data(1 + m, 1 + m, 1)), expect, 1.0e-13_dp, 'd: same contraction for every m')
        end do
    end subroutine test_d_shell

    subroutine test_noncollinear_2atoms()
        ! 2 átomos com espinores [up|dn] e acoplamento spin-flip ENTRE átomos;
        ! a mesma T espacial contrai os quatro setores.
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: h, s, hr, sred
        character(len=8) :: labels(2)
        character(len=512) :: errmsg
        real(dp) :: hu(2, 2), x(2, 2), hd(2, 2), smin, kmin(3)
        integer :: ierr, i

        labels = ['Cc', 'Cc']
        inp%spin_mode = 'noncollinear'
        inp%basis(1) = 'Cc s2'
        call set_rule(inp, 1, 'Cc', 0, 1, [0.9_dp, 0.2_dp])

        hu = reshape([1.0_dp, 0.2_dp, 0.2_dp, 0.5_dp], [2, 2])
        x = reshape([0.3_dp, 0.1_dp, 0.05_dp, 0.2_dp], [2, 2])
        hd = reshape([2.0_dp, 0.1_dp, 0.1_dp, -1.0_dp], [2, 2])

        call init_hs(h, [4, 4], reshape([0, 0, 0, 1, 0, 0, -1, 0, 0], [3, 3]), 'nc2')
        ! on-site por átomo: [[Hu, X], [X^T, Hd]]
        h%Data(1:2, 1:2, 1) = cmplx(hu, 0.0_dp, kind=dp)
        h%Data(3:4, 3:4, 1) = cmplx(hd, 0.0_dp, kind=dp)
        h%Data(1:2, 3:4, 1) = cmplx(x, 0.0_dp, kind=dp)
        h%Data(3:4, 1:2, 1) = cmplx(transpose(x), 0.0_dp, kind=dp)
        h%Data(5:6, 5:6, 1) = cmplx(hu, 0.0_dp, kind=dp)
        h%Data(7:8, 7:8, 1) = cmplx(hd, 0.0_dp, kind=dp)
        h%Data(5:6, 7:8, 1) = cmplx(x, 0.0_dp, kind=dp)
        h%Data(7:8, 5:6, 1) = cmplx(transpose(x), 0.0_dp, kind=dp)
        ! acoplamento spin-flip entre átomos no R = 0:
        ! bloco (a1, a2) = [[0, x12*I], [x21*I, 0]], (a2, a1) = transposto
        h%Data(1, 7, 1) = (0.15_dp, 0.0_dp); h%Data(2, 8, 1) = (0.15_dp, 0.0_dp)
        h%Data(3, 5, 1) = (0.05_dp, 0.0_dp); h%Data(4, 6, 1) = (0.05_dp, 0.0_dp)
        h%Data(7, 1, 1) = (0.15_dp, 0.0_dp); h%Data(8, 2, 1) = (0.15_dp, 0.0_dp)
        h%Data(5, 3, 1) = (0.05_dp, 0.0_dp); h%Data(6, 4, 1) = (0.05_dp, 0.0_dp)
        ! hopping em +a1 (a1, a2): up-up 0.3, dn-dn 0.5, flip up1-dn2 0.2, dn1-up2 0.1
        h%Data(1, 5, 2) = (0.3_dp, 0.0_dp); h%Data(2, 6, 2) = (0.3_dp, 0.0_dp)
        h%Data(3, 7, 2) = (0.5_dp, 0.0_dp); h%Data(4, 8, 2) = (0.5_dp, 0.0_dp)
        h%Data(1, 7, 2) = (0.2_dp, 0.0_dp); h%Data(2, 8, 2) = (0.2_dp, 0.0_dp)
        h%Data(3, 5, 2) = (0.1_dp, 0.0_dp); h%Data(4, 6, 2) = (0.1_dp, 0.0_dp)
        ! R = -a1: bloco (a2, a1) = transposto do de +a1
        h%Data(5, 1, 3) = (0.3_dp, 0.0_dp); h%Data(6, 2, 3) = (0.3_dp, 0.0_dp)
        h%Data(7, 3, 3) = (0.5_dp, 0.0_dp); h%Data(8, 4, 3) = (0.5_dp, 0.0_dp)
        h%Data(7, 1, 3) = (0.2_dp, 0.0_dp); h%Data(8, 2, 3) = (0.2_dp, 0.0_dp)
        h%Data(5, 3, 3) = (0.1_dp, 0.0_dp); h%Data(6, 4, 3) = (0.1_dp, 0.0_dp)

        call init_hs(s, [4, 4], reshape([0, 0, 0, 1, 0, 0, -1, 0, 0], [3, 3]), 'nc2 S')
        do i = 1, 8
            s%Data(i, i, 1) = (1.0_dp, 0.0_dp)
        end do

        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'nc2: build')
        if (ierr /= 0) return
        call check(all(ct%orbs_red == [2, 2]) .and. ct%norb_red == 4, 'nc2: orbs 4,4 -> 2,2')
        call check_pair_symmetry(h, 1.0e-9_dp, ierr, errmsg)
        call check(ierr == 0, 'nc2: pair symmetry with spin-flip coupling')

        call contract_hs(h, ct, hr, ierr, errmsg)
        call check(ierr == 0, 'nc2: contract')
        if (ierr /= 0) return
        call check_close(real(hr%Data(1, 1, 1)), 0.902_dp, 1.0e-13_dp, 'nc2: onsite up')
        call check_close(real(hr%Data(2, 2, 1)), 1.616_dp, 1.0e-13_dp, 'nc2: onsite dn')
        call check_close(real(hr%Data(1, 2, 1)), 0.278_dp, 1.0e-13_dp, 'nc2: onsite spin-flip')
        call check_close(real(hr%Data(1, 4, 1)), 0.15_dp * 0.85_dp, 1.0e-13_dp, 'nc2: inter-atom flip x12')
        call check_close(real(hr%Data(2, 3, 1)), 0.05_dp * 0.85_dp, 1.0e-13_dp, 'nc2: inter-atom flip x21')
        call check(abs(hr%Data(1, 3, 1)) < 1.0e-15_dp, 'nc2: no up-up inter-atom flip')
        call check_close(real(hr%Data(1, 3, 2)), 0.3_dp * 0.85_dp, 1.0e-13_dp, 'nc2: hopping up-up')
        call check_close(real(hr%Data(2, 4, 2)), 0.5_dp * 0.85_dp, 1.0e-13_dp, 'nc2: hopping dn-dn')
        call check_close(real(hr%Data(1, 4, 2)), 0.2_dp * 0.85_dp, 1.0e-13_dp, 'nc2: hopping flip 0.2')
        call check_close(real(hr%Data(2, 3, 2)), 0.1_dp * 0.85_dp, 1.0e-13_dp, 'nc2: hopping flip 0.1')

        call contract_hs(s, ct, sred, ierr, errmsg)
        call check(ierr == 0, 'nc2: contract S')
        if (ierr /= 0) return
        call lambda_min_grid(sred, [2, 2, 2], 1.0e-10_dp, smin, kmin, ierr, errmsg)
        call check(ierr == 0, 'nc2: lambda_min ok')
        call check_close(smin, 0.85_dp, 1.0e-12_dp, 'nc2: lambda_min = 0.85')
    end subroutine test_noncollinear_2atoms

    subroutine test_title_long()
        ! título original maior que 256 caracteres: a tag hscontract
        ! precisa sobreviver (o original é abreviado com '...')
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: h, hr
        character(len=8) :: labels(1)
        character(len=512) :: errmsg
        character(len=600) :: long
        integer :: ierr

        labels = ['Aa']
        inp%spin_mode = 'none'
        inp%basis(1) = 'Aa s1'
        inp%nrule = 0
        long = repeat('Lorem ipsum dolor sit amet ', 21)   ! 600 caracteres
        call init_hs(h, [1], reshape([0, 0, 0], [3, 1]), trim(long))
        h%Data(1, 1, 1) = (1.0_dp, 0.0_dp)
        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'title: build')
        if (ierr /= 0) return
        call contract_hs(h, ct, hr, ierr, errmsg)
        call check(ierr == 0, 'title: contract')
        if (ierr /= 0) return
        call check(len_trim(hr%Title) <= 256, 'title: fits in 256')
        call check(index(hr%Title, 'hscontract') > 0, 'title: tag preserved')
        call check(index(hr%Title, '...') > 0, 'title: original abbreviated')
        call check(index(hr%Title, 'Lorem ipsum') == 1, 'title: original head kept')
    end subroutine test_title_long

    ! ------------------------------------------------------------------
    ! Erros de construção/aplicação
    ! ------------------------------------------------------------------

    subroutine test_build_errors()
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        character(len=8) :: labels(2)
        character(len=512) :: errmsg
        real(dp) :: nan
        integer :: ierr

        ! espécie sem base declarada
        labels = ['Aa', 'Zz']
        inp%spin_mode = 'none'
        inp%basis(1) = 'Aa s2'
        inp%nrule = 0
        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr /= 0, 'build err: species without explicit basis')
        call check(index(errmsg, 'Zz') > 0, 'build err: message names the species')

        ! posto deficiente: duas funções de saída idênticas
        labels = ['Aa', 'Aa']
        inp%basis(1) = 'Aa s2'
        call set_rule(inp, 1, 'Aa', 0, 2, [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp])
        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr /= 0, 'build err: rank-deficient coefficients')
        call check(index(errmsg, 'rank') > 0, 'build err: message names the rank problem')

        ! coeficiente não finito (NaN bit pattern, sem operação de FP)
        nan = transfer(int(z'7FF8000000000000', int64), 1.0_dp)
        call set_rule(inp, 1, 'Aa', 0, 1, [nan, 1.0_dp])
        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr /= 0, 'build err: non-finite coefficients')

        ! ncoeff inconsistente na API direta
        call set_rule(inp, 1, 'Aa', 0, 1, [0.5_dp, 0.5_dp])
        inp%rule(1)%ncoeff = 1
        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr /= 0, 'build err: ncoeff mismatch (direct API)')

        ! spin_mode inválido na API direta
        inp%rule(1)%ncoeff = 2
        inp%spin_mode = 'diagonal'
        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr /= 0, 'build err: invalid spin_mode (direct API)')
    end subroutine test_build_errors

    subroutine test_apply_errors()
        type(contract_input_t) :: inp
        type(contraction_t) :: ct
        type(HSDataFormat) :: h, hr, h2
        character(len=8) :: labels(1)
        character(len=512) :: errmsg
        integer :: ierr

        labels = ['Aa']
        inp%spin_mode = 'none'
        inp%basis(1) = 'Aa s1'
        inp%nrule = 0
        call build_contraction(inp, labels, ct, ierr, errmsg)
        call check(ierr == 0, 'apply err: build')
        if (ierr /= 0) return

        ! Orbs incompatível com a base declarada
        call init_hs(h, [2], reshape([0, 0, 0], [3, 1]), 'x')
        h%Data(1, 1, 1) = (1.0_dp, 0.0_dp)
        h%Data(2, 2, 1) = (1.0_dp, 0.0_dp)
        call contract_hs(h, ct, hr, ierr, errmsg)
        call check(ierr /= 0, 'apply err: Orbs mismatch')
        call check(index(errmsg, 'basis') > 0, 'apply err: message names the basis')

        ! H/S incompatíveis
        call init_hs(h2, [1], reshape([0, 0, 0], [3, 1]), 'y')
        call check_hs_compatible(h, 'H', h2, 'S', ierr, errmsg)
        call check(ierr /= 0, 'apply err: incompatible H/S')

        ! quebra da estrutura algébrica: H(0) não hermitiana
        call init_hs(h2, [1], reshape([0, 0, 0], [3, 1]), 'z')
        h2%Data(1, 1, 1) = (1.0_dp, 0.5_dp)
        call check_pair_symmetry(h2, 1.0e-9_dp, ierr, errmsg)
        call check(ierr /= 0, 'apply err: non-hermitian H(0)')
        call check(index(errmsg, 'hermiticity') > 0, 'apply err: hermiticity message')

        ! quebra da estrutura algébrica: H(-R) /= H(R)^H
        call init_hs(h2, [1], reshape([0, 0, 0, 1, 0, 0, -1, 0, 0], [3, 3]), 'w')
        h2%Data(1, 1, 1) = (1.0_dp, 0.0_dp)
        h2%Data(1, 1, 2) = (0.5_dp, 0.0_dp)
        h2%Data(1, 1, 3) = (0.3_dp, 0.0_dp)
        call check_pair_symmetry(h2, 1.0e-9_dp, ierr, errmsg)
        call check(ierr /= 0, 'apply err: H(-R) /= H(R)^H')
        call check(index(errmsg, 'R = (') > 0, 'apply err: message names R')
    end subroutine test_apply_errors

    ! ------------------------------------------------------------------
    ! Caso sintético para o teste CLI (test_hscontract_cli.sh)
    ! ------------------------------------------------------------------

    subroutine write_case(dir)
        character(len=*), intent(in) :: dir
        type(HSDataFormat) :: h, s
        integer :: i, u
        real(dp) :: v

        ! 2 átomos: Aa s1p1 (4) e Bb s2 (2); R = 0, +-a1
        call init_hs(h, [4, 2], reshape([0, 0, 0, 1, 0, 0, -1, 0, 0], [3, 3]), 'synthetic H (eV)')
        s = h
        s%Title = 'synthetic S'
        do i = 1, 6
            v = 0.3_dp * i
            h%Data(i, i, 1) = cmplx(v, 0.0_dp, kind=dp)
            s%Data(i, i, 1) = (1.0_dp, 0.0_dp)
            h%Data(i, i, 2) = cmplx(-0.2_dp, 0.0_dp, kind=dp)
            h%Data(i, i, 3) = cmplx(-0.2_dp, 0.0_dp, kind=dp)
            s%Data(i, i, 2) = cmplx(-0.05_dp, 0.0_dp, kind=dp)
            s%Data(i, i, 3) = cmplx(-0.05_dp, 0.0_dp, kind=dp)
        end do
        h%has_meta = .true.
        h%fermi_ev = 0.0_dp
        h%energy_reference = 'absolute'
        call h%WriteFile(trim(dir) // '/case_H.hsio')
        call s%WriteFile(trim(dir) // '/case_S.hsio')

        ! fixtures collinear (mesma estrutura, mesmos sítios): H_up/H_dn
        block
            type(HSDataFormat) :: hu, hd
            hu = h; hd = h
            hu%Title = 'synthetic H_up (eV)'
            hd%Title = 'synthetic H_dn (eV)'
            do i = 1, 6
                hu%Data(i, i, 1) = hu%Data(i, i, 1) + (0.5_dp, 0.0_dp)
                hd%Data(i, i, 1) = hd%Data(i, i, 1) - (0.5_dp, 0.0_dp)
            end do
            call hu%WriteFile(trim(dir) // '/coll_H_up.hsio')
            call hd%WriteFile(trim(dir) // '/coll_H_dn.hsio')
        end block
        call s%WriteFile(trim(dir) // '/coll_S.hsio')

        open(newunit=u, file=trim(dir) // '/case.nml', status='replace', action='write')
        write(u, '(A)') '&STRUCTURE'
        write(u, '(A)') '    LatticeVectors(1,:) = 3.0, 0.0, 0.0'
        write(u, '(A)') '    LatticeVectors(2,:) = 0.0, 3.0, 0.0'
        write(u, '(A)') '    LatticeVectors(3,:) = 0.0, 0.0, 3.0'
        write(u, '(A)') '/'
        write(u, '(A)') '&SPECIES'
        write(u, '(A)') '2'
        write(u, '(A)') 'Aa 4'
        write(u, '(A)') 'Bb 2'
        write(u, '(A)') '/'
        write(u, '(A)') '&ATOMS'
        write(u, '(A)') '2'
        write(u, '(A)') 'Aa  0.0  0.0  0.0'
        write(u, '(A)') 'Bb  1.5  0.0  0.0'
        write(u, '(A)') '/'
        close(u)
    end subroutine write_case


    ! ------------------------------------------------------------------
    ! Caso de bandas (test_hscontract_bands.sh): modelo 1D nao ortogonal,
    ! multi-zeta, shell d totalmente desacoplada das shells s/p (H e S), com
    ! energias onsite altas. A poda d -> 0 e EXATA para as bandas de s/p,
    ! logo a tolerancia de comparacao e numerica (nao um limite físico
    ! arbitrario). Janela E_F = 0 +- 3 eV: contem as 5 bandas s/p; as 5
    ! bandas d ficam fora (~6.3-8.1 eV).
    ! ------------------------------------------------------------------

    !> Dois modelos 1D nao ortogonais multi-zeta ('Aa s2p1d1', 10 orbitais,
    !> R = 0, +-a1, +-2a1), shell d totalmente desacoplada das shells s/p com
    !> energias onsite altas (~7 eV) e bandas fora da janela E_F +- 3 eV:
    !> - 'model'    : onsite do zeta 2 do s em +0.8 eV (banda dentro da
    !>               janela; usada nos casos EXATOS: T = I e poda d -> 0);
    !> - 'appmodel' : onsite do zeta 2 do s em +5.0 eV (banda fora da
    !>               janela; usada na CONTRACAO APROXIMADA s2 -> 1 com pesos
    !>               [0.9, 0.1] -- nao normalizados). A janela fica com os
    !>               4 estados relevantes (s1 + 3p) nos dois lados.
    subroutine write_bands_case(dir)
        character(len=*), intent(in) :: dir
        integer :: u

        call write_bands_model(trim(dir), 'model', 0.8_dp)
        call write_bands_model(trim(dir), 'appmodel', 5.0_dp)

        open(newunit=u, file=trim(dir) // '/model.nml', status='replace', action='write')
        write(u, '(A)') '&STRUCTURE'
        write(u, '(A)') '    LatticeVectors(1,:) = 3.0, 0.0, 0.0'
        write(u, '(A)') '    LatticeVectors(2,:) = 0.0, 3.0, 0.0'
        write(u, '(A)') '    LatticeVectors(3,:) = 0.0, 0.0, 3.0'
        write(u, '(A)') '/'
        write(u, '(A)') '&SPECIES'
        write(u, '(A)') '1'
        write(u, '(A)') 'Aa 10'
        write(u, '(A)') '/'
        write(u, '(A)') '&ATOMS'
        write(u, '(A)') '1'
        write(u, '(A)') 'Aa  0.0  0.0  0.0'
        write(u, '(A)') '/'
        close(u)

        open(newunit=u, file=trim(dir) // '/contract_id.nml', status='replace', action='write')
        write(u, '(A)') '&HSCONTRACT'
        write(u, '(A)') "    spin_mode = 'none'"
        write(u, '(A)') "    h_file = 'model_H.hsio'"
        write(u, '(A)') "    s_file = 'model_S.hsio'"
        write(u, '(A)') "    structure_file = 'model.nml'"
        write(u, '(A)') "    output_prefix = 'rid'"
        write(u, '(A)') "    basis = 'Aa s2p1d1'"
        write(u, '(A)') '    check_kgrid = 2, 2, 2'
        write(u, '(A)') '    tol_s_min = 1.0e-10'
        write(u, '(A)') '/'
        close(u)

        open(newunit=u, file=trim(dir) // '/contract_prune.nml', status='replace', action='write')
        write(u, '(A)') '&HSCONTRACT'
        write(u, '(A)') "    spin_mode = 'none'"
        write(u, '(A)') "    h_file = 'model_H.hsio'"
        write(u, '(A)') "    s_file = 'model_S.hsio'"
        write(u, '(A)') "    structure_file = 'model.nml'"
        write(u, '(A)') "    output_prefix = 'rpr'"
        write(u, '(A)') "    basis = 'Aa s2p1d1'"
        write(u, '(A)') '    check_kgrid = 2, 2, 2'
        write(u, '(A)') '    tol_s_min = 1.0e-10'
        write(u, '(A)') '/'
        write(u, '(A)') '&RADIAL_RULE'
        write(u, '(A)') "    species = 'Aa'"
        write(u, '(A)') '    l = 2'
        write(u, '(A)') '    n_out = 0'
        write(u, '(A)') '/'
        close(u)

        ! contração aproximada: poda d exata + s2 -> 1 com pesos nao triviais
        ! nos DOIS zetas (nao normalizados: c.c = 0.82)
        open(newunit=u, file=trim(dir) // '/contract_appr.nml', status='replace', action='write')
        write(u, '(A)') '&HSCONTRACT'
        write(u, '(A)') "    spin_mode = 'none'"
        write(u, '(A)') "    h_file = 'appmodel_H.hsio'"
        write(u, '(A)') "    s_file = 'appmodel_S.hsio'"
        write(u, '(A)') "    structure_file = 'model.nml'"
        write(u, '(A)') "    output_prefix = 'rap'"
        write(u, '(A)') "    basis = 'Aa s2p1d1'"
        write(u, '(A)') '    check_kgrid = 2, 2, 2'
        write(u, '(A)') '    tol_s_min = 1.0e-10'
        write(u, '(A)') '/'
        write(u, '(A)') '&RADIAL_RULE'
        write(u, '(A)') "    species = 'Aa'"
        write(u, '(A)') '    l = 0'
        write(u, '(A)') '    n_out = 1'
        write(u, '(A)') '    coefficients = 0.9, 0.1'
        write(u, '(A)') '/'
        write(u, '(A)') '&RADIAL_RULE'
        write(u, '(A)') "    species = 'Aa'"
        write(u, '(A)') '    l = 2'
        write(u, '(A)') '    n_out = 0'
        write(u, '(A)') '/'
        close(u)
    end subroutine write_bands_case

    !> Escreve <prefix>_H.hsio / <prefix>_S.hsio do modelo 1D; e2 é a
    !> energia onsite do zeta 2 do s (0.8 eV no modelo exato, 5.0 eV no
    !> modelo da contração aproximada).
    subroutine write_bands_model(dir, prefix, e2)
        character(len=*), intent(in) :: dir, prefix
        real(dp), intent(in) :: e2
        type(HSDataFormat) :: h, s
        integer :: i

        call init_hs(h, [10], reshape([0, 0, 0, 1, 0, 0, -1, 0, 0, 2, 0, 0, -2, 0, 0], [3, 5]), &
                     'synthetic bands H (eV)')
        s = h
        s%Title = 'synthetic bands S'
        h%has_meta = .true.
        h%fermi_ev = 0.0_dp
        h%energy_reference = 'absolute'

        ! R = 0
        s%Data(1, 1, 1) = (1.0_dp, 0.0_dp);  s%Data(2, 2, 1) = (0.9_dp, 0.0_dp)
        s%Data(1, 2, 1) = (0.3_dp, 0.0_dp);  s%Data(2, 1, 1) = (0.3_dp, 0.0_dp)
        do i = 3, 5
            s%Data(i, i, 1) = (1.0_dp, 0.0_dp)
            s%Data(1, i, 1) = (0.1_dp, 0.0_dp);  s%Data(i, 1, 1) = (0.1_dp, 0.0_dp)
            s%Data(2, i, 1) = (0.05_dp, 0.0_dp); s%Data(i, 2, 1) = (0.05_dp, 0.0_dp)
        end do
        do i = 6, 10
            s%Data(i, i, 1) = (1.0_dp, 0.0_dp)
        end do
        h%Data(1, 1, 1) = (-1.0_dp, 0.0_dp); h%Data(2, 2, 1) = cmplx(e2, 0.0_dp, kind=dp)
        h%Data(1, 2, 1) = (0.05_dp, 0.0_dp); h%Data(2, 1, 1) = (0.05_dp, 0.0_dp)
        h%Data(3, 3, 1) = (0.4_dp, 0.0_dp);  h%Data(4, 4, 1) = (0.5_dp, 0.0_dp)
        h%Data(5, 5, 1) = (0.6_dp, 0.0_dp)
        do i = 3, 5
            h%Data(1, i, 1) = (0.1_dp, 0.0_dp);  h%Data(i, 1, 1) = (0.1_dp, 0.0_dp)
            h%Data(2, i, 1) = (-0.05_dp, 0.0_dp); h%Data(i, 2, 1) = (-0.05_dp, 0.0_dp)
        end do
        do i = 6, 10
            h%Data(i, i, 1) = cmplx(7.0_dp + 0.1_dp * (i - 6), 0.0_dp, kind=dp)
        end do

        ! R = +-a1 (indices 2, 3) e R = +-2a1 (indices 4, 5)
        call set_hopping(h, s, 2, -0.4_dp, -0.15_dp, 0.02_dp, -0.2_dp, 0.05_dp, -0.02_dp, -0.3_dp, &
                         0.05_dp, 0.02_dp, 0.03_dp, 0.01_dp, 0.01_dp)
        s%Data(:, :, 3) = s%Data(:, :, 2)
        h%Data(:, :, 3) = h%Data(:, :, 2)
        call set_hopping(h, s, 4, -0.08_dp, -0.05_dp, 0.0_dp, -0.04_dp, 0.0_dp, 0.0_dp, -0.05_dp, &
                         0.01_dp, 0.008_dp, 0.01_dp, 0.0_dp, 0.005_dp)
        s%Data(:, :, 5) = s%Data(:, :, 4)
        h%Data(:, :, 5) = h%Data(:, :, 4)

        call h%WriteFile(trim(dir) // '/' // trim(prefix) // '_H.hsio')
        call s%WriteFile(trim(dir) // '/' // trim(prefix) // '_S.hsio')
    end subroutine write_bands_model

    ! Preenche o bloco H/S do vizinho ir (hopping 1D, diagonal em cada setor)
    subroutine set_hopping(h, s, ir, t_s1, t_s2, t_s12, t_p, t_sp_s1, t_sp_s2, t_d, &
                           ov_s1, ov_s2, ov_p, ov_sp, ov_d)
        type(HSDataFormat), intent(inout) :: h, s
        integer, intent(in) :: ir
        real(dp), intent(in) :: t_s1, t_s2, t_s12, t_p, t_sp_s1, t_sp_s2, t_d
        real(dp), intent(in) :: ov_s1, ov_s2, ov_p, ov_sp, ov_d
        integer :: i
        h%Data(1, 1, ir) = cmplx(t_s1, 0.0_dp, kind=dp)
        h%Data(2, 2, ir) = cmplx(t_s2, 0.0_dp, kind=dp)
        h%Data(1, 2, ir) = cmplx(t_s12, 0.0_dp, kind=dp)
        h%Data(2, 1, ir) = cmplx(t_s12, 0.0_dp, kind=dp)
        do i = 3, 5
            h%Data(i, i, ir) = cmplx(t_p, 0.0_dp, kind=dp)
            h%Data(1, i, ir) = cmplx(t_sp_s1, 0.0_dp, kind=dp)
            h%Data(i, 1, ir) = cmplx(t_sp_s1, 0.0_dp, kind=dp)
            h%Data(2, i, ir) = cmplx(t_sp_s2, 0.0_dp, kind=dp)
            h%Data(i, 2, ir) = cmplx(t_sp_s2, 0.0_dp, kind=dp)
        end do
        do i = 6, 10
            h%Data(i, i, ir) = cmplx(t_d, 0.0_dp, kind=dp)
        end do
        s%Data(1, 1, ir) = cmplx(ov_s1, 0.0_dp, kind=dp)
        s%Data(2, 2, ir) = cmplx(ov_s2, 0.0_dp, kind=dp)
        do i = 3, 5
            s%Data(i, i, ir) = cmplx(ov_p, 0.0_dp, kind=dp)
            s%Data(1, i, ir) = cmplx(ov_sp, 0.0_dp, kind=dp)
            s%Data(i, 1, ir) = cmplx(ov_sp, 0.0_dp, kind=dp)
            s%Data(2, i, ir) = cmplx(ov_sp, 0.0_dp, kind=dp)
            s%Data(i, 2, ir) = cmplx(ov_sp, 0.0_dp, kind=dp)
        end do
        do i = 6, 10
            s%Data(i, i, ir) = cmplx(ov_d, 0.0_dp, kind=dp)
        end do
    end subroutine set_hopping

    ! ------------------------------------------------------------------
    ! bandcmp: compara dois arquivos .bands.bin do LowFoldIn (binário bruto:
    ! int32 nOrbs, int32 nk, real64 x(nk), real64 E(nOrbs, nk)). Por k, reúne
    ! as energias na janela [emin, emax], ordena os dois conjuntos e compara
    ! posição a posição; métricas max |dE| e RMS; falha se max > tol.
    ! ------------------------------------------------------------------

    subroutine run_bandcmp()
        character(len=256) :: f1, f2
        character(len=64) :: cbuf
        real(dp) :: emin, emax, tol
        integer :: ios, u1, u2, nb1, nb2, nk1, nk2, ik, n1, n2, ntot, i
        real(dp), allocatable :: x(:), e1(:, :), e2(:, :), w1(:), w2(:)
        real(dp) :: dmax, drms, sumsq, dtmp, lo1, hi1, lo2, hi2

        call get_command_argument(2, f1)
        call get_command_argument(3, f2)
        call get_command_argument(4, cbuf)
        read(cbuf, *, iostat=ios) emin
        if (ios == 0) call get_command_argument(5, cbuf)
        if (ios == 0) read(cbuf, *, iostat=ios) emax
        if (ios == 0) call get_command_argument(6, cbuf)
        if (ios == 0) read(cbuf, *, iostat=ios) tol
        if (len_trim(f1) == 0 .or. len_trim(f2) == 0 .or. ios /= 0 .or. emax <= emin .or. tol <= 0.0_dp) then
            print '(A)', 'usage: test_basis_contract bandcmp F1.bands.bin F2.bands.bin EMIN EMAX TOL'
            stop 1
        end if

        open(newunit=u1, file=trim(f1), form='unformatted', access='stream', status='old', iostat=ios)
        if (ios /= 0) then; print '(A)', 'bandcmp: cannot open ' // trim(f1); stop 1; end if
        open(newunit=u2, file=trim(f2), form='unformatted', access='stream', status='old', iostat=ios)
        if (ios /= 0) then; print '(A)', 'bandcmp: cannot open ' // trim(f2); stop 1; end if
        read(u1, iostat=ios) nb1, nk1
        if (ios /= 0) then; print '(A)', 'bandcmp: truncated ' // trim(f1); stop 1; end if
        read(u2, iostat=ios) nb2, nk2
        if (ios /= 0) then; print '(A)', 'bandcmp: truncated ' // trim(f2); stop 1; end if
        if (nk1 /= nk2) then
            print '(A,I0,A,I0)', 'bandcmp: k-point count mismatch ', nk1, ' vs ', nk2
            stop 1
        end if
        allocate(x(nk1), e1(nb1, nk1), e2(nb2, nk2), w1(max(nb1, nb2)), w2(max(nb1, nb2)))
        read(u1, iostat=ios) x
        if (ios == 0) read(u1, iostat=ios) e1
        read(u2, iostat=ios) x
        if (ios == 0) read(u2, iostat=ios) e2
        if (ios /= 0) then; print '(A)', 'bandcmp: truncated data'; stop 1; end if
        close(u1); close(u2)

        dmax = 0.0_dp
        sumsq = 0.0_dp
        ntot = 0
        lo1 = huge(1.0_dp); hi1 = -huge(1.0_dp); lo2 = huge(1.0_dp); hi2 = -huge(1.0_dp)
        do ik = 1, nk1
            call window_sort(e1(:, ik), emin, emax, w1, n1)
            call window_sort(e2(:, ik), emin, emax, w2, n2)
            if (n1 /= n2) then
                print '(A,I0,A,I0,A,I0)', 'bandcmp: band count in window at k #', ik, ': ', n1, ' vs ', n2
                stop 1
            end if
            do i = 1, n1
                dtmp = abs(w1(i) - w2(i))
                dmax = max(dmax, dtmp)
                sumsq = sumsq + dtmp * dtmp
            end do
            ntot = ntot + n1
            if (n1 > 0) then
                lo1 = min(lo1, w1(1)); hi1 = max(hi1, w1(n1))
                lo2 = min(lo2, w2(1)); hi2 = max(hi2, w2(n1))
            end if
        end do
        if (ntot == 0) then
            print '(A)', 'bandcmp: no bands in window'
            stop 1
        end if
        print '(A,I0,A,F8.3,A,F8.3,A)', '[bandcmp] k-points: ', nk1, ', window: [', emin, ', ', emax, '] eV'
        print '(A,I0)', '[bandcmp] bands in window per k: ', n1
        if (nb1 == nb2) print '(A)', '[bandcmp] (same band count in both files)'
        print '(A,ES10.3,A,ES10.3,A)', '[bandcmp] max |dE| = ', dmax, ' eV, RMS = ', &
            sqrt(sumsq / real(ntot, dp)), ' eV'
        print '(A,2F10.6)', '[bandcmp] extrema full : ', lo1, hi1
        print '(A,2F10.6)', '[bandcmp] extrema red  : ', lo2, hi2
        if (dmax > tol) then
            print '(A,ES10.3,A,ES10.3)', 'bandcmp: FAIL max |dE| = ', dmax, ' > tol = ', tol
            stop 1
        end if
        print '(A)', 'bandcmp: PASS'
    end subroutine run_bandcmp

    !> Energias de e(:) na janela [emin, emax], ordenadas (inserção), em w(1:n).
    subroutine window_sort(e, emin, emax, w, n)
        real(dp), intent(in) :: e(:), emin, emax
        real(dp), intent(out) :: w(:)
        integer, intent(out) :: n
        integer :: i, j
        real(dp) :: v
        n = 0
        do i = 1, size(e)
            if (e(i) < emin .or. e(i) > emax) cycle
            v = e(i)
            j = n
            do while (j >= 1)
                if (w(j) <= v) exit
                w(j + 1) = w(j)
                j = j - 1
            end do
            w(j + 1) = v
            n = n + 1
        end do
    end subroutine window_sort

end program test_basis_contract
