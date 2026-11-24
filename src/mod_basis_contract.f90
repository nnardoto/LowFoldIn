! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Contração da base LCAO de H(R)/S(R) em arquivos .hsio (utilitário hscontract).
!>
!> Dada uma transformação T bloco-diagonal por átomo com coeficientes reais,
!> a base contratida é
!>     H_r(R) = T† H(R) T,   S_r(R) = T† S(R) T    para todo R.
!> A base contratida NÃO é ortogonal em geral (T† S T ≠ I mesmo com pesos
!> ortonormais); nenhum peso é normalizado automaticamente.
!>
!> Regras radiais por espécie (&RADIAL_RULE: species, l, n_out, coefficients):
!>   - shell l sem regra: preservada (identidade);
!>   - n_out = 0: a shell l inteira é removida;
!>   - 1 <= n_out <= nzeta(l): nzeta -> n_out funções, com coefficients
!>     ordenados por função de saída (grupos de nzeta pesos); a mesma
!>     contração radial vale para todos os m — uma regra sempre se restringe
!>     a um único l e nunca seleciona ou mistura m isoladamente.
!>
!> A ordem dos orbitais é a do OpenMX (src/mod_sk_table.f90:19-26): l
!> crescente, zeta crescente, m crescente; não colinear: átomo-major
!> [orbitais espaciais up | orbitais espaciais down] por átomo. A base de
!> cada espécie deve ser declarada explicitamente (nada é inferido de nomes
!> PAO); nspatial é validado contra Orbs do HSIO (2*nspatial no não colinear).
!>
!> A biblioteca nunca para o programa: erros voltam por ierr/errmsg.
!> A leitura/escrita densa do HSIO (HSDataFormat) é aceita e documentada.
module basis_contract
    use types, only: dp
    use sk_table, only: species_basis_t, shell_t, parse_basis, basis_shell_index
    use HS_io, only: HSDataFormat
    use core, only: compute_Hk_at_kpoint
    use zlinalg, only: herm_check, herm_eigvals
    use, intrinsic :: iso_fortran_env, only: iostat_end, int64
    implicit none
    private

    ! ==================================================================
    ! CONSTANTES PÚBLICAS
    ! ==================================================================
    integer, parameter, public :: BC_PATH_LEN     = 512   ! buffers de caminho
    integer, parameter, public :: BC_MAX_SPECIES = 32     ! entradas de basis
    integer, parameter, public :: BC_MAX_COEFF   = 64     ! coefficients por regra
    integer, parameter, public :: BC_SPEC_STR_LEN = 64    ! 'Especie s3p2d2'
    integer, parameter, public :: BC_LABEL_LEN   = 16     ! species de &RADIAL_RULE

    !> Tolerância relativa padrão do cheque H(-R) = H(R)† (escala = max|H|).
    real(dp), parameter, public :: BC_TOL_PAIR = 1.0e-9_dp
    !> Posto completo: w_min(C C^T) > BC_RANK_RTOL * w_max(C C^T).
    real(dp), parameter :: BC_RANK_RTOL = 1.0e-12_dp
    !> Marca posições não usadas no buffer coefficients do namelist.
    real(dp), parameter :: COEFF_SENTINEL = -huge(1.0_dp)

    ! ==================================================================
    ! TIPOS
    ! ==================================================================

    !> Uma regra radial: contração das nzeta(l) funções radiais da espécie
    !> para n_out (0 = remover a shell). coeff(1:ncoeff) com ncoeff =
    !> n_out*nzeta, ordenado por função de saída (grupos de nzeta pesos).
    type, public :: radial_rule_t
        character(len=BC_LABEL_LEN) :: species = ''
        integer          :: l      = -1
        integer          :: n_out  = -1
        integer          :: ncoeff = 0
        real(dp)         :: coeff(BC_MAX_COEFF) = 0.0_dp
    end type radial_rule_t

    !> Conteúdo de contract.nml (&HSCONTRACT + blocos &RADIAL_RULE).
    type, public :: contract_input_t
        character(len=16)          :: spin_mode = 'none'  ! 'none'|'collinear'|'noncollinear'
        character(len=BC_PATH_LEN) :: h_file = '', h_up_file = '', h_dn_file = ''
        character(len=BC_PATH_LEN) :: s_file = '', structure_file = '', output_prefix = ''
        character(len=BC_SPEC_STR_LEN) :: basis(BC_MAX_SPECIES) = ''
        integer          :: check_kgrid(3) = 2
        real(dp)         :: tol_s_min = 1.0e-10_dp
        integer          :: nrule = 0
        type(radial_rule_t), allocatable :: rule(:)
    end type contract_input_t

    !> Bloco por átomo: T espacial (real) e a cópia complexa usada na
    !> aplicação, já em dobro no não colinear (blkdiag(T, T), [up|dn]).
    type :: atom_ct_t
        real(dp), allocatable :: t(:, :)      ! (nspatial x nspatial_red)
        complex(dp), allocatable :: tc(:, :) ! (orbs_full x orbs_red), com spinor
    end type atom_ct_t

    !> Contração completa de um sistema: T bloco-diagonal por átomo.
    type, public :: contraction_t
        integer :: nat = 0
        integer :: spinor = 1                       ! 2 no não colinear
        integer :: norb_full = 0, norb_red = 0
        integer, allocatable :: orbs_full(:), orbs_red(:)   ! (nat)
        integer, allocatable :: off_full(:), off_red(:)     ! (nat+1)
        type(atom_ct_t), allocatable :: atom(:)
    end type contraction_t

    ! ==================================================================
    ! API PÚBLICA
    ! ==================================================================
    public :: read_contract_input
    public :: build_contraction, contract_hs
    public :: check_hs_compatible, check_pair_symmetry
    public :: lambda_min_grid

contains

    ! ==================================================================
    ! LEITURA DE contract.nml
    ! ==================================================================

    !> Lê &HSCONTRACT (primeiro e único grupo) seguido de qualquer número de
    !> blocos &RADIAL_RULE. Valida: spin_mode, arquivos exigidos por modo de
    !> spin, basis (uma entrada explícita por espécie, sem duplicatas),
    !> check_kgrid >= 1, tol_s_min > 0 e, por regra, species/l/n_out/
    !> coefficients contra a base declarada (nzeta, n_out <= nzeta,
    !> ncoeff = n_out*nzeta, coeficientes finitos, sem regra duplicada).
    !> Nunca para o programa; ierr = 1 e errmsg descreve o problema.
    subroutine read_contract_input(filename, inp, ierr, errmsg)
        character(len=*),       intent(in)  :: filename
        type(contract_input_t), intent(out) :: inp
        integer,                intent(out) :: ierr
        character(len=*),        intent(out) :: errmsg

        ! &HSCONTRACT
        character(len=16)          :: spin_mode
        character(len=BC_PATH_LEN) :: h_file, h_up_file, h_dn_file, s_file, structure_file, output_prefix
        character(len=BC_SPEC_STR_LEN) :: basis(BC_MAX_SPECIES)
        integer :: check_kgrid(3)
        real(dp) :: tol_s_min
        namelist /HSCONTRACT/ spin_mode, h_file, h_up_file, h_dn_file, s_file, &
                              structure_file, output_prefix, basis, check_kgrid, tol_s_min

        ! &RADIAL_RULE
        character(len=BC_LABEL_LEN) :: species
        integer :: l, n_out
        real(dp) :: coefficients(BC_MAX_COEFF)
        namelist /RADIAL_RULE/ species, l, n_out, coefficients

        type(species_basis_t) :: b(BC_MAX_SPECIES)
        type(radial_rule_t), allocatable :: rules(:)
        character(len=512) :: msg
        character(len=:), allocatable :: msg2
        logical :: ok
        integer :: iu, ios, iblock, ir, nb, ib, nzl, nc, i

        ierr = 0
        errmsg = ''
        call check_header_position(filename, ok, ierr, errmsg)
        if (ierr /= 0) return

        open(newunit=iu, file=filename, status='old', action='read', iostat=ios, iomsg=msg)
        if (ios /= 0) then
            ierr = 1
            errmsg = 'cannot open contract file "' // trim(filename) // '": ' // trim(msg)
            return
        end if

        spin_mode = ''
        h_file = ''; h_up_file = ''; h_dn_file = ''
        s_file = ''; structure_file = ''; output_prefix = ''
        basis = ''
        check_kgrid = 2
        tol_s_min = 1.0e-10_dp
        read(iu, nml=HSCONTRACT, iostat=ios, iomsg=msg)
        if (ios /= 0) then
            close(iu)
            ierr = 1
            errmsg = 'contract file "' // trim(filename) // '": error reading &HSCONTRACT: ' // trim(msg)
            return
        end if

        ! -- campos escalares ------------------------------------------
        inp%spin_mode = lower(adjustl(spin_mode))
        inp%h_file = h_file; inp%h_up_file = h_up_file; inp%h_dn_file = h_dn_file
        inp%s_file = s_file; inp%structure_file = structure_file
        inp%output_prefix = output_prefix
        inp%check_kgrid = check_kgrid
        if (tol_s_min <= 0.0_dp) then
            close(iu)
            ierr = 1
            errmsg = 'contract file: tol_s_min must be > 0'
            return
        end if
        inp%tol_s_min = tol_s_min

        select case (trim(inp%spin_mode))
        case ('none', 'noncollinear')
            if (len_trim(inp%h_file) == 0) then
                ierr = 1
                errmsg = "contract file: h_file is required for spin_mode='" // trim(inp%spin_mode) // "'"
            end if
        case ('collinear')
            if (len_trim(inp%h_up_file) == 0 .or. len_trim(inp%h_dn_file) == 0) then
                ierr = 1
                errmsg = "contract file: h_up_file and h_dn_file are required for spin_mode='collinear'"
            end if
        case default
            ierr = 1
            errmsg = "contract file: invalid spin_mode '" // trim(spin_mode) // &
                     "' (expected none, collinear or noncollinear)"
        end select
        if (ierr /= 0) then
            close(iu); return
        end if
        if (len_trim(inp%s_file) == 0) then
            ierr = 1; errmsg = 'contract file: s_file is required'
        else if (len_trim(inp%structure_file) == 0) then
            ierr = 1; errmsg = 'contract file: structure_file is required'
        else if (len_trim(inp%output_prefix) == 0) then
            ierr = 1; errmsg = 'contract file: output_prefix is required'
        else if (any(inp%check_kgrid < 1)) then
            ierr = 1; errmsg = 'contract file: check_kgrid components must be >= 1'
        else if (inp%tol_s_min <= 0.0_dp) then
            ierr = 1; errmsg = 'contract file: tol_s_min must be > 0'
        end if
        if (ierr /= 0) then
            close(iu); return
        end if

        ! -- basis: uma entrada explícita por espécie -------------------
        nb = count_leading_str(basis)
        if (nb < 0) then
            close(iu); ierr = 1; errmsg = 'contract file: gap in basis entries'; return
        end if
        if (nb == 0) then
            close(iu)
            ierr = 1
            errmsg = "contract file: basis is required (e.g. basis = 'Bi s3p2d2', 'Sb s3p2d2')"
            return
        end if
        do ib = 1, nb
            call parse_basis(basis(ib), b(ib), ierr, msg2)
            if (ierr /= 0) then
                close(iu); errmsg = 'contract file: basis entry ' // itoa(ib) // ': ' // msg2; return
            end if
            do i = 1, ib - 1
                if (trim(b(i)%name) == trim(b(ib)%name)) then
                    close(iu)
                    ierr = 1
                    errmsg = 'contract file: duplicate basis entry for species ''' // trim(b(ib)%name) // ''''
                    return
                end if
            end do
        end do
        inp%basis = basis

        ! -- blocos &RADIAL_RULE ----------------------------------------
        allocate(rules(0))
        iblock = 0
        do
            iblock = iblock + 1
            species = ''
            l = -huge(0)
            n_out = -huge(0)
            coefficients = COEFF_SENTINEL
            read(iu, nml=RADIAL_RULE, iostat=ios, iomsg=msg)
            if (ios == iostat_end) exit
            if (ios /= 0) then
                close(iu)
                ierr = 1
                errmsg = 'contract file: error reading &RADIAL_RULE block ' // itoa(iblock) // ': ' // trim(msg)
                return
            end if
            call validate_rule_block(species, l, n_out, coefficients, b, nb, rules, ir, ierr, errmsg)
            if (ierr /= 0) then
                close(iu)
                errmsg = 'contract file: &RADIAL_RULE block ' // itoa(iblock) // ': ' // errmsg
                return
            end if
        end do
        close(iu)

        inp%nrule = size(rules)
        if (inp%nrule > 0) inp%rule = rules
    end subroutine read_contract_input

    !> Valida um bloco &RADIAL_RULE contra as bases já analisadas e contra as
    !> regras anteriores; devolve em ir o índice (1-based) da regra na lista.
    subroutine validate_rule_block(species, l, n_out, coefficients, b, nb, rules, ir, ierr, errmsg)
        character(len=BC_LABEL_LEN), intent(in)    :: species
        integer,                     intent(in)    :: l, n_out
        real(dp),                    intent(in)    :: coefficients(:)
        type(species_basis_t),       intent(in)    :: b(:)
        integer,                     intent(in)    :: nb
        type(radial_rule_t), allocatable, intent(inout) :: rules(:)
        integer,                     intent(out)   :: ir
        integer,                     intent(out)   :: ierr
        character(len=*),            intent(out)   :: errmsg

        integer :: ib, nzl, nc, k

        ierr = 0
        errmsg = ''
        ir = 0
        if (len_trim(species) == 0) then
            ierr = 1; errmsg = 'species is required'; return
        end if
        if (l < 0 .or. l > 3) then
            ierr = 1; errmsg = 'l must be 0..3 (got ' // itoa(l) // ')'; return
        end if
        if (n_out < 0) then
            ierr = 1; errmsg = 'n_out must be >= 0 (got ' // itoa(n_out) // ')'; return
        end if
        ib = 0
        do k = 1, nb
            if (trim(b(k)%name) == trim(adjustl(species))) ib = k
        end do
        if (ib == 0) then
            ierr = 1
            errmsg = "species '" // trim(species) // "' has no basis entry"
            return
        end if
        nzl = 0
        do k = 1, b(ib)%nshell
            if (b(ib)%shell(k)%l == l) nzl = nzl + 1
        end do
        if (nzl == 0) then
            ierr = 1
            errmsg = "species '" // trim(species) // "' has no " // lchar(l) // ' shell'
            return
        end if
        if (n_out > nzl) then
            ierr = 1
            errmsg = 'n_out = ' // itoa(n_out) // ' exceeds nzeta(' // lchar(l) // ') = ' // itoa(nzl) &
                     // " of species '" // trim(species) // "'"
            return
        end if
        nc = count_leading_coeff(coefficients)
        if (nc < 0) then
            ierr = 1; errmsg = 'gap in coefficients list'; return
        end if
        if (nc /= n_out * nzl) then
            ierr = 1
            errmsg = 'expected n_out*nzeta = ' // itoa(n_out * nzl) // ' coefficients, got ' // itoa(nc)
            return
        end if
        if (nc > 0) then
            if (.not. reals_finite(coefficients(1:nc))) then
                ierr = 1; errmsg = 'coefficients must be finite'
                return
            end if
        end if
        do k = 1, size(rules)
            if (trim(rules(k)%species) == trim(adjustl(species)) .and. rules(k)%l == l) then
                ierr = 1
                errmsg = 'duplicate rule for species ''' // trim(species) // ''' and l = ' // itoa(l)
                return
            end if
        end do

        ir = size(rules) + 1
        rules = [rules, radial_rule_t(trim(adjustl(species)), l, n_out)]
        rules(ir)%ncoeff = nc
        if (nc > 0) rules(ir)%coeff(1:nc) = coefficients(1:nc)
    end subroutine validate_rule_block

    !> Pré-varredura de texto: &HSCONTRACT deve ser o primeiro grupo do
    !> namelist e aparecer exatamente uma vez; os únicos grupos aceitos são
    !> &HSCONTRACT e &RADIAL_RULE. Qualquer outro grupo (p. ex. um
    !> &RADIAL_RULES com erro de digitação, que a leitura de namelist sozinha
    !> pularia em silêncio, produzindo contração identidade) é rejeitado aqui.
    subroutine check_header_position(filename, ok, ierr, errmsg)
        character(len=*), intent(in)  :: filename
        logical,          intent(out) :: ok
        integer,          intent(out) :: ierr
        character(len=*), intent(out) :: errmsg
        integer :: iu, ios, ngroup, nheader, i
        character(len=1024) :: line, tok
        character(len=512) :: msg

        ok = .false.
        ierr = 0
        errmsg = ''
        open(newunit=iu, file=filename, status='old', action='read', iostat=ios, iomsg=msg)
        if (ios /= 0) then
            ierr = 1
            errmsg = 'cannot open contract file "' // trim(filename) // '": ' // trim(msg)
            return
        end if
        ngroup = 0
        nheader = 0
        do
            read(iu, '(A)', iostat=ios) line
            if (ios /= 0) exit
            do i = 1, len(line)
                if (line(i:i) == char(9)) line(i:i) = ' '
            end do
            line = adjustl(line)
            if (line(1:1) /= '&') cycle
            ngroup = ngroup + 1
            tok = upper(first_token(line))
            select case (trim(tok))
            case ('&HSCONTRACT')
                nheader = nheader + 1
                if (ngroup /= 1) then
                    close(iu)
                    ierr = 1
                    errmsg = 'contract file: &HSCONTRACT must be the first namelist group'
                    return
                end if
            case ('&RADIAL_RULE')
                ! aceito em qualquer número (lidos em sequência depois do header)
            case default
                close(iu)
                ierr = 1
                errmsg = 'contract file: unexpected namelist group ''' // trim(tok) // &
                         ''' (allowed: &HSCONTRACT, &RADIAL_RULE)'
                return
            end select
        end do
        close(iu)
        if (nheader /= 1) then
            ierr = 1
            errmsg = 'contract file: expected exactly one &HSCONTRACT group, found ' // itoa(nheader)
            return
        end if
        ok = .true.
    end subroutine check_header_position

    ! ==================================================================
    ! CONSTRUÇÃO DA CONTRAÇÃO
    ! ==================================================================

    !> Constrói a contração (T bloco-diagonal por átomo) a partir da entrada
    !> e dos rótulos de espécie por átomo (ordem do arquivo de estrutura =
    !> ordem HSIO). Exige base explícita para toda espécie; valida n_out,
    !> ncoeff, coeficientes finitos e posto completo (BC_RANK_RTOL relativo).
    !> Nunca para o programa.
    subroutine build_contraction(inp, labels, ct, ierr, errmsg)
        type(contract_input_t), intent(in)  :: inp
        character(len=*),       intent(in)  :: labels(:)
        type(contraction_t),    intent(out) :: ct
        integer,                intent(out) :: ierr
        character(len=*),       intent(out) :: errmsg

        type(species_basis_t) :: b(BC_MAX_SPECIES)
        character(len=512) :: msg
        character(len=:), allocatable :: msg2
        real(dp), allocatable :: t(:, :)
        integer :: nat, nb, ib, a, l, k, nzl, nout, nsp_red, nsp_full, spinor, i
        integer :: off_f, off_r
        type(radial_rule_t) :: rule

        ierr = 0
        errmsg = ''
        ct%nat = 0
        ct%spinor = 1
        ct%norb_full = 0
        ct%norb_red = 0

        nat = size(labels)
        if (nat < 1) then
            ierr = 1; errmsg = 'build_contraction: structure has no atoms'; return
        end if

        select case (trim(inp%spin_mode))
        case ('none', 'collinear')
            spinor = 1
        case ('noncollinear')
            spinor = 2
        case default
            ierr = 1
            errmsg = "build_contraction: invalid spin_mode '" // trim(inp%spin_mode) // "'"
            return
        end select

        nb = count_leading_str(inp%basis)
        if (nb <= 0) then
            ierr = 1; errmsg = 'build_contraction: basis is required for every species'; return
        end if
        do ib = 1, nb
            call parse_basis(inp%basis(ib), b(ib), ierr, msg2)
            if (ierr /= 0) then
                errmsg = msg2
                return
            end if
        end do

        ct%nat = nat
        ct%spinor = spinor
        allocate(ct%orbs_full(nat), ct%orbs_red(nat), ct%off_full(nat + 1), ct%off_red(nat + 1))
        allocate(ct%atom(nat))
        ct%off_full(1) = 0
        ct%off_red(1) = 0

        do a = 1, nat
            ib = 0
            do k = 1, nb
                if (trim(b(k)%name) == trim(adjustl(labels(a)))) ib = k
            end do
            if (ib == 0) then
                ierr = 1
                errmsg = "atom " // itoa(a) // ": species '" // trim(adjustl(labels(a))) // &
                        "' has no explicit basis entry"
                return
            end if
            nsp_full = b(ib)%nspatial

            ! contagem de funções de saída e validação das regras
            nsp_red = 0
            do l = 0, 3
                nzl = 0
                do k = 1, b(ib)%nshell
                    if (b(ib)%shell(k)%l == l) nzl = nzl + 1
                end do
                if (nzl == 0) cycle
                nout = nzl
                rule%n_out = -1
                do k = 1, inp%nrule
                    if (trim(inp%rule(k)%species) == trim(b(ib)%name) .and. inp%rule(k)%l == l) then
                        rule = inp%rule(k)
                        nout = rule%n_out
                        exit
                    end if
                end do
                if (rule%n_out >= 0) then
                    call check_rule_numerics(rule, b(ib), nzl, a, ierr, msg)
                    if (ierr /= 0) then
                        errmsg = msg
                        return
                    end if
                end if
                nsp_red = nsp_red + nout * (2 * l + 1)
            end do
            if (nsp_red < 1) then
                ierr = 1
                errmsg = 'atom ' // itoa(a) // ": the contraction removed all orbitals of species '" // &
                        trim(b(ib)%name) // "'"
                return
            end if

            ! montagem de T (nsp_full x nsp_red)
            allocate(t(nsp_full, nsp_red))
            t = 0.0_dp
            off_r = 0
            do l = 0, 3
                nzl = 0
                do k = 1, b(ib)%nshell
                    if (b(ib)%shell(k)%l == l) nzl = nzl + 1
                end do
                if (nzl == 0) cycle
                rule%n_out = -1
                do k = 1, inp%nrule
                    if (trim(inp%rule(k)%species) == trim(b(ib)%name) .and. inp%rule(k)%l == l) then
                        rule = inp%rule(k)
                        exit
                    end if
                end do
                if (rule%n_out < 0) then
                    ! shell preservada: identidade
                    do k = 1, nzl
                        do i = 1, 2 * l + 1
                            t(b(ib)%shell(basis_shell_index(b(ib), l, k))%offset + i, off_r + i) = 1.0_dp
                        end do
                        off_r = off_r + 2 * l + 1
                    end do
                else
                    do nout = 1, rule%n_out
                        do k = 1, nzl
                            do i = 1, 2 * l + 1
                                t(b(ib)%shell(basis_shell_index(b(ib), l, k))%offset + i, off_r + i) &
                                    = rule%coeff((nout - 1) * nzl + k)
                            end do
                        end do
                        off_r = off_r + 2 * l + 1
                    end do
                end if
            end do

            ct%orbs_full(a) = spinor * nsp_full
            ct%orbs_red(a) = spinor * nsp_red
            ct%atom(a)%t = t
            deallocate(t)

            ! cópia complexa com o layout do HSIO (não colinear: [up|dn])
            allocate(ct%atom(a)%tc(ct%orbs_full(a), ct%orbs_red(a)))
            ct%atom(a)%tc = (0.0_dp, 0.0_dp)
            ct%atom(a)%tc(1:nsp_full, 1:nsp_red) = cmplx(ct%atom(a)%t, 0.0_dp, kind=dp)
            if (spinor == 2) then
                ct%atom(a)%tc(nsp_full + 1:2 * nsp_full, nsp_red + 1:2 * nsp_red) &
                    = cmplx(ct%atom(a)%t, 0.0_dp, kind=dp)
            end if
        end do

        ct%off_full(1) = 0
        ct%off_red(1) = 0
        do a = 1, nat
            ct%off_full(a + 1) = ct%off_full(a) + ct%orbs_full(a)
            ct%off_red(a + 1) = ct%off_red(a) + ct%orbs_red(a)
        end do
        ct%norb_full = ct%off_full(nat + 1)
        ct%norb_red = ct%off_red(nat + 1)
        if (ct%norb_red < 1) then
            ierr = 1; errmsg = 'build_contraction: the contraction removed all orbitals'; return
        end if
    end subroutine build_contraction

    !> Checagens numéricas de uma regra contra a base: ncoeff = n_out*nzl,
    !> coeficientes finitos e posto completo de C (n_out x nzl) via os
    !> autovalores de C C^T (BC_RANK_RTOL relativo ao maior).
    subroutine check_rule_numerics(rule, b, nzl, iatom, ierr, errmsg)
        type(radial_rule_t),   intent(in)  :: rule
        type(species_basis_t), intent(in)  :: b
        integer,               intent(in)  :: nzl, iatom
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        complex(dp), allocatable :: cc(:, :)
        real(dp), allocatable :: w(:)
        character(len=:), allocatable :: msg
        character(len=64) :: where
        integer :: n, j, z, i
        real(dp) :: dot

        ierr = 0
        errmsg = ''
        where = "atom " // itoa(iatom) // ", species '" // trim(b%name) // "': "
        n = rule%n_out
        if (n < 0 .or. n > nzl) then
            ierr = 1
            errmsg = trim(where) // 'n_out must be <= nzeta'
            return
        end if
        if (n == 0) return
        if (rule%ncoeff /= n * nzl) then
            ierr = 1
            errmsg = trim(where) // 'expected ' // itoa(n * nzl) // ' coefficients, got ' // itoa(rule%ncoeff)
            return
        end if
        if (.not. reals_finite(rule%coeff(1:rule%ncoeff))) then
            ierr = 1
            errmsg = trim(where) // 'coefficients must be finite'
            return
        end if

        allocate(cc(n, n), w(n))
        cc = (0.0_dp, 0.0_dp)
        do j = 1, n
            do i = 1, n
                dot = 0.0_dp
                do z = 1, nzl
                    dot = dot + rule%coeff((j - 1) * nzl + z) * rule%coeff((i - 1) * nzl + z)
                end do
                cc(j, i) = cmplx(dot, 0.0_dp, kind=dp)
            end do
        end do
        call herm_eigvals(cc, w, ierr, msg)
        if (ierr /= 0) then
            errmsg = trim(where) // 'rank check failed'
            return
        end if
        if (w(1) <= max(BC_RANK_RTOL * w(n), tiny(1.0_dp))) then
            ierr = 1
            errmsg = trim(where) // 'coefficient matrix (' // itoa(n) // 'x' // itoa(nzl) // &
                     ', l = ' // itoa(rule%l) // ') is not full rank'
            return
        end if
    end subroutine check_rule_numerics

    ! ==================================================================
    ! APLICAÇÃO
    ! ==================================================================

    !> Contrai um HSDataFormat: Data -> T† Data T em todo R, preservando
    ! Neighbors/IndexList/Degenerecency, meta (E_F/referência), a ordem
    !> atômica e o Title com proveniência; atualiza Orbs/NOrbs. Valida
    ! Orbs/NOrbs do arquivo contra a contração (nspatial, ou 2*nspatial no
    !> não colinear). Nunca para o programa.
    subroutine contract_hs(hs, ct, hs_red, ierr, errmsg)
        type(HSDataFormat),  intent(in)  :: hs
        type(contraction_t), intent(in)  :: ct
        type(HSDataFormat),  intent(out) :: hs_red
        integer,             intent(out) :: ierr
        character(len=*),    intent(out) :: errmsg

        complex(dp), allocatable :: tmp(:, :)
        integer :: a1, a2, ir, n1, n2, r1, r2, o1, o2, p1, p2

        ierr = 0
        errmsg = ''
        if (.not. allocated(hs%Data) .or. .not. allocated(hs%IndexList) .or. .not. allocated(hs%Orbs)) then
            ierr = 1; errmsg = 'contract_hs: Data/Orbs/IndexList not allocated'; return
        end if
        if (hs%nAtoms /= ct%nat .or. size(hs%Orbs) /= ct%nat) then
            ierr = 1
            errmsg = 'contract_hs: nAtoms mismatch (file ' // itoa(hs%nAtoms) // ', contraction ' // itoa(ct%nat) // ')'
            return
        end if
        if (hs%NOrbs /= ct%norb_full) then
            ierr = 1
            errmsg = 'contract_hs: NOrbs mismatch (file ' // itoa(hs%NOrbs) // &
                     ', basis expects ' // itoa(ct%norb_full) // ')'
            return
        end if
        if (any(hs%Orbs /= ct%orbs_full)) then
            do a1 = 1, ct%nat
                if (hs%Orbs(a1) /= ct%orbs_full(a1)) exit
            end do
            ierr = 1
            errmsg = 'contract_hs: Orbs(' // itoa(a1) // ') = ' // itoa(hs%Orbs(a1)) // &
                     ' does not match the declared basis (' // itoa(ct%orbs_full(a1)) // ')'
            return
        end if
        if (any(shape(hs%Data) /= [ct%norb_full, ct%norb_full, hs%Neighbors])) then
            ierr = 1; errmsg = 'contract_hs: Data shape inconsistent with NOrbs/Neighbors'; return
        end if

        hs_red%Title = provenance_title(hs%Title, ct)
        hs_red%Neighbors = hs%Neighbors
        hs_red%nAtoms = hs%nAtoms
        hs_red%NOrbs = ct%norb_red
        hs_red%Orbs = ct%orbs_red
        hs_red%IndexList = hs%IndexList
        if (allocated(hs%Degenerecency)) hs_red%Degenerecency = hs%Degenerecency
        hs_red%has_meta = hs%has_meta
        hs_red%fermi_ev = hs%fermi_ev
        hs_red%fermi_shift_applied = hs%fermi_shift_applied
        hs_red%energy_reference = hs%energy_reference
        allocate(hs_red%Data(ct%norb_red, ct%norb_red, hs%Neighbors))
        hs_red%Data = (0.0_dp, 0.0_dp)

        allocate(tmp(maxval(ct%orbs_full), maxval(ct%orbs_red)))
        do ir = 1, hs%Neighbors
            do a1 = 1, ct%nat
                n1 = ct%orbs_full(a1); o1 = ct%off_full(a1)
                r1 = ct%orbs_red(a1);  p1 = ct%off_red(a1)
                do a2 = 1, ct%nat
                    n2 = ct%orbs_full(a2); o2 = ct%off_full(a2)
                    r2 = ct%orbs_red(a2);  p2 = ct%off_red(a2)
                    if (all(hs%Data(o1 + 1:o1 + n1, o2 + 1:o2 + n2, ir) == (0.0_dp, 0.0_dp))) cycle
                    tmp(1:n1, 1:r2) = matmul(hs%Data(o1 + 1:o1 + n1, o2 + 1:o2 + n2, ir), ct%atom(a2)%tc)
                    hs_red%Data(p1 + 1:p1 + r1, p2 + 1:p2 + r2, ir) &
                        = matmul(conjg(transpose(ct%atom(a1)%tc)), tmp(1:n1, 1:r2))
                end do
            end do
        end do
        deallocate(tmp)
    end subroutine contract_hs

    !> Title de saída: título original + proveniência da contração. O título
    !> original é abreviado se necessário -- a tag 'hscontract' nunca é
    !> perdida por truncamento da atribuição.
    function provenance_title(orig, ct) result(t)
        character(len=*),   intent(in) :: orig
        type(contraction_t), intent(in) :: ct
        character(len=256) :: t, o
        character(len=64) :: tag
        integer :: nkeep
        write(tag, '(A,I0,A,I0,A)') 'hscontract: NOrbs ', ct%norb_full, ' -> ', ct%norb_red, &
            merge(' (noncollinear)', '               ', ct%spinor == 2)
        o = trim(adjustl(orig))
        if (len_trim(o) == 0) then
            t = trim(tag)
            return
        end if
        nkeep = len(t) - 3 - len_trim(tag)   ! espaço para ' | ' + tag
        if (len_trim(o) <= nkeep) then
            t = trim(o) // ' | ' // trim(tag)
        else
            ! abrevia o original: mantém o início e marca o corte com '...'
            t = o(1:max(nkeep - 3, 1)) // '... | ' // trim(tag)
        end if
    end function provenance_title

    ! ==================================================================
    ! CHECAGENS
    ! ==================================================================

    !> Compatibilidade estrutural de dois HSIO (nAtoms, NOrbs, Orbs,
    !> Neighbors, IndexList, Degenerecency). errmsg nomeia os dois arquivos.
    subroutine check_hs_compatible(h1, name1, h2, name2, ierr, errmsg)
        type(HSDataFormat), intent(in) :: h1, h2
        character(len=*),   intent(in) :: name1, name2
        integer,            intent(out) :: ierr
        character(len=*),   intent(out) :: errmsg

        ierr = 0
        errmsg = ''
        if (h1%nAtoms /= h2%nAtoms .or. h1%NOrbs /= h2%NOrbs) then
            ierr = 1
            errmsg = '"' // trim(name1) // '" and "' // trim(name2) // &
                     '": nAtoms/NOrbs differ (HSIO bases must be identical)'
            return
        end if
        if (allocated(h1%Orbs) .and. allocated(h2%Orbs)) then
            if (any(h1%Orbs /= h2%Orbs)) then
                ierr = 1
                errmsg = '"' // trim(name1) // '" and "' // trim(name2) // '": Orbs differ'
                return
            end if
        end if
        if (h1%Neighbors /= h2%Neighbors) then
            ierr = 1
            errmsg = '"' // trim(name1) // '" and "' // trim(name2) // '": Neighbors differ'
            return
        end if
        if (allocated(h1%IndexList) .and. allocated(h2%IndexList)) then
            if (any(h1%IndexList /= h2%IndexList)) then
                ierr = 1
                errmsg = '"' // trim(name1) // '" and "' // trim(name2) // '": IndexList differs'
                return
            end if
        end if
        if (allocated(h1%Degenerecency) .and. allocated(h2%Degenerecency)) then
            if (any(h1%Degenerecency /= h2%Degenerecency)) then
                ierr = 1
                errmsg = '"' // trim(name1) // '" and "' // trim(name2) // '": Degenerecency differs'
                return
            end if
        end if
    end subroutine check_hs_compatible

    !> Estrutura algébrica de H(R)/S(R): H(0) hermitiana e H(-R) = H(R)†
    !> para todo par R/-R presente na lista (pares sem -R são aceitos, p. ex.
    !> arquivos já reduzidos a um semi-espaço). Tolerância relativa ao
    !> maior |elemento| (tol_rel; padrão BC_TOL_PAIR). Falha nomeando R.
    subroutine check_pair_symmetry(hs, tol_rel, ierr, errmsg)
        type(HSDataFormat), intent(in) :: hs
        real(dp),           intent(in) :: tol_rel
        integer,            intent(out) :: ierr
        character(len=*),   intent(out) :: errmsg

        real(dp) :: scale, tol, dev
        integer :: i, j

        ierr = 0
        errmsg = ''
        if (.not. allocated(hs%Data) .or. .not. allocated(hs%IndexList)) then
            ierr = 1; errmsg = 'check_pair_symmetry: Data/IndexList not allocated'; return
        end if
        scale = maxval(abs(hs%Data))
        if (scale == 0.0_dp) return
        tol = tol_rel * scale
        do i = 1, hs%Neighbors
            if (all(hs%IndexList(:, i) == 0)) then
                dev = herm_check(hs%Data(:, :, i))
                if (dev > tol) then
                    ierr = 1
                    write(errmsg, '(A,ES9.2,A,ES9.2,A)') 'R = (0,0,0): |H - H^H| = ', dev, &
                        ' exceeds tolerance ', tol, ' (hermiticity violated)'
                    return
                end if
            else
                do j = 1, hs%Neighbors
                    if (all(hs%IndexList(:, j) == -hs%IndexList(:, i))) then
                        dev = maxval(abs(hs%Data(:, :, j) - transpose(conjg(hs%Data(:, :, i)))))
                        if (dev > tol) then
                            ierr = 1
                            write(errmsg, '(A,3I4,A,ES9.2,A,ES9.2,A)') 'R = (', hs%IndexList(:, i), &
                                '): |H(-R) - H(R)^H| = ', dev, ' exceeds tolerance ', tol, &
                                ' (algebraic structure violated)'
                            return
                        end if
                        exit
                    end if
                end do
            end if
        end do
    end subroutine check_pair_symmetry

    !> Menor autovalor de S_r(k) numa grade k uniforme que inclui k = 0:
    !> k_i = (n_i - 1)/N_i, n_i = 1..check_kgrid(i), em coordenadas reduzidas.
    !> Sem piso artificial: qualquer lambda_min <= tol_s_min é erro, com o k
    !> identificado. Em sucesso devolve o mínimo da grade e o k onde ocorre.
    !>Nunca para o programa.
    subroutine lambda_min_grid(hs, kgrid, tol_s_min, smin, kmin, ierr, errmsg)
        type(HSDataFormat), intent(in)  :: hs
        integer,            intent(in)  :: kgrid(3)
        real(dp),          intent(in)  :: tol_s_min
        real(dp),          intent(out) :: smin
        real(dp),          intent(out) :: kmin(3)
        integer,           intent(out) :: ierr
        character(len=*),  intent(out) :: errmsg

        complex(dp), allocatable :: sk(:, :)
        real(dp), allocatable :: w(:)
        character(len=:), allocatable :: msg2
        real(dp) :: kpt(3)
        integer :: n, i1, i2, i3

        ierr = 0
        errmsg = ''
        smin = 0.0_dp
        kmin = 0.0_dp
        if (.not. allocated(hs%Data) .or. .not. allocated(hs%IndexList)) then
            ierr = 1; errmsg = 'lambda_min_grid: Data/IndexList not allocated'; return
        end if
        if (any(kgrid < 1)) then
            ierr = 1; errmsg = 'lambda_min_grid: kgrid components must be >= 1'; return
        end if
        n = hs%NOrbs
        if (n < 1) then
            ierr = 1; errmsg = 'lambda_min_grid: NOrbs < 1'; return
        end if
        allocate(sk(n, n), w(n))
        smin = huge(1.0_dp)
        do i3 = 1, kgrid(3)
            do i2 = 1, kgrid(2)
                do i1 = 1, kgrid(1)
                    kpt = [real(i1 - 1, dp) / kgrid(1), real(i2 - 1, dp) / kgrid(2), real(i3 - 1, dp) / kgrid(3)]
                    call compute_Hk_at_kpoint(hs%Data, hs%IndexList, sk, kpt)
                    ! S_r(k) precisa ser hermitiana para o zheev; listas de R sem
                    ! o par -R produzem S(k) não hermitiana -- rejeita nomeando k
                    if (maxval(abs(sk)) > 0.0_dp) then
                        if (herm_check(sk) > BC_TOL_PAIR * maxval(abs(sk))) then
                            ierr = 1
                            write(errmsg, '(A,A,A)') 'S_r(k) is not hermitian at k = (', &
                                trim(kfmt(kpt)), ') (check the R/-R list of the input)'
                            return
                        end if
                    end if
                    call herm_eigvals(sk, w, ierr, msg2)
                    if (ierr /= 0) then
                        errmsg = trim(msg2) // ' at k = (' // trim(kfmt(kpt)) // ')'
                        return
                    end if
                    if (w(1) <= tol_s_min) then
                        smin = w(1)
                        kmin = kpt
                        ierr = 1
                        write(errmsg, '(A,ES12.4,A,ES9.2,A,A,A)') 'lambda_min(S_r(k)) = ', w(1), &
                            ' <= tol_s_min = ', tol_s_min, ' at k = (', trim(kfmt(kpt)), ')'
                        return
                    end if
                    if (w(1) < smin) then
                        smin = w(1)
                        kmin = kpt
                    end if
                end do
            end do
        end do
    end subroutine lambda_min_grid

    ! ==================================================================
    ! AUXILIARES
    ! ==================================================================

    pure function kfmt(k) result(s)
        real(dp), intent(in) :: k(3)
        character(len=48) :: s
        write(s, '(3F10.6)') k
        s = adjustl(s)
    end function kfmt

    !> Número de entradas iniciais não vazias de uma lista de strings;
    !> -1 se houver valor não vazio depois de um vazio (buraco na lista).
    pure function count_leading_str(x) result(n)
        character(len=*), intent(in) :: x(:)
        integer :: n, i
        n = 0
        do i = 1, size(x)
            if (len_trim(x(i)) == 0) exit
            n = i
        end do
        if (n < size(x)) then
            do i = n + 2, size(x)
                if (len_trim(x(i)) > 0) then
                    n = -1
                    return
                end if
            end do
        end if
    end function count_leading_str

    !> .true. se todos os reais são finitos (|expoente| < 2^11-1). Teste bit a
    !> bit: comparação com NaN em -ffpe-trap=invalid causaria SIGFPE.
    pure function reals_finite(x) result(ok)
        real(dp), intent(in) :: x(:)
        integer(int64), parameter :: MASK = int(z'7FF0000000000000', int64)
        integer(int64) :: bits
        integer :: i
        logical :: ok
        ok = .true.
        do i = 1, size(x)
            bits = transfer(x(i), bits)
            if (iand(bits, MASK) == MASK) then
                ok = .false.
                return
            end if
        end do
    end function reals_finite

    !> Análogo para o buffer coefficients (sentinela = COEFF_SENTINEL).
    pure function count_leading_coeff(x) result(n)
        real(dp), intent(in) :: x(:)
        integer :: n, i
        n = 0
        do i = 1, size(x)
            if (x(i) == COEFF_SENTINEL) exit
            n = i
        end do
        if (n < size(x)) then
            do i = n + 2, size(x)
                if (x(i) /= COEFF_SENTINEL) then
                    n = -1
                    return
                end if
            end do
        end if
    end function count_leading_coeff

    pure function lchar(l)
        integer, intent(in) :: l
        character(len=1) :: lchar
        lchar = 'spdf'(l + 1:l + 1)
    end function lchar

    pure function itoa(i)
        integer, intent(in) :: i
        character(len=16) :: itoa
        write(itoa, '(I0)') i
    end function itoa

    pure function lower(s) result(u)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: u
        integer :: i, c
        u = s
        do i = 1, len(u)
            c = iachar(u(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) u(i:i) = achar(c + iachar('a') - iachar('A'))
        end do
    end function lower

    pure function upper(s) result(u)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: u
        integer :: i, c
        u = s
        do i = 1, len(u)
            c = iachar(u(i:i))
            if (c >= iachar('a') .and. c <= iachar('z')) u(i:i) = achar(c - iachar('a') + iachar('A'))
        end do
    end function upper

    !> Primeiro token de uma linha que começa com '&'.
    pure function first_token(line) result(tok)
        character(len=*), intent(in) :: line
        character(len=len(line)) :: tok
        integer :: i
        tok = ''
        do i = 1, len_trim(line)
            if (line(i:i) == ' ' .or. line(i:i) == char(9)) exit
            tok(i:i) = line(i:i)
        end do
    end function first_token

end module basis_contract
