! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Modelos de contato entre dois materiais empilhados (A abaixo, B acima).
!>
!> O bloco de contato é H_{A,B} (linhas = orbitais da PL de A, colunas = PL de
!> B), na orientação "de baixo para cima", e o mesmo para S:
!>   - 'none':      H_AB = S_AB = 0 (materiais desacoplados);
!>   - 'scaled':    H_AB = lambda_H H_ref, S_AB = lambda_S S_ref, com H_ref o
!>                  bloco entre PLs vizinhas (H_{m,m+1}) de um material de
!>                  referência ('substrate' ou 'film'); exige que as dimensões
!>                  das PLs sejam compatíveis;
!>   - 'sk':        acoplamento de Slater–Koster posto à mão entre os átomos
!>                  vizinhos dos dois lados (módulo sk_contact), S_AB = 0;
!>                  montado em green_engine a partir da geometria;
!>   - 'extracted': devolve erro aqui; a fase 4 é o modo &GREEN mode =
!>                  'extracted' (módulo extracted_green), em que a região
!>                  central inteira e a costura vêm de uma supercélula DFT.
!>
!> Potenciais: H_mn -> H_mn + (V_m + V_n)/2 S_mn (forma hermitiana local). Um
!> V constante reproduz o deslocamento rígido H + Delta S. A aplicação é feita
!> em central_green::assemble_central/attach_lead; aqui fica apply_potential
!> para uso direto em blocos isolados.
module interface_model
    use types, only: dp
    implicit none
    private

    public :: contact_t, contact_validate, build_contact, apply_potential

    integer, parameter, public :: CONTACT_NONE = 0, CONTACT_SCALED = 1, CONTACT_EXTRACTED = 2, CONTACT_SK = 3
    integer, parameter, public :: REF_SUBSTRATE = 1, REF_FILM = 2

    type :: contact_t
        integer  :: model = CONTACT_NONE
        integer  :: reference = REF_SUBSTRATE
        real(dp) :: lambda_h = 1.0_dp
        real(dp) :: lambda_s = 1.0_dp
    end type contact_t

contains

    !> Converte os nomes do namelist em contact_t.
    subroutine contact_validate(model, reference, lambda_h, lambda_s, c, ierr, errmsg)
        character(len=*), intent(in)  :: model, reference
        real(dp),         intent(in)  :: lambda_h, lambda_s
        type(contact_t),  intent(out) :: c
        integer,          intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        ierr = 0
        errmsg = ''
        select case (trim(adjustl(model)))
        case ('none')
            c%model = CONTACT_NONE
        case ('scaled')
            c%model = CONTACT_SCALED
        case ('sk')
            c%model = CONTACT_SK
        case ('extracted')
            c%model = CONTACT_EXTRACTED
            ierr = 2
            errmsg = "interface model 'extracted' (phase 4) is not a contact block: use &GREEN " // &
                     "mode = 'extracted' with &EXTRACTED and &LEAD (central region taken from a supercell, " // &
                     "docs/surfgreen.md); here use 'none', 'scaled' or 'sk'"
            return
        case default
            ierr = 1; errmsg = "&INTERFACE model must be 'none', 'scaled', 'sk' or 'extracted'"; return
        end select
        select case (trim(adjustl(reference)))
        case ('substrate')
            c%reference = REF_SUBSTRATE
        case ('film')
            c%reference = REF_FILM
        case default
            ierr = 1; errmsg = "&INTERFACE reference must be 'substrate' or 'film'"; return
        end select
        c%lambda_h = lambda_h
        c%lambda_s = lambda_s
    end subroutine contact_validate

    !> Bloco de contato H_AB, S_AB (nrow x ncol). href, sref: H_{m,m+1} e
    !> S_{m,m+1} do material de referência no mesmo k_par (ignorados em 'none').
    subroutine build_contact(c, nrow, ncol, href, sref, hc, sc, ierr, errmsg)
        type(contact_t), intent(in)  :: c
        integer,         intent(in)  :: nrow, ncol
        complex(dp),     intent(in)  :: href(:, :), sref(:, :)
        complex(dp), allocatable, intent(out) :: hc(:, :), sc(:, :)
        integer,         intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        ierr = 0
        errmsg = ''
        allocate(hc(nrow, ncol), sc(nrow, ncol))
        select case (c%model)
        case (CONTACT_NONE)
            hc = (0.0_dp, 0.0_dp)
            sc = (0.0_dp, 0.0_dp)
        case (CONTACT_SCALED)
            if (size(href, 1) /= nrow .or. size(href, 2) /= ncol) then
                ierr = 1
                errmsg = "interface model 'scaled': the reference coupling block does not match the " // &
                         'principal-layer sizes of the two materials'
                return
            end if
            hc = c%lambda_h * href
            sc = c%lambda_s * sref
        case default
            ierr = 2; errmsg = 'build_contact: unsupported interface model'
        end select
    end subroutine build_contact

    !> h <- h + (vm + vn)/2 s.
    pure subroutine apply_potential(h, s, vm, vn)
        complex(dp), intent(inout) :: h(:, :)
        complex(dp), intent(in)    :: s(:, :)
        real(dp),    intent(in)    :: vm, vn
        h = h + (0.5_dp * (vm + vn)) * s
    end subroutine apply_potential

end module interface_model
