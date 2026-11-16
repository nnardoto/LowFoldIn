! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! examples/surface_spectral/gen_model.f90
!
! MODELO SINTÉTICO -- NÃO SÃO DADOS DE DFT, NEM DE ARPES.
!
! Escreve, no diretório ATUAL, uma slab tight-binding determinística, do tipo
! LCAO e não-ortogonal:
!   model_H.hsio, model_S.hsio   (HSIO v2, via HSDataFormat%WriteFile)
!   model.nml                    (arquivo de estrutura, contracts §1, com &SPECIES)
!
! Geometria: rede quadrada com a = 3.0 A no plano (a1 = (3,0,0), a2 = (0,3,0)),
! a3 = (0,0,40) (vácuo). NL = 8 camadas ao longo de z, um átomo por camada,
! espaçamento de 2.0 A, z_l = 13 + 2(l-1) A. Átomo l = camada l (átomo 1 = base,
! átomo 8 = topo). As espécies se alternam: camadas ímpares 'A' (2 orbitais: s, pz)
! e camadas pares 'B' (1 orbital: s). Camada da base = A (átomo 1), camada do
! topo = B (átomo 8).
!
! Hamiltoniano (eV); Data(:,:,iR) = <i na célula 0| H |j na célula IndexList(:,iR)>,
! H(k) = sum_R Data(:,:,R) exp(+i 2 pi k.R). Células: R = 0, +-a1, +-a2.
!   onsite       : A_s = EPS_AS, A_pz = EPS_AP, B_s = EPS_BS
!   superfície   : B_s do topo deslocado por DTOP, A_s da base deslocado por DBOT
!   intercamadas : (R = 0, somente camadas adjacentes)  s-s = TZ_SS,
!                  A_pz - B_s = +TZ_PS se B está acima de A, -TZ_PS se abaixo (pz ímpar)
!   NN no plano  : (R = +-a1, +-a2, mesmo átomo) s-s = TP_S, pz-pz = TP_P
! Overlap: identidade no sítio; mesma conectividade com SZ_SS, SZ_PS (mesma regra
! de sinal), SP_S, SP_P. Todas as matrizes são reais simétricas no sentido
! Data(:,:,-R) = Data(:,:,R)^T, de modo que H(k) e S(k) são hermitianas.
! E_F = 0 por construção (nenhum deslocamento é aplicado em lugar algum).
! ----------------------------------------------------------------------------------------------------------
program gen_model
    use types, only: dp
    use HS_io, only: HSDataFormat
    implicit none

    integer,  parameter :: NL = 8
    real(dp), parameter :: A_LAT = 3.0_dp, DZ = 2.0_dp, C_LAT = 40.0_dp, Z0 = 13.0_dp

    ! ---- parâmetros do modelo (eV) ----
    real(dp), parameter :: EPS_AS = -1.20_dp, EPS_AP = 4.00_dp, EPS_BS = 1.80_dp
    real(dp), parameter :: DTOP   = -2.00_dp, DBOT = 1.80_dp
    real(dp), parameter :: TZ_SS  = -0.60_dp, TZ_PS = 0.50_dp
    real(dp), parameter :: TP_S   = -0.15_dp, TP_P  = 0.10_dp
    ! ---- overlaps (adimensionais) ----
    real(dp), parameter :: SZ_SS  = 0.10_dp, SZ_PS = 0.05_dp
    real(dp), parameter :: SP_S   = 0.08_dp, SP_P  = 0.05_dp

    integer, parameter :: NR = 5
    integer, parameter :: RLIST(3, NR) = reshape([0,0,0, 1,0,0, -1,0,0, 0,1,0, 0,-1,0], [3, NR])

    type(HSDataFormat) :: H, S
    integer            :: norb(NL), off(NL), ntot, l, iR, u, is_a, is_b, ip_a
    real(dp)           :: sgn

    do l = 1, NL
        norb(l) = merge(2, 1, mod(l, 2) == 1)
    end do
    off(1) = 0
    do l = 2, NL
        off(l) = off(l-1) + norb(l-1)
    end do
    ntot = off(NL) + norb(NL)

    call init(H, 'Synthetic model slab H (NOT DFT) - examples/surface_spectral')
    call init(S, 'Synthetic model slab S (NOT DFT) - examples/surface_spectral')

    ! ---- onsite (R = 0) ----
    do l = 1, NL
        if (mod(l, 2) == 1) then
            call setv(H, off(l)+1, off(l)+1, 1, EPS_AS + merge(DBOT, 0.0_dp, l == 1))
            call setv(H, off(l)+2, off(l)+2, 1, EPS_AP)
            call setv(S, off(l)+1, off(l)+1, 1, 1.0_dp)
            call setv(S, off(l)+2, off(l)+2, 1, 1.0_dp)
        else
            call setv(H, off(l)+1, off(l)+1, 1, EPS_BS + merge(DTOP, 0.0_dp, l == NL))
            call setv(S, off(l)+1, off(l)+1, 1, 1.0_dp)
        end if
    end do

    ! ---- intercamadas (R = 0), camadas l e l+1 ----
    do l = 1, NL - 1
        if (mod(l, 2) == 1) then          ! A em l (abaixo), B em l+1 (acima)
            is_a = off(l) + 1; ip_a = off(l) + 2; is_b = off(l+1) + 1; sgn = 1.0_dp
        else                              ! B em l (abaixo), A em l+1 (acima)
            is_b = off(l) + 1; is_a = off(l+1) + 1; ip_a = off(l+1) + 2; sgn = -1.0_dp
        end if
        call setsym(H, is_a, is_b, TZ_SS)
        call setsym(S, is_a, is_b, SZ_SS)
        call setsym(H, ip_a, is_b, sgn * TZ_PS)
        call setsym(S, ip_a, is_b, sgn * SZ_PS)
    end do

    ! ---- primeiros vizinhos no plano (mesmo átomo, R = +-a1, +-a2) ----
    do iR = 2, NR
        do l = 1, NL
            call setv(H, off(l)+1, off(l)+1, iR, TP_S)
            call setv(S, off(l)+1, off(l)+1, iR, SP_S)
            if (norb(l) == 2) then
                call setv(H, off(l)+2, off(l)+2, iR, TP_P)
                call setv(S, off(l)+2, off(l)+2, iR, SP_P)
            end if
        end do
    end do

    call H%WriteFile('model_H.hsio')
    call S%WriteFile('model_S.hsio')

    ! ---- arquivo de estrutura ----
    open(newunit=u, file='model.nml', status='replace', action='write')
    write(u, '(A)') '! SYNTHETIC MODEL slab (NOT DFT) generated by examples/surface_spectral/gen_model.f90'
    write(u, '(A)') '! 8 layers along z (atom 1 = bottom, atom 8 = top), species A (s,pz) / B (s)'
    write(u, '(A)') '&STRUCTURE'
    write(u, '(A,3(F12.6,:,","))') '    LatticeVectors(1,:) = ', A_LAT, 0.0_dp, 0.0_dp
    write(u, '(A,3(F12.6,:,","))') '    LatticeVectors(2,:) = ', 0.0_dp, A_LAT, 0.0_dp
    write(u, '(A,3(F12.6,:,","))') '    LatticeVectors(3,:) = ', 0.0_dp, 0.0_dp, C_LAT
    write(u, '(A)') '/'
    write(u, '(A)') ''
    write(u, '(A)') '&SPECIES'
    write(u, '(A)') '2'
    write(u, '(A)') 'A   2'
    write(u, '(A)') 'B   1'
    write(u, '(A)') '/'
    write(u, '(A)') ''
    write(u, '(A)') '&ATOMS'
    write(u, '(I0)') NL
    do l = 1, NL
        write(u, '(A,3F12.6)') merge('A ', 'B ', mod(l, 2) == 1), 0.0_dp, 0.0_dp, Z0 + DZ * real(l - 1, dp)
    end do
    write(u, '(A)') '/'
    close(u)

    print '(A,I0,A,I0,A,I0)', '[gen_model] synthetic model written: nAtoms = ', NL, &
        '  NOrbs = ', ntot, '  cells = ', NR
    print '(A)', '[gen_model] files: model_H.hsio model_S.hsio model.nml (SYNTHETIC, NOT DFT)'

contains

    subroutine init(M, title)
        type(HSDataFormat), intent(inout) :: M
        character(len=*),   intent(in)    :: title
        M%Title     = title
        M%Neighbors = NR
        M%nAtoms    = NL
        M%NOrbs     = ntot
        allocate(M%Degenerecency(NR), source=1)
        allocate(M%Orbs(NL), source=norb)
        allocate(M%IndexList(3, NR), source=RLIST)
        allocate(M%Data(ntot, ntot, NR), source=(0.0_dp, 0.0_dp))
    end subroutine init

    subroutine setv(M, i, j, iRn, v)
        type(HSDataFormat), intent(inout) :: M
        integer,            intent(in)    :: i, j, iRn
        real(dp),           intent(in)    :: v
        M%Data(i, j, iRn) = cmplx(v, 0.0_dp, kind=dp)
    end subroutine setv

    !> Par simétrico em R = 0: (i,j) e (j,i).
    subroutine setsym(M, i, j, v)
        type(HSDataFormat), intent(inout) :: M
        integer,            intent(in)    :: i, j
        real(dp),           intent(in)    :: v
        call setv(M, i, j, 1, v)
        call setv(M, j, i, 1, v)
    end subroutine setsym

end program gen_model
