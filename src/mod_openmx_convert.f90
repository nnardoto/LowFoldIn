! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! mod_openmx_convert.f90
!
! Monta as matrizes densas H(R) e S(R) (formato HSDataFormat de
! mod_hs_io.f90) a partir dos dados esparsos lidos por mod_scfout_io.f90,
! porte direto da logica de scfout2hsio.py (assemble_dense,
! assemble_noncollinear_H) para Fortran.
!
! Tambem aplica a conversao Hartree -> eV e o deslocamento do nivel de
! Fermi H'(R) = H(R) - E_F*S(R), fisicamente correto no problema
! generalizado nao-ortogonal (ver discussao no scfout2hsio.py).
!
! Nathanael N. Batista / porte para Fortran, Agosto de 2026
! ----------------------------------------------------------------------------------------------------------
module openmx_convert
    use types,     only: dp
    use scfout_io, only: ScfoutData
    implicit none
    private

    public :: HARTREE_TO_EV
    public :: AssembleDense, AssembleNonCollinearH, DoubleBlockDiagonal
    public :: BuildIndexList, HermiticityReport, ApplyFermiShift
    public :: ReorderSpinorAtomMajor

    ! CODATA 2018 (mesmo valor usado internamente pelo OpenMX). O .scfout
    ! guarda Hks/iHks em Hartree; mod_hs_io/mod_core nao fazem nenhuma
    ! conversao de unidades, entao H precisa sair em eV daqui. Overlap (S)
    ! e' adimensional e NUNCA deve ser escalado por este fator.
    real(dp), parameter :: HARTREE_TO_EV = 27.211386245988_dp

contains

    ! -----------------------------------------------------------------
    ! IndexList (3, Neighbors) a partir de atv_ijk -- mesma ordem de
    ! celulas do .scfout (indice 1 = celula central [0,0,0], em
    ! convencao Fortran 1-based)
    ! -----------------------------------------------------------------
    subroutine BuildIndexList(scfout, IndexList)
        type(ScfoutData), intent(in)  :: scfout
        integer, allocatable, intent(out) :: IndexList(:, :)
        integer :: Rn, ncell

        ncell = scfout%TCpyCell + 1
        allocate(IndexList(3, ncell))
        do Rn = 0, scfout%TCpyCell
            IndexList(:, Rn + 1) = scfout%atv_ijk(:, Rn)
        end do
    end subroutine BuildIndexList

    ! -----------------------------------------------------------------
    ! Calcula os offsets globais de orbital (0-based) de cada atomo:
    ! atom_offset(ct_AN) = indice do primeiro orbital do atomo ct_AN
    ! -----------------------------------------------------------------
    subroutine OrbitalOffsets(scfout, atom_offset, NOrbs)
        type(ScfoutData), intent(in)  :: scfout
        integer, allocatable, intent(out) :: atom_offset(:)
        integer, intent(out) :: NOrbs
        integer :: ct_AN

        allocate(atom_offset(scfout%atomnum + 1))
        atom_offset(1) = 0
        do ct_AN = 1, scfout%atomnum
            atom_offset(ct_AN + 1) = atom_offset(ct_AN) + scfout%TotalNumOrbs(ct_AN)
        end do
        NOrbs = atom_offset(scfout%atomnum + 1)
    end subroutine OrbitalOffsets

    ! -----------------------------------------------------------------
    ! Reordena H/S do caso nao-colinear (SpinP_switch=3) da convencao
    ! "spin-major" produzida por AssembleNonCollinearH/DoubleBlockDiagonal
    ! (todos os N orbitais espaciais "up" primeiro, depois todos os N
    ! "down": Hin(off1+i, off2+j) = bloco up-up do par de atomos, e
    ! Hin(N+off1+i, N+off2+j) = bloco down-down, com N=NOrbsSpatial) para a
    ! convencao "atomo-major, spin-menor" (para cada atomo, seus orbitais
    ! up seguidos pelos seus orbitais down, atomo apos atomo) -- que e' a
    ! UNICA convencao compativel com o campo Orbs(:) do formato
    ! HSDataFormat, que assume blocos CONTIGUOS por atomo (ver cabecalho de
    ! mod_hs_io.f90). A rotina WriteHS (em tools/openmx2hsio.f90) ja'
    ! assume essa contiguidade ao gravar Orbs(iAtom) = 2*TotalNumOrbs(iAtom);
    ! sem este reordenamento essa premissa fica FALSA para os dados
    ! realmente gravados, e qualquer ferramenta que localize orbitais de um
    ! atomo via Orbs(:) (como o `unfold`) leria indices errados no caso
    ! nao-colinear -- ainda que o calculo de bandas em si (que nao usa
    ! Orbs(:) para nada) continue correto de qualquer forma.
    !
    ! E' uma transformacao de similaridade por PERMUTACAO (unitaria):
    ! preserva hermiticidade e todos os autovalores, so' reordena rotulos
    ! de base. Aplicada a cada bloco de vizinho independentemente.
    ! -----------------------------------------------------------------
    subroutine ReorderSpinorAtomMajor(scfout, Hin, Sin, Hout, Sout)
        type(ScfoutData), intent(in)  :: scfout
        complex(dp),      intent(in)  :: Hin(:, :, :), Sin(:, :, :)
        complex(dp), allocatable, intent(out) :: Hout(:, :, :), Sout(:, :, :)

        integer, allocatable :: perm(:), atom_offset(:), new_atom_offset(:)
        integer :: NOrbsSpatial, NOrbsTotal, ncell
        integer :: ct_AN, o, TNO, dummyNOrbs

        call OrbitalOffsets(scfout, atom_offset, dummyNOrbs)
        NOrbsSpatial = dummyNOrbs
        NOrbsTotal   = 2 * NOrbsSpatial
        ncell        = size(Hin, 3)

        if (size(Hin, 1) /= NOrbsTotal .or. size(Hin, 2) /= NOrbsTotal .or. &
            size(Sin, 1) /= NOrbsTotal .or. size(Sin, 2) /= NOrbsTotal .or. &
            size(Sin, 3) /= ncell) then
            print '(A)', "[openmx_convert] Erro: dimensao de Hin/Sin incompativel com 2*NOrbsSpatial " // &
                          "em ReorderSpinorAtomMajor."
            stop 1
        end if

        ! Deslocamentos na NOVA ordem (atomo-major): cada atomo ocupa
        ! 2*TotalNumOrbs(ct_AN) posicoes contiguas, up seguido de down.
        allocate(new_atom_offset(scfout%atomnum + 1))
        new_atom_offset(1) = 0
        do ct_AN = 1, scfout%atomnum
            new_atom_offset(ct_AN + 1) = new_atom_offset(ct_AN) + 2 * scfout%TotalNumOrbs(ct_AN)
        end do

        allocate(perm(NOrbsTotal))
        do ct_AN = 1, scfout%atomnum
            TNO = scfout%TotalNumOrbs(ct_AN)
            do o = 1, TNO
                ! perm(nova_posicao) = posicao_antiga
                perm(new_atom_offset(ct_AN) + o)       = atom_offset(ct_AN) + o                  ! up
                perm(new_atom_offset(ct_AN) + TNO + o) = NOrbsSpatial + atom_offset(ct_AN) + o    ! down
            end do
        end do

        allocate(Hout(NOrbsTotal, NOrbsTotal, ncell))
        allocate(Sout(NOrbsTotal, NOrbsTotal, ncell))
        Hout = Hin(perm, perm, :)
        Sout = Sin(perm, perm, :)
    end subroutine ReorderSpinorAtomMajor

    ! -----------------------------------------------------------------
    ! Monta um bloco denso complex(dp) (NOrbs,NOrbs,Neighbors) a partir
    ! de um canal "flat" (Hks_flat(:,ispin) ou OLP_flat), posicionando
    ! cada bloco atomo-vizinho na celula certa (via ncn_flat) e nos
    ! indices de orbital globais certos (via atom_offset). Port direto
    ! de assemble_dense() do scfout2hsio.py.
    ! -----------------------------------------------------------------
    subroutine AssembleDense(scfout, flat_channel, Hout, NOrbsOut)
        type(ScfoutData), intent(in)  :: scfout
        real(dp),          intent(in)  :: flat_channel(:)
        complex(dp), allocatable, intent(out) :: Hout(:, :, :)
        integer, intent(out) :: NOrbsOut

        integer, allocatable :: atom_offset(:)
        integer :: ct_AN, h_AN, p, off1, off2, TNO1, TNO2, Gh_AN, Rn
        integer :: i, j, blk, ncell

        call OrbitalOffsets(scfout, atom_offset, NOrbsOut)
        ncell = scfout%TCpyCell + 1

        allocate(Hout(NOrbsOut, NOrbsOut, ncell))
        Hout = cmplx(0.0_dp, 0.0_dp, kind=dp)

        p = 0
        do ct_AN = 1, scfout%atomnum
            TNO1 = scfout%TotalNumOrbs(ct_AN)
            off1 = atom_offset(ct_AN)
            do h_AN = 0, scfout%FNAN(ct_AN)
                p = p + 1
                Gh_AN = scfout%natn_flat(scfout%neighbor_offset(ct_AN) + h_AN + 1)
                Rn    = scfout%ncn_flat(scfout%neighbor_offset(ct_AN) + h_AN + 1)
                TNO2  = scfout%TotalNumOrbs(Gh_AN)
                off2  = atom_offset(Gh_AN)
                blk   = scfout%block_offset(p)

                do i = 1, TNO1
                    do j = 1, TNO2
                        Hout(off1 + i, off2 + j, Rn + 1) = Hout(off1 + i, off2 + j, Rn + 1) &
                            + cmplx(flat_channel(blk + (i - 1) * TNO2 + j), 0.0_dp, kind=dp)
                    end do
                end do
            end do
        end do
    end subroutine AssembleDense

    ! -----------------------------------------------------------------
    ! Retorna o bloco (TNO1,TNO2) do par (ct_AN,h_AN) de um canal flat,
    ! ja' reorganizado como matriz Fortran normal (i,j).
    ! -----------------------------------------------------------------
    subroutine ExtractBlock(scfout, flat_channel, ct_AN, h_AN_idx_global, TNO1, TNO2, Bout)
        type(ScfoutData), intent(in) :: scfout
        real(dp),          intent(in) :: flat_channel(:)
        integer,           intent(in) :: ct_AN            ! nao usado diretamente, mantido p/ clareza
        integer,           intent(in) :: h_AN_idx_global   ! indice global do par p (1..n_pairs)
        integer,           intent(in) :: TNO1, TNO2
        real(dp), allocatable, intent(out) :: Bout(:, :)
        integer :: i, j, blk

        blk = scfout%block_offset(h_AN_idx_global)
        allocate(Bout(TNO1, TNO2))
        do i = 1, TNO1
            do j = 1, TNO2
                Bout(i, j) = flat_channel(blk + (i - 1) * TNO2 + j)
            end do
        end do
    end subroutine ExtractBlock

    ! -----------------------------------------------------------------
    ! Encontra, para cada celula Rn, a celula Rn_neg tal que
    ! atv_ijk(:,Rn_neg) == -atv_ijk(:,Rn). Necessario para montar o
    ! bloco spin down-up (H_du) por hermiticidade no caso nao-colinear.
    ! -----------------------------------------------------------------
    subroutine NegCellMap(scfout, neg_cell)
        type(ScfoutData), intent(in) :: scfout
        integer, allocatable, intent(out) :: neg_cell(:)   ! 0:TCpyCell
        integer :: Rn, Rn2

        allocate(neg_cell(0:scfout%TCpyCell))
        neg_cell = -1
        do Rn = 0, scfout%TCpyCell
            do Rn2 = 0, scfout%TCpyCell
                if (all(scfout%atv_ijk(:, Rn2) == -scfout%atv_ijk(:, Rn))) then
                    neg_cell(Rn) = Rn2
                    exit
                end if
            end do
            if (neg_cell(Rn) < 0) then
                print '(A,I0,A)', "[openmx_convert] Erro: nao foi possivel achar a " // &
                    "celula oposta a atv_ijk(:,", Rn, "). A montagem do caso nao-colinear " // &
                    "depende dessa simetria (R e -R ambos presentes no .scfout)."
                stop 1
            end if
        end do
    end subroutine NegCellMap

    ! -----------------------------------------------------------------
    ! Para um atomo Gh_AN, encontra o indice de vizinho (global, 1..n_pairs)
    ! tal que (natn(Gh_AN,h'), ncn(Gh_AN,h')) == (target_atom, target_cell).
    ! Busca linear dentro da lista de vizinhos de Gh_AN (tipicamente
    ! pequena, dezenas de elementos).
    ! -----------------------------------------------------------------
    function FindMirrorNeighbor(scfout, Gh_AN, target_atom, target_cell) result(p_mirror)
        type(ScfoutData), intent(in) :: scfout
        integer,           intent(in) :: Gh_AN, target_atom, target_cell
        integer :: p_mirror
        integer :: h, off, natn_val, ncn_val

        p_mirror = -1
        off = scfout%neighbor_offset(Gh_AN)
        do h = 0, scfout%FNAN(Gh_AN)
            natn_val = scfout%natn_flat(off + h + 1)
            ncn_val  = scfout%ncn_flat(off + h + 1)
            if (natn_val == target_atom .and. ncn_val == target_cell) then
                ! indice global do par (Gh_AN, h): soma dos pares de todos
                ! os atomos anteriores + h + 1
                p_mirror = GlobalPairIndex(scfout, Gh_AN, h)
                return
            end if
        end do
    end function FindMirrorNeighbor

    ! Indice global (1-based, na mesma ordem de leitura ct_AN/h_AN) do
    ! par (ct_AN, h_AN [0-based])
    function GlobalPairIndex(scfout, ct_AN, h_AN) result(p)
        type(ScfoutData), intent(in) :: scfout
        integer,           intent(in) :: ct_AN, h_AN
        integer :: p, k

        ! numero de pares de todos os atomos 1..ct_AN-1
        p = 0
        do k = 1, ct_AN - 1
            p = p + scfout%FNAN(k) + 1
        end do
        p = p + h_AN + 1
    end function GlobalPairIndex

    ! -----------------------------------------------------------------
    ! Monta o Hamiltoniano spinor completo (caso SpinP_switch==3), shape
    ! (2*NOrbs,2*NOrbs,Neighbors), convencao de bloco:
    !
    !   Data(:,:,n) = [ H_uu  H_ud ]
    !                 [ H_du  H_dd ]
    !
    ! H_uu, H_dd, H_ud vem direto de Hks+i*iHks (Hks[2],Hks[3],iHks[2]
    ! para H_ud). H_du e' obtido por hermiticidade via bloco espelhado:
    ! H_du(R)_{ct_AN,Gh_AN} = conj( H_ud(-R)_{Gh_AN,ct_AN} )^T.
    ! Port direto de assemble_noncollinear_H() do scfout2hsio.py.
    ! -----------------------------------------------------------------
    subroutine AssembleNonCollinearH(scfout, Hout, NOrbsOut2)
        type(ScfoutData), intent(in)  :: scfout
        complex(dp), allocatable, intent(out) :: Hout(:, :, :)
        integer, intent(out) :: NOrbsOut2

        integer, allocatable :: atom_offset(:), neg_cell(:)
        integer :: N, ncell
        integer :: ct_AN, h_AN, p, off1, off2, TNO1, TNO2, Gh_AN, Rn, Rn_neg
        integer :: i, j, p_mirror
        real(dp), allocatable :: Huu_r(:,:), Hdd_r(:,:), Hud_re(:,:), Hud_im1(:,:), Hud_im2(:,:)
        real(dp), allocatable :: Hud_mirror_re(:,:), Hud_mirror_im1(:,:), Hud_mirror_im2(:,:)
        complex(dp), allocatable :: Hud_mirror(:,:), Hdu(:,:)
        integer :: NOrbs

        call OrbitalOffsets(scfout, atom_offset, NOrbs)
        NOrbsOut2 = 2 * NOrbs
        ncell = scfout%TCpyCell + 1

        allocate(Hout(NOrbsOut2, NOrbsOut2, ncell))
        Hout = cmplx(0.0_dp, 0.0_dp, kind=dp)

        call NegCellMap(scfout, neg_cell)

        N = NOrbs
        p = 0
        do ct_AN = 1, scfout%atomnum
            TNO1 = scfout%TotalNumOrbs(ct_AN)
            off1 = atom_offset(ct_AN)
            do h_AN = 0, scfout%FNAN(ct_AN)
                p = p + 1
                Gh_AN = scfout%natn_flat(scfout%neighbor_offset(ct_AN) + h_AN + 1)
                Rn    = scfout%ncn_flat(scfout%neighbor_offset(ct_AN) + h_AN + 1)
                TNO2  = scfout%TotalNumOrbs(Gh_AN)
                off2  = atom_offset(Gh_AN)

                call ExtractBlock(scfout, scfout%Hks_flat(:, 0), ct_AN, p, TNO1, TNO2, Huu_r)
                call ExtractBlock(scfout, scfout%Hks_flat(:, 1), ct_AN, p, TNO1, TNO2, Hdd_r)
                call ExtractBlock(scfout, scfout%Hks_flat(:, 2), ct_AN, p, TNO1, TNO2, Hud_re)
                call ExtractBlock(scfout, scfout%Hks_flat(:, 3), ct_AN, p, TNO1, TNO2, Hud_im1)
                call ExtractBlock(scfout, scfout%iHks_flat(:, 2), ct_AN, p, TNO1, TNO2, Hud_im2)

                ! Huu = Hks[0] + i*iHks[0]; Hdd = Hks[1] + i*iHks[1]
                block
                    real(dp), allocatable :: iHuu(:,:), iHdd(:,:)
                    call ExtractBlock(scfout, scfout%iHks_flat(:, 1), ct_AN, p, TNO1, TNO2, iHdd)
                    call ExtractBlock(scfout, scfout%iHks_flat(:, 0), ct_AN, p, TNO1, TNO2, iHuu)
                    do i = 1, TNO1
                        do j = 1, TNO2
                            Hout(off1 + i, off2 + j, Rn + 1) = Hout(off1 + i, off2 + j, Rn + 1) &
                                + cmplx(Huu_r(i, j), iHuu(i, j), kind=dp)
                            Hout(N + off1 + i, N + off2 + j, Rn + 1) = Hout(N + off1 + i, N + off2 + j, Rn + 1) &
                                + cmplx(Hdd_r(i, j), iHdd(i, j), kind=dp)
                            ! H_ud = Hks[2] + i*(Hks[3] + iHks[2])
                            Hout(off1 + i, N + off2 + j, Rn + 1) = Hout(off1 + i, N + off2 + j, Rn + 1) &
                                + cmplx(Hud_re(i, j), Hud_im1(i, j) + Hud_im2(i, j), kind=dp)
                        end do
                    end do
                end block

                ! H_du via bloco espelhado (Gh_AN -> ct_AN na celula -Rn)
                Rn_neg = neg_cell(Rn)
                p_mirror = FindMirrorNeighbor(scfout, Gh_AN, ct_AN, Rn_neg)
                if (p_mirror < 0) then
                    print '(A,I0,A,I0,A)', "[openmx_convert] Erro: bloco espelhado (atomo ", &
                        Gh_AN, " -> ", ct_AN, ") nao encontrado na lista de vizinhos do .scfout. " // &
                        "Nao e' possivel montar H_du."
                    stop 1
                end if

                call ExtractBlock(scfout, scfout%Hks_flat(:, 2), Gh_AN, p_mirror, TNO2, TNO1, Hud_mirror_re)
                call ExtractBlock(scfout, scfout%Hks_flat(:, 3), Gh_AN, p_mirror, TNO2, TNO1, Hud_mirror_im1)
                call ExtractBlock(scfout, scfout%iHks_flat(:, 2), Gh_AN, p_mirror, TNO2, TNO1, Hud_mirror_im2)

                allocate(Hud_mirror(TNO2, TNO1))
                Hud_mirror = cmplx(Hud_mirror_re, Hud_mirror_im1 + Hud_mirror_im2, kind=dp)

                allocate(Hdu(TNO1, TNO2))
                Hdu = transpose(conjg(Hud_mirror))

                do i = 1, TNO1
                    do j = 1, TNO2
                        Hout(N + off1 + i, off2 + j, Rn + 1) = Hout(N + off1 + i, off2 + j, Rn + 1) + Hdu(i, j)
                    end do
                end do

                deallocate(Hud_mirror, Hdu)
                deallocate(Huu_r, Hdd_r, Hud_re, Hud_im1, Hud_im2)
            end do
        end do
    end subroutine AssembleNonCollinearH

    ! -----------------------------------------------------------------
    ! Dado um bloco espacial (N,N,Neighbors), devolve (2N,2N,Neighbors)
    ! com o mesmo bloco copiado na diagonal (up-up e down-down) e zero
    ! nos quadrantes fora da diagonal. Usado para dobrar o Overlap S no
    ! caso nao-colinear.
    ! -----------------------------------------------------------------
    subroutine DoubleBlockDiagonal(Sin, Sout)
        complex(dp), intent(in)  :: Sin(:, :, :)
        complex(dp), allocatable, intent(out) :: Sout(:, :, :)
        integer :: N, ncell

        N = size(Sin, 1)
        ncell = size(Sin, 3)
        allocate(Sout(2 * N, 2 * N, ncell))
        Sout = cmplx(0.0_dp, 0.0_dp, kind=dp)
        Sout(1:N, 1:N, :) = Sin
        Sout(N + 1:2 * N, N + 1:2 * N, :) = Sin
    end subroutine DoubleBlockDiagonal

    ! -----------------------------------------------------------------
    ! Diagnostico: soma Data sobre todas as celulas (== H(k=Gamma)) e
    ! mede o quanto isso desvia de ser hermitiano. Deve ser ~0 (ordem de
    ! erro de arredondamento) se a montagem estiver fisicamente correta.
    ! -----------------------------------------------------------------
    subroutine HermiticityReport(data, label)
        complex(dp), intent(in) :: data(:, :, :)
        character(len=*), intent(in) :: label
        complex(dp), allocatable :: gamma(:, :)
        real(dp) :: dev, scale

        allocate(gamma(size(data,1), size(data,2)))
        gamma = sum(data, dim=3)
        dev = maxval(abs(gamma - transpose(conjg(gamma))))
        scale = maxval(abs(gamma)) + 1.0e-30_dp
        print '(A, A, ES12.3E3, A, ES12.3E3, A)', &
            "  [diagnostico ", trim(label), dev, " ] desvio de hermiticidade em Gamma (relativo: ", &
            dev/scale, ") -- deve ser proximo de zero"
    end subroutine HermiticityReport

    ! -----------------------------------------------------------------
    ! Desloca H para que o nivel de Fermi caia em 0 eV: H'(R) = H(R) -
    ! E_F*S(R) para TODO vizinho R (forma correta no problema
    ! generalizado nao-ortogonal -- nao basta subtrair so' do termo
    ! on-site, ja que S tem overlaps entre vizinhos).
    ! -----------------------------------------------------------------
    subroutine ApplyFermiShift(Hdata, Sdata, fermi_ev)
        complex(dp), intent(inout) :: Hdata(:, :, :)
        complex(dp), intent(in)    :: Sdata(:, :, :)
        real(dp),    intent(in)    :: fermi_ev

        Hdata = Hdata - fermi_ev * Sdata
    end subroutine ApplyFermiShift

end module openmx_convert