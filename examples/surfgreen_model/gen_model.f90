! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! examples/surfgreen_model/gen_model.f90
!
! MODELOS SINTÉTICOS -- NÃO SÃO DADOS DE DFT.
!
! Escreve, no diretório ATUAL, dois bulks tight-binding não ortogonais numa
! célula adaptada à superfície (a1, a2 no plano, a3 inclinado ligando as camadas):
!   sub_H.hsio, sub_S.hsio, sub_struct.nml     substrato A (espécie 'A', orbitais s e pz)
!   film_H.hsio, film_S.hsio, film_struct.nml  filme B     (espécie 'B', orbitais s e pz)
! As redes no plano são iguais (filme "deformado" para o substrato).
!
! H(R) (eV), R = 0, ±a1, ±a2, ±a3, ordem dos orbitais (s, pz):
!   sítio : diag(EPS_S, EPS_P)
!   ±a1,±a2: s-s = TP_SS, pz-pz = TP_PP
!   +a3   : [[TZ_SS, TZ_SP], [-TZ_SP, TZ_PP]]; -a3 = transposta (pz ímpar)
! S(R): identidade no sítio, mesma conectividade com amplitudes menores.
! O dataset meta do .hsio registra E_F = 0 e energy_reference = 'absolute'.
! ----------------------------------------------------------------------------------------------------------
program gen_surfgreen_model
    use types, only: dp
    use HS_io, only: HSDataFormat
    implicit none

    call write_material('sub', 'A', [0.8_dp, 0.5_dp, 3.2_dp], &
                        eps=[-1.0_dp, 1.0_dp], tp=[-0.30_dp, 0.10_dp], tz=[-0.50_dp, 0.40_dp, 0.60_dp], &
                        sscale=0.08_dp)
    call write_material('film', 'B', [0.6_dp, 0.4_dp, 3.0_dp], &
                        eps=[-0.4_dp, 1.6_dp], tp=[-0.24_dp, 0.08_dp], tz=[-0.40_dp, 0.32_dp, 0.48_dp], &
                        sscale=0.06_dp)
    print '(A)', 'Wrote sub_H.hsio, sub_S.hsio, sub_struct.nml, film_H.hsio, film_S.hsio, film_struct.nml (synthetic models)'

contains

    subroutine write_material(prefix, species, a3, eps, tp, tz, sscale)
        character(len=*), intent(in) :: prefix, species
        real(dp), intent(in) :: a3(3), eps(2), tp(2), tz(3), sscale
        type(HSDataFormat) :: h, s
        integer, parameter :: nr = 7
        integer :: rl(3, nr), ir, u
        rl = reshape([0,0,0, 1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1], [3, nr])
        h%Title = 'synthetic ' // prefix // ' H (eV)'
        h%Neighbors = nr; h%nAtoms = 1; h%NOrbs = 2
        allocate(h%Orbs(1), h%IndexList(3, nr), h%Degenerecency(nr), h%Data(2, 2, nr))
        h%Orbs = 2; h%IndexList = rl; h%Degenerecency = 1
        h%Data = (0.0_dp, 0.0_dp)
        s = h
        s%Title = 'synthetic ' // prefix // ' S'
        do ir = 1, nr
            select case (ir)
            case (1)
                h%Data(:, :, ir) = cmplx(reshape([eps(1), 0.0_dp, 0.0_dp, eps(2)], [2, 2]), 0.0_dp, dp)
                s%Data(:, :, ir) = cmplx(reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), 0.0_dp, dp)
            case (2:5)
                h%Data(:, :, ir) = cmplx(reshape([tp(1), 0.0_dp, 0.0_dp, tp(2)], [2, 2]), 0.0_dp, dp)
                s%Data(:, :, ir) = -sscale * h%Data(:, :, ir)
            case (6)
                h%Data(:, :, ir) = cmplx(reshape([tz(1), -tz(2), tz(2), tz(3)], [2, 2]), 0.0_dp, dp)
                s%Data(:, :, ir) = -sscale * h%Data(:, :, ir)
            case (7)
                h%Data(:, :, ir) = transpose(h%Data(:, :, 6))
                s%Data(:, :, ir) = transpose(s%Data(:, :, 6))
            end select
        end do
        h%has_meta = .true.
        h%fermi_ev = 0.0_dp
        h%energy_reference = 'absolute'
        call h%WriteFile(prefix // '_H.hsio')
        call s%WriteFile(prefix // '_S.hsio')
        open(newunit=u, file=prefix // '_struct.nml', status='replace', action='write')
        write(u, '(A)') '! Synthetic bulk in a surface-adapted cell (a3 tilted), NOT DFT.'
        write(u, '(A)') '&STRUCTURE'
        write(u, '(A)') '    LatticeVectors(1,:) = 3.0, 0.0, 0.0'
        write(u, '(A)') '    LatticeVectors(2,:) = 0.0, 3.0, 0.0'
        write(u, '(A,3(F8.4,:,","))') '    LatticeVectors(3,:) = ', a3
        write(u, '(A)') '/'
        write(u, '(A)') '&SPECIES'
        write(u, '(A)') '1'
        write(u, '(A,A)') species, ' 2'
        write(u, '(A)') '/'
        write(u, '(A)') '&ATOMS'
        write(u, '(A)') '1'
        write(u, '(A,A)') species, '  0.0  0.0  1.0'
        write(u, '(A)') '/'
        close(u)
    end subroutine write_material

end program gen_surfgreen_model
