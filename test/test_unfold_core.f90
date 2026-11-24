! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Testes do núcleo de unfolding (module unfold), com dados sintéticos em memória.
!>
!> Cadeia diatômica 1D (PC: A em x=0, B em x=0.5, a=1; hopping t1 A->B, t2 B->A+1),
!> SC de período 2 ao longo de x. Compara o motor de M NÃO-diagonal
!> M = [[2,1,0],[0,1,0],[0,0,1]] (det 2) contra o caminho M diagonal diag(2,1,1)
!> (mesma SC física ao longo da cadeia; y é desacoplado), com e sem
!> representantes de célula DESLOCADOS por um vetor da SC em um dos orbitais-PC.
!> A geometria (e portanto H(R)) é construída a partir dos rótulos reais.
program test_unfold_core
    use types,   only: dp
    use HS_io,   only: HSDataFormat
    use unfold,  only: UnfoldConfig, UnfoldResultFormat, UnfoldPerformance, RunUnfolding
    use testing, only: check, check_close, test_summary
    implicit none

    integer,  parameter :: NA = 4, NK = 5
    real(dp), parameter :: T1 = -1.0_dp, T2 = -0.6_dp, EA = 0.4_dp, EB = -0.4_dp
    ! átomos da SC: tipo (1=A, 2=B) e rótulo de célula n (posição x = n + 0.5*(tipo-1))
    integer,  parameter :: KIND_(NA)     = [1, 2, 1, 2]
    integer,  parameter :: LAB_REF(NA)   = [0, 0, 1, 1]
    integer,  parameter :: LAB_SHIFT(NA) = [0, 0, 1, -1]   ! B da 2ª célula: 1 - 2 = -1
    type(UnfoldResultFormat) :: ref, nd_plain, nd_shift, dg_shift
    real(dp)                 :: tol

    call run_case(.true.,  LAB_REF,   ref)
    call run_case(.false., LAB_REF,   nd_plain)
    call run_case(.false., LAB_SHIFT, nd_shift)
    call run_case(.true.,  LAB_SHIFT, dg_shift)

    tol = 1.0e-10_dp
    call compare(ref, nd_plain, tol, "M nao-diagonal, rotulos canonicos == M diagonal")
    call compare(ref, dg_shift, tol, "M diagonal, rotulos deslocados == rotulos canonicos")
    call compare(ref, nd_shift, tol, "M nao-diagonal, rotulos deslocados == M diagonal")

    call test_summary('test_unfold_core')

contains

    subroutine compare(a, b, tolerance, msg)
        type(UnfoldResultFormat), intent(in) :: a, b
        real(dp),                 intent(in) :: tolerance
        character(len=*),         intent(in) :: msg
        call check_close(maxval(abs(a%Energies - b%Energies)), 0.0_dp, tolerance, msg // " (E)")
        call check_close(maxval(abs(a%Weights  - b%Weights)),  0.0_dp, tolerance, msg // " (W)")
        call check_close(maxval(abs(b%SumCheck - 2.0_dp)),     0.0_dp, tolerance, msg // " (soma)")
    end subroutine compare

    subroutine run_case(diag, labels, res)
        logical,                  intent(in)  :: diag
        integer,                  intent(in)  :: labels(NA)
        type(UnfoldResultFormat), intent(out) :: res
        type(UnfoldConfig)       :: cfg
        type(HSDataFormat)       :: H, S
        type(UnfoldPerformance)  :: perf
        integer                  :: i, j, ir, ik, R
        real(dp)                 :: xi, xj, d

        H%nAtoms = NA;  H%NOrbs = NA;  H%Neighbors = 5
        allocate(H%Orbs(NA), source=1)
        allocate(H%Degenerecency(5), source=1)
        allocate(H%IndexList(3, 5), source=0)
        do ir = 1, 5
            H%IndexList(1, ir) = ir - 3          ! R = -2..2 ao longo de A1
        end do
        allocate(H%Data(NA, NA, 5), source=(0.0_dp, 0.0_dp))
        S = H
        S%Data = (0.0_dp, 0.0_dp)
        do i = 1, NA
            S%Data(i, i, 3) = (1.0_dp, 0.0_dp)
            H%Data(i, i, 3) = cmplx(merge(EA, EB, KIND_(i) == 1), 0.0_dp, dp)
        end do
        do ir = 1, 5
            R = ir - 3
            do i = 1, NA
                do j = 1, NA
                    xi = real(labels(i), dp) + 0.5_dp * real(KIND_(i) - 1, dp)
                    xj = real(labels(j), dp) + 0.5_dp * real(KIND_(j) - 1, dp)
                    d = xj + 2.0_dp * real(R, dp) - xi
                    if (abs(abs(d) - 0.5_dp) > 1.0e-9_dp) cycle
                    if ((KIND_(i) == 1) .eqv. (d > 0.0_dp)) then
                        H%Data(i, j, ir) = cmplx(T1, 0.0_dp, dp)
                    else
                        H%Data(i, j, ir) = cmplx(T2, 0.0_dp, dp)
                    end if
                end do
            end do
        end do

        if (diag) then
            cfg%M = reshape([2,0,0, 0,1,0, 0,0,1], [3,3])
            cfg%MDiag = [2, 1, 1]
            cfg%DiagonalM = .true.
        else
            cfg%M = transpose(reshape([2,1,0, 0,1,0, 0,0,1], [3,3]))
            cfg%MDiag = 0
            cfg%DiagonalM = .false.
        end if
        cfg%nCells = 2
        allocate(cfg%PCOrbit(NA), source=KIND_)
        allocate(cfg%CellIndex(3, NA), source=0)
        cfg%CellIndex(1, :) = labels
        cfg%PCPath%npaths = 2
        cfg%PCPath%total_kpoints = NK
        allocate(cfg%PCPath%FullPath(3, NK), source=0.0_dp)
        allocate(cfg%PCPath%xpoints(NK), source=0.0_dp)
        do ik = 1, NK
            cfg%PCPath%FullPath(1, ik) = 0.04_dp + 0.09_dp * real(ik - 1, dp)
        end do

        call RunUnfolding(cfg, H, S, res, perf)
    end subroutine run_case

end program test_unfold_core
