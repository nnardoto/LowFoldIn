! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Regressão do unfolding com M NÃO-diagonal e NÃO-simétrica, em rede 2D real.
!>
!> Modelo tight-binding tipo hBN (rede hexagonal, 2 sítios com on-site ±DELTA,
!> hopping T e sobreposição SOV entre primeiros vizinhos, base não ortogonal).
!> A supercélula √3×√3 R30° tem M = [[1,1,0],[-1,2,0],[0,0,1]] (linhas), que
!> NÃO é simétrica: erros de transposição no mapeamento k -> K ou na tabela de
!> grupo Z^3/M.Z^3 passam despercebidos com M simétrica ou cadeias 1D.
!>
!> Como a SC é perfeita, para cada k da PC os estados da SC com peso não nulo
!> reproduzem exatamente as bandas da PC. Comparamos somas invariantes sob a
!> escolha de base em subespaços degenerados:
!>   sum_N W_N = 2,  sum_N W_N E_N = E_1 + E_2,  sum_N W_N E_N^2 = E_1^2 + E_2^2,
!> com E_1, E_2 do problema generalizado 2x2 da PC (resolvido analiticamente).
program test_unfold_nonsym
    use types,   only: dp
    use HS_io,   only: HSDataFormat
    use core,    only: compute_Hk_at_kpoint
    use unfold,  only: UnfoldConfig, UnfoldResultFormat, UnfoldPerformance, RunUnfolding
    use testing, only: check_close, test_summary
    implicit none

    real(dp), parameter :: A = 2.5_dp, T = -2.6_dp, SOV = 0.1_dp, DELTA = 2.3_dp
    integer,  parameter :: NK = 13, MAXR = 2, NRMAX = (2*MAXR + 1)**2

    real(dp) :: a1(3), a2(3), pcpos(3,2)
    integer  :: M(3,3)

    type(HSDataFormat)       :: Hpc, Spc, Hsc, Ssc
    type(UnfoldConfig)       :: cfg
    type(UnfoldResultFormat) :: res
    type(UnfoldPerformance)  :: perf
    integer,  allocatable    :: kindsc(:), cellsc(:,:)
    complex(dp)              :: hk(2,2), sk(2,2)
    real(dp)                 :: kpt(3), qa, qb, qc, s1, s2, w1, w2, err0, err1, err2
    integer                  :: ik, nsc

    a1 = [A, 0.0_dp, 0.0_dp]
    a2 = [A/2, A*sqrt(3.0_dp)/2, 0.0_dp]
    pcpos(:,1) = 0.0_dp
    pcpos(:,2) = (a1 + a2) / 3
    M = transpose(reshape([1,1,0, -1,2,0, 0,0,1], [3,3]))   ! linhas de M

    ! --- PC (identidade) e SC √3 a partir das mesmas regras de vizinhança ---
    call build([1,0,0, 0,1,0, 0,0,1], Hpc, Spc)
    call build([1,1,0, -1,2,0, 0,0,1], Hsc, Ssc)
    call check_close(real(nsc, dp), 6.0_dp, 0.5_dp, "SC sqrt3 tem 6 atomos")

    ! --- configuração do unfolding (rótulos conhecidos por construção) ---
    cfg%M = M
    cfg%MDiag = 0
    cfg%DiagonalM = .false.
    cfg%nCells = 3
    allocate(cfg%PCOrbit(nsc), source=kindsc)
    allocate(cfg%CellIndex(3, nsc), source=cellsc)
    cfg%PCPath%npaths = 2
    cfg%PCPath%total_kpoints = NK
    allocate(cfg%PCPath%FullPath(3, NK), cfg%PCPath%xpoints(NK))
    cfg%PCPath%xpoints = 0.0_dp
    do ik = 1, NK          ! segmento genérico (fora de linhas de alta simetria)
        cfg%PCPath%FullPath(:, ik) = [0.03_dp, 0.11_dp, 0.0_dp] + &
                                     real(ik - 1, dp) / real(NK - 1, dp) * [0.61_dp, 0.17_dp, 0.0_dp]
    end do

    call RunUnfolding(cfg, Hsc, Ssc, res, perf)

    err0 = 0; err1 = 0; err2 = 0
    do ik = 1, NK
        kpt = cfg%PCPath%FullPath(:, ik)
        call compute_Hk_at_kpoint(Hpc%Data, Hpc%IndexList, hk, kpt)
        call compute_Hk_at_kpoint(Spc%Data, Spc%IndexList, sk, kpt)
        ! det(H - E S) = 0, com S_11 = S_22 = 1:  qa E^2 - qb E + qc = 0
        qa = 1.0_dp - abs(sk(1,2))**2
        qb = real(hk(1,1) + hk(2,2), dp) - 2.0_dp * real(hk(1,2) * conjg(sk(1,2)), dp)
        qc = real(hk(1,1) * hk(2,2), dp) - abs(hk(1,2))**2
        s1 = qb / qa                       ! E_1 + E_2
        s2 = s1**2 - 2.0_dp * qc / qa      ! E_1^2 + E_2^2
        w1 = sum(res%Weights(:, ik) * res%Energies(:, ik))
        w2 = sum(res%Weights(:, ik) * res%Energies(:, ik)**2)
        err0 = max(err0, abs(sum(res%Weights(:, ik)) - 2.0_dp))
        err1 = max(err1, abs(w1 - s1))
        err2 = max(err2, abs(w2 - s2))
    end do
    call check_close(err0, 0.0_dp, 1.0e-10_dp, "regra de soma sum W = p (M nao simetrica)")
    call check_close(err1, 0.0_dp, 1.0e-8_dp,  "sum W E = soma das bandas da PC (M nao simetrica)")
    call check_close(err2, 0.0_dp, 1.0e-7_dp,  "sum W E^2 = soma dos quadrados da PC (M nao simetrica)")

    call test_summary('test_unfold_nonsym')

contains

    !> Monta H(R), S(R) da supercélula com linhas de M dadas em `mrows`.
    !> Para a SC, também preenche kindsc/cellsc/possc/nsc (rótulos exatos).
    subroutine build(mrows, H, S)
        integer,            intent(in)  :: mrows(9)
        type(HSDataFormat), intent(out) :: H, S
        integer     :: Ms(3,3), i, j, k, nR, R(3), n1, n2, p, q
        real(dp)    :: A1v(3), A2v(3), Minv(2,2), det, f(2), d(3), dist
        integer,     allocatable :: kk(:), cc(:,:), Rl(:,:)
        real(dp),    allocatable :: pp(:,:)
        complex(dp), allocatable :: Hd(:,:,:), Sd(:,:,:)
        logical     :: used

        Ms = transpose(reshape(mrows, [3,3]))
        A1v = Ms(1,1)*a1 + Ms(1,2)*a2
        A2v = Ms(2,1)*a1 + Ms(2,2)*a2
        det = real(Ms(1,1)*Ms(2,2) - Ms(1,2)*Ms(2,1), dp)
        Minv = reshape([real(Ms(2,2),dp), -real(Ms(2,1),dp), -real(Ms(1,2),dp), real(Ms(1,1),dp)], [2,2]) / det

        allocate(kk(64), cc(3,64), pp(3,64))
        k = 0
        do n2 = -6, 6
            do n1 = -6, 6
                f(1) = n1*Minv(1,1) + n2*Minv(2,1)
                f(2) = n1*Minv(1,2) + n2*Minv(2,2)
                if (any(f < -1.0e-9_dp) .or. any(f >= 1.0_dp - 1.0e-9_dp)) cycle
                do i = 1, 2
                    k = k + 1
                    kk(k) = i; cc(:,k) = [n1, n2, 0]
                    pp(:,k) = pcpos(:,i) + n1*a1 + n2*a2
                end do
            end do
        end do

        allocate(Hd(k,k,NRMAX), Sd(k,k,NRMAX), Rl(3,NRMAX))
        Hd = 0; Sd = 0; nR = 0
        do j = -MAXR, MAXR
            do i = -MAXR, MAXR
                R = [i, j, 0]
                used = all(R == 0)
                do p = 1, k
                    do q = 1, k
                        d = pp(:,q) + R(1)*A1v + R(2)*A2v - pp(:,p)
                        dist = norm2(d)
                        if (all(R == 0) .and. p == q) then
                            Hd(p,q,nR+1) = merge(DELTA, -DELTA, kk(p) == 1)
                            Sd(p,q,nR+1) = 1.0_dp
                            used = .true.
                        else if (abs(dist - A/sqrt(3.0_dp)) < 1.0e-6_dp) then
                            Hd(p,q,nR+1) = T
                            Sd(p,q,nR+1) = SOV
                            used = .true.
                        end if
                    end do
                end do
                if (used) then
                    nR = nR + 1; Rl(:, nR) = R
                end if
            end do
        end do

        call init(H, Hd, Rl, nR, k); call init(S, Sd, Rl, nR, k)
        if (k > 2) then
            nsc = k
            allocate(kindsc(k), source=kk(1:k))
            allocate(cellsc(3,k), source=cc(:,1:k))
        end if

    end subroutine build

    subroutine init(X, D, Rl, nR, n)
        type(HSDataFormat), intent(out) :: X
        complex(dp),        intent(in)  :: D(:,:,:)
        integer,            intent(in)  :: Rl(:,:), nR, n
        X%Title = "test_unfold_nonsym"
        X%Neighbors = nR; X%nAtoms = n; X%NOrbs = n
        allocate(X%Degenerecency(nR), source=1)
        allocate(X%Orbs(n), source=1)
        allocate(X%IndexList(3,nR), source=Rl(:,1:nR))
        allocate(X%Data(n,n,nR), source=D(:,:,1:nR))
    end subroutine init

end program test_unfold_nonsym
