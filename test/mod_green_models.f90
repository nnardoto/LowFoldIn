! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Modelos sintéticos para os testes do módulo de Green: H(R), S(R) hermitianos
!> (H(-R) = H(R)^dagger) com valores pseudoaleatórios determinísticos, e
!> estruturas correspondentes.
module green_models
    use types, only: dp
    use HS_io, only: HSDataFormat
    use structure_io, only: structure_t
    implicit none
    private

    public :: make_model, make_structure, rnd, set_seed, read_cmat, slab_eigen

    integer, save :: state = 12345

contains

    subroutine set_seed(s)
        integer, intent(in) :: s
        state = s
    end subroutine set_seed

    !> Número pseudoaleatório em [-1, 1) (LCG de 31 bits).
    real(dp) function rnd()
        state = modulo(1103515245 * state + 12345, 2147483647)
        rnd = 2.0_dp * real(state, dp) / 2147483647.0_dp - 1.0_dp
    end function rnd

    !> H(R) e S(R) para R em [-r(1), r(1)] x [-r(2), r(2)] x [-r(3), r(3)].
    !> amp_h, amp_s: amplitudes dos termos fora do sítio (decaem como e^{-|R|_1});
    !> amp_s = 0 dá S = I. cplx: valores complexos. orbs: orbitais por átomo.
    subroutine make_model(orbs, r, amp_h, amp_s, cplx, seed, hh, hs)
        integer,  intent(in) :: orbs(:), r(3), seed
        real(dp), intent(in) :: amp_h, amp_s
        logical,  intent(in) :: cplx
        type(HSDataFormat), intent(out) :: hh, hs
        complex(dp), allocatable :: xh(:, :, :), xs(:, :, :)
        integer :: n, nr, i, j, k, ir, jr, i1, i2, i3
        real(dp) :: dec

        call set_seed(seed)
        n = sum(orbs)
        nr = product(2 * r + 1)
        allocate(xh(n, n, nr), xs(n, n, nr))
        hh%nAtoms = size(orbs); hh%NOrbs = n; hh%Neighbors = nr
        allocate(hh%Orbs(size(orbs)), hh%IndexList(3, nr), hh%Degenerecency(nr), hh%Data(n, n, nr))
        hh%Orbs = orbs
        hh%Degenerecency = 1
        hh%Title = 'synthetic H'
        ir = 0
        do i3 = -r(3), r(3)
            do i2 = -r(2), r(2)
                do i1 = -r(1), r(1)
                    ir = ir + 1
                    hh%IndexList(:, ir) = [i1, i2, i3]
                    dec = exp(-real(abs(i1) + abs(i2) + abs(i3), dp))
                    do j = 1, n
                        do i = 1, n
                            xh(i, j, ir) = amp_h * dec * cval()
                            xs(i, j, ir) = amp_s * dec * cval()
                        end do
                    end do
                end do
            end do
        end do
        ! H(R) = (X(R) + X(-R)^dagger)/2
        do ir = 1, nr
            jr = find_r(-hh%IndexList(:, ir))
            hh%Data(:, :, ir) = 0.5_dp * (xh(:, :, ir) + conjg(transpose(xh(:, :, jr))))
        end do
        hs = hh
        hs%Title = 'synthetic S'
        do ir = 1, nr
            jr = find_r(-hh%IndexList(:, ir))
            hs%Data(:, :, ir) = 0.5_dp * (xs(:, :, ir) + conjg(transpose(xs(:, :, jr))))
        end do
        ir = find_r([0, 0, 0])
        do k = 1, n
            hh%Data(k, k, ir) = hh%Data(k, k, ir) + cmplx(-1.0_dp + 2.0_dp * real(k - 1, dp) / &
                                                       real(max(n - 1, 1), dp), 0.0_dp, dp)
            hs%Data(k, k, ir) = hs%Data(k, k, ir) + (1.0_dp, 0.0_dp)
        end do

    contains

        complex(dp) function cval()
            real(dp) :: re, im
            re = rnd()
            im = 0.0_dp
            if (cplx) im = rnd()
            cval = cmplx(re, im, dp)
        end function cval

        integer function find_r(rv)
            integer, intent(in) :: rv(3)
            integer :: q
            find_r = 0
            do q = 1, nr
                if (all(hh%IndexList(:, q) == rv)) then
                    find_r = q
                    return
                end if
            end do
        end function find_r
    end subroutine make_model

    !> Estrutura com rede `lattice` (linhas a_i) e nat átomos em posições
    !> pseudoaleatórias dentro da célula. Rótulos prefix//índice.
    subroutine make_structure(lattice, nat, prefix, seed, st)
        real(dp),          intent(in)  :: lattice(3, 3)
        integer,           intent(in)  :: nat, seed
        character(len=*),  intent(in)  :: prefix
        type(structure_t), intent(out) :: st
        integer :: i
        real(dp) :: f(3)
        call set_seed(seed)
        st%nat = nat
        st%lattice = lattice
        allocate(st%cart(3, nat), st%label(nat))
        do i = 1, nat
            f = 0.5_dp + 0.4_dp * [rnd(), rnd(), rnd()]
            st%cart(:, i) = matmul(f, lattice)
            write(st%label(i), '(A,I0)') prefix, i
        end do
    end subroutine make_structure

    !> Lê n*n números "re im" (ordem de colunas) de uma unidade aberta.
    subroutine read_cmat(u, a)
        integer,     intent(in)  :: u
        complex(dp), intent(out) :: a(:, :)
        integer :: i, j
        real(dp) :: re, im
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                read(u, *) re, im
                a(i, j) = cmplx(re, im, dp)
            end do
        end do
    end subroutine read_cmat

    !> Autovalores e autovetores S-ortonormais de H c = E S c (zhegv).
    subroutine slab_eigen(h, s, e, c)
        complex(dp), intent(in)  :: h(:, :), s(:, :)
        real(dp),    intent(out) :: e(:)
        complex(dp), intent(out) :: c(:, :)
        interface
            subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
                use types, only: dp
                character   :: jobz, uplo
                integer     :: itype, n, lda, ldb, lwork, info
                complex(dp) :: a(lda, *), b(ldb, *), work(*)
                real(dp)    :: w(*), rwork(*)
            end subroutine zhegv
        end interface
        complex(dp), allocatable :: b(:, :), work(:)
        real(dp), allocatable :: rwork(:)
        integer :: n, info
        n = size(h, 1)
        allocate(b(n, n), work(4 * n), rwork(3 * n))
        c = h
        b = s
        call zhegv(1, 'V', 'U', n, c, n, b, n, e, work, 4 * n, rwork, info)
        if (info /= 0) error stop 'slab_eigen: zhegv failed'
    end subroutine slab_eigen

end module green_models
