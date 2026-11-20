! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Álgebra linear complexa densa usada pelo módulo de funções de Green.
!>
!> Interfaces finas sobre LAPACK/BLAS que devolvem status em vez de parar:
!>   - lu_factor / lu_solve: fatoração LU (zgetrf) e solução com vários lados
!>     direitos (zgetrs), com transposição opcional;
!>   - lu_solve_identity: X = A^{-1} pela solução A X = I (sem inversa explícita
!>     fora do LAPACK; usada só para blocos pequenos de G);
!>   - zmm: C = alpha*op(A)*op(B) + beta*C via zgemm;
!>   - herm_check: max |A - A^H|;
!>   - herm_eigvals: autovalores de uma matriz hermitiana (zheev), para
!>     diagnósticos como o menor autovalor de S.
!>
!> Nenhuma rotina usa `stop`; erros são devolvidos em ierr/errmsg.
module zlinalg
    use types, only: dp
    implicit none
    private

    public :: lu_factor, lu_solve, lu_solve_identity, zmm, herm_check, herm_eigvals
    public :: identity_matrix, herm_sqrt, diag_of_product

    complex(dp), parameter :: ZONE = (1.0_dp, 0.0_dp), ZZERO = (0.0_dp, 0.0_dp)

    interface
        subroutine zgetrf(m, n, a, lda, ipiv, info)
            use types, only: dp
            integer     :: m, n, lda, info
            integer     :: ipiv(*)
            complex(dp) :: a(lda, *)
        end subroutine zgetrf

        subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            use types, only: dp
            character   :: trans
            integer     :: n, nrhs, lda, ldb, info
            integer     :: ipiv(*)
            complex(dp) :: a(lda, *), b(ldb, *)
        end subroutine zgetrs

        subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            use types, only: dp
            character   :: transa, transb
            integer     :: m, n, k, lda, ldb, ldc
            complex(dp) :: alpha, beta
            complex(dp) :: a(lda, *), b(ldb, *), c(ldc, *)
        end subroutine zgemm

        subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
            use types, only: dp
            character   :: jobz, uplo
            integer     :: n, lda, lwork, info
            complex(dp) :: a(lda, *), work(*)
            real(dp)    :: w(*), rwork(*)
        end subroutine zheev
    end interface

contains

    !> Fatoração LU in-place: A = P L U. ipiv é alocado aqui.
    !> ierr = 1 (dimensões), 2 (matriz singular: U(info,info) = 0), 3 (argumento inválido).
    subroutine lu_factor(a, ipiv, ierr, errmsg)
        complex(dp), contiguous, intent(inout) :: a(:, :)
        integer, allocatable,    intent(out)   :: ipiv(:)
        integer,                 intent(out)   :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: n, info

        ierr = 0
        errmsg = ''
        n = size(a, 1)
        allocate(ipiv(max(n, 1)))
        if (size(a, 2) /= n) then
            ierr = 1; errmsg = 'lu_factor: matrix must be square'; return
        end if
        if (n == 0) return
        call zgetrf(n, n, a, n, ipiv, info)
        if (info > 0) then
            ierr = 2; errmsg = 'lu_factor: singular matrix'
        else if (info < 0) then
            ierr = 3; errmsg = 'lu_factor: invalid argument to zgetrf'
        end if
    end subroutine lu_factor

    !> Resolve op(A) X = B com os fatores de lu_factor; B é sobrescrito por X.
    !> trans = 'N' (padrão), 'T' ou 'C'.
    subroutine lu_solve(a_lu, ipiv, b, ierr, errmsg, trans)
        complex(dp), contiguous, intent(in)    :: a_lu(:, :)
        integer,                 intent(in)    :: ipiv(:)
        complex(dp), contiguous, intent(inout) :: b(:, :)
        integer,                 intent(out)   :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=1), optional, intent(in) :: trans
        character(len=1) :: t
        integer :: n, nrhs, info

        ierr = 0
        errmsg = ''
        t = 'N'
        if (present(trans)) t = trans
        n = size(a_lu, 1)
        nrhs = size(b, 2)
        if (size(a_lu, 2) /= n .or. size(b, 1) /= n .or. size(ipiv) < n) then
            ierr = 1; errmsg = 'lu_solve: inconsistent dimensions'; return
        end if
        if (n == 0 .or. nrhs == 0) return
        call zgetrs(t, n, nrhs, a_lu, n, ipiv, b, n, info)
        if (info /= 0) then
            ierr = 3; errmsg = 'lu_solve: invalid argument to zgetrs'
        end if
    end subroutine lu_solve

    !> X = A^{-1} a partir dos fatores LU (solução com o lado direito identidade).
    subroutine lu_solve_identity(a_lu, ipiv, x, ierr, errmsg)
        complex(dp), contiguous, intent(in)  :: a_lu(:, :)
        integer,                 intent(in)  :: ipiv(:)
        complex(dp), contiguous, intent(out) :: x(:, :)
        integer,                 intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        call identity_matrix(x)
        call lu_solve(a_lu, ipiv, x, ierr, errmsg)
    end subroutine lu_solve_identity

    !> C = alpha * op(A) * op(B) + beta * C, op = 'N', 'T' ou 'C'.
    !> alpha = 1 e beta = 0 por padrão. As dimensões são deduzidas de A, B e C;
    !> não há verificação além da que o zgemm faz (uso interno).
    subroutine zmm(transa, transb, a, b, c, alpha, beta)
        character(len=1),        intent(in)    :: transa, transb
        complex(dp), contiguous, intent(in)    :: a(:, :), b(:, :)
        complex(dp), contiguous, intent(inout) :: c(:, :)
        complex(dp), optional,   intent(in)    :: alpha, beta
        complex(dp) :: al, be
        integer :: m, n, k
        al = ZONE
        be = ZZERO
        if (present(alpha)) al = alpha
        if (present(beta))  be = beta
        m = size(c, 1)
        n = size(c, 2)
        if (transa == 'N' .or. transa == 'n') then
            k = size(a, 2)
        else
            k = size(a, 1)
        end if
        if (m == 0 .or. n == 0) return
        call zgemm(transa, transb, m, n, k, al, a, size(a, 1), b, size(b, 1), be, c, m)
    end subroutine zmm

    !> max_ij |A_ij - conj(A_ji)|.
    pure real(dp) function herm_check(a) result(dev)
        complex(dp), intent(in) :: a(:, :)
        integer :: i, j
        dev = 0.0_dp
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                dev = max(dev, abs(a(i, j) - conjg(a(j, i))))
            end do
        end do
    end function herm_check

    !> Autovalores (ordem crescente) de uma matriz hermitiana; a não é alterada.
    subroutine herm_eigvals(a, w, ierr, errmsg)
        complex(dp), intent(in)  :: a(:, :)
        real(dp),    intent(out) :: w(:)
        integer,     intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        complex(dp), allocatable :: tmp(:, :), work(:)
        real(dp), allocatable :: rwork(:)
        integer :: n, info, lwork

        ierr = 0
        errmsg = ''
        n = size(a, 1)
        if (size(a, 2) /= n .or. size(w) /= n) then
            ierr = 1; errmsg = 'herm_eigvals: inconsistent dimensions'; return
        end if
        if (n == 0) return
        lwork = max(1, 2 * n)
        allocate(tmp(n, n), work(lwork), rwork(max(1, 3 * n - 2)))
        tmp = a
        call zheev('N', 'U', n, tmp, n, w, work, lwork, rwork, info)
        if (info /= 0) then
            ierr = 2; errmsg = 'herm_eigvals: zheev failed'
        end if
    end subroutine herm_eigvals

    !> S^{1/2} de uma matriz hermitiana definida positiva (zheev com vetores).
    !> smin = menor autovalor; ierr = 2 se smin <= 0.
    subroutine herm_sqrt(a, ahalf, smin, ierr, errmsg)
        complex(dp), intent(in)  :: a(:, :)
        complex(dp), contiguous, intent(out) :: ahalf(:, :)
        real(dp),    intent(out) :: smin
        integer,     intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        complex(dp), allocatable :: v(:, :), vs(:, :), work(:)
        real(dp), allocatable :: w(:), rwork(:)
        integer :: n, info, lwork, j

        ierr = 0
        errmsg = ''
        smin = 0.0_dp
        n = size(a, 1)
        if (size(a, 2) /= n .or. any(shape(ahalf) /= [n, n])) then
            ierr = 1; errmsg = 'herm_sqrt: inconsistent dimensions'; return
        end if
        if (n == 0) return
        lwork = max(1, 2 * n)
        allocate(v(n, n), vs(n, n), w(n), work(lwork), rwork(max(1, 3 * n - 2)))
        v = a
        call zheev('V', 'U', n, v, n, w, work, lwork, rwork, info)
        if (info /= 0) then
            ierr = 2; errmsg = 'herm_sqrt: zheev failed'; return
        end if
        smin = w(1)
        if (smin <= 0.0_dp) then
            ierr = 2; errmsg = 'herm_sqrt: matrix is not positive definite'; return
        end if
        do j = 1, n
            vs(:, j) = v(:, j) * sqrt(w(j))
        end do
        call zmm('N', 'C', vs, v, ahalf)
    end subroutine herm_sqrt

    !> d(i) = d(i) + sum_j A(i,j) B(j,i): diagonal de A*B acumulada, sem formar o produto.
    pure subroutine diag_of_product(a, b, d)
        complex(dp), intent(in)    :: a(:, :), b(:, :)
        complex(dp), intent(inout) :: d(:)
        integer :: i, j
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                d(i) = d(i) + a(i, j) * b(j, i)
            end do
        end do
    end subroutine diag_of_product

    pure subroutine identity_matrix(x)
        complex(dp), intent(out) :: x(:, :)
        integer :: i
        x = ZZERO
        do i = 1, min(size(x, 1), size(x, 2))
            x(i, i) = ZONE
        end do
    end subroutine identity_matrix

end module zlinalg
