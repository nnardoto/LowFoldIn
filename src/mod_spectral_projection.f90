! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Blocos numéricos puros para pesos espectrais resolvidos por átomo/região
!> (docs/design/contracts.md §4): grupos de orbitais a partir de listas de átomos,
!> produto S·C, pesos de Mulliken e de Löwdin e um autoresolvedor hermitiano
!> generalizado que nunca para o programa.
!>
!> Posição no pipeline: chamado por ponto k pelo spectral_engine (`process_kpoint`)
!> depois de construídos H(k), S(k): autoproblema -> (S·C) -> pesos W(banda, grupo).
!> Sem E/S, sem estado de módulo; as rotinas da biblioteca nunca param (erros via ierr).
!>
!> Convenções:
!>   - As matrizes são hermitianas complex(dp) em armazenamento completo, ordem de
!>     colunas; aqui o LAPACK/BLAS referencia o triângulo superior.
!>   - Os autovetores são as colunas de C: C(mu, b) = coeficiente do orbital mu
!>     na banda b; as bandas crescem em energia e satisfazem C^† S C = I.
!>   - Índices de orbitais e de átomos são base 1; o átomo a possui os orbitais
!>     contíguos offset(a)+1 .. offset(a+1) (o próprio offset é base 0).
!>   - Os pesos têm sinal: nada é truncado nem passado por abs() (veja cada rotina).
!>   - O S passado às rotinas de peso deve ser o S(k) original, nunca o espaço de
!>     trabalho (destruído) do zhegv.
module spectral_projection
    use types, only: dp
    implicit none
    private

    public :: orbital_groups_t
    public :: atom_orbital_offsets, build_orbital_groups
    public :: overlap_times_vectors, mulliken_weights, lowdin_weights
    public :: generalized_eigensolve

    ! ==================================================================
    ! Tipos
    ! ==================================================================

    !> Grupos de índices de orbitais em layout CSR:
    !> grupo g = orb(start(g):start(g+1)-1). norb = total de orbitais do
    !> sistema (não o número armazenado em orb: grupos sobrepostos podem repetir
    !> orbitais, e os grupos não precisam cobrir todos eles).
    type :: orbital_groups_t
        integer              :: ngroup = 0, norb = 0
        integer, allocatable :: start(:)             ! (ngroup+1)
        integer, allocatable :: orb(:)
    end type orbital_groups_t

    ! ==================================================================
    ! INTERFACES EXPLÍCITAS PARA BLAS/LAPACK (double complex)
    ! (as assinaturas externas não têm intent, daí as cópias locais abaixo)
    ! ==================================================================
    interface
        subroutine zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            use types, only: dp
            character   :: side, uplo
            integer     :: m, n, lda, ldb, ldc
            complex(dp) :: alpha, beta
            complex(dp) :: a(lda, *), b(ldb, *), c(ldc, *)
        end subroutine zhemm

        subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
            use types, only: dp
            character   :: jobz, uplo
            integer     :: itype, n, lda, ldb, lwork, info
            complex(dp) :: a(lda, *), b(ldb, *), work(*)
            real(dp)    :: w(*), rwork(*)
        end subroutine zhegv

        subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
            use types, only: dp
            character   :: jobz, uplo
            integer     :: n, lda, lwork, info
            complex(dp) :: a(lda, *), work(*)
            real(dp)    :: w(*), rwork(*)
        end subroutine zheev

        subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            use types, only: dp
            character   :: transa, transb
            integer     :: m, n, k, lda, ldb, ldc
            complex(dp) :: alpha, beta
            complex(dp) :: a(lda, *), b(ldb, *), c(ldc, *)
        end subroutine zgemm
    end interface

contains

    ! ==================================================================
    ! Contabilidade de orbitais
    ! ==================================================================

    !> Deslocamentos de orbitais por átomo: offset(1) = 0, offset(a+1) = offset(a) + orbs(a).
    !> Os orbitais do átomo a são offset(a)+1 .. offset(a+1).
    !> offset deve ter tamanho nat+1 (nat = size(orbs_per_atom)).
    pure subroutine atom_orbital_offsets(orbs_per_atom, offset)
        integer, intent(in)  :: orbs_per_atom(:)
        integer, intent(out) :: offset(:)           ! (nat+1)
        integer              :: a

        offset(1) = 0
        do a = 1, size(orbs_per_atom)
            offset(a + 1) = offset(a) + orbs_per_atom(a)
        end do
    end subroutine atom_orbital_offsets

    !> Listas CSR genéricas de átomos -> listas de orbitais. O grupo g contém, para cada
    !> átomo a em list_atoms(list_start(g):list_start(g+1)-1) (nessa ordem), os
    !> orbitais offset(a)+1 .. offset(a)+orbs(a). Um átomo pode aparecer em vários
    !> grupos (grupos sobrepostos são permitidos).
    !> list_start(1:nlist+1) deve ser base 1 e monótono (list_start(nlist+1)-1
    !> = número de entradas usadas em list_atoms). Os índices de átomos não têm o
    !> intervalo verificado aqui: o chamador deve passar 1 <= a <= size(orbs_per_atom).
    !> groups%norb = total de orbitais do sistema. Nunca para o programa.
    subroutine build_orbital_groups(orbs_per_atom, nlist, list_start, list_atoms, groups)
        integer,                intent(in)  :: orbs_per_atom(:)
        integer,                intent(in)  :: nlist
        integer,                intent(in)  :: list_start(:)     ! (nlist+1)
        integer,                intent(in)  :: list_atoms(:)
        type(orbital_groups_t), intent(out) :: groups

        integer, allocatable :: offset(:)
        integer              :: nat, g, i, a, mu, pos, ntot

        nat = size(orbs_per_atom)
        allocate(offset(nat + 1))
        call atom_orbital_offsets(orbs_per_atom, offset)

        groups%ngroup = nlist
        groups%norb   = offset(nat + 1)
        allocate(groups%start(nlist + 1))

        ! 1. Primeira passada: tamanho de cada grupo -> ponteiros de início CSR.
        groups%start(1) = 1
        do g = 1, nlist
            ntot = 0
            do i = list_start(g), list_start(g + 1) - 1
                ntot = ntot + orbs_per_atom(list_atoms(i))
            end do
            groups%start(g + 1) = groups%start(g) + ntot
        end do

        ! 2. Segunda passada: expande cada átomo em seu intervalo contíguo de orbitais.
        allocate(groups%orb(groups%start(nlist + 1) - 1))
        pos = 0
        do g = 1, nlist
            do i = list_start(g), list_start(g + 1) - 1
                a = list_atoms(i)
                do mu = offset(a) + 1, offset(a) + orbs_per_atom(a)
                    pos = pos + 1
                    groups%orb(pos) = mu
                end do
            end do
        end do
    end subroutine build_orbital_groups

    ! ==================================================================
    ! Pesos
    ! ==================================================================

    !> SC = S·C com S hermitiana (armazenamento completo, triângulo superior referenciado).
    !> S e C não são modificados. S(n,n), C(n,nband), SC(n,nband); nband = n
    !> no motor spectral. Se n <= 0 ou nband <= 0, retorna com SC intacto.
    subroutine overlap_times_vectors(n, nband, S, C, SC)
        integer,     intent(in)  :: n, nband
        complex(dp), intent(in)  :: S(n, n)
        complex(dp), intent(in)  :: C(n, nband)
        complex(dp), intent(out) :: SC(n, nband)

        complex(dp), allocatable :: Swork(:, :), Cwork(:, :)

        if (n <= 0 .or. nband <= 0) return
        ! As cópias locais mantêm os argumentos intent(in) formalmente intactos na
        ! chamada externa (zhemm não modifica A nem B, mas sua interface não tem intent).
        allocate(Swork(n, n), Cwork(n, nband))
        Swork = S
        Cwork = C
        call zhemm('L', 'U', n, nband, (1.0_dp, 0.0_dp), Swork, n, Cwork, n, &
                   (0.0_dp, 0.0_dp), SC, n)
    end subroutine overlap_times_vectors

    !> Mulliken weights W(b,g) = Re Σ_{μ∈g} conj(C(μ,b)) · SC(μ,b).
    !> Valores negativos são mantidos como estão (sem truncar, sem abs): os pesos de
    !> Mulliken são populações, não probabilidades. SC deve ser S·C com o mesmo S
    !> usado no autoproblema (veja overlap_times_vectors). O grupo formado por
    !> todos os orbitais soma 1 por banda para C normalizado (C^† S C = I).
    subroutine mulliken_weights(n, nband, C, SC, groups, W)
        integer,                intent(in)  :: n, nband
        complex(dp),            intent(in)  :: C(n, nband), SC(n, nband)
        type(orbital_groups_t), intent(in)  :: groups
        real(dp),               intent(out) :: W(nband, groups%ngroup)

        integer  :: g, b, i, mu
        real(dp) :: acc

        do g = 1, groups%ngroup
            do b = 1, nband
                acc = 0.0_dp
                do i = groups%start(g), groups%start(g + 1) - 1
                    mu = groups%orb(i)
                    acc = acc + real(conjg(C(mu, b)) * SC(mu, b), dp)
                end do
                W(b, g) = acc
            end do
        end do
    end subroutine mulliken_weights

    !> Löwdin weights W(b,g) = Σ_{μ∈g} |(S^{1/2} C)(μ,b)|².
    !> Não negativos por construção (base ortogonalizada).
    !> S^{1/2} = U diag(sqrt(w)) U^† obtido de zheev('V','U') sobre uma CÓPIA de S.
    !> smin = menor autovalor de S. Se smin <= tol_s_min: ierr = 1 e W fica
    !> indefinido (os autovalores nunca são limitados nem truncados).
    !> Falha do zheev: ierr = info (< 0, argumento ilegal) ou 1000 + info (> 0,
    !> sem convergência), ambos distintos de 1. S e C não são modificados.
    subroutine lowdin_weights(n, nband, S, C, groups, W, tol_s_min, smin, ierr)
        integer,                intent(in)  :: n, nband
        complex(dp),            intent(in)  :: S(n, n)
        complex(dp),            intent(in)  :: C(n, nband)
        type(orbital_groups_t), intent(in)  :: groups
        real(dp),               intent(out) :: W(nband, groups%ngroup)
        real(dp),               intent(in)  :: tol_s_min
        real(dp),               intent(out) :: smin
        integer,                intent(out) :: ierr

        complex(dp), allocatable :: U(:, :), US(:, :), Shalf(:, :), Cwork(:, :), Ct(:, :), work(:)
        complex(dp)              :: wquery(1)
        real(dp),    allocatable :: ev(:), rwork(:)
        real(dp)                 :: acc
        integer                  :: lwork, info, i, g, b, mu

        ierr = 0
        smin = 0.0_dp
        if (n <= 0 .or. nband <= 0) return

        allocate(U(n, n), ev(n), rwork(max(1, 3 * n - 2)))
        U = S

        ! 1. Diagonaliza S (cópia em U): S = U diag(ev) U^†.
        ! Consulta de espaço de trabalho do LAPACK (lwork = -1 devolve o tamanho ótimo em wquery).
        call zheev('V', 'U', n, U, n, ev, wquery, -1, rwork, info)
        if (info /= 0) then
            ierr = merge(info, 1000 + info, info < 0)
            return
        end if
        lwork = max(1, 2 * n - 1, int(real(wquery(1), dp)))
        allocate(work(lwork))
        call zheev('V', 'U', n, U, n, ev, work, lwork, rwork, info)
        if (info /= 0) then
            ierr = merge(info, 1000 + info, info < 0)
            return
        end if

        ! 2. Guarda de positividade definida: os autovalores crescem, mas usa-se
        ! minval para que a verificação não dependa da ordenação. O chamador aborta.
        smin = minval(ev)
        if (smin <= tol_s_min) then
            ierr = 1
            return
        end if

        ! 3. S^{1/2} = (U diag(sqrt(ev))) U^H  (escala das colunas, depois zgemm 'N','C')
        allocate(US(n, n), Shalf(n, n))
        do i = 1, n
            US(:, i) = U(:, i) * sqrt(ev(i))
        end do
        call zgemm('N', 'C', n, n, n, (1.0_dp, 0.0_dp), US, n, U, n, (0.0_dp, 0.0_dp), Shalf, n)

        ! 4. Ct = S^{1/2} C (a cópia local de C mantém o argumento intent(in) intacto)
        allocate(Cwork(n, nband), Ct(n, nband))
        Cwork = C
        call zgemm('N', 'N', n, nband, n, (1.0_dp, 0.0_dp), Shalf, n, Cwork, n, &
                   (0.0_dp, 0.0_dp), Ct, n)

        ! 5. Soma |Ct|^2 sobre os orbitais de cada grupo.
        do g = 1, groups%ngroup
            do b = 1, nband
                acc = 0.0_dp
                do i = groups%start(g), groups%start(g + 1) - 1
                    mu = groups%orb(i)
                    acc = acc + real(Ct(mu, b), dp)**2 + aimag(Ct(mu, b))**2
                end do
                W(b, g) = acc
            end do
        end do
    end subroutine lowdin_weights

    ! ==================================================================
    ! Autoresolvedor
    ! ==================================================================

    !> Autoproblema hermitiano generalizado H C = S C diag(E) via
    !> zhegv(1,'V','U') sobre cópias internas de H e S. Autovalores crescentes,
    !> C^† S C = I. Nunca para o programa: ierr = info do zhegv (ou da consulta de espaço):
    !> info < 0 argumento ilegal, 0 < info <= n o autoresolvedor não convergiu,
    !> info > n S não é definida positiva (convenção do LAPACK). Com ierr /= 0, E e C
    !> ficam indefinidos. H e S não são modificados.
    subroutine generalized_eigensolve(n, H, S, E, C, ierr)
        integer,     intent(in)  :: n
        complex(dp), intent(in)  :: H(n, n), S(n, n)
        real(dp),    intent(out) :: E(n)
        complex(dp), intent(out) :: C(n, n)
        integer,     intent(out) :: ierr

        complex(dp), allocatable :: Bwork(:, :), work(:)
        complex(dp)              :: wquery(1)
        real(dp),    allocatable :: rwork(:)
        integer                  :: lwork, info

        ierr = 0
        if (n <= 0) return

        ! O zhegv sobrescreve A com os autovetores e B com seu fator de Cholesky,
        ! então ambos são passados como cópias (C serve também de cópia de A).
        C = H
        allocate(Bwork(n, n))
        Bwork = S
        allocate(rwork(max(1, 3 * n - 2)))

        ! Consulta de espaço de trabalho do LAPACK (lwork = -1)
        call zhegv(1, 'V', 'U', n, C, n, Bwork, n, E, wquery, -1, rwork, info)
        if (info /= 0) then
            ierr = info
            return
        end if
        lwork = max(1, 2 * n - 1, int(real(wquery(1), dp)))
        allocate(work(lwork))

        call zhegv(1, 'V', 'U', n, C, n, Bwork, n, E, work, lwork, rwork, info)
        ierr = info
    end subroutine generalized_eigensolve

end module spectral_projection
