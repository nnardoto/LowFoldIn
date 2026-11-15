! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Pure numerical building blocks for atom/region-resolved spectral weights
!> (contracts §4): orbital groups from atom lists, S·C product, Mulliken
!> weights and a non-stopping generalized Hermitian eigensolver.
module spectral_projection
    use types, only: dp
    implicit none
    private

    public :: orbital_groups_t
    public :: atom_orbital_offsets, build_orbital_groups
    public :: overlap_times_vectors, mulliken_weights
    public :: generalized_eigensolve

    !> Groups of orbital indices in CSR layout:
    !> group g = orb(start(g):start(g+1)-1).
    type :: orbital_groups_t
        integer :: ngroup = 0, norb = 0
        integer, allocatable :: start(:)   ! (ngroup+1)
        integer, allocatable :: orb(:)
    end type orbital_groups_t

    ! ==================================================================
    ! EXPLICIT INTERFACES FOR BLAS/LAPACK (double complex)
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
    end interface

contains

    !> Orbital offsets per atom: offset(1) = 0, offset(a+1) = offset(a) + orbs(a).
    !> Orbitals of atom a are offset(a)+1 .. offset(a+1).
    pure subroutine atom_orbital_offsets(orbs_per_atom, offset)
        integer, intent(in)  :: orbs_per_atom(:)
        integer, intent(out) :: offset(:)          ! (nat+1)
        integer :: a

        offset(1) = 0
        do a = 1, size(orbs_per_atom)
            offset(a + 1) = offset(a) + orbs_per_atom(a)
        end do
    end subroutine atom_orbital_offsets

    !> Generic CSR atom lists -> orbital lists. Group g contains, for each atom
    !> a in list_atoms(list_start(g):list_start(g+1)-1) (in that order), the
    !> orbitals offset(a)+1 .. offset(a)+orbs(a). Atoms may appear in several
    !> groups (overlapping groups are allowed).
    subroutine build_orbital_groups(orbs_per_atom, nlist, list_start, list_atoms, groups)
        integer,                intent(in)  :: orbs_per_atom(:)
        integer,                intent(in)  :: nlist
        integer,                intent(in)  :: list_start(:)     ! (nlist+1)
        integer,                intent(in)  :: list_atoms(:)
        type(orbital_groups_t), intent(out) :: groups

        integer, allocatable :: offset(:)
        integer :: nat, g, i, a, mu, pos, ntot

        nat = size(orbs_per_atom)
        allocate(offset(nat + 1))
        call atom_orbital_offsets(orbs_per_atom, offset)

        groups%ngroup = nlist
        groups%norb   = offset(nat + 1)
        allocate(groups%start(nlist + 1))

        groups%start(1) = 1
        do g = 1, nlist
            ntot = 0
            do i = list_start(g), list_start(g + 1) - 1
                ntot = ntot + orbs_per_atom(list_atoms(i))
            end do
            groups%start(g + 1) = groups%start(g) + ntot
        end do

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

    !> SC = S·C with S Hermitian (full storage, upper triangle referenced).
    !> S and C are not modified.
    subroutine overlap_times_vectors(n, nband, S, C, SC)
        integer,     intent(in)  :: n, nband
        complex(dp), intent(in)  :: S(n, n)
        complex(dp), intent(in)  :: C(n, nband)
        complex(dp), intent(out) :: SC(n, nband)

        complex(dp), allocatable :: Swork(:, :), Cwork(:, :)

        if (n <= 0 .or. nband <= 0) return
        ! Local copies keep the intent(in) arguments formally untouched for the
        ! external call (zhemm does not modify A or B, but its interface has no intent).
        allocate(Swork(n, n), Cwork(n, nband))
        Swork = S
        Cwork = C
        call zhemm('L', 'U', n, nband, (1.0_dp, 0.0_dp), Swork, n, Cwork, n, &
                   (0.0_dp, 0.0_dp), SC, n)
    end subroutine overlap_times_vectors

    !> Mulliken weights W(b,g) = Re Σ_{μ∈g} conj(C(μ,b)) · SC(μ,b).
    !> Negative values are kept as-is (no clipping, no abs).
    subroutine mulliken_weights(n, nband, C, SC, groups, W)
        integer,                intent(in)  :: n, nband
        complex(dp),            intent(in)  :: C(n, nband), SC(n, nband)
        type(orbital_groups_t), intent(in)  :: groups
        real(dp),               intent(out) :: W(nband, groups%ngroup)

        integer :: g, b, i, mu
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

    !> Generalized Hermitian eigenproblem H C = S C diag(E) via
    !> zhegv(1,'V','U') on internal copies of H and S. Eigenvalues ascending,
    !> C^† S C = I. Never stops: ierr = info from zhegv (or the workspace query).
    subroutine generalized_eigensolve(n, H, S, E, C, ierr)
        integer,     intent(in)  :: n
        complex(dp), intent(in)  :: H(n, n), S(n, n)
        real(dp),    intent(out) :: E(n)
        complex(dp), intent(out) :: C(n, n)
        integer,     intent(out) :: ierr

        complex(dp), allocatable :: Bwork(:, :), work(:)
        complex(dp) :: wquery(1)
        real(dp),    allocatable :: rwork(:)
        integer :: lwork, info

        ierr = 0
        if (n <= 0) return

        C = H
        allocate(Bwork(n, n))
        Bwork = S
        allocate(rwork(max(1, 3 * n - 2)))

        ! Workspace query
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
