! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

program test_structure_io
    use types, only: dp
    use testing, only: check, check_close, test_summary
    use structure_io, only: structure_t, read_structure, structure_fingerprint, fnv1a64_hex, &
                            reciprocal_lattice, cart_to_frac, species_norb_of
    implicit none

    character(len=*), parameter :: TMP = 'build/test_struct_tmp.nml'
    character(len=*), parameter :: NL = achar(10)
    real(dp), parameter :: TWO_PI = 2.0_dp * acos(-1.0_dp)
    character(len=*), parameter :: HEX_LAT = &
        "&STRUCTURE" // NL // &
        "  LatticeVectors(1,:) = 3.0, 0.0, 0.0" // NL // &
        "  LatticeVectors(2,:) = -1.5, 2.598076211353316, 0.0" // NL // &
        "  LatticeVectors(3,:) = 0.0, 0.0, 20.0" // NL // &
        "/" // NL

    type(structure_t) :: st, st2
    integer :: ierr, i, j
    character(len=256) :: errmsg
    character(len=16) :: fp0
    real(dp) :: a(3,3), b(3,3), r(3), f(3)

    ! 1. FNV-1a 64 test vectors
    call check(fnv1a64_hex("") == "cbf29ce484222325", "fnv1a64('')")
    call check(fnv1a64_hex("a") == "af63dc4c8601ec8c", "fnv1a64('a')")
    call check(fnv1a64_hex("foobar") == "85944171f73967e8", "fnv1a64('foobar')")

    ! 2. Read a small non-orthogonal structure with comments and blank lines
    call write_file(HEX_LAT // NL // &
        "&SPECIES" // NL // "! species table" // NL // NL // "2" // NL // &
        "Mo  9" // NL // "# comment" // NL // "S   4" // NL // "/" // NL // NL // &
        "&ATOMS" // NL // "  3" // NL // &
        "Mo  0.0  0.0  10.0" // NL // &
        "! the two S" // NL // NL // &
        "S   1.5  0.8660254037844386  11.5   ! trailing comment" // NL // &
        "S   1.5  0.8660254037844386  8.5" // NL // "/" // NL)
    call read_structure(TMP, st, ierr, errmsg)
    call check(ierr == 0, "read hex structure: " // trim(errmsg))
    call check(st%nat == 3, "nat == 3")
    call check(st%has_species .and. st%nspecies == 2, "species table present")
    if (st%nat == 3) then
        call check(st%label(1) == "Mo" .and. st%label(2) == "S" .and. st%label(3) == "S", "labels")
        call check_close(st%cart(3, 1), 10.0_dp, 1.0e-14_dp, "cart z atom 1")
        call check_close(st%cart(1, 2), 1.5_dp, 1.0e-14_dp, "cart x atom 2")
        call check_close(st%cart(2, 3), 0.8660254037844386_dp, 1.0e-14_dp, "cart y atom 3")
        call check(species_norb_of(st, 1) == 9, "norb Mo")
        call check(species_norb_of(st, 3) == 4, "norb S")
        call check_close(st%lattice(2, 2), 2.598076211353316_dp, 1.0e-14_dp, "lattice a2y")
    end if

    ! 3. Reciprocal lattice: b_i . a_j = 2π δ_ij
    a = st%lattice
    b = reciprocal_lattice(st)
    do i = 1, 3
        do j = 1, 3
            call check_close(dot_product(b(i, :), a(j, :)), merge(TWO_PI, 0.0_dp, i == j), &
                1.0e-12_dp, "b_i . a_j")
        end do
    end do

    ! 4. cart_to_frac round trip
    r = [0.7_dp, -1.3_dp, 4.2_dp]
    f = cart_to_frac(st, r)
    do i = 1, 3
        call check_close(f(1) * a(1, i) + f(2) * a(2, i) + f(3) * a(3, i), r(i), 1.0e-12_dp, &
            "cart_to_frac round trip")
    end do
    r = a(2, :)
    f = cart_to_frac(st, r)
    call check(all(abs(f - [0.0_dp, 1.0_dp, 0.0_dp]) < 1.0e-12_dp), "cart_to_frac(a2) = (0,1,0)")

    ! 5. Fingerprint
    fp0 = structure_fingerprint(st)
    call check(len_trim(fp0) == 16 .and. verify(fp0, "0123456789abcdef") == 0, "fingerprint is 16 hex")
    st2 = st
    st2%cart(1, 2) = st2%cart(1, 2) + 1.0e-7_dp
    call check(structure_fingerprint(st2) == fp0, "fingerprint stable under 1e-7 shift")
    st2 = st
    st2%cart(3, 3) = st2%cart(3, 3) + 1.0e-3_dp
    call check(structure_fingerprint(st2) /= fp0, "fingerprint changes under 1e-3 shift")
    st2 = st
    st2%label(2) = "Se"
    call check(structure_fingerprint(st2) /= fp0, "fingerprint changes with label")
    st2 = st
    st2%cart(:, 1) = st%cart(:, 2)
    st2%label(1) = st%label(2)
    st2%cart(:, 2) = st%cart(:, 1)
    st2%label(2) = st%label(1)
    call check(structure_fingerprint(st2) /= fp0, "fingerprint changes with atom order")
    st2 = st
    st2%lattice(3, 3) = 30.0_dp
    st2%lattice(2, 1) = 0.0_dp
    call check(structure_fingerprint(st2) == fp0, "fingerprint independent of lattice")
    st2 = st
    st2%cart(1, 1) = -0.00001_dp
    fp0 = structure_fingerprint(st2)
    st2%cart(1, 1) = 0.00001_dp
    call check(structure_fingerprint(st2) == fp0, "fingerprint: -1e-5 and +1e-5 equal")

    ! 6. Error cases (no stop)
    call read_structure('build/does_not_exist_struct.nml', st2, ierr, errmsg)
    call check(ierr /= 0 .and. index(errmsg, 'does_not_exist_struct.nml') > 0, "missing file -> ierr")

    call write_file(HEX_LAT // "&SPECIES" // NL // "1" // NL // "Mo 9" // NL // "/" // NL)
    call read_structure(TMP, st2, ierr, errmsg)
    call check(ierr /= 0, "missing &ATOMS -> ierr")

    call write_file(HEX_LAT // "&SPECIES" // NL // "1" // NL // "Mo 9" // NL // "/" // NL // &
        "&ATOMS" // NL // "2" // NL // "Mo 0 0 0" // NL // "W 1 0 0" // NL // "/" // NL)
    call read_structure(TMP, st2, ierr, errmsg)
    call check(ierr /= 0 .and. index(errmsg, 'W') > 0, "label not in &SPECIES -> ierr")

    call write_file("&STRUCTURE" // NL // &
        "  LatticeVectors(1,:) = 1.0, 0.0, 0.0" // NL // &
        "  LatticeVectors(2,:) = 2.0, 0.0, 0.0" // NL // &
        "  LatticeVectors(3,:) = 0.0, 0.0, 1.0" // NL // "/" // NL // &
        "&ATOMS" // NL // "1" // NL // "X 0 0 0" // NL // "/" // NL)
    call read_structure(TMP, st2, ierr, errmsg)
    call check(ierr /= 0, "singular lattice -> ierr")

    call write_file(HEX_LAT // "&ATOMS" // NL // "3" // NL // "X 0 0 0" // NL // "X 1 0 0" // NL)
    call read_structure(TMP, st2, ierr, errmsg)
    call check(ierr /= 0, "fewer atom lines than nAtoms -> ierr")

    ! No &SPECIES block: accepted, species_norb_of = -1
    call write_file(HEX_LAT // "&ATOMS" // NL // "1" // NL // "X 0 0 0" // NL // "/" // NL)
    call read_structure(TMP, st2, ierr, errmsg)
    call check(ierr == 0 .and. st2%nat == 1 .and. .not. st2%has_species, "no &SPECIES accepted")
    call check(species_norb_of(st2, 1) == -1, "species_norb_of = -1 without table")

    call delete_file(TMP)

    ! 7. Real example
    call read_structure('examples/UnfoldDiatomic/sc.nml', st2, ierr, errmsg)
    call check(ierr == 0 .and. st2%nat == 6, "read examples/UnfoldDiatomic/sc.nml: " // trim(errmsg))
    if (st2%nat == 6) call check(species_norb_of(st2, 2) == 1, "UnfoldDiatomic norb(B) == 1")

    call test_summary("test_structure_io")

contains

    subroutine write_file(text)
        character(len=*), intent(in) :: text
        integer :: iu
        open(newunit=iu, file=TMP, status='replace', action='write', form='formatted')
        write(iu, '(A)', advance='no') text
        close(iu)
    end subroutine write_file

    subroutine delete_file(name)
        character(len=*), intent(in) :: name
        integer :: iu, ios
        open(newunit=iu, file=name, status='old', iostat=ios)
        if (ios == 0) close(iu, status='delete')
    end subroutine delete_file

end program test_structure_io
