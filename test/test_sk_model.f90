! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Modelo SK (sk_model):
!>  - leitura dos JSON do ajuste em Python (te_skg_{pd,spd}.json);
!>  - slab TeLead com a camada de cima (3 átomos), periódico em a1 e c: os
!>    autovalores em 8 pontos k coincidem com os do sk_gen.py
!>    (test/data/sk/sk_slab_n1_bands.txt) a 1e-10 eV;
!>  - H(k) hermitiana; S = I; escrita e releitura do JSON.
program test_sk_model
    use types, only: dp
    use testing, only: check, check_close, test_summary
    use HS_io, only: HSBlockSparse
    use structure_io, only: structure_t, read_structure
    use sk_model, only: sk_model_t, read_sk_model_json, write_sk_model_json, sk_build_blocksparse, &
                        sk_identity_blocksparse, sk_hk_dense
    use zlinalg, only: herm_eigvals
    implicit none

    call json_read()
    call slab_bands('pd', 'test/data/sk/te_skg_pd.json')
    call slab_bands('spd', 'test/data/sk/te_skg_spd.json')
    call json_roundtrip()
    call test_summary('test_sk_model')

contains

    subroutine json_read()
        type(sk_model_t) :: m
        integer :: ierr
        character(len=:), allocatable :: errmsg
        call read_sk_model_json('test/data/sk/te_skg_pd.json', m, ierr, errmsg)
        call check(ierr == 0, 'read te_skg_pd.json: ' // errmsg)
        if (ierr /= 0) return
        call check(m%norb == 8 .and. m%nl == 2 .and. m%nch == 18, 'pd: 8 orbitals, 2 on-site, 18 channels')
        call check(size(m%shells) == 4 .and. m%nso_per_atom == 2, 'pd: 4 shells, 2 occupied bands per atom')
        call check_close(m%p(1), -5.190013321705335_dp, 1.0e-15_dp, 'pd: E_p')
        call check_close(m%p(20), -0.34812466106983236_dp, 1.0e-15_dp, 'pd: dd delta_2')
        call check_close(m%rc, 4.8_dp, 1.0e-15_dp, 'pd: rc')
        call check(all(m%ch(:, 11) == [1, 2, 1, 0]), 'pd: channel 11 = pd sigma shell 2')
        call read_sk_model_json('test/data/sk/te_skg_spd.json', m, ierr, errmsg)
        call check(ierr == 0 .and. m%norb == 9 .and. m%loff(0) == 0 .and. m%loff(1) == 1 .and. m%loff(2) == 4, &
                   'spd: offsets s, p, d = 0, 1, 4')
    end subroutine json_read

    subroutine slab_bands(basis, jfile)
        character(len=*), intent(in) :: basis, jfile
        type(sk_model_t) :: m
        type(structure_t) :: st
        type(HSBlockSparse) :: h, s
        integer :: ierr, u, nk, nb, ik, i, n, j
        character(len=:), allocatable :: errmsg
        character(len=256) :: emsg
        character(len=8) :: tag
        real(dp) :: k(3), dmax, herm
        real(dp), allocatable :: ref(:), w(:), cart(:, :)
        complex(dp), allocatable :: hk(:, :), sk(:, :)
        integer :: top(3)

        call read_sk_model_json(jfile, m, ierr, errmsg)
        call check(ierr == 0, 'read ' // jfile)
        if (ierr /= 0) return
        call read_structure('test/data/sk/TeLead.structure.nml', st, ierr, emsg)
        call check(ierr == 0, 'read TeLead.structure.nml')
        if (ierr /= 0) return
        ! os 3 átomos de maior y, na ordem original
        do i = 1, 3
            top(i) = maxloc(st%cart(2, :), dim=1, mask=[(all(top(1:i - 1) /= j), j = 1, st%nat)])
        end do
        call sort3(top)
        cart = st%cart(:, top)
        call sk_build_blocksparse(m, st%lattice, cart, [.true., .false., .true.], h, ierr, errmsg)
        call check(ierr == 0, basis // ': build slab N=1: ' // errmsg)
        if (ierr /= 0) return
        call check(all(h%IndexList(:, 1) == 0), basis // ': R = 0 is the first cell')
        call check(all(h%IndexList(2, :) == 0), basis // ': no cell along the vacuum axis')
        call sk_identity_blocksparse(h, s)

        open(newunit=u, file='test/data/sk/sk_slab_n1_bands.txt', status='old', action='read')
        dmax = 0.0_dp
        herm = 0.0_dp
        n = 0
        do
            read(u, *, iostat=ierr) tag, nk, nb
            if (ierr /= 0) exit
            allocate(ref(nb))
            do ik = 1, nk
                read(u, *) k
                read(u, *) ref
                if (trim(tag) /= basis) cycle
                n = n + 1
                call sk_hk_dense(h, k, hk)
                herm = max(herm, maxval(abs(hk - conjg(transpose(hk)))))
                allocate(w(size(hk, 1)))
                call herm_eigvals(hk, w, ierr, errmsg)
                call check(size(w) == nb, basis // ': number of bands')
                if (size(w) == nb) dmax = max(dmax, maxval(abs(w - ref)))
                if (ik == nk) then
                    call sk_hk_dense(s, k, sk)
                    do i = 1, size(sk, 1)
                        sk(i, i) = sk(i, i) - 1.0_dp
                    end do
                    call check(maxval(abs(sk)) == 0.0_dp, basis // ': S(k) = I')
                end if
                deallocate(w)
            end do
            deallocate(ref)
        end do
        close(u)
        call check(n == 8, basis // ': 8 reference k-points')
        call check_close(herm, 0.0_dp, 1.0e-13_dp, basis // ': H(k) Hermitian')
        call check_close(dmax, 0.0_dp, 1.0e-10_dp, basis // ': slab N=1 bands = sk_gen.py (1e-10 eV)')
    end subroutine slab_bands

    subroutine json_roundtrip()
        type(sk_model_t) :: m, m2
        integer :: ierr
        character(len=:), allocatable :: errmsg
        call read_sk_model_json('test/data/sk/te_skg_spd.json', m, ierr, errmsg)
        call write_sk_model_json('test_sk_model_tmp.json', m, ierr, errmsg, note='teste')
        call check(ierr == 0, 'write JSON')
        call read_sk_model_json('test_sk_model_tmp.json', m2, ierr, errmsg)
        call check(ierr == 0, 'reread JSON: ' // errmsg)
        if (ierr /= 0) return
        call check(m2%nch == m%nch .and. all(m2%ch == m%ch) .and. m2%basis == m%basis, 'JSON round trip: channels')
        call check_close(maxval(abs(m2%p - m%p)), 0.0_dp, 1.0e-15_dp, 'JSON round trip: parameters')
        open(newunit=ierr, file='test_sk_model_tmp.json', status='old')
        close(ierr, status='delete')
    end subroutine json_roundtrip

    pure subroutine sort3(a)
        integer, intent(inout) :: a(3)
        integer :: i, j, t
        do i = 1, 2
            do j = i + 1, 3
                if (a(j) < a(i)) then
                    t = a(i); a(i) = a(j); a(j) = t
                end if
            end do
        end do
    end subroutine sort3

end program test_sk_model
