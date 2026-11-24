! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> skgen: Hamiltoniano Slater–Koster ortogonal de um cristal em .hsio (S = I).
!>
!> Lê os parâmetros de um JSON (sk_model), a geometria de um structure.nml
!> (openmx2hsio --write-structure) e grava PREFIX_H.hsio, PREFIX_S.hsio
!> (bloco-esparsos) e PREFIX.structure.nml, para que surfgreen, spectral e
!> sktrench consumam o SK como consomem a DFT.
!>
!> Uso: skgen entrada.nml
!>
!>   &SKGEN
!>     model_file      = 'te_skg_pd.json'
!>     structure_file  = 'TeLead.structure.nml'
!>     periodic        = .true., .false., .true.   ! eixos de rede periódicos
!>     keep_layers     = 5        ! 0 = todos os átomos
!>     atoms_per_layer = 3        ! mantém os keep_layers*atoms_per_layer átomos
!>     normal_axis     = 2        ! de maior coordenada cartesiana normal_axis (ordem original)
!>     output_prefix   = 'TeSK_l5'
!>   /
!>   &SKBANDS                     ! opcional: autovalores em pontos k (frações)
!>     nk = 2
!>     kpoints(:,1) = 0.0, 0.0, 0.0
!>     kpoints(:,2) = 0.5, 0.0, 0.5
!>   /
!> Saída opcional: PREFIX.bands.dat (ik, k1, k2, k3, autovalores em eV).
program skgen_main
    use, intrinsic :: iso_fortran_env, only: error_unit
    use types, only: dp
    use HS_io, only: HSBlockSparse, write_hs_blocksparse
    use structure_io, only: structure_t, read_structure, write_structure
    use sk_model, only: sk_model_t, read_sk_model_json, sk_build_blocksparse, sk_identity_blocksparse, sk_hk_dense
    use zlinalg, only: herm_eigvals
    implicit none

    integer, parameter :: MAXK = 4096
    character(len=512) :: model_file = '', structure_file = '', output_prefix = 'sk'
    logical :: periodic(3) = [.true., .true., .true.]
    integer :: keep_layers = 0, atoms_per_layer = 3, normal_axis = 2
    integer :: nk = 0
    real(dp) :: kpoints(3, MAXK) = 0.0_dp
    namelist /SKGEN/ model_file, structure_file, periodic, keep_layers, atoms_per_layer, normal_axis, output_prefix
    namelist /SKBANDS/ nk, kpoints

    character(len=512) :: nml, emsg
    character(len=:), allocatable :: errmsg
    type(sk_model_t) :: m
    type(structure_t) :: st, sk
    type(HSBlockSparse) :: h, s
    integer :: u, ios, ierr, nkeep, i, ik
    integer, allocatable :: order(:), keep(:)
    logical, allocatable :: mask(:)
    complex(dp), allocatable :: hk(:, :)
    real(dp), allocatable :: w(:)
    logical :: have_bands

    if (command_argument_count() < 1) then
        write(error_unit, '(a)') 'uso: skgen entrada.nml'
        stop 1
    end if
    call get_command_argument(1, nml)
    open(newunit=u, file=nml, status='old', action='read', iostat=ios)
    if (ios /= 0) call die('cannot open ' // trim(nml))
    read(u, nml=SKGEN, iostat=ios)
    if (ios /= 0) call die('error reading &SKGEN')
    rewind(u)
    read(u, nml=SKBANDS, iostat=ios)
    have_bands = (ios == 0 .and. nk > 0)
    close(u)
    if (nk > MAXK) call die('&SKBANDS: nk too large')

    call read_sk_model_json(trim(model_file), m, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    call read_structure(trim(structure_file), st, ierr, emsg)
    if (ierr /= 0) call die(trim(emsg))

    ! --- átomos mantidos (ordem original) ---
    allocate(mask(st%nat))
    mask = .true.
    if (keep_layers > 0) then
        nkeep = keep_layers * atoms_per_layer
        if (nkeep > st%nat .or. normal_axis < 1 .or. normal_axis > 3) call die('&SKGEN: invalid keep_layers/normal_axis')
        order = argsort_desc(st%cart(normal_axis, :))
        mask = .false.
        mask(order(1:nkeep)) = .true.
    end if
    keep = pack([(i, i = 1, st%nat)], mask)
    sk%nat = size(keep)
    sk%lattice = st%lattice
    sk%cart = st%cart(:, keep)
    sk%label = st%label(keep)
    sk%has_species = .true.
    sk%nspecies = 1
    sk%species_name = [st%label(keep(1))]
    sk%species_norb = [m%norb]
    do i = 2, sk%nat
        if (sk%label(i) /= sk%label(1)) call die('skgen: only one species is supported')
    end do

    call sk_build_blocksparse(m, sk%lattice, sk%cart, periodic, h, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    call sk_identity_blocksparse(h, s)
    call write_hs_blocksparse(trim(output_prefix) // '_H.hsio', h, 'H', ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    call write_hs_blocksparse(trim(output_prefix) // '_S.hsio', s, 'S', ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    call write_structure(trim(output_prefix) // '.structure.nml', sk, ierr, emsg)
    if (ierr /= 0) call die(trim(emsg))
    print '(a,i0,a,i0,a,i0,a,i0,a)', 'skgen: ', sk%nat, ' atoms (', m%norb, ' orbitals each), ', &
        h%Neighbors, ' cells, ', h%nBlocks, ' blocks'
    print '(a)', '  written: ' // trim(output_prefix) // '_H.hsio, _S.hsio (S = I), .structure.nml'
    print '(a)', '  kept atoms (indices in ' // trim(structure_file) // '):'
    print '(16(1x,i0))', keep

    if (have_bands) then
        open(newunit=u, file=trim(output_prefix) // '.bands.dat', status='replace', action='write')
        write(u, '(a)') '# skgen: ik k1 k2 k3 (fractions) then eigenvalues (eV), ascending'
        allocate(w(h%NOrbs))
        do ik = 1, nk
            call sk_hk_dense(h, kpoints(:, ik), hk)
            hk = (0.5_dp, 0.0_dp) * (hk + conjg(transpose(hk)))
            call herm_eigvals(hk, w, ierr, errmsg)
            if (ierr /= 0) call die(errmsg)
            write(u, '(i6,3(1x,es24.16e3))', advance='no') ik, kpoints(:, ik)
            write(u, '(*(1x,es24.16e3))') w
        end do
        close(u)
        print '(a)', '  written: ' // trim(output_prefix) // '.bands.dat'
    end if

contains

    subroutine die(msg)
        character(len=*), intent(in) :: msg
        write(error_unit, '(a)') 'skgen: ' // trim(msg)
        stop 1
    end subroutine die

    !> Índices em ordem decrescente de x (estável).
    function argsort_desc(x) result(idx)
        real(dp), intent(in) :: x(:)
        integer :: idx(size(x))
        integer :: i, j, t
        idx = [(i, i = 1, size(x))]
        do i = 2, size(x)
            t = idx(i)
            j = i - 1
            do while (j >= 1)
                if (x(idx(j)) >= x(t)) exit
                idx(j + 1) = idx(j)
                j = j - 1
            end do
            idx(j + 1) = t
        end do
    end function argsort_desc

end program skgen_main
