! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> hsiocell: supercélula diagonal de um .hsio (H e S) com eixos cortados e
!> átomos removidos (módulo hs_supercell). Usos:
!>  - tile lateral do volume (supercell = 11, 1, 1, tudo periódico), para o
!>    lead/substrato com a mesma célula lateral de uma supercélula DFT;
!>  - filme a partir do volume (supercell = N, n_cam, 1, periodic = T, F, T,
!>    vacuum_ang > 0), com a cadeia de cima de uma coluna removida
!>    (remove_cells), pronto para surfgreen mode = 'extracted' com o
!>    volume tile N x 1 x 1 como lead (embed_axis = 2, right = 'vacuum').
!>
!> Uso: hsiocell entrada.nml
!>
!>   &HSIOCELL
!>     h_file = 'TeBulk_H.hsio', s_file = 'TeBulk_S.hsio', structure_file = 'TeBulk.structure.nml'
!>     supercell = 11, 8, 1
!>     periodic = .true., .false., .true.
!>     vacuum_ang = 20.0                 ! somado ao comprimento de cada eixo cortado
!>     nremove_cells = 1
!>     remove_cells(:,1) = 5, 7, 0       ! cópia t = (t1, t2, t3), 0-based: todos os seus átomos
!>     nremove = 0                       ! remove_atoms = índices I(b,t) antes da remoção
!>     output_prefix = 'TeFilm_N11_L8'
!>   /
!> Índice antes da remoção: I(b, t) = b + nat*(t1 + n1*(t2 + n2*t3)).
!> Saídas: PREFIX_H.hsio, PREFIX_S.hsio (bloco-esparsos, meta de E_F copiada)
!> e PREFIX.structure.nml.
program hsiocell_main
    use, intrinsic :: iso_fortran_env, only: error_unit
    use types, only: dp
    use HS_io, only: HSBlockSparse, read_hs_blocksparse, write_hs_blocksparse
    use structure_io, only: structure_t, read_structure, write_structure
    use hs_supercell, only: hs_supercell_build, hs_supercell_index
    implicit none

    integer, parameter :: MAXR = 100000
    character(len=512) :: h_file = '', s_file = '', structure_file = '', output_prefix = 'supercell'
    integer :: supercell(3) = 1
    logical :: periodic(3) = .true.
    real(dp) :: vacuum_ang = 0.0_dp
    integer :: nremove = 0, nremove_cells = 0
    integer :: remove_atoms(MAXR) = 0, remove_cells(3, MAXR) = 0
    namelist /HSIOCELL/ h_file, s_file, structure_file, supercell, periodic, vacuum_ang, nremove, remove_atoms, &
        nremove_cells, remove_cells, output_prefix

    character(len=512) :: nml, emsg
    character(len=:), allocatable :: errmsg
    type(HSBlockSparse) :: hp, sp, hs, ss
    type(structure_t) :: st, sc
    integer :: u, ios, ierr, nat, ntot, k, b, t1, t2, t3, i, a
    integer, allocatable :: newidx(:), newidx_s(:)
    logical, allocatable :: keep(:)

    if (command_argument_count() < 1) then
        write(error_unit, '(a)') 'uso: hsiocell entrada.nml'
        stop 1
    end if
    call get_command_argument(1, nml)
    open(newunit=u, file=nml, status='old', action='read', iostat=ios)
    if (ios /= 0) call die('cannot open ' // trim(nml))
    read(u, nml=HSIOCELL, iostat=ios)
    close(u)
    if (ios /= 0) call die('error reading &HSIOCELL')
    if (nremove < 0 .or. nremove > MAXR .or. nremove_cells < 0 .or. nremove_cells > MAXR) call die('bad nremove')

    call read_hs_blocksparse(trim(h_file), hp, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    call read_hs_blocksparse(trim(s_file), sp, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    call read_structure(trim(structure_file), st, ierr, emsg)
    if (ierr /= 0) call die(trim(emsg))
    nat = hp%nAtoms
    if (st%nat /= nat .or. sp%nAtoms /= nat) call die('structure, H and S have different numbers of atoms')

    ntot = nat * product(supercell)
    allocate(keep(ntot))
    keep = .true.
    do k = 1, nremove
        if (remove_atoms(k) < 1 .or. remove_atoms(k) > ntot) call die('remove_atoms out of range')
        keep(remove_atoms(k)) = .false.
    end do
    do k = 1, nremove_cells
        if (any(remove_cells(:, k) < 0) .or. any(remove_cells(:, k) >= supercell)) call die('remove_cells out of range')
        do b = 1, nat
            keep(hs_supercell_index(nat, supercell, b, remove_cells(:, k))) = .false.
        end do
    end do

    call hs_supercell_build(hp, supercell, periodic, keep, hs, newidx, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    call hs_supercell_build(sp, supercell, periodic, keep, ss, newidx_s, ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    if (hs%Neighbors /= ss%Neighbors) call die('H and S give different cell lists (different block sets)')
    if (any(hs%IndexList /= ss%IndexList)) call die('H and S give different cell lists (different block sets)')

    ! estrutura
    sc%nat = count(keep)
    do a = 1, 3
        sc%lattice(a, :) = real(supercell(a), dp) * st%lattice(a, :)
        if (.not. periodic(a)) sc%lattice(a, :) = sc%lattice(a, :) + vacuum_ang * st%lattice(a, :) / norm2(st%lattice(a, :))
    end do
    allocate(sc%cart(3, sc%nat), sc%label(sc%nat))
    do t3 = 0, supercell(3) - 1
        do t2 = 0, supercell(2) - 1
            do t1 = 0, supercell(1) - 1
                do b = 1, nat
                    i = hs_supercell_index(nat, supercell, b, [t1, t2, t3])
                    if (newidx(i) == 0) cycle
                    sc%cart(:, newidx(i)) = st%cart(:, b) + real(t1, dp) * st%lattice(1, :) + &
                        real(t2, dp) * st%lattice(2, :) + real(t3, dp) * st%lattice(3, :)
                    sc%label(newidx(i)) = st%label(b)
                end do
            end do
        end do
    end do
    sc%has_species = st%has_species
    if (st%has_species) then
        sc%nspecies = st%nspecies
        sc%species_name = st%species_name
        sc%species_norb = st%species_norb
    end if

    call write_hs_blocksparse(trim(output_prefix) // '_H.hsio', hs, 'H', ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    call write_hs_blocksparse(trim(output_prefix) // '_S.hsio', ss, 'S', ierr, errmsg)
    if (ierr /= 0) call die(errmsg)
    call write_structure(trim(output_prefix) // '.structure.nml', sc, ierr, emsg)
    if (ierr /= 0) call die(trim(emsg))
    print '(a,3(1x,i0),a,i0,a,i0,a,i0,a,i0,a)', 'hsiocell: supercell', supercell, ': ', sc%nat, ' atoms (', &
        ntot - sc%nat, ' removed), ', hs%Neighbors, ' cells, ', hs%nBlocks, ' blocks'
    print '(a)', '  written: ' // trim(output_prefix) // '_H.hsio, _S.hsio, .structure.nml'

contains

    subroutine die(msg)
        character(len=*), intent(in) :: msg
        write(error_unit, '(a)') 'hsiocell: ' // trim(msg)
        stop 1
    end subroutine die

end program hsiocell_main
