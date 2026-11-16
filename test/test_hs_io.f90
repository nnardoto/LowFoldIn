! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! Teste de ida e volta WriteFile -> ReadFile do formato .hsio (blocos esparsos),
! incluindo blocos átomo-átomo inteiramente nulos.
program test_hs_io
    use types, only: dp
    use HS_io, only: HSDataFormat
    implicit none

    type(HSDataFormat) :: A, B
    character(len=*), parameter :: fname = "test_hs_io_roundtrip.hsio"
    integer                     :: i, j, R, nfail, unit_del

    nfail = 0

    ! 2 átomos (2 e 3 orbitais), 3 células R
    A%Title     = "roundtrip"
    A%Neighbors = 3
    A%nAtoms    = 2
    A%NOrbs     = 5
    allocate(A%Orbs(2), A%IndexList(3, 3), A%Data(5, 5, 3), A%Degenerecency(3))
    A%Orbs          = [2, 3]
    A%IndexList     = reshape([0, 0, 0,  1, 0, 0,  -1, 0, 0], [3, 3])
    A%Degenerecency = [1, 1, 1]

    A%Data = (0.0_dp, 0.0_dp)
    do R = 1, 2                      ! R = 3 fica totalmente nulo
        do j = 1, 5
            do i = 1, 5
                ! Zera o bloco (átomo 2, átomo 2) em R = 2 para testar blocos nulos
                if (R == 2 .and. i > 2 .and. j > 2) cycle
                A%Data(i, j, R) = cmplx(real(100*R + 10*i + j, dp), real(i - j, dp) / 7.0_dp, kind=dp)
            end do
        end do
    end do

    call A%WriteFile(fname)
    call B%ReadFile(fname)

    if (B%Neighbors /= A%Neighbors .or. B%nAtoms /= A%nAtoms .or. B%NOrbs /= A%NOrbs) then
        print '(A)', "FALHA: Neighbors/nAtoms/NOrbs"; nfail = nfail + 1
    end if
    if (.not. allocated(B%Orbs) .or. .not. allocated(B%IndexList) .or. .not. allocated(B%Data)) then
        print '(A)', "FALHA: arrays nao alocados"; stop 1
    end if
    if (any(B%Orbs /= A%Orbs)) then
        print '(A)', "FALHA: Orbs"; nfail = nfail + 1
    end if
    if (any(B%IndexList /= A%IndexList)) then
        print '(A)', "FALHA: IndexList"; nfail = nfail + 1
    end if
    if (any(shape(B%Data) /= shape(A%Data))) then
        print '(A)', "FALHA: shape(Data)"; stop 1
    end if
    if (any(B%Data /= A%Data)) then
        print '(A)', "FALHA: Data"; nfail = nfail + 1
    end if

    open(newunit=unit_del, file=fname, status='old')
    close(unit_del, status='delete')

    if (nfail > 0) then
        print '(A,I0,A)', "test_hs_io: ", nfail, " falha(s)"
        stop 1
    end if
    print '(A)', "test_hs_io: OK"
end program test_hs_io
