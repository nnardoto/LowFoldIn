! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Leitura da Hamiltoniana real-space do Wannier90 (`seedname_hr.dat`).
!>
!> O formato armazena H(R) acompanhado do peso de degenerescencia de cada
!> vetor R. A representacao interna do LowFoldIn usa diretamente
!> H(k)=sum_R H(R) exp(i 2 pi k.R); por isso os blocos sao divididos pelo
!> peso durante a leitura e `Degenerecency` fica unitario. A base de Wannier
!> e ortonormal; `build_identity_overlap` constroi S(R)=I delta(R,0) com a
!> mesma lista de R, para interoperar com os motores LCAO existentes.
module wannier90_io
    use types, only: dp
    use hs_io, only: HSDataFormat
    implicit none
    private

    public :: read_wannier90_hr, build_identity_overlap

contains

    subroutine read_wannier90_hr(filename, hs, ierr, errmsg)
        character(len=*),    intent(in)  :: filename
        type(HSDataFormat),  intent(out) :: hs
        integer,             intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg

        integer :: u, ios, nwann, nrpts, ir, i, j, m, n
        integer :: r1, r2, r3
        integer, allocatable :: deg(:)
        real(dp) :: hre, him
        character(len=256) :: title

        call hs%Clear()
        ierr = 0
        errmsg = ''

        open(newunit=u, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            call fail(1, "cannot open '" // trim(filename) // "'")
            return
        end if

        read(u, '(A)', iostat=ios) title
        if (ios /= 0) then
            call fail_close(2, 'missing title line')
            return
        end if
        read(u, *, iostat=ios) nwann
        if (ios /= 0 .or. nwann < 1) then
            call fail_close(2, 'invalid number of Wannier functions')
            return
        end if
        read(u, *, iostat=ios) nrpts
        if (ios /= 0 .or. nrpts < 1) then
            call fail_close(2, 'invalid number of R vectors')
            return
        end if

        allocate(deg(nrpts))
        read(u, *, iostat=ios) deg
        if (ios /= 0 .or. any(deg < 1)) then
            call fail_close(2, 'invalid or truncated degeneracy list')
            return
        end if

        hs%Title = trim(title) // ' [Wannier90; degeneracy-normalized]'
        hs%Neighbors = nrpts
        hs%nAtoms = nwann
        hs%NOrbs = nwann
        allocate(hs%Degenerecency(nrpts), hs%Orbs(nwann))
        allocate(hs%IndexList(3, nrpts), hs%Data(nwann, nwann, nrpts))
        hs%Degenerecency = 1
        hs%Orbs = 1
        hs%Data = cmplx(0.0_dp, 0.0_dp, dp)

        do ir = 1, nrpts
            do j = 1, nwann
                do i = 1, nwann
                    read(u, *, iostat=ios) r1, r2, r3, m, n, hre, him
                    if (ios /= 0) then
                        call fail_close(2, 'truncated H(R) matrix data')
                        return
                    end if
                    if (m /= i .or. n /= j) then
                        call fail_close(2, 'unexpected orbital ordering in H(R)')
                        return
                    end if
                    if (i == 1 .and. j == 1) then
                        hs%IndexList(:, ir) = [r1, r2, r3]
                    else if (any(hs%IndexList(:, ir) /= [r1, r2, r3])) then
                        call fail_close(2, 'R vector changes inside a matrix block')
                        return
                    end if
                    hs%Data(i, j, ir) = cmplx(hre, him, dp) / real(deg(ir), dp)
                end do
            end do
        end do

        close(u, iostat=ios)
        if (ios /= 0) call fail(3, 'error while closing input file')

    contains

        subroutine fail(code, message)
            integer,          intent(in) :: code
            character(len=*), intent(in) :: message
            ierr = code
            errmsg = 'read_wannier90_hr: ' // message
            call hs%Clear()
        end subroutine fail

        subroutine fail_close(code, message)
            integer,          intent(in) :: code
            character(len=*), intent(in) :: message
            integer :: close_ios
            close(u, iostat=close_ios)
            call fail(code, message)
        end subroutine fail_close
    end subroutine read_wannier90_hr

    subroutine build_identity_overlap(hr, sr, ierr, errmsg)
        type(HSDataFormat), intent(in)  :: hr
        type(HSDataFormat), intent(out) :: sr
        integer,            intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg

        integer :: ir, izero

        call sr%Clear()
        ierr = 0
        errmsg = ''
        izero = 0
        do ir = 1, hr%Neighbors
            if (all(hr%IndexList(:, ir) == 0)) then
                if (izero /= 0) then
                    ierr = 1
                    errmsg = 'build_identity_overlap: duplicate R=(0,0,0) block'
                    return
                end if
                izero = ir
            end if
        end do
        if (izero == 0) then
            ierr = 1
            errmsg = 'build_identity_overlap: missing R=(0,0,0) block'
            return
        end if

        sr%Title = 'Identity overlap for ' // trim(hr%Title)
        sr%Neighbors = hr%Neighbors
        sr%nAtoms = hr%nAtoms
        sr%NOrbs = hr%NOrbs
        allocate(sr%Degenerecency(hr%Neighbors), sr%Orbs(hr%nAtoms))
        allocate(sr%IndexList(3, hr%Neighbors))
        allocate(sr%Data(hr%NOrbs, hr%NOrbs, hr%Neighbors))
        sr%Degenerecency = 1
        sr%Orbs = hr%Orbs
        sr%IndexList = hr%IndexList
        sr%Data = cmplx(0.0_dp, 0.0_dp, dp)
        do ir = 1, hr%NOrbs
            sr%Data(ir, ir, izero) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
    end subroutine build_identity_overlap

end module wannier90_io
