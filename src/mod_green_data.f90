! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Formato .green.bin (v1): resultados do surfgreen.
!>
!> Container binário (stream, little/big-endian nativo com canário) no mesmo
!> estilo do .hsio v2:
!>
!>   [header] [dataset 1] ... [dataset n] [TOC]
!>
!>   header: magic char(8) = 'LFGREEN'//NUL, header_bytes, ver_major, ver_minor,
!>           canary (int32), toc_offset (int64), toc_count (int32)
!>   TOC:    por dataset, name char(64), kind, schema (int32), offset, nbytes (int64)
!>
!> Datasets (schema 1):
!>   meta        (kind 1) nlines int32; linhas char(256) 'chave = valor'
!>   axes        (kind 2) kmode, nk, ne, nq, n1, n2, nvertex (int32);
!>                        kfrac(3,nk), kcart(3,nk), kdist(nk), energy(ne), qz(nq) (real64);
!>                        vertex_k(nvertex) int32; label(nvertex) char(16)
!>   maps        (kind 3) nreg int32; name(nreg) char(64);
!>                        total(ne,nk), A(ne,nk,nreg) (real64)
!>   coherent    (kind 4) icoh(ne,nk,nq), iinc(ne,nk,nq) (real64)   [opcional]
!>   diagnostics (kind 5) niter(ne,nk) int32; residual(ne,nk), min_neg_im_gii(ne,nk),
!>                        smin(nk) (real64)
!>
!> Um leitor ignora datasets de kind desconhecido. read_green_file devolve
!> ierr/errmsg e não para a execução.
module green_data
    use, intrinsic :: iso_fortran_env, only: int64
    use types, only: dp
    implicit none
    private

    public :: green_data_t, write_green_file, read_green_file, meta_add, meta_get
    public :: GREEN_KPATH, GREEN_KGRID, GREEN_NAME_LEN, GREEN_META_LEN, GREEN_LABEL_LEN

    integer, parameter :: GREEN_KPATH = 1, GREEN_KGRID = 2
    integer, parameter :: GREEN_NAME_LEN = 64, GREEN_META_LEN = 256, GREEN_LABEL_LEN = 16

    character(len=8), parameter :: MAGIC = 'LFGREEN' // achar(0)
    integer, parameter :: VER_MAJOR = 1, VER_MINOR = 0
    integer, parameter :: CANARY = 168496141
    integer, parameter :: HEADER_BYTES = 8 + 4 + 4 + 4 + 4 + 8 + 4
    integer, parameter :: TOC_ENTRY_BYTES = 64 + 4 + 4 + 8 + 8
    integer, parameter :: KIND_META = 1, KIND_AXES = 2, KIND_MAPS = 3, KIND_COH = 4, KIND_DIAG = 5
    integer, parameter :: SCHEMA_V1 = 1

    type :: green_data_t
        ! meta
        integer :: nmeta = 0
        character(len=GREEN_META_LEN), allocatable :: meta(:)
        ! axes
        integer :: kmode = GREEN_KPATH
        integer :: nk = 0, ne = 0, nq = 0, n1 = 0, n2 = 0, nvertex = 0
        real(dp), allocatable :: kfrac(:, :), kcart(:, :), kdist(:), energy(:), qz(:)
        integer,  allocatable :: vertex_k(:)
        character(len=GREEN_LABEL_LEN), allocatable :: label(:)
        ! maps
        integer :: nreg = 0
        character(len=GREEN_NAME_LEN), allocatable :: region(:)
        real(dp), allocatable :: total(:, :)          !< (ne, nk)
        real(dp), allocatable :: amap(:, :, :)        !< (ne, nk, nreg)
        ! coherent
        logical :: has_coherent = .false.
        real(dp), allocatable :: icoh(:, :, :), iinc(:, :, :)   !< (ne, nk, nq)
        ! diagnostics
        integer,  allocatable :: niter(:, :)
        real(dp), allocatable :: residual(:, :), min_neg_im_gii(:, :), smin(:)
    end type green_data_t

contains

    !> Acrescenta uma linha 'chave = valor' aos metadados.
    subroutine meta_add(gd, key, value)
        type(green_data_t), intent(inout) :: gd
        character(len=*),   intent(in)    :: key, value
        character(len=GREEN_META_LEN), allocatable :: tmp(:)
        if (.not. allocated(gd%meta)) allocate(gd%meta(0))
        allocate(tmp(gd%nmeta + 1))
        tmp(1:gd%nmeta) = gd%meta(1:gd%nmeta)
        tmp(gd%nmeta + 1) = trim(key) // ' = ' // trim(value)
        call move_alloc(tmp, gd%meta)
        gd%nmeta = gd%nmeta + 1
    end subroutine meta_add

    !> Valor da chave (string vazia se ausente).
    function meta_get(gd, key) result(value)
        type(green_data_t), intent(in) :: gd
        character(len=*),   intent(in) :: key
        character(len=:), allocatable :: value
        integer :: i, p
        value = ''
        do i = 1, gd%nmeta
            p = index(gd%meta(i), ' = ')
            if (p > 0) then
                if (gd%meta(i)(1:p - 1) == trim(key)) then
                    value = trim(gd%meta(i)(p + 3:))
                    return
                end if
            end if
        end do
    end function meta_get

    subroutine write_green_file(filename, gd, ierr, errmsg)
        character(len=*),   intent(in)  :: filename
        type(green_data_t), intent(in)  :: gd
        integer,            intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: u, ios, ntoc, i
        integer(int64) :: off(5), nby(5), pos, toc_offset
        character(len=GREEN_NAME_LEN) :: names(5)
        integer :: kinds(5)

        ierr = 0
        errmsg = ''
        open(newunit=u, file=filename, status='replace', action='write', form='unformatted', &
             access='stream', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = 'cannot open green file for writing: ' // trim(filename); return
        end if
        write(u, pos=1, iostat=ios) MAGIC, HEADER_BYTES, VER_MAJOR, VER_MINOR, CANARY, 0_int64, 0
        ntoc = 0

        ! meta
        call begin_ds('meta', KIND_META)
        write(u, iostat=ios) gd%nmeta
        do i = 1, gd%nmeta
            if (ios == 0) write(u, iostat=ios) gd%meta(i)
        end do
        call end_ds()

        ! axes
        call begin_ds('axes', KIND_AXES)
        if (ios == 0) write(u, iostat=ios) gd%kmode, gd%nk, gd%ne, gd%nq, gd%n1, gd%n2, gd%nvertex
        if (ios == 0) write(u, iostat=ios) gd%kfrac, gd%kcart, gd%kdist, gd%energy
        if (ios == 0 .and. gd%nq > 0) write(u, iostat=ios) gd%qz
        if (ios == 0 .and. gd%nvertex > 0) write(u, iostat=ios) gd%vertex_k, gd%label
        call end_ds()

        ! maps
        call begin_ds('maps', KIND_MAPS)
        if (ios == 0) write(u, iostat=ios) gd%nreg
        if (ios == 0 .and. gd%nreg > 0) write(u, iostat=ios) gd%region
        if (ios == 0) write(u, iostat=ios) gd%total
        if (ios == 0 .and. gd%nreg > 0) write(u, iostat=ios) gd%amap
        call end_ds()

        ! coherent
        if (gd%has_coherent) then
            call begin_ds('coherent', KIND_COH)
            if (ios == 0) write(u, iostat=ios) gd%icoh, gd%iinc
            call end_ds()
        end if

        ! diagnostics
        call begin_ds('diagnostics', KIND_DIAG)
        if (ios == 0) write(u, iostat=ios) gd%niter, gd%residual, gd%min_neg_im_gii, gd%smin
        call end_ds()

        ! TOC
        if (ios == 0) inquire(u, pos=toc_offset)
        do i = 1, ntoc
            if (ios == 0) write(u, iostat=ios) names(i), kinds(i), SCHEMA_V1, off(i), nby(i)
        end do
        if (ios == 0) write(u, pos=1 + 8 + 4 + 4 + 4 + 4, iostat=ios) toc_offset, ntoc
        close(u)
        if (ios /= 0) then
            ierr = 2; errmsg = 'error writing green file: ' // trim(filename)
        end if

    contains

        subroutine begin_ds(name, kind)
            character(len=*), intent(in) :: name
            integer,          intent(in) :: kind
            ntoc = ntoc + 1
            names(ntoc) = name
            kinds(ntoc) = kind
            if (ios == 0) inquire(u, pos=off(ntoc))
        end subroutine begin_ds

        subroutine end_ds()
            if (ios == 0) inquire(u, pos=pos)
            nby(ntoc) = pos - off(ntoc)
        end subroutine end_ds

    end subroutine write_green_file

    subroutine read_green_file(filename, gd, ierr, errmsg)
        character(len=*),   intent(in)  :: filename
        type(green_data_t), intent(out) :: gd
        integer,            intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=8) :: magic_read
        integer :: u, ios, hb, vmaj, vmin, can, ntoc, i, kind, schema
        integer(int64) :: toc_offset, off, nby
        character(len=GREEN_NAME_LEN) :: name
        logical :: got_axes, got_maps, got_diag

        ierr = 0
        errmsg = ''
        open(newunit=u, file=filename, status='old', action='read', form='unformatted', &
             access='stream', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = 'cannot open green file: ' // trim(filename); return
        end if
        read(u, pos=1, iostat=ios) magic_read, hb, vmaj, vmin, can, toc_offset, ntoc
        if (ios /= 0 .or. magic_read /= MAGIC) then
            close(u); ierr = 2; errmsg = trim(filename) // ' is not a .green.bin file'; return
        end if
        if (can /= CANARY) then
            close(u); ierr = 2; errmsg = 'green file written with a different byte order'; return
        end if
        if (vmaj /= VER_MAJOR) then
            close(u); ierr = 2; errmsg = 'unsupported .green.bin major version'; return
        end if
        got_axes = .false.; got_maps = .false.; got_diag = .false.
        do i = 1, ntoc
            read(u, pos=toc_offset + int(i - 1, int64) * TOC_ENTRY_BYTES, iostat=ios) name, kind, schema, off, nby
            if (ios /= 0) exit
            if (schema /= SCHEMA_V1) cycle
            select case (kind)
            case (KIND_META)
                read(u, pos=off, iostat=ios) gd%nmeta
                if (ios == 0) then
                    allocate(gd%meta(gd%nmeta))
                    if (gd%nmeta > 0) read(u, iostat=ios) gd%meta
                end if
            case (KIND_AXES)
                read(u, pos=off, iostat=ios) gd%kmode, gd%nk, gd%ne, gd%nq, gd%n1, gd%n2, gd%nvertex
                if (ios /= 0) exit
                allocate(gd%kfrac(3, gd%nk), gd%kcart(3, gd%nk), gd%kdist(gd%nk), gd%energy(gd%ne))
                allocate(gd%qz(gd%nq), gd%vertex_k(gd%nvertex), gd%label(gd%nvertex))
                read(u, iostat=ios) gd%kfrac, gd%kcart, gd%kdist, gd%energy
                if (ios == 0 .and. gd%nq > 0) read(u, iostat=ios) gd%qz
                if (ios == 0 .and. gd%nvertex > 0) read(u, iostat=ios) gd%vertex_k, gd%label
                got_axes = (ios == 0)
            case (KIND_MAPS)
                if (.not. got_axes) cycle
                read(u, pos=off, iostat=ios) gd%nreg
                if (ios /= 0) exit
                allocate(gd%region(gd%nreg), gd%total(gd%ne, gd%nk), gd%amap(gd%ne, gd%nk, gd%nreg))
                if (gd%nreg > 0) read(u, iostat=ios) gd%region
                if (ios == 0) read(u, iostat=ios) gd%total
                if (ios == 0 .and. gd%nreg > 0) read(u, iostat=ios) gd%amap
                got_maps = (ios == 0)
            case (KIND_COH)
                if (.not. got_axes) cycle
                allocate(gd%icoh(gd%ne, gd%nk, gd%nq), gd%iinc(gd%ne, gd%nk, gd%nq))
                read(u, pos=off, iostat=ios) gd%icoh, gd%iinc
                gd%has_coherent = (ios == 0)
            case (KIND_DIAG)
                if (.not. got_axes) cycle
                allocate(gd%niter(gd%ne, gd%nk), gd%residual(gd%ne, gd%nk))
                allocate(gd%min_neg_im_gii(gd%ne, gd%nk), gd%smin(gd%nk))
                read(u, pos=off, iostat=ios) gd%niter, gd%residual, gd%min_neg_im_gii, gd%smin
                got_diag = (ios == 0)
            end select
            if (ios /= 0) exit
        end do
        close(u)
        if (ios /= 0) then
            ierr = 3; errmsg = 'truncated or unreadable green file: ' // trim(filename); return
        end if
        if (.not. (got_axes .and. got_maps .and. got_diag)) then
            ierr = 3; errmsg = 'green file is missing axes, maps or diagnostics: ' // trim(filename)
        end if
        if (.not. allocated(gd%meta)) allocate(gd%meta(0))
    end subroutine read_green_file

end module green_data
