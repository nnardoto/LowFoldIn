! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Discrete spectral data (band energies + projection weights per atomic
!> region along a k-path) and its intermediate on-disk formats
!> (docs/design/contracts.md §5).
!>
!> Pipeline position: the spectral engine produces a spectral_data_t and
!> write_spectral_file stores it; downstream tools (specmap) read it back
!> with read_spectral_file, so they need neither the Hamiltonian nor the
!> structure.
!>
!>   - binary stream file (write_spectral_file / read_spectral_file)
!>   - text export for debugging (write_spectral_text, never read back)
!>
!> Conventions
!>   - Units are fixed by the format: energies in eV, k distances in 1/Angstrom.
!>   - Array shapes follow Fortran order, band index first:
!>     energy(nband,nk,nspin), weight(nband,nk,nregion,nspin).
!>   - All indices are 1-based. Region 1 is always named 'total'.
!>   - nspin = 2 only for collinear spin (spin_mode = SPIN_COLLINEAR).
!>   - Binary files use native byte order; an endian canary detects a
!>     mismatch. Integers are int32 and reals real64 on disk, independent of
!>     the default kinds used in memory.
!>
!> Library routines never stop: errors are returned through ierr/errmsg
!> (ierr = 0 means success).
module spectral_data
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use types, only: dp
    implicit none
    private

    ! ==================================================================
    ! Enumerations
    ! ==================================================================
    integer, parameter, public :: SPIN_NONE = 0, SPIN_COLLINEAR = 1, SPIN_NONCOLLINEAR = 2
    integer, parameter, public :: METHOD_MULLIKEN = 1, METHOD_LOWDIN = 2   ! projection scheme
    integer, parameter, public :: EREF_UNKNOWN = 0, EREF_FERMI = 1, EREF_ABSOLUTE = 2  ! energy zero

    ! ==================================================================
    ! Data type
    ! ==================================================================

    !> In-memory image of a spectral data file. The allocatable arrays must
    !> be allocated with the shapes given below before writing; read
    !> allocates them. vertex_k/vertex_label are optional when nvertex = 0.
    !> fermi_energy_ev is meaningful only if fermi_known; energy_ref tells
    !> whether energies are already shifted by E_F (EREF_FERMI) or absolute.
    type, public :: spectral_data_t
        integer :: version_major = 1, version_minor = 0
        integer :: nk = 0, nband = 0, nspin = 1, nregion = 0, nvertex = 0
        integer :: spin_mode = SPIN_NONE, method = METHOD_MULLIKEN, energy_ref = EREF_UNKNOWN
        logical :: fermi_known = .false.
        real(dp) :: fermi_energy_ev = 0
        character(len=16)  :: fingerprint = ''
        character(len=512) :: source = ''          ! free text: input file names
        real(dp), allocatable :: kfrac(:,:)        ! (3,nk)
        real(dp), allocatable :: kdist(:)          ! (nk) 1/Angstrom
        integer,  allocatable :: vertex_k(:)       ! (nvertex)
        character(len=16), allocatable :: vertex_label(:)    ! (nvertex)
        character(len=64), allocatable :: region_name(:)     ! (nregion); region 1 = 'total'
        integer,  allocatable :: region_natoms(:)  ! (nregion)
        real(dp), allocatable :: energy(:,:,:)     ! (nband, nk, nspin) eV
        real(dp), allocatable :: weight(:,:,:,:)   ! (nband, nk, nregion, nspin)
    end type spectral_data_t

    public :: write_spectral_file, read_spectral_file, write_spectral_text, find_region

    ! ==================================================================
    ! Binary format constants (layout table in contracts.md §5)
    ! ==================================================================
    character(len=8), parameter :: SD_MAGIC   = "LFSPECT" // achar(0)
    character(len=8), parameter :: SD_TRAILER = "LFSPEND" // achar(0)
    integer(int32),   parameter :: SD_VER_MAJOR = 1_int32, SD_VER_MINOR = 0_int32
    integer(int32),   parameter :: SD_CANARY = 168496141_int32   ! 0x0A0B0C0D, byte-order probe

    !> Fixed header size in bytes (items 1..9 of the layout table):
    !> magic + versions + canary + type sizes + 5 dims + 4 flags + E_F +
    !> fingerprint + source.
    integer(int64), parameter :: SD_HEADER_BYTES = 8 + 2*4 + 4 + 3*4 + 5*4 + 4*4 + 8 + 16 + 512

contains

    ! ==================================================================
    ! Public API
    ! ==================================================================

    !> Index of the region whose trimmed name equals trim(name); 0 if absent.
    integer function find_region(sd, name) result(idx)
        type(spectral_data_t), intent(in) :: sd
        character(len=*),      intent(in) :: name
        integer :: i
        idx = 0
        if (.not. allocated(sd%region_name)) return
        do i = 1, size(sd%region_name)
            if (trim(adjustl(sd%region_name(i))) == trim(adjustl(name))) then
                idx = i
                return
            end if
        end do
    end function find_region

    !> Write sd to a binary stream file (layout of contracts.md §5).
    !> The file is overwritten if it exists. sd is validated first (shapes,
    !> enumerations, region 1 = 'total'); the same checks apply to the text
    !> writer.
    !> ierr: 0 ok; 1 invalid sd or cannot open file; 2 I/O error while
    !> writing; 3 error on close. errmsg describes the failure.
    subroutine write_spectral_file(filename, sd, ierr, errmsg)
        character(len=*),      intent(in)  :: filename
        type(spectral_data_t), intent(in)  :: sd
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        integer :: u, ios, i
        integer(int32) :: fk

        call validate(sd, ierr, errmsg)
        if (ierr /= 0) return

        open(newunit=u, file=filename, access='stream', form='unformatted', &
             status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            errmsg = "write_spectral_file: cannot open '" // trim(filename) // "' for writing"
            return
        end if

        fk = 0_int32
        if (sd%fermi_known) fk = 1_int32

        ! Records follow the layout table order. Any write error leaves ios
        ! non-zero, which short-circuits the remaining writes.
        ! Header (items 1..9)
        write(u, iostat=ios) SD_MAGIC, SD_VER_MAJOR, SD_VER_MINOR, SD_CANARY, &
            int(storage_size(0_int32)/8, int32), int(storage_size(0_int64)/8, int32), &
            int(storage_size(0.0_real64)/8, int32), &
            int(sd%nk, int32), int(sd%nband, int32), int(sd%nspin, int32), &
            int(sd%nregion, int32), int(sd%nvertex, int32), &
            int(sd%spin_mode, int32), int(sd%method, int32), int(sd%energy_ref, int32), fk, &
            real(sd%fermi_energy_ev, real64), sd%fingerprint, sd%source
        ! k-points (item 10), vertices (11-12), regions (13-14)
        if (ios == 0) write(u, iostat=ios) real(sd%kfrac, real64), real(sd%kdist, real64)
        if (ios == 0 .and. sd%nvertex > 0) then
            write(u, iostat=ios) int(sd%vertex_k, int32)
            do i = 1, sd%nvertex
                if (ios == 0) write(u, iostat=ios) sd%vertex_label(i)
            end do
        end if
        do i = 1, sd%nregion
            if (ios == 0) write(u, iostat=ios) sd%region_name(i)
        end do
        if (ios == 0) write(u, iostat=ios) int(sd%region_natoms, int32)
        ! Bulk arrays (15-16), then n_ext = 0 (no extension blocks in v1.0) and trailer
        if (ios == 0) write(u, iostat=ios) real(sd%energy, real64)
        if (ios == 0) write(u, iostat=ios) real(sd%weight, real64)
        if (ios == 0) write(u, iostat=ios) 0_int32, SD_TRAILER
        if (ios /= 0) then
            close(u, iostat=i)
            ierr = 2
            errmsg = "write_spectral_file: I/O error while writing '" // trim(filename) // "'"
            return
        end if
        close(u, iostat=ios)
        if (ios /= 0) then
            ierr = 3
            errmsg = "write_spectral_file: error closing '" // trim(filename) // "'"
            return
        end if
        ierr = 0
        errmsg = ''
    end subroutine write_spectral_file

    !> Read a binary spectral file, validating the full layout.
    !> Every array is allocated by this routine. Before each read the number
    !> of bytes still available is checked against what is needed, so a
    !> truncated or corrupted file yields an error rather than a crash.
    !> On any error sd is reset to its default (empty) state.
    !> ierr: 0 ok; 1 missing/unopenable file; 2 bad magic; 3 unsupported
    !> version, endianness or type sizes; 4 invalid dimensions/flags;
    !> 5 invalid vertex/region content; 6 bad extension block;
    !> 7 trailer missing or trailing data; 8 file truncated; 9 read I/O error.
    !> Only version_major = 1 is accepted; any minor version is read, and
    !> unknown extension blocks are skipped using their size field.
    subroutine read_spectral_file(filename, sd, ierr, errmsg)
        character(len=*),      intent(in)  :: filename
        type(spectral_data_t), intent(out) :: sd
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        integer :: u, ios, i
        logical :: exists
        integer(int64) :: fsize, pos, nbytes_ext
        character(len=8)  :: magic, trailer
        character(len=16) :: tag
        integer(int32) :: vmaj, vmin, canary, sz(3), dims(5), flags(4), n_ext
        integer(int32), allocatable :: itmp(:)
        real(real64) :: ef

        ierr = 0
        errmsg = ''

        inquire(file=filename, exist=exists, size=fsize)
        if (.not. exists) then
            call fail(1, "file '" // trim(filename) // "' does not exist")
            return
        end if
        if (fsize < 0) then
            call fail(1, "cannot determine size of '" // trim(filename) // "'")
            return
        end if

        open(newunit=u, file=filename, access='stream', form='unformatted', &
             status='old', action='read', iostat=ios)
        if (ios /= 0) then
            call fail(1, "cannot open '" // trim(filename) // "' for reading")
            return
        end if

        ! `pos` is the 1-based byte offset of the next unread byte. It is
        ! advanced by hand after each read so that have() can check the
        ! remaining size; the reads without pos= continue sequentially from
        ! the previous one.
        pos = 1
        ! ---- fixed header --------------------------------------------
        if (.not. have(8_int64, "magic")) return
        read(u, pos=pos, iostat=ios) magic
        if (.not. ok_read("magic")) return
        if (magic /= SD_MAGIC) then
            call fail_close(2, "bad magic: not a LowFoldIn spectral file")
            return
        end if
        pos = pos + 8

        if (.not. have(SD_HEADER_BYTES - 8, "header")) return
        read(u, iostat=ios) vmaj, vmin, canary
        if (.not. ok_read("version/canary")) return
        if (canary /= SD_CANARY) then
            call fail_close(3, "endian canary mismatch (file written with different byte order?)")
            return
        end if
        if (vmaj /= SD_VER_MAJOR) then
            call fail_close(3, "unsupported version_major (expected 1)")
            return
        end if
        read(u, iostat=ios) sz
        if (.not. ok_read("type sizes")) return
        if (sz(1) /= 4 .or. sz(2) /= 8 .or. sz(3) /= 8) then
            call fail_close(3, "unexpected on-disk type sizes (expected int32=4, int64=8, real64=8)")
            return
        end if
        read(u, iostat=ios) dims, flags, ef
        if (.not. ok_read("dimensions")) return
        read(u, iostat=ios) sd%fingerprint, sd%source
        if (.not. ok_read("fingerprint/source")) return
        pos = pos + (SD_HEADER_BYTES - 8)

        sd%version_major = int(vmaj)
        sd%version_minor = int(vmin)
        sd%nk = int(dims(1)); sd%nband = int(dims(2)); sd%nspin = int(dims(3))
        sd%nregion = int(dims(4)); sd%nvertex = int(dims(5))
        sd%spin_mode = int(flags(1)); sd%method = int(flags(2)); sd%energy_ref = int(flags(3))
        sd%fermi_energy_ev = real(ef, dp)

        if (sd%nk < 1 .or. sd%nband < 1 .or. sd%nregion < 1 .or. sd%nvertex < 0) then
            call fail_close(4, "invalid dimensions (need nk>=1, nband>=1, nregion>=1, nvertex>=0)")
            return
        end if
        if (sd%nspin /= 1 .and. sd%nspin /= 2) then
            call fail_close(4, "invalid nspin (must be 1 or 2)")
            return
        end if
        if (sd%spin_mode < SPIN_NONE .or. sd%spin_mode > SPIN_NONCOLLINEAR) then
            call fail_close(4, "invalid spin_mode")
            return
        end if
        if ((sd%nspin == 2) .neqv. (sd%spin_mode == SPIN_COLLINEAR)) then
            call fail_close(4, "inconsistent nspin and spin_mode (nspin=2 iff collinear)")
            return
        end if
        if (sd%method /= METHOD_MULLIKEN .and. sd%method /= METHOD_LOWDIN) then
            call fail_close(4, "invalid method")
            return
        end if
        if (sd%energy_ref < EREF_UNKNOWN .or. sd%energy_ref > EREF_ABSOLUTE) then
            call fail_close(4, "invalid energy_ref")
            return
        end if
        if (flags(4) /= 0 .and. flags(4) /= 1) then
            call fail_close(4, "invalid fermi_known flag (must be 0 or 1)")
            return
        end if
        sd%fermi_known = (flags(4) == 1)

        ! ---- k-points ------------------------------------------------
        ! 8 bytes * (3 kfrac + 1 kdist) per k-point
        if (.not. have(8_int64 * 4_int64 * int(sd%nk, int64), "kfrac/kdist")) return
        allocate(sd%kfrac(3, sd%nk), sd%kdist(sd%nk))
        read(u, iostat=ios) sd%kfrac, sd%kdist
        if (.not. ok_read("kfrac/kdist")) return
        pos = pos + 32_int64 * int(sd%nk, int64)

        ! ---- vertices ------------------------------------------------
        ! 20 bytes per vertex: int32 index + char(16) label
        if (.not. have(20_int64 * int(sd%nvertex, int64), "vertices")) return
        allocate(itmp(sd%nvertex), sd%vertex_k(sd%nvertex), sd%vertex_label(sd%nvertex))
        if (sd%nvertex > 0) then
            read(u, iostat=ios) itmp
            if (.not. ok_read("vertex_k")) return
            sd%vertex_k = int(itmp)
            do i = 1, sd%nvertex
                read(u, iostat=ios) sd%vertex_label(i)
                if (.not. ok_read("vertex_label")) return
            end do
        end if
        deallocate(itmp)
        pos = pos + 20_int64 * int(sd%nvertex, int64)
        if (sd%nvertex > 0) then
            if (any(sd%vertex_k < 1) .or. any(sd%vertex_k > sd%nk)) then
                call fail_close(5, "vertex_k out of range 1..nk")
                return
            end if
        end if

        ! ---- regions -------------------------------------------------
        ! 68 bytes per region: char(64) name + int32 natoms
        if (.not. have(68_int64 * int(sd%nregion, int64), "regions")) return
        allocate(sd%region_name(sd%nregion), sd%region_natoms(sd%nregion), itmp(sd%nregion))
        do i = 1, sd%nregion
            read(u, iostat=ios) sd%region_name(i)
            if (.not. ok_read("region_name")) return
        end do
        read(u, iostat=ios) itmp
        if (.not. ok_read("region_natoms")) return
        sd%region_natoms = int(itmp)
        deallocate(itmp)
        pos = pos + 68_int64 * int(sd%nregion, int64)
        if (trim(adjustl(sd%region_name(1))) /= 'total') then
            call fail_close(5, "region 1 must be named 'total'")
            return
        end if

        ! ---- energies and weights --------------------------------------
        nbytes_ext = 8_int64 * int(sd%nband, int64) * int(sd%nk, int64) * int(sd%nspin, int64)
        if (.not. have(nbytes_ext, "energy")) return
        allocate(sd%energy(sd%nband, sd%nk, sd%nspin))
        read(u, iostat=ios) sd%energy
        if (.not. ok_read("energy")) return
        pos = pos + nbytes_ext

        nbytes_ext = nbytes_ext * int(sd%nregion, int64)
        if (.not. have(nbytes_ext, "weight")) return
        allocate(sd%weight(sd%nband, sd%nk, sd%nregion, sd%nspin))
        read(u, iostat=ios) sd%weight
        if (.not. ok_read("weight")) return
        pos = pos + nbytes_ext

        ! ---- extension blocks ------------------------------------------
        ! Forward-compatibility hook: each block is {tag char(16), nbytes
        ! int64, payload}. Payloads are skipped, never interpreted, here.
        if (.not. have(4_int64, "n_ext")) return
        read(u, iostat=ios) n_ext
        if (.not. ok_read("n_ext")) return
        pos = pos + 4
        if (n_ext < 0) then
            call fail_close(6, "negative n_ext")
            return
        end if
        do i = 1, int(n_ext)
            if (.not. have(24_int64, "extension block header")) return
            read(u, pos=pos, iostat=ios) tag, nbytes_ext
            if (.not. ok_read("extension block header")) return
            pos = pos + 24
            if (nbytes_ext < 0) then
                call fail_close(6, "extension block '" // trim(tag) // "' has negative size")
                return
            end if
            if (.not. have(nbytes_ext, "extension block payload")) return
            pos = pos + nbytes_ext
        end do

        ! ---- trailer ---------------------------------------------------
        if (.not. have(8_int64, "trailer")) return
        read(u, pos=pos, iostat=ios) trailer
        if (.not. ok_read("trailer")) return
        if (trailer /= SD_TRAILER) then
            call fail_close(7, "trailer not found at expected position (corrupted file)")
            return
        end if
        pos = pos + 8
        if (pos /= fsize + 1) then
            call fail_close(7, "unexpected trailing data after trailer")
            return
        end if

        close(u, iostat=ios)

    contains

        !> True if at least n bytes remain from pos; otherwise closes the
        !> file and sets ierr = 8.
        logical function have(n, what)
            integer(int64),   intent(in) :: n
            character(len=*), intent(in) :: what
            have = (fsize - pos + 1 >= n)
            if (.not. have) call fail_close(8, "file truncated: not enough bytes for " // what)
        end function have

        !> True if the last read succeeded (ios == 0); otherwise closes the
        !> file and sets ierr = 9.
        logical function ok_read(what)
            character(len=*), intent(in) :: what
            ok_read = (ios == 0)
            if (.not. ok_read) call fail_close(9, "I/O error reading " // what)
        end function ok_read

        subroutine fail_close(code, msg)
            integer,          intent(in) :: code
            character(len=*), intent(in) :: msg
            integer :: ios2
            close(u, iostat=ios2)
            call fail(code, msg)
        end subroutine fail_close

        subroutine fail(code, msg)
            integer,          intent(in) :: code
            character(len=*), intent(in) :: msg
            ierr = code
            errmsg = "read_spectral_file: " // msg
            call reset_sd(sd)
        end subroutine fail

    end subroutine read_spectral_file

    !> Human-readable export (debugging only, never read back).
    !> Layout: '#' metadata header, a '#' line naming the columns, then one
    !> row per (spin, k, band):
    !>   ispin ik kdist k1 k2 k3 iband energy w(1..nregion)
    !> kdist in 1/Angstrom, k1..k3 fractional, energy in eV.
    !> Same validation and ierr convention as write_spectral_file
    !> (0 ok; 1 invalid sd or cannot open; 2 I/O error; 3 error closing).
    subroutine write_spectral_text(filename, sd, ierr, errmsg)
        character(len=*),      intent(in)  :: filename
        type(spectral_data_t), intent(in)  :: sd
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        integer :: u, ios, i, is, ik, ib
        character(len=:), allocatable :: line

        call validate(sd, ierr, errmsg)
        if (ierr /= 0) return

        open(newunit=u, file=filename, access='sequential', form='formatted', &
             status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            errmsg = "write_spectral_text: cannot open '" // trim(filename) // "' for writing"
            return
        end if

        write(u, '(A)', iostat=ios) "# LowFoldIn spectral data (text export, not read back)"
        if (ios == 0) write(u, '(A,I0,A,I0)', iostat=ios) "# version ", &
            sd%version_major, ".", sd%version_minor
        if (ios == 0) write(u, '(A,I0,A,I0,A,I0,A,I0,A,I0)', iostat=ios) &
            "# nk ", sd%nk, "  nband ", sd%nband, "  nspin ", sd%nspin, &
            "  nregion ", sd%nregion, "  nvertex ", sd%nvertex
        if (ios == 0) write(u, '(A,A)', iostat=ios) "# spin_mode ", spin_mode_name(sd%spin_mode)
        if (ios == 0) write(u, '(A,A)', iostat=ios) "# method ", method_name(sd%method)
        if (ios == 0) write(u, '(A,A)', iostat=ios) "# energy_reference ", eref_name(sd%energy_ref)
        if (ios == 0) write(u, '(A,L1,A,ES23.15)', iostat=ios) "# fermi_known ", sd%fermi_known, &
            "  fermi_energy_ev ", sd%fermi_energy_ev
        if (ios == 0) write(u, '(A,A)', iostat=ios) "# fingerprint ", trim(sd%fingerprint)
        if (ios == 0) write(u, '(A,A)', iostat=ios) "# source ", trim(sd%source)
        if (ios == 0) write(u, '(A)', iostat=ios) "# units energy eV, k_distance 1/Angstrom"
        do i = 1, sd%nvertex
            if (ios == 0) write(u, '(A,I0,A,A,A,I0,A,ES23.15)', iostat=ios) "# vertex ", i, &
                " label ", trim(sd%vertex_label(i)), " ik ", sd%vertex_k(i), &
                " kdist ", sd%kdist(sd%vertex_k(i))
        end do
        do i = 1, sd%nregion
            if (ios == 0) write(u, '(A,I0,A,A,A,I0)', iostat=ios) "# region ", i, &
                " name ", trim(sd%region_name(i)), " natoms ", sd%region_natoms(i)
        end do
        line = "# ispin ik kdist k1 k2 k3 iband energy"
        do i = 1, sd%nregion
            line = line // " " // trim(adjustl(sd%region_name(i)))
        end do
        if (ios == 0) write(u, '(A)', iostat=ios) line

        do is = 1, sd%nspin
            do ik = 1, sd%nk
                do ib = 1, sd%nband
                    if (ios /= 0) exit
                    write(u, '(I0,1X,I0,4(1X,ES23.15),1X,I0,1X,ES23.15,*(1X,ES23.15))', &
                        iostat=ios) is, ik, sd%kdist(ik), sd%kfrac(:, ik), ib, &
                        sd%energy(ib, ik, is), sd%weight(ib, ik, :, is)
                end do
            end do
        end do

        if (ios /= 0) then
            close(u, iostat=i)
            ierr = 2
            errmsg = "write_spectral_text: I/O error while writing '" // trim(filename) // "'"
            return
        end if
        close(u, iostat=ios)
        if (ios /= 0) then
            ierr = 3
            errmsg = "write_spectral_text: error closing '" // trim(filename) // "'"
            return
        end if
        ierr = 0
        errmsg = ''
    end subroutine write_spectral_text

    ! ==================================================================
    ! Internals
    ! ==================================================================

    !> Deallocate every array of sd and restore the default scalars, so a
    !> failed read leaves nothing allocated.
    subroutine reset_sd(sd)
        type(spectral_data_t), intent(inout) :: sd
        sd = spectral_data_t()
    end subroutine reset_sd

    !> Consistency checks shared by the writers: dimensions, enumerations,
    !> nspin/spin_mode agreement, array allocation and shapes, vertex range
    !> and region 1 = 'total'. Mirrors the checks read_spectral_file applies
    !> to file contents. Sets ierr = 1 and errmsg on the first failure.
    subroutine validate(sd, ierr, errmsg)
        type(spectral_data_t), intent(in)  :: sd
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        ierr = 1
        errmsg = ''
        if (sd%nk < 1 .or. sd%nband < 1 .or. sd%nregion < 1 .or. sd%nvertex < 0) then
            errmsg = "spectral_data: invalid dimensions (need nk>=1, nband>=1, nregion>=1, nvertex>=0)"
            return
        end if
        if (sd%nspin /= 1 .and. sd%nspin /= 2) then
            errmsg = "spectral_data: nspin must be 1 or 2"
            return
        end if
        if (sd%spin_mode < SPIN_NONE .or. sd%spin_mode > SPIN_NONCOLLINEAR) then
            errmsg = "spectral_data: invalid spin_mode"
            return
        end if
        if ((sd%nspin == 2) .neqv. (sd%spin_mode == SPIN_COLLINEAR)) then
            errmsg = "spectral_data: nspin=2 iff spin_mode=SPIN_COLLINEAR"
            return
        end if
        if (sd%method /= METHOD_MULLIKEN .and. sd%method /= METHOD_LOWDIN) then
            errmsg = "spectral_data: invalid method"
            return
        end if
        if (sd%energy_ref < EREF_UNKNOWN .or. sd%energy_ref > EREF_ABSOLUTE) then
            errmsg = "spectral_data: invalid energy_ref"
            return
        end if
        if (.not. (allocated(sd%kfrac) .and. allocated(sd%kdist) .and. &
                   allocated(sd%region_name) .and. allocated(sd%region_natoms) .and. &
                   allocated(sd%energy) .and. allocated(sd%weight))) then
            errmsg = "spectral_data: required arrays not allocated"
            return
        end if
        if (size(sd%kfrac, 1) /= 3 .or. size(sd%kfrac, 2) /= sd%nk .or. size(sd%kdist) /= sd%nk) then
            errmsg = "spectral_data: kfrac/kdist shape mismatch"
            return
        end if
        if (size(sd%region_name) /= sd%nregion .or. size(sd%region_natoms) /= sd%nregion) then
            errmsg = "spectral_data: region arrays shape mismatch"
            return
        end if
        if (any(shape(sd%energy) /= [sd%nband, sd%nk, sd%nspin])) then
            errmsg = "spectral_data: energy shape mismatch"
            return
        end if
        if (any(shape(sd%weight) /= [sd%nband, sd%nk, sd%nregion, sd%nspin])) then
            errmsg = "spectral_data: weight shape mismatch"
            return
        end if
        if (sd%nvertex > 0) then
            if (.not. (allocated(sd%vertex_k) .and. allocated(sd%vertex_label))) then
                errmsg = "spectral_data: vertex arrays not allocated"
                return
            end if
        end if
        if (allocated(sd%vertex_k)) then
            if (size(sd%vertex_k) /= sd%nvertex) then
                errmsg = "spectral_data: vertex_k shape mismatch"
                return
            end if
            if (any(sd%vertex_k < 1) .or. any(sd%vertex_k > sd%nk)) then
                errmsg = "spectral_data: vertex_k out of range 1..nk"
                return
            end if
        end if
        if (allocated(sd%vertex_label)) then
            if (size(sd%vertex_label) /= sd%nvertex) then
                errmsg = "spectral_data: vertex_label shape mismatch"
                return
            end if
        end if
        if (trim(adjustl(sd%region_name(1))) /= 'total') then
            errmsg = "spectral_data: region 1 must be named 'total'"
            return
        end if
        ierr = 0
    end subroutine validate

    ! Enumeration -> text name, used by the text export ("invalid" if unknown).
    pure function spin_mode_name(m) result(s)
        integer, intent(in) :: m
        character(len=:), allocatable :: s
        select case (m)
        case (SPIN_NONE);         s = "none"
        case (SPIN_COLLINEAR);    s = "collinear"
        case (SPIN_NONCOLLINEAR); s = "noncollinear"
        case default;             s = "invalid"
        end select
    end function spin_mode_name

    pure function method_name(m) result(s)
        integer, intent(in) :: m
        character(len=:), allocatable :: s
        select case (m)
        case (METHOD_MULLIKEN); s = "mulliken"
        case (METHOD_LOWDIN);   s = "lowdin"
        case default;           s = "invalid"
        end select
    end function method_name

    pure function eref_name(m) result(s)
        integer, intent(in) :: m
        character(len=:), allocatable :: s
        select case (m)
        case (EREF_UNKNOWN);  s = "unknown"
        case (EREF_FERMI);    s = "fermi"
        case (EREF_ABSOLUTE); s = "absolute"
        case default;         s = "invalid"
        end select
    end function eref_name

end module spectral_data
