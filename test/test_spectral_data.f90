! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Testes do módulo spectral_data (contracts.md §5).
program test_spectral_data
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use types,         only: dp
    use testing,       only: check, test_summary
    use spectral_data, only: spectral_data_t, write_spectral_file, read_spectral_file, &
                             write_spectral_text, find_region, SPIN_NONE, SPIN_COLLINEAR, SPIN_NONCOLLINEAR, &
                             METHOD_MULLIKEN, METHOD_LOWDIN, EREF_FERMI, EREF_ABSOLUTE
    implicit none

    character(len=*), parameter :: F1 = "build/test_spectral_tmp.bin"
    character(len=*), parameter :: F2 = "build/test_spectral_tmp2.bin"
    character(len=*), parameter :: FT = "build/test_spectral_tmp.txt"

    type(spectral_data_t)         :: a, b
    integer                       :: ierr, i, n, sizes(4)
    integer(int64)                :: fsize
    character(len=512)            :: errmsg
    character(len=:), allocatable :: buf

    ! ---- 1. ida e volta, colinear ----------------------------------------
    call make_sd(a, nspin=2, spin_mode=SPIN_COLLINEAR, nvertex=2)
    call write_spectral_file(F1, a, ierr, errmsg)
    call check(ierr == 0, "write collinear: " // trim(errmsg))
    call read_spectral_file(F1, b, ierr, errmsg)
    call check(ierr == 0, "read collinear: " // trim(errmsg))
    if (ierr == 0) call check(sd_equal(a, b), "round trip collinear identical")

    ! ---- 2. ida e volta, não colinear, nvertex = 0 ------------------------
    call make_sd(a, nspin=1, spin_mode=SPIN_NONCOLLINEAR, nvertex=0)
    call write_spectral_file(F2, a, ierr, errmsg)
    call check(ierr == 0, "write noncollinear: " // trim(errmsg))
    call read_spectral_file(F2, b, ierr, errmsg)
    call check(ierr == 0, "read noncollinear: " // trim(errmsg))
    if (ierr == 0) call check(sd_equal(a, b), "round trip noncollinear identical")

    ! ---- 3. arquivos truncados ----------------------------------------------
    call make_sd(a, nspin=2, spin_mode=SPIN_COLLINEAR, nvertex=2)
    call write_spectral_file(F1, a, ierr, errmsg)
    call slurp(F1, buf)
    n = len(buf)
    sizes = [4, 40, n/2, n - 1]
    do i = 1, size(sizes)
        call dump(F2, buf(1:sizes(i)))
        call read_spectral_file(F2, b, ierr, errmsg)
        call check(ierr /= 0, "truncated file rejected")
    end do
    ! Também alguns outros pontos de corte (região de fronteira de cada seção)
    do i = 8, n - 1, 37
        call dump(F2, buf(1:i))
        call read_spectral_file(F2, b, ierr, errmsg)
        call check(ierr /= 0, "truncated file rejected (scan)")
    end do

    ! ---- 4. corrupções --------------------------------------------------
    call dump(F2, buf)
    call patch_bytes(F2, 1_int64, "XXXXXXXX")
    call read_spectral_file(F2, b, ierr, errmsg)
    call check(ierr /= 0, "bad magic rejected")

    call dump(F2, buf)
    call patch_int32(F2, 9_int64, 2_int32)
    call read_spectral_file(F2, b, ierr, errmsg)
    call check(ierr /= 0, "wrong major version rejected")

    call dump(F2, buf)
    call patch_int32(F2, 13_int64, 1_int32)   ! versão menor: aceita
    call read_spectral_file(F2, b, ierr, errmsg)
    call check(ierr == 0, "other minor version accepted: " // trim(errmsg))
    if (ierr == 0) call check(b%version_minor == 1, "minor version read back")

    call dump(F2, buf // achar(0))
    call read_spectral_file(F2, b, ierr, errmsg)
    call check(ierr /= 0, "trailing byte rejected")
    call check(index(errmsg, "trailing") > 0, "trailing byte message")

    ! ---- 5. um bloco de extensão -----------------------------------------
    ! Substitui o final [n_ext=0][trailer] (12 bytes) por
    ! [n_ext=1][tag char16][nbytes=8 int64][8 bytes de carga][trailer].
    call dump(F2, buf(1:n - 12))
    fsize = int(n - 12, int64)
    call append_ext(F2, fsize + 1)
    call read_spectral_file(F2, b, ierr, errmsg)
    call check(ierr == 0, "extension block accepted: " // trim(errmsg))
    if (ierr == 0) call check(sd_equal(a, b), "data identical with extension block")

    ! ---- 6. o escritor rejeita sd inconsistente --------------------------------
    call make_sd(a, nspin=2, spin_mode=SPIN_NONE, nvertex=2)
    call write_spectral_file(F2, a, ierr, errmsg)
    call check(ierr /= 0, "writer rejects nspin=2 with spin_mode none")

    call make_sd(a, nspin=2, spin_mode=SPIN_COLLINEAR, nvertex=2)
    a%region_name(1) = 'all'
    call write_spectral_file(F2, a, ierr, errmsg)
    call check(ierr /= 0, "writer rejects region 1 not 'total'")

    call make_sd(a, nspin=2, spin_mode=SPIN_COLLINEAR, nvertex=2)
    a%vertex_k(2) = 6
    call write_spectral_file(F2, a, ierr, errmsg)
    call check(ierr /= 0, "writer rejects vertex_k out of range")

    call make_sd(a, nspin=2, spin_mode=SPIN_COLLINEAR, nvertex=2)
    deallocate(a%weight)
    allocate(a%weight(3, 5, 2, 2))
    a%weight = 0
    call write_spectral_file(F2, a, ierr, errmsg)
    call check(ierr /= 0, "writer rejects weight shape mismatch")

    ! ---- 7. find_region e exportação em texto -----------------------------------
    call make_sd(a, nspin=2, spin_mode=SPIN_COLLINEAR, nvertex=2)
    call check(find_region(a, 'total') == 1, "find_region total")
    call check(find_region(a, 'bottom') == 3, "find_region bottom")
    call check(find_region(a, 'top   ') == 2, "find_region trimmed")
    call check(find_region(a, 'middle') == 0, "find_region absent")

    call write_spectral_text(FT, a, ierr, errmsg)
    call check(ierr == 0, "write text: " // trim(errmsg))
    call check(count_data_lines(FT) == a%nspin * a%nk * a%nband, "text data line count")

    ! ---- limpeza ------------------------------------------------------------
    call remove_file(F1)
    call remove_file(F2)
    call remove_file(FT)

    call test_summary("test_spectral_data")

contains

    subroutine make_sd(sd, nspin, spin_mode, nvertex)
        type(spectral_data_t), intent(out) :: sd
        integer, intent(in)                :: nspin, spin_mode, nvertex
        integer                            :: ik, ib, ir, is

        sd%nk = 5; sd%nband = 3; sd%nspin = nspin; sd%nregion = 3; sd%nvertex = nvertex
        sd%spin_mode = spin_mode
        if (nspin == 2) then
            sd%method = METHOD_MULLIKEN
            sd%energy_ref = EREF_FERMI
            sd%fermi_known = .true.
            sd%fermi_energy_ev = -3.141592653589793_dp
        else
            sd%method = METHOD_LOWDIN
            sd%energy_ref = EREF_ABSOLUTE
            sd%fermi_known = .false.
            sd%fermi_energy_ev = 0.0_dp
        end if
        sd%fingerprint = '0123456789abcdef'
        sd%source = 'structure.xyz test_H.hsio test_S.hsio'

        allocate(sd%kfrac(3, sd%nk), sd%kdist(sd%nk))
        do ik = 1, sd%nk
            sd%kfrac(:, ik) = [0.1_dp * ik, -0.2_dp * ik / 3.0_dp, 1.0_dp / (ik + 6)]
            sd%kdist(ik) = 0.123456789_dp * (ik - 1)
        end do
        allocate(sd%vertex_k(nvertex), sd%vertex_label(nvertex))
        if (nvertex == 2) then
            sd%vertex_k = [1, 5]
            sd%vertex_label = [character(len=16) :: 'G', 'M']
        end if
        allocate(sd%region_name(3), sd%region_natoms(3))
        sd%region_name = [character(len=64) :: 'total', 'top', 'bottom']
        sd%region_natoms = [10, 4, 6]

        allocate(sd%energy(sd%nband, sd%nk, sd%nspin))
        allocate(sd%weight(sd%nband, sd%nk, sd%nregion, sd%nspin))
        do is = 1, sd%nspin
            do ik = 1, sd%nk
                do ib = 1, sd%nband
                    sd%energy(ib, ik, is) = -5.0_dp + ib * 1.7_dp + ik / 7.0_dp + is * 0.01_dp
                    do ir = 1, sd%nregion
                        sd%weight(ib, ik, ir, is) = sin(real(ib + 3*ik + 7*ir + 11*is, dp))
                    end do
                end do
            end do
        end do
    end subroutine make_sd

    logical function sd_equal(x, y)
        type(spectral_data_t), intent(in) :: x, y
        sd_equal = .false.
        if (x%version_major /= y%version_major .or. x%version_minor /= y%version_minor) return
        if (x%nk /= y%nk .or. x%nband /= y%nband .or. x%nspin /= y%nspin) return
        if (x%nregion /= y%nregion .or. x%nvertex /= y%nvertex) return
        if (x%spin_mode /= y%spin_mode .or. x%method /= y%method) return
        if (x%energy_ref /= y%energy_ref) return
        if (x%fermi_known .neqv. y%fermi_known) return
        if (x%fermi_energy_ev /= y%fermi_energy_ev) return
        if (x%fingerprint /= y%fingerprint .or. x%source /= y%source) return
        if (any(shape(x%kfrac) /= shape(y%kfrac))) return
        if (any(x%kfrac /= y%kfrac) .or. any(x%kdist /= y%kdist)) return
        if (size(x%vertex_k) /= size(y%vertex_k)) return
        if (any(x%vertex_k /= y%vertex_k)) return
        if (any(x%vertex_label /= y%vertex_label)) return
        if (any(x%region_name /= y%region_name)) return
        if (any(x%region_natoms /= y%region_natoms)) return
        if (any(shape(x%energy) /= shape(y%energy))) return
        if (any(x%energy /= y%energy)) return
        if (any(shape(x%weight) /= shape(y%weight))) return
        if (any(x%weight /= y%weight)) return
        sd_equal = .true.
    end function sd_equal

    subroutine slurp(fname, s)
        character(len=*), intent(in)               :: fname
        character(len=:), allocatable, intent(out) :: s
        integer                                    :: u
        integer(int64)                             :: sz
        inquire(file=fname, size=sz)
        allocate(character(len=sz) :: s)
        open(newunit=u, file=fname, access='stream', form='unformatted', status='old', action='read')
        read(u) s
        close(u)
    end subroutine slurp

    subroutine dump(fname, s)
        character(len=*), intent(in) :: fname, s
        integer                      :: u
        open(newunit=u, file=fname, access='stream', form='unformatted', status='replace', &
             action='write')
        write(u) s
        close(u)
    end subroutine dump

    subroutine patch_bytes(fname, pos, s)
        character(len=*), intent(in) :: fname, s
        integer(int64),   intent(in) :: pos
        integer                      :: u
        open(newunit=u, file=fname, access='stream', form='unformatted', status='old', &
             action='readwrite')
        write(u, pos=pos) s
        close(u)
    end subroutine patch_bytes

    subroutine patch_int32(fname, pos, v)
        character(len=*), intent(in) :: fname
        integer(int64),   intent(in) :: pos
        integer(int32),   intent(in) :: v
        integer                      :: u
        open(newunit=u, file=fname, access='stream', form='unformatted', status='old', &
             action='readwrite')
        write(u, pos=pos) v
        close(u)
    end subroutine patch_int32

    subroutine append_ext(fname, pos)
        character(len=*), intent(in) :: fname
        integer(int64),   intent(in) :: pos
        integer                      :: u
        character(len=16)            :: tag
        tag = 'TESTEXT'
        open(newunit=u, file=fname, access='stream', form='unformatted', status='old', &
             action='readwrite')
        write(u, pos=pos) 1_int32, tag, 8_int64, 2.5_real64, "LFSPEND" // achar(0)
        close(u)
    end subroutine append_ext

    integer function count_data_lines(fname) result(nl)
        character(len=*), intent(in) :: fname
        integer                      :: u, ios
        character(len=4096)          :: line
        nl = 0
        open(newunit=u, file=fname, status='old', action='read')
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            nl = nl + 1
        end do
        close(u)
    end function count_data_lines

    subroutine remove_file(fname)
        character(len=*), intent(in) :: fname
        integer                      :: u, ios
        logical                      :: ex
        inquire(file=fname, exist=ex)
        if (.not. ex) return
        open(newunit=u, file=fname, status='old', iostat=ios)
        if (ios == 0) close(u, status='delete')
    end subroutine remove_file

end program test_spectral_data
