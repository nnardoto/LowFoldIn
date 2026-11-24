! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Potencial on-site por átomo (escalar ou bloco), independente da geometria
!> de quem o aplica.
!>
!> Cada entrada é identificada por um rótulo inteiro `tag` (na trincheira: o
!> deslocamento de coluna em relação à coluna da trincheira) e pela posição
!> cartesiana de referência do átomo (Å; na trincheira, a posição no lead, sem
!> o deslocamento de coluna). `site` é só informativo (índice do átomo na
!> estrutura de origem). O bloco V (norb x norb, eV, real simétrico) entra
!> como H_ii <- H_ii + V na base ortogonal; um escalar v é gravado como v*I.
!> Quem aplica procura a entrada por (tag, posição), com as posições
!> comparadas módulo translações de rede dadas (site_potential_find). Átomos
!> sem entrada ficam sem potencial.
!>
!> Formato de texto (v1); linhas que começam com '#' são comentários:
!>   LFSITEPOT 1
!>   norb <n>
!>   nentry <m>
!>   entry <tag> <site> <x> <y> <z> scalar <v>
!>   entry <tag> <site> <x> <y> <z> block
!>   <n linhas com n números>
!>
!> Rotinas sem stop; erros por ierr/errmsg.
module site_potential
    use types, only: dp
    implicit none
    private

    public :: site_potential_t, read_site_potential, write_site_potential, site_potential_find
    public :: site_potential_add

    type :: site_potential_t
        integer :: norb = 0
        integer :: n = 0
        integer, allocatable :: tag(:), site(:)
        real(dp), allocatable :: pos(:, :)          !< (3, n), Å
        real(dp), allocatable :: v(:, :, :)         !< (norb, norb, n), eV
        logical, allocatable :: is_scalar(:)
    end type site_potential_t

contains

    !> Acrescenta uma entrada (bloco v ou escalar, se v tem forma (1,1)).
    subroutine site_potential_add(sp, tag, site, pos, v)
        type(site_potential_t), intent(inout) :: sp
        integer,  intent(in) :: tag, site
        real(dp), intent(in) :: pos(3), v(:, :)
        integer, allocatable :: ti(:), si(:)
        real(dp), allocatable :: pi(:, :), vi(:, :, :)
        logical, allocatable :: li(:)
        integer :: n, i
        n = sp%n
        if (.not. allocated(sp%tag)) then
            allocate(sp%tag(0), sp%site(0), sp%pos(3, 0), sp%v(sp%norb, sp%norb, 0), sp%is_scalar(0))
        end if
        allocate(ti(n + 1), si(n + 1), pi(3, n + 1), vi(sp%norb, sp%norb, n + 1), li(n + 1))
        ti(1:n) = sp%tag; si(1:n) = sp%site; pi(:, 1:n) = sp%pos; vi(:, :, 1:n) = sp%v; li(1:n) = sp%is_scalar
        ti(n + 1) = tag; si(n + 1) = site; pi(:, n + 1) = pos
        vi(:, :, n + 1) = 0.0_dp
        if (size(v, 1) == 1 .and. size(v, 2) == 1) then
            do i = 1, sp%norb
                vi(i, i, n + 1) = v(1, 1)
            end do
            li(n + 1) = .true.
        else
            vi(:, :, n + 1) = v
            li(n + 1) = .false.
        end if
        call move_alloc(ti, sp%tag); call move_alloc(si, sp%site); call move_alloc(pi, sp%pos)
        call move_alloc(vi, sp%v); call move_alloc(li, sp%is_scalar)
        sp%n = n + 1
    end subroutine site_potential_add

    subroutine read_site_potential(filename, sp, ierr, errmsg)
        character(len=*),       intent(in)  :: filename
        type(site_potential_t), intent(out) :: sp
        integer,                intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=1024) :: line
        character(len=16) :: key, kind
        integer :: u, ios, ver, nent, tag, site, i, k
        real(dp) :: pos(3), sv(1, 1)
        real(dp), allocatable :: b(:, :)

        ierr = 0
        errmsg = ''
        open(newunit=u, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = "cannot open site potential file '" // trim(filename) // "'"; return
        end if
        call next_line(ios)
        read(line, *, iostat=ios) key, ver
        if (ios /= 0 .or. key /= 'LFSITEPOT' .or. ver /= 1) then
            call fail("not a LFSITEPOT 1 file"); return
        end if
        call next_line(ios)
        read(line, *, iostat=ios) key, sp%norb
        if (ios /= 0 .or. key /= 'norb' .or. sp%norb < 1) then
            call fail("expected 'norb <n>'"); return
        end if
        call next_line(ios)
        read(line, *, iostat=ios) key, nent
        if (ios /= 0 .or. key /= 'nentry' .or. nent < 0) then
            call fail("expected 'nentry <m>'"); return
        end if
        allocate(b(sp%norb, sp%norb))
        allocate(sp%tag(0), sp%site(0), sp%pos(3, 0), sp%v(sp%norb, sp%norb, 0), sp%is_scalar(0))
        do k = 1, nent
            call next_line(ios)
            if (ios /= 0) then
                call fail('unexpected end of file'); return
            end if
            read(line, *, iostat=ios) key, tag, site, pos, kind
            if (ios /= 0 .or. key /= 'entry') then
                call fail("bad entry line: " // trim(line)); return
            end if
            if (kind == 'scalar') then
                read(line, *, iostat=ios) key, tag, site, pos, kind, sv(1, 1)
                if (ios /= 0) then
                    call fail("bad scalar entry: " // trim(line)); return
                end if
                call site_potential_add(sp, tag, site, pos, sv)
            else if (kind == 'block') then
                do i = 1, sp%norb
                    call next_line(ios)
                    if (ios == 0) read(line, *, iostat=ios) b(i, :)
                    if (ios /= 0) then
                        call fail('bad block row'); return
                    end if
                end do
                call site_potential_add(sp, tag, site, pos, b)
            else
                call fail("entry kind must be 'scalar' or 'block'"); return
            end if
        end do
        close(u)

    contains

        subroutine next_line(st)
            integer, intent(out) :: st
            do
                read(u, '(a)', iostat=st) line
                if (st /= 0) return
                line = adjustl(line)
                if (len_trim(line) == 0) cycle
                if (line(1:1) == '#' .or. line(1:1) == '!') cycle
                return
            end do
        end subroutine next_line

        subroutine fail(msg)
            character(len=*), intent(in) :: msg
            ierr = 1
            errmsg = 'site potential ' // trim(filename) // ': ' // msg
            close(u)
        end subroutine fail

    end subroutine read_site_potential

    subroutine write_site_potential(filename, sp, ierr, errmsg, comment)
        character(len=*),       intent(in)  :: filename
        type(site_potential_t), intent(in)  :: sp
        integer,                intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=*), optional, intent(in) :: comment(:)
        integer :: u, ios, k, i
        ierr = 0
        errmsg = ''
        open(newunit=u, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = "cannot write site potential file '" // trim(filename) // "'"; return
        end if
        write(u, '(a)') 'LFSITEPOT 1'
        if (present(comment)) then
            do i = 1, size(comment)
                write(u, '(a)') '# ' // trim(comment(i))
            end do
        end if
        write(u, '(a,i0)') 'norb ', sp%norb
        write(u, '(a,i0)') 'nentry ', sp%n
        do k = 1, sp%n
            if (sp%is_scalar(k)) then
                write(u, '(a,i0,1x,i0,3(1x,es24.16e3),a,es24.16e3)') 'entry ', sp%tag(k), sp%site(k), sp%pos(:, k), &
                    ' scalar ', sp%v(1, 1, k)
            else
                write(u, '(a,i0,1x,i0,3(1x,es24.16e3),a)') 'entry ', sp%tag(k), sp%site(k), sp%pos(:, k), ' block'
                do i = 1, sp%norb
                    write(u, '(*(1x,es24.16e3))') sp%v(i, :, k)
                end do
            end if
        end do
        close(u)
    end subroutine write_site_potential

    !> Índice da primeira entrada com rótulo `tag` e posição igual a `pos` módulo
    !> combinações inteiras (|n| <= 2) das translações `trans` (3, nt; pode ser
    !> nt = 0), a no máximo tol (Å). 0 se não há.
    function site_potential_find(sp, tag, pos, trans, tol, match_tag) result(idx)
        type(site_potential_t), intent(in) :: sp
        integer,  intent(in) :: tag
        real(dp), intent(in) :: pos(3), trans(:, :), tol
        logical,  intent(in), optional :: match_tag   !< .false.: ignora o rótulo (padrão .true.)
        integer :: idx
        logical :: usetag
        integer :: k, n1, n2, nt, m1, m2
        real(dp) :: d(3)
        idx = 0
        usetag = .true.
        if (present(match_tag)) usetag = match_tag
        nt = size(trans, 2)
        m1 = merge(2, 0, nt >= 1)
        m2 = merge(2, 0, nt >= 2)
        do k = 1, sp%n
            if (usetag .and. sp%tag(k) /= tag) cycle
            do n1 = -m1, m1
                do n2 = -m2, m2
                    d = sp%pos(:, k) - pos
                    if (nt >= 1) d = d + n1 * trans(:, 1)
                    if (nt >= 2) d = d + n2 * trans(:, 2)
                    if (norm2(d) <= tol) then
                        idx = k
                        return
                    end if
                end do
            end do
        end do
    end function site_potential_find

end module site_potential
