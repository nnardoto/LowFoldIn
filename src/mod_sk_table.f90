! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Integrais de dois centros de Slater–Koster (s, p, d) e descrição da base
!> PAO de cada espécie.
!>
!> Convenção da tabela (Slater e Koster, Phys. Rev. 94, 1498 (1954)): para o
!> orbital a (momento angular la) no átomo i e b (lb) no átomo j, com
!> n = (l, m, n) o vetor unitário de i para j,
!>
!>   E_ab(n) = tabela SK com (V_sigma, V_pi, V_delta) do par,      la <= lb
!>   E_ab(n) = (-1)^(la+lb) [E_ba(n)]^T (mesmos V),                 la >  lb
!>
!> ou seja, o parâmetro de um par sempre se refere à tabela com o menor l
!> primeiro; a espécie de cada orbital diz em que átomo ele está. Com essa
!> convenção os parâmetros universais de Harrison (V_ss_sigma < 0,
!> V_sp_sigma > 0, V_pp_sigma > 0, V_pp_pi < 0) têm os sinais usuais.
!>
!> Ordem das componentes m (a do OpenMX, conferida no manual: análise de
!> população): s; px, py, pz; d3z^2-r^2, dx^2-y^2, dxy, dxz, dyz.
!>
!> Base de uma espécie: string do OpenMX como 's3p2d2f1' (número de funções
!> radiais por l). Orbitais de um átomo, na ordem do OpenMX:
!>   para l = 0, 1, 2, 3: para zeta = 1..n_l: para m = 1..2l+1.
!> Com spin não colinear o átomo tem 2x esse número, na ordem
!> [todos os orbitais spin up | todos spin down] (openmx2hsio, atom-major).
module sk_table
    use types, only: dp
    implicit none
    private

    public :: sk_block, sk_channel_block, shell_t, species_basis_t, parse_basis, basis_shell_index, l_letter

    real(dp), parameter :: S3 = 1.7320508075688772935274463415059_dp

    !> Uma camada radial (l, zeta) de uma espécie e o deslocamento do seu
    !> primeiro orbital dentro do átomo (0-based, parte espacial).
    type :: shell_t
        integer :: l = 0
        integer :: zeta = 1
        integer :: offset = 0
    end type shell_t

    type :: species_basis_t
        character(len=16) :: name = ''
        character(len=32) :: spec = ''
        integer :: nspatial = 0
        integer :: nshell = 0
        type(shell_t), allocatable :: shell(:)
    end type species_basis_t

contains

    pure character(len=1) function l_letter(l)
        integer, intent(in) :: l
        character(len=4), parameter :: letters = 'spdf'
        l_letter = letters(l + 1:l + 1)
    end function l_letter

    !> Lê 'nome spec' (ex.: 'Te s3p2d2' ou 'Te s3p2d2f1').
    subroutine parse_basis(text, b, ierr, errmsg)
        character(len=*),      intent(in)  :: text
        type(species_basis_t), intent(out) :: b
        integer,               intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=64) :: nm, spec
        integer :: ios, i, l, nz, nsh, z, off
        integer :: count(0:3)
        character :: c

        ierr = 0
        errmsg = ''
        read(text, *, iostat=ios) nm, spec
        if (ios /= 0) then
            ierr = 1; errmsg = "basis entry must be 'Species s3p2d2', got '" // trim(text) // "'"; return
        end if
        count = 0
        i = 1
        spec = adjustl(spec)
        do while (i <= len_trim(spec))
            c = spec(i:i)
            l = index('spdf', c) - 1
            if (l < 0) then
                ierr = 1; errmsg = "invalid basis string '" // trim(spec) // "' (expected e.g. s3p2d2)"; return
            end if
            i = i + 1
            nz = 0
            do while (i <= len_trim(spec))
                if (index('0123456789', spec(i:i)) == 0) exit
                nz = 10 * nz + (iachar(spec(i:i)) - iachar('0'))
                i = i + 1
            end do
            if (nz < 1 .or. count(l) /= 0) then
                ierr = 1; errmsg = "invalid basis string '" // trim(spec) // "' (expected e.g. s3p2d2)"; return
            end if
            count(l) = nz
        end do
        b%name = nm
        b%spec = spec
        nsh = sum(count)
        b%nshell = nsh
        allocate(b%shell(nsh))
        nsh = 0
        off = 0
        do l = 0, 3
            do z = 1, count(l)
                nsh = nsh + 1
                b%shell(nsh) = shell_t(l, z, off)
                off = off + 2 * l + 1
            end do
        end do
        b%nspatial = off
    end subroutine parse_basis

    !> Índice da camada (l, zeta) na base, 0 se não existe.
    pure integer function basis_shell_index(b, l, zeta)
        type(species_basis_t), intent(in) :: b
        integer, intent(in) :: l, zeta
        integer :: k
        basis_shell_index = 0
        do k = 1, b%nshell
            if (b%shell(k)%l == l .and. b%shell(k)%zeta == zeta) then
                basis_shell_index = k
                return
            end if
        end do
    end function basis_shell_index

    !> Bloco (2la+1) x (2lb+1) para o canal ch (0 = sigma, 1 = pi, 2 = delta)
    !> com valor unitário.
    pure function sk_channel_block(la, lb, dir, ch) result(e)
        integer,  intent(in) :: la, lb, ch
        real(dp), intent(in) :: dir(3)
        real(dp) :: e(2 * la + 1, 2 * lb + 1)
        real(dp) :: v(3)
        v = 0.0_dp
        v(ch + 1) = 1.0_dp
        e = sk_block(la, lb, dir, v)
    end function sk_channel_block

    !> Bloco de dois centros E_ab(n), la, lb em 0..2, v = (V_sigma, V_pi,
    !> V_delta). dir não precisa estar normalizado (é normalizado aqui).
    pure function sk_block(la, lb, dir, v) result(e)
        integer,  intent(in) :: la, lb
        real(dp), intent(in) :: dir(3), v(3)
        real(dp) :: e(2 * la + 1, 2 * lb + 1)
        real(dp) :: n(3)
        n = dir / norm2(dir)
        if (la <= lb) then
            e = sk_ordered(la, lb, n, v)
        else
            e = real((-1)**(la + lb), dp) * transpose(sk_ordered(lb, la, n, v))
        end if
    end function sk_block

    !> Tabela SK para la <= lb, na ordem de componentes do OpenMX.
    pure function sk_ordered(la, lb, n, v) result(e)
        integer,  intent(in) :: la, lb
        real(dp), intent(in) :: n(3), v(3)
        real(dp) :: e(2 * la + 1, 2 * lb + 1)
        ! ordem padrão das d nas fórmulas: xy, yz, zx, x2-y2, 3z2-r2
        ! ordem OpenMX: 3z2-r2, x2-y2, xy, xz, yz  -> padrão (5, 4, 1, 3, 2)
        integer, parameter :: dperm(5) = [5, 4, 1, 3, 2]
        real(dp) :: l, m, nn, l2, m2, n2, sg, pi, dl
        real(dp) :: sd(5), pd(3, 5), dd(5, 5)
        integer :: a, b

        l = n(1); m = n(2); nn = n(3)
        l2 = l * l; m2 = m * m; n2 = nn * nn
        sg = v(1); pi = v(2); dl = v(3)

        if (la == 0 .and. lb == 0) then
            e(1, 1) = sg
        else if (la == 0 .and. lb == 1) then
            e(1, :) = [l, m, nn] * sg
        else if (la == 1 .and. lb == 1) then
            e(1, 1) = l2 * sg + (1.0_dp - l2) * pi
            e(2, 2) = m2 * sg + (1.0_dp - m2) * pi
            e(3, 3) = n2 * sg + (1.0_dp - n2) * pi
            e(1, 2) = l * m * (sg - pi);  e(2, 1) = e(1, 2)
            e(1, 3) = l * nn * (sg - pi); e(3, 1) = e(1, 3)
            e(2, 3) = m * nn * (sg - pi); e(3, 2) = e(2, 3)
        else if (la == 0 .and. lb == 2) then
            sd(1) = S3 * l * m * sg
            sd(2) = S3 * m * nn * sg
            sd(3) = S3 * nn * l * sg
            sd(4) = 0.5_dp * S3 * (l2 - m2) * sg
            sd(5) = (n2 - 0.5_dp * (l2 + m2)) * sg
            do b = 1, 5
                e(1, b) = sd(dperm(b))
            end do
        else if (la == 1 .and. lb == 2) then
            ! linhas x, y, z; colunas xy, yz, zx, x2-y2, 3z2-r2
            pd(1, 1) = S3 * l2 * m * sg + m * (1.0_dp - 2.0_dp * l2) * pi
            pd(1, 2) = S3 * l * m * nn * sg - 2.0_dp * l * m * nn * pi
            pd(1, 3) = S3 * l2 * nn * sg + nn * (1.0_dp - 2.0_dp * l2) * pi
            pd(2, 1) = S3 * m2 * l * sg + l * (1.0_dp - 2.0_dp * m2) * pi
            pd(2, 2) = S3 * m2 * nn * sg + nn * (1.0_dp - 2.0_dp * m2) * pi
            pd(2, 3) = S3 * l * m * nn * sg - 2.0_dp * l * m * nn * pi
            pd(3, 1) = S3 * l * m * nn * sg - 2.0_dp * l * m * nn * pi
            pd(3, 2) = S3 * n2 * m * sg + m * (1.0_dp - 2.0_dp * n2) * pi
            pd(3, 3) = S3 * n2 * l * sg + l * (1.0_dp - 2.0_dp * n2) * pi
            pd(1, 4) = 0.5_dp * S3 * l * (l2 - m2) * sg + l * (1.0_dp - l2 + m2) * pi
            pd(2, 4) = 0.5_dp * S3 * m * (l2 - m2) * sg - m * (1.0_dp + l2 - m2) * pi
            pd(3, 4) = 0.5_dp * S3 * nn * (l2 - m2) * sg - nn * (l2 - m2) * pi
            pd(1, 5) = l * (n2 - 0.5_dp * (l2 + m2)) * sg - S3 * l * n2 * pi
            pd(2, 5) = m * (n2 - 0.5_dp * (l2 + m2)) * sg - S3 * m * n2 * pi
            pd(3, 5) = nn * (n2 - 0.5_dp * (l2 + m2)) * sg + S3 * nn * (l2 + m2) * pi
            do a = 1, 3
                do b = 1, 5
                    e(a, b) = pd(a, dperm(b))
                end do
            end do
        else
            ! d-d (simétrico)
            dd(1, 1) = 3.0_dp * l2 * m2 * sg + (l2 + m2 - 4.0_dp * l2 * m2) * pi + (n2 + l2 * m2) * dl
            dd(2, 2) = 3.0_dp * m2 * n2 * sg + (m2 + n2 - 4.0_dp * m2 * n2) * pi + (l2 + m2 * n2) * dl
            dd(3, 3) = 3.0_dp * n2 * l2 * sg + (n2 + l2 - 4.0_dp * n2 * l2) * pi + (m2 + n2 * l2) * dl
            dd(1, 2) = 3.0_dp * l * m2 * nn * sg + l * nn * (1.0_dp - 4.0_dp * m2) * pi + l * nn * (m2 - 1.0_dp) * dl
            dd(1, 3) = 3.0_dp * l2 * m * nn * sg + m * nn * (1.0_dp - 4.0_dp * l2) * pi + m * nn * (l2 - 1.0_dp) * dl
            dd(2, 3) = 3.0_dp * m * n2 * l * sg + m * l * (1.0_dp - 4.0_dp * n2) * pi + m * l * (n2 - 1.0_dp) * dl
            dd(1, 4) = 1.5_dp * l * m * (l2 - m2) * sg + 2.0_dp * l * m * (m2 - l2) * pi + 0.5_dp * l * m * (l2 - m2) * dl
            dd(2, 4) = 1.5_dp * m * nn * (l2 - m2) * sg - m * nn * (1.0_dp + 2.0_dp * (l2 - m2)) * pi + &
                       m * nn * (1.0_dp + 0.5_dp * (l2 - m2)) * dl
            dd(3, 4) = 1.5_dp * nn * l * (l2 - m2) * sg + nn * l * (1.0_dp - 2.0_dp * (l2 - m2)) * pi - &
                       nn * l * (1.0_dp - 0.5_dp * (l2 - m2)) * dl
            dd(1, 5) = S3 * l * m * (n2 - 0.5_dp * (l2 + m2)) * sg - 2.0_dp * S3 * l * m * n2 * pi + &
                       0.5_dp * S3 * l * m * (1.0_dp + n2) * dl
            dd(2, 5) = S3 * m * nn * (n2 - 0.5_dp * (l2 + m2)) * sg + S3 * m * nn * (l2 + m2 - n2) * pi - &
                       0.5_dp * S3 * m * nn * (l2 + m2) * dl
            dd(3, 5) = S3 * l * nn * (n2 - 0.5_dp * (l2 + m2)) * sg + S3 * l * nn * (l2 + m2 - n2) * pi - &
                       0.5_dp * S3 * l * nn * (l2 + m2) * dl
            dd(4, 4) = 0.75_dp * (l2 - m2)**2 * sg + (l2 + m2 - (l2 - m2)**2) * pi + (n2 + 0.25_dp * (l2 - m2)**2) * dl
            dd(4, 5) = 0.5_dp * S3 * (l2 - m2) * (n2 - 0.5_dp * (l2 + m2)) * sg + S3 * n2 * (m2 - l2) * pi + &
                       0.25_dp * S3 * (1.0_dp + n2) * (l2 - m2) * dl
            dd(5, 5) = (n2 - 0.5_dp * (l2 + m2))**2 * sg + 3.0_dp * n2 * (l2 + m2) * pi + 0.75_dp * (l2 + m2)**2 * dl
            do a = 1, 5
                do b = a + 1, 5
                    dd(b, a) = dd(a, b)
                end do
            end do
            do a = 1, 5
                do b = 1, 5
                    e(a, b) = dd(dperm(a), dperm(b))
                end do
            end do
        end if
    end function sk_ordered

end module sk_table
