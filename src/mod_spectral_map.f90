! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Expressões de mapa sobre regiões nomeadas e saída em texto dos mapas espectrais
!> (docs/design/contracts.md §7).
!>
!> Posição no pipeline: usado pela ferramenta specmap depois da etapa de alargamento
!> (spectral_broadening), que transforma pares discretos (energia, peso) em um
!> mapa de intensidade A(ne,nk); este módulo seleciona/combina regiões e escreve
!> as saídas prontas para o gnuplot. Não faz física por si só.
!>
!>   - parse_map_expression: 'top-bottom', 'a+b-c', '-top' -> índices de região
!>     e coeficientes +/-1 (uma região repetida tem os coeficientes somados)
!>   - map_file_tag:         'top-bottom' -> 'top_minus_bottom'
!>   - write_map_file:       mapa em texto pronto para pm3d (um bloco por k)
!>   - write_kticks_file:    include do gnuplot com xtics/xrange/EF/ELABEL
!>
!> Convenções: energias em eV, distância k em 1/Angstrom (como em
!> spectral_data); os índices de região são de base 1 em sd%region_name.
!> As rotinas da biblioteca nunca param: os erros são devolvidos via ierr/errmsg
!> (ierr = 0 significa sucesso).
module spectral_map
    use types, only: dp
    use spectral_data, only: spectral_data_t, find_region, EREF_FERMI, EREF_ABSOLUTE
    implicit none
    private

    public :: parse_map_expression, write_map_file, write_kticks_file, map_file_tag

contains

    ! ==================================================================
    ! Expressões de mapa
    ! ==================================================================

    !> Analisa uma expressão de mapa: nomes de região unidos por '+' ou '-'; um
    !> '-' inicial é permitido; espaços são ignorados. Devolve os índices de região
    !> distintos (na ordem da primeira aparição) e seus coeficientes (soma dos +/-1
    !> de cada ocorrência; uma região repetida é combinada, p. ex. 'top+top' -> 2).
    !> Erros: expressão vazia, termo vazio ('top+', '++', '+top'), região
    !> desconhecida (a mensagem lista os nomes disponíveis).
    !> Nota: uma região cujo nome contém '+' ou '-' não pode ser endereçada.
    !> Em caso de erro (ierr = 1), region_idx e coeff são devolvidos com tamanho 0.
    !> O chamador combina as intensidades por região como sum(coeff(i)*I(region_idx(i))).
    subroutine parse_map_expression(expr, sd, region_idx, coeff, ierr, errmsg)
        character(len=*),      intent(in)  :: expr
        type(spectral_data_t), intent(in)  :: sd
        integer,  allocatable, intent(out) :: region_idx(:)
        real(dp), allocatable, intent(out) :: coeff(:)
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        character(len=:), allocatable :: s, name
        integer                       :: i, n, nterm, idx, j
        real(dp)                      :: sgn
        integer,  allocatable         :: tidx(:)
        real(dp), allocatable         :: tcoef(:)
        logical                       :: found, at_end

        ierr = 0
        errmsg = ''
        allocate(region_idx(0), coeff(0))

        s = strip_spaces(expr)
        n = len(s)
        if (n == 0) then
            call fail("empty map expression")
            return
        end if

        ! Varredura única da esquerda para a direita. 'name' acumula o termo atual;
        ! '+'/'-' (ou o fim da string, i = n+1) o fecha, e o operador que o
        ! fechou define o sinal do PRÓXIMO termo. Existem no máximo n termos,
        ! daí o tamanho dos arrays auxiliares.
        allocate(tidx(n), tcoef(n))
        nterm = 0
        sgn = 1.0_dp
        name = ''
        i = 1
        ! Um '-' inicial nega o primeiro termo ('+' no início é um erro de
        ! termo vazio, capturado abaixo).
        if (s(1:1) == '-') then
            sgn = -1.0_dp
            i = 2
        end if

        do while (i <= n + 1)
            if (i > n) then
                at_end = .true.
            else
                at_end = (s(i:i) == '+' .or. s(i:i) == '-')
            end if
            if (at_end) then
                ! Fim de um termo.
                if (len(name) == 0) then
                    call fail("empty term in map expression '" // trim(expr) // "'")
                    return
                end if
                idx = find_region(sd, name)
                if (idx == 0) then
                    call fail("unknown region '" // name // "' in map expression '" // &
                              trim(expr) // "'; available: " // available_names(sd))
                    return
                end if
                ! Une regiões repetidas somando seus coeficientes com sinal
                ! (assim 'top-top' resulta em coeficiente 0, não em erro).
                found = .false.
                do j = 1, nterm
                    if (tidx(j) == idx) then
                        tcoef(j) = tcoef(j) + sgn
                        found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    nterm = nterm + 1
                    tidx(nterm) = idx
                    tcoef(nterm) = sgn
                end if
                ! Sinal do termo seguinte (nenhum resta no fim da string).
                if (i <= n) then
                    if (s(i:i) == '+') then
                        sgn = 1.0_dp
                    else
                        sgn = -1.0_dp
                    end if
                end if
                name = ''
            else
                name = name // s(i:i)
            end if
            i = i + 1
        end do

        region_idx = tidx(1:nterm)
        coeff = tcoef(1:nterm)

    contains

        subroutine fail(msg)
            character(len=*), intent(in) :: msg
            ierr = 1
            errmsg = msg
            if (allocated(region_idx)) deallocate(region_idx)
            if (allocated(coeff)) deallocate(coeff)
            allocate(region_idx(0), coeff(0))
        end subroutine fail

    end subroutine parse_map_expression

    !> Tag de nome de arquivo de uma expressão de mapa: espaços removidos, um '-'
    !> inicial vira 'minus_', depois '+' -> '_plus_' e '-' -> '_minus_'.
    !> Usada para montar '<prefix>.<tag>.map.dat'. A expressão não é
    !> validada aqui.
    function map_file_tag(expr) result(tag)
        character(len=*), intent(in)  :: expr
        character(len=:), allocatable :: tag
        character(len=:), allocatable :: s
        integer                       :: i, i0

        s = strip_spaces(expr)
        tag = ''
        i0 = 1
        if (len(s) > 0) then
            if (s(1:1) == '-') then
                tag = 'minus_'
                i0 = 2
            end if
        end if
        do i = i0, len(s)
            select case (s(i:i))
            case ('+'); tag = tag // '_plus_'
            case ('-'); tag = tag // '_minus_'
            case default; tag = tag // s(i:i)
            end select
        end do
    end function map_file_tag

    ! ==================================================================
    ! Arquivos de saída
    ! ==================================================================

    !> Escreve um mapa A(ne,nk) na grade de energia grid(ne) ao longo de sd%kdist(nk).
    !> grid em eV, kdist em 1/Angstrom; A é a intensidade (unidades arbitrárias).
    !> Cada linha de cabeçalho é escrita com prefixo '# ' (a menos que já
    !> comece com '#'). Depois, para cada k, ne linhas
    !> 'k_distance energy_eV intensity' seguidas de uma linha em branco, de modo que o
    !> arquivo contém nk*ne linhas de dados e exatamente nk linhas em branco (a última
    !> encerra o bloco final; inofensiva para o pm3d do gnuplot).
    !> ierr: 0 ok; 1 tamanhos incompatíveis ou sem kdist; 2 não abre; 3 erro de escrita/fechamento.
    subroutine write_map_file(filename, sd, A, grid, header_lines, ierr, errmsg)
        character(len=*),      intent(in)  :: filename
        type(spectral_data_t), intent(in)  :: sd
        real(dp),              intent(in)  :: A(:,:), grid(:)
        character(len=*),      intent(in)  :: header_lines(:)
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        integer                     :: u, ios, ios2, ik, ie, i
        character(len=*), parameter :: ROWFMT = '(ES17.9E3, 1X, ES17.9E3, 1X, ES18.10E3)'

        ierr = 0
        errmsg = ''
        if (size(A, 1) /= size(grid)) then
            ierr = 1
            errmsg = "write_map_file: size(A,1) /= size(grid)"
            return
        end if
        if (.not. allocated(sd%kdist)) then
            ierr = 1
            errmsg = "write_map_file: spectral data has no kdist"
            return
        end if
        if (size(A, 2) /= size(sd%kdist)) then
            ierr = 1
            errmsg = "write_map_file: size(A,2) /= nk"
            return
        end if

        open(newunit=u, file=filename, status='replace', action='write', &
             form='formatted', iostat=ios)
        if (ios /= 0) then
            ierr = 2
            errmsg = "cannot open map file for writing: " // trim(filename)
            return
        end if

        ! Cabeçalho: cada linha vira um comentário do gnuplot.
        do i = 1, size(header_lines)
            if (len_trim(header_lines(i)) > 0) then
                if (header_lines(i)(1:1) == '#') then
                    write(u, '(A)', iostat=ios) trim(header_lines(i))
                else
                    write(u, '(A)', iostat=ios) '# ' // trim(header_lines(i))
                end if
            else
                write(u, '(A)', iostat=ios) '#'
            end if
            if (ios /= 0) exit
        end do
        if (ios == 0) write(u, '(A)', iostat=ios) '# columns: k_distance energy_eV intensity'

        ! Dados: blocos por k (a energia varia mais rápido) separados por linhas
        ! em branco, o layout que o 'splot ... with pm3d' do gnuplot espera.
        do ik = 1, size(A, 2)
            if (ios /= 0) exit
            do ie = 1, size(grid)
                write(u, ROWFMT, iostat=ios) sd%kdist(ik), grid(ie), A(ie, ik)
                if (ios /= 0) exit
            end do
            if (ios == 0) write(u, '(A)', iostat=ios) ''
        end do

        if (ios /= 0) then
            ierr = 3
            errmsg = "error writing map file: " // trim(filename)
            close(u, iostat=ios2)
            return
        end if
        close(u, iostat=ios)
        if (ios /= 0) then
            ierr = 3
            errmsg = "error closing map file: " // trim(filename)
        end if
    end subroutine write_map_file

    !> Escreve um arquivo de include do gnuplot com:
    !>   set xrange [0:kmax]
    !>   set xtics ("label" kdist, ...)     ('G'/'Gamma' -> {/Symbol G})
    !>   linhas verticais tracejadas nos vértices interiores (cor VLINE_COLOR,
    !>   padrão "white"; um script pode definir VLINE_COLOR antes de carregar)
    !>   EF = 0 (referência fermi) | E_F (absoluta e conhecida) | NaN
    !>   ELABEL = "E - E_F (eV)" | "E (eV)" | "E (eV, reference unknown)"
    !> A escolha de EF/ELABEL segue sd%energy_ref, de modo que os gráficos podem
    !> marcar corretamente o nível de Fermi qualquer que seja o zero das energias.
    !> As marcas de vértice são omitidas se sd não tiver vértices.
    !> ierr: 0 ok; 1 sem pontos k; 2 não abre; 3 erro no fechamento.
    subroutine write_kticks_file(filename, sd, ierr, errmsg)
        character(len=*),      intent(in)  :: filename
        type(spectral_data_t), intent(in)  :: sd
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg

        integer                       :: u, ios, ios2, iv, ik
        logical                       :: first
        real(dp)                      :: kmax, kd
        character(len=:), allocatable :: line
        character(len=32)             :: num

        ierr = 0
        errmsg = ''
        if (.not. allocated(sd%kdist) .or. sd%nk < 1) then
            ierr = 1
            errmsg = "write_kticks_file: spectral data has no k-points"
            return
        end if
        kmax = sd%kdist(sd%nk)

        open(newunit=u, file=filename, status='replace', action='write', &
             form='formatted', iostat=ios)
        if (ios /= 0) then
            ierr = 2
            errmsg = "cannot open kticks file for writing: " // trim(filename)
            return
        end if

        if (ios == 0) write(u, '(A)', iostat=ios) '# gnuplot include written by specmap: k-path ticks, Fermi level, energy label'
        if (ios == 0) write(u, '(A)', iostat=ios) 'set xrange [0:' // fmt_real(kmax) // ']'

        if (sd%nvertex > 0 .and. allocated(sd%vertex_k) .and. allocated(sd%vertex_label)) then
            line = 'set xtics ('
            first = .true.
            do iv = 1, sd%nvertex
                ik = sd%vertex_k(iv)
                if (ik < 1 .or. ik > sd%nk) cycle
                if (.not. first) line = line // ', '
                first = .false.
                line = line // '"' // gp_label(sd%vertex_label(iv)) // '" ' // fmt_real(sd%kdist(ik))
            end do
            line = line // ')'
            if (ios == 0) write(u, '(A)', iostat=ios) line
            if (ios == 0) write(u, '(A)', iostat=ios) 'if (!exists("VLINE_COLOR")) VLINE_COLOR = "white"'
            do iv = 1, sd%nvertex
                ik = sd%vertex_k(iv)
                if (ik < 1 .or. ik > sd%nk) cycle
                kd = sd%kdist(ik)
                if (kd <= 0.0_dp .or. kd >= kmax) cycle
                num = fmt_real(kd)
                if (ios == 0) write(u, '(A)', iostat=ios) 'set arrow from ' // trim(num) // ', graph 0 to ' // trim(num) // &
                    ', graph 1 nohead dt 2 lw 1 lc rgb VLINE_COLOR front'
            end do
        end if

        ! Referência de energia: com EREF_FERMI as energias já são
        ! relativas a E_F, então o nível de Fermi fica em 0.
        if (sd%energy_ref == EREF_FERMI) then
            if (ios == 0) write(u, '(A)', iostat=ios) 'EF = 0.0'
            if (ios == 0) write(u, '(A)', iostat=ios) 'ELABEL = "E - E_F (eV)"'
        else if (sd%energy_ref == EREF_ABSOLUTE) then
            if (sd%fermi_known) then
                if (ios == 0) write(u, '(A)', iostat=ios) 'EF = ' // fmt_real(sd%fermi_energy_ev)
            else
                if (ios == 0) write(u, '(A)', iostat=ios) 'EF = NaN'
            end if
            if (ios == 0) write(u, '(A)', iostat=ios) 'ELABEL = "E (eV)"'
        else
            if (ios == 0) write(u, '(A)', iostat=ios) 'EF = NaN'
            if (ios == 0) write(u, '(A)', iostat=ios) 'ELABEL = "E (eV, reference unknown)"'
        end if

        if (ios /= 0) then
            ierr = 3
            errmsg = "error writing kticks file: " // trim(filename)
            close(u, iostat=ios2)
            return
        end if

        close(u, iostat=ios)
        if (ios /= 0) then
            ierr = 3
            errmsg = "error closing kticks file: " // trim(filename)
        end if
    end subroutine write_kticks_file

    ! ==================================================================
    ! Rotinas auxiliares privadas
    ! ==================================================================

    !> Remove espaços e tabulações.
    pure function strip_spaces(s) result(r)
        character(len=*), intent(in)  :: s
        character(len=:), allocatable :: r
        integer                       :: i
        r = ''
        do i = 1, len(s)
            if (s(i:i) /= ' ' .and. s(i:i) /= achar(9)) r = r // s(i:i)
        end do
    end function strip_spaces

    !> Nomes de região separados por vírgula, para as mensagens de erro "unknown region".
    function available_names(sd) result(r)
        type(spectral_data_t), intent(in) :: sd
        character(len=:), allocatable     :: r
        integer                           :: i
        r = ''
        if (.not. allocated(sd%region_name)) return
        do i = 1, size(sd%region_name)
            if (i > 1) r = r // ', '
            r = r // trim(sd%region_name(i))
        end do
    end function available_names

    !> Rótulo do gnuplot (enhanced text): 'G'/'Gamma' -> '{/Symbol G}';
    !> aspas duplas são removidas.
    function gp_label(lab) result(r)
        character(len=*), intent(in)  :: lab
        character(len=:), allocatable :: r
        character(len=:), allocatable :: t
        integer                       :: i
        t = trim(adjustl(lab))
        if (t == 'G' .or. t == 'g' .or. t == 'Gamma' .or. t == 'GAMMA' .or. t == 'gamma') then
            r = '{/Symbol G}'
            return
        end if
        r = ''
        do i = 1, len(t)
            if (t(i:i) /= '"') r = r // t(i:i)
        end do
    end function gp_label

    !> Notação científica compacta, sem preenchimento, para números do gnuplot.
    function fmt_real(x) result(r)
        real(dp), intent(in)          :: x
        character(len=:), allocatable :: r
        character(len=32)             :: buf
        write(buf, '(ES16.8E3)') x
        r = trim(adjustl(buf))
    end function fmt_real

end module spectral_map
