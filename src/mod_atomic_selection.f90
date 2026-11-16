! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Auxiliares de seleção de átomos (docs/design/contracts.md §2 para
!> parse_index_list, §6 e §6.1 para as regras geométricas).
!>
!> Papel no pipeline: usadas pelo CLI atomselect para transformar pedidos do
!> usuário (listas de índices, fatias ao longo de uma normal, caixas, camadas
!> topo/base) nas listas de átomos guardadas como regiões em regions.nml
!> (ver atomic_regions).
!>
!> Convenções: índices de átomos começam em 1, na ordem do arquivo de
!> estrutura; posições são cartesianas em Å (structure_t%cart); as listas de
!> índices devolvidas são ascendentes. As rotinas da biblioteca nunca param o
!> programa: erros voltam como ierr /= 0 mais errmsg.
!>
!>  - listas numéricas de índices ("1-12", "15, 17, 22", "1-10, 15-20");
!>  - seleção geométrica: projeção sobre uma normal arbitrária, s_i = r_i . n_hat
!>    (Å cartesianos; a normal é dada por suas componentes CARTESIANAS e nunca
!>    se assume que seja z), seleção por faixa e por caixa cartesiana;
!>  - detecção de camadas ao longo da normal (ligação simples sobre s
!>    ordenado), camadas de referência manuais e regiões de superfície
!>    topo/base.
!>
!> Periodicidade: as coordenadas projetadas NÃO são reembrulhadas. A detecção de
!> camadas supõe que a fatia não é cortada pela borda da célula ao longo de
!> n_hat; `check_periodic_split` detecta o caso oposto para que os chamadores
!> avisem (nada é reembrulhado).
module atomic_selection
    use types,        only: dp
    use structure_io, only: structure_t
    implicit none
    private

    public :: parse_index_list
    public :: project_positions, select_by_range, select_by_box
    public :: detect_layers, check_periodic_split, assign_to_references
    public :: surface_layers

    real(dp),         parameter :: NORMAL_EPS = 1.0e-12_dp   ! |normal| abaixo disto é tratada como zero
    character(len=*), parameter :: DIGITS     = '0123456789'
    integer,          parameter :: MAX_DIGITS = 9            ! mantém os valores bem dentro do integer padrão

contains

    ! ==================================================================
    ! Listas de índices
    ! ==================================================================

    !> Interpreta uma lista, separada por vírgulas, de índices (base 1) e faixas
    !> inclusivas a-b (p. ex. "1-10, 15, 17-20"; brancos em volta dos itens são
    !> ignorados) numa `list` ascendente e ordenada. Os índices devem estar em
    !> 1..nmax e não podem se repetir (inclusive entre itens); duplicatas são
    !> reportadas, não removidas. Só são aceitos inteiros decimais não negativos
    !> (sem sinal). Nunca para o programa: em erro ierr = 1, errmsg diz o motivo
    !> e `list` tem tamanho 0.
    subroutine parse_index_list(text, nmax, list, ierr, errmsg)
        character(len=*),     intent(in)  :: text
        integer,              intent(in)  :: nmax
        integer, allocatable, intent(out) :: list(:)
        integer,              intent(out) :: ierr
        character(len=*),     intent(out) :: errmsg

        logical,                      allocatable :: seen(:)
        integer                                   :: pos, next, ia, ib, k, i, dash
        character(len=:),             allocatable :: item

        ierr = 0
        errmsg = ''
        allocate(list(0))

        if (nmax < 1) then
            call fail('index list: nmax must be >= 1')
            return
        end if
        if (len_trim(text) == 0) then
            call fail('index list: empty text')
            return
        end if

        ! Percorre o texto item a item (itens separados por vírgulas). `seen` marca
        ! os índices já tomados, o que detecta duplicatas e também produz o
        ! resultado ordenado no final.
        allocate(seen(nmax), source=.false.)
        pos = 1
        k = 0
        do
            next = index(text(pos:), ',')
            if (next == 0) then
                item = trim(adjustl(text(pos:)))
            else
                item = trim(adjustl(text(pos:pos+next-2)))
            end if
            k = k + 1

            if (len(item) == 0) then
                call fail('index list: empty item at position ' // itoa(k))
                return
            end if

            ! índice único "n" ou faixa "a-b"
            dash = index(item, '-')
            if (dash == 0) then
                call to_int(item, ia)
                if (ierr /= 0) return
                ib = ia
            else
                call to_int(trim(adjustl(item(1:dash-1))), ia)
                if (ierr /= 0) return
                call to_int(trim(adjustl(item(dash+1:))), ib)
                if (ierr /= 0) return
                if (ia > ib) then
                    call fail('index list: reversed range "' // item // '"')
                    return
                end if
            end if

            if (ia < 1 .or. ib > nmax) then
                if (ia < 1) then
                    call fail('index list: index ' // itoa(ia) // ' < 1')
                else
                    call fail('index list: index ' // itoa(ib) // ' > ' // itoa(nmax))
                end if
                return
            end if

            do i = ia, ib
                if (seen(i)) then
                    call fail('index list: duplicated index ' // itoa(i))
                    return
                end if
                seen(i) = .true.
            end do

            if (next == 0) exit
            pos = pos + next
            ! uma vírgula final não deixa nada depois dela
            if (pos > len(text)) then
                call fail('index list: empty item at position ' // itoa(k + 1))
                return
            end if
        end do

        list = pack([(i, i = 1, nmax)], seen)

    contains

        !> Define o estado de erro e reinicia `list` com tamanho 0.
        subroutine fail(msg)
            character(len=*), intent(in) :: msg
            ierr = 1
            errmsg = msg
            if (allocated(list)) deallocate(list)
            allocate(list(0))
        end subroutine fail

        !> Converte um token só de dígitos; senão chama fail (e assim ierr é definido).
        subroutine to_int(tok, val)
            character(len=*), intent(in)  :: tok
            integer,          intent(out) :: val
            integer                       :: ios
            val = 0
            if (len(tok) == 0 .or. len(tok) > MAX_DIGITS .or. verify(tok, DIGITS) /= 0) then
                call fail('index list: not a non-negative integer: "' // tok // '"')
                return
            end if
            read(tok, *, iostat=ios) val
            if (ios /= 0) call fail('index list: not an integer: "' // tok // '"')
        end subroutine to_int

    end subroutine parse_index_list

    ! ==================================================================
    ! Seleção geométrica
    ! ==================================================================

    !> Projeta cada átomo sobre a normal: s(i) = r_i . n_hat (Å, com sinal, medido
    !> a partir da origem cartesiana), n_hat = normal/|normal|. `normal` contém
    !> componentes CARTESIANAS (não fracionárias) e não precisa ter norma 1.
    !> s tem tamanho st%nat. Erro (ierr = 1) se |normal| < 1e-12; s fica zerado.
    subroutine project_positions(st, normal, s, ierr, errmsg)
        type(structure_t),     intent(in)  :: st
        real(dp),              intent(in)  :: normal(3)
        real(dp), allocatable, intent(out) :: s(:)
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg
        real(dp)                           :: nhat(3), nrm
        integer                            :: i

        ierr = 0
        errmsg = ''
        allocate(s(st%nat))
        s = 0.0_dp
        nrm = norm2(normal)
        if (nrm < NORMAL_EPS) then
            ierr = 1
            errmsg = 'projection: normal vector has zero length'
            return
        end if
        nhat = normal / nrm
        do i = 1, st%nat
            s(i) = dot_product(st%cart(:, i), nhat)
        end do
    end subroutine project_positions

    !> Átomos com smin <= s(i) <= smax (inclusivo, mesmas unidades de s), índices
    !> ascendentes, possivelmente vazia. Nunca falha.
    subroutine select_by_range(s, smin, smax, list)
        real(dp),             intent(in)  :: s(:)
        real(dp),             intent(in)  :: smin, smax
        integer, allocatable, intent(out) :: list(:)
        integer :: i
        list = pack([(i, i = 1, size(s))], s >= smin .and. s <= smax)
    end subroutine select_by_range

    !> Átomos com lo(k) <= r_i(k) <= hi(k) para k = x,y,z (Å cartesianos); use
    !> -huge/huge para limites abertos. Índices ascendentes (possivelmente vazia).
    subroutine select_by_box(st, lo, hi, list)
        type(structure_t),    intent(in)  :: st
        real(dp),             intent(in)  :: lo(3), hi(3)
        integer, allocatable, intent(out) :: list(:)
        logical,              allocatable :: inside(:)
        integer                           :: i
        allocate(inside(st%nat))
        do i = 1, st%nat
            inside(i) = all(st%cart(:, i) >= lo) .and. all(st%cart(:, i) <= hi)
        end do
        list = pack([(i, i = 1, st%nat)], inside)
    end subroutine select_by_box

    ! ==================================================================
    ! Camadas
    ! ==================================================================

    !> Detecção de camadas por ligação simples sobre coordenadas projetadas (Å):
    !> os átomos são ordenados por s e uma nova camada começa quando
    !> s(próximo) - s(anterior) > tol, de modo que uma camada pode ser mais larga
    !> que tol se seus átomos se encadearem. As camadas são numeradas 1..nlayer
    !> em s ASCENDENTE; layer_of(i) é a camada do átomo i e layer_center(l) = s
    !> médio da camada l. Exige tol > 0 (senão ierr = 1; as saídas ficam
    !> vazias/zeradas). s vazio dá nlayer = 0. s não é reembrulhado (ver notas do
    !> módulo).
    subroutine detect_layers(s, tol, layer_of, nlayer, layer_center, ierr, errmsg)
        real(dp),              intent(in)  :: s(:)
        real(dp),              intent(in)  :: tol
        integer,  allocatable, intent(out) :: layer_of(:)
        integer,               intent(out) :: nlayer
        real(dp), allocatable, intent(out) :: layer_center(:)
        integer,               intent(out) :: ierr
        character(len=*),      intent(out) :: errmsg
        integer,               allocatable :: idx(:), cnt(:)
        integer                            :: n, k, l

        ierr = 0
        errmsg = ''
        n = size(s)
        nlayer = 0
        allocate(layer_of(n))
        layer_of = 0
        allocate(layer_center(0))
        if (.not. (tol > 0.0_dp)) then
            ierr = 1
            errmsg = 'layers: tolerance must be > 0'
            return
        end if
        if (n == 0) return

        ! 1. Ordena os átomos por s e corta em camadas nos intervalos maiores que tol
        idx = sort_index(s)
        nlayer = 1
        layer_of(idx(1)) = 1
        do k = 2, n
            if (s(idx(k)) - s(idx(k - 1)) > tol) nlayer = nlayer + 1
            layer_of(idx(k)) = nlayer
        end do

        ! 2. Centros das camadas = s médio dos membros
        deallocate(layer_center)
        allocate(layer_center(nlayer), cnt(nlayer))
        layer_center = 0.0_dp
        cnt = 0
        do k = 1, n
            l = layer_of(k)
            layer_center(l) = layer_center(l) + s(k)
            cnt(l) = cnt(l) + 1
        end do
        layer_center = layer_center / real(cnt, dp)
    end subroutine detect_layers

    !> warn = .true. se a fatia parece contínua através de uma fronteira periódica
    !> ao longo de n_hat: para todo vetor de rede a_i com p = |a_i . n_hat| > tol,
    !> a imagem periódica do átomo mais baixo fica em min(s) + p; se
    !> min(s) + p - max(s) <= tol a imagem se junta ao topo da fatia (ou se
    !> intercala com ele) e a detecção de camadas em s sem reembrulhar não é
    !> confiável. Vetores de rede com p <= tol (no plano) levam as camadas sobre
    !> si mesmas dentro da tolerância e são ignorados. Normal nula ou s vazio ->
    !> .false.
    subroutine check_periodic_split(st, normal, s, tol, warn)
        type(structure_t), intent(in)  :: st
        real(dp),          intent(in)  :: normal(3)
        real(dp),          intent(in)  :: s(:)
        real(dp),          intent(in)  :: tol
        logical,           intent(out) :: warn
        real(dp)                       :: nhat(3), nrm, p, smin, smax
        integer                        :: i

        warn = .false.
        nrm = norm2(normal)
        if (nrm < NORMAL_EPS .or. size(s) == 0) return
        nhat = normal / nrm
        smin = minval(s)
        smax = maxval(s)
        do i = 1, 3
            p = abs(dot_product(st%lattice(i, :), nhat))
            if (p <= tol) cycle
            if (smin + p - smax <= tol) warn = .true.
        end do
    end subroutine check_periodic_split

    !> Atribui cada átomo ao centro de referência manual mais próximo (mesmas
    !> unidades de s, Å). As camadas são numeradas por centro ASCENDENTE
    !> (layer_of(i) = posição do centro escolhido), logo a numeração coincide com
    !> a de detect_layers; um centro sem átomos ainda ocupa sua posição. Erros
    !> (ierr = 1, layer_of = 0): tol <= 0, nenhum centro, centros duplicados ou
    !> átomo mais distante que tol de todos os centros (geometria ambígua nunca
    !> é forçada).
    subroutine assign_to_references(s, ref_centers, tol, layer_of, ierr, errmsg)
        real(dp),             intent(in)  :: s(:)
        real(dp),             intent(in)  :: ref_centers(:)
        real(dp),             intent(in)  :: tol
        integer, allocatable, intent(out) :: layer_of(:)
        integer,              intent(out) :: ierr
        character(len=*),     intent(out) :: errmsg
        integer,  allocatable :: idx(:)
        real(dp), allocatable :: c(:)
        real(dp)              :: d, dbest
        integer               :: i, k, kbest, nref
        character(len=32)     :: b1, b2

        ierr = 0
        errmsg = ''
        allocate(layer_of(size(s)))
        layer_of = 0
        nref = size(ref_centers)
        if (.not. (tol > 0.0_dp)) then
            ierr = 1
            errmsg = 'layers: tolerance must be > 0'
            return
        end if
        if (nref == 0) then
            ierr = 1
            errmsg = 'layers: no reference centers given'
            return
        end if
        idx = sort_index(ref_centers)
        c = ref_centers(idx)
        do k = 2, nref
            if (c(k) == c(k - 1)) then
                write(b1, '(F12.6)') c(k)
                ierr = 1
                errmsg = 'layers: duplicated reference center ' // trim(adjustl(b1))
                return
            end if
        end do

        ! Centro mais próximo de cada átomo, rejeitando átomos longe demais de todos
        do i = 1, size(s)
            kbest = 0
            dbest = huge(1.0_dp)
            do k = 1, nref
                d = abs(s(i) - c(k))
                if (d < dbest) then
                    dbest = d
                    kbest = k
                end if
            end do
            if (dbest > tol) then
                write(b1, '(F12.6)') s(i)
                write(b2, '(F12.6)') tol
                ierr = 1
                errmsg = 'layers: atom ' // itoa(i) // ' (projected ' // trim(adjustl(b1)) // &
                         ') is farther than tol = ' // trim(adjustl(b2)) // ' from every reference center'
                layer_of = 0
                return
            end if
            layer_of(i) = kbest
        end do
    end subroutine assign_to_references

    !> Regiões de superfície a partir de uma atribuição de camadas (numeradas por s
    !> ascendente): top_list = átomos das ntop camadas mais altas
    !> (nlayer-ntop+1..nlayer), bottom_list = átomos das nbottom camadas mais
    !> baixas (1..nbottom). ntop ou nbottom pode ser 0 (não pedido), mas não
    !> ambos; valores negativos e ntop + nbottom > nlayer são erros (ierr = 1,
    !> listas vazias), logo topo e base nunca se sobrepõem. Listas ordenadas
    !> ascendentemente; átomos com layer_of < 1 (não atribuídos) são ignorados.
    subroutine surface_layers(layer_of, nlayer, ntop, nbottom, top_list, bottom_list, ierr, errmsg)
        integer,              intent(in)  :: layer_of(:)
        integer,              intent(in)  :: nlayer, ntop, nbottom
        integer, allocatable, intent(out) :: top_list(:), bottom_list(:)
        integer,              intent(out) :: ierr
        character(len=*),     intent(out) :: errmsg
        integer :: i

        ierr = 0
        errmsg = ''
        allocate(top_list(0), bottom_list(0))
        if (ntop < 0 .or. nbottom < 0) then
            ierr = 1
            errmsg = 'surface: number of layers must be >= 0'
            return
        end if
        if (ntop == 0 .and. nbottom == 0) then
            ierr = 1
            errmsg = 'surface: request at least one top or bottom layer'
            return
        end if
        if (ntop + nbottom > nlayer) then
            ierr = 1
            errmsg = 'surface: top (' // itoa(ntop) // ') + bottom (' // itoa(nbottom) // &
                     ') layers exceed the ' // itoa(nlayer) // ' layers found'
            return
        end if
        top_list = pack([(i, i = 1, size(layer_of))], layer_of > nlayer - ntop .and. layer_of >= 1)
        bottom_list = pack([(i, i = 1, size(layer_of))], layer_of >= 1 .and. layer_of <= nbottom)
    end subroutine surface_layers

    ! ==================================================================
    ! AUXILIARES PRIVADAS
    ! ==================================================================

    !> Permutação de ordenação ascendente estável (merge sort de baixo para cima):
    !> idx tal que v(idx) é ascendente; empates mantêm a ordem original. Usada no
    !> lugar de um intrínseco porque o Fortran não tem um, e a estabilidade deixa
    !> os resultados determinísticos para coordenadas iguais.
    function sort_index(v) result(idx)
        real(dp), intent(in)  :: v(:)
        integer,  allocatable :: idx(:)
        integer,  allocatable :: tmp(:)
        integer               :: n, i, width, lo, mid, hi, a, b, k

        n = size(v)
        idx = [(i, i = 1, n)]
        allocate(tmp(n))
        width = 1
        do while (width < n)
            lo = 1
            do while (lo <= n)
                mid = min(lo + width - 1, n)
                hi = min(lo + 2 * width - 1, n)
                a = lo
                b = mid + 1
                k = lo
                do while (k <= hi)
                    if (a <= mid .and. (b > hi .or. v(idx(a)) <= v(idx(min(b, n))))) then
                        tmp(k) = idx(a)
                        a = a + 1
                    else
                        tmp(k) = idx(b)
                        b = b + 1
                    end if
                    k = k + 1
                end do
                lo = lo + 2 * width
            end do
            idx = tmp
            width = 2 * width
        end do
    end function sort_index

    !> Inteiro para string sem brancos (para mensagens de erro).
    pure function itoa(i) result(s)
        integer,          intent(in)  :: i
        character(len=:), allocatable :: s
        character(len=16)             :: buf
        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

end module atomic_selection
