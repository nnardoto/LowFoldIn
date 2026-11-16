! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Caminho k linear por partes (e grade k 2D, build_kgrid) em coordenadas fracionárias recíprocas, com
!> distância acumulada em 1/Å cartesiano (métrica recíproca,
!> docs/design/contracts.md §3.1).
!>
!> Posição no pipeline: o motor espectral constrói o caminho uma única vez,
!> avalia o Hamiltoniano em cada ponto kfrac e copia os dados kfrac/kdist/vértices
!> para spectral_data_t para o gráfico.
!>
!> Convenções
!>   - kfrac são coordenadas fracionárias (múltiplos dos vetores da rede
!>     recíproca b_i); kdist é cartesiano, em 1/Å, com kdist(1) = 0.
!>   - recip(3,3) guarda b_i como LINHAS e já inclui o fator 2*pi.
!>   - os índices vertex_k são posições (base 1) na lista de pontos k.
!> As rotinas da biblioteca nunca param: os erros são devolvidos via ierr/errmsg.
module kpath
    use types, only: dp
    implicit none
    private

    integer, parameter, public :: KLABEL_LEN = 16

    !> Um caminho k gerado. Todos os arrays só são alocados em caso de sucesso.
    type, public :: kpath_t
        integer                                :: nk = 0, nvertex = 0
        real(dp), allocatable                  :: kfrac(:,:)      ! (3,nk)
        real(dp), allocatable                  :: kdist(:)        ! (nk) 1/Å, kdist(1)=0
        integer,  allocatable                  :: vertex_k(:)     ! (nvertex) índice do vértice na lista k
        character(len=KLABEL_LEN), allocatable :: label(:)        ! (nvertex)
    end type kpath_t

    !> Grade k 2D (plano gerado por span1 e span2) para mapas de energia constante.
    !> Ponto (i,j): kfrac = origin + (i-1)/(n1-1)*span1 + (j-1)/(n2-1)*span2,
    !> i = 1..n1, j = 1..n2, ambos os extremos incluídos. Índice plano
    !> ik = i + (j-1)*n1 (i varia mais rápido).
    !> kcart = kfrac·recip (1/Å). kuv são as coordenadas no plano da grade ao
    !> longo da base ortonormal (e_u, e_v): e_u ∥ span1 cartesiano e e_v é a
    !> parte de span2 ortogonal a e_u (Gram-Schmidt). Γ projeta em (0,0).
    type, public :: kgrid_t
        integer               :: n1 = 0, n2 = 0, nk = 0
        real(dp)              :: origin(3) = 0.0_dp, span1(3) = 0.0_dp, span2(3) = 0.0_dp
        real(dp)              :: e_u(3) = 0.0_dp, e_v(3) = 0.0_dp      ! base no plano (cartesiana)
        real(dp), allocatable :: kfrac(:,:)   ! (3,nk)
        real(dp), allocatable :: kcart(:,:)   ! (3,nk) 1/Å
        real(dp), allocatable :: kuv(:,:)     ! (2,nk) 1/Å
    end type kgrid_t

    public :: build_kpath, build_kgrid, kgrid_inplane_axes

contains

    ! ==================================================================
    ! Public API
    ! ==================================================================

    !> Constrói o caminho k.
    !> vertex(3,nvertex) fracionário; npoints(nvertex-1) >= 1 intervalos por
    !> segmento; label(nvertex); recip(3,3) com linhas b_i (1/Å, inclui 2π).
    !> Segmento s: vertex_s + j/npoints(s)*(vertex_{s+1}-vertex_s),
    !> j = 0..npoints(s)-1; o vértice final é acrescentado uma única vez.
    !> nk = sum(npoints) + 1. Incrementos de kdist = norm2(matmul(dk_frac, recip)).
    !> Os pontos são uniformemente espaçados dentro de cada segmento e os vértices
    !> de junção não são duplicados, de modo que kdist é contínuo entre segmentos.
    !> Códigos de ierr: 1 vertex não é (3,n); 2 nvertex < 2; 3 tamanho de npoints;
    !> 4 tamanho de label; 5 npoints < 1; 6 segmento de comprimento zero;
    !> 7 label maior que KLABEL_LEN. Em caso de erro errmsg é preenchida e path
    !> fica vazio.
    subroutine build_kpath(vertex, npoints, label, recip, path, ierr, errmsg)
        real(dp),         intent(in)  :: vertex(:,:)
        integer,          intent(in)  :: npoints(:)
        character(len=*), intent(in)  :: label(:)
        real(dp),         intent(in)  :: recip(3,3)
        type(kpath_t),    intent(out) :: path
        integer,          intent(out) :: ierr
        character(len=*), intent(out) :: errmsg

        integer  :: nvertex, nk, s, j, k
        real(dp) :: dv(3), seglen, step, dist0

        ierr = 0
        errmsg = ''
        nvertex = size(vertex, 2)

        ! 1. Valida todas as entradas antes de alocar qualquer coisa.
        if (size(vertex, 1) /= 3) then
            ierr = 1
            errmsg = 'build_kpath: vertex must have shape (3,nvertex)'
            return
        end if
        if (nvertex < 2) then
            ierr = 2
            errmsg = 'build_kpath: at least 2 vertices are required'
            return
        end if
        if (size(npoints) /= nvertex - 1) then
            ierr = 3
            write(errmsg, '(A,I0,A,I0)') 'build_kpath: npoints must have nvertex-1 = ', &
                nvertex - 1, ' values, got ', size(npoints)
            return
        end if
        if (size(label) /= nvertex) then
            ierr = 4
            write(errmsg, '(A,I0,A,I0)') 'build_kpath: label must have nvertex = ', &
                nvertex, ' values, got ', size(label)
            return
        end if
        do s = 1, nvertex - 1
            if (npoints(s) < 1) then
                ierr = 5
                write(errmsg, '(A,I0,A,I0)') 'build_kpath: npoints(', s, ') must be >= 1, got ', &
                    npoints(s)
                return
            end if
            ! Comprimento zero cobre vértices iguais e um recip degenerado.
            seglen = norm2(matmul(vertex(:, s + 1) - vertex(:, s), recip))
            if (seglen <= 0.0_dp) then
                ierr = 6
                write(errmsg, '(A,I0,A,I0)') 'build_kpath: zero-length segment between vertices ', &
                    s, ' and ', s + 1
                return
            end if
        end do
        do s = 1, nvertex
            if (len_trim(label(s)) > KLABEL_LEN) then
                ierr = 7
                write(errmsg, '(A,I0,A,I0,A)') 'build_kpath: label(', s, ') longer than ', &
                    KLABEL_LEN, ' characters'
                return
            end if
        end do

        ! 2. Aloca: cada segmento contribui com npoints(s) pontos (início
        !    incluído, fim excluído); o último vértice acrescenta mais um.
        nk = sum(npoints) + 1
        path%nk = nk
        path%nvertex = nvertex
        allocate(path%kfrac(3, nk), path%kdist(nk), path%vertex_k(nvertex), path%label(nvertex))
        do s = 1, nvertex
            path%label(s) = label(s)
        end do

        ! 3. Preenche os pontos segmento a segmento. dist0 é a distância
        !    acumulada no início do segmento s. O comprimento é medido no
        !    espaço cartesiano (dv_frac * recip, linhas = b_i), pois passos
        !    fracionários ao longo de b_i diferentes não são equivalentes.
        k = 0
        dist0 = 0.0_dp
        do s = 1, nvertex - 1
            dv = vertex(:, s + 1) - vertex(:, s)
            seglen = norm2(matmul(dv, recip))
            step = seglen / real(npoints(s), dp)
            path%vertex_k(s) = k + 1
            do j = 0, npoints(s) - 1
                k = k + 1
                path%kfrac(:, k) = vertex(:, s) + (real(j, dp) / real(npoints(s), dp)) * dv
                path%kdist(k) = dist0 + real(j, dp) * step
            end do
            dist0 = dist0 + seglen
        end do
        ! 4. Vértice final, acrescentado uma vez (fecha o último segmento).
        k = k + 1
        path%kfrac(:, k) = vertex(:, nvertex)
        path%kdist(k) = dist0
        path%vertex_k(nvertex) = k
    end subroutine build_kpath

    !> Constrói uma grade k 2D (veja kgrid_t).
    !> origin, span1, span2 fracionários (3); n1, n2 >= 2 pontos por direção;
    !> recip(3,3) com linhas b_i (1/Å, inclui 2π).
    !> Códigos de ierr: 1 n1 ou n2 < 2; 2 span1 ou span2 de comprimento zero;
    !> 3 span1 e span2 colineares. Em caso de erro errmsg é preenchida e grid
    !> fica vazio.
    subroutine build_kgrid(origin, span1, span2, n1, n2, recip, grid, ierr, errmsg)
        real(dp),         intent(in)  :: origin(3), span1(3), span2(3)
        integer,          intent(in)  :: n1, n2
        real(dp),         intent(in)  :: recip(3,3)
        type(kgrid_t),    intent(out) :: grid
        integer,          intent(out) :: ierr
        character(len=*), intent(out) :: errmsg

        integer  :: i, j, k
        real(dp) :: eu(3), ev(3)

        ierr = 0
        errmsg = ''
        if (n1 < 2 .or. n2 < 2) then
            ierr = 1
            write(errmsg, '(A,I0,A,I0,A)') 'build_kgrid: n1 and n2 must be >= 2 (got n1=', n1, &
                ', n2=', n2, ')'
            return
        end if
        call kgrid_inplane_axes(span1, span2, recip, eu, ev, ierr, errmsg)
        if (ierr /= 0) return

        grid%n1 = n1
        grid%n2 = n2
        grid%nk = n1 * n2
        grid%origin = origin
        grid%span1 = span1
        grid%span2 = span2
        grid%e_u = eu
        grid%e_v = ev
        allocate(grid%kfrac(3, grid%nk), grid%kcart(3, grid%nk), grid%kuv(2, grid%nk))
        do j = 1, n2
            do i = 1, n1
                k = i + (j - 1) * n1
                grid%kfrac(:, k) = origin + (real(i - 1, dp) / real(n1 - 1, dp)) * span1 &
                                          + (real(j - 1, dp) / real(n2 - 1, dp)) * span2
                grid%kcart(:, k) = matmul(grid%kfrac(:, k), recip)
                grid%kuv(1, k) = dot_product(grid%kcart(:, k), eu)
                grid%kuv(2, k) = dot_product(grid%kcart(:, k), ev)
            end do
        end do
    end subroutine build_kgrid

    !> Base ortonormal (e_u, e_v) do plano gerado por span1 e span2 (fracionários):
    !> e_u ∥ span1·recip; e_v = componente de span2·recip ortogonal a e_u.
    !> ierr: 0 ok; 2 span de comprimento zero; 3 spans colineares (tolerância
    !> relativa 1e-10 em |e_v| antes da normalização).
    subroutine kgrid_inplane_axes(span1, span2, recip, eu, ev, ierr, errmsg)
        real(dp),         intent(in)  :: span1(3), span2(3), recip(3,3)
        real(dp),         intent(out) :: eu(3), ev(3)
        integer,          intent(out) :: ierr
        character(len=*), intent(out) :: errmsg

        real(dp) :: c1(3), c2(3), l1, l2, lv

        ierr = 0
        errmsg = ''
        eu = 0.0_dp
        ev = 0.0_dp
        c1 = matmul(span1, recip)
        c2 = matmul(span2, recip)
        l1 = norm2(c1)
        l2 = norm2(c2)
        if (.not. (l1 > 0.0_dp) .or. .not. (l2 > 0.0_dp)) then
            ierr = 2
            errmsg = 'k-grid: span1 and span2 must have non-zero length'
            return
        end if
        eu = c1 / l1
        ev = c2 - dot_product(c2, eu) * eu
        lv = norm2(ev)
        if (lv <= 1.0e-10_dp * l2) then
            ierr = 3
            errmsg = 'k-grid: span1 and span2 are collinear'
            eu = 0.0_dp
            ev = 0.0_dp
            return
        end if
        ev = ev / lv
    end subroutine kgrid_inplane_axes

end module kpath
