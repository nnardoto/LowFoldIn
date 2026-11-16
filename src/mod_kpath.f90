! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Caminho k linear por partes em coordenadas fracionárias recíprocas, com
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

    public :: build_kpath

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

end module kpath
