! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Blocos de camada principal (PL) a partir de H(R) e S(R).
!>
!> O bulk é descrito numa célula adaptada à superfície: dois vetores de rede no
!> plano e um vetor de empilhamento a_s (eixo `stack_axis`). Uma camada
!> principal agrupa P células consecutivas ao longo de a_s (P = pl_cells). Com
!> R = n_s a_s + R_par, o bloco entre a célula a e a célula b = a + n_s (a,b em
!> 0..P-1 dentro da PL) cai na PL vizinha p = floor(b/P) ∈ {-1, 0, 1}, desde que
!> |n_s| <= P. Essa condição é verificada: hoppings nunca são truncados.
!>
!> No gauge de rede do LowFoldIn (fase só com R):
!>
!>   H_{m,m+p}(k_par) = sum_{R_par} H_p(R_par) exp(2 pi i (k_1 n_1 + k_2 n_2)),
!>
!> em que (k_1, k_2) são coordenadas fracionárias ao longo dos dois eixos no
!> plano (na ordem crescente de eixo). Para P = 1:
!>   H00 + H01 e^{2 pi i k_s} + H10 e^{-2 pi i k_s} = H(k_1, k_2, k_s).
!>
!> Orbitais dentro de uma PL: célula a = 0..P-1 (lenta), depois a ordem do
!> .hsio (átomo-major) dentro da célula.
!>
!> Gauge atômico (opcional, exige a estrutura): U_ii = exp(i k_par . tau_i),
!> com tau_i a posição cartesiana do orbital i incluindo o deslocamento m*P*a_s
!> da PL m e k_par a projeção no plano de k_1 b_1 + k_2 b_2. G~ = U^dagger G U.
module layer_blocks
    use types, only: dp
    use HS_io, only: HSDataFormat
    use structure_io, only: structure_t, reciprocal_lattice
    implicit none
    private

    public :: pl_blocks_t, build_pl_blocks, pl_blocks_at_k, attach_pl_geometry
    public :: layer_gauge_phases, kpar_cartesian, pl_inplane_axes

    real(dp), parameter :: TWO_PI = 6.283185307179586476925286766559_dp

    type :: pl_blocks_t
        integer :: nbasis = 0        !< orbitais por célula (NOrbs do .hsio)
        integer :: ncell  = 1        !< P: células por PL
        integer :: norb   = 0        !< nbasis * P
        integer :: nat_cell = 0      !< átomos por célula
        integer :: stack_axis = 3
        integer :: inplane(2) = [1, 2]
        integer :: max_ns = 0        !< max |n_s| no .hsio original
        integer :: nr = 0            !< número de R_par distintos
        logical :: has_overlap = .false.
        integer, allocatable :: orbs_per_atom(:)     !< (nat_cell)
        integer, allocatable :: rpar(:, :)           !< (2, nr)
        complex(dp), allocatable :: h(:, :, :, :)    !< (norb, norb, nr, -1:1)
        complex(dp), allocatable :: s(:, :, :, :)    !< idem; ausente => S = I
        ! -- geometria (attach_pl_geometry) --
        logical :: has_geometry = .false.
        real(dp) :: lattice(3, 3) = 0.0_dp           !< lattice(i,:) = a_i, Å
        real(dp) :: recip(3, 3) = 0.0_dp             !< recip(i,:) = b_i, 1/Å
        real(dp) :: normal(3) = 0.0_dp               !< normal unitária ao plano
        real(dp) :: pl_shift(3) = 0.0_dp             !< P * a_s, Å
        real(dp), allocatable :: tau(:, :)           !< (3, norb) posição dos orbitais na PL 0
        integer,  allocatable :: orb_atom(:)         !< (norb) átomo (na célula) do orbital
        integer,  allocatable :: orb_cell(:)         !< (norb) célula a do orbital
    end type pl_blocks_t

contains

    !> Eixos no plano (os dois diferentes de stack_axis, em ordem crescente).
    pure function pl_inplane_axes(stack_axis) result(ax)
        integer, intent(in) :: stack_axis
        integer :: ax(2)
        select case (stack_axis)
        case (1);    ax = [2, 3]
        case (2);    ax = [1, 3]
        case default; ax = [1, 2]
        end select
    end function pl_inplane_axes

    !> Constrói os blocos de PL. hs é opcional (ausente => base ortonormal).
    !> pl_cells = 0 escolhe P = max(1, max|n_s|). Com pl_cells > 0 e
    !> max|n_s| > pl_cells devolve erro (ierr = 2) em vez de truncar.
    subroutine build_pl_blocks(hh, stack_axis, pl_cells, blocks, ierr, errmsg, hs)
        type(HSDataFormat), intent(in) :: hh
        integer,            intent(in) :: stack_axis, pl_cells
        type(pl_blocks_t),  intent(out) :: blocks
        integer,            intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        type(HSDataFormat), intent(in), optional :: hs

        integer :: nb, nrh, ir, js, ns, p, a, b, bb, irp, P_, row0, col0, k
        integer, allocatable :: rp_tmp(:, :), map_h(:), map_s(:)
        logical, allocatable :: nzh(:), nzs(:)
        integer :: ax(2)

        ierr = 0
        errmsg = ''
        if (stack_axis < 1 .or. stack_axis > 3) then
            ierr = 1; errmsg = 'build_pl_blocks: stacking_axis must be 1, 2 or 3'; return
        end if
        if (pl_cells < 0) then
            ierr = 1; errmsg = 'build_pl_blocks: pl_cells must be >= 0 (0 = automatic)'; return
        end if
        if (.not. allocated(hh%Data) .or. .not. allocated(hh%IndexList)) then
            ierr = 1; errmsg = 'build_pl_blocks: H data not loaded'; return
        end if
        nb = hh%NOrbs
        nrh = size(hh%Data, 3)
        if (size(hh%Data, 1) /= nb .or. size(hh%Data, 2) /= nb .or. size(hh%IndexList, 2) /= nrh) then
            ierr = 1; errmsg = 'build_pl_blocks: inconsistent H(R) dimensions'; return
        end if
        if (present(hs)) then
            if (.not. allocated(hs%Data)) then
                ierr = 1; errmsg = 'build_pl_blocks: S data not loaded'; return
            end if
            if (hs%NOrbs /= nb) then
                ierr = 1; errmsg = 'build_pl_blocks: H and S have different numbers of orbitals'; return
            end if
            if (allocated(hh%Orbs) .and. allocated(hs%Orbs)) then
                if (size(hh%Orbs) /= size(hs%Orbs)) then
                    ierr = 1; errmsg = 'build_pl_blocks: H and S have different numbers of atoms'; return
                end if
                if (any(hh%Orbs /= hs%Orbs)) then
                    ierr = 1; errmsg = 'build_pl_blocks: H and S have different orbitals per atom'; return
                end if
            end if
        end if

        ax = pl_inplane_axes(stack_axis)
        blocks%stack_axis = stack_axis
        blocks%inplane = ax
        blocks%nbasis = nb
        ! Alcance ao longo de a_s medido só sobre os R com blocos não nulos: o
        ! .scfout lista todas as cópias de célula do OpenMX (inclusive as que
        ! ficam além do raio de corte dos PAOs, com dados nulos).
        allocate(nzh(nrh))
        blocks%max_ns = 0
        do ir = 1, nrh
            nzh(ir) = any(hh%Data(:, :, ir) /= (0.0_dp, 0.0_dp))
            if (nzh(ir)) blocks%max_ns = max(blocks%max_ns, abs(hh%IndexList(stack_axis, ir)))
        end do
        if (present(hs)) then
            allocate(nzs(size(hs%IndexList, 2)))
            do ir = 1, size(hs%IndexList, 2)
                nzs(ir) = any(hs%Data(:, :, ir) /= (0.0_dp, 0.0_dp))
                if (nzs(ir)) blocks%max_ns = max(blocks%max_ns, abs(hs%IndexList(stack_axis, ir)))
            end do
        end if

        if (pl_cells == 0) then
            P_ = max(1, blocks%max_ns)
        else
            P_ = pl_cells
            if (blocks%max_ns > P_) then
                ierr = 2
                errmsg = 'build_pl_blocks: H(R)/S(R) couple cells up to |n_s| = ' // itoa(blocks%max_ns) // &
                         ' along the stacking axis, beyond pl_cells = ' // itoa(P_) // &
                         '; use pl_cells >= ' // itoa(blocks%max_ns) // ' (or 0 = automatic)'
                return
            end if
        end if
        blocks%ncell = P_
        blocks%norb = nb * P_

        if (allocated(hh%Orbs)) then
            blocks%nat_cell = size(hh%Orbs)
            blocks%orbs_per_atom = hh%Orbs
        else
            blocks%nat_cell = 0
            allocate(blocks%orbs_per_atom(0))
        end if

        ! Lista de R_par distintos (união de H e S)
        k = nrh
        if (present(hs)) k = k + size(hs%IndexList, 2)
        allocate(rp_tmp(2, k))
        blocks%nr = 0
        allocate(map_h(nrh))
        map_h = 0
        do ir = 1, nrh
            if (nzh(ir)) call find_or_add(hh%IndexList(ax, ir), map_h(ir))
        end do
        if (present(hs)) then
            allocate(map_s(size(hs%IndexList, 2)))
            map_s = 0
            do ir = 1, size(hs%IndexList, 2)
                if (nzs(ir)) call find_or_add(hs%IndexList(ax, ir), map_s(ir))
            end do
        end if
        if (blocks%nr == 0) call find_or_add([0, 0], k)   ! H(R) inteiramente nulo
        blocks%rpar = rp_tmp(:, 1:blocks%nr)

        allocate(blocks%h(blocks%norb, blocks%norb, blocks%nr, -1:1))
        blocks%h = (0.0_dp, 0.0_dp)
        do ir = 1, nrh
            if (.not. nzh(ir)) cycle
            ns = hh%IndexList(stack_axis, ir)
            irp = map_h(ir)
            do a = 0, P_ - 1
                b = a + ns
                p = floor_div(b, P_)
                bb = b - p * P_
                row0 = a * nb
                col0 = bb * nb
                blocks%h(row0 + 1:row0 + nb, col0 + 1:col0 + nb, irp, p) = &
                    blocks%h(row0 + 1:row0 + nb, col0 + 1:col0 + nb, irp, p) + hh%Data(:, :, ir)
            end do
        end do

        blocks%has_overlap = present(hs)
        if (present(hs)) then
            allocate(blocks%s(blocks%norb, blocks%norb, blocks%nr, -1:1))
            blocks%s = (0.0_dp, 0.0_dp)
            do js = 1, size(hs%IndexList, 2)
                if (.not. nzs(js)) cycle
                ns = hs%IndexList(stack_axis, js)
                irp = map_s(js)
                do a = 0, P_ - 1
                    b = a + ns
                    p = floor_div(b, P_)
                    bb = b - p * P_
                    row0 = a * nb
                    col0 = bb * nb
                    blocks%s(row0 + 1:row0 + nb, col0 + 1:col0 + nb, irp, p) = &
                        blocks%s(row0 + 1:row0 + nb, col0 + 1:col0 + nb, irp, p) + hs%Data(:, :, js)
                end do
            end do
        end if

        ! Índices célula/átomo de cada orbital da PL (independem da geometria)
        allocate(blocks%orb_cell(blocks%norb), blocks%orb_atom(blocks%norb))
        do a = 0, P_ - 1
            k = a * nb
            if (blocks%nat_cell > 0) then
                do ir = 1, blocks%nat_cell
                    blocks%orb_cell(k + 1:k + blocks%orbs_per_atom(ir)) = a
                    blocks%orb_atom(k + 1:k + blocks%orbs_per_atom(ir)) = ir
                    k = k + blocks%orbs_per_atom(ir)
                end do
            else
                blocks%orb_cell(k + 1:k + nb) = a
                blocks%orb_atom(k + 1:k + nb) = 0
            end if
        end do

    contains

        subroutine find_or_add(r2, idx)
            integer, intent(in)  :: r2(2)
            integer, intent(out) :: idx
            integer :: i
            do i = 1, blocks%nr
                if (all(rp_tmp(:, i) == r2)) then
                    idx = i
                    return
                end if
            end do
            blocks%nr = blocks%nr + 1
            rp_tmp(:, blocks%nr) = r2
            idx = blocks%nr
        end subroutine find_or_add

    end subroutine build_pl_blocks

    !> Blocos de PL em k_par = (k_1, k_2) (fracionário ao longo de blocks%inplane).
    !> h(:,:,p) = H_{m,m+p}, s(:,:,p) idem, p = -1, 0, 1. Sem S, s(:,:,0) = I e
    !> s(:,:,±1) = 0.
    subroutine pl_blocks_at_k(blocks, kpar, h, s)
        type(pl_blocks_t), intent(in)  :: blocks
        real(dp),          intent(in)  :: kpar(2)
        complex(dp),       intent(out) :: h(:, :, -1:)
        complex(dp),       intent(out) :: s(:, :, -1:)
        integer :: ir, p, i
        complex(dp) :: ph

        h(:, :, -1:1) = (0.0_dp, 0.0_dp)
        s(:, :, -1:1) = (0.0_dp, 0.0_dp)
        do ir = 1, blocks%nr
            ph = exp(cmplx(0.0_dp, TWO_PI * (kpar(1) * real(blocks%rpar(1, ir), dp) + &
                                             kpar(2) * real(blocks%rpar(2, ir), dp)), dp))
            do p = -1, 1
                h(:, :, p) = h(:, :, p) + ph * blocks%h(:, :, ir, p)
                if (blocks%has_overlap) s(:, :, p) = s(:, :, p) + ph * blocks%s(:, :, ir, p)
            end do
        end do
        if (.not. blocks%has_overlap) then
            do i = 1, blocks%norb
                s(i, i, 0) = (1.0_dp, 0.0_dp)
            end do
        end if
    end subroutine pl_blocks_at_k

    !> Associa a geometria (estrutura da célula do bulk) aos blocos: posições
    !> dos orbitais da PL 0, normal ao plano, deslocamento entre PLs.
    subroutine attach_pl_geometry(blocks, st, ierr, errmsg)
        type(pl_blocks_t), intent(inout) :: blocks
        type(structure_t), intent(in)    :: st
        integer,           intent(out)   :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        integer :: i, at, a
        real(dp) :: nrm(3), nn

        ierr = 0
        errmsg = ''
        if (blocks%nat_cell <= 0) then
            ierr = 1; errmsg = 'attach_pl_geometry: H(R) has no per-atom orbital counts'; return
        end if
        if (st%nat /= blocks%nat_cell) then
            ierr = 1
            errmsg = 'attach_pl_geometry: structure has ' // itoa(st%nat) // ' atoms but H(R) has ' // &
                     itoa(blocks%nat_cell)
            return
        end if
        blocks%lattice = st%lattice
        blocks%recip = reciprocal_lattice(st)
        nrm = cross(st%lattice(blocks%inplane(1), :), st%lattice(blocks%inplane(2), :))
        nn = norm2(nrm)
        if (nn <= 1.0e-12_dp) then
            ierr = 1; errmsg = 'attach_pl_geometry: in-plane lattice vectors are collinear'; return
        end if
        nrm = nrm / nn
        ! Orienta a normal no sentido de a_s
        if (dot_product(nrm, st%lattice(blocks%stack_axis, :)) < 0.0_dp) nrm = -nrm
        blocks%normal = nrm
        blocks%pl_shift = real(blocks%ncell, dp) * st%lattice(blocks%stack_axis, :)
        if (allocated(blocks%tau)) deallocate(blocks%tau)
        allocate(blocks%tau(3, blocks%norb))
        do i = 1, blocks%norb
            at = blocks%orb_atom(i)
            a = blocks%orb_cell(i)
            blocks%tau(:, i) = st%cart(:, at) + real(a, dp) * st%lattice(blocks%stack_axis, :)
        end do
        blocks%has_geometry = .true.
    end subroutine attach_pl_geometry

    !> k_par cartesiano (1/Å): projeção no plano de k_1 b_1 + k_2 b_2.
    pure function kpar_cartesian(blocks, kpar) result(kc)
        type(pl_blocks_t), intent(in) :: blocks
        real(dp),          intent(in) :: kpar(2)
        real(dp) :: kc(3)
        kc = kpar(1) * blocks%recip(blocks%inplane(1), :) + kpar(2) * blocks%recip(blocks%inplane(2), :)
        kc = kc - dot_product(kc, blocks%normal) * blocks%normal
    end function kpar_cartesian

    !> Fases do gauge atômico para os orbitais da PL m (m pode ser negativo):
    !> phase(i) = exp(i k_par . (tau_i + m * pl_shift)). origin (Å) é subtraído
    !> das posições (opcional), para fixar a referência de fase global.
    subroutine layer_gauge_phases(blocks, kpar, m, phase, origin)
        type(pl_blocks_t), intent(in)  :: blocks
        real(dp),          intent(in)  :: kpar(2)
        integer,           intent(in)  :: m
        complex(dp),       intent(out) :: phase(:)
        real(dp), optional, intent(in) :: origin(3)
        real(dp) :: kc(3), r0(3)
        integer :: i
        kc = kpar_cartesian(blocks, kpar)
        r0 = 0.0_dp
        if (present(origin)) r0 = origin
        do i = 1, blocks%norb
            phase(i) = exp(cmplx(0.0_dp, dot_product(kc, blocks%tau(:, i) + &
                           real(m, dp) * blocks%pl_shift - r0), dp))
        end do
    end subroutine layer_gauge_phases

    ! ------------------------------------------------------------------
    pure integer function floor_div(a, b)
        integer, intent(in) :: a, b
        floor_div = (a - modulo(a, b)) / b
    end function floor_div

    pure function cross(u, v) result(w)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: w(3)
        w = [u(2) * v(3) - u(3) * v(2), u(3) * v(1) - u(1) * v(3), u(1) * v(2) - u(2) * v(1)]
    end function cross

    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=16) :: buf
        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

end module layer_blocks
