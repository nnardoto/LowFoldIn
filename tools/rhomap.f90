! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> rhomap: densidade eletrônica de janela em espaço real (base PAO do OpenMX).
!>
!> Fontes da matriz densidade D_ij (por pares de átomos, gauge atômico):
!>   source = 'scfout': autovetores de H(k), S(k) do .scfout (supercélula),
!>       pesos do partial.charge do OpenMX 4.0 (Band_DFT_Dosout.c):
!>       |f((e-mu)/kT) - f((e-mu-W)/kT)|, x cortado em +-50, k = i/n da grade
!>       &WINDOW kgrid com peso 1/N, D = 2 PDM (sem spin). W = e_window_ev.
!>   source = 'pairs': PREFIX.rho_pairs.dat do surfgreen extracted (&RHO),
!>       janela pairs_window, D = spin_factor * D_pairs.
!> n(r) = sum D phi phi (density_grid) numa grade &GRID; mode 'integrate'
!> integra ao longo de dw (nw pontos; e/Å^2), mode 'plane' dá n no plano
!> (e/Å^3), mode 'cube' dá n em toda a grade 3D (A8): PREFIX.cube.bin (stream real64,
!> n(iu,iv,iw) em e/Bohr^3, iu mais rápido) e PREFIX.cube.hdr (texto, origem e passos em Bohr).
!> mode 'elf' (A9) grava, além do PREFIX.cube.bin/.hdr, PREFIX.grad.bin (grad n (3,nu,nv,nw), e/Bohr^4,
!> gradiente analítico das PAO), PREFIX.tau.bin (tau = 1/2 sum D grad phi . grad phi, Ha/Bohr^3) e
!> PREFIX.elf.bin (ELF de camada fechada, 1/(1 + (Dp/Dh)^2), Dp = tau - |grad n|^2/(8n),
!> Dh = (3/10)(3 pi^2)^(2/3) n^(5/3); -1 onde n < 1e-5 e/Bohr^3, ELF indefinida no vácuo).
!> Saídas: PREFIX.map.dat, PREFIX.atoms.dat.
!>
!> Exemplo:
!>   &RHOMAP source='scfout', scfout_file='loneChain.scfout', output_prefix='m05' /
!>   &SPECIES nspec=2, label='Te','E', pao_file='Te7.0.pao','Te7.0.pao', basis='s3p2d2','s3p2d2',
!>            rcut_bohr=7.0,7.0 /
!>   &WINDOW e_window_ev=-0.5, temperature_k=300, kgrid=1,1,3 /
!>   &GRID mode='integrate', units='bohr', origin=..., du=..., dv=..., dw=..., nu=405, nv=648, nw=54 /
program rhomap_cli
    use types, only: dp
    use scfout_io, only: ScfoutData, ReadScfoutFile
    use pao_basis, only: pao_species_t, read_pao_species
    use density_grid, only: dens_system_t, density_map, dens_set_pairs
    !$ use omp_lib
    implicit none

    interface
        subroutine zhegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info)
            use types, only: dp
            integer :: itype, n, lda, ldb, lwork, lrwork, liwork, info
            character :: jobz, uplo
            complex(dp) :: a(lda, *), b(ldb, *), work(*)
            real(dp) :: w(*), rwork(*)
            integer :: iwork(*)
        end subroutine zhegvd
    end interface

    real(dp), parameter :: BOHR = 0.529177249_dp, HA_EV = 27.2113845_dp, KB_EV = 0.00008617251324_dp
    real(dp), parameter :: PI = acos(-1.0_dp)
    integer, parameter :: MAXSP = 8

    character(len=16) :: source = 'scfout', mode = 'integrate', units = 'ang'
    character(len=512) :: scfout_file = '', pairs_file = '', output_prefix = 'rhomap', inp_file
    character(len=16) :: pairs_window = ''
    real(dp) :: spin_factor = 2.0_dp
    integer :: nspec = 0
    character(len=16) :: label(MAXSP) = ''
    character(len=512) :: pao_file(MAXSP) = ''
    character(len=16) :: basis(MAXSP) = ''
    real(dp) :: rcut_bohr(MAXSP) = 0.0_dp
    real(dp) :: e_window_ev = 0.0_dp, temperature_k = 300.0_dp
    integer :: kgrid(3) = [1, 1, 1]
    real(dp) :: origin(3) = 0.0_dp, du(3) = 0.0_dp, dv(3) = 0.0_dp, dw(3) = 0.0_dp
    integer :: nu = 0, nv = 0, nw = 1
    integer :: nplanes = 0, plane_iw(16) = 0      ! planos extras (n em iw dado), saída PREFIX.plane<iw>.dat

    namelist /RHOMAP/ source, scfout_file, pairs_file, pairs_window, spin_factor, output_prefix
    namelist /SPECIES/ nspec, label, pao_file, basis, rcut_bohr
    namelist /WINDOW/ e_window_ev, temperature_k, kgrid
    namelist /GRID/ mode, units, origin, du, dv, dw, nu, nv, nw, nplanes, plane_iw

    type(dens_system_t) :: ds
    real(dp), allocatable :: map(:, :)
    integer :: u, ios, ierr, i, j, s
    character(len=:), allocatable :: msg
    real(dp) :: fac, t0, t1, area, trs
    logical :: integ, cubem
    real(dp), allocatable :: cub(:, :, :), grd(:, :, :, :), tau(:, :, :)
    character(len=64) :: line_

    if (command_argument_count() < 1) then
        write(*, '(A)') 'usage: rhomap input.nml'
        stop 1
    end if
    call get_command_argument(1, inp_file)
    open(newunit=u, file=trim(inp_file), status='old', action='read', iostat=ios)
    if (ios /= 0) call die('cannot open ' // trim(inp_file))
    read(u, nml=RHOMAP, iostat=ios); if (ios /= 0) call die('bad &RHOMAP')
    rewind(u); read(u, nml=SPECIES, iostat=ios); if (ios /= 0) call die('bad &SPECIES')
    rewind(u); read(u, nml=GRID, iostat=ios); if (ios /= 0) call die('bad &GRID')
    if (trim(source) == 'scfout') then
        rewind(u); read(u, nml=WINDOW, iostat=ios); if (ios /= 0) call die('bad &WINDOW')
    end if
    close(u)
    if (nspec < 1 .or. nspec > MAXSP) call die('&SPECIES nspec out of range')
    if (nu < 1 .or. nv < 1 .or. nw < 1) call die('&GRID nu, nv, nw must be >= 1')
    integ = trim(mode) == 'integrate'
    cubem = trim(mode) == 'cube' .or. trim(mode) == 'elf'
    if (.not. integ .and. .not. cubem .and. trim(mode) /= 'plane') &
        call die('&GRID mode must be integrate, plane, cube or elf')
    if (trim(units) == 'ang') then
        origin = origin / BOHR; du = du / BOHR; dv = dv / BOHR; dw = dw / BOHR
    else if (trim(units) /= 'bohr') then
        call die('&GRID units must be ang or bohr')
    end if
    if (.not. integ .and. .not. cubem) nw = 1

    allocate(ds%sp(nspec))
    do s = 1, nspec
        call read_pao_species(trim(pao_file(s)), trim(basis(s)), rcut_bohr(s), ds%sp(s), ierr, msg)
        if (ierr /= 0) call die(msg)
        write(*, '(A,A,A,I0,A,F6.2,A)') 'species ', trim(label(s)), ': ', ds%sp(s)%norb, ' orbitals, rcut ', &
            rcut_bohr(s), ' bohr'
    end do

    call cpu_wall(t0)
    if (trim(source) == 'scfout') then
        call dm_from_scfout()
    else if (trim(source) == 'pairs') then
        call dm_from_pairs()
    else
        call die('&RHOMAP source must be scfout or pairs')
    end if
    call dens_set_pairs(ds)
    call cpu_wall(t1)
    write(*, '(A,I0,A,I0,A,F10.1,A)') 'density matrix: ', ds%nsite, ' sites, ', ds%npair, ' pairs (', t1 - t0, ' s)'

    allocate(map(nu, nv))
    call cpu_wall(t0)
    if (cubem) then
        allocate(cub(nu, nv, nw))
        if (trim(mode) == 'elf') then
            allocate(grd(3, nu, nv, nw), tau(nu, nv, nw))
            call density_map(ds, origin, du, dv, dw, nu, nv, nw, .true., map, ierr, msg, cube=cub, grad=grd, tau=tau)
        else
            call density_map(ds, origin, du, dv, dw, nu, nv, nw, .true., map, ierr, msg, cube=cub)
        end if
        if (ierr /= 0) call die(msg)
        call cpu_wall(t1)
        write(*, '(A,F10.1,A)') 'grid: ', t1 - t0, ' s'
        write(*, '(A,ES24.16,A)') 'sum cube * dV = ', sum(cub) * abs(dot_product(cross(du, dv), dw)), ' e'
        open(newunit=u, file=trim(output_prefix) // '.cube.hdr', status='replace', action='write')
        write(u, '(A)') '# rhomap cube: n(iu,iv,iw) real64 stream in PREFIX.cube.bin, iu fastest, e/Bohr^3'
        write(u, '(A,3ES24.16)') 'origin_bohr ', origin
        write(u, '(A,3ES24.16)') 'du_bohr ', du
        write(u, '(A,3ES24.16)') 'dv_bohr ', dv
        write(u, '(A,3ES24.16)') 'dw_bohr ', dw
        write(u, '(A,3(1X,I0))') 'n ', nu, nv, nw
        close(u)
        open(newunit=u, file=trim(output_prefix) // '.cube.bin', status='replace', action='write', &
             access='stream', form='unformatted')
        write(u) cub
        close(u)
        if (trim(mode) == 'elf') call write_elf()
        stop
    end if
    call density_map(ds, origin, du, dv, dw, nu, nv, nw, integ, map, ierr, msg)
    if (ierr /= 0) call die(msg)
    call cpu_wall(t1)
    area = norm2(cross(du, dv))
    write(*, '(A,F10.1,A)') 'grid: ', t1 - t0, ' s'
    if (integ) then
        write(*, '(A,ES24.16,A)') 'sum map * dA = ', sum(map) * area, ' e (whole cell if the grid spans it)'
        fac = 1.0_dp / BOHR**2
    else
        fac = 1.0_dp / BOHR**3
    end if
    call write_outputs(trim(output_prefix) // '.map.dat', integ, origin)
    do i = 1, nplanes
        call density_map(ds, origin + plane_iw(i) * dw, du, dv, dw, nu, nv, 1, .false., map, ierr, msg)
        if (ierr /= 0) call die(msg)
        fac = 1.0_dp / BOHR**3
        write(line_, '(A,I0,A)') '.plane', plane_iw(i), '.dat'
        call write_outputs(trim(output_prefix) // trim(line_), .false., origin + plane_iw(i) * dw)
    end do

contains

    !> A9: grava grad n, tau e a ELF (mode 'elf'); confere tau >= tau_W e 0 <= ELF <= 1.
    subroutine write_elf()
        real(dp), parameter :: RHOMIN = 1.0e-5_dp
        real(dp), allocatable :: elf(:, :, :)
        real(dp) :: ck, n, tw, dp_, dh, rmin
        integer :: iu, iv, iw, nmask
        ck = 0.3_dp * (3.0_dp * PI * PI)**(2.0_dp / 3.0_dp)
        allocate(elf(nu, nv, nw))
        nmask = 0
        rmin = huge(1.0_dp)
        do iw = 1, nw
            do iv = 1, nv
                do iu = 1, nu
                    n = cub(iu, iv, iw)
                    if (n < RHOMIN) then
                        elf(iu, iv, iw) = -1.0_dp
                        nmask = nmask + 1
                        cycle
                    end if
                    tw = sum(grd(:, iu, iv, iw)**2) / (8.0_dp * n)
                    dp_ = tau(iu, iv, iw) - tw
                    rmin = min(rmin, dp_ / tau(iu, iv, iw))
                    dh = ck * n**(5.0_dp / 3.0_dp)
                    elf(iu, iv, iw) = 1.0_dp / (1.0_dp + (dp_ / dh)**2)
                end do
            end do
        end do
        write(*, '(A,I0,A,ES10.3,A)') 'elf: ', nmask, ' points masked (n < ', RHOMIN, ' e/Bohr^3)'
        write(*, '(A,ES12.4)') 'elf: min (tau - tau_W)/tau over unmasked points = ', rmin
        write(*, '(A,2F10.6)') 'elf: min, max ELF (unmasked) = ', minval(elf, mask=elf >= 0.0_dp), maxval(elf)
        open(newunit=u, file=trim(output_prefix) // '.grad.bin', status='replace', action='write', &
             access='stream', form='unformatted')
        write(u) grd
        close(u)
        open(newunit=u, file=trim(output_prefix) // '.tau.bin', status='replace', action='write', &
             access='stream', form='unformatted')
        write(u) tau
        close(u)
        open(newunit=u, file=trim(output_prefix) // '.elf.bin', status='replace', action='write', &
             access='stream', form='unformatted')
        write(u) elf
        close(u)
    end subroutine write_elf

    subroutine die(m)
        character(len=*), intent(in) :: m
        write(*, '(A)') 'rhomap: ERROR: ' // trim(m)
        stop 1
    end subroutine die

    subroutine cpu_wall(t)
        real(dp), intent(out) :: t
        integer(8) :: c, r
        call system_clock(c, r)
        t = real(c, dp) / real(r, dp)
    end subroutine cpu_wall

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
    end function cross

    integer function spec_of(lab) result(k)
        character(len=*), intent(in) :: lab
        do k = 1, nspec
            if (trim(label(k)) == trim(lab)) return
        end do
        call die('no &SPECIES entry for label ' // trim(lab))
    end function spec_of

    !> Matriz densidade de janela pelos autovetores (partial.charge do OpenMX).
    subroutine dm_from_scfout()
        type(ScfoutData) :: sd
        integer :: nat, n, ik, nk, a, h, p, gh, ni, nj, k1, k2, k3, mu, nu_, ns, l, info, lwork, lrwork, liwork
        integer, allocatable :: off(:), sel(:), iwork(:)
        integer(8) :: o, nd
        complex(dp), allocatable :: hk(:, :), sk(:, :), work(:), cw(:, :), m(:, :)
        real(dp), allocatable :: ev(:), rwork(:), wgt(:), sblk(:)
        real(dp) :: kf(3), th, beta, x1, x2, mu_ch, wwin, tr_w, tr_ds, f1, f2
        complex(dp) :: ph, wq(1)
        real(dp) :: rq(1)
        integer :: iq(1)
        call ReadScfoutFile(sd, trim(scfout_file))
        if (sd%SpinP_switch /= 0) call die('scfout: only SpinP_switch = 0 is supported')
        if (.not. sd%has_atom_labels) call die('scfout: no atom labels in the embedded input')
        nat = sd%atomnum
        ds%nsite = nat
        allocate(ds%pos(3, nat), ds%spec(nat), off(nat + 1))
        ds%pos = sd%cart_bohr
        ds%nper = 3
        ds%lat = transpose(sd%lattice_bohr)
        off(1) = 0
        do a = 1, nat
            ds%spec(a) = spec_of(sd%atom_label(a))
            if (ds%sp(ds%spec(a))%norb /= sd%TotalNumOrbs(a)) call die('scfout: basis size differs from &SPECIES')
            off(a + 1) = off(a) + sd%TotalNumOrbs(a)
        end do
        n = off(nat + 1)
        ds%npair = sd%n_pairs
        allocate(ds%pi(ds%npair), ds%pj(ds%npair), ds%pn(3, ds%npair), ds%poff(ds%npair))
        p = 0
        nd = 0
        do a = 1, nat
            do h = 0, sd%FNAN(a)
                p = p + 1
                gh = sd%natn_flat(sd%neighbor_offset(a) + h + 1)
                ds%pi(p) = a; ds%pj(p) = gh
                ds%pn(:, p) = sd%atv_ijk(:, sd%ncn_flat(sd%neighbor_offset(a) + h + 1))
                ds%poff(p) = nd
                nd = nd + int(sd%TotalNumOrbs(a), 8) * sd%TotalNumOrbs(gh)
            end do
        end do
        allocate(ds%dmat(nd))
        ds%dmat = 0.0_dp
        beta = 1.0_dp / (KB_EV * (temperature_k / HA_EV))
        mu_ch = sd%ChemP
        wwin = e_window_ev / HA_EV
        nk = product(kgrid)
        write(*, '(A,I0,A,F14.10,A,F8.3,A,I0)') 'scfout: N = ', n, ', ChemP = ', mu_ch, ' Ha, window ', &
            e_window_ev, ' eV, nk = ', nk
        allocate(hk(n, n), sk(n, n), ev(n), wgt(n), sel(n))
        tr_w = 0.0_dp
        do ik = 0, nk - 1
            k1 = ik / (kgrid(2) * kgrid(3)); k2 = mod(ik / kgrid(3), kgrid(2)); k3 = mod(ik, kgrid(3))
            ! Band_DFT_Dosout.c: k = i/n (+ Shift_K_Point = 1e-12, desprezível)
            kf = [real(k1, dp) / kgrid(1), real(k2, dp) / kgrid(2), real(k3, dp) / kgrid(3)]
            hk = (0.0_dp, 0.0_dp); sk = (0.0_dp, 0.0_dp)
            do p = 1, ds%npair
                a = ds%pi(p); gh = ds%pj(p)
                ni = sd%TotalNumOrbs(a); nj = sd%TotalNumOrbs(gh)
                th = 2.0_dp * PI * dot_product(kf, real(ds%pn(:, p), dp))
                ph = cmplx(cos(th), sin(th), dp)
                o = sd%block_offset(p)
                do mu = 1, ni
                    do nu_ = 1, nj
                        hk(off(a) + mu, off(gh) + nu_) = hk(off(a) + mu, off(gh) + nu_) + &
                            ph * sd%Hks_flat(o + (mu - 1) * nj + nu_, 0)
                        sk(off(a) + mu, off(gh) + nu_) = sk(off(a) + mu, off(gh) + nu_) + &
                            ph * sd%OLP_flat(o + (mu - 1) * nj + nu_)
                    end do
                end do
            end do
            call zhegvd(1, 'V', 'U', n, hk, n, sk, n, ev, wq, -1, rq, -1, iq, -1, info)
            lwork = int(real(wq(1))); lrwork = int(rq(1)); liwork = iq(1)
            allocate(work(lwork), rwork(lrwork), iwork(liwork))
            call zhegvd(1, 'V', 'U', n, hk, n, sk, n, ev, work, lwork, rwork, lrwork, iwork, liwork, info)
            deallocate(work, rwork, iwork)
            if (info /= 0) call die('zhegvd failed')
            ns = 0
            do l = 1, n
                x1 = max(-50.0_dp, min(50.0_dp, (ev(l) - mu_ch) * beta))
                x2 = max(-50.0_dp, min(50.0_dp, (ev(l) - (mu_ch + wwin)) * beta))
                f1 = 1.0_dp / (1.0_dp + exp(x1)); f2 = 1.0_dp / (1.0_dp + exp(x2))
                if (abs(f1 - f2) > 1.0e-15_dp) then
                    ns = ns + 1; sel(ns) = l; wgt(ns) = abs(f1 - f2)
                end if
            end do
            tr_w = tr_w + sum(wgt(1:ns)) / nk
            write(*, '(A,3F8.4,A,I0,A,F12.6)') 'k ', kf, ': states in window ', ns, ', sum w = ', sum(wgt(1:ns))
            allocate(cw(n, ns))
            do l = 1, ns
                cw(:, l) = conjg(hk(:, sel(l))) * wgt(l)
            end do
            ! PDM_ij(R) += Re[ e^{i 2 pi k.R} sum_l w_l conj(c_il) c_jl ] / nk
            !$omp parallel do default(shared) private(p, a, gh, ni, nj, th, ph, m, l, mu, nu_) schedule(dynamic, 64)
            do p = 1, ds%npair
                a = ds%pi(p); gh = ds%pj(p)
                ni = sd%TotalNumOrbs(a); nj = sd%TotalNumOrbs(gh)
                th = 2.0_dp * PI * dot_product(kf, real(ds%pn(:, p), dp))
                ph = cmplx(cos(th), sin(th), dp)
                allocate(m(ni, nj))
                m = matmul(cw(off(a) + 1:off(a) + ni, :), transpose(hk(off(gh) + 1:off(gh) + nj, sel(1:ns))))
                do nu_ = 1, nj
                    do mu = 1, ni
                        ds%dmat(ds%poff(p) + (nu_ - 1) * ni + mu) = ds%dmat(ds%poff(p) + (nu_ - 1) * ni + mu) + &
                            real(ph * m(mu, nu_), dp) / nk
                    end do
                end do
                deallocate(m)
            end do
            !$omp end parallel do
            deallocate(cw)
        end do
        ! tr(PDM S) contra sum_k w_k sum_l w_l
        tr_ds = 0.0_dp
        do p = 1, ds%npair
            a = ds%pi(p); gh = ds%pj(p)
            ni = sd%TotalNumOrbs(a); nj = sd%TotalNumOrbs(gh)
            o = sd%block_offset(p)
            do mu = 1, ni
                do nu_ = 1, nj
                    tr_ds = tr_ds + ds%dmat(ds%poff(p) + (nu_ - 1) * ni + mu) * sd%OLP_flat(o + (mu - 1) * nj + nu_)
                end do
            end do
        end do
        write(*, '(A,ES24.16,A,ES24.16,A,ES10.3)') 'tr(PDM S) = ', tr_ds, ' ; sum of weights = ', tr_w, &
            ' ; rel. diff ', abs(tr_ds - tr_w) / max(abs(tr_w), tiny(1.0_dp))
        ds%dmat = 2.0_dp * ds%dmat
        write(*, '(A,ES24.16)') 'electrons in window (x2 spin) = ', 2.0_dp * tr_ds
    end subroutine dm_from_scfout

    !> Matriz densidade dos pares do surfgreen extracted (&RHO).
    subroutine dm_from_pairs()
        integer :: v, q, npr, kind_, ia, icol, ja, ni, nj, nper, a, per(3), k, si, sj, ninv(3), nj_(3), mu
        real(dp) :: lat(3, 3), d(3), dn, ri(3), x(3)
        real(dp), allocatable :: pos(:, :), blk(:)
        integer, allocatable :: pspec(:), qi(:), qj(:), qn(:, :), qni(:), qnj(:)
        character(len=1024) :: line
        character(len=16) :: li, lj, word, wname
        integer(8) :: nd
        real(dp) :: pm(3, 3), g(3, 3), gi(3, 3)
        logical :: inwin
        open(newunit=v, file=trim(pairs_file), status='old', action='read', iostat=ios)
        if (ios /= 0) call die('cannot open ' // trim(pairs_file))
        npr = 0; nper = 0
        do
            read(v, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:8) == 'lattice ') then
                read(line(9:), *) a, lat(a, :), per(a)
            else if (line(1:5) == 'pair ') then
                npr = npr + 1
            end if
        end do
        if (npr == 0) call die('no pairs in ' // trim(pairs_file))
        ds%nper = 0
        do a = 1, 3
            if (per(a) == 1) then
                ds%nper = ds%nper + 1
                ds%lat(:, ds%nper) = lat(a, :) / BOHR
            end if
        end do
        nper = ds%nper
        pm = 0.0_dp
        pm(:, 1:nper) = ds%lat(:, 1:nper)
        g = matmul(transpose(pm), pm)
        do a = nper + 1, 3
            g(a, a) = 1.0_dp
        end do
        call inv3(g, gi)
        allocate(pos(3, 2 * npr), pspec(2 * npr), qi(npr), qj(npr), qn(3, npr), qni(npr), qnj(npr))
        ds%nsite = 0
        rewind(v)
        q = 0
        do
            read(v, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:5) /= 'pair ') cycle
            q = q + 1
            read(line(6:), *) k, kind_, ia, icol, ja, ni, nj, d, dn, ri, li, lj
            call site_of(ri / BOHR, spec_of(li), pos, pspec, pm, gi, nper, si, ninv)
            call site_of((ri + d) / BOHR, spec_of(lj), pos, pspec, pm, gi, nper, sj, nj_)
            qi(q) = si; qj(q) = sj; qn(:, q) = nj_ - ninv; qni(q) = ni; qnj(q) = nj
        end do
        ds%npair = npr
        allocate(ds%pos(3, ds%nsite), ds%spec(ds%nsite))
        ds%pos = pos(:, 1:ds%nsite); ds%spec = pspec(1:ds%nsite)
        allocate(ds%pi(npr), ds%pj(npr), ds%pn(3, npr), ds%poff(npr))
        ds%pi = qi; ds%pj = qj; ds%pn = qn
        nd = 0
        do q = 1, npr
            ds%poff(q) = nd
            nd = nd + int(qni(q), 8) * qnj(q)
        end do
        allocate(ds%dmat(nd))
        ds%dmat = 0.0_dp
        ! blocos da janela pedida
        rewind(v)
        inwin = .false.
        q = 0
        do
            read(v, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:7) == 'window ') then
                read(line(8:), *) wname
                inwin = trim(wname) == trim(pairs_window)
                cycle
            end if
            if (.not. inwin .or. line(1:6) /= 'block ') cycle
            read(line(7:), *) q
            allocate(blk(qnj(q)))
            do mu = 1, qni(q)
                read(v, *) blk
                ds%dmat(ds%poff(q) + mu:ds%poff(q) + int(qni(q), 8) * qnj(q):qni(q)) = spin_factor * blk
            end do
            deallocate(blk)
        end do
        close(v)
        if (q == 0) call die('window ' // trim(pairs_window) // ' not found in ' // trim(pairs_file))
        write(*, '(A,I0,A)') 'pairs: ', nper, ' periodic axes'
    end subroutine dm_from_pairs

    !> Sítio de r (Bohr): procura módulo os eixos periódicos; cria se novo. nimg: coeficientes da imagem.
    subroutine site_of(r, sp_, pos, pspec, pm, gi, nper, s_, nimg)
        real(dp), intent(in) :: r(3), pm(3, 3), gi(3, 3)
        integer, intent(in) :: sp_, nper
        real(dp), intent(inout) :: pos(:, :)
        integer, intent(inout) :: pspec(:)
        integer, intent(out) :: s_, nimg(3)
        real(dp) :: dd(3), c(3), res(3)
        integer :: t
        nimg = 0
        do t = 1, ds%nsite
            if (pspec(t) /= sp_) cycle
            dd = r - pos(:, t)
            c = 0.0_dp
            if (nper > 0) c = matmul(gi, matmul(transpose(pm), dd))
            c(nper + 1:3) = 0.0_dp
            res = dd - matmul(pm, real(nint(c), dp))
            if (norm2(res) < 1.0e-4_dp) then
                s_ = t; nimg = nint(c); return
            end if
        end do
        ds%nsite = ds%nsite + 1
        s_ = ds%nsite
        pos(:, s_) = r; pspec(s_) = sp_
    end subroutine site_of

    subroutine inv3(a, b)
        real(dp), intent(in) :: a(3, 3)
        real(dp), intent(out) :: b(3, 3)
        real(dp) :: det
        det = a(1, 1) * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)) - a(1, 2) * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1)) &
            + a(1, 3) * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
        b(1, 1) = a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2); b(1, 2) = a(1, 3) * a(3, 2) - a(1, 2) * a(3, 3)
        b(1, 3) = a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2); b(2, 1) = a(2, 3) * a(3, 1) - a(2, 1) * a(3, 3)
        b(2, 2) = a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1); b(2, 3) = a(1, 3) * a(2, 1) - a(1, 1) * a(2, 3)
        b(3, 1) = a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1); b(3, 2) = a(1, 2) * a(3, 1) - a(1, 1) * a(3, 2)
        b(3, 3) = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
        b = b / det
    end subroutine inv3

    subroutine write_outputs(fname, integ_, org)
        character(len=*), intent(in) :: fname
        logical, intent(in) :: integ_
        real(dp), intent(in) :: org(3)
        integer :: w, iu, iv, n1, n2, n3, i
        real(dp) :: r(3), bmin(3), bmax(3), p(3)
        open(newunit=w, file=fname, status='replace', action='write')
        write(w, '(A,A,A,L1)') '# rhomap source=', trim(source), ' integrate=', integ_
        if (integ_) then
            write(w, '(A)') '# value: n integrated along dw (e/Angstrom^2)'
        else
            write(w, '(A)') '# value: n on the plane (e/Angstrom^3)'
        end if
        write(w, '(A,3ES24.16)') '# origin_ang ', org * BOHR
        write(w, '(A,3ES24.16)') '# du_ang ', du * BOHR
        write(w, '(A,3ES24.16)') '# dv_ang ', dv * BOHR
        write(w, '(A,3ES24.16,1X,I0)') '# dw_ang nw ', dw * BOHR, nw
        write(w, '(A,I0,1X,I0)') '# nu nv ', nu, nv
        write(w, '(A)') '# iu iv x y z value (x, y, z: Angstrom, point iw = 0)'
        do iv = 1, nv
            do iu = 1, nu
                r = (org + (iu - 1) * du + (iv - 1) * dv) * BOHR
                write(w, '(I0,1X,I0,3(1X,F12.6),1X,ES24.16)') iu - 1, iv - 1, r, map(iu, iv) * fac
            end do
        end do
        close(w)
        if (.not. integ_ .and. integ) return
        ! átomos (e imagens) na caixa da grade
        bmin = huge(1.0_dp); bmax = -huge(1.0_dp)
        do iu = 0, 1
            do iv = 0, 1
                do n3 = 0, 1
                    r = origin + iu * (nu - 1) * du + iv * (nv - 1) * dv + n3 * (nw - 1) * dw
                    bmin = min(bmin, r); bmax = max(bmax, r)
                end do
            end do
        end do
        open(newunit=w, file=trim(output_prefix) // '.atoms.dat', status='replace', action='write')
        write(w, '(A)') '# site label x y z (Angstrom), atoms and periodic images inside the grid box'
        do i = 1, ds%nsite
            do n1 = -3, 3
                do n2 = -3, 3
                    do n3 = -3, 3
                        if (ds%nper < 1 .and. n1 /= 0) cycle
                        if (ds%nper < 2 .and. n2 /= 0) cycle
                        if (ds%nper < 3 .and. n3 /= 0) cycle
                        p = ds%pos(:, i) + n1 * ds%lat(:, 1) + n2 * ds%lat(:, 2) + n3 * ds%lat(:, 3)
                        if (any(p < bmin - 1.0e-6_dp) .or. any(p > bmax + 1.0e-6_dp)) cycle
                        write(w, '(I0,1X,A,3(1X,F12.6))') i, trim(label(ds%spec(i))), p * BOHR
                    end do
                end do
            end do
        end do
        close(w)
    end subroutine write_outputs

end program rhomap_cli
