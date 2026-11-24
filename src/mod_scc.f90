! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Ciclo eletrostático autoconsistente do modo 'extracted' (A5, "full LowFoldIn").
!>
!> Modelo: H(V) = H_0 + 1/2 {S, V} (V escalar por átomo central, forma do &SITEPOT),
!>   V_i = sum_j Gamma_ij dq_j,   dq_j = s * q_j - q_ref(j),
!> com q_j a carga de Mulliken com termos cruzados (q_DS do &RHO) integrada por contorno
!> de e_lo a e_fermi (E_F fixo do volume, grande-canônico) e s o fator de spin.
!> Gamma sem ajuste: interação coulombiana entre densidades esféricas
!>   rho(r) = sum_m |phi_{l,mul,m}(r)|^2 / (2l+1)   (radial do .pao, normalizada a 1),
!> on-site U_H = Gamma(0) (Hartree, sem xc, sem epsilon), Gamma(R) = 1/R para R >= 2 rcut.
!> Longo alcance: soma de Ewald 2D nos eixos periódicos (eixos k do extracted), com o
!> termo G = 0 referido ao lado do substrato (V -> 0 e campo nulo no fundo; cargas só
!> na região central). Mistura de Anderson em V; parada em max|V_out - V_in| < tol_ev.
!>
!> Unidades: Å e eV fora; a radial em Bohr/Hartree internamente.
!>
!> Dipolos (A6, multipole = 'dipole', sem parâmetro novo): densidade de elétrons de cada sítio
!>   n_i(r) = dq_i rho_s(r - R_i) + dp_i . grad rho_s(r - R_i)      (p = dipolo de carga, elétron = -e)
!> com a mesma rho_s do Gamma; p_i = -s sum_{mu in i, nu} D_{mu nu} <nu|r - R_i|mu> (Mulliken simétrico, termos
!> cruzados, integrais nas PAO). Com Phi(d) o núcleo de Gamma (Ewald 2D + curto alcance, mesma referência):
!>   V_i = sum_j [Phi(d_ij) dq_j + G_ij . dp_j],   F_i = sum_j [G_ij dq_j - T_ij dp_j],
!>   G_ij = grad Phi(d_ij) (eV/Bohr), T_ij = -grad grad Phi(d_ij) (eV/Bohr^2), d_ij = R_i - R_j,
!>   T_ii = (4 pi/3) int rho_s^2 d^3r (= (1/3) int int grad rho_s . grad' rho_s / |r - r'|) + rede.
!> V_i = dE/dq_i e F_i = -dE/dp_i (campo elétrico no sítio, eV/Bohr) da energia
!>   E = 1/2 dq Gamma dq + dq G dp + 1/2 dp T dp  (parte simétrica do núcleo),
!> e o Hamiltoniano variacional é H_{mu nu} += 1/2 (V_i + V_j) S_{mu nu} + 1/2 [F_i . X^(i) + F_j . X^(j)],
!> X^(i) = <mu|r - R_i|nu> (sinal +: energia do elétron num campo F é +F.r).
!>
!> Cargas de ligação (A10, multipole = 'bond', sem parâmetro novo): sítios = átomos centrais e pares b = (i, j, R)
!> de &RHO (todos os pares até rmax, d /= 0), o par no ponto médio, com a mesma forma rho_s:
!>   q_i = s sum_{mu,nu in i} D_{mu nu} S_{nu mu} (só on-site),  p_i = -s sum_{mu,nu in i} D_{mu nu} X^(i)_{nu mu},
!>   q_b = s sum_{mu in i, nu in j} (D_{mu nu} S_{nu mu} + D_{nu mu} S_{mu nu})  (sem dipolo de par);
!> Gamma entre todos os sítios (on-site de par = U_H), carga (átomo ou par)-dipolo (átomo) e dipolo-dipolo como no
!> A6; H variacional: bloco on-site V_i S + F_i . X^(i), bloco do par V_b S (substitui 1/2 (V_i + V_j) S).
!> Referência q_ref/p_ref do volume para átomos e pares (qref_write/qref_file, linhas 'bond').
!> Rotinas sem stop; erros por ierr/errmsg.
module scc
    use types, only: dp
    use pao_basis, only: pao_species_t, read_pao_species, pao_eval, pao_dip_cache_t
    use extracted_green, only: extracted_input_t, run_extracted
    use green_data, only: green_data_t
    use structure_io, only: structure_t, read_structure
    use scc_planar, only: planar_t, planar_setup, planar_charge, planar_dh, planar_mult, planar_pot, planar_read_ref
    implicit none
    private

    public :: scc_input_t, input_has_scc, read_scc_input, run_scc
    public :: radial_density_t, build_radial_density, finish_radial_density, gamma_au, ewald2d_ref, scc_gamma_matrix
    public :: anderson_t, anderson_step
    public :: gamma_derivs_au, onsite_dipole_au, ewald2d_derivs, scc_dipole_kernels, scc_point_field
    public :: scc_cross_kernel, scc_image_kernels, scc_phi_full

    real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp
    real(dp), parameter :: HARTREE_EV = 27.211386245988_dp
    real(dp), parameter :: BOHR_ANG = 0.529177210903_dp
    real(dp), parameter :: KE_EVANG = HARTREE_EV * BOHR_ANG     !< e^2/(4 pi eps0) em eV·Å
    integer, parameter :: MAXK = 1024, PLEN = 512

    !> Espelho de &SCC.
    type :: scc_input_t
        character(len=PLEN) :: pao_file = ''
        character(len=32)   :: basis = 's3p3d2f1'
        real(dp)            :: rcut_bohr = 7.0_dp
        integer             :: gamma_l = 1, gamma_mul = 1
        real(dp)            :: dr_bohr = 5.0e-4_dp
        real(dp)            :: e_lo = -25.0_dp, e_fermi = 0.0_dp      !< eV, na origem de energia do extracted
        integer             :: npoints = 64
        real(dp)            :: rmax = 7.6_dp                          !< Å (pares de Mulliken)
        real(dp)            :: spin_factor = 2.0_dp
        character(len=8)    :: qref_mode = 'value'                    !< 'value' | 'file'
        real(dp)            :: qref = 6.0_dp                          !< elétrons por átomo (com spin)
        character(len=PLEN) :: qref_file = ''
        integer             :: qref_nbasis = 0                        !< átomo de base b = mod(atom-1, nbasis)+1
        character(len=PLEN) :: qref_write = ''                        !< grava q médio por átomo de base
        character(len=PLEN) :: v_init_file = ''                       !< perfil anterior (coluna V) para reinício
        integer             :: k_avg_axis = 0, k_navg = 1, k_list_axis = 3, k_nlist = 0
        real(dp)            :: k_list(MAXK) = 0.0_dp, k_w(MAXK) = 1.0_dp
        real(dp)            :: mix_beta = 0.1_dp
        integer             :: mix_hist = 6, maxiter = 60
        real(dp)            :: tol_ev = 1.0e-4_dp
        real(dp)            :: ewald_alpha = 0.0_dp                   !< 1/Å; 0 = sqrt(pi/A)
        character(len=16)   :: substrate = 'auto'                     !< 'left' | 'right' | 'auto' | 'none' (A7) | 'dielectric' (A14)
        real(dp)            :: eps_sub = 1.0_dp                       !< A14: epsilon do meio abaixo da região central
        logical             :: final_map = .true.
        logical             :: one_shot = .false.                     !< só V_out(V_init), sem ciclo
        character(len=8)    :: multipole = 'monopole'                 !< 'monopole' (A5) | 'dipole' (A6) | 'bond' (A10)
        real(dp)            :: dip_grid_bohr = 0.1_dp                 !< passo da grade das integrais <mu|r|nu>
        real(dp)            :: tol_f = 1.0e-4_dp                      !< eV/Bohr (dipole)
        logical             :: field_out = .false.                    !< A7: campo pontual nos núcleos (dipole)
        real(dp)            :: field_qatom = 6.0_dp                   !< A7: elétrons do átomo neutro (referência)
        logical             :: rho_out = .false.                      !< A8: grava PREFIX_scc.rho_pairs.dat (dipole)
        logical             :: planar_full = .false.                  !< A19: G_par = 0 pela densidade completa
        character(len=PLEN) :: planar_ref_file = ''                   !< A19: rho_pairs.dat do volume (mesma janela)
        character(len=PLEN) :: planar_init_file = ''                  !< A19: PREFIX.scc_planar.dat anterior (dV_in)
        real(dp)            :: planar_dy = 0.1_dp                     !< A19: passo da malha normal (Å)
        real(dp)            :: planar_vac = 8.0_dp                    !< A19: malha até planar_vac acima do topo (Å)
    end type scc_input_t

    !> Densidade esférica tabelada (Bohr, normalizada a 1) e potenciais acumulados.
    type :: radial_density_t
        integer :: n = 0
        real(dp) :: dr = 0.0_dp, rc = 0.0_dp, norm = 0.0_dp
        real(dp), allocatable :: rho(:), w(:), v(:)   !< rho(r_k); V(r_k); W(s) = int_0^s V(t) t dt
        real(dp), allocatable :: qc(:)                !< carga acumulada Q(r_k) (A6)
    end type radial_density_t

    !> Histórico da mistura de Anderson.
    type :: anderson_t
        integer :: m = 0, hist = 6
        real(dp) :: beta = 0.1_dp
        real(dp), allocatable :: x(:, :), f(:, :)   !< (n, hist+1): entradas e resíduos
    end type anderson_t

contains

    logical function input_has_scc(filename) result(yes)
        character(len=*), intent(in) :: filename
        integer :: u, ios
        character(len=512) :: line
        yes = .false.
        open(newunit=u, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) >= 4) then
                if (line(1:1) == '&' .and. (line(2:4) == 'SCC' .or. line(2:4) == 'scc')) then
                    if (len_trim(line) == 4 .or. line(5:5) == ' ') then
                        yes = .true.; exit
                    end if
                end if
            end if
        end do
        close(u)
    end function input_has_scc

    subroutine read_scc_input(filename, si, ierr, errmsg)
        character(len=*), intent(in) :: filename
        type(scc_input_t), intent(out) :: si
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=PLEN) :: pao_file, qref_file, qref_write, v_init_file, planar_ref_file, planar_init_file
        logical :: planar_full
        real(dp) :: planar_dy, planar_vac
        character(len=32) :: basis
        character(len=8) :: qref_mode, multipole
        character(len=16) :: substrate
        real(dp) :: rcut_bohr, dr_bohr, e_lo, e_fermi, rmax, spin_factor, qref, mix_beta, tol_ev, ewald_alpha
        real(dp) :: eps_sub
        real(dp) :: dip_grid_bohr, tol_f, field_qatom
        real(dp) :: k_list(MAXK), k_w(MAXK)
        integer :: gamma_l, gamma_mul, npoints, qref_nbasis, k_avg_axis, k_navg, k_list_axis, k_nlist
        integer :: mix_hist, maxiter, u, ios
        logical :: final_map, one_shot, field_out, rho_out
        character(len=256) :: iomsg
        namelist /SCC/ pao_file, basis, rcut_bohr, gamma_l, gamma_mul, dr_bohr, e_lo, e_fermi, npoints, rmax, &
                       spin_factor, qref_mode, qref, qref_file, qref_nbasis, qref_write, v_init_file, &
                       k_avg_axis, k_navg, k_list_axis, k_nlist, k_list, k_w, mix_beta, mix_hist, maxiter, &
                       tol_ev, ewald_alpha, substrate, final_map, one_shot, multipole, dip_grid_bohr, tol_f, &
                       field_out, field_qatom, rho_out, eps_sub, planar_full, planar_ref_file, planar_init_file, &
                       planar_dy, planar_vac
        ierr = 0
        errmsg = ''
        pao_file = si%pao_file; basis = si%basis; rcut_bohr = si%rcut_bohr; gamma_l = si%gamma_l
        gamma_mul = si%gamma_mul; dr_bohr = si%dr_bohr; e_lo = si%e_lo; e_fermi = si%e_fermi
        npoints = si%npoints; rmax = si%rmax; spin_factor = si%spin_factor; qref_mode = si%qref_mode
        qref = si%qref; qref_file = si%qref_file; qref_nbasis = si%qref_nbasis; qref_write = si%qref_write
        v_init_file = si%v_init_file; k_avg_axis = si%k_avg_axis; k_navg = si%k_navg
        k_list_axis = si%k_list_axis; k_nlist = si%k_nlist; k_list = si%k_list; k_w = si%k_w
        mix_beta = si%mix_beta; mix_hist = si%mix_hist; maxiter = si%maxiter; tol_ev = si%tol_ev
        ewald_alpha = si%ewald_alpha; substrate = si%substrate; final_map = si%final_map; one_shot = si%one_shot
        multipole = si%multipole; dip_grid_bohr = si%dip_grid_bohr; tol_f = si%tol_f
        field_out = si%field_out; field_qatom = si%field_qatom; rho_out = si%rho_out; eps_sub = si%eps_sub
        planar_full = si%planar_full; planar_ref_file = si%planar_ref_file; planar_init_file = si%planar_init_file
        planar_dy = si%planar_dy; planar_vac = si%planar_vac
        open(newunit=u, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1; errmsg = 'scc: cannot open ' // trim(filename); return
        end if
        read(u, nml=SCC, iostat=ios, iomsg=iomsg)
        close(u)
        if (ios /= 0) then
            ierr = 2; errmsg = 'scc: error reading &SCC: ' // trim(iomsg); return
        end if
        si%pao_file = pao_file; si%basis = basis; si%rcut_bohr = rcut_bohr; si%gamma_l = gamma_l
        si%gamma_mul = gamma_mul; si%dr_bohr = dr_bohr; si%e_lo = e_lo; si%e_fermi = e_fermi
        si%npoints = npoints; si%rmax = rmax; si%spin_factor = spin_factor; si%qref_mode = adjustl(qref_mode)
        si%qref = qref; si%qref_file = qref_file; si%qref_nbasis = qref_nbasis; si%qref_write = qref_write
        si%v_init_file = v_init_file; si%k_avg_axis = k_avg_axis; si%k_navg = k_navg
        si%k_list_axis = k_list_axis; si%k_nlist = k_nlist; si%k_list = k_list; si%k_w = k_w
        si%mix_beta = mix_beta; si%mix_hist = mix_hist; si%maxiter = maxiter; si%tol_ev = tol_ev
        si%ewald_alpha = ewald_alpha; si%substrate = adjustl(substrate); si%final_map = final_map
        si%one_shot = one_shot
        si%multipole = adjustl(multipole); si%dip_grid_bohr = dip_grid_bohr; si%tol_f = tol_f
        si%field_out = field_out; si%field_qatom = field_qatom; si%rho_out = rho_out; si%eps_sub = eps_sub
        si%planar_full = planar_full; si%planar_ref_file = planar_ref_file; si%planar_init_file = planar_init_file
        si%planar_dy = planar_dy; si%planar_vac = planar_vac
        if (si%planar_full) then
            if (trim(si%multipole) /= 'dipole' .or. si%field_out) then
                ierr = 2; errmsg = "&SCC planar_full needs multipole = 'dipole' (without field_out)"; return
            end if
            if (trim(si%substrate) == 'dielectric' .or. trim(si%substrate) == 'none') then
                ierr = 2; errmsg = "&SCC planar_full needs substrate = 'auto', 'left' or 'right'"; return
            end if
            if (len_trim(si%planar_ref_file) == 0 .or. si%planar_dy <= 0.0_dp .or. si%planar_vac < 0.0_dp) then
                ierr = 2; errmsg = '&SCC planar_full needs planar_ref_file, planar_dy > 0 and planar_vac >= 0'; return
            end if
            if (trim(si%qref_mode) /= 'file') then
                ierr = 2; errmsg = "&SCC planar_full needs qref_mode = 'file' (qref_nbasis keys the reference pairs)"
                return
            end if
        end if
        if (trim(si%substrate) == 'dielectric') then
            if (si%eps_sub < 1.0_dp) then
                ierr = 2; errmsg = "&SCC substrate = 'dielectric' needs eps_sub >= 1"; return
            end if
            if (trim(si%multipole) == 'bond' .or. si%field_out) then
                ierr = 2; errmsg = "&SCC substrate = 'dielectric' is implemented for multipole = 'monopole'/'dipole' " // &
                                   "without field_out"; return
            end if
        end if
        if (si%field_out .and. trim(si%multipole) /= 'dipole') then
            ierr = 2; errmsg = "&SCC field_out needs multipole = 'dipole'"; return
        end if
        if (trim(si%multipole) /= 'monopole' .and. trim(si%multipole) /= 'dipole' .and. &
            trim(si%multipole) /= 'bond') then
            ierr = 2; errmsg = "&SCC multipole must be 'monopole', 'dipole' or 'bond'"; return
        end if
        if (si%dip_grid_bohr <= 0.0_dp .or. si%tol_f <= 0.0_dp) then
            ierr = 2; errmsg = '&SCC needs dip_grid_bohr > 0 and tol_f > 0'; return
        end if
        if (len_trim(si%pao_file) == 0) then
            ierr = 2; errmsg = '&SCC pao_file is required'; return
        end if
        if (si%k_nlist < 1 .or. si%k_nlist > MAXK .or. si%k_navg < 1) then
            ierr = 2; errmsg = '&SCC needs k_nlist >= 1 (and k_navg >= 1)'; return
        end if
        if (si%k_list_axis < 1 .or. si%k_list_axis > 3 .or. si%k_avg_axis < 0 .or. si%k_avg_axis > 3) then
            ierr = 2; errmsg = '&SCC k_list_axis must be 1..3 and k_avg_axis 0..3'; return
        end if
        if (si%e_lo >= si%e_fermi .or. si%npoints < 2) then
            ierr = 2; errmsg = '&SCC needs e_lo < e_fermi and npoints >= 2'; return
        end if
        if (trim(si%qref_mode) /= 'value' .and. trim(si%qref_mode) /= 'file') then
            ierr = 2; errmsg = "&SCC qref_mode must be 'value' or 'file'"; return
        end if
        if (trim(si%qref_mode) == 'file' .and. (si%qref_nbasis < 1 .or. len_trim(si%qref_file) == 0)) then
            ierr = 2; errmsg = "&SCC qref_mode = 'file' needs qref_file and qref_nbasis"; return
        end if
        if (si%mix_hist < 0 .or. si%mix_beta <= 0.0_dp .or. si%maxiter < 1) then
            ierr = 2; errmsg = '&SCC needs mix_hist >= 0, mix_beta > 0, maxiter >= 1'; return
        end if
    end subroutine read_scc_input

    ! ------------------------------------------------------------------
    ! Gamma: densidade esférica e Coulomb entre duas delas
    ! ------------------------------------------------------------------

    !> rho(r) = sum_m phi_{l,mul,m}(r z)^2 / (2l+1) numa malha uniforme de passo dr (Bohr) até rcut,
    !> normalizada a 1, e W(s) = int_0^s V(t) t dt do potencial de Hartree V da densidade.
    subroutine build_radial_density(sp, l, mul, dr, rd, ierr, errmsg)
        type(pao_species_t), intent(in) :: sp
        integer, intent(in) :: l, mul
        real(dp), intent(in) :: dr
        type(radial_density_t), intent(out) :: rd
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        real(dp), allocatable :: phi(:)
        integer :: o0, k, n, ll, m
        real(dp) :: r
        ierr = 0
        errmsg = ''
        if (l < 0 .or. l > sp%lmax) then
            ierr = 2; errmsg = 'scc: gamma_l outside the basis'; return
        end if
        if (mul < 1 .or. mul > sp%nmul(l)) then
            ierr = 2; errmsg = 'scc: gamma_mul outside the basis'; return
        end if
        o0 = 0
        do ll = 0, l - 1
            o0 = o0 + sp%nmul(ll) * (2 * ll + 1)
        end do
        o0 = o0 + (mul - 1) * (2 * l + 1)
        n = nint(sp%rcut / dr)
        rd%n = n
        rd%dr = sp%rcut / real(n, dp)
        rd%rc = sp%rcut
        allocate(rd%rho(0:n), rd%w(0:n), rd%v(0:n), phi(sp%norb))
        do k = 0, n
            r = real(k, dp) * rd%dr
            ! direção genérica (não alinhada a eixos); a soma em m é invariante
            call pao_eval(sp, r * [0.48_dp, 0.36_dp, 0.8_dp], phi)
            rd%rho(k) = 0.0_dp
            do m = 1, 2 * l + 1
                rd%rho(k) = rd%rho(k) + phi(o0 + m)**2
            end do
            rd%rho(k) = rd%rho(k) / real(2 * l + 1, dp)
        end do
        rd%rho(n) = 0.0_dp
        call finish_radial_density(rd, ierr, errmsg)
    end subroutine build_radial_density

    !> Normaliza rd%rho (já tabelada em 0..n, passo dr, raio rc) e monta V e W.
    subroutine finish_radial_density(rd, ierr, errmsg)
        type(radial_density_t), intent(inout) :: rd
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        real(dp), allocatable :: q(:), p(:), v(:)
        integer :: k, n
        ierr = 0
        errmsg = ''
        n = rd%n
        allocate(q(0:n), p(0:n), v(0:n))
        if (.not. allocated(rd%w)) allocate(rd%w(0:n))
        if (.not. allocated(rd%v)) allocate(rd%v(0:n))
        ! norma e potenciais (trapézio)
        q(0) = 0.0_dp
        do k = 1, n
            q(k) = q(k - 1) + 0.5_dp * rd%dr * 4.0_dp * PI * (rd%rho(k - 1) * (real(k - 1, dp) * rd%dr)**2 + &
                                                            rd%rho(k) * (real(k, dp) * rd%dr)**2)
        end do
        rd%norm = q(n)
        if (rd%norm <= 0.0_dp) then
            ierr = 3; errmsg = 'scc: radial density with zero norm'; return
        end if
        rd%rho = rd%rho / rd%norm
        q = q / rd%norm
        p(n) = 0.0_dp
        do k = n - 1, 0, -1
            p(k) = p(k + 1) + 0.5_dp * rd%dr * 4.0_dp * PI * (rd%rho(k) * real(k, dp) * rd%dr + &
                                                            rd%rho(k + 1) * real(k + 1, dp) * rd%dr)
        end do
        v(0) = p(0)
        do k = 1, n
            v(k) = q(k) / (real(k, dp) * rd%dr) + p(k)
        end do
        rd%w(0) = 0.0_dp
        do k = 1, n
            rd%w(k) = rd%w(k - 1) + 0.5_dp * rd%dr * (v(k - 1) * real(k - 1, dp) * rd%dr + v(k) * real(k, dp) * rd%dr)
        end do
        rd%v = v
        rd%qc = q
    end subroutine finish_radial_density

    !> W(s) com a cauda exata V = 1/s fora da esfera.
    pure real(dp) function wcum(rd, s) result(w)
        type(radial_density_t), intent(in) :: rd
        real(dp), intent(in) :: s
        integer :: k
        real(dp) :: t
        if (s >= rd%rc) then
            w = rd%w(rd%n) + (s - rd%rc)
        else
            t = s / rd%dr
            k = min(int(t), rd%n - 1)
            t = t - real(k, dp)
            w = (1.0_dp - t) * rd%w(k) + t * rd%w(k + 1)
        end if
    end function wcum

    !> Gamma(R) em Hartree para R em Bohr: int rho(r) <V>(r, R) d^3r, <V> = [W(R+r) - W(|R-r|)] / (2 r R).
    pure real(dp) function gamma_au(rd, rr) result(g)
        type(radial_density_t), intent(in) :: rd
        real(dp), intent(in) :: rr
        integer :: k
        real(dp) :: r, f, vk
        if (rr >= 2.0_dp * rd%rc) then
            g = 1.0_dp / rr; return
        end if
        g = 0.0_dp
        do k = 1, rd%n
            r = real(k, dp) * rd%dr
            if (rr < 1.0e-12_dp) then
                vk = rd%v(k)
            else
                vk = (wcum(rd, rr + r) - wcum(rd, abs(rr - r))) / (2.0_dp * r * rr)
            end if
            f = 4.0_dp * PI * rd%rho(k) * r * r * vk
            if (k == rd%n) f = 0.5_dp * f
            g = g + f * rd%dr
        end do
    end function gamma_au

    ! ------------------------------------------------------------------
    ! Ewald 2D
    ! ------------------------------------------------------------------

    !> sum_L 1/|d + L| (Å^-1) na rede 2D {a1, a2} (L = 0 excluído se self), com o termo G = 0
    !> referido a z -> -inf: phi_0 = -(2 pi / A) [z erf(alpha z) + e^{-alpha^2 z^2}/(alpha sqrt(pi)) + z],
    !> z = d . nhat. Para cargas acima de i (z < 0) o termo some longe; abaixo, -(4 pi/A) z.
    !> sym (A7, opcional, padrão .false.): sem substrato, termo G = 0 simétrico
    !> phi_0 = -(2 pi / A) [z erf(alpha z) + e^{-alpha^2 z^2}/(alpha sqrt(pi))] (potencial a menos de constante;
    !> para slab neutro o campo é o físico, nulo dos dois lados longe).
    real(dp) function ewald2d_ref(d, a1, a2, nhat, alpha, self, sym) result(phi)
        real(dp), intent(in) :: d(3), a1(3), a2(3), nhat(3), alpha
        logical, intent(in) :: self
        logical, intent(in), optional :: sym
        real(dp), parameter :: TOLX = 6.0_dp      ! erfc(6) ~ 2e-17
        real(dp) :: area, z, dp_(3), b1(3), b2(3), g(3), gn, rv(3), r, c(3), x
        integer :: n1, n2, m1, m2, nmax1, nmax2
        c = cross(a1, a2)
        area = norm2(c)
        z = dot_product(d, nhat)
        dp_ = d - z * nhat
        phi = 0.0_dp
        ! espaço real
        nmax1 = ceiling((TOLX / alpha + abs(z)) / (area / norm2(a2))) + 1
        nmax2 = ceiling((TOLX / alpha + abs(z)) / (area / norm2(a1))) + 1
        do n1 = -nmax1, nmax1
            do n2 = -nmax2, nmax2
                if (self .and. n1 == 0 .and. n2 == 0) cycle
                rv = d + real(n1, dp) * a1 + real(n2, dp) * a2
                r = norm2(rv)
                if (alpha * r > TOLX) cycle
                phi = phi + erfc(alpha * r) / r
            end do
        end do
        if (self) phi = phi - 2.0_dp * alpha / sqrt(PI)
        ! recíproco (G /= 0)
        b1 = 2.0_dp * PI * cross(a2, nhat) / dot_product(a1, cross(a2, nhat))
        b2 = 2.0_dp * PI * cross(nhat, a1) / dot_product(a2, cross(nhat, a1))
        nmax1 = ceiling(2.0_dp * alpha * TOLX / (2.0_dp * PI / norm2(a1))) + 2
        nmax2 = ceiling(2.0_dp * alpha * TOLX / (2.0_dp * PI / norm2(a2))) + 2
        do m1 = -nmax1, nmax1
            do m2 = -nmax2, nmax2
                if (m1 == 0 .and. m2 == 0) cycle
                g = real(m1, dp) * b1 + real(m2, dp) * b2
                gn = norm2(g)
                if (gn / (2.0_dp * alpha) > TOLX + alpha * abs(z) + 1.0_dp) cycle
                phi = phi + (PI / area) * cos(dot_product(g, dp_)) / gn * (ef(gn, z) + ef(gn, -z))
            end do
        end do
        ! G = 0 referido ao fundo (padrão) ou simétrico (sym)
        x = alpha * z
        if (present(sym)) then
            if (sym) then
                phi = phi - (2.0_dp * PI / area) * (z * erf(x) + exp(-x * x) / (alpha * sqrt(PI)))
                return
            end if
        end if
        phi = phi - (2.0_dp * PI / area) * (z * erf(x) + exp(-x * x) / (alpha * sqrt(PI)) + z)
    contains
        pure real(dp) function ef(gg, s)
            real(dp), intent(in) :: gg, s
            real(dp) :: xx
            xx = gg / (2.0_dp * alpha) + alpha * s
            if (xx > 0.0_dp) then
                ef = exp(-gg * gg / (4.0_dp * alpha * alpha) - alpha * alpha * s * s) * erfc_scaled(xx)
            else
                ef = exp(gg * s) * erfc(xx)
            end if
        end function ef
    end function ewald2d_ref

    !> A10: d = 0 módulo a rede 2D {a1, a2} (tolerância 1e-5 Å): dois sítios no mesmo ponto (ou imagens).
    pure logical function coincident(d, a1, a2, nhat) result(yes)
        real(dp), intent(in) :: d(3), a1(3), a2(3), nhat(3)
        real(dp), parameter :: TOL = 1.0e-5_dp
        real(dp) :: dp_(3), g11, g12, g22, b1, b2, det, c1, c2
        yes = .false.
        if (abs(dot_product(d, nhat)) > TOL) return
        dp_ = d - dot_product(d, nhat) * nhat
        g11 = dot_product(a1, a1); g12 = dot_product(a1, a2); g22 = dot_product(a2, a2)
        b1 = dot_product(dp_, a1); b2 = dot_product(dp_, a2)
        det = g11 * g22 - g12 * g12
        c1 = anint((g22 * b1 - g12 * b2) / det)
        c2 = anint((g11 * b2 - g12 * b1) / det)
        yes = norm2(dp_ - c1 * a1 - c2 * a2) < TOL
    end function coincident

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
    end function cross

    !> Gamma_ij (eV) dos átomos em pos(3, n) (Å), rede 2D {a1, a2}, normal nhat (para o vácuo):
    !> Ewald de 1/r (fundo como referência) + sum_L [gamma(|d+L|) - 1/|d+L|] a curto alcance + U_H no on-site.
    subroutine scc_gamma_matrix(rd, pos, a1, a2, nhat, alpha, gam, uh_ev, sym)
        type(radial_density_t), intent(in) :: rd
        real(dp), intent(in) :: pos(:, :), a1(3), a2(3), nhat(3), alpha
        real(dp), intent(out) :: gam(:, :), uh_ev
        logical, intent(in), optional :: sym
        logical :: sy, same
        integer :: i, j, n, n1, n2, nm1, nm2
        real(dp) :: d(3), rv(3), r, rcut_ang, area, corr
        n = size(pos, 2)
        sy = .false.
        if (present(sym)) sy = sym
        uh_ev = HARTREE_EV * gamma_au(rd, 0.0_dp)
        rcut_ang = 2.0_dp * rd%rc * BOHR_ANG
        area = norm2(cross(a1, a2))
        nm1 = ceiling(rcut_ang / (area / norm2(a2))) + 1
        nm2 = ceiling(rcut_ang / (area / norm2(a1))) + 1
        !$omp parallel do default(shared) private(i, j, d, rv, r, corr, n1, n2, same) schedule(dynamic)
        do i = 1, n
            do j = 1, n
                d = pos(:, i) - pos(:, j)
                same = i == j
                if (.not. same) same = coincident(d, a1, a2, nhat)   ! A10: sítios distintos no mesmo ponto (mod rede)
                if (same) d = 0.0_dp
                gam(i, j) = KE_EVANG * ewald2d_ref(d, a1, a2, nhat, alpha, same, sy)
                corr = 0.0_dp
                do n1 = -nm1, nm1
                    do n2 = -nm2, nm2
                        rv = d + real(n1, dp) * a1 + real(n2, dp) * a2
                        r = norm2(rv)
                        if (r < 1.0e-8_dp .or. r >= rcut_ang) cycle
                        corr = corr + HARTREE_EV * gamma_au(rd, r / BOHR_ANG) - KE_EVANG / r
                    end do
                end do
                gam(i, j) = gam(i, j) + corr
                if (same) gam(i, j) = gam(i, j) + uh_ev
            end do
        end do
        !$omp end parallel do
    end subroutine scc_gamma_matrix

    ! ------------------------------------------------------------------
    ! Dipolos (A6): derivadas de Gamma e núcleos carga-dipolo e dipolo-dipolo
    ! ------------------------------------------------------------------

    !> V(s) e Q(s) da densidade (interpolação linear; cauda exata 1/s e 1 fora da esfera).
    pure subroutine vq_at(rd, s, v, q)
        type(radial_density_t), intent(in) :: rd
        real(dp), intent(in) :: s
        real(dp), intent(out) :: v, q
        integer :: k
        real(dp) :: t
        if (s >= rd%rc) then
            v = 1.0_dp / s; q = 1.0_dp
        else
            t = s / rd%dr
            k = min(int(t), rd%n - 1)
            t = t - real(k, dp)
            v = (1.0_dp - t) * rd%v(k) + t * rd%v(k + 1)
            q = (1.0_dp - t) * rd%qc(k) + t * rd%qc(k + 1)
        end if
    end subroutine vq_at

    !> Derivadas radiais g' (Hartree/Bohr) e g'' (Hartree/Bohr^2) de Gamma(R), R em Bohr (R > 0),
    !> pela derivada analítica de <V>(r, R) = [W(R+r) - W(|R-r|)] / (2 r R), com W' = s V, W'' = V - Q/s.
    pure subroutine gamma_derivs_au(rd, rr, g1, g2)
        type(radial_density_t), intent(in) :: rd
        real(dp), intent(in) :: rr
        real(dp), intent(out) :: g1, g2
        integer :: k
        real(dp) :: r, sp, sm, vp, qp, vm, qm, a0, a1, a2, w2p, w2m, f
        if (rr >= 2.0_dp * rd%rc) then
            g1 = -1.0_dp / rr**2; g2 = 2.0_dp / rr**3; return
        end if
        g1 = 0.0_dp
        g2 = 0.0_dp
        do k = 1, rd%n
            r = real(k, dp) * rd%dr
            sp = rr + r
            sm = abs(rr - r)
            call vq_at(rd, sp, vp, qp)
            call vq_at(rd, sm, vm, qm)
            a0 = wcum(rd, sp) - wcum(rd, sm)
            a1 = sp * vp - (rr - r) * vm
            w2p = vp - qp / sp
            if (sm > 1.0e-12_dp) then
                w2m = vm - qm / sm
            else
                w2m = rd%v(0)
            end if
            a2 = w2p - w2m
            f = 4.0_dp * PI * rd%rho(k) * r * r * rd%dr
            if (k == rd%n) f = 0.5_dp * f
            g1 = g1 + f * (a1 / (2.0_dp * r * rr) - a0 / (2.0_dp * r * rr * rr))
            g2 = g2 + f * (a2 / (2.0_dp * r * rr) - a1 / (r * rr * rr) + a0 / (r * rr**3))
        end do
    end subroutine gamma_derivs_au

    !> T on-site (Hartree/Bohr^2): -grad grad Gamma(0) = (4 pi / 3) int rho_s^2 d^3r.
    pure real(dp) function onsite_dipole_au(rd) result(t0)
        type(radial_density_t), intent(in) :: rd
        integer :: k
        real(dp) :: r, f
        t0 = 0.0_dp
        do k = 1, rd%n
            r = real(k, dp) * rd%dr
            f = 4.0_dp * PI * r * r * rd%rho(k)**2
            if (k == rd%n) f = 0.5_dp * f
            t0 = t0 + f * rd%dr
        end do
        t0 = 4.0_dp * PI / 3.0_dp * t0
    end function onsite_dipole_au

    !> Gradiente (Å^-2) e hessiana (Å^-3) de ewald2d_ref em d, por diferenças centrais com extrapolação de
    !> Richardson (passos h e h/2, h = 0,02 Å; erro O(h^4)). self: a soma exclui L = 0 (o termo liso
    !> erf(alpha r)/r - 2 alpha/sqrt(pi) que o Ewald deixa no lugar dele é descontado da hessiana).
    subroutine ewald2d_derivs(d, a1, a2, nhat, alpha, self, grad, hess, sym)
        real(dp), intent(in) :: d(3), a1(3), a2(3), nhat(3), alpha
        logical, intent(in) :: self
        real(dp), intent(out) :: grad(3), hess(3, 3)
        logical, intent(in), optional :: sym
        real(dp), parameter :: H0 = 0.02_dp
        real(dp) :: gh(3, 2), hh(3, 3, 2), f0, fp, fm, fpp, fpm, fmp, fmm, h, ea(3), eb(3)
        integer :: lev, a, b
        f0 = ewald2d_ref(d, a1, a2, nhat, alpha, self, sym)
        do lev = 1, 2
            h = H0 / real(lev, dp)
            do a = 1, 3
                ea = 0.0_dp; ea(a) = h
                fp = ewald2d_ref(d + ea, a1, a2, nhat, alpha, self, sym)
                fm = ewald2d_ref(d - ea, a1, a2, nhat, alpha, self, sym)
                gh(a, lev) = (fp - fm) / (2.0_dp * h)
                hh(a, a, lev) = (fp - 2.0_dp * f0 + fm) / (h * h)
                do b = a + 1, 3
                    eb = 0.0_dp; eb(b) = h
                    fpp = ewald2d_ref(d + ea + eb, a1, a2, nhat, alpha, self, sym)
                    fpm = ewald2d_ref(d + ea - eb, a1, a2, nhat, alpha, self, sym)
                    fmp = ewald2d_ref(d - ea + eb, a1, a2, nhat, alpha, self, sym)
                    fmm = ewald2d_ref(d - ea - eb, a1, a2, nhat, alpha, self, sym)
                    hh(a, b, lev) = (fpp - fpm - fmp + fmm) / (4.0_dp * h * h)
                    hh(b, a, lev) = hh(a, b, lev)
                end do
            end do
        end do
        grad = (4.0_dp * gh(:, 2) - gh(:, 1)) / 3.0_dp
        hess = (4.0_dp * hh(:, :, 2) - hh(:, :, 1)) / 3.0_dp
        if (self) then
            do a = 1, 3
                hess(a, a) = hess(a, a) + 4.0_dp * alpha**3 / (3.0_dp * sqrt(PI))
            end do
        end if
    end subroutine ewald2d_derivs

    !> Núcleos dos dipolos para os átomos em pos(3, n) (Å), mesma rede/referência de scc_gamma_matrix:
    !>   gk(:, i, j) = grad Phi(d_ij) (eV/Bohr),  tk(:, :, i, j) = -grad grad Phi(d_ij) (eV/Bohr^2),
    !> Phi = Ewald de 1/r + sum_L [gamma - 1/r] (curto alcance); t0_ev = T on-site (eV/Bohr^2), já em tk(:,:,i,i).
    subroutine scc_dipole_kernels(rd, pos, a1, a2, nhat, alpha, gk, tk, t0_ev, sym)
        type(radial_density_t), intent(in) :: rd
        real(dp), intent(in) :: pos(:, :), a1(3), a2(3), nhat(3), alpha
        real(dp), intent(out) :: gk(:, :, :), tk(:, :, :, :), t0_ev
        logical, intent(in), optional :: sym
        logical :: sy
        integer :: i, j, n, n1, n2, nm1, nm2, a, b
        real(dp) :: d(3), rv(3), r, rb, rcut_ang, area, ge(3), he(3, 3), g1, g2, f1, f2, u(3), gb(3), hb(3, 3)
        n = size(pos, 2)
        sy = .false.
        if (present(sym)) sy = sym
        t0_ev = HARTREE_EV * onsite_dipole_au(rd)
        rcut_ang = 2.0_dp * rd%rc * BOHR_ANG
        area = norm2(cross(a1, a2))
        nm1 = ceiling(rcut_ang / (area / norm2(a2))) + 1
        nm2 = ceiling(rcut_ang / (area / norm2(a1))) + 1
        !$omp parallel do default(shared) private(i, j, d, rv, r, rb, ge, he, g1, g2, f1, f2, u, gb, hb, n1, n2, a, b) &
        !$omp schedule(dynamic)
        do i = 1, n
            do j = 1, n
                d = pos(:, i) - pos(:, j)
                call ewald2d_derivs(d, a1, a2, nhat, alpha, i == j, ge, he, sy)
                gb = 0.0_dp
                hb = 0.0_dp
                do n1 = -nm1, nm1
                    do n2 = -nm2, nm2
                        rv = d + real(n1, dp) * a1 + real(n2, dp) * a2
                        r = norm2(rv)
                        if (r < 1.0e-8_dp .or. r >= rcut_ang) cycle
                        rb = r / BOHR_ANG
                        call gamma_derivs_au(rd, rb, g1, g2)
                        f1 = g1 + 1.0_dp / rb**2
                        f2 = g2 - 2.0_dp / rb**3
                        u = rv / r
                        gb = gb + f1 * u
                        do b = 1, 3
                            do a = 1, 3
                                hb(a, b) = hb(a, b) + f2 * u(a) * u(b)
                            end do
                            hb(b, b) = hb(b, b) + f1 / rb
                        end do
                        do b = 1, 3
                            do a = 1, 3
                                hb(a, b) = hb(a, b) - f1 / rb * u(a) * u(b)
                            end do
                        end do
                    end do
                end do
                gk(:, i, j) = KE_EVANG * BOHR_ANG * ge + HARTREE_EV * gb
                tk(:, :, i, j) = -(KE_EVANG * BOHR_ANG**2 * he + HARTREE_EV * hb)
                if (i == j) then
                    do a = 1, 3
                        tk(a, a, i, j) = tk(a, a, i, j) + t0_ev
                    end do
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine scc_dipole_kernels

    ! ------------------------------------------------------------------
    ! A14: substrato dielétrico (imagens eletrostáticas)
    ! ------------------------------------------------------------------

    !> Phi(d) completo (eV) = Ewald 2D de 1/r (referência no fundo, sem L = 0 excluído) + sum_L [gamma - 1/r]
    !> a curto alcance, e (se want_d) grad Phi (eV/Bohr) e hess Phi (eV/Bohr^2), d em Å. É o núcleo de
    !> scc_gamma_matrix / scc_dipole_kernels para um par de sítios distintos (tk = -hess).
    subroutine scc_phi_full(rd, d, a1, a2, nhat, alpha, want_d, phi, grad, hess)
        type(radial_density_t), intent(in) :: rd
        real(dp), intent(in) :: d(3), a1(3), a2(3), nhat(3), alpha
        logical, intent(in) :: want_d
        real(dp), intent(out) :: phi, grad(3), hess(3, 3)
        integer :: n1, n2, nm1, nm2, a, b
        real(dp) :: rv(3), r, rb, rcut_ang, area, ge(3), he(3, 3), g1, g2, f1, f2, u(3), gb(3), hb(3, 3), corr
        rcut_ang = 2.0_dp * rd%rc * BOHR_ANG
        area = norm2(cross(a1, a2))
        nm1 = ceiling(rcut_ang / (area / norm2(a2))) + 1
        nm2 = ceiling(rcut_ang / (area / norm2(a1))) + 1
        phi = KE_EVANG * ewald2d_ref(d, a1, a2, nhat, alpha, .false.)
        corr = 0.0_dp
        gb = 0.0_dp
        hb = 0.0_dp
        do n1 = -nm1, nm1
            do n2 = -nm2, nm2
                rv = d + real(n1, dp) * a1 + real(n2, dp) * a2
                r = norm2(rv)
                if (r < 1.0e-8_dp .or. r >= rcut_ang) cycle
                corr = corr + HARTREE_EV * gamma_au(rd, r / BOHR_ANG) - KE_EVANG / r
                if (.not. want_d) cycle
                rb = r / BOHR_ANG
                call gamma_derivs_au(rd, rb, g1, g2)
                f1 = g1 + 1.0_dp / rb**2
                f2 = g2 - 2.0_dp / rb**3
                u = rv / r
                gb = gb + f1 * u
                do b = 1, 3
                    do a = 1, 3
                        hb(a, b) = hb(a, b) + (f2 - f1 / rb) * u(a) * u(b)
                    end do
                    hb(b, b) = hb(b, b) + f1 / rb
                end do
            end do
        end do
        phi = phi + corr
        grad = 0.0_dp
        hess = 0.0_dp
        if (.not. want_d) return
        call ewald2d_derivs(d, a1, a2, nhat, alpha, .false., ge, he)
        grad = KE_EVANG * BOHR_ANG * ge + HARTREE_EV * gb
        hess = KE_EVANG * BOHR_ANG**2 * he + HARTREE_EV * hb
    end subroutine scc_phi_full

    !> A14: soma aos núcleos o meio dielétrico semi-infinito (epsilon) abaixo do plano h = zplane (h = r . nhat,
    !> nhat para o vácuo). Cada fonte j (carga dq_j com a forma rho_s, dipolo dp_j) em R_j tem a imagem em
    !> M R_j (M R = R - 2 (R . n - zplane) n, parte linear Ml = I - 2 n n^T) com carga -beta dq_j e dipolo
    !> -beta Ml dp_j, beta = (eps - 1)/(eps + 1) (componentes paralelas trocam de sinal; a normal não). Com
    !> K(R_i, R_j) = -beta Phi(d'), d' = R_i - M R_j (Phi de scc_phi_full, a mesma referência G = 0 do 'auto'):
    !>   gam_ij += -beta Phi(d'),  gk_ij (V de dp_j) += -beta Ml grad Phi(d'),
    !>   gkf_ij (F de dq_j) += -beta grad Phi(d'),  tk_ij (F = -tk dp_j) += beta hess Phi(d') Ml.
    !> Na forma direta gk = gkf; com a imagem diferem (V_i = sum_j [gam dq + gk . dp], F_i = sum_j [gkf dq - tk dp]).
    !> Válido para fontes acima do plano; ponto e imagem nunca coincidem (|d' . n| = h_i + h_j - 2 zplane > 0).
    subroutine scc_image_kernels(rd, pos, a1, a2, nhat, alpha, zplane, beta, gam, gk, gkf, tk)
        type(radial_density_t), intent(in) :: rd
        real(dp), intent(in) :: pos(:, :), a1(3), a2(3), nhat(3), alpha, zplane, beta
        real(dp), intent(inout) :: gam(:, :)
        real(dp), intent(inout), optional :: gk(:, :, :), gkf(:, :, :), tk(:, :, :, :)
        integer :: i, j, n, a
        real(dp) :: ml(3, 3), rm(3), d(3), phi, g(3), h(3, 3)
        logical :: wd
        n = size(pos, 2)
        wd = present(gk) .and. present(gkf) .and. present(tk)
        ml = 0.0_dp
        do a = 1, 3
            ml(a, a) = 1.0_dp
            ml(:, a) = ml(:, a) - 2.0_dp * nhat * nhat(a)
        end do
        !$omp parallel do default(shared) private(i, j, rm, d, phi, g, h) schedule(dynamic)
        do i = 1, n
            do j = 1, n
                rm = pos(:, j) - 2.0_dp * (dot_product(pos(:, j), nhat) - zplane) * nhat
                d = pos(:, i) - rm
                call scc_phi_full(rd, d, a1, a2, nhat, alpha, wd, phi, g, h)
                gam(i, j) = gam(i, j) - beta * phi
                if (.not. wd) cycle
                gk(:, i, j) = gk(:, i, j) - beta * matmul(ml, g)
                gkf(:, i, j) = gkf(:, i, j) - beta * g
                tk(:, :, i, j) = tk(:, :, i, j) + beta * matmul(h, ml)
            end do
        end do
        !$omp end parallel do
    end subroutine scc_image_kernels

    !> A10: núcleo carga-dipolo entre sítios spos(3, ns) e átomos apos(3, na) (Å), mesma rede/referência:
    !>   gk(:, s, a) = grad Phi(R_s - R_a) (eV/Bohr), como gk de scc_dipole_kernels (só o gradiente;
    !> diferenças centrais com Richardson, h = 0,02 Å). Sítio no mesmo ponto do átomo: soma sem L = 0.
    subroutine scc_cross_kernel(rd, spos, apos, a1, a2, nhat, alpha, gk, sym)
        type(radial_density_t), intent(in) :: rd
        real(dp), intent(in) :: spos(:, :), apos(:, :), a1(3), a2(3), nhat(3), alpha
        real(dp), intent(out) :: gk(:, :, :)
        logical, intent(in), optional :: sym
        real(dp), parameter :: H0 = 0.02_dp
        logical :: sy, same
        integer :: s, a, n1, n2, nm1, nm2, lev, c
        real(dp) :: d(3), rv(3), r, rb, rcut_ang, area, gh(3, 2), ea(3), h, g1, g2, gb(3)
        sy = .false.
        if (present(sym)) sy = sym
        rcut_ang = 2.0_dp * rd%rc * BOHR_ANG
        area = norm2(cross(a1, a2))
        nm1 = ceiling(rcut_ang / (area / norm2(a2))) + 1
        nm2 = ceiling(rcut_ang / (area / norm2(a1))) + 1
        !$omp parallel do default(shared) private(s, a, d, same, lev, h, c, ea, gh, gb, n1, n2, rv, r, rb, g1, g2) &
        !$omp schedule(dynamic) collapse(2)
        do a = 1, size(apos, 2)
            do s = 1, size(spos, 2)
                d = spos(:, s) - apos(:, a)
                same = coincident(d, a1, a2, nhat)
                if (same) d = 0.0_dp
                do lev = 1, 2
                    h = H0 / real(lev, dp)
                    do c = 1, 3
                        ea = 0.0_dp; ea(c) = h
                        gh(c, lev) = (ewald2d_ref(d + ea, a1, a2, nhat, alpha, same, sy) - &
                                      ewald2d_ref(d - ea, a1, a2, nhat, alpha, same, sy)) / (2.0_dp * h)
                    end do
                end do
                gb = 0.0_dp
                do n1 = -nm1, nm1
                    do n2 = -nm2, nm2
                        rv = d + real(n1, dp) * a1 + real(n2, dp) * a2
                        r = norm2(rv)
                        if (r < 1.0e-8_dp .or. r >= rcut_ang) cycle
                        rb = r / BOHR_ANG
                        call gamma_derivs_au(rd, rb, g1, g2)
                        gb = gb + (g1 + 1.0_dp / rb**2) * rv / r
                    end do
                end do
                gk(:, s, a) = KE_EVANG * BOHR_ANG * (4.0_dp * gh(:, 2) - gh(:, 1)) / 3.0_dp + HARTREE_EV * gb
            end do
        end do
        !$omp end parallel do
    end subroutine scc_cross_kernel

    !> Campo elétrico pontual (A7) nos núcleos pos(3, n) (Å) de cargas/dipolos com a forma rho_s do Gamma:
    !>   U(r) = sum_{j,L} [qn_j v_s(r - R_j - L) + p_j . grad v_s(r - R_j - L)]   (energia potencial do elétron),
    !>   E(R_i) = grad U(R_i)  (= -grad phi, campo elétrico usual; V/Å),
    !> qn_j = excesso de elétrons (e), p_j = dipolo de carga (e Bohr, elétron = -e), v_s = potencial de rho_s
    !> (1/r fora de rc). Longo alcance pelo Ewald 2D de 1/r (mesma referência de G = 0 de scc_gamma_matrix,
    !> sym opcional), curto alcance pela diferença v_s - 1/r dentro de rc. Termo próprio (j = i, L = 0):
    !> grad v_s(0) = 0 e grad grad v_s(0) = -(4 pi/3) rho_s(0) I, devolvido à parte em eself (V/Å; já em efield).
    subroutine scc_point_field(rd, pos, a1, a2, nhat, alpha, qn, p, efield, eself, sym)
        type(radial_density_t), intent(in) :: rd
        real(dp), intent(in) :: pos(:, :), a1(3), a2(3), nhat(3), alpha, qn(:), p(:, :)
        real(dp), intent(out) :: efield(:, :), eself(:, :)
        logical, intent(in), optional :: sym
        logical :: sy
        integer :: i, j, n, n1, n2, nm1, nm2, a
        real(dp) :: d(3), rv(3), r, rb, rc_ang, area, ge(3), he(3, 3), v, q, f1, f2, rho, u(3), hb(3, 3), t, &
                    eau(3), pa(3)
        integer :: k
        n = size(pos, 2)
        sy = .false.
        if (present(sym)) sy = sym
        rc_ang = rd%rc * BOHR_ANG
        area = norm2(cross(a1, a2))
        nm1 = ceiling(rc_ang / (area / norm2(a2))) + 1
        nm2 = ceiling(rc_ang / (area / norm2(a1))) + 1
        !$omp parallel do default(shared) private(i, j, d, rv, r, rb, ge, he, v, q, f1, f2, rho, u, hb, t, eau, pa, &
        !$omp n1, n2, a, k) schedule(dynamic)
        do i = 1, n
            efield(:, i) = 0.0_dp
            eau = 0.0_dp                                    ! Hartree/Bohr
            do j = 1, n
                d = pos(:, i) - pos(:, j)
                pa = p(:, j) * BOHR_ANG                     ! e Å
                call ewald2d_derivs(d, a1, a2, nhat, alpha, i == j, ge, he, sy)
                ! Ewald (e^2 = KE eV Å): campo em V/Å
                efield(:, i) = efield(:, i) + KE_EVANG * (qn(j) * ge + matmul(he, pa))
                ! curto alcance: imagens com |d + L| < rc (para j = i, L = 0 fica no termo próprio)
                do n1 = -nm1, nm1
                    do n2 = -nm2, nm2
                        rv = d + real(n1, dp) * a1 + real(n2, dp) * a2
                        r = norm2(rv)
                        if (r < 1.0e-8_dp .or. r >= rc_ang) cycle
                        rb = r / BOHR_ANG
                        call vq_at(rd, rb, v, q)
                        t = rb / rd%dr
                        k = min(int(t), rd%n - 1)
                        t = t - real(k, dp)
                        rho = (1.0_dp - t) * rd%rho(k) + t * rd%rho(k + 1)
                        f1 = (1.0_dp - q) / rb**2                       ! v_s' - (1/r)'
                        f2 = 2.0_dp * (q - 1.0_dp) / rb**3 - 4.0_dp * PI * rho   ! v_s'' - (1/r)''
                        u = rv / r
                        do a = 1, 3
                            hb(:, a) = (f2 - f1 / rb) * u * u(a)
                            hb(a, a) = hb(a, a) + f1 / rb
                        end do
                        eau = eau + qn(j) * f1 * u + matmul(hb, p(:, j))
                    end do
                end do
            end do
            eself(:, i) = -(4.0_dp * PI / 3.0_dp) * rd%rho(0) * p(:, i) * HARTREE_EV / BOHR_ANG
            efield(:, i) = efield(:, i) + eau * HARTREE_EV / BOHR_ANG + eself(:, i)
        end do
        !$omp end parallel do
    end subroutine scc_point_field

    ! ------------------------------------------------------------------
    ! Mistura de Anderson
    ! ------------------------------------------------------------------

    !> Dado x (entrada) e f = g(x) - x, devolve a próxima entrada em x.
    subroutine anderson_step(st, x, f)
        type(anderson_t), intent(inout) :: st
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: f(:)
        integer :: n, m, i, j, k, piv
        real(dp), allocatable :: dx(:, :), df(:, :), a(:, :), b(:), th(:), xb(:), fb(:)
        real(dp) :: t
        n = size(x)
        if (.not. allocated(st%x)) then
            allocate(st%x(n, st%hist + 1), st%f(n, st%hist + 1))
            st%m = 0
        end if
        ! empilha (x, f) no histórico, descartando o mais antigo
        if (st%m == st%hist + 1) then
            st%x(:, 1:st%hist) = st%x(:, 2:st%hist + 1)
            st%f(:, 1:st%hist) = st%f(:, 2:st%hist + 1)
            st%m = st%hist
        end if
        st%m = st%m + 1
        st%x(:, st%m) = x
        st%f(:, st%m) = f
        m = st%m - 1
        allocate(xb(n), fb(n))
        xb = x
        fb = f
        if (m > 0) then
            allocate(dx(n, m), df(n, m), a(m, m), b(m), th(m))
            do k = 1, m
                dx(:, k) = st%x(:, st%m) - st%x(:, st%m - k)
                df(:, k) = st%f(:, st%m) - st%f(:, st%m - k)
            end do
            a = matmul(transpose(df), df)
            b = matmul(transpose(df), f)
            do k = 1, m
                a(k, k) = a(k, k) * (1.0_dp + 1.0e-10_dp) + 1.0e-30_dp
            end do
            ! eliminação de Gauss com pivoteamento parcial
            do k = 1, m
                piv = k - 1 + maxloc(abs(a(k:m, k)), 1)
                if (piv /= k) then
                    a([k, piv], :) = a([piv, k], :)
                    b([k, piv]) = b([piv, k])
                end if
                do i = k + 1, m
                    t = a(i, k) / a(k, k)
                    a(i, k:m) = a(i, k:m) - t * a(k, k:m)
                    b(i) = b(i) - t * b(k)
                end do
            end do
            do i = m, 1, -1
                th(i) = b(i)
                do j = i + 1, m
                    th(i) = th(i) - a(i, j) * th(j)
                end do
                th(i) = th(i) / a(i, i)
            end do
            xb = x - matmul(dx, th)
            fb = f - matmul(df, th)
        end if
        x = xb + st%beta * fb
    end subroutine anderson_step

    ! ------------------------------------------------------------------
    ! Ciclo
    ! ------------------------------------------------------------------

    !> Executa o ciclo autoconsistente e, se final_map, o mapa com o V convergido.
    !> Saídas: PREFIX.scc_log.dat (uma linha por iteração), PREFIX.scc_profile.dat (perfil final)
    !> e PREFIX.green.bin (via gd, gravado por quem chama).
    subroutine run_scc(inp, si, gd, ierr, errmsg, notes, converged)
        type(extracted_input_t), intent(in) :: inp
        type(scc_input_t), intent(in) :: si
        type(green_data_t), intent(out) :: gd
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: errmsg
        character(len=256), allocatable, intent(out) :: notes(:)
        logical, intent(out) :: converged
        type(extracted_input_t) :: ic
        type(pao_species_t) :: pao
        type(radial_density_t) :: rd
        type(structure_t) :: st
        type(green_data_t) :: gtmp
        type(anderson_t) :: mix
        real(dp), allocatable :: pos(:, :), gam(:, :), q(:, :), qref(:), dq(:), vin(:), vout(:), f(:), qb(:)
        integer, allocatable :: col(:), atom(:), cntb(:)
        ! dipolos (A6)
        type(pao_dip_cache_t) :: dc
        real(dp), allocatable :: gk(:, :, :), tk(:, :, :, :), pp(:, :, :), pref(:, :), dpv(:, :), fin(:, :), fout(:, :)
        real(dp), allocatable :: gkf(:, :, :)      ! A14: núcleo do campo de dq (= gk sem imagens)
        real(dp) :: zpl, hmin, hbel, beta
        logical :: diel
        real(dp), allocatable :: xv(:), xf(:)
        real(dp) :: t0d, dvmax, dfmax
        logical :: dip, sy
        character(len=256), allocatable :: nt(:)
        character(len=:), allocatable :: perr
        character(len=512) :: msg
        integer :: n, it, i, ul, a, kax(2), na, side, nbad
        real(dp) :: a1(3), a2(3), nhat(3), alpha, uh, t0, t1, ae(3)
        logical :: sl
        ! cargas de ligação (A10)
        logical :: bnd
        integer :: nbd, ns, nsv, ngh, nall
        real(dp), allocatable :: gpos(:, :), gdq(:)
        integer, allocatable :: bij(:, :)
        real(dp), allocatable :: bdd(:, :), bpos(:, :), spos(:, :), gks(:, :, :), gka(:, :, :), qon(:, :), qbd(:, :), &
                                 pon(:, :, :)
        ! componente planar pela densidade completa (A19)
        logical :: pla
        integer :: npl
        type(planar_t) :: pl
        integer, allocatable :: pic(:), pja(:), pkind(:), ric(:), nrow(:)
        real(dp), allocatable :: rvec(:, :)
        real(dp), allocatable :: pvec(:, :), pss(:, :, :), dmat(:, :, :), dref(:, :, :), dhp(:, :, :), cfull(:), &
                                 cref(:), cmul(:), vfull(:), vmul(:), pvin(:), pvout(:), dpn(:)
        real(dp) :: dvp
        integer :: upl

        ierr = 0
        errmsg = ''
        converged = .false.
        allocate(notes(0))
        dip = trim(si%multipole) == 'dipole'
        bnd = trim(si%multipole) == 'bond'
        nbd = 0
        ! --- Gamma ---
        call read_pao_species(trim(si%pao_file), trim(si%basis), si%rcut_bohr, pao, ierr, perr)
        if (ierr /= 0) then
            errmsg = perr; return
        end if
        call build_radial_density(pao, si%gamma_l, si%gamma_mul, si%dr_bohr, rd, ierr, perr)
        if (ierr /= 0) then
            errmsg = perr; return
        end if
        ! --- geometria: átomos centrais e eixos periódicos ---
        ic = inp
        ic%have_rho = .false.
        if (bnd) then
            ic%rho_rmax = si%rmax
            call run_extracted(ic, gtmp, ierr, perr, structure_only=.true., notes=nt, cent_pos=pos, &
                               cent_col=col, cent_atom=atom, bond_ij=bij, bond_d=bdd, bond_pos=bpos)
            nbd = size(bij, 2)
        else
            call run_extracted(ic, gtmp, ierr, perr, structure_only=.true., notes=nt, cent_pos=pos, &
                               cent_col=col, cent_atom=atom)
        end if
        if (ierr /= 0) then
            errmsg = 'scc (structure): ' // perr; return
        end if
        n = size(pos, 2)
        nsv = n
        call read_structure(trim(inp%structure_file), st, ierr, msg)
        if (ierr /= 0) then
            errmsg = 'scc: ' // trim(msg); return
        end if
        na = 0
        do a = 1, 3
            if (a == inp%embed_axis .or. any(inp%vacuum_axes == a)) cycle
            na = na + 1
            if (na <= 2) kax(na) = a
        end do
        if (na /= 2) then
            ierr = 2; errmsg = 'scc: needs exactly two periodic (k) axes'; return
        end if
        a1 = st%lattice(kax(1), :)
        a2 = st%lattice(kax(2), :)
        nhat = cross(a1, a2)
        nhat = nhat / norm2(nhat)
        ae = st%lattice(inp%embed_axis, :)
        ! substrato: lado do lead; normal aponta para longe dele
        select case (trim(si%substrate))
        case ('left')
            side = 1
        case ('right')
            side = 2
        case default
            side = merge(1, 2, trim(inp%left) == 'lead')
        end select
        if ((side == 1 .and. dot_product(nhat, ae) < 0.0_dp) .or. (side == 2 .and. dot_product(nhat, ae) > 0.0_dp)) &
            nhat = -nhat
        alpha = si%ewald_alpha
        if (alpha <= 0.0_dp) alpha = sqrt(PI / norm2(cross(a1, a2)))
        sy = trim(si%substrate) == 'none'
        if (sy) call add("scc: substrate = 'none': symmetric 2D Ewald G = 0 term (no substrate reference); " // &
                         "dipole cycle keeps mean V fixed (gauge)")
        if (bnd) then
            call bond_run()
            return
        end if
        allocate(gam(n, n))
        call scc_gamma_matrix(rd, pos, a1, a2, nhat, alpha, gam, uh, sy)
        write(msg, '(A,F10.5,A,F9.6,A,I0,A)') 'scc: U_H = ', uh, ' eV (Hartree, p density norm ', rd%norm, &
            ' before normalization); ', n, ' central atoms'
        call add(msg)
        if (dip) then
            allocate(gk(3, n, n), tk(3, 3, n, n))
            call scc_dipole_kernels(rd, pos, a1, a2, nhat, alpha, gk, tk, t0d, sy)
            write(msg, '(A,F12.6,A)') 'scc: dipoles on; T_ii (on-site) = ', t0d, ' eV/Bohr^2 = (4 pi/3) int rho_s^2'
            call add(msg)
            dc%h = si%dip_grid_bohr
            allocate(gkf(3, n, n))
            gkf = gk
        end if
        ! --- A14: meio dielétrico abaixo da região central (imagens) ---
        diel = trim(si%substrate) == 'dielectric'
        if (diel) then
            hmin = minval(matmul(nhat, pos))
            hbel = -huge(1.0_dp)
            do a = 1, st%nat
                if (dot_product(st%cart(:, a), nhat) < hmin - 1.0e-3_dp) &
                    hbel = max(hbel, dot_product(st%cart(:, a), nhat))
            end do
            if (hbel <= -huge(1.0_dp)) then
                ierr = 2; errmsg = "scc: substrate = 'dielectric' needs atoms below the central region " // &
                                   "(structure_file) to place the plane"; return
            end if
            zpl = 0.5_dp * (hmin + hbel)
            beta = (si%eps_sub - 1.0_dp) / (si%eps_sub + 1.0_dp)
            if (dip) then
                call scc_image_kernels(rd, pos, a1, a2, nhat, alpha, zpl, beta, gam, gk, gkf, tk)
            else
                call scc_image_kernels(rd, pos, a1, a2, nhat, alpha, zpl, beta, gam)
            end if
            write(msg, '(A,F9.3,A,F8.5,A,F9.4,A,F8.4,A)') "scc: substrate = 'dielectric': eps = ", si%eps_sub, &
                ', beta = (eps-1)/(eps+1) = ', beta, '; plane at h = ', zpl, ' A (normal), ', hmin - zpl, &
                ' A below the lowest central atom (midway to the atom below)'
            call add(msg)
        end if
        ! --- referência ---
        allocate(qref(n))
        if (dip) then
            allocate(pref(3, n))
            pref = 0.0_dp
        end if
        if (trim(si%qref_mode) == 'value') then
            qref = si%qref
        else
            call read_qref(ierr, errmsg)
            if (ierr /= 0) return
        end if
        ! --- entrada do ciclo ---
        ic = inp
        ic%kmode = 'list'
        ic%klist_avg_axis = si%k_avg_axis
        ic%klist_navg = si%k_navg
        ic%klist_axis = si%k_list_axis
        ic%klist_n = si%k_nlist
        ic%klist_k(1:si%k_nlist) = si%k_list(1:si%k_nlist)
        ic%klist_w(1:si%k_nlist) = si%k_w(1:si%k_nlist)
        ic%ne = 1
        ic%emin = si%e_fermi
        ic%emax = si%e_fermi
        ic%overlap_diagnostic = .false.
        ic%have_rho = .true.
        ic%rho_contour = .true.
        ic%rho_nwin = 1
        ic%rho_name(1) = 'occ'
        ic%rho_lo(1) = si%e_lo
        ic%rho_hi(1) = si%e_fermi
        ic%rho_ncontour = si%npoints
        ic%rho_contour_eta = 0.0_dp
        ic%rho_rmax = si%rmax
        ic%output_prefix = trim(inp%output_prefix) // '_scc'
        allocate(vin(n), vout(n), dq(n), f(n))
        vin = 0.0_dp
        pla = dip .and. si%planar_full
        npl = 0
        if (pla) then
            call planar_init(ierr, errmsg)
            if (ierr /= 0) return
            npl = pl%ny
        end if
        if (dip) then
            allocate(fin(3, n), fout(3, n), dpv(3, n), xv(4 * n + npl), xf(4 * n + npl))
            fin = 0.0_dp
        end if
        if (len_trim(si%v_init_file) > 0) then
            call read_vinit(ierr, errmsg)
            if (ierr /= 0) return
        end if
        mix%hist = si%mix_hist
        mix%beta = si%mix_beta
        open(newunit=ul, file=trim(inp%output_prefix) // '.scc_log.dat', status='replace', action='write')
        if (dip) then
            write(ul, '(A)') '# it  max|Vout-Vin| (eV)  max|Fout-Fin| (eV/Bohr)  sum dq (e)  max|dq| (e)  ' // &
                'max|dp| (e Bohr)  wall (s)  V_i (eV) ...'
            if (pla) then
                open(newunit=upl, file=trim(inp%output_prefix) // '.scc_planar_log.dat', status='replace', &
                     action='write')
                write(upl, '(A)') '# it  max|dVpl_out - dVpl_in| (eV)  int dn (e, full - bulk)  sum dq (e)  ' // &
                    'int dn_mult (e)  dVpl_out(top) (eV)  V_full(top) (eV)  V_mult(top) (eV)  max|dhp| (eV)'
            end if
            call dipole_cycle()
            close(ul)
            if (pla) close(upl)
            if (ierr /= 0) return
            call write_profile()
            if (pla) call write_planar()
            if (si%field_out) call write_field()
            if (.not. si%final_map) return
            ic = inp
            if (pla) then
                ic%rho_rmax = si%rmax
                call run_extracted(ic, gd, ierr, perr, nonconverged=nbad, notes=nt, scc_v=vin, scc_f=fin, &
                                   scc_pao=pao, scc_dip=dc, scc_dhp=dhp)
            else
                call run_extracted(ic, gd, ierr, perr, nonconverged=nbad, notes=nt, scc_v=vin, scc_f=fin, &
                                   scc_pao=pao, scc_dip=dc)
            end if
            if (ierr /= 0) then
                errmsg = 'scc (final map): ' // perr; return
            end if
            notes = [notes, nt]
            return
        end if
        write(ul, '(A)') '# it  max|Vout-Vin| (eV)  sum dq (e)  max|dq| (e)  wall (s)  V_i (eV) ...'
        do it = 1, si%maxiter
            t0 = wtime()
            call run_extracted(ic, gtmp, ierr, perr, nonconverged=nbad, scc_v=vin, rho_q=q, rho_files=.false.)
            if (ierr /= 0) then
                close(ul); errmsg = 'scc (iteration ' // itoa(it) // '): ' // perr; return
            end if
            dq = si%spin_factor * q(:, 1) - qref
            vout = matmul(gam, dq)
            f = vout - vin
            t1 = wtime()
            write(ul, '(I4,1X,ES12.4,1X,ES14.6,1X,ES12.4,1X,F9.2,*(1X,F11.6))') it, maxval(abs(f)), sum(dq), &
                maxval(abs(dq)), t1 - t0, vin
            flush(ul)
            write(msg, '(A,I0,A,ES10.3,A,ES11.3,A,F8.1,A,I0)') 'scc it ', it, ': max|dV| = ', maxval(abs(f)), &
                ' eV, sum dq = ', sum(dq), ' e, ', t1 - t0, ' s, nonconverged decimations ', nbad
            call add(msg)
            if (it == 1 .and. len_trim(si%qref_write) > 0) call write_qref()
            if (maxval(abs(f)) < si%tol_ev .or. si%one_shot) then
                converged = maxval(abs(f)) < si%tol_ev
                exit
            end if
            if (it == si%maxiter) exit
            call anderson_step(mix, vin, f)
        end do
        close(ul)
        call write_profile()
        if (.not. si%final_map) return
        ! --- mapa com o V da última entrada (convergida) ---
        ic = inp
        call run_extracted(ic, gd, ierr, perr, nonconverged=nbad, notes=nt, scc_v=vin)
        if (ierr /= 0) then
            errmsg = 'scc (final map): ' // perr; return
        end if
        notes = [notes, nt]

    contains

        subroutine add(s)
            character(len=*), intent(in) :: s
            character(len=256) :: l
            l = s
            notes = [notes, l]
        end subroutine add

        !> Ciclo com dipolos: Anderson no vetor (V, F); parada em max|dV| < tol_ev e max|dF| < tol_f.
        subroutine dipole_cycle()
            integer :: j
            do it = 1, si%maxiter
                t0 = wtime()
                if (pla) then
                    call run_extracted(ic, gtmp, ierr, perr, nonconverged=nbad, scc_v=vin, rho_q=q, &
                                       rho_files=si%rho_out, scc_f=fin, scc_pao=pao, scc_dip=dc, rho_p=pp, notes=nt, &
                                       rho_d=dmat, scc_dhp=dhp)
                else
                    call run_extracted(ic, gtmp, ierr, perr, nonconverged=nbad, scc_v=vin, rho_q=q, &
                                       rho_files=si%rho_out, scc_f=fin, scc_pao=pao, scc_dip=dc, rho_p=pp, notes=nt)
                end if
                if (ierr /= 0) then
                    errmsg = 'scc (iteration ' // itoa(it) // '): ' // perr; return
                end if
                if (it == 1) then
                    do j = 1, size(nt)
                        if (index(nt(j), 'dipoles:') > 0) call add(trim(nt(j)))
                    end do
                end if
                dq = si%spin_factor * q(:, 1) - qref
                dpv = si%spin_factor * pp(:, :, 1) - pref
                vout = matmul(gam, dq)
                do i = 1, n
                    do j = 1, n
                        vout(i) = vout(i) + dot_product(gk(:, i, j), dpv(:, j))
                    end do
                    fout(:, i) = 0.0_dp
                    do j = 1, n
                        fout(:, i) = fout(:, i) + gkf(:, i, j) * dq(j) - matmul(tk(:, :, i, j), dpv(:, j))
                    end do
                end do
                ! substrate = 'none' (A7): V só é definido a menos de constante; calibre fixo pela média
                ! (mean V_out = mean V_in, = 0 com V_init nulo). Neutralidade conferida por sum dq no log.
                if (sy) vout = vout - (sum(vout) - sum(vin)) / real(n, dp)
                f = vout - vin
                dvmax = maxval(abs(f))
                dfmax = maxval(abs(fout - fin))
                if (pla) then
                    call planar_out()
                    dvmax = max(dvmax, dvp)
                end if
                t1 = wtime()
                write(ul, '(I4,2(1X,ES12.4),1X,ES14.6,2(1X,ES12.4),1X,F9.2,*(1X,F11.6))') it, dvmax, dfmax, sum(dq), &
                    maxval(abs(dq)), maxval(abs(dpv)), t1 - t0, vin
                flush(ul)
                write(msg, '(A,I0,A,ES10.3,A,ES10.3,A,ES11.3,A,F8.1,A,I0)') 'scc it ', it, ': max|dV| = ', dvmax, &
                    ' eV, max|dF| = ', dfmax, ' eV/Bohr, sum dq = ', sum(dq), ' e, ', t1 - t0, &
                    ' s, nonconverged decimations ', nbad
                call add(msg)
                if (it == 1 .and. len_trim(si%qref_write) > 0) call write_qref()
                if ((dvmax < si%tol_ev .and. dfmax < si%tol_f) .or. si%one_shot) then
                    converged = dvmax < si%tol_ev .and. dfmax < si%tol_f
                    exit
                end if
                if (it == si%maxiter) exit
                xv(1:n) = vin
                xv(n + 1:4 * n) = reshape(fin, [3 * n])
                xf(1:n) = f
                xf(n + 1:4 * n) = reshape(fout - fin, [3 * n])
                if (pla) then
                    xv(4 * n + 1:) = pvin
                    xf(4 * n + 1:) = pvout - pvin
                end if
                call anderson_step(mix, xv, xf)
                vin = xv(1:n)
                fin = reshape(xv(n + 1:4 * n), [3, n])
                if (pla) then
                    pvin = xv(4 * n + 1:)
                    call planar_dh(pl, pvin, dhp)
                end if
            end do
        end subroutine dipole_cycle

        !> A19: malha, perfis, pares e referência; dV_in inicial (zero ou planar_init_file).
        subroutine planar_init(lerr, lmsg)
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            type(extracted_input_t) :: icp
            integer :: ax, c
            real(dp) :: tp0
            character(len=:), allocatable :: pm
            lerr = 0; lmsg = ''
            tp0 = wtime()
            ax = 0
            do c = 1, 3
                if (abs(abs(nhat(c)) - 1.0_dp) < 1.0e-10_dp) ax = c
            end do
            if (ax == 0) then
                lerr = 2; lmsg = 'scc: planar_full needs the surface normal along a Cartesian axis'; return
            end if
            icp = ic
            icp%rho_rmax = si%rmax
            call run_extracted(icp, gtmp, lerr, pm, structure_only=.true., rho_pair_ic=pic, rho_pair_ja=pja, &
                               rho_pair_kind=pkind, rho_pair_vec=pvec, rho_pair_s=pss)
            if (lerr /= 0) then
                lmsg = 'scc (planar pairs): ' // pm; return
            end if
            call planar_read_ref(trim(si%planar_ref_file), si%qref_nbasis, pao%norb, atom, ric, rvec, dref, nrow, &
                                 lerr, pm)
            if (lerr /= 0) then
                lmsg = 'scc: ' // pm; return
            end if
            call planar_setup(pl, pao, si%dip_grid_bohr, ax, sign(1.0_dp, nhat(ax)), matmul(nhat, pos), [pic, ric], &
                              reshape([pvec, rvec], [3, size(pic) + size(ric)]), [pkind, (0, c = 1, size(ric))], &
                              size(ric), si%planar_dy, si%planar_vac, norm2(cross(a1, a2)), rd%rho, rd%dr, lerr, pm)
            if (lerr /= 0) then
                lmsg = 'scc: ' // pm; return
            end if
            allocate(cfull(pl%ny), cref(pl%ny), cmul(pl%ny), vfull(pl%ny), vmul(pl%ny), pvin(pl%ny), pvout(pl%ny), &
                     dpn(n), dhp(pao%norb, pao%norb, pl%npair))
            call planar_charge(pl, dref, cref, ref=.true.)
            cref = si%spin_factor * cref
            pvin = 0.0_dp
            pvout = 0.0_dp
            cfull = 0.0_dp; cmul = 0.0_dp; vfull = 0.0_dp; vmul = 0.0_dp
            if (len_trim(si%planar_init_file) > 0) then
                call read_planar_init(lerr, lmsg)
                if (lerr /= 0) return
            end if
            call planar_dh(pl, pvin, dhp)
            write(msg, '(A,I0,A,F7.3,A,F10.4,A,F10.4,A,I0,A,I0,A,F7.1,A)') 'scc: planar_full (A19): ', pl%ny, &
                ' nodes, dy = ', pl%dy, ' A, h = ', pl%y(1), ' .. ', pl%y(pl%ny), ' A; ', pl%nprof, ' profiles, ', &
                pl%npair, ' pairs (', wtime() - tp0, ' s)'
            call add(msg)
            write(msg, '(A,ES14.6,A,I0,A,I0,A,*(1X,I0))') 'scc: planar reference (bulk rows): int n_ref = ', &
                sum(cref), ' e; ', size(ric), ' bulk pairs (system rows: ', size(pic), '); pairs per basis atom:', nrow
            call add(msg)
            write(msg, '(A,ES12.4)') 'scc: planar max|dV_in| = ', maxval(abs(pvin))
            call add(msg)
        end subroutine planar_init

        subroutine read_planar_init(lerr, lmsg)
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            integer :: u, ios, k
            real(dp) :: vals(9)
            character(len=512) :: ln
            lerr = 0; lmsg = ''
            open(newunit=u, file=trim(si%planar_init_file), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                lerr = 1; lmsg = 'scc: cannot open planar_init_file'; return
            end if
            k = 0
            do
                read(u, '(A)', iostat=ios) ln
                if (ios /= 0) exit
                if (ln(1:1) == '#' .or. len_trim(ln) == 0) cycle
                read(ln, *, iostat=ios) vals
                if (ios /= 0) cycle
                k = k + 1
                if (k > pl%ny) exit
                if (abs(vals(1) - pl%y(k)) > 1.0e-6_dp) exit
                pvin(k) = vals(9)
            end do
            close(u)
            if (k /= pl%ny) then
                lerr = 2; lmsg = 'scc: planar_init_file does not match the planar grid'; pvin = 0.0_dp
            end if
        end subroutine read_planar_init

        !> A19: densidade planar completa, multipolos e dV_out = phi[dn] - phi[dn_mult] desta iteração.
        subroutine planar_out()
            integer :: j, qq
            real(dp) :: tds(2)
            call planar_charge(pl, dmat, cfull)
            cfull = si%spin_factor * cfull
            do j = 1, n
                dpn(j) = dot_product(nhat, dpv(:, j))
            end do
            call planar_mult(pl, dq, dpn, cmul)
            call planar_pot(pl, cfull - cref, vfull)
            call planar_pot(pl, cmul, vmul)
            pvout = vfull - vmul
            dvp = maxval(abs(pvout - pvin))
            write(upl, '(I4,1X,ES12.4,3(1X,ES14.6),4(1X,ES13.5))') it, dvp, sum(cfull - cref), sum(dq), sum(cmul), &
                pvout(pl%ny), vfull(pl%ny), vmul(pl%ny), maxval(abs(dhp))
            flush(upl)
            if (it == 1) then
                ! teste 3: sum_k c_k = s sum_q tr(D_q S_q) = s sum_i q_DS (S da grade: exato; S do arquivo: erro da grade)
                tds = 0.0_dp
                do qq = 1, pl%npair
                    tds(1) = tds(1) + sum(dmat(:, :, qq) * pss(:, :, qq))
                    tds(2) = tds(2) + (sum(dmat(:, :, qq) * merge(transpose(sum(pl%prof(:, :, :, &
                        pl%pidx(qq)), dim=3)), sum(pl%prof(:, :, :, pl%pidx(qq)), dim=3), pl%ptr(qq))))
                end do
                tds = si%spin_factor * tds
                write(msg, '(A,ES22.14,A,ES10.3,A,ES10.3,A,ES22.14)') 'scc: planar sum_k n_k dy A = ', sum(cfull), &
                    ' e; - s tr(D S_grid) = ', sum(cfull) - tds(2), '; - s tr(D S_file) = ', &
                    sum(cfull) - tds(1), '; s sum_i q_DS = ', si%spin_factor * sum(q(:, 1))
                call add(msg)
            end if
        end subroutine planar_out

        subroutine write_planar()
            integer :: u, k
            real(dp) :: fa
            fa = 1.0_dp / (pl%area * pl%dy)
            open(newunit=u, file=trim(inp%output_prefix) // '.scc_planar.dat', status='replace', action='write')
            write(u, '(A)') '# scc planar G_par = 0 term (A19, last iteration): dV = phi[n - n_ref] - phi[n_mult]'
            write(u, '(A,F12.6,A,F10.6,A,I0,A,3F10.6)') '# area (A^2) = ', pl%area, ' ; dy (A) = ', pl%dy, &
                ' ; nodes = ', pl%ny, ' ; normal = ', nhat
            write(u, '(A)') '# h (A, along normal)  n  n_ref  dn  dn_mult (e/A^3, planar mean, spin summed)  ' // &
                'V_full = phi[dn]  V_mult = phi[dn_mult]  dV_out  dV_in (eV, electron energy, bottom reference)'
            do k = 1, pl%ny
                write(u, '(F12.6,4(1X,ES15.7),4(1X,ES15.7))') pl%y(k), cfull(k) * fa, cref(k) * fa, &
                    (cfull(k) - cref(k)) * fa, cmul(k) * fa, vfull(k), vmul(k), pvout(k), pvin(k)
            end do
            close(u)
        end subroutine write_planar

        !> A10: ciclo com cargas de ligação. Sítios 1..n = átomos (carga on-site e dipolo on-site), n+1..ns = pares
        !> (carga no ponto médio). Anderson em (V(ns), F(3n)); parada como no ciclo de dipolos.
        subroutine bond_run()
            integer :: j
            ns = n + nbd
            nsv = ns
            ngh = 0
            allocate(gpos(3, 0), gdq(0))
            allocate(qref(ns), pref(3, n))
            qref = 0.0_dp
            pref = 0.0_dp
            if (trim(si%qref_mode) == 'value') then
                qref(1:n) = si%qref
                call add("scc: bond mode with qref_mode = 'value': bond references are zero (use qref_write)")
            else
                call read_qref(ierr, errmsg)
                if (ierr /= 0) return
                call read_qref_bond(ierr, errmsg)
                if (ierr /= 0) return
            end if
            ! sítios: átomos, pares e (fixos) pares cortados da superfície, com dq = -1/2 q_ref (metade do lado do
            ! átomo do filme; o átomo ausente leva a outra metade, unidade neutra como os átomos ausentes do A5/A6)
            nall = ns + ngh
            allocate(spos(3, nall))
            spos(:, 1:n) = pos
            if (nbd > 0) spos(:, n + 1:ns) = bpos
            if (ngh > 0) spos(:, ns + 1:) = gpos
            allocate(gam(nall, nall))
            call scc_gamma_matrix(rd, spos, a1, a2, nhat, alpha, gam, uh, sy)
            write(msg, '(A,F10.5,A,I0,A,I0,A,I0,A,F9.5,A)') 'scc: U_H = ', uh, ' eV; bond charges on: ', n, &
                ' atoms + ', nbd, ' bond sites (midpoints) + ', ngh, ' cut bonds (fixed, sum dq = ', sum(gdq), ' e)'
            call add(msg)
            allocate(gk(3, n, n), tk(3, 3, n, n), gks(3, nall, n), gka(3, n, nall))
            call scc_dipole_kernels(rd, pos, a1, a2, nhat, alpha, gk, tk, t0d, sy)
            call scc_cross_kernel(rd, spos, pos, a1, a2, nhat, alpha, gks, sy)
            call scc_cross_kernel(rd, pos, spos, a1, a2, nhat, alpha, gka, sy)
            write(msg, '(A,F12.6,A)') 'scc: on-site dipoles; T_ii = ', t0d, ' eV/Bohr^2'
            call add(msg)
            dc%h = si%dip_grid_bohr
            ic = inp
            ic%kmode = 'list'
            ic%klist_avg_axis = si%k_avg_axis
            ic%klist_navg = si%k_navg
            ic%klist_axis = si%k_list_axis
            ic%klist_n = si%k_nlist
            ic%klist_k(1:si%k_nlist) = si%k_list(1:si%k_nlist)
            ic%klist_w(1:si%k_nlist) = si%k_w(1:si%k_nlist)
            ic%ne = 1
            ic%emin = si%e_fermi
            ic%emax = si%e_fermi
            ic%overlap_diagnostic = .false.
            ic%have_rho = .true.
            ic%rho_contour = .true.
            ic%rho_nwin = 1
            ic%rho_name(1) = 'occ'
            ic%rho_lo(1) = si%e_lo
            ic%rho_hi(1) = si%e_fermi
            ic%rho_ncontour = si%npoints
            ic%rho_contour_eta = 0.0_dp
            ic%rho_rmax = si%rmax
            ic%output_prefix = trim(inp%output_prefix) // '_scc'
            allocate(vin(ns), vout(ns), dq(nall), f(ns), fin(3, n), fout(3, n), dpv(3, n), xv(ns + 3 * n), &
                     xf(ns + 3 * n))
            vin = 0.0_dp
            fin = 0.0_dp
            if (len_trim(si%v_init_file) > 0) then
                call read_vinit(ierr, errmsg)
                if (ierr /= 0) return
            end if
            mix%hist = si%mix_hist
            mix%beta = si%mix_beta
            open(newunit=ul, file=trim(inp%output_prefix) // '.scc_log.dat', status='replace', action='write')
            write(ul, '(A)') '# it  max|Vout-Vin| (eV)  max|Fout-Fin| (eV/Bohr)  sum dq (e, with the fixed cut bonds)  ' // &
                'max|dq| atoms  ' // &
                'max|dq| bonds (e)  max|dp| (e Bohr)  wall (s)  V_i atoms (eV) ...'
            do it = 1, si%maxiter
                t0 = wtime()
                call run_extracted(ic, gtmp, ierr, perr, nonconverged=nbad, rho_files=si%rho_out, scc_f=fin, &
                                   scc_pao=pao, scc_dip=dc, scc_vb=vin, rho_qb=qbd, rho_qon=qon, rho_pon=pon, notes=nt)
                if (ierr /= 0) then
                    close(ul); errmsg = 'scc (iteration ' // itoa(it) // '): ' // perr; return
                end if
                if (it == 1) then
                    do j = 1, size(nt)
                        if (index(nt(j), 'dipoles:') > 0 .or. index(nt(j), 'bond') > 0) call add(trim(nt(j)))
                    end do
                end if
                dq(1:n) = si%spin_factor * qon(:, 1) - qref(1:n)
                if (nbd > 0) dq(n + 1:ns) = si%spin_factor * qbd(:, 1) - qref(n + 1:ns)
                if (ngh > 0) dq(ns + 1:) = gdq
                dpv = si%spin_factor * pon(:, :, 1) - pref
                vout = matmul(gam(1:ns, :), dq)
                do i = 1, ns
                    do j = 1, n
                        vout(i) = vout(i) + dot_product(gks(:, i, j), dpv(:, j))
                    end do
                end do
                do i = 1, n
                    fout(:, i) = 0.0_dp
                    do j = 1, nall
                        fout(:, i) = fout(:, i) + gka(:, i, j) * dq(j)
                    end do
                    do j = 1, n
                        fout(:, i) = fout(:, i) - matmul(tk(:, :, i, j), dpv(:, j))
                    end do
                end do
                if (sy) vout = vout - (sum(vout) - sum(vin)) / real(ns, dp)
                f = vout - vin
                dvmax = maxval(abs(f))
                dfmax = maxval(abs(fout - fin))
                t1 = wtime()
                write(ul, '(I4,2(1X,ES12.4),1X,ES14.6,3(1X,ES12.4),1X,F9.2,*(1X,F11.6))') it, dvmax, dfmax, sum(dq), &
                    maxval(abs(dq(1:n))), maxval(abs(dq(n + 1:ns))), maxval(abs(dpv)), t1 - t0, vin(1:n)
                flush(ul)
                write(msg, '(A,I0,A,ES10.3,A,ES10.3,A,ES11.3,A,F8.1,A,I0)') 'scc it ', it, ': max|dV| = ', dvmax, &
                    ' eV, max|dF| = ', dfmax, ' eV/Bohr, sum dq = ', sum(dq), ' e, ', t1 - t0, &
                    ' s, nonconverged decimations ', nbad
                call add(msg)
                if (it == 1 .and. len_trim(si%qref_write) > 0) call write_qref_bond()
                if ((dvmax < si%tol_ev .and. dfmax < si%tol_f) .or. si%one_shot) then
                    converged = dvmax < si%tol_ev .and. dfmax < si%tol_f
                    exit
                end if
                if (it == si%maxiter) exit
                xv(1:ns) = vin
                xv(ns + 1:) = reshape(fin, [3 * n])
                xf(1:ns) = f
                xf(ns + 1:) = reshape(fout - fin, [3 * n])
                call anderson_step(mix, xv, xf)
                vin = xv(1:ns)
                fin = reshape(xv(ns + 1:), [3, n])
            end do
            close(ul)
            call write_profile_bond()
            if (.not. si%final_map) return
            ic = inp
            ic%rho_rmax = si%rmax
            call run_extracted(ic, gd, ierr, perr, nonconverged=nbad, notes=nt, scc_vb=vin, scc_f=fin, &
                               scc_pao=pao, scc_dip=dc)
            if (ierr /= 0) then
                errmsg = 'scc (final map): ' // perr; return
            end if
            notes = [notes, nt]
        end subroutine bond_run

        !> Chave do par no volume: átomos de base (b_i <= b_j) e d de i para j; b_i = b_j: d com a primeira
        !> componente não nula positiva.
        subroutine bkey(ib, nbas, bi, bj, dv)
            integer, intent(in) :: ib, nbas
            integer, intent(out) :: bi, bj
            real(dp), intent(out) :: dv(3)
            integer :: c, t
            logical :: neg
            bi = modulo(bij(3, ib) - 1, nbas) + 1
            bj = modulo(bij(4, ib) - 1, nbas) + 1
            dv = bdd(:, ib)
            neg = .false.
            do c = 1, 3
                if (abs(dv(c)) > 1.0e-6_dp) then
                    neg = dv(c) < 0.0_dp
                    exit
                end if
            end do
            if (bj < bi .or. (bi == bj .and. neg)) then
                t = bi; bi = bj; bj = t
                dv = -dv
            end if
        end subroutine bkey

        subroutine read_qref_bond(lerr, lmsg)
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            integer :: u, ios, m, k, bi, bj, ki, kj, ib
            real(dp) :: dv(3), qq, kd(3)
            real(dp), allocatable :: rq(:), rd_(:, :)
            integer, allocatable :: rb(:, :)
            character(len=512) :: ln
            lerr = 0; lmsg = ''
            allocate(rq(0), rd_(3, 0), rb(2, 0))
            open(newunit=u, file=trim(si%qref_file), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                lerr = 1; lmsg = 'scc: cannot open qref_file'; return
            end if
            do
                read(u, '(A)', iostat=ios) ln
                if (ios /= 0) exit
                if (ln(1:4) /= 'bond') cycle
                read(ln(5:), *, iostat=ios) bi, bj, dv, qq
                if (ios /= 0) cycle
                rq = [rq, qq]
                rd_ = reshape([rd_, dv], [3, size(rq)])
                rb = reshape([rb, bi, bj], [2, size(rq)])
            end do
            close(u)
            m = size(rq)
            do ib = 1, nbd
                call bkey(ib, si%qref_nbasis, ki, kj, kd)
                do k = 1, m
                    if (rb(1, k) == ki .and. rb(2, k) == kj .and. norm2(rd_(:, k) - kd) < 1.0e-3_dp) exit
                end do
                if (k > m) then
                    lerr = 2
                    write(ln, '(A,I0,A,I0,1X,I0,3F10.4)') 'scc: qref_file lacks bond ', ib, ' (basis pair, d): ', &
                        ki, kj, kd
                    lmsg = trim(ln); return
                end if
                qref(n + ib) = rq(k)
            end do
            ! pares cortados: pares do volume a partir de um átomo central que não existem no sistema
            do i = 1, n
                ki = modulo(atom(i) - 1, si%qref_nbasis) + 1
                do k = 1, m
                    do kj = 1, 2
                        if (rb(kj, k) /= ki) cycle
                        kd = merge(1.0_dp, -1.0_dp, kj == 1) * rd_(:, k)     ! i como b_i (+d) ou como b_j (-d)
                        call cut_check(i, kd, rq(k))
                    end do
                end do
            end do
        end subroutine read_qref_bond

        !> A10: par (atom i, d) do volume ausente da lista de pares -> par cortado fixo, dq = -1/2 q_ref.
        subroutine cut_check(ia, dv, qr)
            integer, intent(in) :: ia
            real(dp), intent(in) :: dv(3), qr
            integer :: ib
            do ib = 1, nbd
                if (bij(1, ib) == ia .and. norm2(bdd(:, ib) - dv) < 1.0e-3_dp) return
                if (bij(2, ib) == ia .and. norm2(bdd(:, ib) + dv) < 1.0e-3_dp) return
            end do
            ngh = ngh + 1
            gpos = reshape([gpos, pos(:, ia) + 0.5_dp * dv], [3, ngh])
            gdq = [gdq, -0.5_dp * qr]
        end subroutine cut_check

        subroutine write_qref_bond()
            integer :: u, b, nb_, c, ib, k, m, ki, kj
            real(dp), allocatable :: pb(:, :), kq(:), kdv(:, :), kmx(:)
            integer, allocatable :: kb(:, :), kc(:), kof(:)
            real(dp) :: pdev, kd(3)
            nb_ = max(si%qref_nbasis, 1)
            allocate(qb(nb_), cntb(nb_), pb(3, nb_))
            qb = 0.0_dp; cntb = 0; pb = 0.0_dp
            do i = 1, n
                b = modulo(atom(i) - 1, nb_) + 1
                qb(b) = qb(b) + si%spin_factor * qon(i, 1)
                pb(:, b) = pb(:, b) + si%spin_factor * pon(:, i, 1)
                cntb(b) = cntb(b) + 1
            end do
            open(newunit=u, file=trim(si%qref_write), status='replace', action='write')
            write(u, '(A)') '# scc q_ref (bond mode, A10): basis atom b, mean on-site s*q over central copies, copies, ' // &
                'max |q - mean|, mean on-site p_ref (x y z, e Bohr; electron = -e), max |p - mean|'
            write(u, '(A)') '# bond lines: bond b_i b_j dx dy dz (Å, i -> j) mean s*q_b copies max |q_b - mean| |d|'
            do b = 1, nb_
                qb(b) = qb(b) / real(max(cntb(b), 1), dp)
                pb(:, b) = pb(:, b) / real(max(cntb(b), 1), dp)
                pdev = 0.0_dp
                do i = 1, n
                    if (modulo(atom(i) - 1, nb_) + 1 /= b) cycle
                    do c = 1, 3
                        pdev = max(pdev, abs(si%spin_factor * pon(c, i, 1) - pb(c, b)))
                    end do
                end do
                write(u, '(I0,1X,ES24.16,1X,I0,1X,ES12.4,3(1X,ES24.16),1X,ES12.4)') b, qb(b), cntb(b), &
                    maxval(abs(si%spin_factor * qon(:, 1) - qb(b)), mask=modulo(atom - 1, nb_) + 1 == b), pb(:, b), pdev
            end do
            ! pares: agrupa por chave
            allocate(kq(nbd), kdv(3, nbd), kmx(nbd), kb(2, nbd), kc(nbd), kof(nbd))
            m = 0
            do ib = 1, nbd
                call bkey(ib, nb_, ki, kj, kd)
                do k = 1, m
                    if (kb(1, k) == ki .and. kb(2, k) == kj .and. norm2(kdv(:, k) - kd) < 1.0e-4_dp) exit
                end do
                if (k > m) then
                    m = m + 1
                    kb(:, m) = [ki, kj]; kdv(:, m) = kd; kq(m) = 0.0_dp; kc(m) = 0
                end if
                kof(ib) = k
                kq(k) = kq(k) + si%spin_factor * qbd(ib, 1)
                kc(k) = kc(k) + 1
            end do
            kq(1:m) = kq(1:m) / real(kc(1:m), dp)
            kmx(1:m) = 0.0_dp
            do ib = 1, nbd
                k = kof(ib)
                kmx(k) = max(kmx(k), abs(si%spin_factor * qbd(ib, 1) - kq(k)))
            end do
            do k = 1, m
                write(u, '(A,I0,1X,I0,3(1X,ES24.16),1X,ES24.16,1X,I0,1X,ES12.4,1X,F10.5)') 'bond ', kb(:, k), &
                    kdv(:, k), kq(k), kc(k), kmx(k), norm2(kdv(:, k))
            end do
            close(u)
        end subroutine write_qref_bond

        subroutine write_profile_bond()
            integer :: u, ib
            open(newunit=u, file=trim(inp%output_prefix) // '.scc_profile.dat', status='replace', action='write')
            write(u, '(A)') '# scc profile (bond mode, A10; last iteration): V_in is the potential used'
            write(u, '(A,F12.6,A,I0,A,L1)') '# U_H (eV) = ', uh, ' ; iterations = ', it, ' ; converged = ', converged
            write(u, '(A,3F12.6)') '# normal (to vacuum) = ', nhat
            write(u, '(A,F12.6,A,I0,A,I0)') '# T_ii = ', t0d, ' eV/Bohr^2; atoms ', n, ' ; bond sites ', nbd
            write(u, '(A)') '# atoms: i atom col x y z height q_on(spin) q_ref dq V_out V_in p(3) p_ref(3) dp(3) ' // &
                'F_out(3) F_in(3)  (on-site charge and dipole)'
            write(u, '(A)') '# bonds: i=n+b atom_i atom_j xm ym zm height q_b(spin) q_ref dq V_out V_in ic jc |d| ' // &
                'dx dy dz  (midpoint; jc = -1: j in the seam)'
            do i = 1, n
                write(u, '(I0,1X,I0,1X,I0,4(1X,F12.6),1X,F14.9,1X,F12.8,1X,ES14.6,2(1X,F12.7),15(1X,ES15.7))') &
                    i, atom(i), col(i), pos(:, i), dot_product(pos(:, i), nhat), si%spin_factor * qon(i, 1), &
                    qref(i), dq(i), vout(i), vin(i), si%spin_factor * pon(:, i, 1), pref(:, i), dpv(:, i), &
                    fout(:, i), fin(:, i)
            end do
            do ib = 1, nbd
                i = n + ib
                write(u, '(I0,1X,I0,1X,I0,4(1X,F12.6),1X,F14.9,1X,F12.8,1X,ES14.6,2(1X,F12.7),2(1X,I0),4(1X,F10.5))') &
                    i, bij(3, ib), bij(4, ib), spos(:, i), dot_product(spos(:, i), nhat), si%spin_factor * qbd(ib, 1), &
                    qref(i), dq(i), vout(i), vin(i), bij(1:2, ib), norm2(bdd(:, ib)), bdd(:, ib)
            end do
            close(u)
        end subroutine write_profile_bond

        real(dp) function wtime()
            integer(8) :: c, r
            call system_clock(c, r)
            wtime = real(c, dp) / real(r, dp)
        end function wtime

        subroutine read_qref(lerr, lmsg)
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            integer :: u, ios, b, nc
            real(dp) :: qq, md, pv(3)
            real(dp), allocatable :: qbv(:), pbv(:, :)
            logical, allocatable :: got(:)
            character(len=256) :: ln
            lerr = 0; lmsg = ''
            allocate(qbv(si%qref_nbasis), got(si%qref_nbasis), pbv(3, si%qref_nbasis))
            got = .false.
            open(newunit=u, file=trim(si%qref_file), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                lerr = 1; lmsg = 'scc: cannot open qref_file'; return
            end if
            do
                read(u, '(A)', iostat=ios) ln
                if (ios /= 0) exit
                if (ln(1:1) == '#' .or. len_trim(ln) == 0 .or. ln(1:4) == 'bond') cycle
                if (dip .or. bnd) then
                    read(ln, *, iostat=ios) b, qq, nc, md, pv
                    if (ios /= 0) then
                        close(u); lerr = 2; lmsg = 'scc: dipole mode needs p_ref columns in qref_file'; return
                    end if
                    if (b >= 1 .and. b <= si%qref_nbasis) pbv(:, b) = pv
                else
                    read(ln, *, iostat=ios) b, qq
                end if
                if (ios /= 0 .or. b < 1 .or. b > si%qref_nbasis) cycle
                qbv(b) = qq; got(b) = .true.
            end do
            close(u)
            if (.not. all(got)) then
                lerr = 2; lmsg = 'scc: qref_file lacks some basis atoms'; return
            end if
            do i = 1, n
                qref(i) = qbv(modulo(atom(i) - 1, si%qref_nbasis) + 1)
                if (dip .or. bnd) pref(:, i) = pbv(:, modulo(atom(i) - 1, si%qref_nbasis) + 1)
            end do
        end subroutine read_qref

        subroutine write_qref()
            integer :: u, b, nb_, c
            real(dp), allocatable :: pb(:, :)
            real(dp) :: pdev
            nb_ = max(si%qref_nbasis, 1)
            allocate(qb(nb_), cntb(nb_), pb(3, nb_))
            qb = 0.0_dp; cntb = 0; pb = 0.0_dp
            do i = 1, n
                b = modulo(atom(i) - 1, nb_) + 1
                qb(b) = qb(b) + si%spin_factor * q(i, 1)
                if (dip) pb(:, b) = pb(:, b) + si%spin_factor * pp(:, i, 1)
                cntb(b) = cntb(b) + 1
            end do
            open(newunit=u, file=trim(si%qref_write), status='replace', action='write')
            if (dip) then
                write(u, '(A)') '# scc q_ref: basis atom b, mean s*q_DS over central copies, copies, max |q - mean|, ' // &
                    'mean p_ref (x y z, e Bohr; electron = -e), max |p - mean|'
            else
                write(u, '(A)') '# scc q_ref: basis atom b, mean s*q_DS over central copies, copies, max |q - mean|'
            end if
            do b = 1, nb_
                qb(b) = qb(b) / real(max(cntb(b), 1), dp)
                if (dip) then
                    pb(:, b) = pb(:, b) / real(max(cntb(b), 1), dp)
                    pdev = 0.0_dp
                    do i = 1, n
                        if (modulo(atom(i) - 1, nb_) + 1 /= b) cycle
                        do c = 1, 3
                            pdev = max(pdev, abs(si%spin_factor * pp(c, i, 1) - pb(c, b)))
                        end do
                    end do
                    write(u, '(I0,1X,ES24.16,1X,I0,1X,ES12.4,3(1X,ES24.16),1X,ES12.4)') b, qb(b), cntb(b), &
                        maxval(abs(si%spin_factor * q(:, 1) - qb(b)), mask=modulo(atom - 1, nb_) + 1 == b), pb(:, b), pdev
                else
                    write(u, '(I0,1X,ES24.16,1X,I0,1X,ES12.4)') b, qb(b), cntb(b), &
                        maxval(abs(si%spin_factor * q(:, 1) - qb(b)), mask=modulo(atom - 1, nb_) + 1 == b)
                end if
            end do
            close(u)
        end subroutine write_qref

        subroutine read_vinit(lerr, lmsg)
            integer, intent(out) :: lerr
            character(len=:), allocatable, intent(out) :: lmsg
            integer :: u, ios, k, ii, kf
            real(dp) :: vals(11), vald(26)
            character(len=512) :: ln
            lerr = 0; lmsg = ''
            open(newunit=u, file=trim(si%v_init_file), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                lerr = 1; lmsg = 'scc: cannot open v_init_file'; return
            end if
            k = 0
            kf = 0
            do
                read(u, '(A)', iostat=ios) ln
                if (ios /= 0) exit
                if (ln(1:1) == '#' .or. len_trim(ln) == 0) cycle
                read(ln, *, iostat=ios) ii, vals
                if (ios /= 0) cycle
                if (ii < 1 .or. ii > nsv) cycle
                vin(ii) = vals(11)
                k = k + 1
                if (dip .or. bnd) then
                    read(ln, *, iostat=ios) ii, vald
                    if (ios == 0) then
                        fin(:, ii) = vald(24:26)
                        kf = kf + 1
                    end if
                end if
            end do
            if ((dip .or. bnd) .and. kf /= n) call add('scc: v_init_file without F_in columns; F starts at zero')
            close(u)
            if (k /= nsv) then
                lerr = 2; lmsg = 'scc: v_init_file does not cover all central atoms'
            end if
        end subroutine read_vinit

        !> A7: campo pontual nos núcleos com referência no átomo neutro (q_i - field_qatom, p_i total).
        subroutine write_field()
            integer :: u
            real(dp), allocatable :: qn(:), pt(:, :), ef(:, :), es(:, :), er(:, :), ed(:, :), ws(:, :), eq(:, :), &
                                     zp(:, :)
            real(dp) :: ev2v
            allocate(qn(n), pt(3, n), ef(3, n), es(3, n), er(3, n), ed(3, n), ws(3, n), eq(3, n), zp(3, n))
            qn = si%spin_factor * q(:, 1) - si%field_qatom
            pt = si%spin_factor * pp(:, :, 1)
            call scc_point_field(rd, pos, a1, a2, nhat, alpha, qn, pt, ef, es, sy)
            ! partes (lineares): referência do volume (q_ref - q_atom, p_ref), desvio de carga dq e de dipolo dp
            zp = 0.0_dp
            call scc_point_field(rd, pos, a1, a2, nhat, alpha, qref - si%field_qatom, pref, er, ws, sy)
            call scc_point_field(rd, pos, a1, a2, nhat, alpha, dq, zp, eq, ws, sy)
            call scc_point_field(rd, pos, a1, a2, nhat, alpha, 0.0_dp * dq, dpv, ed, ws, sy)
            ev2v = 1.0_dp / BOHR_ANG          ! eV/Bohr -> V/Å
            open(newunit=u, file=trim(inp%output_prefix) // '.scc_field.dat', status='replace', action='write')
            write(u, '(A)') '# scc point field at the nuclei (A7): reference = neutral atom (q - q_atom, total p)'
            write(u, '(A,F10.5,A,A)') '# q_atom = ', si%field_qatom, ' ; substrate = ', trim(si%substrate)
            write(u, '(A,ES12.4,A)') '# rho_s(0) = ', rd%rho(0), ' Bohr^-3 (self term -(4 pi/3) rho_s(0) p)'
            write(u, '(A,3F12.6)') '# normal = ', nhat
            write(u, '(A)') '# i atom col x y z (Å) q-q_atom (e) p(3) (e Bohr) E(3) (V/Å, with self) ' // &
                'E_self(3) (V/Å) F_out(3) F_in(3) (V/Å, rho_s-averaged, bulk reference) ' // &
                'E_ref(3) (from q_ref, p_ref) E_dq(3) (from dq) E_dp(3) (from dp) (V/Å, point, with self)'
            do i = 1, n
                write(u, '(I0,1X,I0,1X,I0,3(1X,F12.6),1X,ES15.7,24(1X,ES15.7))') i, atom(i), col(i), pos(:, i), &
                    qn(i), pt(:, i), ef(:, i), es(:, i), fout(:, i) * ev2v, fin(:, i) * ev2v, er(:, i), eq(:, i), &
                    ed(:, i)
            end do
            close(u)
        end subroutine write_field

        subroutine write_profile()
            integer :: u
            open(newunit=u, file=trim(inp%output_prefix) // '.scc_profile.dat', status='replace', action='write')
            write(u, '(A)') '# scc profile (last iteration): V_in is the potential used, V_out = Gamma dq'
            write(u, '(A,F12.6,A,I0,A,L1)') '# U_H (eV) = ', uh, ' ; iterations = ', it, ' ; converged = ', converged
            write(u, '(A,3F12.6)') '# normal (to vacuum) = ', nhat
            if (dip) then
                write(u, '(A,F12.6,A)') '# dipoles: T_ii (on-site) = ', t0d, ' eV/Bohr^2; p in e Bohr (electron = -e), ' // &
                    'F in eV/Bohr (electric field)'
                write(u, '(A)') '# i atom col x y z height(Å, along normal) q_DS(spin) q_ref dq V_out V_in' // &
                    ' p(3) p_ref(3) dp(3) F_out(3) F_in(3)'
            else
                write(u, '(A)') '# i atom col x y z height(Å, along normal) q_DS(spin) q_ref dq V_out V_in'
            end if
            do i = 1, n
                if (dip) then
                    write(u, '(I0,1X,I0,1X,I0,4(1X,F12.6),1X,F14.9,1X,F12.8,1X,ES14.6,2(1X,F12.7),15(1X,ES15.7))') &
                        i, atom(i), col(i), pos(:, i), dot_product(pos(:, i), nhat), si%spin_factor * q(i, 1), &
                        qref(i), dq(i), vout(i), vin(i), si%spin_factor * pp(:, i, 1), pref(:, i), dpv(:, i), &
                        fout(:, i), fin(:, i)
                else
                    write(u, '(I0,1X,I0,1X,I0,4(1X,F12.6),1X,F14.9,1X,F12.8,1X,ES14.6,2(1X,F12.7))') i, atom(i), &
                        col(i), pos(:, i), dot_product(pos(:, i), nhat), si%spin_factor * q(i, 1), qref(i), dq(i), &
                        vout(i), vin(i)
                end if
            end do
            close(u)
        end subroutine write_profile

    end subroutine run_scc

    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=16) :: b
        write(b, '(I0)') i
        s = trim(b)
    end function itoa

end module scc
