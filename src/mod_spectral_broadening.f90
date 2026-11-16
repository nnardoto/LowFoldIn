! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

!> Alargamento de dados espectrais discretos (energia, peso) em mapas contínuos
!> A(k,E) sobre uma grade de energia uniforme, mais a ocupação de Fermi-Dirac
!> (docs/design/contracts.md §7, specmap).
!>
!> Posição no pipeline: pós-processamento. O specmap lê os pares discretos
!> (energia de banda, peso de projeção) produzidos por `spectral`
!> (módulo spectral_data, contracts §5) e usa este módulo para transformá-los
!> em uma imagem A(E,k) suave; nada aqui toca autovetores nem E/S.
!>
!> Convenções:
!>   - Energias, larguras e a grade estão todas em eV; temperatura em K.
!>   - Ordem dos arrays: energy(nband, nk), weight(nband, nk), A(ne, nk), isto é,
!>     bandas primeiro e índice da grade de energia primeiro na saída.
!>   - Lorentziana: largura = HWHM eta,           L(x) = (1/pi) eta/(x^2+eta^2)
!>   - Gaussiana:   largura = desvio padrão sigma (NÃO HWHM),
!>                  G(x) = exp(-x^2/(2 sigma^2)) / (sigma sqrt(2 pi))
!>   - Ambos os núcleos integram 1 sobre a reta real; os pesos são usados com
!>     o sinal (sem abs, sem truncamento).
!> As rotinas da biblioteca nunca param o programa. Os chamadores devem validar
!> width > 0, ne >= 2 e emax > emin antes (o specmap faz isso na validação de entrada).
module spectral_broadening
    use types, only: dp
    implicit none
    private

    ! ==================================================================
    ! Constantes
    ! ==================================================================

    !> Seletor de núcleo para `broaden`.
    integer, parameter, public :: KERNEL_LORENTZIAN = 1, KERNEL_GAUSSIAN = 2

    !> Constante de Boltzmann em eV/K.
    real(dp), parameter, public :: KB_EV_PER_K = 8.617333262e-5_dp

    !> Contribuições gaussianas além deste número de sigma são ignoradas
    !> (exp(-72) ~ 5e-32, muito abaixo da relevância da dupla precisão).
    real(dp), parameter :: GAUSS_CUTOFF_NSIGMA = 12.0_dp

    real(dp), parameter :: PI = acos(-1.0_dp)

    public :: lorentzian, gaussian, fermi_dirac, energy_grid, broaden

contains

    ! ==================================================================
    ! Núcleos e ocupação (elementais, utilizáveis em arrays)
    ! ==================================================================

    !> Lorentziana normalizada com meia largura a meia altura eta.
    elemental function lorentzian(x, eta) result(y)
        real(dp), intent(in) :: x, eta
        real(dp) :: y
        y = (eta / PI) / (x*x + eta*eta)
    end function lorentzian

    !> Gaussiana normalizada com desvio padrão sigma.
    elemental function gaussian(x, sigma) result(y)
        real(dp), intent(in) :: x, sigma
        real(dp) :: y
        y = exp(-x*x / (2.0_dp*sigma*sigma)) / (sigma * sqrt(2.0_dp*PI))
    end function gaussian

    !> Ocupação de Fermi-Dirac f = 1/(1+exp((e-mu)/(kB T))), segura contra overflow.
    !> T == 0 dá a função degrau com f(mu) = 1/2. T < 0 não tem significado
    !> e deve ser rejeitado pelo chamador; é tratado como T = 0.
    !> e, mu em eV; temperature_k em K.
    elemental function fermi_dirac(e, mu, temperature_k) result(f)
        real(dp), intent(in) :: e, mu, temperature_k
        real(dp) :: f
        real(dp) :: z, ez

        ! T <= 0: função degrau (f = 1/2 exatamente em e == mu).
        if (temperature_k <= 0.0_dp) then
            if (e < mu) then
                f = 1.0_dp
            else if (e > mu) then
                f = 0.0_dp
            else
                f = 0.5_dp
            end if
            return
        end if

        ! z = (e-mu)/(kB T). Segurança contra overflow: além de |z| = 700 exp()
        ! estouraria, então satura; para 0 < z usa-se exp(-z), de modo que o
        ! argumento de exp nunca seja positivo.
        z = (e - mu) / (KB_EV_PER_K * temperature_k)
        if (z > 700.0_dp) then
            f = 0.0_dp
        else if (z < -700.0_dp) then
            f = 1.0_dp
        else if (z > 0.0_dp) then
            ez = exp(-z)
            f = ez / (1.0_dp + ez)
        else
            f = 1.0_dp / (1.0_dp + exp(z))
        end if
    end function fermi_dirac

    ! ==================================================================
    ! Grade e alargamento
    ! ==================================================================

    !> Grade uniforme E_i = emin + (i-1)(emax-emin)/(ne-1), i = 1..ne (ne >= 2).
    !> O último ponto é fixado exatamente em emax (sem deriva de arredondamento).
    !> grid deve ter tamanho >= ne; arrays menores são preenchidos só até size(grid).
    subroutine energy_grid(emin, emax, ne, grid)
        real(dp), intent(in)  :: emin, emax
        integer,  intent(in)  :: ne
        real(dp), intent(out) :: grid(:)
        integer :: i
        real(dp) :: de

        de = (emax - emin) / real(max(ne - 1, 1), dp)
        do i = 1, min(ne, size(grid))
            grid(i) = emin + real(i - 1, dp) * de
        end do
        if (ne >= 2 .and. size(grid) >= ne) grid(ne) = emax
    end subroutine energy_grid

    !> A(:,k) = sum_n weight(n,k) * K(grid - energy(n,k)).
    !> energy(nband,nk), weight(nband,nk), grid(ne) -> A(ne,nk).
    !> As caudas lorentzianas são sempre calculadas por inteiro (sem corte); as
    !> contribuições gaussianas com |x| > 12 sigma são ignoradas. Núcleo desconhecido -> A = 0.
    !> width é eta (HWHM) para a Lorentziana e sigma para a Gaussiana, em eV.
    !> ne é obtido de size(grid); A deve ter shape (ne, nk). Bandas com
    !> peso exatamente 0 são ignoradas. Custo O(nk * nband * ne).
    subroutine broaden(energy, weight, grid, kernel, width, A)
        real(dp), intent(in)  :: energy(:,:), weight(:,:), grid(:)
        integer,  intent(in)  :: kernel
        real(dp), intent(in)  :: width
        real(dp), intent(out) :: A(:,:)

        integer  :: ik, ib, ie, nband, nk, ne
        real(dp) :: x, w, cut

        nband = size(energy, 1)
        nk    = size(energy, 2)
        ne    = size(grid)
        A = 0.0_dp
        cut = GAUSS_CUTOFF_NSIGMA * width

        do ik = 1, nk
            do ib = 1, nband
                w = weight(ib, ik)
                if (w == 0.0_dp) cycle
                select case (kernel)
                case (KERNEL_LORENTZIAN)
                    ! Caudas lentas 1/x^2: sem corte, então toda banda contribui a todo ponto da grade.
                    do ie = 1, ne
                        x = grid(ie) - energy(ib, ik)
                        A(ie, ik) = A(ie, ik) + w * lorentzian(x, width)
                    end do
                case (KERNEL_GAUSSIAN)
                    do ie = 1, ne
                        x = grid(ie) - energy(ib, ik)
                        if (abs(x) > cut) cycle
                        A(ie, ik) = A(ie, ik) + w * gaussian(x, width)
                    end do
                end select
            end do
        end do
    end subroutine broaden

end module spectral_broadening
