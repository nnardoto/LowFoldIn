! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

! ----------------------------------------------------------------------------------------------------------
! tools/openmx2hsio.f90
!
! Converte um arquivo binario .scfout do OpenMX (SCFOUT_VERSION = 3,
! OpenMX 3.9/4.0) para o formato binario HSDataFormat (mod_hs_io.f90)
! usado pelo LowFoldIn, em Fortran puro -- sem depender de Python.
!
! Uso:
!   openmx2hsio <arquivo.scfout> [--out-prefix PREFIXO] [--fermi E_F_HARTREE]
!               [--no-fermi-shift] [--write-structure] [--structure-file ARQ]
!
! Gera:
!   PREFIXO_S.hsio
!   PREFIXO_H.hsio                       (SpinP_switch == 0)
! ou
!   PREFIXO_H_up.hsio, PREFIXO_H_dn.hsio (SpinP_switch == 1, spin colinear;
!                                          um S.hsio so', compartilhado pelos
!                                          dois canais de spin -- rode o
!                                          `unfold` duas vezes, uma por canal,
!                                          com o MESMO mapeamento atomo->PC)
! ou
!   PREFIXO_H.hsio, PREFIXO_S.hsio       (SpinP_switch == 3, spin nao-colinear
!                                          / spinor, NOrbs dobrado. Escrito em
!                                          ordem ATOMO-MAJOR -- up e down
!                                          contiguos por atomo -- via
!                                          ReorderSpinorAtomMajor, a UNICA
!                                          ordem compativel com Orbs(:) do
!                                          formato HSDataFormat. A ordem
!                                          "spin-major" bruta que sai de
!                                          AssembleNonCollinearH/scfout NAO e'
!                                          usavel diretamente por ferramentas
!                                          que indexam orbitais por atomo via
!                                          Orbs(:), como o `unfold` -- ver
!                                          comentario de ReorderSpinorAtomMajor
!                                          em mod_openmx_convert.f90.)
!
! Por padrao, o nivel de Fermi (ChemP) e' lido automaticamente do
! proprio .scfout e H e' deslocado para que E_F caia em 0 eV, via
! H'(R) = H(R) - E_F*S(R) para todo vizinho R.
!
! Opcionalmente (--write-structure, ou --structure-file ARQ que tambem o
! ativa) escreve tambem PREFIXO.structure.nml (ou ARQ): rede (Angstrom),
! especies com numero de orbitais e posicoes atomicas, extraidos da geometria
! e do input OpenMX embutido no .scfout. Esse arquivo e' o structure_file
! usado por `spectral` e `atomselect` (ver docs/spectral.md e
! docs/atomselect.md). O caminho nao pode coincidir com o .scfout nem com os
! .hsio gerados.
!
! Codigo de saida: 0 em sucesso (ou --help); 1 em qualquer erro.
!
! Nathanael N. Batista, Agosto de 2026
! ----------------------------------------------------------------------------------------------------------
program openmx2hsio
    use types,          only: dp
    use scfout_io,      only: ScfoutData, ReadScfoutFile
    use openmx_convert, only: HARTREE_TO_EV, AssembleDense, AssembleNonCollinearH, &
                               DoubleBlockDiagonal, BuildIndexList, HermiticityReport, &
                               ApplyFermiShift, ReorderSpinorAtomMajor
    use HS_io,          only: HSDataFormat
    implicit none

    ! O .scfout guarda geometria em Bohr; o structure.nml usa Angstrom.
    real(dp), parameter :: BOHR_TO_ANG = 0.529177210903_dp

    character(len=1024) :: scfout_path, out_prefix, structure_path, arg
    real(dp) :: fermi_hartree
    logical  :: has_fermi_override, no_fermi_shift, write_structure, custom_structure_path
    integer  :: nargs, iarg, ios

    type(ScfoutData) :: scfout
    integer, allocatable :: IndexList(:, :)
    complex(dp), allocatable :: Sdata(:, :, :), Hdata(:, :, :)
    integer :: NOrbs, NOrbs2, neighbors
    real(dp) :: fermi_hartree_final, fermi_ev
    logical  :: apply_shift
    character(len=64) :: fermi_source

    ! -------------------------------------------------------------
    ! 1. Parsing da linha de comando
    ! -------------------------------------------------------------
    nargs = command_argument_count()
    if (nargs < 1) then
        call print_usage()
        stop 1
    end if

    call get_command_argument(1, scfout_path)
    out_prefix = default_prefix(trim(scfout_path))
    has_fermi_override = .false.
    no_fermi_shift = .false.
    write_structure = .false.
    custom_structure_path = .false.
    fermi_hartree = 0.0_dp

    iarg = 2
    do while (iarg <= nargs)
        call get_command_argument(iarg, arg)
        select case (trim(arg))
        case ("--out-prefix")
            iarg = iarg + 1
            if (iarg > nargs) then
                print *, "Erro: --out-prefix precisa de um valor"
                stop 1
            end if
            call get_command_argument(iarg, out_prefix)
        case ("--fermi")
            iarg = iarg + 1
            if (iarg > nargs) then
                print *, "Erro: --fermi precisa de um valor (Hartree)"
                stop 1
            end if
            call get_command_argument(iarg, arg)
            read(arg, *, iostat=ios) fermi_hartree
            if (ios /= 0) then
                print *, "Erro: valor invalido para --fermi: ", trim(arg)
                stop 1
            end if
            has_fermi_override = .true.
        case ("--no-fermi-shift")
            no_fermi_shift = .true.
        case ("--write-structure")
            write_structure = .true.
        case ("--structure-file")
            iarg = iarg + 1
            if (iarg > nargs) then
                print *, "Erro: --structure-file precisa de um caminho"
                stop 1
            end if
            call get_command_argument(iarg, structure_path)
            write_structure = .true.
            custom_structure_path = .true.
        case ("-h", "--help")
            call print_usage()
            stop 0
        case default
            print *, "Aviso: argumento desconhecido ignorado: ", trim(arg)
        end select
        iarg = iarg + 1
    end do

    ! Caminho padrao da estrutura (a menos que --structure-file tenha sido dado)
    if (write_structure .and. .not. custom_structure_path) then
        structure_path = trim(out_prefix) // ".structure.nml"
    end if

    ! -------------------------------------------------------------
    ! 2. Le o .scfout
    ! -------------------------------------------------------------
    call ReadScfoutFile(scfout, trim(scfout_path))
    ! Falha cedo (antes de qualquer escrita) se a estrutura sobrescreveria um arquivo em uso
    if (write_structure) call ValidateStructurePath(trim(structure_path), trim(scfout_path), &
                                                    trim(out_prefix))

    ! -------------------------------------------------------------
    ! 3. Decide o nivel de Fermi: --fermi (override) > ChemP do
    !    .scfout > nenhum shift (--no-fermi-shift)
    ! -------------------------------------------------------------
    apply_shift = .not. no_fermi_shift
    if (no_fermi_shift) then
        fermi_hartree_final = 0.0_dp
        fermi_source = "(desabilitado)"
    else if (has_fermi_override) then
        fermi_hartree_final = fermi_hartree
        fermi_source = "fornecido via --fermi"
    else
        fermi_hartree_final = scfout%ChemP
        fermi_source = "lido do proprio .scfout (ChemP)"
    end if
    fermi_ev = fermi_hartree_final * HARTREE_TO_EV

    call BuildIndexList(scfout, IndexList)
    neighbors = size(IndexList, 2)

    ! -------------------------------------------------------------
    ! 4. Monta e escreve S e H conforme SpinP_switch
    ! -------------------------------------------------------------
    select case (scfout%SpinP_switch)
    case (0, 1)
        call AssembleDense(scfout, scfout%OLP_flat, Sdata, NOrbs)
        call WriteHS(trim(out_prefix) // "_S.hsio", "S", scfout, NOrbs, neighbors, IndexList, Sdata)

        if (scfout%SpinP_switch == 0) then
            call AssembleDense(scfout, scfout%Hks_flat(:, 0), Hdata, NOrbs)
            Hdata = Hdata * HARTREE_TO_EV
            if (apply_shift) call ApplyFermiShift(Hdata, Sdata, fermi_ev)
            call WriteHS(trim(out_prefix) // "_H.hsio", "H, nao-polarizado, eV", &
                         scfout, NOrbs, neighbors, IndexList, Hdata)
        else
            call AssembleDense(scfout, scfout%Hks_flat(:, 0), Hdata, NOrbs)
            Hdata = Hdata * HARTREE_TO_EV
            if (apply_shift) call ApplyFermiShift(Hdata, Sdata, fermi_ev)
            call WriteHS(trim(out_prefix) // "_H_up.hsio", "H, spin up, eV", &
                         scfout, NOrbs, neighbors, IndexList, Hdata)

            call AssembleDense(scfout, scfout%Hks_flat(:, 1), Hdata, NOrbs)
            Hdata = Hdata * HARTREE_TO_EV
            if (apply_shift) call ApplyFermiShift(Hdata, Sdata, fermi_ev)
            call WriteHS(trim(out_prefix) // "_H_dn.hsio", "H, spin dn, eV", &
                         scfout, NOrbs, neighbors, IndexList, Hdata)
        end if

    case (3)
        block
            complex(dp), allocatable :: Sspatial(:, :, :), Sspinor(:, :, :)
            complex(dp), allocatable :: HdataReordered(:, :, :), SspinorReordered(:, :, :)
            integer :: NOrbsSpatial

            call AssembleDense(scfout, scfout%OLP_flat, Sspatial, NOrbsSpatial)
            call DoubleBlockDiagonal(Sspatial, Sspinor)
            NOrbs2 = 2 * NOrbsSpatial

            call AssembleNonCollinearH(scfout, Hdata, NOrbs2)
            Hdata = Hdata * HARTREE_TO_EV
            if (apply_shift) call ApplyFermiShift(Hdata, Sspinor, fermi_ev)

            call HermiticityReport(Hdata, "H (spinor, ordem original spin-major)")
            call HermiticityReport(Sspinor, "S (spinor, ordem original spin-major)")

            ! Reordena de "spin-major" (todos os orbitais up, depois todos
            ! os down) para "atomo-major" (up e down contiguos por atomo),
            ! a UNICA ordem compativel com Orbs(:) do formato HSDataFormat.
            ! Pura permutacao de base: autovalores e hermiticidade
            ! preservados, so' a ordem das linhas/colunas muda.
            call ReorderSpinorAtomMajor(scfout, Hdata, Sspinor, HdataReordered, SspinorReordered)

            call WriteHS(trim(out_prefix) // "_H.hsio", &
                         "H, nao-colinear, eV, ordem atomo-major (up,dn contiguos por atomo)", &
                         scfout, NOrbs2, neighbors, IndexList, HdataReordered)
            call WriteHS(trim(out_prefix) // "_S.hsio", &
                         "S, nao-colinear, ordem atomo-major (up,dn contiguos por atomo)", &
                         scfout, NOrbs2, neighbors, IndexList, SspinorReordered)
        end block

    case default
        print '(A,I0)', "Erro: SpinP_switch inesperado: ", scfout%SpinP_switch
        stop 1
    end select

    ! Estrutura escrita apos os .hsio, a partir dos mesmos dados do .scfout
    if (write_structure) call WriteStructureFile(trim(structure_path), scfout)

    ! -------------------------------------------------------------
    ! 5. Resumo no console
    ! -------------------------------------------------------------
    print *
    print '(A)',    "Arquivo lido: " // trim(scfout_path)
    print '(A,I0)', "  atomnum        = ", scfout%atomnum
    print '(A,I0)', "  SpinP_switch   = ", scfout%SpinP_switch
    print '(A,I0,A,I0,A)', "  TCpyCell       = ", scfout%TCpyCell, " (Neighbors = ", neighbors, ")"
    print '(A,I0)', "  Solver         = ", scfout%Solver
    print '(A,F16.9,A)', "  ChemP (bruto)  = ", scfout%ChemP, " Hartree (lido do .scfout)"
    print '(A,F12.4,A)', "  E_Temp         = ", scfout%E_Temp, " K"
    print '(A,F14.9)', "  Unidades de H  = eV (fator ", HARTREE_TO_EV
    print '(A)',    "  Unidades de S  = adimensional (Overlap nao e' escalado)"
    if (apply_shift) then
        print '(A,F12.6,A,A,A)', "  Nivel de Fermi = deslocado para 0 eV (E_F = ", fermi_ev, &
                                  " eV, ", trim(fermi_source), ")"
    else
        print '(A)', "  Nivel de Fermi = NAO deslocado (--no-fermi-shift)."
    end if

contains

    subroutine print_usage()
        print '(A)', "Uso: openmx2hsio <arquivo.scfout> [--out-prefix PREFIXO] " // &
                      "[--fermi E_F_HARTREE] [--no-fermi-shift] [--write-structure]"
        print '(A)', ""
        print '(A)', "  --out-prefix PREFIXO   prefixo dos arquivos de saida " // &
                      "(padrao: nome do .scfout sem extensao)"
        print '(A)', "  --fermi E_F_HARTREE    sobrescreve o nivel de Fermi (por padrao, " // &
                      "lido automaticamente do ChemP do .scfout)"
        print '(A)', "  --no-fermi-shift       desliga o deslocamento do nivel de Fermi"
        print '(A)', "  --write-structure      escreve PREFIXO.structure.nml a partir da geometria e do input embutido"
        print '(A)', "  --structure-file FILE  escreve a estrutura em FILE (tambem ativa --write-structure)"
    end subroutine print_usage

    function default_prefix(path) result(prefix)
        character(len=*), intent(in) :: path
        character(len=1024) :: prefix
        integer :: dot_pos, slash_pos, i

        slash_pos = 0
        do i = len_trim(path), 1, -1
            if (path(i:i) == '/' ) then
                slash_pos = i
                exit
            end if
        end do

        dot_pos = 0
        do i = len_trim(path), slash_pos + 1, -1
            if (path(i:i) == '.') then
                dot_pos = i
                exit
            end if
        end do

        if (dot_pos > slash_pos + 1) then
            prefix = path(slash_pos + 1 : dot_pos - 1)
        else
            prefix = path(slash_pos + 1 : len_trim(path))
        end if
    end function default_prefix

    !> Aborta (stop 1) se `path` coincidir com o .scfout de entrada ou com algum
    !> .hsio de saida (_S, _H, _H_up, _H_dn), para nao sobrescreve-los.
    subroutine ValidateStructurePath(path, source_path, prefix)
        character(len=*), intent(in) :: path, source_path, prefix
        character(len=1024) :: reserved(4)
        integer :: i

        reserved(1) = source_path
        reserved(2) = trim(prefix) // '_S.hsio'
        reserved(3) = trim(prefix) // '_H.hsio'
        reserved(4) = trim(prefix) // '_H_up.hsio'
        do i = 1, size(reserved)
            if (trim(path) == trim(reserved(i))) then
                print '(A,A)', 'Erro: --structure-file colide com um arquivo de entrada/saida: ', trim(path)
                stop 1
            end if
        end do
        if (trim(path) == trim(prefix) // '_H_dn.hsio') then
            print '(A,A)', 'Erro: --structure-file colide com um arquivo de saida: ', trim(path)
            stop 1
        end if
    end subroutine ValidateStructurePath

    !> Escreve o structure.nml (&STRUCTURE, &SPECIES, &ATOMS) em Angstrom.
    !> Requer os rotulos atomicos recuperados do input embutido no .scfout.
    !> O numero de orbitais por especie e' TotalNumOrbs (x2 em spinor,
    !> SpinP_switch == 3, para casar com o NOrbs dobrado do .hsio) e deve ser
    !> o mesmo para todos os atomos da mesma especie.
    subroutine WriteStructureFile(filename, scf)
        character(len=*), intent(in) :: filename
        type(ScfoutData), intent(in) :: scf

        character(len=8), allocatable :: species_label(:)
        integer, allocatable :: species_norb(:)
        integer :: iu, ios_local, iatom, ispecies, nspecies, norb_factor
        logical :: known_species

        if (.not. scf%has_atom_labels) then
            print '(A)', "Erro: nao foi possivel recuperar Atoms.SpeciesAndCoordinates do input embutido."
            print '(A)', "       Use um .scfout com o input OpenMX completo ou forneca a estrutura manualmente."
            stop 1
        end if
        if (.not. allocated(scf%cart_bohr)) then
            print '(A)', "Erro: coordenadas atomicas ausentes no .scfout."
            stop 1
        end if

        ! Fator de spinor e lista de especies unicas (ordem de primeira aparicao)
        norb_factor = merge(2, 1, scf%SpinP_switch == 3)
        allocate(species_label(scf%atomnum), species_norb(scf%atomnum))
        nspecies = 0
        do iatom = 1, scf%atomnum
            known_species = .false.
            do ispecies = 1, nspecies
                if (trim(species_label(ispecies)) == trim(scf%atom_label(iatom))) then
                    known_species = .true.
                    if (species_norb(ispecies) /= norb_factor * scf%TotalNumOrbs(iatom)) then
                        print '(A,A)', "Erro: contagem de orbitais inconsistente para a especie ", &
                                      trim(scf%atom_label(iatom))
                        stop 1
                    end if
                    exit
                end if
            end do
            if (.not. known_species) then
                nspecies = nspecies + 1
                species_label(nspecies) = scf%atom_label(iatom)
                species_norb(nspecies) = norb_factor * scf%TotalNumOrbs(iatom)
            end if
        end do

        open(newunit=iu, file=filename, status='replace', action='write', iostat=ios_local)
        if (ios_local /= 0) then
            print '(A,A)', "Erro: nao foi possivel escrever a estrutura: ", trim(filename)
            stop 1
        end if

        write(iu, '(A)') '&STRUCTURE'
        do iatom = 1, 3
            write(iu, '(A,I0,A,3(1X,ES24.16))') '    LatticeVectors(', iatom, ',:) =', &
                                                scf%lattice_bohr(iatom, :) * BOHR_TO_ANG
        end do
        write(iu, '(A)') '/'
        write(iu, '(A)') ''
        write(iu, '(A)') '&SPECIES'
        write(iu, '(I0)') nspecies
        do ispecies = 1, nspecies
            write(iu, '(A,1X,I0)') trim(species_label(ispecies)), species_norb(ispecies)
        end do
        write(iu, '(A)') '/'
        write(iu, '(A)') ''
        write(iu, '(A)') '&ATOMS'
        write(iu, '(I0)') scf%atomnum
        do iatom = 1, scf%atomnum
            write(iu, '(A,3(1X,ES24.16))') trim(scf%atom_label(iatom)), &
                                           scf%cart_bohr(:, iatom) * BOHR_TO_ANG
        end do
        write(iu, '(A)') '/'
        close(iu)
        print '(A,A)', "  Escrito: ", trim(filename)
    end subroutine WriteStructureFile

    ! Monta um objeto HSDataFormat e chama WriteFile (o MESMO codigo que
    ! o Fortran usa para ler depois -- garante compatibilidade binaria
    ! sem precisar reimplementar o layout do arquivo).
    subroutine WriteHS(filename, tag, scf, nOrbsLocal, nNeighLocal, IdxList, Data)
        character(len=*), intent(in) :: filename, tag
        type(ScfoutData), intent(in) :: scf
        integer, intent(in) :: nOrbsLocal, nNeighLocal
        integer, intent(in) :: IdxList(:, :)
        complex(dp), intent(in) :: Data(:, :, :)

        type(HSDataFormat) :: out_data

        out_data%Title      = "scfout2hsio (Fortran): " // trim(scfout_path) // " (" // trim(tag) // ")"
        out_data%Neighbors  = nNeighLocal
        out_data%nAtoms     = scf%atomnum
        out_data%NOrbs      = nOrbsLocal

        allocate(out_data%Degenerecency(nNeighLocal))
        out_data%Degenerecency = 1

        allocate(out_data%Orbs(scf%atomnum))
        if (nOrbsLocal == sum(scf%TotalNumOrbs)) then
            out_data%Orbs = scf%TotalNumOrbs
        else
            ! caso spinor (NOrbs dobrado): reporta 2x orbitais por atomo
            out_data%Orbs = 2 * scf%TotalNumOrbs
        end if

        allocate(out_data%IndexList(3, nNeighLocal))
        out_data%IndexList = IdxList

        allocate(out_data%Data(nOrbsLocal, nOrbsLocal, nNeighLocal))
        out_data%Data = Data

        call out_data%WriteFile(filename)
        print '(A,A)', "  Escrito: ", trim(filename)
    end subroutine WriteHS

end program openmx2hsio
