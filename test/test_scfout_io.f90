! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

program test_scfout_io
    use scfout_io, only: ParseEmbeddedAtomLabels
    use testing, only: check, test_summary
    implicit none

    character(len=256) :: lines(8)
    character(len=8) :: labels(3)
    logical :: found

    lines = ''
    lines(1) = '# representative OpenMX input'
    lines(2) = '<Atoms.SpeciesAndCoordinates' // achar(0) // 'ignored bytes'
    lines(3) = '  1 Te 0.0 0.0 0.0 8.0 8.0'
    lines(4) = '  2 Bi 0.0 0.0 1.0 7.5 7.5'
    lines(5) = '  # comment inside the block'
    lines(6) = '  3 Sb 0.0 0.0 2.0 7.5 7.5'
    lines(7) = 'Atoms.SpeciesAndCoordinates>'
    call ParseEmbeddedAtomLabels(lines, 3, labels, found)
    call check(found, 'complete atom block found')
    call check(all(labels == ['Te      ', 'Bi      ', 'Sb      ']), 'labels preserve atom order')

    lines(6) = '  2 Sb 0.0 0.0 2.0'
    call ParseEmbeddedAtomLabels(lines, 3, labels, found)
    call check(.not. found, 'duplicate atom index rejected')

    lines = ''
    lines(1) = '<Atoms.SpeciesAndCoordinates'
    lines(2) = '1 SpeciesNameTooLong 0.0 0.0 0.0'
    lines(3) = 'Atoms.SpeciesAndCoordinates>'
    call ParseEmbeddedAtomLabels(lines, 3, labels, found)
    call check(.not. found, 'label longer than structure format rejected')

    lines = ''
    call ParseEmbeddedAtomLabels(lines, 3, labels, found)
    call check(.not. found, 'missing atom block reported')

    lines = ''
    lines(1) = '# <Atoms.SpeciesAndCoordinates'
    lines(2) = '! 1 Wrong 0.0 0.0 0.0'
    lines(3) = '<Atoms.SpeciesAndCoordinates'
    lines(4) = '1 Te 0.0 0.0 0.0'
    lines(5) = '2 Bi 0.0 0.0 1.0'
    lines(6) = '3 Sb 0.0 0.0 2.0'
    lines(7) = 'Atoms.SpeciesAndCoordinates>'
    call ParseEmbeddedAtomLabels(lines, 3, labels, found)
    call check(found .and. labels(1) == 'Te', 'commented delimiter ignored')

    lines(7) = ''
    call ParseEmbeddedAtomLabels(lines, 3, labels, found)
    call check(.not. found, 'missing closing delimiter rejected')

    lines(4) = '999 Te 0.0 0.0 0.0'
    call ParseEmbeddedAtomLabels(lines, 3, labels, found)
    call check(.not. found, 'out-of-range atom index rejected safely')

    lines(4) = 'not_an_index Te 0.0 0.0 0.0'
    call ParseEmbeddedAtomLabels(lines, 3, labels, found)
    call check(.not. found, 'malformed atom index rejected safely')

    call test_summary('test_scfout_io')
end program test_scfout_io
