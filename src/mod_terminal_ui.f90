! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

module terminal_ui
    implicit none
    
    ! Largura fixa padrão para toda a interface do software
    integer, parameter :: UI_WIDTH = 65
    integer, parameter :: dp       = kind(1.0d0)

contains

    ! 1. Desenha uma linha divisória horizontal simples
    subroutine print_line(char)
        character(len=1), intent(in) :: char
        print '(A)', repeat(char, UI_WIDTH)
    end subroutine print_line

    ! 2. Desenha o topo/base de uma caixa fechada
    subroutine print_box_close()
        print '(A)', "+" // repeat("-", UI_WIDTH - 2) // "+"
    end subroutine print_box_close

    ! 3. Desenha o TOPO de uma caixa com texto centralizado
    subroutine print_box(text)
        character(len=*), intent(in)  :: text
        character(len=:), allocatable :: clean_text
        integer                       :: text_len, total_spaces, left_spaces, right_spaces

        clean_text = trim(adjustl(text))
        text_len   = len_trim(clean_text)

        total_spaces = UI_WIDTH - 2 - text_len
        left_spaces  = total_spaces / 2
        right_spaces = total_spaces - left_spaces

        call print_box_close()
        print '(A, A, A, A, A)', "|", repeat(" ", left_spaces), clean_text, repeat(" ", right_spaces), "|"
        call print_box_close()
    end subroutine print_box

    ! 4. NOVA SUBROTINA: Imprime uma linha interna alinhada à esquerda com fechamento automático
    subroutine print_row(label, value)
        character(len=*), intent(in) :: label
        character(len=*), intent(in) :: value
        character(len=:), allocatable :: full_line
        integer                      :: content_len, remaining_spaces

        ! Monta o texto da linha combinando a chave e o valor
        full_line = "  " // trim(label) // trim(adjustl(value))
        content_len = len_trim(full_line)
        
        ! Calcula o preenchimento para bater na borda direita fixa
        remaining_spaces = UI_WIDTH - 2 - content_len
        
        print '(A, A, A, A)', "|", trim(full_line), repeat(" ", remaining_spaces), "|"
    end subroutine print_row

    ! 5. Desenha um CABEÇALHO DE PAINEL (Linhas duplas com texto centralizado)
    subroutine print_header(title)
        character(len=*), intent(in) :: title
        character(len=:), allocatable :: clean_title
        integer                      :: title_len, total_spaces, left_spaces, right_spaces

        clean_title = trim(adjustl(title))
        title_len   = len_trim(clean_title)

        total_spaces = UI_WIDTH - 2 - title_len
        left_spaces  = total_spaces / 2
        right_spaces = total_spaces - left_spaces

        call print_line("=")
        print '(A, A, A, A, A)', "|", repeat(" ", left_spaces), clean_title, repeat(" ", right_spaces), "|"
        call print_line("=")
    end subroutine print_header

    ! 6. Cabeçalho específico para a tabela de K-Points
    subroutine print_kpath_table_header()
        call print_header("CAMINHO DE PONTOS K")
        print '(A6, 3A14, A12)', "  Path", "Kx", "Ky", "Kz", "NKPTS"
        call print_line("-")
    end subroutine print_kpath_table_header

    ! 7. Linha de dados da tabela de K-Points
    subroutine print_kpath_row(idx, kpoints, nkpts)
        integer,      intent(in) :: idx
        real(dp),     intent(in) :: kpoints(3)
        integer,      intent(in) :: nkpts
        
        print '(I6, " : ", 3F14.8, I11)', idx, kpoints, nkpts
    end subroutine print_kpath_row

        ! Imprime uma linha com a Chave alinhada à esquerda e o Valor alinhado à direita
    subroutine print_row_split(label, value)
        character(len=*), intent(in) :: label
        character(len=*), intent(in) :: value
        character(len=:), allocatable :: clean_label, clean_value
        integer                      :: label_len, value_len, middle_spaces

        ! Limpa os espaços em branco das pontas das strings
        clean_label = trim(adjustl(label))
        clean_value = trim(adjustl(value))
        
        label_len = len_trim(clean_label)
        value_len = len_trim(clean_value)

        ! Calcula o espaço exato que deve ficar no MEIO entre a chave e o valor
        ! Descontamos as 2 barras laterais '|' e os 2 espaços de respiro nas bordas ("  ")
        middle_spaces = UI_WIDTH - 2 - 4 - label_len - value_len

        ! Se o texto for grande demais e o espaço virar negativo, garante pelo menos 1 espaço
        if (middle_spaces < 1) middle_spaces = 1

        ! Imprime com alinhamento perfeito nas duas extremidades
        print '(A, A, A, A, A, A, A)', "|  ", clean_label, &
            repeat(" ", middle_spaces), clean_value, "  |"
    end subroutine print_row_split
end module terminal_ui
