! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2026 Nathanael Nardoto Batista

module types
    use, intrinsic :: iso_fortran_env, only: real32, real64, real128
    implicit none
    private
    
    ! Parâmetros públicos de Kind de precisão
    public :: sp, dp, qp

    ! 1. Definições de Kinds de Precisão
    integer, parameter :: sp = real32  ! Single Precision (Precisão Simples - 4 bytes)
    integer, parameter :: dp = real64  ! Double Precision (Precisão Dupla - 8 bytes)
    integer, parameter :: qp = real128 ! Quad Precision (Precisão Quádrupla - 16 bytes)

end module types
