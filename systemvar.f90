! ==============================================================================
! Module: SYSTEMVAR
! 
! Purpose: Contains the global variables needed in Vassembly.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            09/01/2013   Original Code
!
! Contains:
!
! Modules -
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      MODULE SYSTEMVAR

        IMPLICIT NONE


        !=== ALLOCATABLE ARRAYS ===!

        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: gb
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: psum

        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: gbond

        INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: lbond
        INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: mapnn
        INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: maphp
        INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: maptwo

        INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nncp
        INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: mapauto

        !=== VARIABLES ===!

        !=== Gas Constant in kcal / (mol*K) ===!
        DOUBLE PRECISION, PARAMETER :: gcons = 1.987206d-3

        DOUBLE PRECISION, SAVE :: ratep1,ratep2,rkap1,rkap2
        DOUBLE PRECISION, SAVE :: vol,temp,beta,time,tfinal,tnext

        INTEGER, SAVE :: nbond,ncp,nps,nnps,nnmax
        INTEGER, SAVE :: npro,nrna,nsum,ITOT

      END MODULE
