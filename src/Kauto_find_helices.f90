!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
! SUBROUTINE FIND_HELIX
!
! Purpose: scan a structure for segments of a molecule that are a helix
!
!
! Algorithm: 
! (1) get coordinates array and NR scalar from main program
! (2) get helix-type constraints array (alpha, 310, B-DNA,...) from main
! (3) for NR total atoms in molecule, and KR atoms from helix-type 
!      constraint, do the following NR - KR + 1 times:
!      (i) 
!
! Arguments:
!
!  input===========
!
!  natoms:             
!  oldZ(natoms):       same as oldX(natoms), but the z-components...
!
!  output==========
!   
!--------------------------------------------------------------------------

SUBROUTINE find_helix()
IMPLICIT NONE
!--declare passenger variables to/fro MAIN..........
   integer, intent(in)              :: natoms
   double precision, intent(in)     :: phi_idx, theta_idx
   double precision, dimension(:),   &
                     allocatable,    &
                     intent(in)     :: oldX, oldY, oldZ
   double precision, dimension(:,:), &
                     allocatable,    &
                     intent(out)    :: new_xyz
!--declare internal variables that MAIN don-t need...
   double precision                  &
        :: coef1, coef2, coef3, coef4, coef5, coef6, coef7, coef8, coef9
   double precision                 &
        :: o1, o2, o3, r
   integer                          :: iatom

!--ALLOCATE the arrays whose dimensions needed user input variable natoms..
   ALLOCATE(new_xyz     (natoms,3))


! =========================================================================
! =========================================================================

END SUBROUTINE find_helix

