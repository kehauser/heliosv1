!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
! SUBROUTINE ROTATE
!
! Purpose: given phi AND theta, rotate X, Y, Z coordinates over a unit sphere
!
! Arguments:
!
!  input===========
!  natoms:             number of atoms, needed to allocate new_XYZ
!
!  phi_idx:          grid counter in \phi
!                   = loop_step(iPhi) * grid_size(dPhi) + phi_start(phi_beg)
!  theta_idx:        grid counter in \theta
!                   = loop_step(iTheta) * grd_siz(dTheta) + the_strt(the_beg)
!  oldX(natoms):       array of the x-components of the input coordinates
!                   whose dimension is the number of atoms (natoms)
!  oldY(natoms):       same as oldX(natoms), but the y-components...
!  oldZ(natoms):       same as oldX(natoms), but the z-components...
!
!  output==========
!  new_xyz(natoms,3):   array of rotated coordinates 
!                  whose dimensions are the number of atoms (natoms) 
!                   by 3 (for the X, the Y, and the Z components)
!   
!--------------------------------------------------------------------------

SUBROUTINE rotate (natoms, phi_idx, theta_idx, oldX, oldY, oldZ, new_xyz)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ROTATE COORDINATES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Setup factors of elements of the transformation matix
  o1  = SIN(phi_idx) * COS(theta_idx)
  o2  = SIN(phi_idx) * SIN(theta_idx)
  o3  = COS(phi_idx)
  r   = SQRT(1-(o1*o1))

  !Setup elements of the transformation matrix. 
  ! See EQN(k) in accompanying PDF documentation.
  !
  !     ( x )   ( coef1  coef2  coef3 ) ( X )
  !     ( y ) = ( coef4  coef5  coef6 ) ( Y )
  !     ( z )   ( coef7  coef8  coef9 ) ( Z )
  !
  coef1 =   SQRT(1.0-(o1*o1))
  coef2 =   0
  coef3 =   o1
  coef4 =  -o1*o2 / r
  coef5 =   o3    / r
  coef6 =   o2
  coef7 =  -o1*o3 / r
  coef8 =  -o2    / r
  coef9 =   o3

! Linear Algebra: loop over iatoms, we transform the coordinates
  DO iatom = 1,natoms         ! Loop-Atoms
     new_xyz(iatom,1) = coef1*oldX(iatom) + coef2*oldY(iatom) + coef3*oldZ(iatom)
     new_xyz(iatom,2) = coef4*oldX(iatom) + coef5*oldY(iatom) + coef6*oldZ(iatom)
     new_xyz(iatom,3) = coef7*oldX(iatom) + coef8*oldY(iatom) + coef9*oldZ(iatom)
  ENDDO

! =========================================================================
! =========================================================================

END SUBROUTINE rotate
