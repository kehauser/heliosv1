!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
! SUBROUTINE PITCHANG:
!
! Purpose: given rotated coordinates, calculate the pitch and sweep of helix
!
! Arguments:
!
!  input===========
!  nhelix:         number of atoms in the helix
!                  =natoms if ssDNA, =natoms/2 if dsDNA
!  strandBstart:   =1, if ssDNA (dsDNA = 0)
!                  =(natoms/2), if dsDNA (dsDNA = 0)
!  X0, Y0:         center of the circle, needed to calculate angles
!  newX(nhelix):       rotated X coordinates (see ROTATE subroutine)
!                   whose dimension is natoms, the number of atoms
!  newY(nhelix):       rotated Y coordinates, together with newX to calculate SWEEP
!  newZ_max:       the highest atom of the helix 
!  newZ_min:       the lowest  atom of the helix
! 
!  output==========
!  thetahelix:    scalar of the angle swept by the complete helix
!                  in degrees...
!  thedstep(nhelix-1):array of the angles between successive steps along 
!                  the helix, with dimension natoms-1 (10 atoms have 9 steps)
!                   in degrees...
!  pit:           scalar value of the pitch (for first strand, or only strand)
!--------------------------------------------------------------------------
SUBROUTINE pitchang(natoms, nhelix,strandbstart,x0,y0,newX,newY,newZ_min,newZ_max,theta_helix,thed_step,pit)
IMPLICIT NONE
!--declare passenger variables from MAIN..........
   integer         , intent(in)     :: nhelix, strandbstart,natoms
   double precision, intent(in)     :: x0, y0, newZ_max, newZ_min
   double precision, dimension(:),   &
                     allocatable,    &
                     intent(in)     :: newX, newY
   !-declare passenger variables to MAIN
   double precision, dimension(:),   &
                     allocatable,    &
                     intent(out)     :: thed_step
   double precision, intent(out)     :: theta_helix, pit
!--------------------------------------------------------------------------
!--declare internal variables that MAIN don-t need...
   double precision                  &
        :: vecix,veciy, vecip1x,vecip1y, dotiip1, magip1, &
           dotovermags, zetapitch, kappa, magi, thetastep
   integer                          :: iatom
   DOUBLE PRECISION, PARAMETER  :: pi  = 3.141592653589793
   DOUBLE PRECISION, PARAMETER  :: rtd = 180.0/pi

!==========================================================================

!--allocate the arrays whose dimensions needed user input variable natoms..
   ALLOCATE(thed_step    (natoms-1))
!==========================================================================

thetastep=0.0D+00
theta_helix=0.0D+00
zetapitch=0.0D+00
kappa=0.0D+00
pit=0.0D+00


   DO iatom = strandbstart, nhelix-1
      vecix    = newX(iatom)   - x0                !ith x crd
      veciy    = newY(iatom)   - y0                !ith y crd
      vecip1x  = newX(iatom+1) - x0                !i+1th x crd
      vecip1y  = newY(iatom+1) - y0                !i+1th y crd
       dotiip1 = vecix*vecip1x + veciy*vecip1y     !dot i,i+1
       magi    = SQRT(vecix*vecix + &              !mag of ith vector
                      veciy*veciy)
       magip1  = SQRT(vecip1x*vecip1x + &          !mag of i+1th vector
                      vecip1y*vecip1y)

       ! ALPHA = ARCCOS(A*B/|A||B|), A is ith vec, B is i+1th vec...
       dotovermags        = dotiip1/(magi*magip1) !value holder
       thetastep          = ACOS(dotovermags)     !angle, in rad., bet steps
       thed_step(iatom)   = thetastep*rtd

       theta_helix = theta_helix + thed_step(iatom) !sum over steps, total angle

       !write(*,*) 'in-routine_thed_step(', iatom, ')', thed_step(iatom)

   ENDDO

!===========================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   now calculate the pitch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       zetapitch = ABS (newZ_max - newZ_min)   !height, e_z, of helix
       kappa     = zetapitch * rtd / theta_helix   !this needs to be in rad.
       pit       = ABS ( kappa * 2*pi )

END SUBROUTINE pitchang
