!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! 
! SUBROUTINE PRINTER:
!
! Purpose: print output from kHelix, parameters and helical axis PDB file
!

SUBROUTINE printer(what_to_print, best_params, thed_step, newZ, helixout_name)
IMPLICIT NONE
integer, dimension(17), intent(in) :: what_to_print
character(len=100), intent(in) :: helixout_name
double precision, dimension(:), allocatable, intent(in) :: thed_step, newZ
double precision, dimension(12), intent(out) :: best_params
!-------local variables
double precision, dimension(:), allocatable :: rise
double precision :: best_phi, best_the, best_res, best_rad, pit, &
                    theta_helix, pit2, theta_helix2, sing0, sing1, &
                    sing2, sing3, sqrt_residual
double precision :: sum_rise,sum_twis,sd_rise,sd_twis,va,va2,vb,vb2,&
                    avg_rise,avg_twis,sd_rise2,sd_twis2,avg_rise2,avg_twis2,&
                    avg_rise3,avg_twis3,sd_rise3,sd_twis3,sum_rise3,sum_twis3
integer          :: dsDNA, iatom, natoms, nhelix, print_sing, print_step, &
                    print_to_plot, oradian
integer, parameter            :: helixout_io = 42
DOUBLE PRECISION, PARAMETER   :: pi = 3.141592653589793
DOUBLE PRECISION, PARAMETER   :: rtd = 180.0/pi

! what_to_print(index) = best_params(index); 
!  if what_to_print(index) = 1, print corresponding best_params(index)
best_phi        = best_params(1)
best_the        = best_params(2)
best_res        = best_params(3)
best_rad        = best_params(4)
pit             = best_params(5)
theta_helix     = best_params(6)
pit2            = best_params(7)
theta_helix2    = best_params(8)
sing0           = best_params(9)
sing1           = best_params(10)
sing2           = best_params(11)
sing3           = best_params(12)
! get user parameters (integer) control switches, or just pass simple integers
dsDNA           = what_to_print(7)  
print_sing      = what_to_print(9)  
print_step      = what_to_print(14)
natoms          = what_to_print(15)
print_to_plot   = what_to_print(16)
oradian         = what_to_print(17)

ALLOCATE(rise(natoms))

sqrt_residual = SQRT ( best_res )

IF ( oradian == 1 ) THEN
   best_phi     = best_params(1) / rtd
   best_the     = best_params(2) / rtd
   theta_helix  = best_params(6) / rtd
   theta_helix2 = best_params(8) / rtd
ENDIF

   CALL kopen(helixout_io, helixout_name)

!---------------------------------------------------------------------------
!
!   PRINT WITH PRETTY FORMAT, HARDER TO PLOT THO..
!
!---------------------------------------------------------------------------
IF (  ( dsDNA == 1 ) .AND. ( print_sing == 1 ) .AND. ( print_to_plot == 0 )  ) THEN
   WRITE( helixout_io, 701 ) best_phi, best_the, sqrt_residual, best_rad, pit, & 
                             theta_helix, pit2, theta_helix2, &
                             sing0, sing1, sing2, sing3
ENDIF
IF (  ( dsDNA == 1 ) .AND. ( print_sing == 0 ) .AND. ( print_to_plot == 0 ) ) THEN
   WRITE( helixout_io, 702 ) best_phi, best_the, sqrt_residual, best_rad, pit, & 
                             theta_helix, pit2, theta_helix2
ENDIF
IF (  ( dsDNA == 0 ) .AND. ( print_sing == 1 ) .AND. ( print_to_plot == 0 )  ) THEN
   WRITE( helixout_io, 721 ) best_phi, best_the, sqrt_residual, best_rad, pit, & 
                             theta_helix, &
                             sing0, sing1, sing2, sing3
ENDIF
IF (  ( dsDNA == 0 ) .AND. ( print_sing == 0 ) .AND. ( print_to_plot == 0 )  ) THEN
   WRITE( helixout_io, 722 ) best_phi, best_the, sqrt_residual, best_rad, pit, &
                             theta_helix
ENDIF
IF (  ( dsDNA == 1 ) .AND. ( print_step == 1 ) .AND. ( print_to_plot == 0 )  ) THEN
!
!
!    add oradian for thed_step output....   
!
!
   nhelix = (natoms/2)
   WRITE( helixout_io, 1111 )
   DO iatom = 1, nhelix-1
      WRITE( helixout_io, 730 ) iatom, thed_step(iatom), newZ(iatom)
   ENDDO
   WRITE( helixout_io, *   ) 'Strand B steps......................'
   DO iatom = nhelix+1, natoms-1
      WRITE( helixout_io, 730 ) iatom, thed_step(iatom), newZ(iatom)
   ENDDO
!===============================================================
!do iatom = 1, natoms-1
!write(*,*) 'thed_step(',iatom,')', ' = ', thed_step(iatom)
!end do
!===============================================================

ENDIF
IF (  ( dsDNA == 0 ) .AND. ( print_step == 1 ) .AND. ( print_to_plot == 0 )  ) THEN
   WRITE( helixout_io, 1111 )
   DO iatom = 1, natoms-1
      WRITE( helixout_io, 730 ) iatom, thed_step(iatom), newZ(iatom)
   ENDDO
ENDIF


!---------------------------------------------------------------------------
!
!   Calculate RISE between steps...
!
!---------------------------------------------------------------------------
IF ( dsDNA == 1 ) THEN
   nhelix=(natoms)/2
   DO iatom = 1, nhelix-1
      rise(iatom) = ABS ( newZ(iatom+1) -  newZ(iatom) )
   ENDDO
   DO iatom = nhelix, natoms-1
      rise(iatom) = ABS ( newZ(iatom+1) -  newZ(iatom) )
   ENDDO
ENDIF

IF ( dsDNA == 0 ) THEN
   DO iatom = 1, natoms-1
      rise(iatom) = ABS ( newZ(iatom+1) -  newZ(iatom) )
   ENDDO
ENDIF


!---------------------------------------------------------------------------
!
!   Calculate Statistics of RISE and TWIST for each STEP...
!
!    Sample standard deviation BECAUSE we assume the statistics are general.
!
!---------------------------------------------------------------------------
IF ( dsDNA == 0 ) THEN
   sum_rise = 0
   sum_twis = 0
   DO iatom = 1, natoms-1
      sum_rise = sum_rise + rise(iatom)
      sum_twis = sum_twis + thed_step(iatom)
   ENDDO
   avg_rise = sum_rise / ( natoms-1 )
   avg_twis = sum_twis / ( natoms-1 )
   sum_rise = 0
   sum_twis = 0
   DO iatom = 1, natoms-1
      va = avg_rise - rise(iatom)
      vb = avg_twis - thed_step(iatom)
      va2 = va * va
      vb2 = vb * vb
      sum_rise = sum_rise + va2
      sum_twis = sum_twis + vb2
   ENDDO
   sd_rise = SQRT ( sum_rise / ( natoms-1 ) )
   sd_twis = SQRT ( sum_twis / ( natoms-1 ) ) 
ENDIF

! calculate average and s.d. for the strands separately, and together
IF ( dsDNA == 1 ) THEN
   !average first stand:
   nhelix = (natoms/2)
   sum_rise = 0
   sum_twis = 0
   DO iatom = 1, nhelix-1
      sum_rise = sum_rise + rise(iatom)
      sum_twis = sum_twis + thed_step(iatom)
   ENDDO
   avg_rise = sum_rise / ( nhelix-1 ) 
   avg_twis = sum_twis / ( nhelix-1 )

   sum_rise3=sum_rise !store strand one for average of both, below
   sum_twis3=sum_twis !idem
   sum_rise = 0 !reset for second strand
   sum_twis = 0 !idem

   !average second strand:
   DO iatom = nhelix+1, natoms-1
      sum_rise = sum_rise + rise(iatom)
      sum_twis = sum_twis + thed_step(iatom)
      !get sum for both strands, too:
      sum_rise3= sum_rise3+ rise(iatom)
      sum_twis3= sum_twis3+ thed_step(iatom)
   ENDDO
   avg_rise2= sum_rise / ( nhelix-1 )
   avg_twis2= sum_twis / ( nhelix-1 )
   !average for both strands:
   avg_rise3= sum_rise3/ ( natoms-2 )
   avg_twis3= sum_twis3/ ( natoms-2 )

   !s.d. first strand:
   sum_rise = 0
   sum_twis = 0 
   DO iatom = 1, nhelix-1
      va = avg_rise - rise(iatom)
      vb = avg_twis - thed_step(iatom)
      va2 = va * va
      vb2 = vb * vb
      sum_rise = sum_rise + va2
      sum_twis = sum_twis + vb2
   ENDDO
   sd_rise = SQRT ( sum_rise / ( nhelix-1 ) )
   sd_twis = SQRT ( sum_twis / ( nhelix-1 ) )

   !s.d. second strand:
   sum_rise3= sum_rise
   sum_twis3= sum_twis
   sum_rise = 0
   sum_twis = 0
   DO iatom = nhelix+1, natoms-1
      va = avg_rise2 - rise(iatom)
      vb = avg_twis2 - thed_step(iatom)
      va2 = va * va
      vb2 = vb * vb
      sum_rise = sum_rise + va2
      sum_twis = sum_twis + vb2
      !s.d. for both strands, sum:
      sum_rise3= sum_rise3+ va2
      sum_twis3= sum_twis3+ vb2
   ENDDO
   sd_rise2 = SQRT ( sum_rise / ( nhelix-1 ) )
   sd_twis2 = SQRT ( sum_twis / ( nhelix-1 ) )
   !s.d. both strands:
   sd_rise3 = SQRT ( sum_rise3/ ( natoms-2 ) )
   sd_twis3 = SQRT ( sum_twis3/ ( natoms-2 ) ) 

ENDIF            
!---------------------------------------------------------------------------
!
!   PRINT WITH simple FORMAT, EASIER TO PLOT...
!
!---------------------------------------------------------------------------
IF (  ( dsDNA == 1 ) .AND. ( print_sing == 1 ) .AND. ( print_to_plot == 1 )  ) THEN
   WRITE( helixout_io, 501 ) best_phi, best_the, sqrt_residual, best_rad, pit, &
                             theta_helix, pit2, theta_helix2, &
                             sing0, sing1, sing2, sing3
ENDIF
IF (  ( dsDNA == 1 ) .AND. ( print_sing == 0 ) .AND. ( print_to_plot == 1 ) ) THEN
   WRITE( helixout_io, 502 ) best_phi, best_the, sqrt_residual, best_rad, pit, &
                             theta_helix, pit2, theta_helix2
ENDIF
IF (  ( dsDNA == 0 ) .AND. ( print_sing == 1 ) .AND. ( print_to_plot == 1 )  ) THEN
   WRITE( helixout_io, 521 ) best_phi, best_the, sqrt_residual, best_rad, pit, &
                             theta_helix, &
                             sing0, sing1, sing2, sing3
ENDIF
IF (  ( dsDNA == 0 ) .AND. ( print_sing == 0 ) .AND. ( print_to_plot == 1 )  ) THEN
   WRITE( helixout_io, 522 ) best_phi, best_the, sqrt_residual, best_rad, pit, &
                             theta_helix
ENDIF
IF (  ( dsDNA == 1 ) .AND. ( print_step == 1 ) .AND. ( print_to_plot == 1 )  ) THEN
   nhelix = (natoms/2)
   WRITE( helixout_io, 1111 )
   DO iatom = 1, nhelix-1
      WRITE( helixout_io, 730 ) iatom, thed_step(iatom), rise(iatom), newZ(iatom)
   ENDDO
   WRITE( helixout_io, 740 ) avg_rise, sd_rise
   WRITE( helixout_io, 741 ) avg_twis, sd_twis
      
   WRITE( helixout_io, 1111 )
   DO iatom = nhelix+1, natoms-1
      WRITE( helixout_io, 730 ) iatom, thed_step(iatom), rise(iatom), newZ(iatom)
   ENDDO
   WRITE( helixout_io, 740 ) avg_rise2, sd_rise2
   WRITE( helixout_io, 741 ) avg_twis2, sd_twis2
   WRITE( helixout_io, 752 )
   WRITE( helixout_io, 753 )
   WRITE( helixout_io, 750 ) avg_rise3, sd_rise3
   WRITE( helixout_io, 751 ) avg_twis3, sd_twis3
ENDIF
IF (  ( dsDNA == 0 ) .AND. ( print_step == 1 ) .AND. ( print_to_plot == 1 )  ) THEN
   WRITE( helixout_io, 1111 )
   DO iatom = 1, natoms-1
      WRITE( helixout_io, 730 ) iatom, thed_step(iatom), rise(iatom), newZ(iatom)
   ENDDO
   WRITE( helixout_io, 740 ) avg_rise, sd_rise
   WRITE( helixout_io, 741 ) avg_twis, sd_twis
ENDIF
!--------------------------------------------------------------------------
!
!   DEAL WITH THE INPUT AND OUTPUT FILES
!
!--------------------------------------------------------------------------

!---format output style for helical parameters...
701 format (3x, 'Phi:                ', f10.4, /, &
            3x, 'Theta:              ', f10.4, /, &
            3x, 'Residual:           ', f10.4, /, &
            3x, 'Radius:             ', f10.4, /, &
            3x, 'Pitch of A:         ', f10.4, /, &
            3x, 'Sweep of A:         ', f10.4, /, &
            3x, 'Pitch of B:         ', f10.4, /, &
            3x, 'Sweep of B:         ', f10.4, /, &
            3x, 'Singular Value (0): ', f10.4, /, &
            3x, 'Singular Value (1): ', f10.4, /, &
            3x, 'Singular Value (2): ', f10.4, /, &
            3x, 'Singular Value (3): ', f10.4, /    )
501 format (3x,12f10.4)

702 format (3x, 'Phi:                ', f10.4, /, &
            3x, 'Theta:              ', f10.4, /, &
            3x, 'Residual:           ', f10.4, /, &
            3x, 'Radius:             ', f10.4, /, &
            3x, 'Pitch of A:         ', f10.4, /, &
            3x, 'Sweep of A:         ', f10.4, /, &
            3x, 'Pitch of B:         ', f10.4, /, &
            3x, 'Sweep of B:         ', f10.4, /     )
502 format (3x,8f10.4) 
            
721 format (3x, 'Phi:                ', f10.4, /, &
            3x, 'Theta:              ', f10.4, /, &
            3x, 'Residual:           ', f10.4, /, &
            3x, 'Radius:             ', f10.4, /, &
            3x, 'Pitch of A:         ', f10.4, /, &
            3x, 'Sweep of A:         ', f10.4, /, &
            3x, 'Singular Value (0): ', f10.4, /, &
            3x, 'Singular Value (1): ', f10.4, /, &
            3x, 'Singular Value (2): ', f10.4, /, &
            3x, 'Singular Value (3): ', f10.4, /    )
521 format (3x,10f10.4)

722 format (3x, 'Phi:                ', f10.4, /, &
            3x, 'Theta:              ', f10.4, /, &
            3x, 'Residual:           ', f10.4, /, &
            3x, 'Radius:             ', f10.4, /, &
            3x, 'Pitch of A:         ', f10.4, /, &
            3x, 'Sweep of A:         ', f10.4, /    )
522 format (3x,6f10.4)

730 format (3x, 'Step: ', i4, '    Twist (degrees): ', f10.4, '    Rise: ', f10.4, '    z-coord: ', f10.4, / )

740 format ('Rise:  ', f10.4, ' (avg.) ', f10.4, ' (s.d.)', / )
741 format ('Twist: ', f10.4, ' (avg.) ', f10.4, ' (s.d.)', / )
750 format ('Rise:  ', f10.4, ' (avg.) ', f10.4, ' (s.d.)', / )                                                             
751 format ('Twist: ', f10.4, ' (avg.) ', f10.4, ' (s.d.)', / )    
752 format ('---------------------------------------------------')
753 format ('Average and standard deviation for both DNA strands')

1111 format(/1x, 64('.'), /10x, &
             'Angles swept between successive steps               :::', &
             /1x, 64('.')/)

END SUBROUTINE printer
