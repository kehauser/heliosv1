!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Helios v1.0. GNU stuff here. 2015. HTH, kevin.
!
! Kevin Hauser, Mosavverul Hassan, Yiging He, Carlos Simmerling
!  and  Evangelos Coutsias
!
! The Louis and Beatrice Laufer Center for Physical and Quantitative Biology
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM HELIOS
IMPLICIT NONE

!---------------------------------------------------------------------------
!
!      INTERFACE BLOCK FOR SUBROUTINES VARIABLE DECLARATIONS...
!
!---------------------------------------------------------------------------
INTERFACE
    SUBROUTINE rotate(natoms, phi_idx, theta_idx, oldX,oldY,oldZ,new_xyz)
      INTEGER,               INTENT(IN) :: natoms
      DOUBLE PRECISION,      INTENT(IN) :: phi_idx, theta_idx
      DOUBLE PRECISION, DIMENSION(:),                                  &
        ALLOCATABLE,         INTENT(IN) :: oldX, oldY, oldZ
      DOUBLE PRECISION, DIMENSION(:,:),                                &
        ALLOCATABLE,        INTENT(OUT) :: new_xyz
    END SUBROUTINE
    SUBROUTINE pitchang(natoms, nhelix, strandbstart, x0,y0,newX,newY, &
                        newZ_min, newZ_max, theta_helix, thed_step, pit)
      INTEGER,               INTENT(IN) :: nhelix, strandbstart, natoms
      DOUBLE PRECISION,      INTENT(IN) :: x0, y0, newZ_min, newZ_max
      DOUBLE PRECISION, DIMENSION(:),                                  &
        ALLOCATABLE,         INTENT(IN) :: newX, newY
      DOUBLE PRECISION, DIMENSION(:),                                  &
        ALLOCATABLE,        INTENT(OUT) :: thed_step
      DOUBLE PRECISION,     INTENT(OUT) :: theta_helix, pit
    END SUBROUTINE
    SUBROUTINE read_control_file(control_file,inputhelix,num_grid,     &
                 grid_phi_beg,grid_phi_end,grid_theta_beg,             &
                 grid_theta_end, natoms, nframes, dsDNA, coord_type,   &
              oradian,all_helix_out,opt_axis_out,print_sing,print_step,&
                 print_to_plot,helixout_name,helix_atom_names,         &
                 opt_Axis_out_name,ktest,autohelix)
      character(LEN=100),intent(in)     :: control_file
      character(LEN=100),intent(out)    :: inputhelix, helixout_name,  &
                                           opt_Axis_out_name
      character(LEN=4),  intent(out)    :: helix_atom_names !len=2; CA
      double precision, intent(out)     :: grid_phi_beg,grid_phi_end,  &
                                           grid_theta_beg,grid_theta_end
      integer, intent(out)              :: num_grid, natoms, nframes,  &
                                           dsDNA, coord_type, oradian, &
                                           all_helix_out,opt_axis_out, &
                    print_sing,print_step, print_to_plot,ktest,autohelix
    END SUBROUTINE

    SUBROUTINE readpdb(autohelix, kautocut, inputhelix, helix_atom_names, natoms, &
                 & oldX, oldY, oldZ, ca_oldX, ca_oldY, ca_oldZ, resName, resN)
      CHARACTER(LEN=100), INTENT(IN)               :: inputhelix
      CHARACTER(LEN=4), INTENT(IN)                 :: helix_atom_names
      INTEGER, INTENT(IN)                          :: natoms, autohelix
      DOUBLE PRECISION,             INTENT(IN)     :: kautocut
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE,   &
                                      INTENT(OUT)  :: oldX, oldY, oldZ
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE,   &
                                      INTENT(OUT)  :: ca_oldX, ca_oldY, ca_oldZ
      CHARACTER(LEN=3),DIMENSION (:),ALLOCATABLE,   &
                                      INTENT(OUT)  :: resName
      INTEGER,         DIMENSION (:),ALLOCATABLE,   & 
                                      INTENT(OUT)  :: resN

    END SUBROUTINE
    SUBROUTINE readtraj(inputhelix, nframes, natoms, all_xyz)
      character(len=100), intent(in)    :: inputhelix
      integer, intent(in)               :: natoms, nframes
      double precision, dimension(:),                                  &
        allocatable, intent(out)        :: all_xyz
    END SUBROUTINE
    SUBROUTINE printer(what_to_print, best_params, thed_step, newZ,    &
                 helixout_name)
      integer, dimension(17), intent(in):: what_to_print
      character(len=100), intent(in)    :: helixout_name
      double precision, dimension(:),                                  &
        allocatable, intent(in)         :: thed_step, newZ
      double precision, dimension(12),                                 &
        intent(out)                     :: best_params
    END SUBROUTINE
    SUBROUTINE kopen(helixout, helixout_name)
      integer, intent(in)               :: helixout
      character(100),intent(in)         :: helixout_name
    END SUBROUTINE
    SUBROUTINE kerr(khelix_error)
      integer, intent(in)               :: khelix_error
    END SUBROUTINE
    SUBROUTINE printer_optHelix(opt_Axis_out_name, natoms, new_xyz, &
                                resName, resN)
      double precision,dimension(:,:),   &
                          intent(in)    :: new_xyz
      character(len=100), intent(in)    :: opt_Axis_out_name
      integer, intent(in)               :: natoms
      CHARACTER(LEN=3),DIMENSION (:),ALLOCATABLE,   &
                                      INTENT(IN)   :: resName
      INTEGER,         DIMENSION (:),ALLOCATABLE,   &
                                      INTENT(IN)   :: resN

    END SUBROUTINE
END INTERFACE
!===========================================================================

!---------------------------------------------------------------------------
!
!      MAIN PROGRAM VARIABLE-S DECLARATION
!
!---------------------------------------------------------------------------
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::                      &
     all_xyz, oldX, oldY, oldZ, newX, newY, newZ, newZA, newZB,        &
     newxyz_BMAT, sing, thed_step, thed_step2, ca_oldX, ca_oldY, ca_oldZ
   INTEGER,         DIMENSION (:), ALLOCATABLE :: lowestRes
   INTEGER,         DIMENSION (17)             :: what_to_print
   DOUBLE PRECISION,DIMENSION (12)             :: best_params
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE ::                      &
     newxyz_AMAT, new_xyz, resmat
   CHARACTER(LEN=100)                          ::                      &
     inputhelix, control_file, helixout_name, opt_Axis_out_name
   INTEGER                                     ::                      & 
     I, min_iphi, min_itheta, all_helix_out, opt_axis_out,             &
     print_sing, print_step, print_to_plot, oradian, dsDNA,            &
     strandbstart, coord_type, num_grid,                               &
     iatom, iPhi, iTheta, jatom, iFrame, FrameShift, nframes, natoms,  &
     nhelix, ierr, khelix_error, helixout_io, ktest, autohelix,resolhelix
   DOUBLE PRECISION                            ::                      &
     dPhi,dTheta, phi_idx,theta_idx, phideg, thedeg, START, FINISH,    &
     x0, y0, x02, y02, rad_svd, residual, best_res,                    &
     best_phi, best_the, grid_phi_beg,grid_phi_end,grid_theta_beg,     &
     grid_theta_end, pit,pit2,theta_helix, theta_helix2, newZ_min,     &
     newZ_max, phi_beg,the_beg, phi_end,the_end,phi_range,the_range,   &
     kautocut
   CHARACTER(LEN=3),DIMENSION (:),ALLOCATABLE :: resName
   INTEGER,         DIMENSION (:),ALLOCATABLE :: resN

   character(8)     :: date
   character(10)    :: time
   character(512)   :: char_tmp_512

   CHARACTER(LEN=4)                            :: helix_atom_names !len=2;CA
   DOUBLE PRECISION, PARAMETER                 :: pi = 3.141592653589793
   DOUBLE PRECISION, PARAMETER                 :: rtd = 180.0/pi
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! for SVD-Lawson.f                                ON  PRBLOCK  U  W
   INTEGER,SAVE,DIMENSION(4)            :: kpvec = (/ 1, 000000, -1,69 /)
   INTEGER,PARAMETER                    :: mx = 3
   CHARACTER(LEN=8), DIMENSION(1)       :: names = (/' '/)
   DOUBLE PRECISION,DIMENSION(mx)       :: d(mx)
   DOUBLE PRECISION,DIMENSION(mx)       :: work(2*mx)
!===========================================================================

!---------------------------------------------------------------------------
!
!      Start timing the program
!
!---------------------------------------------------------------------------
   CALL CPU_TIME(START)
!===========================================================================

!---------------------------------------------------------------------------
!
!      Open the control file. Expected program usage: ./kHelix.o input
!
!---------------------------------------------------------------------------
   CALL getarg(1,control_file)
   OPEN(15,FILE=control_file,ACTION='READ',STATUS='OLD', IOSTAT=ierr)
    IF ( ierr == 1 ) THEN
         khelix_error =1
    ENDIF
!===========================================================================

!---------------------------------------------------------------------------
!
!      Set default control variables. Get overwritten if provided by user.
!
!---------------------------------------------------------------------------
num_grid=90
nframes=1                 !default = just one frame
grid_phi_beg=0.0
grid_phi_end=180.0
grid_theta_beg=0.0
grid_theta_end=360.0
autohelix=0               !default = NO, you know what you-re doing
kautocut=10
dsDNA=0                   !default = NO, just one helix
coord_type=0              !=0 (just x y z info), nframes > 1 Ok
opt_axis_out=0            !default = NO

!!! HK -- add this to read_control_file
resolhelix=100            !number of points used to represent the helical axis
!!!

all_helix_out=0           !default = NO
print_sing=0              !default = NO
print_step=0              !default = NO
print_to_plot=0           !default = NO
oradian=0                 !default = NO, we want degrees
helix_atom_names='P'      !default = phosphorous atom (P)
helixout_name='helios.dat'!default OUTPUT file_name
opt_Axis_out_name='Helix_along_z_now.pdb' !default OUTPUT file_name
ktest=0              !run after build =1, stop program; =0, normal operation
!===========================================================================
   !get the control variables from the user via subroutine read_control_file
   CALL read_control_file(control_file,inputhelix,num_grid,grid_phi_beg,    &
                          grid_phi_end,grid_theta_beg, grid_theta_end,      &
                          natoms, nframes, dsDNA, coord_type, oradian,      &
                        all_helix_out, opt_axis_out, print_sing, print_step,&
        print_to_plot,helixout_name,helix_atom_names,opt_Axis_out_name,ktest,&
        autohelix)
   !now write the control variables into the output file
!----------------------------------------------------------------------------
! this bit is just to test that the program has compiled
!----------------------------------------------------------------------------
IF ( ktest == 1 ) THEN
   WRITE(*,*) '                                                             '
   WRITE(*,*) 'Helios was built.                                            '
   WRITE(*,*) '                                                             '

   CALL getcwd(char_tmp_512)
   WRITE(*,'(a,a)') 'Helios was just run here:', TRIM(char_tmp_512)

   WRITE(*,*) '                                         '
   WRITE(*,*) ' Please try the tests: test_1, test_2, test_3                     '
   WRITE(*,*) '                                         '

!SSSTTTOOOPPP THE PROGRAM HERE....
stop
ENDIF
!===========================================================================
   helixout_io = 42
   CALL kopen(helixout_io,helixout_name)
   WRITE(helixout_io,1000)
   WRITE(helixout_io, 988) ' inputhelix',       inputhelix
   WRITE(helixout_io, 977) ' coord_type',       coord_type
   WRITE(helixout_io, 977) ' num_grd',          num_grid
   WRITE(helixout_io, 977) ' natoms',           natoms
   WRITE(helixout_io, 977) ' nframes',          nframes
   WRITE(helixout_io, 966) ' grid_phi_beg',     grid_phi_beg
   WRITE(helixout_io, 966) ' grid_phi_end',     grid_phi_end
   WRITE(helixout_io, 966) ' grid_theta_beg',   grid_theta_beg
   WRITE(helixout_io, 966) ' grid_theta_end',   grid_theta_end
   WRITE(helixout_io, 977) ' dsDNA',            dsDNA
   WRITE(helixout_io, 977) ' opt_axis_out',     opt_axis_out
   WRITE(helixout_io, 977) ' all_helix_out',    all_helix_out
   WRITE(helixout_io, 977) ' oradian',          oradian
   WRITE(helixout_io, 988) ' helix_atom_names', helix_atom_names
   WRITE(helixout_io, 988) ' helixout_name',    helixout_name 
   WRITE(helixout_io, 977) ' print_sing',       print_sing
   WRITE(helixout_io, 977) ' print_step',       print_step
   WRITE(helixout_io, 977) ' print_to_plot',    print_to_plot
!---------------------------------------------------------------------------
   !PRINT visual buffer
   WRITE(helixout_io,1011)
   !PRINT run details
   CALL date_and_time(DATE=date, TIME=time)
   WRITE(helixout_io,'(12(a),/)') '| Run on ', date(5:6), '/', date(7:8), '/',  &
        date(1:4), ' at ', time(1:2), ':', time(3:4), ':', time(5:6)
   CALL get_command_argument(0, char_tmp_512)
   WRITE(helixout_io,'(a,a)') '| Executable path:   ', TRIM(char_tmp_512)
   CALL getcwd(char_tmp_512)
   WRITE(helixout_io,'(a,a)') '| Working directory: ', TRIM(char_tmp_512)
   !PRINT visual buffer
   WRITE(helixout_io,1012)
IF ( ( print_to_plot == 1 ) .AND. ( dsDNA == 0) )  THEN
   WRITE(helixout_io,5112)
ENDIF
IF ( ( print_to_plot == 1 ) .AND. ( dsDNA == 1) )  THEN
   WRITE(helixout_io,5113)
ENDIF
!---------------------------------------------------------------------------

!============================================================================

!----------------------------------------------------------------------------
!
!  Allocate arrays based on input parameters
!
!----------------------------------------------------------------------------
   ALLOCATE(all_xyz     (3*natoms*nframes)  )
   ALLOCATE(oldX        (natoms)            )
   ALLOCATE(oldY        (natoms)            )
   ALLOCATE(oldZ        (natoms)            )
   ALLOCATE(newX        (natoms)            )
   ALLOCATE(newY        (natoms)            )
   ALLOCATE(newZ        (natoms)            )
   ALLOCATE(newZA       (natoms/2)          )
   ALLOCATE(newZB       (natoms/2)          )
   ALLOCATE(thed_step   (natoms)            )
   ALLOCATE(thed_step2  (natoms)            )
   ALLOCATE(newxyz_AMAT (natoms,3)          )
   ALLOCATE(newxyz_BMAT (natoms)            )
   ALLOCATE(new_xyz     (natoms,3)          )
   ALLOCATE(sing        (3*mx)              )
   ALLOCATE(resmat      (num_grid,num_grid) )   !COSTLIEST PART OF CODE......
   ALLOCATE(lowestRes   (2)                 )
!   ALLOCATE(resN        (natoms)            )
!   ALLOCATE(resName     (natoms)            )
!============================================================================

!----------------------------------------------------------------------------
!
!  GET COORDINATES FROM USER-S FILE...
!
!----------------------------------------------------------------------------
! file of coordinates, and only coordinates, separated by white space only
IF ( coord_type == 0 ) THEN
    OPEN(51,FILE=inputhelix,ACTION='READ',STATUS='OLD', IOSTAT=ierr)
      IF ( ierr .NE. 0 ) THEN
           khelix_error = 2
           CALL kerr(khelix_error)
      ENDIF
    CLOSE(1)
    READ(51,*) all_xyz
ENDIF
! file is in PDB format, only a single frame is currently supported...
IF ( coord_type == 1 ) THEN
     IF ( nframes > 1 ) THEN 
          khelix_error = 3
          CALL kerr(khelix_error)
     ENDIF
     IF ( autohelix == 1 ) THEN
        natoms=9999 !HK
        helix_atom_names='P' !HK
     ENDIF
    CALL readpdb(autohelix, kautocut, inputhelix, helix_atom_names, natoms, &
                 & oldX, oldY, oldZ, ca_oldX, ca_oldY, ca_oldZ, resName, resN)
ENDIF
! file is in amber coordinate format, good for big MD data!
IF ( coord_type == 2 ) THEN
    CALL readtraj(inputhelix, nframes, natoms, all_xyz)
ENDIF
!============================================================================

!----------------------------------------------------------------------------
!
!   MOD-KHELIX
!   Quick checks of user inputs. THIS IS WHERE USER MODIFIES DISK-SAVERS
!
!----------------------------------------------------------------------------
!comment this out if you-re crazy, or REALLY know what you-re doing..
IF ( ( all_helix_out == 1 ) .AND. ( nframes > 1 ) ) THEN
   khelix_error=69
   CALL kerr(khelix_error)
ENDIF

!---dummy-up PIT2 and THETA_HELIX2 __IF__ dsDNA == 0
IF ( dsDNA == 0 ) THEN
    pit2=0.0D+00
    theta_helix2=0.0D+00
ENDIF

!============================================================================

!----------------------------------------------------------------------------
!
!  MAIN PROGRAM LOOP: ONE FRAME AT A TIME
!
!----------------------------------------------------------------------------

DO iFrame = 1,nframes   ! For each frame, for each atom, load XYZs

  IF ( coord_type == 0 .or. coord_type == 2) THEN
    ! For each structure/frame in the input coordinates, build arrays
    FrameShift = ( iFrame -1 ) * ( 3 ) * ( natoms )
    DO iatom = 1,natoms
     oldX(iatom) = all_xyz ( 3*(iatom-1) + 1 + FrameShift )
     oldY(iatom) = all_xyz ( 3*(iatom-1) + 2 + FrameShift )
     oldZ(iatom) = all_xyz ( 3*(iatom-1) + 3 + FrameShift )
    ENDDO
  ENDIF

!  Initialize the rotated components arrays
   DO iatom = 1,natoms
      newX(iatom)=oldX(iatom)
      newY(iatom)=oldY(iatom)
      newZ(iatom)=oldZ(iatom)
   ENDDO
!============================================================================

!----------------------------------------------------------------------------
!                                                                           
!  START LOOP FOR ROTATIONS OVER THE UNIT SPHERE ALONG PHI AND THETA
!
!----------------------------------------------------------------------------
   
   !Convert input values from degrees to radians
   phi_beg     = grid_phi_beg   / rtd
   the_beg     = grid_theta_beg / rtd
   phi_end     = grid_phi_end     / rtd
   the_end    = grid_theta_end   / rtd
   
   !Find the range in Phi and Theta over which to search and for loop control
   phi_range = phi_end  - phi_beg
   the_range = the_end - the_beg

   !Find the step size through Phi and Theta, and for loop control
   !  Remind: Spherical coordinates: 0 < phi < PI ; 0 < theta < 2*PI
   dTheta      = the_range / num_grid
   dPhi        = phi_range / num_grid

 !Now launch the loops over the appropriate ranges in Phi and Theta, 
 ! incremented based on the resolution (step size) based on num_grid
 DO iPhi    = 2,num_grid-1   ! Loop-PHI
  DO iTheta = 2,num_grid-1   ! Loop-THETA

     !Get the Phi and Theta angles, in radians, back from the loop counters
     ! so we can apply the rotation matrices to the input coordinates
      phi_idx   = ( iPhi   * dPhi )   + phi_beg
      theta_idx = ( iTheta * dTheta ) + the_beg
      !Convert phi_idx and theta_idx to degrees for writing
      phideg = phi_idx   * rtd
      thedeg = theta_idx * rtd

!----------------------------------------------------------------------------
!
! ROTATE COORDINATES
!
!----------------------------------------------------------------------------

      CALL rotate (natoms, phi_idx, theta_idx, oldX, oldY, oldZ, new_xyz)

!============================================================================

!----------------------------------------------------------------------------
!
! SOLVE THE LINEAR PROBLEM, A * X = B
!
!----------------------------------------------------------------------------
   
   !B matrix              B = ( x**2 + y**2 )   
   DO iatom = 1,natoms
      newxyz_BMAT(iatom) = new_xyz(iatom,1)*new_xyz(iatom,1)                &
                           + new_xyz(iatom,2)*new_xyz(iatom,2)  
   ENDDO
   !A matrix              A = ( 1 x y )                 
   newxyz_AMAT(:,1) = 1                    
   DO iatom = 1,natoms                     
      newxyz_AMAT(iatom,2) = new_xyz(iatom,1)   
      newxyz_AMAT(iatom,3) = new_xyz(iatom,2)   

      newX(iatom)          = new_xyz(iatom,1)
      newY(iatom)          = new_xyz(iatom,2)
   ENDDO

! FIND THE BEST ESTIMATE OF X USING SVD:
! Call the SVA routine from lawson.f90 to get our parameters
!  See accompanying documentation for the variables below.

   CALL SVA (newxyz_AMAT,natoms,natoms,3,natoms,newxyz_BMAT,sing,kpvec,&
             & names,1,D,work)

! Calculate the radius and circle center from the output  parameters >>>
!  Back-out the actual parameters, R, X0, and Y0 from the fit results
    x0  = newxyz_AMAT(2,3) /2  
    y0  = newxyz_AMAT(3,3) /2 
    x02 = x0*x0
    y02 = y0*y0
    rad_svd = SQRT( x02 + y02 + (newxyz_AMAT(1,3)) )

!----------------------------------------------------------------------------
!
!   CALCULATE THE RESIDUAL: HELIX
!
!----------------------------------------------------------------------------
!  Calculate the residual as per Vageli >>>
   residual = 0.0D+00
   DO iatom = 1,natoms
    residual = residual + ABS( &
                     & (new_xyz(iatom,1)-x0)*(new_xyz(iatom,1)-x0) +                &
                     & (new_xyz(iatom,2)-y0)*(new_xyz(iatom,2)-y0) +                &
                     & (rad_svd)*(rad_svd)                 -                &
                     & (2*rad_svd)*( SQRT (                                 &
                     & (new_xyz(iatom,1)-x0)*(new_xyz(iatom,1)-x0) +                &
                     & (new_xyz(iatom,2)-y0)*(new_xyz(iatom,2)-y0) ) ) )
   ENDDO
   resmat(iPhi:iPhi,iTheta:iTheta) = residual

!----------------------------------------------------------------------------
!
!   PRINT THE GRID OF INFORMATION
!
!----------------------------------------------------------------------------

IF ( all_helix_out == 1 ) THEN
    !--------get the pitch and angle for ONE helix
    IF ( dsDNA == 0 ) THEN
         strandbstart = 1
         newZ_max      = MAXVAL (new_xyz(:,3), natoms)  
         newZ_min      = MINVAL (new_xyz(:,3), natoms)  
         CALL pitchang(natoms, natoms,strandbstart,x0,y0,newX,newY,&
                       & newZ_min,newZ_max,theta_helix,thed_step,pit)
    ENDIF
    !----IF dsDNA, then need pitch and angle for EACH helix
    IF ( dsDNA == 1 ) THEN
         !dsDNA needs the newZ_max and newZ_min for each strand, 
         !  so we first break newZ into TWO...
         DO iatom = 1, natoms/2
            jatom = (natoms/2) + iatom
            newZA(iatom) = new_xyz(iatom,3) 
            newZB(iatom) = new_xyz(jatom,3) 
         ENDDO
         strandbstart  = 1
         nhelix        = natoms/2
         newZ_max      = MAXVAL ( newZA(:), natoms )
         newZ_min      = MINVAL ( newZA(:), natoms )
         CALL pitchang(natoms, nhelix,strandbstart,x0,y0,newX,newY,&
                       & newZ_min,newZ_max,theta_helix,thed_step2,pit)
           DO iatom = 1, natoms/2
            thed_step(iatom) = thed_step2(iatom)
           ENDDO
         strandbstart  = (natoms/2)+1
         nhelix        = natoms
         newZ_max      = MAXVAL ( newZB(:), natoms )
         newZ_min      = MINVAL ( newZB(:), natoms )
         CALL pitchang(natoms, nhelix,strandbstart,x0,y0,newX,newY,&
                       & newZ_min,newZ_max,theta_helix2,thed_step2,pit2)
           DO iatom = 1, natoms/2
            jatom=iatom + (natoms/2)
            thed_step(jatom) = thed_step2(iatom)
           ENDDO
    ENDIF
!===Setup the printable values to send to PRINTER
        best_params(1) = phideg
      what_to_print(1) = 1
        best_params(2) = thedeg
      what_to_print(2) = 1
        best_params(3) = residual
      what_to_print(3) = 1
        best_params(4) = rad_svd
      what_to_print(4) = 1
        best_params(5) = pit
      what_to_print(5) = 1
        best_params(6) = theta_helix
      what_to_print(6) = 1
        best_params(7) = pit2
      what_to_print(7) = dsDNA
        best_params(8) = theta_helix2
      what_to_print(8) = dsDNA
        best_params(9) = sing(4)
      what_to_print(9) = print_sing
        best_params(10)= sing(1)
      what_to_print(10)= print_sing
        best_params(11)= sing(2)
      what_to_print(11)= print_sing
        best_params(12)= sing(3)
      what_to_print(12)= print_sing
      what_to_print(13)= iframe
      what_to_print(14)= print_step
      what_to_print(15)= natoms
      what_to_print(16)= print_to_plot
      what_to_print(17)= oradian

      CALL printer(what_to_print, best_params, thed_step2, newZ, helixout_name)
           helixout_io = 42
      CALL kopen(helixout_io, helixout_name)
ENDIF !! all_helix_out == 1

!----------------------------------------------------------------------------
!
! FINISH LOOP OVER UNIT SPHERE
!
!----------------------------------------------------------------------------

  ENDDO ! iTheta_complete the loop of the phi-s   !
 ENDDO ! iPhi_complete the loop of the theta-s    !

!============================================================================
!============================================================================

!----------------------------------------------------------------------------
!
! FIND THE BEST HELIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!----------------------------------------------------------------------------

   !Find the combination of Phi and Theta (indices of *MAT below) that
   ! yielded the lowest residual, ie the best fit
   lowestRes   = MINLOC ( resmat , MASK = resmat .GT. 0.1D-11)

   !Set the index values as variables for which to use to select parameters
   ! corresponding to the "Best" Phi and Theta rotations.
   min_iphi    = lowestRes(1)
   min_itheta  = lowestRes(2)

   IF ( min_iphi == 0 ) THEN
      khelix_error = 4
      CALL kerr(khelix_error)
   ENDIF
   IF ( min_itheta == 0 ) THEN
      khelix_error = 5
      CALL kerr(khelix_error)
   ENDIF

   !The element_id of resmat with the smallest value (min_iphi,min_itheta)
   !  needs to now be used to return the actual value of that residual
   best_res = SQRT ( resmat (min_iphi,min_itheta) ) !do sqrt only at end...

!----------------------------------------------------------------------------
!
! OBTAIN THE PARAMETERS AND PROPERTIES FROM THE kHELIXXXXXXXXXXXXXXXXXXXXXXX
!
!----------------------------------------------------------------------------

!Now we use the min_iphi and min_itheta to calculate the helical parameters
! Originally, we had calculated ALL helical parameters WITHIN the above 
!  loop, stored EACH parameter into a HUGE matrix (like resmat!), then used
!  min_iphi,min_itheta to pluck the parameters from those matrices. Here, we
!  only build one HUGE matrix, resmat, then use min_iphi and min_itheta to 
!  recover the parameters to rotate the helix to the Optimal orientation..

! first, rotate the coordinate to (min_iphi,min_itheta)
   !convert element_id to actual angle
   phi_idx   = ( min_iphi   * dPhi   ) + phi_beg  
   theta_idx = ( min_itheta * dTheta ) + the_beg

   CALL rotate(natoms, phi_idx, theta_idx, oldX, oldY, oldZ, new_xyz)

!now build the new{x,y,z} matices (from new_xyz) to send over to sva (again)
DO iatom = 1,natoms
   newX(iatom) = new_xyz(iatom,1)
   newY(iatom) = new_xyz(iatom,2)
   newZ(iatom) = new_xyz(iatom,3)
ENDDO

! second, rerun SVA to get the RADIUS, X0 and y0 
   !B matrix
   DO iatom = 1,natoms
      newxyz_BMAT(iatom) = new_xyz(iatom,1)*new_xyz(iatom,1) &
                        & + new_xyz(iatom,2)*new_xyz(iatom,2)
   ENDDO
   !A matrix
   newxyz_AMAT(:,1) = 1
   DO iatom = 1,natoms
      newxyz_AMAT(iatom,2) = new_xyz(iatom,1)
      newxyz_AMAT(iatom,3) = new_xyz(iatom,2)

      newX(iatom)          = new_xyz(iatom,1)
      newY(iatom)          = new_xyz(iatom,2)
   ENDDO

   CALL SVA (newxyz_AMAT,natoms,natoms,3,natoms,newxyz_BMAT,sing,kpvec,&
             & names,1,d,work)

    x0  = newxyz_AMAT(2,3) /2
    y0  = newxyz_AMAT(3,3) /2
    x02 = x0*x0
    y02 = y0*y0
    rad_svd = SQRT( x02 + y02 + (newxyz_AMAT(1,3)) )

! third, we calculate the pitch and angles of the best_rotated helix
IF ( dsDNA == 0 ) THEN
    nhelix = natoms
    strandbstart = 1
    !find the iatom with largest z-component value
    newZ_max      = MAXVAL ( new_xyz(:,3), natoms )
    !find the iatom with smallest z-component value
    newZ_min      = MINVAL ( new_xyz(:,3), natoms )
    CALL pitchang(natoms, nhelix,strandbstart,x0,y0,newX,newY,&
                  &newZ_min,newZ_max,theta_helix,thed_step,pit)
ENDIF

IF ( dsDNA == 1 ) THEN
    !dsDNA needs the newZ_max and newZ_min for each strand, 
    !  so we first break newZ into TWO...
    DO iatom = 1, natoms/2
       jatom = (natoms/2) + iatom 
       newZA(iatom) = new_xyz(iatom,3) ! get the first helix
       newZB(iatom) = new_xyz(jatom,3) ! get the second helix
    ENDDO

    newZ(:) = new_xyz(:,3)

    !now CALL pitchang to get the pitch and angles for the first helix
    strandbstart  = 1
    nhelix        = (natoms/2)
    !find the iatom with largest z-component value
    newZ_max      = MAXVAL ( newZA(:), natoms/2 )
    !find the iatom with largest z-component value
    newZ_min      = MINVAL ( newZA(:), natoms/2 )
    CALL pitchang(natoms,nhelix,strandbstart,x0,y0,newX,newY,&
                  & newZ_min,newZ_max,theta_helix,thed_step2,pit)
         DO iatom = 1, natoms/2-1
          thed_step(iatom) = thed_step2(iatom)
         ENDDO
    !now CALL pitchang to get the pitch and angles for the second helix 
    !  NOTE: strandBstart!
    strandbstart  = (natoms/2) + 1
    nhelix        = natoms
    !find the iatom with largest z-component value
    newZ_max      = MAXVAL ( newZB(:), natoms/2 )
    !find the iatom with largest z-component value
    newZ_min      = MINVAL ( newZB(:), natoms/2 )
    CALL pitchang(natoms,nhelix,strandbstart,x0,y0,newX,newY,&
                  & newZ_min,newZ_max,theta_helix2,thed_step2,pit2)
         DO iatom = 2, natoms/2
            jatom = iatom + natoms/2 -1
            thed_step(jatom) = thed_step2(jatom)
         ENDDO
ENDIF
 
!----------------------------------------------------------------------------
!
! WRITING CONTROL for final (default) output
!
!----------------------------------------------------------------------------

   phi_idx   = ( min_iphi   * dPhi   ) + phi_beg
   theta_idx = ( min_itheta * dTheta ) + the_beg
   best_phi = 0.0D+00  !initialize
   best_the = 0.0D+00  !initialize 
   IF ( oradian == 0 ) THEN
      best_phi = phi_idx
      best_the = theta_idx
   ENDIF

! what_to_print(index) = best_params(index); 
!  if what_to_print(index) = 1, print corresponding best_params(index)
IF ( opt_axis_out == 0 ) THEN
        best_params(1) = best_phi*rtd
      what_to_print(1) = 1
        best_params(2) = best_the*rtd
      what_to_print(2) = 1
        best_params(3) = best_res
      what_to_print(3) = 1
        best_params(4) = rad_svd
      what_to_print(4) = 1
        best_params(5) = pit
      what_to_print(5) = 1
        best_params(6) = theta_helix
      what_to_print(6) = 1
        best_params(7) = pit2
      what_to_print(7) = dsDNA
        best_params(8) = theta_helix2
      what_to_print(8) = dsDNA
        best_params(9) = sing(4)
      what_to_print(9) = print_sing
        best_params(10)= sing(1)
      what_to_print(10)= print_sing
        best_params(11)= sing(2)
      what_to_print(11)= print_sing
        best_params(12)= sing(3)
      what_to_print(12)= print_sing
      what_to_print(13)= iframe
      what_to_print(14)= print_step
      what_to_print(15)= natoms
      what_to_print(16)= print_to_plot
      what_to_print(17)= oradian
      CALL printer(what_to_print, best_params, thed_step, newZ, helixout_name)
          helixout_io = 42
      CALL kopen(helixout_io, helixout_name)
ENDIF
!===========================================================================

!----------------------------------------------------------------------------
!
!  PRINT the points in the NEW helical frame (so helical axis along Z)
!
!----------------------------------------------------------------------------
IF ( opt_axis_out == 1 ) THEN
!
!!!!!!!
      CALL printer_optHelix(opt_Axis_out_name, natoms, new_xyz, resName, resN)
!!!!!!!
!
ENDIF
!===========================================================================

!----------------------------------------------------------------------------
!
!  ZERO-OUT ARRAYS AND VARIABLES
!
!----------------------------------------------------------------------------
!  Clear variables that are summed within an iPhi iTheta loop
   min_iphi   = 0  !V -1
   min_itheta = 0  !V -1

!  Clear the arrays
   DO I = 1, num_grid-2
      resmat(I,:)  = 0.0D+00
      resmat(:,I)  = 0.0D+00
   ENDDO

!  Clean up the coordinate arrays
   DO iatom = 1,natoms
      newX(iatom)          = 0.0D+00
      newY(iatom)          = 0.0D+00
      newZ(iatom)          = 0.0D+00
      newxyz_AMAT(iatom,3) = 0.0D+00
      newxyz_AMAT(iatom,2) = 0.0D+00
      newxyz_AMAT(iatom,1) = 0.0D+00
      newxyz_BMAT(iatom)   = 0.0D+00
   ENDDO

!  End of the iFrame loop, so we-re done with the trajectory!
ENDDO

!----------------------------------------------------------------------------
!
!  FORMAT CENTER
!
!----------------------------------------------------------------------------

!  Format for output of control variables...
988 format(4x,a16,12x,a,4x)
977 format(4x,a16,12x,i16,4x)
966 format(4x,a16,20x,f8.4,4x)
1000 format(/10x, 55('~'), /10x, &
             'kHelix  v1.0                                     2015', &
             /10x, 55('~')/)
1011 format(/1x, 64('~'), /10x, &
             'kHelix  execution info                              :::', &
             /1x, 64('~')/)
1012 format(/1x, 64('~'), /10x, &
             'kHelix  results                                     :::', &
             /1x, 64('~')/)
1200 format(/10x, 55('~'), /10x, &
             'kHelix  execution timing:                         :::', &
             /10x, 55('~')/)
1201 format(3x, a, 4x, f16.6, 4x)

5112 format(/1x, 64('.'), /3x, &
'       Phi     Theta  Residual    Radius    Pitch1    Sweep1', &
             /1x, 64('.')/)

5113 format(/1x, 84('.'), /3x, &
'       Phi     Theta  Residual    Radius    Pitch1    Sweep1    Pitch2    Sweep2', &
             /1x, 84('.')/)


!===========================================================================

! CPU TIMING
   helixout_io=42
   CALL CPU_TIME(FINISH)
   WRITE(helixout_io,1200)
   WRITE(helixout_io,1201) 'Duration of Helios execution (seconds):',&
                            & finish-start

! Close the output file...
   CLOSE(helixout_io)

! End program..
END PROGRAM HELIOS
