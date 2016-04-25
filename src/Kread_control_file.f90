!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! 
! SUBROUTINE READ_CONTROL_FILE:
!
! Purpose: read in all the user-supplied input parameters
!

SUBROUTINE read_control_file(control_file,inputhelix,num_grid,grid_phi_beg,&
                          grid_phi_end,grid_theta_beg, grid_theta_end, &
                          natoms, nframes, dsDNA, coord_type, oradian,      &
                          all_helix_out, opt_axis_out, print_sing, print_step, &
         print_to_plot,helixout_name, helix_atom_names, opt_Axis_out_name,ktest,&
         autohelix)
IMPLICIT NONE
   character(LEN=100),intent(in)  :: control_file

   !=====  variables for the program to run, read from user input file..
   double precision, intent(out)  :: grid_phi_beg, grid_phi_end,   &
                                     grid_theta_beg, grid_theta_end
   integer, intent(out)           :: natoms, nframes, dsDNA, coord_type, oradian, &
                                     all_helix_out, opt_axis_out,print_sing,&
                                     num_grid,print_step,print_to_plot,ktest,&
                                     autohelix
   character(LEN=100),intent(out) :: inputhelix, helixout_name,opt_Axis_out_name
   character(LEN=4),  intent(out) :: helix_atom_names !len=2; CA

   !---local variables that operate on the control file that need not be passed
   integer                :: ios=0, line=0, pos, kint_test, khelix_error
   character(LEN=100)     :: label, buffer
!-------------------------------------------------------------------------------
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!    namelist variables....
!----------------------------------------
!   input parameters:
!----------------------------------------
!   inpFILE            file of coordinates (or trajectory)
!   num_grid           grid resolution (how many grids phi/theta limits)
!===>num_grid=180
!   natoms             number of atoms (or any Cartesian coordinates...)
!   nframes            how many coordinates in the trajectory
!====>nframes=1                 !default = just one frame
!   grid_phi_beg     where does phi start? (0 deg, probably)
!====>grid_phi_beg=0.0
!   grid_phi_end       where does phi end? (180 deg, probably)
!====>grid_phi_end=180.0
!   grid_theta_beg   where does theta start? (0 deg, probably)
!====>grid_theta_beg=0.0
!   grid_theta_end     where does theta end? (360 deg, probably)
!====>grid_theta_end=360.0
!   dsDNA              do we have double-stranded DNA? =1 if so
!===>dsDNA=0                   !default = NO, just one helix
!   coord_type         input coordinates file format type: 3 are supported:
!===>coord_type=0              !=0 (just x y z info), nframes > 1 Ok
                          !=1 (PDB format), nframes > 1 NOT supported
                          !=2 (Amber crd format), nframes > 1 Ok
                          !    see cpptraj (ambermd.org) for more info :-)
!   opt_axis_out    do we want the helical axis to be printed as a PDB?
!====>opt_axis_out=0         !default = NO
!   opt_Axis_out_name  what name shall we give to the optimal helical axis?
!====>opt_Axis_out_name='khelix.optAxis.pdb'   !default name...
!   all_helix_out      do we want output FOR EACH ROTATION?
!====>all_helix_out=0           !default = NO
!   print_sing         do we want the singular values for the SVD solution?
!====>print_sing=0              !default = NO
!   oradian            do we want output to be in radians rather than deg?
!====>oradian=0                 !default = NO, we want degrees
!   helix_atom_names    what ATOM_NAMES out of the pdb file do we want?
!====>helix_atom_names='P'       !default = phosphorous atom (P)
!   helixout_name      what name shall we give to the kHelix output file?
!helixout_name='khelix.dat'   !default name..
!----------------------------------------
!    output files
!----------------------------------------
!    main_output
!    all_grid_output
!    opt_helix_pdb
!    opt_helAxis_pdb
!==========================================================================

   CALL getarg(1,control_file)
   OPEN(15,FILE=control_file,ACTION='READ',STATUS='OLD')

   DO WHILE ( ios == 0 )
      READ(15, '(A)', IOSTAT=ios) buffer

      IF ( ios == 0 ) THEN
         line = line + 1
         pos = SCAN(buffer, ' ')
         label = buffer(1:pos)
         buffer = buffer(pos+1:)
         select CASE (label)
  
         CASE('inputhelix')
           READ(buffer, *, IOSTAT=ios) inputhelix
           kint_test = LEN_TRIM(inputhelix)
           IF ( kint_test > 100 ) THEN
              khelix_error = 6
              CALL kerr(khelix_error)
           ENDIF

         CASE('helix_atom_names')
           READ(buffer, *, IOSTAT=ios)  helix_atom_names
           kint_test = LEN_TRIM(helix_atom_names)
           IF ( kint_test > 2 ) THEN
              khelix_error = 7
              CALL kerr(khelix_error)
           ENDIF
 
         CASE('num_grid')
           READ(buffer, *, IOSTAT=ios) num_grid

         CASE('natoms')
           READ(buffer, *, IOSTAT=ios) natoms

         CASE('nframes')
           READ(buffer, *, IOSTAT=ios) nframes

         CASE('grid_phi_beg')
           READ(buffer, *, IOSTAT=ios) grid_phi_beg
           IF ( grid_phi_beg < 0 ) THEN
              khelix_error = 8
              CALL kerr(khelix_error)
           ENDIF

         CASE('grid_phi_end')
           READ(buffer, *, IOSTAT=ios) grid_phi_end
           IF ( grid_phi_end > 180 ) THEN
              khelix_error = 9 
              CALL kerr(khelix_error)
           ENDIF

         CASE('grid_theta_beg')
           READ(buffer, *, IOSTAT=ios)   grid_theta_beg
           IF ( grid_theta_beg < 0 ) THEN
              khelix_error = 10
              CALL kerr(khelix_error)
           ENDIF

         CASE('grid_theta_end')
           READ(buffer, *, IOSTAT=ios) grid_theta_end
           IF ( grid_theta_end > 360 ) THEN
              khelix_error = 11
              CALL kerr(khelix_error)
           ENDIF

         CASE('dsDNA')
           READ(buffer, *, IOSTAT=ios) dsDNA
           IF (  ( dsDNA .NE. 0 ) &
           .AND. ( dsDNA .NE. 1 ) ) THEN
              khelix_error = 12
              CALL kerr(khelix_error)
           ENDIF

         CASE('coord_type')
           READ(buffer, *, IOSTAT=ios) coord_type
           IF (  ( coord_type .NE. 0 ) &
           .AND. ( coord_type .NE. 1 ) &
           .AND. ( coord_type .NE. 2 ) ) THEN
              khelix_error = 13
              CALL kerr(khelix_error)
           ENDIF

         CASE('opt_axis_out')
           READ(buffer, *, IOSTAT=ios)  opt_axis_out
           IF (  ( opt_axis_out .NE. 0 ) &
           .AND. ( opt_axis_out .NE. 1 ) ) THEN
              khelix_error = 14
              CALL kerr(khelix_error)
           ENDIF

         CASE('all_helix_out')
           READ(buffer, *, IOSTAT=ios) all_helix_out
           IF (  ( all_helix_out .NE. 0 ) &
           .AND. ( all_helix_out .NE. 1 ) ) THEN
              khelix_error = 15
              CALL kerr(khelix_error)
           ENDIF

         CASE('oradian')
           READ(buffer, *, IOSTAT=ios) oradian
           IF (  ( oradian .NE. 0 ) &
           .AND. ( oradian .NE. 1 ) ) THEN
              khelix_error = 16
              CALL kerr(khelix_error)
           ENDIF

         CASE('helixout_name')
           READ(buffer, *, IOSTAT=ios) helixout_name
           kint_test = LEN_TRIM(helixout_name)
           IF ( kint_test > 100 ) THEN
               khelix_error = 17
               CALL kerr(khelix_error)
           ENDIF

         CASE('opt_Axis_out_name')
           READ(buffer, *, IOSTAT=ios) opt_Axis_out_name
           kint_test = LEN_TRIM(opt_Axis_out_name)
           IF (  ( opt_axis_out .EQ. 1 ) &
           .AND. ( kint_test .EQ. 0 ) ) THEN
               khelix_error=23
               CALL kerr(khelix_error)
           ENDIF

         CASE('print_sing')
           READ(buffer, *, IOSTAT=ios) print_sing
           IF (  ( print_sing .NE. 1 ) &
           .AND. ( print_sing .NE. 0 ) ) THEN
               khelix_error=24
               CALL kerr(khelix_error)
           ENDIF

         CASE('print_step')
           READ(buffer, *, IOSTAT=ios) print_step
           IF (  ( print_step .NE. 1 ) &
           .AND. ( print_step .NE. 0 ) ) THEN
               khelix_error=25
               CALL kerr(khelix_error)
           ENDIF

         CASE('print_to_plot')
           READ(buffer, *, IOSTAT=ios) print_to_plot
           IF (  ( print_to_plot .NE. 1 ) &
           .AND. ( print_to_plot .NE. 0 ) ) THEN
               khelix_error=26
               CALL kerr(khelix_error)
           ENDIF

         CASE('ktest')
           READ(buffer, *, IOSTAT=ios) ktest
           IF (  ( ktest .NE. 1 ) &
           .AND. ( ktest .NE. 0 ) ) THEN
               khelix_error=27
               CALL kerr(khelix_error)
           ENDIF

         CASE('autohelix')
           READ(buffer, *, IOSTAT=ios) autohelix
           IF (  ( autohelix .NE. 1 ) &
           .AND. ( autohelix .NE. 0 ) ) THEN
               khelix_error=28
               CALL kerr(khelix_error)
           ENDIF

         CASE DEFAULT

         END SELECT

      ENDIF

   ENDDO

END SUBROUTINE read_control_file
