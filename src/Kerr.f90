!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
! SUBROUTINE KERR:
!
! Purpose: Handling the "most obvious" potential errors the user may encounter
!           with "helpful and not too rude" messages. See Adam Smith-s definition
!            of "rude", or that of the OED...
!

SUBROUTINE kerr(khelix_error)
IMPLICIT NONE
 integer, intent(in) :: khelix_error

!===alright, this sucks, let the user know something broke...
IF ( khelix_error .ne. 0 ) THEN
   WRITE(*,*) '===================================================================='
   WRITE(*,*) '=                                                                  ='
   WRITE(*,*) '=             kHelix exiting due to error(s)...                    ='
   WRITE(*,*) '=                                                                  ='
   WRITE(*,*) '=             Please check the error messages.                     ='
   WRITE(*,*) '=             Please check the documentation to ensure you set     ='
   WRITE(*,*) '=              all the input variables correctly, and that you     ='
   WRITE(*,*) '=               have provided files with the proper formats...     ='
   WRITE(*,*) '=                                                                  ='
   WRITE(*,*) '===================================================================='
ENDIF

!---khelix_error =1 problem opening the control file
IF ( khelix_error == 1 ) THEN
   WRITE(*,*) 'There was an error opening the input file! &
             & Make sure it exists, it follows the correct form, &
             & and that you-ve executed the program correctly:   &
             &  $KHELIXHOME/bin/kHelix.o input'
!killer error
STOP
ENDIF

!---khelix_error = 2 problem opening the coordinate file (coord_type = 0)
IF ( khelix_error == 2 ) THEN
  WRITE(*,*) 'There was an error opening the coordinate file! &
            &  coord_type is set to 0, meaning the file MUST contain &
            &  only x, y, and z double precision info... &
            &  If you have a complex file, try using awk (bash) to &
            &  get your data into a simple form like this: &
            &  1.1 1.22 1.333 2.11 2.222 2.333 3.1415 3.141592 3.1415926 &
            &  where there are THREE atoms in this example... YOU CAN HAVE &
            &  as many atoms as you want, but they must be in this format!'
!killer error
STOP
ENDIF
!---khelix_error = 3 input parameter mismatch..
IF ( khelix_error == 3 ) THEN
   WRITE(*,*) 'You-ve selected incompatible options coord_type=1 and nframes > 1&
             &  ...sorry, but kHelix only supports PDB format input coordinates&
             &  with one frame.. set nframes = 1, please!'
!killer error
STOP
ENDIF
!---khelix_error = 4 min_iphi == 0
IF ( khelix_error == 4 ) THEN
   WRITE(*,*) 'An error occured in the grid search. The variable min_iphi is 0. &
             &  This probably means your grid bounds are not valid'
!killer error
STOP
ENDIF
!---khelix_error = 5 min_itheta == 0
IF ( khelix_error == 5 ) THEN
   WRITE(*,*) 'An error occured in the grid search. The variable min_itheta is 0. &
             &  This probably means your grid bounds are not valid'
!killer error
STOP
ENDIF
!---khelix_error = 6
IF ( khelix_error == 6 ) THEN
   WRITE(*,*) 'The file name of your input coordinate file isn-t valid... it-s too long. &
             &  Try simply changing the file name to something with less than 100 characters,&
             &  please'
!killer error
STOP
ENDIF
!---khelix_error = 7 atom names with more than 2 characters were supplied
IF ( khelix_error == 7 ) THEN
   WRITE(*,*) 'Sorry user, kHelix only supports atom names with two characters:&
             &  try CA for protein, P or C1 for DNA, or X for LHC stuff...'
!killer error
STOP
ENDIF
!---khelix_error = 8 grid_phi_beg < 0...
IF ( khelix_error == 8 ) THEN
   WRITE(*,*) 'Whoops, you just set grid_phi_beg to &
             &  a value less than zero. spherical coordinates don-t work like that...'
!killer error
STOP
ENDIF
!---khelix_error = 9 grid_phi_end > 180...
IF ( khelix_error == 9 ) THEN
   WRITE(*,*) 'Whoops, you just set grid_phi_end to &
             &  a value greater than 180. spherical coordinates don-t work like that...'
!killer error
STOP
ENDIF

!---khelix_error = 10 grid_theta_beg < 0 ...
IF ( khelix_error == 10 ) THEN
   WRITE(*,*) 'Whoops, you just set grid_theta_beg to &
             &  a value less than 0. spherical coordinates don-t work like that...'
!killer error
STOP
ENDIF

!---khelix_error = 11 grid_theta_end > 360
IF ( khelix_error == 11 ) THEN
   WRITE(*,*) 'Whoops, you just set grid_theta_end to &
             &  a value greater than 360. spherical coordinates don-t work like that...'
!killer error
STOP
ENDIF

!---khelix_error = 12 dsDNA is not equal to 0 (ssDNA) or 1 (dsDNA)
IF ( khelix_error == 12 ) THEN
   WRITE(*,*) 'Sorry user, but kHelix only supports single-stranded DNA &
             &  dsDNA = 0, or double-straneded DNA, dsDNA = 1...'
!killer error
STOP
ENDIF

!---khelix_error = 13 coord_type not equal to 0 or 1
IF ( khelix_error == 13 ) THEN
   WRITE(*,*) 'Sorry user, but kHelix only supports 3 types of input files -- &
             &  coord_type must be equal to 0 (simple x y z), to 1 (PDB file), &
             &  or 2 (amber trajectory file WITH NOT BOX INFO!)...'
!killer error
STOP
ENDIF

!---khelix_error = 14 opt_axis_out not equal to 0 or 1
IF ( khelix_error == 14 ) THEN
   WRITE(*,*) 'There are only two valid values of opt_axis_out, 0 (no print) &
             & or 1 (yes print)'
!killer error
STOP
ENDIF

!---khelix_error = 15 all_helix_out not equal to 0 or 1
IF ( khelix_error == 15 ) THEN
   WRITE(*,*) 'there-s only Two valid values of all_helix_out, 0 (print only the one&
              & solution, or 1 (print the solutions for the whole grid!!!!)'
!killer error
STOP
ENDIF

!---khelix_error = 16 oradian not equal to 0 or 1
IF ( khelix_error == 16 ) THEN
   WRITE(*,*) 'There-s only two valid values of oradian, 0 (print degrees), or 1&
              & (print radians)'
!killer error
STOP
ENDIF

!---khelix_error = 17 helixout_name is longer than 100 characters..
IF ( khelix_error == 17 ) THEN
   WRITE(*,*) 'The file name you gave for helixout_name is greater than, &
              & 100 characters. Make it shorter..'
!killer error
STOP
ENDIF

!---khelix_error == 27 more atoms are present in the PDB file than predicted by natoms...
IF ( khelix_error == 27 ) THEN
   WRITE(*,*) 'The value for natoms MUST exactly match the number of atoms matching &
              & helix_atom_name in the PDB file. &
              & Try this: &
              & grep " P " myfile.pdb | wc  and of course change " P " to whatever atom'
!killer error
STOP
ENDIF

!---khelix_error=23, opt_axis_out was request without giving a filename, so default will be used
IF ( khelix_error == 23 ) THEN
   WRITE(*,*) 'You have requested that the optimal helical axis be printed without having&
              & provided a filename for it, so the default name will be used: khelix.optAxis.pdb'
ENDIF

!---khelix_error=24, print_sing provided by user was value other than 0 or 1
IF ( khelix_error == 24 ) THEN
   WRITE(*,*) 'Invalid value of print_sing provided. Using default of 0, do not print them &
              & To print the singular values, set print_sing to 1.'
ENDIF

!---khelix_error=25, print_step provided by user was value other than 0 or 1
IF ( khelix_error == 25 ) THEN
   WRITE(*,*) 'Invalid value of print_step provided. Using default of 0, do not print them &
              & To print the angles of each helix step, set print_step to 1.'
ENDIF

!---khelix_error=26, print_to_plot provided by user was value other than 0 or 1
IF ( khelix_error == 26 ) THEN
   WRITE(*,*) 'Invalid value of print_to_plot provided. Using default of 0, print pretty &
              & To print output data in format easier to plot, set print_to_plot to 1. &
              & I-m defaulting to 0, print pretty.'
ENDIF

!---khelix_error=27, print_to_plot provided by user was value other than 0 or 1
IF ( khelix_error == 27 ) THEN
   WRITE(*,*) 'Invalid value of ktest provided. =0, kHelix runs normally, =1,&
              & kHelix runs only to validate the build. '
ENDIF

!---khelix_error=50, readpdb error, more than 1 strandbreak - missing residues
IF ( khelix_error == 50 ) THEN
   WRITE(*,*) 'Error reading DNA coordinates: seems that there-s more than one&
              & strand break. Either you-re missing a nucleotide in one of your&
              & strands of DNA, or the index of nucleotides are not contiguous.&
              & The nucleotide index of BOTH strands must be contiguous!       &
              & e.g. strand A, 300 to 330 and strand B, 331 to 362             &
              & This is a fatal error. Sorry...'
!killer error
STOP
ENDIF

!---khelix_error = 60, amber coordinate file error - bad format
IF ( khelix_error == 60 ) THEN
   WRITE(*,*) 'The amber coordinate file provided has an unexpected format. The expected format&
              & of the file is that ONLY the first line of the file has header info. The remainder&
              & of the file MUST have only real values ordered by atom index-x, y, and z.'
ENDIF

!---khelix_error=69, user requested all_helix_out with nframe > 1. 
IF ( khelix_error == 69 ) THEN
   WRITE(*,*) 'kHelix disk-saver: all_helix_out = 1 .AND. nframe > 1. See MOD-KHELIX section of &
             & code and modify it appropriately, IF YOU KNOW WHAT YOU-RE DOING.'
!killer error
stop
ENDIF

END SUBROUTINE kerr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@ end kHelix  @@@@@@ @@@ @@ @ @@@@@ @@@ @@ @ @@@@@ @@@ @@ @ @@@@@  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
