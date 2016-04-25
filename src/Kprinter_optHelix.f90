!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
! SUBROUTINE PRINTER_OPTHELIX
!
! Purpose: make a PDB file that we can use to visualize the helical axis
!          that kHelix thinks it has found..
!

SUBROUTINE printer_optHelix(opt_Axis_out_name, natoms, new_xyz, resName, resN)
IMPLICIT NONE
double precision,dimension(:,:),allocatable,intent(in) :: new_xyz
character(len=100), intent(in) :: opt_Axis_out_name
integer, intent(in) :: natoms
CHARACTER(LEN=3),DIMENSION (:),ALLOCATABLE,   &
                                INTENT(IN)  :: resName
INTEGER,         DIMENSION (:),ALLOCATABLE,   &
                                INTENT(IN)  :: resN
character(len=3) :: racecar
!---local variables
   integer            :: iatom, resNumb, iframe
   integer, parameter :: opt_Axis_out_io = 43

   iframe = 1 ! HK -- just in case some Helios SuperUser wishes to do GPL stuff...

CALL kopen(opt_Axis_out_io, opt_Axis_out_name)

    WRITE(opt_Axis_out_io,FMT=44) 'MODEL ', iframe  !Print the MODEL Number

    DO iatom = 1, natoms
       racecar = resName(iatom)
       resNumb = resN   (iatom)

       WRITE(opt_Axis_out_io,FMT=43)  'HETATM',          & ! A6
                                       iatom,            & ! I5   , 1X
                                      '  HX',            & ! A4
                                      ' ',               & ! A1
                                      'kHx',          & ! A3
                                      ' ',               & ! A1   , 1X
                                      '111',          &
                                       new_xyz(iatom,1), & ! F8.3 , 1X
                                       new_xyz(iatom,2), & ! F8.3 , 1X
                                       new_xyz(iatom,3), & ! F8.3 , 1X
                                      '1.000',           & ! A4
                                      '      ',          & ! A6   , 10X
                                      '    ',            & ! A2
                                      '  '                 ! A2
    END DO
    WRITE (opt_Axis_out_io,FMT=45) 'ENDMDL'
  !Format for the ATOMic properties, including coordinates and RES as BETA 
43 format ( A6,I5,1X,A4,A1,A3,1X,A1,A3,2X,F8.3,1X,F8.3,1X,F8.3,1X,A4,A6,10X,A2,A2 )
    !Append at the end of rotation frame ENDMDL to signify end of frame to VMD

    !Format for the MODEL number information
44  format(A6,1X,I4)
    !Format for the ENDMDL append information
45  format(A6)

    CLOSE(opt_Axis_out_io)

END SUBROUTINE printer_optHelix
