!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
! SUBROUTINE KOPEN:
!
! Purpose: Handle output file io, both the parameters output, which includes
!           the user-supplied input parameters (see KHELIXINPUTS above), and
!            the optimal helical axis PDB file, if the user requests it
!

SUBROUTINE kopen(helixout_io, helixout_name)
IMPLICIT NONE
   integer, intent(in)       :: helixout_io  ! logical unit number
   character(len=100),intent(in) :: helixout_name! file name
!----local variables
   integer                :: ios          ! i/o status variable

  OPEN(UNIT = helixout_io, &
       FILE = helixout_name, &
       IOSTAT = ios)

  IF (ios .ne. 0) THEN

      WRITE(helixout_io, '(/,2x,a,i4,a,a)') 'Unit ', helixout_io, ' Error on OPEN: ', helixout_name
  ENDIF

  return

END SUBROUTINE kopen
