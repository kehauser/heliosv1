!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! 
! SUBROUTINE READTRAJ:
!
! Purpose: read in an amber coordinate file formatted trajectory
! 

SUBROUTINE readtraj(inputhelix, nframes, natoms, all_xyz)
IMPLICIT NONE
     character(len=100), intent(in)               :: inputhelix
     integer, intent(in)                          :: natoms, nframes
     double precision, dimension(:),allocatable,   &
              intent(out)                         :: all_xyz

     integer    :: iskip, ierr, khelix_error
ALLOCATE(all_xyz  ( 3* natoms * nframes ) )

OPEN(UNIT=7,file=inputhelix)

DO 50 iskip = 1,1 !skip first line, then go to 50
     READ(7,*)
50 CONTINUE
     READ(7,*,IOSTAT=ierr) all_xyz
     IF ( ierr .ne. 0) THEN
     WRITE(*,*) 'test of IOSTAT:', ierr
     khelix_error = 60
     CALL kerr(khelix_error)
     ENDIF
CLOSE(7)

END SUBROUTINE readtraj
