!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! 
! SUBROUTINE READPDB:
!
! Purpose: read in an Protein Data Bank file format of coordinates
!

SUBROUTINE readpdb(inputhelix, helix_atom_names, natoms, oldX, oldY, oldZ)
  ! PDB ATOMIC COORDINATE ENTRY FORMAT VERSION 3.3 www.wwpdb.org/documentation
IMPLICIT NONE
    CHARACTER(LEN=100), INTENT(IN)               :: inputhelix  !LEN=72
    CHARACTER(LEN=4), INTENT(IN)                 :: helix_atom_names !LEN=2; CA
    INTEGER, INTENT(IN)                          :: natoms
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE,   &
                                    INTENT(OUT)  :: oldX, oldY, oldZ
! --- local_variables
  INTEGER                                     :: iatom, kint_test,khelix_error
  CHARACTER(LEN=6)                            :: pdb_test  !LEN=4; ATOM
  CHARACTER(LEN=3),DIMENSION (:), ALLOCATABLE :: resName
  CHARACTER(LEN=4),DIMENSION (:), ALLOCATABLE :: p_atom_name !LEN=2; CA
  INTEGER,         DIMENSION (:), ALLOCATABLE :: p_atom_serial,resN

ALLOCATE(oldX(natoms), oldY(natoms), oldZ(natoms) )
ALLOCATE(p_atom_name(natoms), resName(natoms), p_atom_serial(natoms), resN(natoms) )

     OPEN(UNIT=7,file=inputhelix)
     iatom=0
      DO WHILE (pdb_test .NE.'END   ')
   IF (iatom == natoms) THEN
   EXIT
   ENDIF
        READ (7,3) pdb_test
        IF ( (pdb_test .EQ.'ATOM  ') .OR. (pdb_test .EQ. 'HETATM') ) THEN
          BACKSPACE (7)
          READ (7,2) pdb_test
            IF (TRIM(ADJUSTL(pdb_test)) .EQ. TRIM(ADJUSTL(helix_atom_names))) THEN
              BACKSPACE(7)
              iatom=iatom+1
              READ (7,1)p_atom_name(iatom),resName(iatom),p_atom_serial(iatom),oldX(iatom),oldY(iatom),oldZ(iatom)
            ENDIF
        ENDIF
      ENDDO
      CLOSE (7)
!---try to catch an error here..
      kint_test = SIZE(PACK(p_atom_serial, p_atom_serial > 0))
      IF ( kint_test .NE. natoms ) THEN
         khelix_error = 27
         CALL kerr( khelix_error )
      ENDIF
!---IF PASS TEST, THEN do while read...
      DO WHILE (pdb_test .NE.'END   ')
   IF (iatom == natoms) THEN
   EXIT
   ENDIF
        READ (7,3) pdb_test
        IF ( (pdb_test .EQ.'ATOM  ') .OR. (pdb_test .EQ. 'HETATM') ) THEN
          BACKSPACE (7)
          READ (7,2) pdb_test
            IF (pdb_test .EQ. helix_atom_names) THEN
              BACKSPACE(7)
              iatom=iatom+1
              READ (7,1)p_atom_serial(iatom),p_atom_name(iatom),resName(iatom),resN(iatom),oldX(iatom),oldY(iatom),oldZ(iatom)
            ENDIF
        ENDIF
      ENDDO
      CLOSE (7)

 1    FORMAT (12x,a4,2x,a3,2x,i4,4x,3f8.3)
 2    FORMAT (11x,a4)
 3    FORMAT (a6)

RETURN
END SUBROUTINE readpdb
