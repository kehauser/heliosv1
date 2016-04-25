!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! 
! SUBROUTINE READPDB:
!
! Purpose: read in an Protein Data Bank file format of coordinates
! PDB ATOMIC COORDINATE ENTRY FORMAT VERSION 3.3 www.wwpdb.org/documentation

SUBROUTINE readpdb(autohelix, inputhelix, helix_atom_names, natoms, &
                 & oldX, oldY, oldZ, ca_oldX, ca_oldY, ca_oldZ)
IMPLICIT NONE
    CHARACTER(LEN=100), INTENT(IN)               :: inputhelix
    CHARACTER(LEN=4), INTENT(IN)                 :: helix_atom_names
    INTEGER, INTENT(IN)                          :: natoms, autohelix
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE,   &
                                    INTENT(OUT)  :: oldX, oldY, oldZ
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE,   &
                                    INTENT(OUT)  :: ca_oldX, ca_oldY, ca_oldZ
! --- local_variables
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE  :: ca_X, ca_Y, ca_Z
  INTEGER                                     :: iatom, kint_test,khelix_error
  INTEGER                                     :: katom, strandbreak, count, tatom
  CHARACTER(LEN=6)                            :: pdb_test
  CHARACTER(LEN=3),DIMENSION (:), ALLOCATABLE :: resName,ca_resName
  CHARACTER(LEN=4),DIMENSION (:), ALLOCATABLE :: p_atom_name, ca_atom_name
  INTEGER,         DIMENSION (:), ALLOCATABLE :: p_atom_serial,ca_atom_serial,resN
  INTEGER                                     :: natomP, natomCA, jatom, ca_natom
  DOUBLE PRECISION                            :: dx, dy, dz, dist

!---run for plain job, one helix (or dsDNA) per PDB. For protein-DNA see below
IF ( autohelix == 0 ) THEN

ALLOCATE(oldX(natoms), oldY(natoms), oldZ(natoms) )
ALLOCATE(p_atom_name(natoms), resName(natoms), p_atom_serial(natoms), resN(natoms) )

     OPEN(UNIT=7,file=inputhelix)
     iatom=0
      DO WHILE (pdb_test .NE.'END   ')
          IF (iatom == natoms) THEN   !IF PDB missing END, EXITs infinite loop
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
              READ (7,1)p_atom_serial(iatom),p_atom_name(iatom),resName(iatom),&
                   &    resN(iatom),oldX(iatom),oldY(iatom),oldZ(iatom)
            ENDIF
        ENDIF
      ENDDO
      CLOSE (7)
ENDIF !----end run for plain job
RETURN


!----AUTOHELIX~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    
!-----------------------------------------------------------------
!    First, gather all the P atoms of the DNA.
!-----------------------------------------------------------------
IF ( autohelix == 1 ) THEN

natoms=9999 !HK
helix_atom_names='P' !HK

ALLOCATE(oldX(natoms), oldY(natoms), oldZ(natoms) )
ALLOCATE(p_atom_name(natoms), resName(natoms), p_atom_serial(natoms), resN(natoms) )
ALLOCATE(ca_oldX(natoms), ca_oldY(natoms), ca_oldZ(natoms) )
ALLOCATE(ca_X(natoms), ca_Y(natoms), ca_Z(natoms) )
ALLOCATE(ca_atom_name(natoms), ca_resName(natoms), ca_atom_serial(natoms) )

     OPEN(UNIT=7,file=inputhelix)
     iatom=0
      DO WHILE (pdb_test .NE.'END   ')
          IF (iatom == natoms) THEN   !IF PDB missing END, EXITs infinite loop
             EXIT
          ENDIF
        READ (7,3) pdb_test
        IF ( (pdb_test .EQ.'ATOM  ') .OR. (pdb_test .EQ. 'HETATM') ) THEN
          BACKSPACE (7)
          READ (7,2) pdb_test
            IF (TRIM(ADJUSTL(pdb_test)) .EQ. TRIM(ADJUSTL(helix_atom_names))) THEN
              BACKSPACE(7)
              iatom=iatom+1
              READ (7,1)p_atom_name(iatom),resName(iatom),p_atom_serial(iatom),&
                   &    oldX(iatom),oldY(iatom),oldZ(iatom)
            ENDIF
        ENDIF
      ENDDO
      CLOSE (7)

!---try to catch simple read error
      kint_test = SIZE(PACK(p_atom_serial, p_atom_serial > 0))
      IF ( kint_test .NE. natoms ) THEN
         khelix_error = 27
         CALL kerr( khelix_error )
      ENDIF
!---try to catch structure problem: If missing nucleotides, assume p_atom_serial will be 1:n, n+1,m
    katom=p_atom_serial(1)
    strandbreak=0
    DO iatom=katom,natom
       count=katom+1
       tatom=katom+p_atom_serial
       IF ( count .NE. tatom ) THEN
          strandbreak=strandbreak+1
       ENDIF
    ENDDO
    IF ( strandbreak .GT. 0 ) THEN
      khelix_error = 50
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
              READ (7,1)p_atom_serial(iatom),p_atom_name(iatom),resName(iatom),&
                   &    resN(iatom),oldX(iatom),oldY(iatom),oldZ(iatom)
            ENDIF                                                                                      
        ENDIF                                                                                          
      ENDDO                                                                                            
      CLOSE (7)           

!-----------------------------------------------------------------
!---now load ALL the protein CA atoms, we-ll deal with them later
!-----------------------------------------------------------------
     helix_atom_names='CA' !HK

     OPEN(UNIT=7,file=inputhelix)
     iatom=0
      DO WHILE (pdb_test .NE.'END   ')
          IF (iatom == natoms) THEN   !IF PDB missing END, EXITs infinite loop
             EXIT
          ENDIF
        READ (7,3) pdb_test
        IF ( (pdb_test .EQ.'ATOM  ') .OR. (pdb_test .EQ. 'HETATM') ) THEN
          BACKSPACE (7)
          READ (7,2) pdb_test
            IF (TRIM(ADJUSTL(pdb_test)) .EQ. TRIM(ADJUSTL(helix_atom_names))) THEN
              BACKSPACE(7)
              iatom=iatom+1
              READ (7,1)ca_atom_name(iatom),ca_resName(iatom),ca_atom_serial(iatom),&
                   &    ca_X(iatom),ca_Y(iatom),ca_Z(iatom)
            ENDIF
        ENDIF
      ENDDO                                                                                            
      CLOSE (7)                                                                                        
                                                                                                       
!---try to catch simple read error                                                                     
      kint_test = SIZE(PACK(p_atom_serial, p_atom_serial > 0))                                         
      IF ( kint_test .NE. natoms ) THEN                                                                
         khelix_error = 27                                                                             
         CALL kerr( khelix_error )                                                                     
      ENDIF                                                                                            
!---try to catch structure problem: If missing nucleotides, assume p_atom_serial will be 1:n, n+1,m    
    katom=p_atom_serial(1)                                                                             
    strandbreak=0                                                                                      
    DO iatom=katom,natom                                                                               
       count=katom+1                                                                                   
       tatom=katom+p_atom_serial                                                                       
       IF ( count .NE. tatom ) THEN                                                                    
          strandbreak=strandbreak+1                                                                    
       ENDIF                                                                                           
    ENDDO                                                                                              
    IF ( strandbreak .GT. 0 ) THEN                                                                     
      khelix_error = 50                                                                                
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
              READ (7,1)ca_atom_serial(iatom),ca_atom_name(iatom),ca_resName(iatom),&
                   &    resN(iatom),ca_X(iatom),ca_Y(iatom),ca_Z(iatom)
            ENDIF
        ENDIF
      ENDDO
      CLOSE (7)


ENDIF!---end gather ALL CA atoms...

!!======================================================================
!!======================================================================
!! now retain ONLY CA atoms if they-re within kautocut to ANY P atom...
!!======================================================================
!!======================================================================

   INTEGER :: natomP, natomCA, jatom, ca_natom
   DOUBLE PRECISION :: dx, dy, dz, dist

   !First, get the number of atoms of DNA(P) and protein(CA):
   natomP=SIZE(PACK(p_atom_serial, p_atom_serial > 0))
   natomCA=SIZE(PACK(ca_atom_serial, ca_atom_serial > 0))

   ca_natom=0

   !Now, for each P atom, and each CA atom...
   DO   iatom = 1, natomP
     DO jatom = 1, natomCA
     !...calculate the distance to check if within cutoff...
     dx = oldX(iatom) - ca_X(jatom)
     dy = oldY(iatom) - ca_Y(jatom)
     dz = oldZ(iatom) - ca_Z(jatom)
     dist = SQRT( dx*dx + dy*dy + dz*dz )
     !...and if it is within the cutoff, save the coordinates to a NEW ARRAY  
     IF ( dist < kautocut ) THEN
        ca_natom = ca_natom+1  !indexing problem when ca_old and ca_ index differ
        ca_oldX(ca_natom) = ca_X(jatom)
        ca_oldY(ca_natom) = ca_Y(jatom)
        ca_oldZ(ca_natom) = ca_Z(jatom)
     ENDIF 

     ENDDO
   ENDDO


 1    FORMAT (12x,a4,2x,a3,2x,i4,4x,3f8.3)
 2    FORMAT (11x,a4)
 3    FORMAT (a6)



RETURN
END SUBROUTINE readpdb
