PROGRAM helang
IMPLICIT NONE
DOUBLE PRECISION   :: phi_one,phi_two,theta_one,theta_two
DOUBLE PRECISION   :: psi, psid
DOUBLE PRECISION   :: p1, p2, t1, t2
DOUBLE PRECISION   :: dp,fac1,fac12,fac2,fac22,fac3,numr,sqnumr,quoti,y
CHARACTER(LEN=100) :: in_phione,in_phitwo,in_thetaone,in_thetatwo
REAL, PARAMETER    :: pi = 3.1415926535
REAL, PARAMETER    :: dtr = pi/180

CALL GETARG(1,in_phione)
!WRITE(*,*) 'what is the PHI (degrees) of FIRST vector?'
READ(in_phione,*) phi_one
!READ(*,*) phi_one
p1 = phi_one * dtr

CALL GETARG(2,in_thetaone)
!WRITE(*,*) 'what is the THETA (degrees) of FIRST vector?'
READ(in_thetaone,*) theta_one
!READ(*,*), theta_one
t1 = theta_one * dtr

CALL GETARG(3,in_phitwo)
!WRITE(*,*) 'what is the PHI (degrees) of SECOND vector?'
READ(in_phitwo,*) phi_two
!READ(*,*) phi_two
p2 = phi_two * dtr

CALL GETARG(4,in_thetatwo)
!WRITE(*,*) 'what is the THETA (degrees) of SECOND vector?'
READ(in_thetatwo,*) theta_two
!READ(*,*) theta_two
t2 = theta_two * dtr

dp    = p2 - p1
fac1  = COS(p2)*SIN(dp)
fac12 = fac1*fac1

fac2  = COS(p1)*SIN(p2)-SIN(p1)*COS(p2)*COS(dp)
fac22 = fac2*fac2

fac3  = SIN(p1)*SIN(p2)+COS(p1)*COS(p2)*COS(dp)

numr  = fac12 + fac22
sqnumr= SQRT(numr)

quoti = sqnumr / fac3
y=0.0d+00
psi = ATAN(quoti) !fix this to be ATAN2???
psid = psi/dtr

WRITE (*,*) psid

END PROGRAM helang
