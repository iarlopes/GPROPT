! Routine for linear approximation
! Igor Lopes, February 2015
SUBROUTINE POLYLIN(NPOINT,X,F)
    IMPLICIT NONE
    ! Arguments
    INTEGER NPOINT
    REAL(8) :: XMIN, FMIN
    REAL(8),DIMENSION(4) :: X , F
    ! Locals
    REAL(8) :: A0,A1
    ! Begin algorithm
    IF(NPOINT.EQ.1)THEN
        A1=F(2)
        A0=F(1)-F(2)*X(1)
    ELSEIF(NPOINT.EQ.2)THEN
        A1=(F(2)-F(1))/(X(2)-X(1))
        A0=F(1)-A1*X(1)
    ELSE
        WRITE(*,*)'Wrong number of input points in POLYLIN'
        GOTO 99
    ENDIF
    WRITE(*,*)'Approximation coefficients:'
    WRITE(*,*)'a0=',A0
    WRITE(*,*)'a1=',A1
    WRITE(11,*)'Approximation coefficients:'
    WRITE(11,*)'a0=',A0
    WRITE(11,*)'a1=',A1
99 CONTINUE    
END SUBROUTINE