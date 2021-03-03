! Routine for quadratic approximation
! Igor Lopes, February 2015
SUBROUTINE POLYQUAD(NPOINT,X,F,XMIN,FMIN)
    IMPLICIT NONE
    REAL(8) R0 /0.0D0/
    REAL(8) R2 /2.0D0/
    ! Arguments
    INTEGER NPOINT
    REAL(8) :: XMIN, FMIN
    REAL(8),DIMENSION(4) :: X , F
    ! Locals
    REAL(8) :: A0,A1,A2
    ! Begin algorithm
    IF(NPOINT.EQ.2)THEN
        A2=((F(2)-F(1))/(X(2)-X(1))-F(3))/(X(2)-X(1))
        A1=F(3)-R2*A2*X(1)
        A0=F(1)-A1*X(1)-A2*X(1)*X(1)
    ELSEIF(NPOINT.EQ.3)THEN
        A2=((F(3)-F(1))/(X(3)-X(1))-(F(2)-F(1))/(X(2)-X(1)))/(X(3)-X(2))
        A1=(F(2)-F(1))/(X(2)-X(1))-A2*(X(1)+X(2))
        A0=F(1)-A1*X(1)-A2*X(1)*X(1)
    ELSE
        WRITE(*,*)'Wrong number of input points in POLYQUAD'
        GOTO 99
    ENDIF
    WRITE(*,*)'Approximation coefficients:'
    WRITE(*,*)'a0=',A0
    WRITE(*,*)'a1=',A1
    WRITE(*,*)'a2=',A2
    WRITE(11,*)'Approximation coefficients:'
    WRITE(11,*)'a0=',A0
    WRITE(11,*)'a1=',A1
    WRITE(11,*)'a2=',A2
    ! Minimum or Maximum
    IF(A2.LT.R0)THEN
        WRITE(11,*)'A maximum is found for'
        WRITE(*,*)'A maximum is found for'
    ELSEIF(A2.GT.R0)THEN
        WRITE(11,*)'A minimum is found for'
        WRITE(*,*)'A minimum is found for'
    ELSE
        WRITE(11,*)'Linear function. Neither minimum nor maximum.'
        GOTO 99
    ENDIF
    XMIN=-A1/(R2*A2)
    FMIN=A0+A1*XMIN+A2*XMIN*XMIN
    WRITE(11,*)'X=',XMIN,'F_approx=',FMIN
    WRITE(*,*)'X=',XMIN,'F_approx=',FMIN
99 CONTINUE    
END SUBROUTINE