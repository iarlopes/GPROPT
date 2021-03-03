! Routine for cubic approximation
! Igor Lopes, February 2015
SUBROUTINE POLYCUBIC(NPOINT,X,F,XMIN,FMIN)
    IMPLICIT NONE
    REAL(8) R0 /0.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) R3 /3.0D0/
    REAL(8) SMALL /1.0D-8/
    ! Arguments
    INTEGER NPOINT
    REAL(8) :: XMIN, FMIN
    REAL(8),DIMENSION(4) :: X , F
    ! Locals
    REAL(8) :: A0,A1,A2,A3,Q1,Q2,Q3,Q4,Q5,Q6,XMAX,FMAX
    ! Begin algorithm
    IF(NPOINT.EQ.3)THEN
        Q1=(F(3)-F(1))/((X(3)-X(2))*(X(3)-X(1))**2)
        Q2=(F(2)-F(1))/((X(3)-X(2))*(X(2)-X(1))**2)
        Q3=F(4)/((X(3)-X(1))*(X(2)-X(1)))
        A3=Q1-Q2+Q3
        A2=((F(2)-F(1))/(X(2)-X(1))-F(4))/(X(2)-X(1))-A3*(R2*X(1)+X(2))
        A1=F(4)-R2*A2*X(1)-R3*A3*X(1)*X(1)
        A0=F(1)-A1*X(1)-A2*X(1)*X(1)-A3*X(1)*X(1)*X(1)
    ELSEIF(NPOINT.EQ.4)THEN
        Q1=X(3)**3*(X(2)-X(1))-X(2)**3*(X(3)-X(1))+X(1)**3*(X(3)-X(2))
        Q2=X(4)**3*(X(2)-X(1))-X(2)**3*(X(4)-X(1))+X(1)**3*(X(4)-X(2))
        Q3=(X(3)-X(2))*(X(2)-X(1))*(X(3)-X(1))
        Q4=(X(4)-X(2))*(X(2)-X(1))*(X(4)-X(1))
        Q5=F(3)*(X(2)-X(1))-F(2)*(X(3)-X(1))+F(1)*(X(3)-X(2))
        Q6=F(4)*(X(2)-X(1))-F(2)*(X(4)-X(1))+F(1)*(X(4)-X(2))
        A3=(Q3*Q6-Q4*Q5)/(Q2*Q3-Q1*Q4)
        A2=(Q5-A3*Q1)/Q3
        A1=(F(2)-F(1))/(X(2)-X(1))-A2*(X(1)+X(2))-A3*(X(2)**3-X(1)**3)/(X(2)-X(1))
        A0=F(1)-A1*X(1)-A2*X(1)*X(1)-A3*X(1)*X(1)*X(1)
    ELSE
        WRITE(*,*)'Wrong number of input points in POLYCUBIC'
        GOTO 99
    ENDIF
    WRITE(*,*)'Approximation coefficients:'
    WRITE(*,*)'a0=',A0
    WRITE(*,*)'a1=',A1
    WRITE(*,*)'a2=',A2
    WRITE(*,*)'a3=',A3
    WRITE(11,*)'Approximation coefficients:'
    WRITE(11,*)'a0=',A0
    WRITE(11,*)'a1=',A1
    WRITE(11,*)'a2=',A2
    WRITE(11,*)'a3=',A3
    ! Minimum or Maximum
    Q1=A2*A2-R3*A1*A3
    IF(Q1.LT.R0)THEN
        WRITE(*,*)'Unreal results'
        GOTO 99
    ELSEIF(Q1.LT.SMALL)THEN
        WRITE(*,*)'Neither minimum nor maximum is found'
        GOTO 99
    ENDIF
    XMIN=(-A2+DSQRT(Q1))/(R3*A3)
    FMIN=A0+A1*XMIN+A2*XMIN*XMIN+A3*XMIN**3
    XMAX=(-A2-DSQRT(Q1))/(R3*A3)
    FMAX=A0+A1*XMAX+A2*XMAX*XMAX+A3*XMAX**3
    WRITE(*,*)'A maximum is found for'
    WRITE(*,*)'X=',XMAX,'F_approx=',FMAX
    WRITE(*,*)'A minimum is found for'
    WRITE(*,*)'X=',XMIN,'F_approx=',FMIN
    WRITE(11,*)'A maximum is found for'
    WRITE(11,*)'X=',XMAX,'F_approx=',FMAX
    WRITE(11,*)'A minimum is found for'
    WRITE(11,*)'X=',XMIN,'F_approx=',FMIN
99 CONTINUE    
END SUBROUTINE