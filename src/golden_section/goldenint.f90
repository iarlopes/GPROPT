! Interface routine for uni-dimensional optimization through
! the golden section method.
! Igor Lopes, February 2015    
SUBROUTINE GOLDENINT
    IMPLICIT NONE
    ! Parameters
    REAL(8) R0 /0.0D0/
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) R5 /5.0D0/
    !Locals
    REAL(8),DIMENSION(2) ::  XO
    REAL(8),DIMENSION(4) ::  XF, F
    REAL(8) ::  TOL, GOLD, XI, XMIN, XMAX, FMIN, EVALFUNC1,FI
    INTEGER ::  I , NITER
    CHARACTER STRING*256,ANS*1
    
    ! Innititalization
    XO=R0; XF=R0
    
    GOLD=(R1+DSQRT(R5))/R2
    STRING='goldensection.res'
    !
    CALL RESULTFILE(STRING)
    ! Begin algorithm
    WRITE(*,*)'Automatic bounds search?(Y/N)'
    READ(*,*)ANS
    IF(ANS.EQ.'Y'.OR.ANS.EQ.'y')THEN
        WRITE(*,*)'Give a lower limit for X...'
        READ(*,*)XMIN
        WRITE(*,*)'and an upper limit...'
        READ(*,*)XMAX
        WRITE(*,*)'and an initial point for the search.'
        READ(*,*)XI
        FI=EVALFUNC1(XI)
        CALL FINDBOUNDS1D(XI,FI,XO(1),XO(2),XMAX,XMIN)
        WRITE(*,*)'Automatic bounds'
        WRITE(11,*)'Automatic bounds'
        WRITE(*,*)'X_l=',XO(1),'X_u=',XO(2)
        WRITE(11,*)'X_l=',XO(1),'X_u=',XO(2)
        GOTO 10
    ENDIF
    
    WRITE(*,*)'Insert lower and upper bounds:'
    WRITE(*,*)'X_lower='
    READ(*,*)XO(1)
    WRITE(*,*)'X_upper='
    READ(*,*)XO(2)
10  WRITE(*,*)'Define relative tolerance:'
    READ(*,*)TOL
    WRITE(11,*)'Tolerance:',TOL
    ! Number of iterations
    NITER=INT(LOG(TOL)/LOG(GOLD-R1))+1
    ! Determine minimum
    CALL GOLDENSECTION(XO,XF,F,NITER)
    XMIN=XF(1)
    FMIN=F(1)
    DO I=2,4
        IF(F(I).LT.FMIN)THEN
            FMIN=F(I)
            XMIN=XF(I)
        ENDIF
    ENDDO
    WRITE(*,*) 'Minimum: X_min=',XMIN,' with F_min=',FMIN
    WRITE(11,*) 'Minimum: X_min=',XMIN,' with F_min=',FMIN
    CLOSE(UNIT=11)
END SUBROUTINE