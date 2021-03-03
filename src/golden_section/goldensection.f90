! Routine for uni-dimensional optimization through
! the golden section method.
! Igor Lopes, February 2015    
SUBROUTINE GOLDENSECTION(XO,XF,F,NITER)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R0 /0.0D0/
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) R5 /5.0D0/
    ! Arguments
    REAL(8),DIMENSION(2) ::  XO
    REAL(8),DIMENSION(4) ::  XF, F
    INTEGER ::  NITER
    ! Locals
    REAL(8)     :: TOL, TAU, EVALFUNC1, GOLD
    INTEGER     :: ITER
    
    ! Initialize
    F=R0
    GOLD=(R1+DSQRT(R5))/R2
    ! Formats  
105 FORMAT('---------------------------------------------------------------------------------------------')
110 FORMAT('ITER.',6X,'X_l',13X,'X_1',13X,'X_2',13X,'X_u',/11X,'F_l',13X,'F_1',13X,'F_2',13X,'F_u',13X,'tol.')
115 FORMAT(I3,4X,F12.6,4X,F12.6,4X,F12.6,4X,F12.6,/11X,G12.6,4X,G12.6,4X,G12.6,4X,G12.6,4X,F12.6)
    WRITE(11,105)
    WRITE(11,110)
    WRITE(11,105)
    TAU=R1/GOLD
    ! Begin algorithm
    XF(1)=XO(1)
    XF(4)=XO(2)
    F(1)=EVALFUNC1(XF(1))
    F(4)=EVALFUNC1(XF(4))
    XF(2)=XF(4)-TAU*(XF(4)-XF(1))
    XF(3)=XF(1)+TAU*(XF(4)-XF(1))
    F(2)=EVALFUNC1(XF(2))
    F(3)=EVALFUNC1(XF(3))
    TOL=(XF(4)-XF(1))/(XO(2)-XO(1))
    WRITE(11,115)0,XF(1),XF(2),XF(3),XF(4),F(1),F(2),F(3),F(4),TOL
    WRITE(11,105)
    DO ITER=1,NITER
        IF(F(2).GT.F(3))THEN
            XF(1)=XF(2)
            F(1)=F(2)
            XF(2)=XF(3)
            F(2)=F(3)
            XF(3)=XF(1)+TAU*(XF(4)-XF(1))
            F(3)=EVALFUNC1(XF(3))
        ELSEIF(F(2).LT.F(3))THEN
            XF(4)=XF(3)
            F(4)=F(3)
            XF(3)=XF(2)
            F(3)=F(2)
            XF(2)=XF(4)-TAU*(XF(4)-XF(1))
            F(2)=EVALFUNC1(XF(2))
        ELSE
            XF(1)=XF(2)
            F(1)=F(2)
            XF(4)=XF(3)
            F(4)=F(3)
            XF(2)=XF(4)-TAU*(XF(4)-XF(1))
            F(2)=EVALFUNC1(XF(2))
            XF(3)=XF(1)+TAU*(XF(4)-XF(1))
            F(3)=EVALFUNC1(XF(3))
        ENDIF
        TOL=(XF(4)-XF(1))/(XO(2)-XO(1))
        WRITE(11,115)ITER,XF(1),XF(2),XF(3),XF(4),F(1),F(2),F(3),F(4),TOL
        WRITE(11,105)
    ENDDO
END SUBROUTINE