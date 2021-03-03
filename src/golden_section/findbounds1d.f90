! Routine for finding initial bounds for the golden 
! section method.
! Igor Lopes, February 2015 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ARGUMENTS
!! XO,FO        >   Initial step and corresponding function
!! XL,XU        <   Lower and upper limits
!! XMAX,XMIN    >   Side constraints: limits for domain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDBOUNDS1D(XO,FO,XL,XU,XMAX,XMIN)
    IMPLICIT NONE
    ! Parameters
    REAL(8) GOLD /1.61803398875/
    REAL(8) R1 /1.0D0/
    INTEGER MITER /50/
    ! Arguments
    REAL(8)     :: XO,XL,XU,XMAX,XMIN,FO
    ! Locals
    REAL(8)     :: XN,X1,FU,FL,F1,EVALFUNC1
    INTEGER     :: I, NITER
    ! Begin algorithm
    XU=XO+R1
    FU=EVALFUNC1(XU)
    IF(FU.GE.FO)THEN
        ! Upper limit is found... but lower limit is not
        GOTO 10
    ENDIF
    XL=XO
    FL=FO
    DO
        X1=XU
        F1=FU
        XU=X1+(X1-XL)*GOLD
        IF(XU.GT.XMAX)THEN
            XU=XMAX
            GOTO 99
        ENDIF
        FU=EVALFUNC1(XU)
        IF(FU.GE.F1)GOTO 99
        XL=X1
        FL=F1
    ENDDO
    ! Find lower limit
10  XL=XO-R1
    FL=EVALFUNC1(XL)
    IF(FL.GE.FO)GOTO 99 !XL is found
    DO
        X1=XL
        F1=FL
        XL=X1-(XU-X1)*GOLD
        IF(XL.LT.XMIN)THEN
            XL=XMIN
            GOTO 99
        ENDIF
        FL=EVALFUNC1(XL)
        IF(FL.GE.F1)GOTO 99
        XU=X1
        FU=F1
    ENDDO           
99  CONTINUE
END SUBROUTINE