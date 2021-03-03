! Routine for determining step length bounds. It is used prior to the
! the golden section method for unidirectional search in n-dimensional
! problem.
! Igor Lopes, February 2015 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ARGUMENTS
!! XO,FO        >   Initial point and corresponding function
!! ALPHL,ALPHU  <   Lower and upper limits
!! XMAX,XMIN    >   Side constraints: limits for domain
!! SDIR         >   Search direction vector
!! NDIM         >   Dimension of the problem (nr design variables)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDBOUNDS(XO,FO,ALPHL,ALPHU,XMAX,XMIN,SDIR,NDIM)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R0 /0.0D0/
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) R5 /5.0D0/
    REAL(8) DELTA /1.0D-2/
    REAL(8) INF /1.0D-8/
    INTEGER MITER /50/
    ! Arguments
    INTEGER NDIM
    REAL(8),DIMENSION(NDIM) :: XO, XMAX, XMIN, SDIR
    REAL(8)     :: ALPHL, ALPHU, FO
    ! Locals
    REAL(8),DIMENSION(NDIM) :: XN
    REAL(8)     :: GOLD, ALPH1, AMAX, AMIN, F1, FU, FL, EVALFUNC, AUX
    INTEGER     :: I, NITER
    ! Initialize
    GOLD=(R1+DSQRT(R5))/R2
    AMAX=1.0E9
    AMIN=-1.0E9
    DO I=1,NDIM
        IF(DABS(SDIR(I)).LE.INF)GOTO 5
        IF(SDIR(I).GT.INF)THEN
            AUX=(XMAX(I)-XO(I))/SDIR(I)
            IF(AUX.LT.AMAX) AMAX=AUX
            AUX=(XMIN(I)-XO(I))/SDIR(I)
            IF(AUX.GT.AMIN) AMIN=AUX
        ELSE
            AUX=(XMIN(I)-XO(I))/SDIR(I)
            IF(AUX.LT.AMAX) AMAX=AUX
            AUX=(XMAX(I)-XO(I))/SDIR(I)
            IF(AUX.GT.AMIN) AMIN=AUX
        ENDIF
5       CONTINUE
    ENDDO
    IF(AMAX.LT.AMIN)THEN
        PAUSE
    ENDIF
    ! Begin algorithm
    ALPHU=DELTA*AMAX
    XN=XO+ALPHU*SDIR
    FU=EVALFUNC(XN,NDIM)
    IF(FU.GE.FO)THEN
        ! Upper limit is found... but lower limit is not
        GOTO 10
    ENDIF
    ALPHL=R0
    FL=FO
    DO
        ALPH1=ALPHU
        F1=FU
        ALPHU=ALPH1+(ALPH1-ALPHL)*GOLD
        IF(ALPHU.GT.AMAX) GOTO 15
        DO I=1,NDIM
            XN(I)=XO(I)+ALPHU*SDIR(I)
        ENDDO
        FU=EVALFUNC(XN,NDIM)
        IF(FU.GE.F1)GOTO 99
        ALPHL=ALPH1
        FL=F1
    ENDDO
    ! Find lower limit
10  ALPHL=DELTA*AMIN
    XN=XO+ALPHL*SDIR
    FL=EVALFUNC(XN,NDIM)
    IF(FL.GE.FO)GOTO 99 !XL is found
    DO
        ALPH1=ALPHL
        F1=FL
        ALPHL=ALPHL-(ALPHU-ALPH1)*GOLD
        IF(ALPHL.LT.AMIN)GOTO 20
        DO I=1,NDIM
            XN(I)=XO(I)+ALPHL*SDIR(I)
        ENDDO
        FL=EVALFUNC(XN,NDIM)
        IF(FL.GE.F1)GOTO 99
        ALPHU=ALPH1
        FU=F1
    ENDDO
    
15  ALPHU=AMAX
    DO I=1,NDIM
        XN(I)=XO(I)+ALPHU*SDIR(I)
    ENDDO
    FU=EVALFUNC(XN,NDIM)
    GOTO 99
    
20  ALPHL=AMIN
    DO I=1,NDIM
        XN(I)=XO(I)+ALPHL*SDIR(I)
    ENDDO
    FL=EVALFUNC(XN,NDIM)
           
99  CONTINUE
END SUBROUTINE