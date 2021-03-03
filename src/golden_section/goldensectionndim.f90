! Routine for n-dimensional optimization through
! the golden section method.
! Igor Lopes, February 2015    
SUBROUTINE GOLDENSECTIONNDIM(XO,AO,AN,FN,NITER,SDIR,NDIM)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R0 /0.0D0/
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) R5 /5.0D0/
    ! Arguments
    INTEGER ::  NITER,NDIM
    REAL(8),DIMENSION(NDIM) ::  XO,SDIR
    REAL(8),DIMENSION(4) ::  AN, FN
    REAL(8),DIMENSION(2) ::  AO
    ! Locals
    INTEGER     :: ITER
    REAL(8)     :: TOL, TAU, EVALFUNC, GOLD
    REAL(8),DIMENSION(NDIM) ::  XN
    ! Initialize
    FN=R0
    GOLD=(R1+DSQRT(R5))/R2
    TAU=R1/GOLD
    ! Begin algorithm
    AN(1)=AO(1)
    AN(4)=AO(2)
    XN=XO+AN(1)*SDIR
    FN(1)=EVALFUNC(XN,NDIM)
    XN=XO+AN(4)*SDIR
    FN(4)=EVALFUNC(XN,NDIM)
    AN(2)=AN(4)-TAU*(AN(4)-AN(1))
    AN(3)=AN(1)+TAU*(AN(4)-AN(1))
    XN=XO+AN(2)*SDIR
    FN(2)=EVALFUNC(XN,NDIM)
    XN=XO+AN(3)*SDIR
    FN(3)=EVALFUNC(XN,NDIM)
    DO ITER=1,NITER
        IF(FN(2).GT.FN(3))THEN
            AN(1)=AN(2)
            FN(1)=FN(2)
            AN(2)=AN(3)
            FN(2)=FN(3)
            AN(3)=AN(1)+TAU*(AN(4)-AN(1))
            XN=XO+AN(3)*SDIR
            FN(3)=EVALFUNC(XN,NDIM)
        ELSEIF(FN(2).LT.FN(3))THEN
            AN(4)=AN(3)
            FN(4)=FN(3)
            AN(3)=AN(2)
            FN(3)=FN(2)
            AN(2)=AN(4)-TAU*(AN(4)-AN(1))
            XN=XO+AN(2)*SDIR
            FN(2)=EVALFUNC(XN,NDIM)
        ELSE
            AN(1)=AN(2)
            FN(1)=FN(2)
            AN(4)=AN(3)
            FN(4)=FN(3)
            AN(2)=AN(4)-TAU*(AN(4)-AN(1))
            XN=XO+AN(2)*SDIR
            FN(2)=EVALFUNC(XN,NDIM)
            AN(3)=AN(1)+TAU*(AN(4)-AN(1))
            XN=XO+AN(3)*SDIR
            FN(3)=EVALFUNC(XN,NDIM)
        ENDIF
    ENDDO
END SUBROUTINE