! Routine for n-dimensional optimization through
! the Newton method.
! Igor Lopes, February 2015    
SUBROUTINE NEWTON(NDIM)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) R5 /5.0D0/
    REAL(8) CONVTOL /1.0D-6/
    REAL(8) INF /1.0D-8/
    ! Arguments
    INTEGER :: NDIM
    ! Locals
    INTEGER                     :: I, ITER, MITER, NITER
    REAL(8)                     :: FO, TOL, EVALFUNC, AMIN, FMIN, GOLD, RELDIF
    REAL(8),DIMENSION(NDIM)     :: XO, XN, GRAD, GRADO, SDIR, XMAX, XMIN
    REAL(8),DIMENSION(NDIM,NDIM):: HMAT, HINV
    REAL(8),DIMENSION(2)        :: ALPHO
    REAL(8),DIMENSION(4)        :: ALPH, FN
    CHARACTER   :: STRING*256
    LOGICAL     :: ERROR
    ! Initialize
    GOLD=(R1+DSQRT(R5))/R2
    STRING='newton.res'
    CALL RESULTFILE(STRING)
    ! Formats
5   FORMAT('-------------------------------------------------------------------------')
10  FORMAT('ITER.',6X,<NDIM>('X_'I3,10X),'F',14X,'alpha')
15  FORMAT(I3,4X,<NDIM>(G12.6,3X),G12.6,3X,G12.6)   
    ! Begin algorithm
    WRITE(*,*)'Define limits for each variable:'
    WRITE(11,*)'Variables limits:'
    WRITE(11,*)'i     Xmin     Xmax'
    DO I=1,NDIM
        WRITE(*,*)'Xmin_',I,' Xmax_',I
        READ(*,*)XMIN(I),XMAX(I)
        WRITE(11,*)I,XMIN(I),XMAX(I)
    ENDDO
    WRITE(*,*)'Give an initial guess for design variables X'
    READ(*,*)(XO(I),I=1,NDIM)
    WRITE(*,*)'Define maximum number of iterations'
    READ(*,*)MITER
    WRITE(*,*)'Define tolerance for step search'
    READ(*,*)TOL
    NITER=INT(LOG(TOL)/LOG(GOLD-1))+1
    !
    FO=EVALFUNC(XO,NDIM)
    WRITE(11,5)
    WRITE(*,5)
    WRITE(11,10)(I,I=1,NDIM)
    WRITE(*,10)(I,I=1,NDIM)
    WRITE(11,5)
    WRITE(*,5)
    WRITE(11,15)0,XO,FO
    WRITE(*,15)0,XO,FO
    !
    CALL GRADIENT(GRAD,XO,NDIM)
    CALL HESSIAN(HMAT,XO,NDIM)
    DO ITER=1,MITER
        CALL RMINVE(HMAT,HINV,NDIM,ERROR)
        CALL MATPRD(HINV,GRAD,SDIR,NDIM,NDIM,1)
        SDIR=-SDIR
        CALL FINDBOUNDS(XO,FO,ALPHO(1),ALPHO(2),XMAX,XMIN,SDIR,NDIM)
        CALL GOLDENSECTIONNDIM(XO,ALPHO,ALPH,FN,NITER,SDIR,NDIM)
        ! Determine minimum
        AMIN=ALPH(1)
        FMIN=FN(1)
        DO I=2,4
            IF(FN(I).LT.FMIN)THEN
                FMIN=FN(I)
                AMIN=ALPH(I)
            ENDIF
        ENDDO
        XN=XO+AMIN*SDIR
        WRITE(11,15)ITER,(XN(I),I=1,NDIM),FMIN,AMIN
        WRITE(*,15)ITER,(XN(I),I=1,NDIM),FMIN,AMIN
        ! Check convergence
        RELDIF=DABS(FMIN-FO)
        IF(DABS(FO).GT.CONVTOL)RELDIF=RELDIF/DABS(FO)
        IF(RELDIF.LT.CONVTOL)THEN
            !FINISH PROCESS
            GOTO 99
        ENDIF
        ! Update HMAT
        GRADO=GRAD
        CALL GRADIENT(GRAD,XN,NDIM)
        CALL HESSIAN(HMAT,XN,NDIM)
        XO=XN
        FO=FMIN
    ENDDO
99  WRITE(11,5)
    WRITE(*,5)
    CLOSE(UNIT=11)
END SUBROUTINE