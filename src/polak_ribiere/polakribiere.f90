! Routine for n-dimensional optimization through
! the Polak-Ribiere conjugate direction method.
! Igor Lopes, February 2015    
SUBROUTINE POLAKRIBIERE(NDIM)
    IMPLICIT NONE
    ! Parameters
    REAL(8) GOLD /1.61803398875D0/
    REAL(8) CONVTOL /1.0D-6/
    REAL(8) R0 /0.0D0/
    ! Arguments
    INTEGER :: NDIM
    ! Locals
    REAL(8)     :: FO,TOL,EVALFUNC,AMIN,FMIN,ANOR,BNOR,CNOR,BETA,SLOPE,RELDIF
    REAL(8),DIMENSION(NDIM) :: XO,XN,GRAD,SDIR,GRADO,XMIN,XMAX
    REAL(8),DIMENSION(2)    :: ALPHO
    REAL(8),DIMENSION(4)    :: ALPH,FN
    INTEGER     :: I,ITER,MITER,NITER
    CHARACTER   :: STRING*256
    ! Initialize
    STRING='polak_ribiere.res'
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
    CALL GRADIENT(GRADO,XO,NDIM)
    ANOR=DOT_PRODUCT(GRADO,GRADO)
    SDIR=-GRADO
    DO ITER=1,MITER
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
        XO=XN
        FO=FMIN
        CALL GRADIENT(GRAD,XO,NDIM)
        BNOR=DOT_PRODUCT(GRAD,GRAD)
        CNOR=DOT_PRODUCT(GRAD,GRADO)
        BETA=(BNOR-CNOR)/ANOR
        SDIR=BETA*SDIR-GRAD
        ANOR=BNOR
        SLOPE=DOT_PRODUCT(SDIR,GRAD)
        GRADO=GRAD
        IF(SLOPE.GE.R0)SDIR=-GRAD
    ENDDO
99  WRITE(11,5)
    WRITE(*,5)
    CLOSE(UNIT=11)
END SUBROUTINE