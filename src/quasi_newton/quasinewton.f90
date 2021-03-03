! Routine for n-dimensional optimization through
! quasi-Newton methods: DFP or BFGS.
! Igor Lopes, February 2015    
SUBROUTINE QUASINEWTON(NDIM)
    IMPLICIT NONE
    ! Parameters
    REAL(8) CONVTOL /1.0D-6/
    REAL(8) R0 /0.0D0/
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) R5 /5.0D0/
    ! Arguments
    INTEGER :: NDIM
    ! Locals
    REAL(8)         :: GOLD, FO, TOL, EVALFUNC, AMIN, FMIN, ANOR, &
                            BNOR, BETA, THETA, SIGMA, TAU, RELDIF
    REAL(8),DIMENSION(NDIM)      :: XO, XN, GRAD, GRADO, SDIR, XMAX, XMIN
    REAL(8),DIMENSION(NDIM,1)    :: DX, DGRAD, AVEC
    REAL(8),DIMENSION(1,NDIM)    :: DXT, AVECT
    REAL(8),DIMENSION(NDIM,NDIM) :: HMAT, DMAT, AMAT
    REAL(8),DIMENSION(2)         :: ALPHO
    REAL(8),DIMENSION(4)         :: ALPH, FN
    INTEGER         :: I, ITER, MITER, NITER
    CHARACTER       :: STRING*256
    ! Initialize
    GOLD=(R1+DSQRT(R5))/R2
    STRING='quasi_newton.res'
    CALL RESULTFILE(STRING)
    HMAT=R0
    DO I=1,NDIM
        HMAT(I,I)=R1
    ENDDO
    ! Formats
5   FORMAT('-------------------------------------------------------------------------')
10  FORMAT('ITER.',6X,<NDIM>('X_'I3,10X),'F',14X,'alpha')
15  FORMAT(I3,4X,<NDIM>(G12.6,3X),G12.6,3X,G12.6)   
    ! Begin algorithm
20  WRITE(*,*)'Define THETA value:'
    WRITE(*,*)'0 - DFP (minimum)'
    WRITE(*,*)'1 - BFGS (maximum)'
    READ(*,*)THETA
    IF(THETA.LT.R0.OR.THETA.GT.R1)THEN
        WRITE(*,*)'THETA out of range.'
        GOTO 20
    ENDIF
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
    WRITE(11,*)'THETA=',THETA
    WRITE(11,*)'Maximum number of iterations=',MITER
    WRITE(11,*)'Tolerance for step search=',TOL
    NITER=INT(LOG(TOL)/LOG(GOLD-R1))+R1
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
    DO ITER=1,MITER
        CALL MATPRD(HMAT,GRAD,SDIR,NDIM,NDIM,1)
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
        DO I=1,NDIM
            DX(I,1)=XN(I)-XO(I)
            DGRAD(I,1)=GRAD(I)-GRADO(I)
        ENDDO
        SIGMA=DOT_PRODUCT(DX(:,1),DGRAD(:,1))
        CALL TRANSPOSE(DX,DXT,NDIM,1)
        CALL MATPRD(HMAT,DGRAD,AVEC,NDIM,NDIM,1)
        TAU=DOT_PRODUCT(DGRAD(:,1),AVEC(:,1))
        CALL MATPRD(DX,DXT,AMAT,NDIM,1,NDIM)
        DMAT=AMAT*(SIGMA+TAU*THETA)/(SIGMA*SIGMA)
        CALL TRANSPOSE(AVEC,AVECT,NDIM,1)
        CALL MATPRD(AVEC,AVECT,AMAT,NDIM,1,NDIM)
        DMAT=DMAT+AMAT*(THETA-R1)/TAU
        CALL MATPRD(AVEC,DXT,AMAT,NDIM,1,NDIM)
        DMAT=DMAT-AMAT*THETA/SIGMA
        CALL MATPRD(DX,AVECT,AMAT,NDIM,1,NDIM)
        DMAT=DMAT-AMAT*THETA/SIGMA
        HMAT=HMAT+DMAT
        XO=XN
        FO=FMIN
    ENDDO
99  WRITE(11,5)
    WRITE(*,5)
    CLOSE(UNIT=11)
END SUBROUTINE