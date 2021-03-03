! Routine for n-dimensional optimization through
! the Powell's method.
! Igor Lopes, February 2015    
SUBROUTINE POWELL(NDIM)
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
    REAL(8)                 :: GOLD, FO, TOL, EVALFUNC, AMIN, FMIN, RELDIF
    REAL(8),DIMENSION(NDIM) :: XO, XN, GRAD, SDIR, SDCJ, XMIN, XMAX
    REAL(8),DIMENSION(NDIM,NDIM) :: HMAT
    REAL(8),DIMENSION(2)    :: ALPHO
    REAL(8),DIMENSION(4)    :: ALPH, FN
    INTEGER                 :: I, IDIM, ITER, MITER, NITER
    CHARACTER               :: STRING*256
    ! Initialize
    STRING='powell.res'
    CALL RESULTFILE(STRING)
    HMAT=R0
    DO IDIM=1,NDIM
        HMAT(IDIM,IDIM)=R1
    ENDDO
    GOLD=(R1+DSQRT(R5))/R2
    ! Formats
105   FORMAT('-------------------------------------------------------------------------')
110  FORMAT('ITER.',6X,<NDIM>('X_'I3,10X),'F',14X,'alpha')
115  FORMAT(I3,4X,<NDIM>(G12.6,3X),G12.6,3X,G12.6)   
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
    WRITE(11,*)'Maximum number of iterations=',MITER
    WRITE(11,*)'Tolerance for step search=',TOL
    NITER=INT(LOG(TOL)/LOG(GOLD-1))+1
    !
    FO=EVALFUNC(XO,NDIM)
    WRITE(11,105)
    WRITE(*,105)
    WRITE(11,110)(I,I=1,NDIM)
    WRITE(*,110)(I,I=1,NDIM)
    WRITE(11,105)
    WRITE(*,105)
    WRITE(11,115)0,XO,FO
    WRITE(*,115)0,XO,FO
    WRITE(11,105)
    WRITE(*,105)
    DO ITER=1,MITER
        SDCJ=R0
        DO IDIM=1,NDIM
            SDIR=HMAT(:,IDIM)
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
            XO=XN
            FO=FMIN
            HMAT(:,IDIM)=AMIN*SDIR
            SDCJ=SDCJ+AMIN*SDIR
            WRITE(11,115)ITER,(XN(I),I=1,NDIM),FMIN,AMIN
            WRITE(*,115)ITER,(XN(I),I=1,NDIM),FMIN,AMIN
        ENDDO
        CALL FINDBOUNDS(XO,FO,ALPHO(1),ALPHO(2),XMAX,XMIN,SDCJ,NDIM)
        CALL GOLDENSECTIONNDIM(XO,ALPHO,ALPH,FN,NITER,SDCJ,NDIM)
        ! Determine minimum
        AMIN=ALPH(1)
        FMIN=FN(1)
        DO I=2,4
            IF(FN(I).LT.FMIN)THEN
                FMIN=FN(I)
                AMIN=ALPH(I)
            ENDIF
        ENDDO
        XN=XO+AMIN*SDCJ
        WRITE(11,115)ITER,(XN(I),I=1,NDIM),FMIN,AMIN
        WRITE(*,115)ITER,(XN(I),I=1,NDIM),FMIN,AMIN
        WRITE(11,105)
        WRITE(*,105)  
        ! Check convergence
        RELDIF=DABS(FMIN-FO)
        IF(DABS(FO).GT.CONVTOL)RELDIF=RELDIF/DABS(FO)
        IF(RELDIF.LT.CONVTOL)THEN
            !FINISH PROCESS
            GOTO 99
        ENDIF
        XO=XN
        FO=FMIN
        DO IDIM=1,NDIM-1
            HMAT(:,IDIM)=HMAT(:,IDIM+1)
        ENDDO
        HMAT(:,NDIM)=SDCJ
    ENDDO
99  WRITE(11,105)
    WRITE(*,105)
    CLOSE(UNIT=11)
END SUBROUTINE