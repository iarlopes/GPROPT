! Output routine for genetic algorithms
! Igor Lopes, February 2015
SUBROUTINE STATISTICS(IGEN,XVEC,VFIT,OBJFUN,NDIM,NPOP,NTOP)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R0 /0.0D0/
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) P10 /0.10D0/
    REAL(8) P90 /0.90D0/
    ! Argument
    INTEGER :: IGEN, NDIM, NPOP, NTOP
    REAL(8),DIMENSION(NPOP) :: VFIT, OBJFUN
    REAL(8),DIMENSION(NDIM,NPOP)    :: XVEC
    ! Locals
    INTEGER     :: IBIT, IDIM, IPOP, IPOS, ISUB, NPOS
    REAL(8)     :: AUX
    ! Initialize
100 FORMAT(I4,2X,I3,2X,<NDIM>(F12.6,3X),F12.6,3X,G12.6)
105 FORMAT('Objective Function max. value:',G12.6)
110 FORMAT('Objective Function min. value:',G12.6)    
115 FORMAT('Objective Function average:   ',G12.6)  
    ! Algorithm
    DO IPOP=1,NTOP
        WRITE(11,100)IGEN,IPOP,(XVEC(IDIM,IPOP),IDIM=1,NDIM),OBJFUN(IPOP),VFIT(IPOP)
        WRITE(*,100)IGEN,IPOP,(XVEC(IDIM,IPOP),IDIM=1,NDIM),OBJFUN(IPOP),VFIT(IPOP)
    ENDDO
    WRITE(11,105)MAXVAL(OBJFUN);WRITE(*,105)MAXVAL(OBJFUN)
    WRITE(11,110)MINVAL(OBJFUN);WRITE(*,110)MINVAL(OBJFUN)
    AUX=R0
    DO IPOP=1,NPOP
        AUX=AUX+OBJFUN(IPOP)
    ENDDO
    AUX=AUX/NPOP
    WRITE(11,115)AUX;WRITE(*,115)AUX
    write(10,*)minval(objfun),maxval(objfun),aux
99  CONTINUE            
END SUBROUTINE