! Routine that computes fitness of each individual,
! and sorts them by decreasing fitness.
! Igor Lopes, February 2015
SUBROUTINE FITNESS(XMIN,XMAX,XVEC,POP,VFIT,OBJF,NBITS,NDIM,NPOP,NTBIT)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    ! Argument
    INTEGER :: NDIM, NPOP, NTBIT
    REAL(8) :: FMAX, FMIN
    INTEGER,DIMENSION(NDIM) :: NBITS
    REAL(8),DIMENSION(NPOP) :: VFIT, OBJF
    REAL(8),DIMENSION(NDIM) :: XMIN, XMAX
    REAL(8),DIMENSION(NDIM,NPOP)    :: XVEC
    INTEGER,DIMENSION(NPOP,NTBIT)   :: POP
    ! Locals
    INTEGER     :: IBIT, IDIM, IPOP, IPOS, ISUB, NPOS
    REAL(8)     :: AUX, BIN2REAL, EVALFUNC
    REAL(8),DIMENSION(NPOP) :: VFITO, OBJFO
    INTEGER,DIMENSION(NPOP,NTBIT)   :: POPO
    REAL(8),DIMENSION(NDIM,NPOP)    :: XVECO
    ! Algorithm
    ! Decode from binary to real
    DO IPOP=1,NPOP
        IPOS=0
        DO IDIM=1,NDIM
            AUX=BIN2REAL(POP(IPOP,IPOS+1:IPOS+NBITS(IDIM)),NBITS(IDIM))
            XVEC(IDIM,IPOP)=XMIN(IDIM)+AUX*(XMAX(IDIM)-XMIN(IDIM))/(R2**NBITS(IDIM)-R1)
            IPOS=IPOS+NBITS(IDIM)
        ENDDO
        OBJF(IPOP)=EVALFUNC(XVEC(:,IPOP),NDIM)
    ENDDO
    ! Minimization of function -> higher fitness for lower function value
    FMAX=MAXVAL(OBJF)
    FMIN=MINVAL(OBJF)
    DO IPOP=1,NPOP
        VFIT(IPOP)=FMAX-OBJF(IPOP)+0.01*(FMAX-FMIN)
        !artifitial increase of fitness in 1% of the best fitness
    ENDDO
    ! Sort by decreasing fitness
    POPO=POP
    VFITO=VFIT
    XVECO=XVEC
    OBJFO=OBJF
    POP=0
    VFIT=MINVAL(VFIT)
    POP(1,:)=POPO(1,:)
    VFIT(1)=VFITO(1)
    XVEC(:,1)=XVECO(:,1)
    OBJF(1)=OBJFO(1)
    DO IPOP=2,NPOP
        DO ISUB=1,IPOP-1
            IF(VFITO(IPOP).GE.VFIT(ISUB))THEN
                DO IPOS=IPOP,ISUB+1,-1
                    VFIT(IPOS)=VFIT(IPOS-1)
                    POP(IPOS,:)=POP(IPOS-1,:)
                    OBJF(IPOS)=OBJF(IPOS-1)
                    XVEC(:,IPOS)=XVEC(:,IPOS-1)
                ENDDO
                VFIT(ISUB)=VFITO(IPOP)
                POP(ISUB,:)=POPO(IPOP,:)
                OBJF(ISUB)=OBJFO(IPOP)
                XVEC(:,ISUB)=XVECO(:,IPOP)
                GOTO 10
            ENDIF
        ENDDO
        VFIT(IPOP)=VFITO(IPOP)
        POP(IPOP,:)=POPO(IPOP,:)
        OBJF(IPOP)=OBJFO(IPOP)
        XVEC(:,IPOP)=XVECO(:,IPOP)
10      CONTINUE        
    ENDDO
END SUBROUTINE