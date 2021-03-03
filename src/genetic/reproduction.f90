! Routine for generating sons through crossover technique.
! Igor Lopes, February 2015
SUBROUTINE REPRODUCTION(POPO,POP,VFITO,VFIT,NSURV,NPOP,NTBIT)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R0 /0.0D0/
    REAL(8) R2 /2.0D0/
    ! Argument
    INTEGER :: NSURV, NPOP, NTBIT
    REAL(8),DIMENSION(NPOP) :: VFIT, VFITO
    INTEGER,DIMENSION(NPOP,NTBIT)   :: POPO, POP
    ! Locals
    INTEGER     :: IPOP, JPOP, IPOS, NPOS
    REAL(8),DIMENSION(NPOP) :: VAFIT
    REAL(8)     :: SUM, RANDN
    INTEGER,DIMENSION(2,NTBIT) :: PAR
    ! Initialize
    VAFIT=R0
    SUM=R0
    IPOS=NSURV
    ! Algorithm
    ! Parents
    VAFIT(1)=VFITO(1)
    SUM=VFITO(1)
    DO IPOP=2,NPOP
        VAFIT(IPOP)=VAFIT(IPOP-1)+VFITO(IPOP)
        SUM=SUM+VFITO(IPOP)
    ENDDO
    VAFIT=VAFIT/SUM
    DO WHILE(IPOS.LT.NPOP)
        ! 1st parent
        CALL RANDOM_NUMBER(RANDN)
        DO IPOP=1,NPOP
            IF(VAFIT(IPOP).GE.RANDN)THEN
                PAR(1,:)=POPO(IPOP,:)
                GOTO 10
            ENDIF
        ENDDO
        ! 2nd parent
10      CALL RANDOM_NUMBER(RANDN)
        DO JPOP=1,NPOP
            IF(VAFIT(JPOP).GE.RANDN)THEN
                IF(JPOP.EQ.IPOP)GOTO 10
                PAR(2,:)=POPO(JPOP,:)
                GOTO 15
            ENDIF
        ENDDO
        ! Position for crossover: NTBIT possible positions
15      CALL RANDOM_NUMBER(RANDN)
        RANDN=RANDN*(NTBIT-1)
        NPOS=CEILING(RANDN)
        ! Sons stored in new population vector
        ! 1st son
        IPOS=IPOS+1
        IF(IPOS.GT.NPOP)GOTO 99
        POP(IPOS,1:NPOS)=PAR(1,1:NPOS)
        POP(IPOS,NPOS+1:NTBIT)=PAR(2,NPOS+1:NTBIT)
        ! 2nd
        IPOS=IPOS+1
        IF(IPOS.GT.NPOP)GOTO 99
        POP(IPOS,1:NPOS)=PAR(2,1:NPOS)
        POP(IPOS,NPOS+1:NTBIT)=PAR(1,NPOS+1:NTBIT)
    ENDDO
99  CONTINUE    
END SUBROUTINE