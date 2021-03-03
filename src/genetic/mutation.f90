! In this routine, some genes of offspring are mutated
! Igor Lopes, February 2015
SUBROUTINE MUTATION(POP,NPOP,PMUT,NSURV,NTBIT)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R0 /0.0D0/
    REAL(8) R2 /2.0D0/
    ! Argument
    INTEGER :: NPOP, NSURV, NTBIT
    INTEGER,DIMENSION(NPOP,NTBIT)   :: POP
    REAL(8) :: PMUT
    ! Locals
    INTEGER     :: IPOP, IBIT, ICOUNT
    REAL(8)     :: RANDN
    ! Initialize
    ICOUNT=0
    ! Algorithm
    DO IBIT=1,NTBIT
        DO IPOP=NSURV+1,NPOP
            CALL RANDOM_NUMBER(RANDN)
            IF(RANDN.LE.PMUT)THEN
                IF(POP(IPOP,IBIT).EQ.0)THEN
                    POP(IPOP,IBIT)=1
                ELSE
                    POP(IPOP,IBIT)=0
                ENDIF
                ICOUNT=ICOUNT+1
            ENDIF
        ENDDO
    ENDDO
    WRITE(11,*)'Nr. mutated bits=',ICOUNT
    WRITE(*,*)'Nr. mutated bits=',ICOUNT
END SUBROUTINE