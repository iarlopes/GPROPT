! Routine that checks termination criteria
! Igor Lopes, February 2015
SUBROUTINE TERMINATION(POPO,POP,NPOP,NTBIT,NTOPT,NGENT,ICOUNT,THEEND)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) P10 /0.10D0/
    REAL(8) P90 /0.90D0/
    ! Argument
    INTEGER :: NPOP, NTBIT, NTOPT, NGENT, ICOUNT
    INTEGER,DIMENSION(NPOP,NTBIT)   :: POP,POPO
    LOGICAL :: THEEND
    ! Locals
    INTEGER     :: IBIT, IPOP
    REAL(8)     :: AUX, ISUM
    ! Initialize
100 FORMAT('TOP ',I3,' has not changed in the last ',I3,' generations!')
105 FORMAT('Reduced diversity in population!')    
    ! Algorithm
    ! Termination criteria
    ! N generations without changing top
    DO IBIT=1,NTBIT
        DO IPOP=1,NTOPT
            IF(POP(IPOP,IBIT).NE.POPO(IPOP,IBIT))THEN
                ICOUNT=0
                GOTO 10
            ENDIF
        ENDDO
    ENDDO
    ICOUNT=ICOUNT+1
    IF(ICOUNT.EQ.NGENT)THEN
        WRITE(11,100)NTOPT,NGENT
        WRITE(*,100)NTOPT,NGENT
        THEEND=.TRUE.
        GOTO 99
    ENDIF
10  CONTINUE    
    ! Check diversity
    DO IBIT=1,NTBIT
        ISUM=0
        DO IPOP=1,NPOP
            ISUM=ISUM+POP(IPOP,IBIT)
        ENDDO
        AUX=ISUM/NPOP
        IF(AUX.GT.P10.AND.AUX.LT.P90)GOTO 99
    ENDDO
    WRITE(11,105)
    WRITE(*,105)
    THEEND=.TRUE.
    !
99  CONTINUE            
END SUBROUTINE