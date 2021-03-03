! Converts the binary string in BIN to its real value
! Igor Lopes, February 2015
DOUBLE PRECISION FUNCTION BIN2REAL(BIN,NBIT)
    IMPLICIT NONE
    REAL(8) R0 /0.0D0/
    REAL(8) R2 /2.0D0/
    ! Arguments
    INTEGER NBIT
    INTEGER,DIMENSION(NBIT) :: BIN
    ! Locals
    INTEGER :: IBIT
    ! Begin algorithm
    BIN2REAL=R0
    DO IBIT=1,NBIT
        BIN2REAL=BIN2REAL+BIN(IBIT)*R2**(NBIT-IBIT)
    ENDDO
    RETURN
END FUNCTION