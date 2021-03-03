! Evaluate function of one variable X
DOUBLE PRECISION FUNCTION EVALFUNC1(X)
    IMPLICIT NONE
    !
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) R3 /3.0D0/
    !
    REAL(8) :: X
    !
    EVALFUNC1=R1-R3*X+EXP(R2*X)
    RETURN
END FUNCTION