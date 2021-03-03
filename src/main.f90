! GPROPT - Generalized PROgram for OPTimization
! 
! The purpose of this program is to be general enough so that several
! optimization algorithms can be implemented
!
! Igor Lopes, February 2015
PROGRAM GPROPT
    IMPLICIT NONE
! Declaration
    INTEGER :: NALGO, NDIM
! Begin algorithm
01  WRITE(*,*)'Introduce number of design variables (0 - QUIT):'
    READ(*,*) NDIM
    IF(NDIM.LT.0)THEN
        WRITE(*,*)'Number of design variables must be larger or equal to 1'
        GOTO 01
    ELSEIF(NDIM.EQ.1)THEN
        WRITE(*,*)'Uni-dimensional problem'
        GOTO 11
    ELSEIF(NDIM.EQ.0)THEN
        GOTO 99
    ENDIF
10  WRITE(*,*)'Select Optimization algorithm:'
    WRITE(*,*)'1 - Powell Method'
    WRITE(*,*)'2 - Steepest Descent'
    WRITE(*,*)'3 - Fletcher-Reeves Conjugate Direction'
    WRITE(*,*)'4 - Polak-Ribiere Conjugate Direction'
    WRITE(*,*)'5 - Quasi-Newton Methods (Variable Metrics)'
    WRITE(*,*)'6 - Newton Method'
    WRITE(*,*)'7 - Genetic Algorithm'
    WRITE(*,*)'0 - QUIT'
    WRITE(*,*)'(Introduce corresponding number)'
    READ(*,*) NALGO
    
    SELECT CASE(NALGO)
    CASE (1)
        WRITE(*,*)'Powell s Method'
        CALL POWELL(NDIM)
    CASE(2)
        WRITE(*,*)'Steepest Descent'
        CALL STEEPDESC(NDIM)
    CASE(3)
        WRITE(*,*)'Fletcher-Reeves Conjugate Direction'
        CALL FLETCHERREEVES(NDIM)
    CASE(4)
        WRITE(*,*)'Polak-Ribiere Conjugate Direction'
        CALL POLAKRIBIERE(NDIM)
    CASE(5)
        WRITE(*,*)'Quasi-Newton Methods'
        CALL QUASINEWTON(NDIM)
    CASE(6)
        WRITE(*,*)'Newton Method'
        CALL NEWTON(NDIM)
    CASE(7)
        WRITE(*,*)'Genetic Algorithm'
        CALL GENETIC(NDIM)  
    CASE DEFAULT
        WRITE(*,*)'Invalid option!'
        GOTO 10
    CASE(0)
        GOTO 99
        
    END SELECT
    GOTO 01
11  CONTINUE
    !1D
20  WRITE(*,*)'Select Optimization algorithm:'
    WRITE(*,*)'1 - Polynomial approximation'
    WRITE(*,*)'2 - Golden Section Method'
    WRITE(*,*)'3 - Genetic Algorithm'
    WRITE(*,*)'0 - QUIT'
    WRITE(*,*)'(Introduce corresponding number)'
    READ(*,*) NALGO
    
    SELECT CASE(NALGO)
    CASE (1)
        WRITE(*,*)'Polynomial approximation'
        CALL POLYINT
    CASE(2)
        WRITE(*,*)'Golden Section Method'
        CALL GOLDENINT
    CASE(3)
        WRITE(*,*)'Genetic Algorithm'
        CALL GENETIC(NDIM)  
    CASE DEFAULT
        WRITE(*,*)'Invalid option!'
        GOTO 20
    CASE(0)
        GOTO 99
    END SELECT
    GOTO 01

99  CONTINUE 
    WRITE(*,*)'CLOSING...'
END PROGRAM
