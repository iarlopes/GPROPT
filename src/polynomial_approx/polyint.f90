! Interface routine for uni-dimensional optimization through
! polynomial approximation.
! Igor Lopes, February 2015    
SUBROUTINE POLYINT
    IMPLICIT NONE
    
    REAL(8) R0 /0.0D0/
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) R5 /5.0D0/
    
    REAL(8),DIMENSION(4) ::  X,F
    REAL(8) ::  XMIN, FMIN
    INTEGER ::  NPOLY
    CHARACTER STRING*256
    
    ! Innititalization
    X=R0; F=R0    
    
101 FORMAT('Linear approximation: F=a0+a1*X')
102 FORMAT('Quadratic approximation: F=a0+a1*X+a2*X^2')
103 FORMAT('Cubic approximation: F=a0+a1*X+a2*X^2+a3*X^3')
110 FORMAT("Insert derivative F'_1")
111 FORMAT('Insert point X_1,F_1')
112 FORMAT('Insert point X_2,F_2')
113 FORMAT('Insert point X_3,F_3')
114 FORMAT('Insert point X_4,F_4')
120 FORMAT("F'_1=",F10.5)
121 FORMAT('X_1=',F10.5,2X,'F_1=',F10.5)
122 FORMAT('X_2=',F10.5,2X,'F_2=',F10.5)
123 FORMAT('X_3=',F10.5,2X,'F_3=',F10.5)
124 FORMAT('X_4=',F10.5,2X,'F_4=',F10.5)
    ! Begin algorithm
10  WRITE(*,*)'Choose kind of polynomial approximation:'
    WRITE(*,*)'1 - Linear 1-point'
    WRITE(*,*)'2 - Linear 2-point'
    WRITE(*,*)'3 - Quadratic 2-point'
    WRITE(*,*)'4 - Quadratic 3-point'
    WRITE(*,*)'5 - Cubic 3-point'
    WRITE(*,*)'6 - Cubic 4-point'
    WRITE(*,*)'0 - QUIT'
    READ(*,*) NPOLY
    
    SELECT CASE(NPOLY)
    CASE(0)
        WRITE(*,*)'CLOSING...'
        STOP
    CASE(1)
        STRING='linear1.res'
        CALL RESULTFILE(STRING)
        WRITE(11,101)
        WRITE(*,101)
        WRITE(*,111);READ(*,*)X(1),F(1)
        WRITE(*,110);READ(*,*)F(2)
        WRITE(11,121)X(1),F(1)
        WRITE(11,120)F(2)
        CALL POLYLIN(1,X,F)
    CASE(2)
        STRING='linear2.res'
        CALL RESULTFILE(STRING)
        WRITE(11,101)
        WRITE(*,101)
        WRITE(*,111);READ(*,*)X(1),F(1)
        WRITE(*,112);READ(*,*)X(2),F(2)
        WRITE(11,121)X(1),F(1)
        WRITE(11,122)X(2),F(2)
        CALL POLYLIN(2,X,F)
    CASE(3)
        STRING='quadratic2.res'
        CALL RESULTFILE(STRING)
        WRITE(11,102)
        WRITE(*,102)
        WRITE(*,111);READ(*,*)X(1),F(1)
        WRITE(*,110);READ(*,*)F(3)
        WRITE(*,112);READ(*,*)X(2),F(2)
        WRITE(11,121)X(1),F(1)
        WRITE(11,120)F(3)
        WRITE(11,122)X(2),F(2)
        CALL POLYQUAD(2,X,F,XMIN,FMIN)
    CASE(4)
        STRING='quadratic3.res'
        CALL RESULTFILE(STRING)
        WRITE(11,102)
        WRITE(*,102)
        WRITE(*,111);READ(*,*)X(1),F(1)
        WRITE(*,112);READ(*,*)X(2),F(2)
        WRITE(*,113);READ(*,*)X(3),F(3)
        WRITE(11,121)X(1),F(1)
        WRITE(11,122)X(2),F(2)
        WRITE(11,123)X(3),F(3)
        CALL POLYQUAD(3,X,F,XMIN,FMIN)
    CASE(5)
        STRING='cubic3.res'
        CALL RESULTFILE(STRING)
        WRITE(11,103)
        WRITE(*,103)
        WRITE(*,110);READ(*,*)F(4)
        WRITE(*,111);READ(*,*)X(1),F(1)
        WRITE(*,112);READ(*,*)X(2),F(2)
        WRITE(*,113);READ(*,*)X(3),F(3)
        WRITE(11,121)X(1),F(1)
        WRITE(11,120)F(4)
        WRITE(11,122)X(2),F(2)
        WRITE(11,123)X(3),F(3)
        CALL POLYCUBIC(3,X,F,XMIN,FMIN)
    CASE(6)
        STRING='cubic4.res'
        CALL RESULTFILE(STRING)
        WRITE(11,103)
        WRITE(*,103)
        WRITE(*,111);READ(*,*)X(1),F(1)
        WRITE(*,112);READ(*,*)X(2),F(2)
        WRITE(*,113);READ(*,*)X(3),F(3)
        WRITE(*,114);READ(*,*)X(4),F(4)
        WRITE(11,121)X(1),F(1)
        WRITE(11,122)X(2),F(2)
        WRITE(11,123)X(3),F(3)
        WRITE(11,124)X(4),F(4)
        CALL POLYCUBIC(4,X,F,XMIN,FMIN)
    CASE DEFAULT
        WRITE(*,*)'Invalid option!'
        GOTO 10    
    END SELECT
    CLOSE(UNIT=11)
END SUBROUTINE