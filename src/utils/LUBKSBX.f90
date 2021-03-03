SUBROUTINE LUBKSBX (A , B , INDX ,  N , NP)
IMPLICIT NONE
!DATA DECLARATION
REAL(8) R0    /0.0D0/
CHARACTER	  NAME*6
DATA NAME/'LUBKSB'/
!SCALAR VARIABLES FROM ARGUMENTS
INTEGER N, NP
!ARRAYS FROM ARGUMENTS
REAL(8), DIMENSION(NP,NP):: A
REAL(8), DIMENSION(N):: B
INTEGER, DIMENSION(N):: INDX
!LOCAL SCALAR VARIABLES
INTEGER II, I, LL, J
REAL(8) SUM
!INITIALIZE LOCAL SCALAR VARIABLES
II=0;   I=0;    LL=0;   J=0;
SUM=R0;
!***********************************************************************
! Routine to solve the set of N linear equations AX=B
! See NUMERICAL RECIPES p36.
!*ACRONYM
! LU_BacK_SuBustitutions
!*DESCRIPTION
!*EXTERNAL
! Arrays
! A      - LU decomposed matrix
!=B      - Right hand side matrix as input and stored solutions as output
! INDX   - Permutation vector
! Variables
! N      - Size of the problem
! NP     - Physical size of A matrix
!***********************************************************************
II=0
! Do forward substitution, equation (2.3.6)
DO 12 I=1,N
	LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
        	DO 11 J=II,I-1
            		SUM=SUM-A(I,J)*B(J)
   		11     CONTINUE
        ELSEIF (SUM.NE.R0) THEN
          	II=I
        ENDIF
        B(I)=SUM
12 CONTINUE
! Do the backsubstitution, equation (2.3.7)
DO 14 I=N,1,-1
	SUM=B(I)
        IF(I.LT.N)THEN
		DO 13 J=I+1,N
            		SUM=SUM-A(I,J)*B(J)
   		13     CONTINUE
        ENDIF
! Store a component of the solution vector Xi
        B(I)=SUM/A(I,I)
14 CONTINUE
RETURN
END
