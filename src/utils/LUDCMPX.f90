SUBROUTINE LUDCMPX (A , INDX , D , N , NP , ERROR )
IMPLICIT NONE
!PARAMETER DECLARATION
INTEGER, PARAMETER:: NMAX=50
!DATA DECLARATION
REAL(8) R0    /0.0D0/
REAL(8) R1    /1.0D0/
REAL(8) TINY  /1.0D-19/
CHARACTER	  NAME*6
DATA NAME/'LUDCMP'/
!SCALAR VARIABLES FROM ARGUMENTS
INTEGER N, NP
LOGICAL ERROR
REAL(8) D
!ARRAYS FROM ARGUMENTS
REAL(8), DIMENSION(NP,NP):: A
INTEGER, DIMENSION(N):: INDX
!LOCAL SCALAR VARIABLES
INTEGER I, J, K, IMAX
REAL(8) AAMAX, SUM, DUM
!LOCAL ARRAYS
REAL(8), DIMENSION(N):: VV
!INTIALIZE LOCAL SCALAR VARIABLES
I=0      ; J=0    ; K=0    ; IMAX=0
AAMAX=R0 ; SUM=R0 ; DUM=R0 ; VV=R0
!***********************************************************************
! Routine to do LU decomposition
! See NUMERICAL RECIPES p35.
!*ACRONYM
! LU_DeCOmPosition
!*DESCRIPTION
!*EXTERNAL
! Arrays
!=A      - LU decomposed matrix
!=INDX   - Permutation vector
! Variables
!=D      - Row interchange indicator
! N      - Size of the problem
! NP     - Physical size of A matrix ( NP >= N )
! ERROR  - Error flag
!*INTERNAL
! Arrays
! VV     - store the implicit scalling of each rows
!***********************************************************************
! Loop over each rows to get the implicit scalling factors
! --------------------------------------------------------
DO 12 I=1,N
	AAMAX=R0
        DO 11 J=1,N
        	IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
   	11   CONTINUE
        IF(AAMAX.EQ.R0)THEN
		ERROR=.TRUE.
          	GOTO 999
        ENDIF
        VV(I)=R1/AAMAX
12 CONTINUE
! Loop over columns of Crout's methods
! ------------------------------------
DO 19 J=1,N
	IF (J.GT.1) THEN
! Calculate upper triangle components Bij see equation(2.3.12), except i=j
        	DO 14 I=1,J-1
            		SUM=A(I,J)
            		IF (I.GT.1)THEN
              			DO 13 K=1,I-1
                			IF(ABS(A(I,K)).LT.TINY.AND.ABS(A(K,J)).LT.TINY)GOTO 13
! If A(I,K) and A(K,J) are very small, skip the following line
                			SUM=SUM-A(I,K)*A(K,J)
   				13         CONTINUE
              			A(I,J)=SUM
            		ENDIF
   		14     CONTINUE
	ENDIF
! Initialize for the search for largest pivot element
        AAMAX=R0
! Calculate diagonal components Bij see equation(2.3.12) ,i=j, and
! lower triangle components aij see equation 2.3.13, i=j+1,..N
        DO 16 I=J,N
        	SUM=A(I,J)
          	IF (J.GT.1)THEN
            		DO 15 K=1,J-1
              			IF(ABS(A(I,K)).LT.TINY.AND.ABS(A(K,J)).LT.TINY)GOTO 15
! If A(I,K) and A(K,J) are very small, skip the following line.
              			SUM=SUM-A(I,K)*A(K,J)
   			15       CONTINUE
            		A(I,J)=SUM
          	ENDIF
! Find largest pivot element
          	DUM=VV(I)*ABS(SUM)
          	IF (DUM.GE.AAMAX) THEN
            		IMAX=I
            		AAMAX=DUM
          	ENDIF
   	16   CONTINUE
! Interchange rows
        IF(J.NE.IMAX)THEN
        	DO 17 K=1,N
            		DUM=A(IMAX,K)
            		A(IMAX,K)=A(J,K)
            		A(J,K)=DUM
   		17     CONTINUE
! Change sign of row channge indicator and interchange scale factor
          	D=-D
          	VV(IMAX)=VV(J)
	ENDIF
! record interchange row number
        INDX(J)=IMAX
! Divided by the pivot element and get final aij values, see equation 2.3.13,
! i=j+1,..N
        IF(J.NE.N)THEN
        	IF(ABS(A(J,J)).LT.TINY)A(J,J)=TINY
          	DUM=R1/A(J,J)
          	DO 18 I=J+1,N
            		A(I,J)=A(I,J)*DUM
   		18     CONTINUE
	ENDIF
! Go back for next column
19 CONTINUE
IF(ABS(A(N,N)).LT.TINY)A(N,N)=TINY
999 CONTINUE
RETURN
END
