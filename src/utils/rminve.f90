SUBROUTINE RMINVE (A , AI , NSIZE , ERROR)
IMPLICIT NONE
!PARAMETER DECLARATION
INTEGER, PARAMETER:: MSIZE=500
!DATA DECLARATION
REAL(8) R0    /0.0D0/
REAL(8) R1    /1.0D0/
CHARACTER*6 NAME
DATA NAME/'RMINVE'/
!SCALAR VARIABLES FROM ARGUMENTS
INTEGER NSIZE
LOGICAL ERROR
!ARRAYS FROM ARGUMENTS
REAL(8), DIMENSION(NSIZE,NSIZE):: A
REAL(8), DIMENSION(NSIZE,NSIZE):: AI
!LOCAL SCALAR VARIABLES
INTEGER ISIZE, JSIZE
REAL(8) DUMMY, DETA
!LOCAL ARRAYS
INTEGER, DIMENSION(NSIZE):: INDX
!INITIALIZE LOCAL VARIABLES
ISIZE=0  ; JSIZE=0
DUMMY=R0 ; DETA=R0
INDX=0
!***********************************************************************
! Evaluate inverse matrix without determinant
!*Arrays
! A     - A matrix
!=AI    - inversed A matrix
!*Variables
! NSIZE - Size of the matrix
!=ERROR - Error flag
!***********************************************************************
! Set up identity matrix
! ----------------------
ERROR=.FALSE.
DO 20 ISIZE=1,NSIZE
	DO 10 JSIZE=1,NSIZE
        	AI(ISIZE,JSIZE)=R0
   	10   CONTINUE
        AI(ISIZE,ISIZE)=R1
20 CONTINUE
! LU decompose the matrix
CALL LUDCMPX(A    ,INDX  ,DUMMY,NSIZE ,NSIZE ,ERROR )
DETA=DUMMY
! Find inverse by columns
DO 30 ISIZE=1,NSIZE
	CALL LUBKSBX(A     ,AI(1,ISIZE),INDX  ,NSIZE ,NSIZE )
! Obtain determinant
        DETA=DETA*A(ISIZE,ISIZE)
30 CONTINUE
999 CONTINUE
RETURN
END