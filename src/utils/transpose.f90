! Subroutine for determine matrix transpose AT=A^T
!
! Igor Lopes, February 2015
SUBROUTINE TRANSPOSE(A,AT,NROWA,NCOLA)
    IMPLICIT NONE
    ! Arguments
    INTEGER :: NROWA , NCOLA
    REAL(8), DIMENSION(NROWA,NCOLA)  :: A
    REAL(8), DIMENSION(NCOLA, NROWA) :: AT
    ! Locals
    INTEGER I, J
    ! Begin algorithm
    DO J=1,NCOLA
	    DO I=1,NROWA
		    AT(J,I)=A(I,J)
	    ENDDO
    ENDDO
END SUBROUTINE
