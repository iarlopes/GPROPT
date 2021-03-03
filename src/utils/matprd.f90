! Function for matrix product C=A*B
! Optimized memory access pattern is considered.
! In Fortran matrices should be accessed column by column.
! Igor Lopes, February 2015
SUBROUTINE MATPRD(A,B,C,NROWA,NCOLA,NCOLB)
    IMPLICIT NONE
    ! Parameter
    REAL(8) R0 /0.0D0/
    ! Arguments
    INTEGER :: NROWA , NCOLA , NCOLB
    REAL(8), DIMENSION(NROWA,NCOLA)  :: A
    REAL(8), DIMENSION(NCOLA, NCOLB) :: B
    REAL(8), DIMENSION(NROWA,NCOLB)  :: C
    ! Locals
    INTEGER I, J, K
    ! Begin algorithm
    C=R0
    DO J=1,NCOLB
	    DO K=1,NCOLA
		    DO I=1,NROWA
			    C(I,J)=C(I,J)+A(I,K)*B(K,J)
		    ENDDO
	    ENDDO
    ENDDO
END SUBROUTINE