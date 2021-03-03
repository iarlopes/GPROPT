! Sub-routine that computes the Hessian matrix of a function,
! through the finite difference method.
! Igor Lopes, February 2015  
SUBROUTINE HESSIAN(H,X,NDIM)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) R4 /4.0D0/
    REAL(8) PERT /1.0D-4/
    ! Arguments
    INTEGER     :: NDIM
    REAL(8),DIMENSION(NDIM) :: X
    REAL(8),DIMENSION(NDIM,NDIM) :: H
    ! Locals
    INTEGER :: I,J
    REAL(8) :: FO,F1,F2,F3,F4,EVALFUNC
    REAL(8),DIMENSION(NDIM) :: X1,X2,X3,X4
    ! Algorithm
    FO=EVALFUNC(X,NDIM)
    DO I=1,NDIM
        DO J=1,NDIM
            IF(I.EQ.J)THEN
                X1=X
                X2=X
                X1(I)=X(I)+PERT
                X2(I)=X(I)-PERT
                F1=EVALFUNC(X1,NDIM)
                F2=EVALFUNC(X2,NDIM)
                H(I,I)=(F1-R2*FO+F2)/(PERT*PERT)
            ELSE
                X1=X
                X2=X
                X3=X
                X4=X
                X1(I)=X(I)+PERT
                X1(J)=X(J)+PERT
                X2(I)=X(I)+PERT
                X2(J)=X(J)-PERT
                X3(I)=X(I)-PERT
                X3(J)=X(J)+PERT
                X4(I)=X(I)-PERT
                X4(J)=X(J)-PERT
                F1=EVALFUNC(X1,NDIM)
                F2=EVALFUNC(X2,NDIM)
                F3=EVALFUNC(X3,NDIM)
                F4=EVALFUNC(X4,NDIM)
                H(I,J)=(F1-F2-F3+F4)/(R4*PERT*PERT)
            ENDIF
        ENDDO
    ENDDO
END SUBROUTINE