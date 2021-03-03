! Sub-routine that computes the gradient of the objective function.
! It can be computed through the finite difference method, or defined
! by its analytical expression.
! Igor Lopes, February 2015
!********************************************************************
! This specific case refers to the equilibrium of a 2 spring system,
! presented on page 72 of "Numerical Optimization Techniques for 
! Engineering Design", by Vanderplaats.    
SUBROUTINE GRADIENT(GRAD,X,NDIM)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    REAL(8) R3 /3.0D0/
    REAL(8) PERT /1.0D-4/
    ! Arguments
    INTEGER     :: NDIM
    REAL(8),DIMENSION(NDIM) :: X,GRAD
    ! Locals
    INTEGER :: I
    REAL(8) :: A1,A2,K1,K2,L1,L2,P1,P2,F1,F,EVALFUNC
    REAL(8),DIMENSION(NDIM) :: X1
    ! Algorithm
    K1=8.0E0    ! Spring constant 1 (N/cm)
    K2=1.0E0    ! Spring constant 2 (N/cm)
    L1=10.0E0   ! Initial length 1 (cm)
    L2=10.0E0   ! Initial length 2 (cm)
    P1=5.0E0    ! Force 1 (N)
    P2=5.0E0    ! Force 2 (N)
    !
    A1=DSQRT(X(1)*X(1)+(L1-X(2))**2)
    A2=DSQRT(X(1)*X(1)+(L2+X(2))**2)
    GRAD(1)=K1*X(1)*(A1-L1)/A1+K2*X(1)*(A2-L2)/A2-P1
    GRAD(2)=K1*(X(2)-L1)*(A1-L1)/A1+K2*(X(2)+L2)*(A2-L2)/A2-P2
    ! Alternative: finite differences method
    !F=EVALFUNC(X,NDIM)
    !DO I=1,NDIM
    !    X1=X
    !    X1(I)=X1(I)+PERT
    !    F1=EVALFUNC(X1,NDIM)
    !    GRAD(I)=(F1-F)/PERT
    !ENDDO
END SUBROUTINE