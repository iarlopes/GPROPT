DOUBLE PRECISION FUNCTION EVALFUNC(X,N)
  ! Interface to evaluate the objective function,
  ! for a design variable X with multiple N dimensions.
  !
  ! Igor A. Rodrigues Lopes
  ! 2015: Initial coding
  ! Mar 2021: Re-organised code with separate routines for different examples
  IMPLICIT NONE
  ! Parameters
  REAL(8) RP5 /0.5D0/
  REAL(8) R0 /0.0D0/
  REAL(8) R1 /1.0D0/
  REAL(8) R2 /2.0D0/
  REAL(8) R3 /3.0D0/
  ! Arguments
  INTEGER :: N
  REAL(8),DIMENSION(N) :: X
  ! ========================================================================================
  ! Select the appropriate problem by comment/uncommenting the appropriate lines
  ! call purelyQuadratic(EVALFUNC, X, N)
  call springs2(EVALFUNC, X, N)
  ! call springs6(EVALFUNC, X, N)
  ! call hyperelasticParameters(EVALFUNC, X, N)
  ! ========================================================================================
CONTAINS
  ! ========================================================================================
  subroutine purelyQuadratic(FValue, x, n)
    ! A purely quadratic function
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: x
    real(8), intent(out) :: FValue
    !
    FValue = X(1)*X(1)-R3*X(1)*X(2)+R2*X(2)*R2*X(2)+X(1)-X(2)
  end subroutine purelyQuadratic
  ! ========================================================================================
  subroutine springs2 (energy, x, ndim)
    ! Minimization of the potential energy in a system with 2 springs (Vanderplaats, 1984).
    !
    ! Igor A. Rodrigues Lopes
    ! 2015: Initial coding
    ! ======================================================================================
    implicit none
    integer, intent(in)  :: ndim
    real(8), dimension(ndim), intent(in)  :: x
    real(8), intent(out) :: energy
    REAL(8) :: A1,A2,K1,K2,L1,L2,P1,P2 ! 2 spring
    ! Algorithm - 2 spring problem
    K1=8.0E0  ! Spring constant 1 (N/cm)
    K2=1.0E0  ! Spring constant 2 (N/cm)
    L1=10.0E0   ! Initial length 1 (cm)
    L2=10.0E0   ! Initial length 2 (cm)
    P1=5.0E0  ! Force 1 (N)
    P2=5.0E0  ! Force 2 (N)
    !
    A1=DSQRT(X(1)*X(1)+(L1-X(2))**2)
    A2=DSQRT(X(1)*X(1)+(L2+X(2))**2)
    energy=RP5*K1*(A1-L1)**2+RP5*K2*(A2-L2)**2-P1*X(1)-P2*X(2)
  end subroutine springs2
  ! ========================================================================================
  subroutine springs6 (energy, x, ndim)
    ! Minimization of the potential energy in a system with 6 springs (Vanderplaats, 1984).
    !
    ! Igor A. Rodrigues Lopes
    ! 2015: Initial coding
    ! ======================================================================================
    implicit none
    integer, intent(in)  :: ndim
    real(8), dimension(ndim), intent(in)  :: x
    real(8), intent(out) :: energy   
    !
    REAL(8) R10 /10.0D0/
    REAL(8) R50 /50.0D0/
    INTEGER,PARAMETER :: NW=5, NS=6
    REAL(8),DIMENSION(NW) :: W
    REAL(8),DIMENSION(NS) :: K,DL,LO
    REAL(8) :: XO, YO, XN, YN
    integer :: i
    ! 6 spring problem
    XO=R0
    XN=60.0E0
    YO=R0
    YN=R0
    LO=10.0E0
    
    DO I=1,NW
       W(I)=R50*I
    ENDDO
    DO I=1,NS
       K(I)=R50*R10+R2*R10*R10*(NW/R3-I)*(NW/R3-I)
    ENDDO
    
    DL(1)=DSQRT((X(1)-XO)**R2+(X(2)-YO)**R2)-LO(1)
    DO I=2,NS-1
       DL(I)=DSQRT((X(2*I-1)-X(2*I-3))**R2+(X(2*I)-X(2*I-2))**R2)-LO(I)
    ENDDO
    DL(NS)=DSQRT((XN-X(2*NS-3))**R2+(YN-X(2*NS-2))**R2)-LO(NS)
    
    energy=R0
    DO I=1,NW
       energy=energy+W(I)*X(2*I)
    ENDDO
    DO I=1,NS
       energy=energy+RP5*K(I)*DL(I)*DL(I)
    ENDDO
  end subroutine springs6
  ! ========================================================================================
  subroutine hyperelasticParameters (FValue, x, ndim)
    ! Parameter identification for Ogden's hyperelastic model, given the stress response
    ! for uniaxial and biaxial deformation.
    !
    ! Igor A. Rodrigues Lopes
    ! 2015: Initial coding
    ! ======================================================================================
    implicit none
    integer, intent(in)  :: ndim
    real(8), dimension(ndim), intent(in)  :: x
    real(8), intent(out) :: FValue
    REAL(8),DIMENSION(10)   :: P1MUA, P2MUA, PMEB, LAMBDA, P1UA, P2UA, PEB
    REAL(8),DIMENSION(4)  :: AP, MUP
    REAL(8)         :: FUA, FEB, FC
    INTEGER         :: K, I
    ! ======================================================================================
    ! Hyperelastic parameter identification
    P1MUA(1)=49.1713758047908     
    P1MUA(2)=95.5020609515927      
    P1MUA(3)=139.194048890133     
    P1MUA(4)=180.435931790392     
    P1MUA(5)=219.403962336700     
    P1MUA(6)=256.263042125065     
    P1MUA(7)=291.167643234921     
    P1MUA(8)=324.262668859696     
    P1MUA(9)=355.684258151587     
    P1MUA(10)=385.560539643794     
    
    P2MUA(1)=24.1807151003576      
    P2MUA(2)=47.3466758278295
    P2MUA(3)=69.5500055375245       
    P2MUA(4)=90.8405113342239      
    P2MUA(5)=111.265817399249        
    P2MUA(6)=130.871490139650       
    P2MUA(7)=149.701155567176       
    P2MUA(8)=167.796609139805       
    P2MUA(9)=185.197918133637       
    P2MUA(10)=201.943516470978       
    
    PMEB(1)=71.7371822572818       
    PMEB(2)=136.769883566410       
    PMEB(3)=195.863728736562      
    PMEB(4)=249.697571912768      
    PMEB(5)=298.874149043269       
    PMEB(6)=343.929420777042      
    PMEB(7)=385.340780268923       
    PMEB(8)=423.534282861733       
    PMEB(9)=458.891010936387     
    PMEB(10)=491.752691255658     
    
    LAMBDA(1)=1.02658363130423
    LAMBDA(2)=1.05387395206178
    LAMBDA(3)=1.08188974864453
    LAMBDA(4)=1.11065030683432
    LAMBDA(5)=1.14017542509914
    LAMBDA(6)=1.17048542822212
    LAMBDA(7)=1.20160118129295
    LAMBDA(8)=1.23354410407117
    LAMBDA(9)=1.26633618573131
    LAMBDA(10)=1.30000000000000
    
    IF(NDIM.NE.8.AND.NDIM.NE.6)THEN
      WRITE(*,*)'nr. of parameters for Ogden model must be 6 or 8, NDIM=',NDIM
      stop 1
    ENDIF
    DO I=1,NDIM/2
      AP(I)=X(I)
      MUP(I)=X(NDIM/2+I)
    ENDDO
    !
    P1UA=R0
    P2UA=R0
    PEB=R0
    FUA=R0
    FEB=R0
    FC=R0
    DO K=1,10
      DO I=1,NDIM/2
        P1UA(K)=P1UA(K)+MUP(I)*(LAMBDA(K)**(AP(I)-R1)-LAMBDA(K)**(-AP(I)-R1))
        P2UA(K)=P2UA(K)+MUP(I)*(R1-LAMBDA(K)**(-AP(I)))
        PEB(K)=PEB(K)+MUP(I)*(LAMBDA(K)**(AP(I)-R1)-LAMBDA(K)**(-R2*AP(I)-R1))
      ENDDO
      FUA=FUA+(P1UA(K)/P1MUA(K)-R1)**2+(P2UA(K)/P2MUA(K)-R1)**2
      FEB=FEB+(PEB(K)/PMEB(K)-R1)**2
    ENDDO
    FValue=FUA+FEB
  end subroutine hyperelasticParameters
  ! ========================================================================================
END FUNCTION