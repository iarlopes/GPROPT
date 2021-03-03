! This routine opens a file whose name is given 
! in STRING variable
!    
! Igor Lopes, February 2015   
SUBROUTINE RESULTFILE(STRING)
    IMPLICIT NONE
    CHARACTER STRING*256
    
    OPEN(UNIT=11,FILE=STRING,STATUS="UNKNOWN")
END SUBROUTINE  