! Main routine for genetic algorithms for optimization.
! Igor Lopes, February 2015
SUBROUTINE GENETIC(NDIM)
    IMPLICIT NONE
    ! Parameters
    REAL(8) R0 /0.0D0/
    REAL(8) R1 /1.0D0/
    REAL(8) R2 /2.0D0/
    ! Argument
    INTEGER NDIM
    ! Locals
    INTEGER     :: IBIT, ICOUNT, IDIM, IPOP, IGEN, NPOP, NTBIT, MGEN,&
                   NGENT, NSONS, NSURV, NTOPOUT, NTOPT
    INTEGER,DIMENSION(NDIM) :: NBITS
    REAL(8),DIMENSION(NDIM) :: XMIN, XMAX, PRECIS
    REAL(8),ALLOCATABLE,DIMENSION(:)    :: VFITO, VFIT, OBJFUN
    INTEGER,ALLOCATABLE,DIMENSION(:,:)  :: POPO, POP
    REAL(8),ALLOCATABLE,DIMENSION(:,:)  :: XVEC
    CHARACTER   :: STRING*256
    REAL(8)     :: RANDN, SRATE, PMUT
    LOGICAL     :: THEEND
    !
100 FORMAT('Define limits for each variable:')
105 FORMAT('Xmin_',I3,', Xmax_',I3)
107 FORMAT(' i ',2X,' Xmin_i ',2X,' Xmax_i ',2X,'nr.bits',2X,'Precision')
108 FORMAT(I3,2X,F8.3,2X,F8.3,3X,I5,3X,F8.6)    
110 FORMAT('Define number of individuals in the population:')
115 FORMAT('Define precision for each variable:')
117 FORMAT('nr. bits per chromossome:',I5,/'nr. chromossomes/individuals in the population:',I5,&
                    /'total nr. of bits in population:',I5)
120 FORMAT('X_',I3,':')
125 FORMAT('Define maximum number of generations:')
127 FORMAT('Maximum number of generations:',I5)
130 FORMAT('Define survival rate (0-1):')
135 FORMAT('No offspring with this survival rate...')
137 FORMAT('Survival rate:',F8.5,/'nr. of surviving individuals:',I5,&
                    /'nr. of sons in new population:',I5)
140 FORMAT('Define mutation probability (0-1):')
143 FORMAT('Mutation probability:',F8.5,/'Expected nr. of mutated bits:',I5)    
145 FORMAT('Value is outside admissible range...')
150 FORMAT('Define nr. of TOP individuals for ouput:')    
153 FORMAT('nr TOP cannot be greater than nr individuals...')
155 FORMAT('Define nr of TOP individuals and nr of generations'/&
                'without changes in TOP for termination criteria:')
157 FORMAT('Termination criteria: i) Maximum number of generations,'/&
            'ii) Reduction of diversity in population - every locus is has more than 90% of the same allele,'/&
            'iii) TOP ',I3,' does not change for ',I3,' generations.')
160 FORMAT('#GEN',2X,'TOP',6X,<NDIM>('X_'I3,10X),'obj.Func',3X,'Fitness')
163 FORMAT('-------------------------------------------------------------------------')    
    ! Initialize
    CALL RANDOM_SEED
    STRING='geneticalgorithm.res'
    CALL RESULTFILE(STRING)
    ICOUNT=0
    POPO=0
    VFITO=R0
    THEEND=.FALSE.
    ! Algorithm
    WRITE(*,100)
    DO IDIM=1,NDIM
        WRITE(*,105)IDIM,IDIM
        READ(*,*)XMIN(IDIM),XMAX(IDIM)
    ENDDO
    WRITE(*,115)
    NTBIT=0
    WRITE(11,107)
    DO IDIM=1,NDIM
        WRITE(*,120)IDIM
        READ(*,*)PRECIS(IDIM)
        NBITS(IDIM)=CEILING(LOG(R1+(XMAX(IDIM)-XMIN(IDIM))/PRECIS(IDIM))/LOG(R2))
        NTBIT=NTBIT+NBITS(IDIM)
        WRITE(11,108)IDIM,XMIN(IDIM),XMAX(IDIM),NBITS(IDIM),PRECIS(IDIM)
    ENDDO
    WRITE(*,110)
    READ(*,*)NPOP
    WRITE(11,117)NTBIT,NPOP,NTBIT*NPOP
    WRITE(*,125)
    READ(*,*)MGEN
    WRITE(11,127)MGEN
10  WRITE(*,130)
    READ(*,*)SRATE
    IF(SRATE.LT.R0.OR.SRATE.GT.R1)THEN
        WRITE(*,145)
        GOTO 10
    ENDIF
    NSURV=INT(SRATE*NPOP)
    NSONS=NPOP-NSURV
    IF(NSONS.EQ.0)THEN
        WRITE(*,135)
        GOTO 10
    ENDIF
    WRITE(11,137)SRATE,NSURV,NSONS
15  WRITE(*,140)
    READ(*,*)PMUT
    IF(PMUT.LT.R0.OR.PMUT.GT.R1)THEN
        WRITE(*,145)
        GOTO 15
    ENDIF
    WRITE(11,143)PMUT,INT(PMUT*NTBIT*NSONS)
20  WRITE(*,150)
    READ(*,*)NTOPOUT
    IF(NTOPOUT.GT.NPOP)THEN
        WRITE(*,153)
        GOTO 20
    ENDIF
25  WRITE(*,155)
    READ(*,*)NTOPT,NGENT
    IF(NTOPOUT.GT.NPOP)THEN
        WRITE(*,153)
        GOTO 25
    ENDIF
    WRITE(*,157)NTOPT,NGENT
    WRITE(11,157)NTOPT,NGENT
    WRITE(11,163);WRITE(*,163)
    WRITE(11,160)(IDIM,IDIM=1,NDIM);WRITE(*,160)(IDIM,IDIM=1,NDIM)
    WRITE(11,163);WRITE(*,163)
    !
    ALLOCATE(POPO(NPOP,NTBIT))
    ALLOCATE(POP(NPOP,NTBIT))
    ALLOCATE(VFIT(NPOP))
    ALLOCATE(VFITO(NPOP))
    ALLOCATE(OBJFUN(NPOP))
    ALLOCATE(XVEC(NDIM,NPOP))
    DO IBIT=1,NTBIT
        DO IPOP=1,NPOP
            CALL RANDOM_NUMBER(RANDN)
            POP(IPOP,IBIT)=NINT(RANDN)
        ENDDO
    ENDDO
    !
    CALL FITNESS(XMIN,XMAX,XVEC,POP,VFIT,OBJFUN,NBITS,NDIM,NPOP,NTBIT)
    CALL STATISTICS(0,XVEC,VFIT,OBJFUN,NDIM,NPOP,NTOPOUT)
    WRITE(11,163);WRITE(*,163)
    DO IGEN=1,MGEN
        POPO=POP
        VFITO=VFIT
        ! Survival of the fittest is implicit
        ! Reproduction
        CALL REPRODUCTION(POPO,POP,VFITO,VFIT,NSURV,NPOP,NTBIT)
        ! Mutation
        IF(PMUT.NE.R0)CALL MUTATION(POP,NPOP,PMUT,NSURV,NTBIT)
        WRITE(11,163);WRITE(*,163)
        !
        CALL FITNESS(XMIN,XMAX,XVEC,POP,VFIT,OBJFUN,NBITS,NDIM,NPOP,NTBIT)
        !
        CALL STATISTICS(IGEN,XVEC,VFIT,OBJFUN,NDIM,NPOP,NTOPOUT)
        WRITE(11,163);WRITE(*,163)
        CALL TERMINATION(POPO,POP,NPOP,NTBIT,NTOPT,NGENT,ICOUNT,THEEND)
        IF(THEEND)GOTO 99
    ENDDO
    
99  CLOSE(UNIT=11)
    
END SUBROUTINE