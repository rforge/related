! Estimating relatedness by using trio genotypes. Other estimators
!   are also included for comparison
! 95% CI estimated by bootstrapping over loci
! Compiling notes:
!   (1) When compiled by FTN95, do not use /OPTIMISE option
!   (2) Better use Ifort:
!       ifort trior.f90 /O2 /Qip /Qipo /Qinline-factor=100 /Qunroll-aggressive /link
!   A missing allele in a genotype is denoted by any integer number <1
!
!Compile command for dianostics:
!g95 *.f90 -Wall -pedantic -fbounds-check -ftrace=full -std=f2003

MODULE Variables
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)              ! Double Precision
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)    ! -2147483648 ~ 2147483647
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)    ! -32768 ~ 32767
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)    ! -128 ~ 127
  INTEGER, PARAMETER :: MUL = 100000000
  REAL(DP), PARAMETER :: Small = 1.E-300_DP
  INTEGER :: NumLoci, SmplSize, MaxNumAllele, BootsSamples, NumTrios, Constrain, &
    KnownAlleleFre, NDIM, EstimateCI95, I_SEED, I025(2), IVar(2), IDUM, NumLoci1
  INTEGER, Target :: MarkVar2(9), MarkVar(66), NumMarkVar2, NumMarkVar
  INTEGER(I4B) :: ThreadID, NumThreads, NumThreadsOMP
  INTEGER, Pointer :: MarkVarPtr(:), NumMarkVarPtr, LocusSelectPtr(:)
  INTEGER, ALLOCATABLE :: ObsNumAllele(:), RawAllele(:, :), &
    ALLELE(:,:), IBDConfig2(:)
  INTEGER(I1B), ALLOCATABLE :: GUnrelated(:, :, :)
  INTEGER, TARGET, ALLOCATABLE :: LocusSelect1(:), LocusSelect2(:)
  INTEGER(I2B), TARGET, ALLOCATABLE :: LocusRepeat1(:), LocusRepeat2(:)
  INTEGER(I2B), Pointer :: LocusRepeat(:)
  INTEGER(I1B),ALLOCATABLE :: IndivGType(:, :, :), RawAlleleNum(:)
  TARGET :: IndivGType
  INTEGER(I2B), ALLOCATABLE :: GConfig(:)
  INTEGER(I4B) :: TrioPointer, RMethod(7)
  LOGICAL :: CI95Q
  TYPE FSProbPrt
    Integer(I4B), Pointer :: R(:)
  END TYPE FSProbPrt
  TYPE(FSProbPrt), Allocatable :: RawAlleleID(:)
  REAL(DP) :: LogMul, Func3_Old, LogLValue
  REAL(DP), ALLOCATABLE :: ObsAlleleFre(:, :), pcom(:), &
    err(:), coe(:, :), Prob_Locus(:),  &
    Prob_Locus2(:), coeXij2(:), RecAlleleFre(:, :)
  REAL, ALLOCATABLE :: REst(:, :), FEst(:, :), RCI95(:, :, :), &
    DCI95(:, :, :, :), FCI95(:, :, :), Delta(:, :, :)
  CHARACTER(Len = 100) :: FName
  CHARACTER(Len = 20), ALLOCATABLE :: IndivID(:)
  LOGICAL :: AccountError
  CHARACTER(Len = 12) :: DATE, TIME
  INTEGER :: AllDyads
  CHARACTER(Len = 999) :: ErrorMessage(5)
END MODULE Variables
!
Module WangVariables
  Use Variables
  IMPLICIT NONE
  REAL(DP), Allocatable :: Fac(:, :), WEIGHT(:)
END Module WangVariables
!
MODULE IBDVariables
  USE VARIABLES
  IMPLICIT NONE
  INTEGER :: IBDValue
  INTEGER(I1B), Allocatable :: G(:, :)
  REAL(DP), Allocatable :: A(:, :)
  REAL(DP) :: FUNC0
END MODULE IBDVariables
!
subroutine related
  USE VARIABLES

  IMPLICIT NONE
!$   include 'mpif.h'
!$  Logical :: ExistQ
  INTEGER :: I, J, K, M, I1, I2, KK, N, ID(2, 99), Idx(2)
  INTEGER(I1B), Pointer :: GPtr(:)
  REAL :: FEst4(2), X, Y(2), RAN1
  REAL, ALLOCATABLE :: RR(:, :), DD(:, :, :), FF(:, :, :)
  Character(Len = 6) :: ThreadString
  EXTERNAL :: INITIATE, DataProcess, InbreedingMT, Wang, LynchRitland, &
    Ritland, Milligan, TrioEst, OUTPUT, Sort, SimuGtype, RAN1, WriteMidR, WriteMidCI

  ! check for user interrupts in the main loop
  external rchkusr

  CALL InitializeThreads
  CALL READ_DATA
  CALL INITIATE
  CALL DataProcess
  IF(RMethod(1) > 0) CALL SimuGtype
  IF(EstimateCI95 == 1) ALLOCATE(RR(BootsSamples, 7), &
    DD(BootsSamples, 2, 4), FF(BootsSamples, 1, 2 : 3))
  DO M = 1, SmplSize
    LocusSelectPtr => LocusSelect1
    GPtr => IndivGType(1, :, M)
    LocusRepeat => LocusRepeat1
    CI95Q = .FALSE.
    CALL InbreedingMT(M, FEst4)
    FEst(2 : 3, M) = FEst4(1 : 2)
    IF(EstimateCI95 == 1)THEN
      LocusSelectPtr => LocusSelect2
      LocusRepeat => LocusRepeat2
      CI95Q = .TRUE.
      LocusRepeat2(:) = 0_I2B
      DO I = 1, BootsSamples
        LoopMain1: DO
          DO J = 1, NumLoci
            LocusSelect2(J) = Floor(RAN1() * NumLoci + 1.0)
            LocusRepeat2(LocusSelect2(J)) = LocusRepeat2(LocusSelect2(J)) + 1_I2B
          END DO
          DO J = 1, NumLoci
            IF(GPtr(LocusSelectPtr(J)) > 0_I1B .AND. &
              ObsNumAllele(LocusSelectPtr(J)) > 1) EXIT LoopMain1
          END DO
        END DO LoopMain1
        CALL InbreedingMT(M, FEst4)
        FF(I, 1, 2 : 3) = FEst4(1 : 2)
      END DO
      DO I = 2, 3
        CALL Sort(BootsSamples, FF(:, 1, I))
        DO J = 1, 2
          FCI95(J,I,M) = FF(I025(J), 1, I)
        END DO
      END DO
    END IF
  END DO
  IF(Constrain /= 1) FORALL(I2 = 1 : SmplSize, I1 = 1 : 4 : 3) FEst(I1, I2) = 0.0
  IF(EstimateCI95 == 1 .AND. Constrain /= 1 .AND. ANY(RMethod(1 : 7 : 6) == 2))THEN
    DEALLOCATE(FF)
    IF(ALL(RMethod(1 : 7 : 6) == 2)) THEN
      Idx(1 : 2) = (/1, 2/)
    ELSE
      Idx(1 : 2) = Merge(1, 2, RMethod(1) == 2)
    END IF
    ALLOCATE(FF(BootsSamples, SmplSize, Idx(1) : Idx(2)))
    FF(:, :, :) = 0.
  END IF
  WRITE(ThreadString, '(I6)') ThreadID
  ThreadString = Trim(AdjustL(ThreadString))
  M = 0
  CALL WriteMidR(ThreadString, M, ID(:, :))
  IF(EstimateCI95 == 1) CALL WriteMidCI(ThreadString, M, ID(:, :))
  IF(NumThreads > 1) N = 0
  KK = SmplSize
  DO I1 = 1, SmplSize, 2 - AllDyads

    ! Check for user interrupts.
    call rchkusr

    IF(AllDyads == 0) KK = I1 + 1
    DO I2 = I1 + 1, KK
      IF(NumThreads > 1) THEN
        N = N + 1
        IF(MOD(N, NumThreads) /= ThreadID) CYCLE
!$      CALL MPI_AllReduce(n, J, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, K)
      END IF
      IF(M == 99) THEN
        CALL WriteMidR(ThreadString, M, ID(:, :))
        IF(EstimateCI95 == 1) CALL WriteMidCI(ThreadString, M, ID(:, :))
        M = 0
      END IF
      M = M + 1
      ID(1 : 2, M) = (/I1, I2/)
      LocusSelectPtr => LocusSelect1
      LocusRepeat => LocusRepeat1
      CI95Q = .FALSE.
      IF(ThreadID == 0 .AND. MOD(M, 9) == 0) WRITE(*,'(A,2I6,A)')'Dyad of individuals of ',I1,I2,&
             '    Calculating point estimate' // &
             Merge('       ', ' & CI95', EstimateCI95 /= 1)
      K = 0
      DO J = 1, NumLoci
        IF(IndivGType(1, J, I1) < 1_I1B .OR. &
           IndivGType(1, J, I2) < 1_I1B) CYCLE
        K = 1
        EXIT
      END DO
      IF(K == 0)THEN  !No locus that has both genotypes available
        REst(1 : 7, M) = 0
        Delta(1 : 2, M, 1 : 4) = 0
        IF(EstimateCI95 == 1)THEN
           RCI95(1 : 2, 1 : 7, M) = 0
           DCI95(1 : 2, 1 : 2, 1 : 4, M) = 0
        END IF
        CYCLE
      END IF
      IF(ALL(RMethod(2 : 3) == 0)) THEN
        REst(2 : 3, M) = 0.0
        Delta(1 : 2, M, 2) = 0.0
      ELSE
        CALL WANG(I1, I2, REst(2, M), REst(3, M), Y)
        Delta(1 : 2, M, 2) = Y(1 : 2)
        IF(RMethod(2) == 0) THEN
          Delta(1 : 2, M, 2) = 0.0
          REst(2, M) = 0.0
        ELSE IF(RMethod(3) == 0) THEN
          REst(3, M) = 0.0
        END IF
      END IF
      IF(RMethod(4) == 0) THEN
        REst(4, M) = 0.0
        Delta(1 : 2, M, 3) = 0.0
      ELSE
        CALL LynchRitland(I1, I2, REst(4, M), Y)
        Delta(1 : 2, M, 3) = Y(1 : 2)
      END IF
      IF(ALL(RMethod(5 : 6) == 0)) THEN
        REst(5 : 6, M) = 0.0
      ELSE
        CALL Ritland(I1, I2, REst(5, M), REst(6, M))
      END IF
      IF(RMethod(7) == 0) THEN
        REst(7, M) = 0.0
        Delta(1 : 2, M, 4) = 0.0
        FEst4(1 : 2) = 0.0
      ELSE
        CALL Milligan(I1, I2, REst(7, M), FEst4, Y)
        Delta(1 : 2, M, 4) = Y(1 : 2)
      END IF
      IF(Constrain /= 1)THEN
        FEst(4, I1) = FEst(4, I1) + FEst4(1)
        FEst(4, I2) = FEst(4, I2) + FEst4(2)
      END IF
      IF(RMethod(1) == 0) THEN
        REst(1, M) = 0.0
        Delta(1 : 2, M, 1) = 0.0
        FEst4(1 : 2) = 0.0
      ELSE
        CALL TrioEst(I1, I2, M, 1, NumTrios, 1, X, FEst4, Y)
        REst(1, M) = X
        Delta(1 : 2, M, 1) = Y(1 : 2)
      END IF
      IF(Constrain /= 1)THEN
        FEst(1, I1) = FEst(1, I1) + FEst4(1)
        FEst(1, I2) = FEst(1, I2) + FEst4(2)
      END IF
      IF(EstimateCI95 == 1)THEN   !CI95 from Bootstrapping
        K = TrioPointer
        LocusSelectPtr => LocusSelect2
        LocusRepeat => LocusRepeat2
        CI95Q = .TRUE.
        DO I = 1, BootsSamples
          LocusRepeat2(:) = 0_I2B
          LoopMain2 : DO
            DO J = 1, NumLoci
              LocusSelect2(J) = Floor(RAN1() * NumLoci + 1.0)
              LocusRepeat2(LocusSelect2(J)) = LocusRepeat2(LocusSelect2(J)) + 1_I2B
            END DO
            DO J = 1, NumLoci
              IF(IndivGType(1, LocusSelect2(J), I1) > 0_I1B .AND. &
                 IndivGType(1, LocusSelect2(J), I2) > 0_I1B) EXIT LoopMain2
            END DO
          END DO LoopMain2
          IF(ANY(RMethod(2 : 3) == 2)) THEN
            CALL WANG(I1, I2, RR(I, 2), RR(I, 3), Y)
            DD(I, 1 : 2, 2) = Y(1 : 2)
          END IF
          IF(RMethod(4) == 2) THEN
            CALL LynchRitland(I1, I2, RR(I, 4), Y)
            DD(I, 1 : 2, 3) = Y(1 : 2)
          END IF
          IF(ANY(RMethod(5 : 6) == 2)) CALL Ritland(I1, I2, RR(I, 5), RR(I, 6))
          IF(RMethod(7) == 2) THEN
            CALL Milligan(I1, I2, RR(I, 7), FEst4, Y)
            DD(I, 1 : 2, 4) = Y(1 : 2)
            IF(Constrain /= 1)THEN
              FF(I, I1, 2) = FF(I, I1, 2) + FEst4(1)
              FF(I, I2, 2) = FF(I, I2, 2) + FEst4(2)
            END IF
          END IF
          IF(RMethod(1) == 2) THEN
            CALL TrioEst(I1, I2, M, K, K, 0, RR(I, 1), FEst4, Y)
            DD(I, 1 : 2, 1) = Y(1 : 2)
            IF(Constrain /= 1)THEN
              FF(I, I1, 1) = FF(I, I1, 1) + FEst4(1)
              FF(I, I2, 1) = FF(I, I2, 1) + FEst4(2)
            END IF
          END IF
        END DO
        DO I = 1, 7
          IF(RMethod(I) == 2) THEN
            CALL Sort(BootsSamples, RR(:, I))
            DO J = 1, 2
              RCI95(J, I, M) = RR(I025(J), I)
            END DO
          ELSE
            RCI95(1 : 2, I, M) = 0.0
          END IF
        END DO
        DO I = 1, 2
          DO K = 1, 4
            J = Merge(K, K + 1 + (K - 3) * 2, K < 3)
            IF(RMethod(J) < 2) THEN
              DCI95(:, I, K, M) = 0.0
            ELSE
              CALL Sort(BootsSamples, DD(:, I, K))
              DO J = 1, 2
                DCI95(J, I, K, M) = DD(I025(J), I, K)
              END DO
            END IF
          END DO
        END DO
      END IF
    END DO
  END DO
  IF(M > 0) THEN
    CALL WriteMidR(ThreadString, M, ID(:, 1 : M))
    IF(EstimateCI95 == 1) CALL WriteMidCI(ThreadString, M, ID(:, 1 : M))
  END IF
  IF(Constrain /= 1)THEN
    X = Merge(1., 1. / (SmplSize - 1.), AllDyads == 0)
    FORALL(M = 1 : SmplSize, J = 1 : 4 : 3) FEst(J, M) = FEst(J, M) * X
    IF(EstimateCI95 == 1)THEN
      DO J = 1, 2
        I = Merge(1, 7, J == 1)
        DO M = 1, SmplSize
          IF(RMethod(I) < 2) THEN
            FCI95(:, J * J, M) = 0.0
          ELSE
            CALL Sort(BootsSamples, FF(:, M, J))
            FORALL(K = 1 : 2) FCI95(K, J * J, M) = FF(I025(K), M, J) * X
          END IF
        END DO
      END DO
    END IF
  END IF
  IF(EstimateCI95 == 1) DEALLOCATE(RR, DD)
  IF(Allocated(FF)) DEALLOCATE(FF)


!$ CALL MPI_Barrier(MPI_COMM_WORLD, J)
  IF(ThreadID == 0) then
    CALL OUTPUT
    ! Clean up global memory; execution no longer ends with related,
    ! so memory leaks matter!
    call clean_mem
  end if
!$  CALL MPI_INITIALIZED(ExistQ, J)
!$  IF(ExistQ) CALL MPI_FINALIZE(J)
END subroutine related
!
subroutine clean_mem
  use Variables
  use WangVariables
  use IBDVariables
  implicit none
  integer :: i

  ! fortunately, this code allocates no pointers, only assigns them.

  if (allocated(Fac))          deallocate(Fac)
  if (allocated(WEIGHT))       deallocate(WEIGHT)
  if (allocated(G))            deallocate(G)
  if (allocated(A))            deallocate(A)
  if (allocated(ObsNumAllele)) deallocate(ObsNumAllele)
  if (allocated(RawAllele))    deallocate(RawAllele)
  if (allocated(ALLELE))       deallocate(ALLELE)
  if (allocated(IBDConfig2))   deallocate(IBDConfig2)
  if (allocated(GUnrelated))   deallocate(GUnrelated)
  if (allocated(LocusSelect1)) deallocate(LocusSelect1)
  if (allocated(LocusSelect2)) deallocate(LocusSelect2)
  if (allocated(LocusRepeat1)) deallocate(LocusRepeat1)
  if (allocated(LocusRepeat2)) deallocate(LocusRepeat2)
  if (allocated(IndivGType))   deallocate(IndivGType)
  if (allocated(RawAlleleNum)) deallocate(RawAlleleNum)
  if (allocated(GConfig))      deallocate(GConfig)
  if (allocated(ObsAlleleFre)) deallocate(ObsAlleleFre)
  if (allocated(pcom))         deallocate(pcom)
  if (allocated(err))          deallocate(err)
  if (allocated(coe))          deallocate(coe)
  if (allocated(Prob_Locus))   deallocate(Prob_Locus)
  if (allocated(Prob_Locus2))  deallocate(Prob_Locus2)
  if (allocated(coeXij2))      deallocate(coeXij2)
  if (allocated(RecAlleleFre)) deallocate(RecAlleleFre)
  if (allocated(REst))         deallocate(REst)
  if (allocated(FEst))         deallocate(FEst)
  if (allocated(RCI95))        deallocate(RCI95)
  if (allocated(DCI95))        deallocate(DCI95)
  if (allocated(FCI95))        deallocate(FCI95)
  if (allocated(Delta))        deallocate(Delta)
  if (allocated(IndivID))      deallocate(IndivID)
  if (allocated(RawAlleleID)) then
    do i = 1, NumLoci
      deallocate(RawAlleleID(i) % R)
    end do
    deallocate(RawAlleleID)
  end if

end subroutine clean_mem
!
Subroutine InitializeThreads
  USE VARIABLES
  IMPLICIT NONE
!$   include 'mpif.h'
!$  INTEGER :: I, J, OMP_GET_NUM_PROCS
!$  EXTERNAL :: StopOnError
  ErrorMessage = ''
  ThreadID = 0
  NumThreads = 1
  NumThreadsOMP = 1
!$  NumThreadsOMP = OMP_GET_NUM_PROCS()
!$  CALL MPI_INIT(I)
!$  IF(I /= MPI_SUCCESS)THEN
!$    ErrorMessage(1) = 'Error starting MPI program. Terminating.'
!$    CALL MPI_ABORT(MPI_COMM_WORLD, J, I)
!$    CALL StopOnError(1, 'InitializeThreads')
!$  ELSE
!$    CALL MPI_COMM_RANK(MPI_COMM_WORLD, ThreadID, I)
!$    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NumThreads, I)
!$  END IF
END Subroutine InitializeThreads
!
SUBROUTINE StopOnError(IFlag, C)
  USE VARIABLES
  IMPLICIT NONE
!$   include 'mpif.h'
  Character(*), INTENT(IN) :: C
  INTEGER :: I, IFlag
!$  Logical :: ExistQ
  EXTERNAL :: DelayedExit
  external rexit
  IF(ThreadID == 0)THEN
    Open(10, File = 'ErrorMessage')
    IF(IFlag == 0) THEN
      WRITE(*, '(A)') 'Program stoped due to insufficient memory!'
      WRITE(10, '(A)') 'Program stoped due to insufficient memory!'
    ELSE
      DO I = 1, 5
        IF(Len_Trim(ErrorMessage(I)) < 1) CYCLE
        WRITE(*, '(A)') Trim(ErrorMessage(I))
        WRITE(10, '(A)') Trim(ErrorMessage(I))
      END DO
    END IF
    WRITE(*, '(A)') 'Program stoped in subroutine ' // Trim(C)
    WRITE(10, '(A)') 'Program stoped in subroutine ' // Trim(C)
    CLOSE(10)
    CALL DelayedExit(2.0)
  END IF
!$  CALL MPI_INITIALIZED(ExistQ, I)
!$  IF(ExistQ) CALL MPI_FINALIZE(I)
  call clean_mem
  call rexit("Related encountered a fatal error.")
END SUBROUTINE StopOnError
!
Subroutine DelayedExit(T)
  IMPLICIT NONE
  REAL, INTENT(IN) :: T
  REAL :: DelayTime(2)
  Call CPU_TIME(DelayTime(1))
  DO
    Call CPU_TIME(DelayTime(2))
    IF(DelayTime(2) - DelayTime(1) > T) EXIT
  END DO
END Subroutine DelayedExit
!
Subroutine WriteMidR(ThreadString, N, ID)
  USE VARIABLES
  IMPLICIT NONE
  Character(Len = 6), INTENT(IN) :: ThreadString
  Character(Len = 1) :: C
  Character(Len = 6) :: Pos, C6
  Character(Len = 999) :: C999
  Character(Len = 12) :: C12
  INTEGER(I4B), INTENT(IN) :: N, ID(2, N)
  INTEGER(I4B) :: I, J, L
  Pos = Merge('Delete', 'keep  ', N == 0)
  OPEN(10, File = 'Relatedness_' // ThreadString, Position = 'APPEND')
  DO J = 1, N
    C999 = ''
    DO L = 1, 2
      WRITE(C6, '(I6)') ID(L, J)
      C999 = Trim(C999) // Trim(AdjustL(C6)) // ','
    END DO
    DO L = 1, 7
      WRITE(C12, '(F12.4)') REst(L, J)
      C999 = Trim(C999) // Trim(AdjustL(C12)) // ','
    END DO
    WRITE(10,'(A)') Trim(C999)
  END DO
  CLOSE(10, Status = Pos)
  DO I = 1, 2
    WRITE(C, '(I1)') I
    OPEN(10, File = 'Delta_' // Trim(ThreadString) // '_' // C, Position = 'APPEND')
    DO J = 1, N
      C999 = ''
      DO L = 1, 2
        WRITE(C6, '(I6)') ID(L, J)
        C999 = Trim(C999) // Trim(AdjustL(C6)) // ','
      END DO
      DO L = 1, 4
        WRITE(C12, '(F12.4)') Delta(I, J, L)
        C999 = Trim(C999) // Trim(AdjustL(C12)) // ','
      END DO
      WRITE(10,'(A)') Trim(C999)
    END DO
    CLOSE(10, Status = Pos)
  END DO
END Subroutine WriteMidR
!
Subroutine WriteMidCI(ThreadString, N, ID)
  USE VARIABLES
  IMPLICIT NONE
  Character(Len = 6), INTENT(IN) :: ThreadString
  Character(Len = 1) :: C
  Character(Len = 6) :: C6, Pos
  Character(Len = 12) :: C12
  Character(Len = 999) :: C999
  INTEGER(I4B), INTENT(IN) :: N, ID(2, N)
  INTEGER(I4B) :: I, J, K, L
  Pos = Merge('Delete', 'keep  ', N == 0)
  OPEN(10, File = 'RCI95_' // Trim(ThreadString), Position = 'APPEND')
  DO J = 1, N
    C999 = ''
    DO L = 1, 2
      WRITE(C6, '(I6)') ID(L, J)
      C999 = Trim(C999) // Trim(AdjustL(C6)) // ','
    END DO
    DO L = 1, 7
      DO I = 1, 2
        WRITE(C12, '(F12.4)') RCI95(I, L, J)
        C999 = Trim(C999) // Trim(AdjustL(C12)) // ','
      END DO
    END DO
    WRITE(10,'(A)') Trim(C999)
  END DO
  CLOSE(10, Status = Pos)
  DO I = 1, 2
    WRITE(C, '(I1)') I
    OPEN(10, File = 'DCI95_' // Trim(ThreadString) // '_' // C, Position = 'APPEND')
    DO J = 1, N
      C999 = ''
      DO L = 1, 2
        WRITE(C6, '(I6)') ID(L, J)
        C999 = Trim(C999) // Trim(AdjustL(C6)) // ','
      END DO
      DO L = 1, 4
        DO K = 1, 2
          WRITE(C12, '(F12.4)') DCI95(K, I, L, J)
          C999 = Trim(C999) // Trim(AdjustL(C12)) // ','
        END DO
      END DO
      WRITE(10,'(A)') Trim(C999)
    END DO
    CLOSE(10, Status = Pos)
  END DO
END Subroutine WriteMidCI
!
SUBROUTINE READ_DATA
  USE Variables
  IMPLICIT NONE
  INTEGER :: I, J, K, IOS
  INTEGER, Allocatable :: GType(:, :), Idx(:)
  CHARACTER(Len = 22) :: ErrLn
  CHARACTER(Len = 8) :: C8
  REAL :: S
  EXTERNAL :: ConvertGenotype, StopOnDataError
  CALL DATE_AND_TIME(DATE,TIME)
  ErrorMessage = ''
  ErrLn = 'Error in reading line '
  ErrorMessage(2) = 'TrioR.DAT'
  OPEN(UNIT=10, FILE='TrioR.DAT', IOSTAT = IOS, ERR = 1, STATUS = 'OLD')
  ErrorMessage(2) = ''
  ErrorMessage(1) = ErrLn // '1, Dataset path & name'
  READ(10, *, IOSTAT = IOS, ERR = 9, End = 9) FName
  FName = Trim(FName)
  ErrorMessage(1) = ErrLn // '2, # loci'
  READ(10, *, IOSTAT = IOS, ERR = 9, End = 9) NumLoci
  IF(NumLoci < 1) then
    if (allocated(GType)) deallocate(GType)
    if (allocated(Idx))   deallocate(Idx)
    CALL StopOnDataError('The # loci should be greater than 0!')
  end if
  ALLOCATE(ObsNumAllele(NumLoci), err(NumLoci))
  ErrorMessage(1) = ErrLn // '3, # alleles per locus'
  READ(10, *, IOSTAT = IOS, ERR = 9, End = 9)ObsNumAllele(1 : NumLoci)
  IF(ANY(ObsNumAllele(1 : NumLoci) < 1) .OR. ANY(ObsNumAllele(1 : NumLoci) > 127)) then
    if (allocated(GType)) deallocate(GType)
    if (allocated(Idx))   deallocate(Idx)
    CALL StopOnDataError('The observed # alleles at each locus should be in range [1,127]!')
  end if
  MaxNumAllele = Maxval(ObsNumAllele(1 : NumLoci)) * Merge(2, 1, KnownAlleleFre == 1)
  ALLOCATE(ObsAlleleFre(MaxNumAllele, NumLoci), Idx(MaxNumAllele))
  ErrorMessage(1) = ErrLn // '4, allele frequency indicator'
  READ(10, *, IOSTAT = IOS, ERR = 9, End = 9)KnownAlleleFre  !1/0=indicating known/unknown allele frequencies
  IF(KnownAlleleFre /= 0 .AND. KnownAlleleFre /= 1) then
    if (allocated(GType)) deallocate(GType)
    if (allocated(Idx))   deallocate(Idx)
    CALL StopOnDataError('Allele frequency indicator should be 0 or 1!')
  end if
  IF(KnownAlleleFre == 1)THEN
    ALLOCATE(ALLELE(NumLoci, MaxNumAllele))
    ErrorMessage(2) = 'FrequencyData.Txt'
    Open(11, FILE = FName(1 : LEN_TRIM(FName) - 10) // 'FrequencyData.Txt', &
      IOSTAT = IOS, ERR = 1, STATUS = 'OLD')
    ErrorMessage(2) = ''
    DO I = 1, NumLoci
      WRITE(C8, '(I8)') I
      ErrorMessage(1) = ErrLn // 'of alleles at locus' // Trim(AdjustL(C8))
      READ(11, *, IOSTAT = IOS, ERR = 9, End = 9) ALLELE(I, 1 : ObsNumAllele(I))
      ErrorMessage(1) = ErrLn // 'of allele frequencies at locus' // Trim(AdjustL(C8))
      READ(11, *, IOSTAT = IOS, ERR = 9, End = 9) ObsAlleleFre(1 : ObsNumAllele(I), I)
      IF(ANY(ObsAlleleFre(1 : ObsNumAllele(I), I) < 0.0) .OR. &
        ANY(ObsAlleleFre(1 : ObsNumAllele(I), I) > 1.0)) then
          if (allocated(GType)) deallocate(GType)
          if (allocated(Idx))   deallocate(Idx)
          CALL StopOnDataError( &
            'Allele frequency should be in the range [0, 1] at locus ' &
            // Trim(AdjustL(C8)))
        end if
      S = Sum(ObsAlleleFre(1 : ObsNumAllele(I), I))
      IF(S < 0.8 .OR. S > 1.2) then
        if (allocated(GType)) deallocate(GType)
        if (allocated(Idx))   deallocate(Idx)
        CALL StopOnDataError( &
          'The sum of allele frequency is far from 1 at locus ' // Trim(AdjustL(C8)))
      end if
    END DO
    CLOSE(11)
  END IF
  ErrorMessage(1) = ErrLn // 'of inbreeding indicator'
  READ(10, *, IOSTAT = IOS, ERR = 9, End = 9)Constrain       !1/0 Constrain/NOT constrain F=0
  IF(Constrain /= 0 .AND. Constrain /= 1) then
    if (allocated(GType)) deallocate(GType)
    if (allocated(Idx))   deallocate(Idx)
    CALL StopOnDataError( &
      'Indicator for inbreeding should be 1 or 0')
  end if
  ErrorMessage(1) = ErrLn // 'of random number seed'
  READ(10, *, IOSTAT = IOS, ERR = 9, End = 9)I_SEED          !Seed of random number in estimation
  ErrorMessage(1) = ErrLn // 'of # bootstrapping samples'
  READ(10, *, IOSTAT = IOS, ERR = 9, End = 9)BootsSamples    !>0/<=0 indicating estimate/notEstimate CI95
  IF(BootsSamples > 100000) then
    if (allocated(GType)) deallocate(GType)
    if (allocated(Idx))   deallocate(Idx)
    CALL StopOnDataError('#bootstrapping samples too big!')
  end if
  EstimateCI95 = Merge(1, 0, BootsSamples > 0)
  ErrorMessage(1) = ErrLn // 'of #reference individuals'
  READ(10, *, IOSTAT = IOS, ERR = 9, End = 9)NumTrios   !Number of reference individuals in Triadic ML method
  IF(NumTrios > 10000) then
    if (allocated(GType)) deallocate(GType)
    if (allocated(Idx))   deallocate(Idx)
    CALL StopOnDataError('#reference individuals too big!')
  end if
  ErrorMessage(1) = ErrLn // 'of sample size'
  READ(10, *, IOSTAT = IOS, ERR = 9, End = 9)SmplSize
  IF(SmplSize < 1) then
    if (allocated(GType)) deallocate(GType)
    if (allocated(Idx))   deallocate(Idx)
    CALL StopOnDataError('# individuals should be > 0!')
  end if
  ALLOCATE(IndivGType(2, NumLoci, SmplSize), IndivID(SmplSize))
  ALLOCATE(RawAlleleNum(NumLoci), RawAlleleID(NumLoci), GType(2, NumLoci))
  DO I = 1, NumLoci
    Allocate(RawAlleleID(I)%R(2))
    RawAlleleNum(I) = 0_I1B
  END DO
  ErrorMessage(2) = 'GenotypeData.Txt'
  Open(11, FILE = FName(1 : LEN_TRIM(FName) - 10) // 'GenotypeData.Txt', &
      IOSTAT = IOS, ERR = 1, STATUS = 'OLD')
  ErrorMessage(2) = ''
  DO I = 1, SmplSize
    WRITE(C8, '(I8)') I
    ErrorMessage(1) = ErrLn // 'of the ID & genotypes of individual ' // Trim(AdjustL(C8))
    READ(11, *, IOSTAT = IOS, ERR = 9, End = 9)IndivID(I),(GType(1:2, J), J = 1, NumLoci)
    CALL ConvertGenotype(IndivGType(:, :, I), GType)
  END DO
  CLOSE(11)
  IF(KnownAlleleFre == 1)THEN
    DO I = 1, NumLoci
      FORALL(K = 1 : RawAlleleNum(I)) Idx(K) = RawAlleleID(I)%R(K)
      DO K = 1, RawAlleleNum(I)
        IF(ANY(RawAlleleID(I)%R(K) == ALLELE(I, 1 : ObsNumAllele(I)))) CYCLE
        ObsNumAllele(I) = ObsNumAllele(I) + 1
        IF(ObsNumAllele(I) > UBound(ALLELE, 2)) THEN
          if (allocated(GType)) deallocate(GType)
          if (allocated(Idx))   deallocate(Idx)
          WRITE(*,'(A, A, I6, A/A)') &
            'Too many alleles which are observed in the sample', &
            ' but have unknown frequencies at locus ', I, ' !', 'Program stops'
          call StopOnDataError('Too many alleles observed + unknown frequencies!')
        END IF
        ALLELE(I, ObsNumAllele(I)) = RawAlleleID(I)%R(K)
        ObsAlleleFre(ObsNumAllele(I), I) = 0.0001
      END DO
      IF(ObsNumAllele(I) > UBound(RawAlleleID(I)%R, 1)) THEN
        Deallocate(RawAlleleID(I)%R)
        Allocate(RawAlleleID(I)%R(ObsNumAllele(I)))
        RawAlleleID(I)%R(1 : RawAlleleNum(I)) = Idx(1 : RawAlleleNum(I))
      END IF
      IF(ObsNumAllele(I) > RawAlleleNum(I)) THEN
        J = RawAlleleNum(I)
        DO K = 1, ObsNumAllele(I)
          IF(ANY(ALLELE(I, K) == Idx(1 : RawAlleleNum(I)))) CYCLE
          J = J + 1
          Idx(J) = K
          RawAlleleID(I)%R(J) = ALLELE(I, K)
        END DO
      END IF
      DO J = 1, RawAlleleNum(I)
        Idx(J) = 0
        DO K = 1, ObsNumAllele(I)
          IF(ALLELE(I, K) /= RawAlleleID(I)%R(J)) CYCLE
          Idx(J) = K
          EXIT
        END DO
      END DO
      FORALL(J = 1 : ObsNumAllele(I), Idx(J) /= 0)
        ALLELE(I, J) = J
        ObsAlleleFre(J, I) = ObsAlleleFre(Idx(J), I)
      END FORALL
    END DO
  END IF
  ErrorMessage(1) = ErrLn // 'of mistyping rate at each locus'
  READ(10, *, IOSTAT = IOS, ERR = 9, End = 9) Err(1 : NumLoci)    !Error rate/locus used in estimation
  ErrorMessage(1) = ErrLn // 'of indicator for dyads selector'
  READ(10, *, IOSTAT = IOS, ERR = 9, End = 9) AllDyads            !1/0=All dyads=Y/N
  IF(AllDyads /= 0 .AND. AllDyads /= 1) then
    if (allocated(GType)) deallocate(GType)
    if (allocated(Idx))   deallocate(Idx)
    CALL StopOnDataError('Dyad selection indicator should be 1 or 0!')
  end if
  READ(10,  *, IOSTAT = IOS, ERR = 9, End = 9) RMethod(1 : 7)    !2/1/0=PointEst+CI95/PointEst/NoEst for methods 1-7
  DO J = 1, 7
    IF(RMethod(J) /= 0 .AND. RMethod(J) /= 1 .AND. RMethod(J) /= 2) then
      if (allocated(GType)) deallocate(GType)
      if (allocated(Idx))   deallocate(Idx)
      CALL StopOnDataError('Estimator selection indicator should be 2, 1 or 0!')
    end if
  END DO
  CLOSE(10)
  IF(EstimateCI95 == 1) THEN
    IF(ALL(RMethod(1 : 7) < 2)) THEN
      EstimateCI95 = 0
      BootsSamples = 0
    END IF
  END IF
  ErrorMessage(1) = ''
  RETURN
1 ErrorMessage(1) = 'Input file ' // Trim(AdjustL(ErrorMessage(2))) // ' does not exist!'
  ErrorMessage(2) = ''
  GOTO 7
9 continue
  if (allocated(GType)) deallocate(GType)
  if (allocated(Idx))   deallocate(Idx)
  CALL StopOnDataError(ErrorMessage(1))
7 continue
  if (allocated(GType)) deallocate(GType)
  if (allocated(Idx))   deallocate(Idx)
  Call StopOnError(1, 'ReadData')
END SUBROUTINE READ_DATA
!
Subroutine StopOnDataError(C)
  Use Variables
  IMPLICIT NONE
  Character(*), INTENT(IN) :: C
  EXTERNAL :: StopOnError
  ErrorMessage(1) = Trim(AdjustL(C))
  ErrorMessage(2) = 'Errors in DATA. Insufficient data or incorrect format.'
  ErrorMessage(3) = 'Please check DATA and format and then re-run the program'
  Call StopOnError(1, 'StopOnDataError')
END Subroutine StopOnDataError
!
Subroutine ConvertGenotype(GType1, GType2)
  Use Variables
  IMPLICIT NONE
  Integer(I4B), INTENT(IN) :: GType2(2, NumLoci)
  Integer(I1B), INTENT(OUT) :: GType1(2, NumLoci)
  Integer(I4B) :: I, J, K, L, M, N, IA(127)
  DO I = 1, NumLoci
    DO J = 1, 2
      IF(GType2(J, I) <= 0) THEN
        GType1(:, I) = 0_I1B
        EXIT
      ELSE
        L = GType2(J, I)
        DO K = 1, RawAlleleNum(I)
          IF(L == RawAlleleID(I)%R(K)) THEN
            GType1(J, I) = INT(K, I1B)
            L = 0
            EXIT
          END IF
        END DO
        IF(L == 0) CYCLE
        IF(RawAlleleNum(I) >= 127_I1B) THEN
          WRITE(*, '(A, I9, A)') &
            'The observed # alleles at locus ', I, ' > 127!'
          call StopOnDataError('More than 127 alleles at a locus!')
        END IF
        RawAlleleNum(I) = RawAlleleNum(I) + 1_I1B
        GType1(J, I) = RawAlleleNum(I)
        N = Size(RawAlleleID(I)%R, 1)
        IF(N < RawAlleleNum(I)) THEN
          FORALL(M = 1 : N) IA(M) = RawAlleleID(I)%R(M)
          Deallocate(RawAlleleID(I)%R)
          Allocate(RawAlleleID(I)%R(N + 10))
          FORALL(M = 1 : N) RawAlleleID(I)%R(M) = IA(M)
        END IF
        RawAlleleID(I)%R(RawAlleleNum(I)) = L
      END IF
    END DO
  END DO
End Subroutine ConvertGenotype
!
SUBROUTINE INITIATE
  USE variables
  USE IBDVariables
  IMPLICIT NONE
  INTEGER I,J
  MaxNumAllele = MAXVAL(ObsNumAllele)
  ALLOCATE(RawAllele(NumLoci, -1  :MaxNumAllele), &
    REst(7, 99), Delta(2, 99, 4), G(6, NumLoci), A(6, NumLoci), &
    LocusSelect1(NumLoci), LocusSelect2(NumLoci), &
    Prob_Locus2(NumLoci), LocusRepeat1(NumLoci), LocusRepeat2(NumLoci))
  IF(Constrain /= 1) THEN
    I = 1
    J = 4
  ELSE
    I = 2
    J = 3
  END IF
  Allocate(FEst(I : J, SmplSize))
  IF(RMethod(1) > 0) THEN
    I = 66
  ELSE IF(RMethod(7) > 0) THEN
    I = 9
  ELSE
    I = 0
  END IF
  IF(I > 0) ALLOCATE(CoeXij2(NumLoci), Coe(NumLoci, I), GConfig(NumLoci), &
    IBDConfig2(NumLoci), Prob_Locus(NumLoci))
  IF(RMethod(1) > 0) ALLOCATE(GUnrelated(2, NumLoci, NumTrios))
  FORALL(J = 1 : NumLoci) LocusSelect1(J) = J
  LocusRepeat1(:) = 1_I2B
  AccountError = (Constrain == 1 .AND. ANY(Err(1 : NumLoci) > 0.0001_DP))
  IF(EstimateCI95 == 1) ALLOCATE(RCI95(2,7,99),DCI95(2,2,4,99),FCI95(2,4,SmplSize))
  LogMul=Log(DBLE(MUL))
  IDUM = -ABS(I_SEED)
  I025(1)=Max(1,NINT(BootsSamples*.025))
  I025(2)=NINT(BootsSamples*.975)
  IF(Constrain == 0)THEN
    NumMarkVar = 66
    FORALL(J = 1 : 66) MarkVar(J) = J
    NumMarkVar2 = 9
    FORALL(J = 1 : 9) MarkVar2(J) = J
  ELSE
    NumMarkVar = 16
    MarkVar(1 : 16) = (/14,19,38,39,40,44,45,46,50,51,52,59,63,64,65,66/)
    NumMarkVar2 = 3
    MarkVar2(1 : 3) = (/7, 8, 9/)
  END IF
END SUBROUTINE INITIATE
!
SUBROUTINE SimuGtype
  use variables
  IMPLICIT NONE
  INTEGER I,J,K,M
  REAL RAN1, RNumber
  REAL(DP) :: CumAlleleFre(NumLoci,MaxNumAllele)
  EXTERNAL :: RAN1
  DO J=1,NumLoci
    CumAlleleFre(J,1)=ObsAlleleFre(1, J)
    DO K=2,ObsNumAllele(J)-1
      CumAlleleFre(J,K)=CumAlleleFre(J,K-1)+ObsAlleleFre(K, J)
    END DO
    CumAlleleFre(J,ObsNumAllele(J))=1.001
  END DO
  DO I=1,NumTrios
    DO J=1,NumLoci
      DO M=1,2
        RNumber = RAN1()
        DO K=1,ObsNumAllele(J)
          IF(RNumber<CumAlleleFre(J,K))THEN
            GUnrelated(M, J, I) = Int(K, I1B)
            EXIT
          END IF
        END DO
      END DO
      IF(GUnrelated(1,J,I)>GUnrelated(2,J,I)) &
        GUnrelated(1:2,J,I)=GUnrelated(2:1:-1,J,I)
    END DO
  END DO
END SUBROUTINE SimuGtype
!
SUBROUTINE DataProcess
  use variables
  use WangVariables
  IMPLICIT NONE
  INTEGER:: LAllele(MaxNumAllele),Tem_Array(1)
  INTEGER:: I,J,K,L,M
  REAL(DP):: X,Y(MaxNumAllele), AA(2 : 4)
  LOGICAL:: FLAG(SmplSize,2)
  DO I=1,NumLoci
    IF(KnownAlleleFre==0)THEN
      ObsNumAllele(I)=0
      ObsAlleleFre(:, I)=0._DP
      DO K=1,SmplSize
        IF(IndivGType(1, I, K)>0)THEN
          DO J=1,2
            M=0
            DO L=1,ObsNumAllele(I)
              IF(IndivGType(J, I, K)==LAllele(L))THEN
                M=L
                EXIT
              END IF
            END DO
            IF(M==0)THEN
              ObsNumAllele(I)=ObsNumAllele(I)+1
              LAllele(ObsNumAllele(I))=IndivGType(J, I, K)
              ObsAlleleFre(ObsNumAllele(I), I)=1._DP
            ELSE
              ObsAlleleFre(M, I)=ObsAlleleFre(M, I)+1._DP
            END IF
          END DO
        END IF
      END DO
      X=1._DP/Sum(ObsAlleleFre(1:ObsNumAllele(I), I))
      DO M=1,ObsNumAllele(I)
        ObsAlleleFre(M, I) = ObsAlleleFre(M, I) * X
      END DO
    ELSE
      FORALL(K=1 : ObsNumAllele(I)) LAllele(K)=Allele(I,K)
    END IF
    FORALL(K=1 : ObsNumAllele(I)) Y(K)=ObsAlleleFre(K, I)
    DO J=1,ObsNumAllele(I)
      Tem_Array=MAXLOC(Y(1:ObsNumAllele(I)))
      L=LAllele(Tem_Array(1))
      RawAllele(I,J)=L
      ObsAlleleFre(J, I)=Y(Tem_Array(1))
      Y(Tem_Array(1))=-2._DP
    END DO
    FLAG=.TRUE.
    DO J=1,ObsNumAllele(I)
      DO K=1,SmplSize
        DO L=1,2
          IF(FLAG(K,L).AND.IndivGType(L, I, K)==RawAllele(I,J))THEN
            IndivGType(L, I, K) = INT(J, I1B)
            FLAG(K,L)=.FALSE.
          END IF
        END DO
      END DO
    END DO
    RawAllele(I,-1)=-1
    DO K=1,SmplSize
      IF(IndivGType(1,I,K) > IndivGType(2,I,K)) &
         IndivGType(1:2,I,K)=IndivGType(2:1:-1,I,K)
    END DO
  END DO
  Allocate(RecAlleleFre(MaxNumAllele, NumLoci))
  DO I=1,NumLoci
    DO J=1,ObsNumAllele(I)
      IF(ObsAlleleFre(J, I) < 0.0001) ObsAlleleFre(J, I) = 0.0001
      RecAlleleFre(J, I) = 1._DP / ObsAlleleFre(J, I)
    END DO
    IF(ObsNumAllele(I) < 2) FORALL(K = 1 : SmplSize) &
      IndivGType(1 : 2, I, K) = 0_I1B
  END DO
! Initialise Wang variables
  Allocate(Fac(6, NumLoci), WEIGHT(NumLoci))
  DO L = 1, NumLoci
    AA(2:4)=0._DP
    DO I=1,ObsNumAllele(L)
      DO K=2,4
        AA(K)=AA(K)+ObsAlleleFre(I, L)**K
      END DO
    END DO
    weight(L) = 1._DP / (2._DP*AA(2)-AA(3))
    IF(WEIGHT(L) <= Small) CYCLE
    Fac(1,L)=(4*(AA(2)-AA(2)**2-2*AA(3)+2*AA(4)))*WEIGHT(L)
    Fac(2,L)=(1-7*AA(2)+4*AA(2)**2+10*AA(3)-8*AA(4))*WEIGHT(L)
    Fac(3,L)=(4*(AA(3)-AA(4)))*WEIGHT(L)
    Fac(4,L)=(2*(AA(2)-3*AA(3)+2*AA(4)))*WEIGHT(L)
    Fac(5,L)=(2*AA(2)**2-AA(4))*WEIGHT(L)
    Fac(6,L)=(AA(2)-2*AA(2)**2+AA(4))*WEIGHT(L)
  END DO
END SUBROUTINE DataProcess
!
SUBROUTINE OUTPUT
  USE variables
  IMPLICIT NONE
  INTEGER I,J,K,L,M,JJ, KK
  Character(Len=1000) TempFileName, Line1, Line2
  Character(Len=1) :: C1
  Character(Len = 6) :: ThreadString
  Character(Len = 5) :: C5
  EXTERNAL :: WriteLine
  WRITE(*,'(A)')'Saving results to files...'
  OPEN(10,File=FNAME)
  WRITE(10,'(A)')'(1)Parameters in simulations'
  WRITE(10,'(/A)')'  SmplSize   NumLoci Constrain KnownAlleleFre'
  WRITE(10,'(4I10)')SmplSize,NumLoci,Constrain,KnownAlleleFre
  WRITE(10,'(/A)')' Number of observed alleles at each locus'
  WRITE(10,'(20I3)')ObsNumAllele(1:NumLoci)
  WRITE(10,'(/A)')' Allele frequency'
  DO I = 1, NumLoci
    WRITE(10,'(/A,I5)')' Locus:',I
    WRITE(10,'(200I8)')(RawAlleleID(I)%R(RawAllele(I,J)), J = 1, ObsNumAllele(I))
    WRITE(10,'(200F8.4)')ObsAlleleFre(1:ObsNumAllele(I), I)
  END DO
  IF(AccountError) THEN
    WRITE(10,'(/A)')' Genotyping error/mutation rate at each locus'
    WRITE(10,'(10(A,I5))')(' Locus', I, I = 1, NumLoci)
    WRITE(10,'(10(F11.4))')(Err(I), I = 1, NumLoci)
  END IF
  WRITE(10,'(//2A)')'(2)Point estimate of R for each dyad'
  WRITE(10,'(/2A)')'Pair,Ind1,Ind2,TrioEst,WEst,LLEst,LREst,REst,QGEst,MEst'
  KK=0
  J=Len_Trim(FName)
  TempFileName=Trim(FName(1:J-10)) // "RelatednessEstimates.Txt"
  OPEN(11,File=Trim(TempFileName))
  DO M = 0, NumThreads - 1
    WRITE(ThreadString, '(I6)') M
    ThreadString = AdjustL(ThreadString)
    OPEN(12, File = 'Relatedness_' // ThreadString)
    DO
      READ(12, *, IOSTAT = K) I, J, REst(1 : 7, 1)
      IF(K /= 0) EXIT
      KK = KK + 1
      CALL WriteLine(C5, Line1, Line2, KK, I, J, 7, REst(1 : 7, 1), 1, 1, Delta(1,1,1))
      WRITE(10,'(A)')Line1(1:I) // Line2(1 : J)
      WRITE(11,'(A)')Line1(1:I) // C5 // Line2(1 : J)
    END DO
    CLOSE(12, STATUS = 'DELETE')
  END DO
  Close(11)
  DO L=1,2
    J=Len_Trim(FName)
    WRITE(C1, '(I1)')Merge(8, 7, L == 1)
    WRITE(10,'(//2A)')'  Point estimate of Delta' // C1 // ' for each dyad'
    TempFileName=Trim(FName(1:J-10)) // 'Delta' // C1 // 'Estimates.Txt'
    OPEN(11,File=Trim(TempFileName))
    WRITE(10,'(/2A)')'Pair,Ind1,Ind2,TrioEst,WEst,LREst,MEst'
    WRITE(C1, '(I1)') L
    KK=0
    DO M = 0, NumThreads - 1
      WRITE(ThreadString, '(I6)') M
      ThreadString = AdjustL(ThreadString)
      OPEN(12, File = 'Delta_' // Trim(ThreadString) // '_' // C1)
      DO
        READ(12, *, IOSTAT = K) I, J, Delta(1, 1, 1:4)
        IF(K /= 0) EXIT
        KK = KK + 1
        CALL WriteLine(C5, Line1, Line2, KK, I, J, 4, Delta(1,1,1:4), 1, 1, Delta(1,1,1))
        WRITE(10,'(A)')Line1(1:I) // Line2(1 : J)
        WRITE(11,'(A)')Line1(1:I) // C5 // Line2(1 : J)
      END DO
      CLOSE(12)
    END DO
    Close(11)
  END DO
  J=Len_Trim(FName)
  TempFileName=Trim(FName(1:J-10)) // "Delta7&8Estimates.Txt"
  OPEN(11,File=Trim(TempFileName))
  KK=0
  DO M = 0, NumThreads - 1
    WRITE(ThreadString, '(I6)') M
    ThreadString = AdjustL(ThreadString)
    OPEN(12, File = 'Delta_' // Trim(ThreadString) // '_1')
    OPEN(13, File = 'Delta_' // Trim(ThreadString) // '_2')
    DO
      READ(12, *, IOSTAT = K) I, J, Delta(1, 1, 1:4)
      IF(K /= 0) EXIT
      READ(13, *, IOSTAT = K) I, J, Delta(2, 1, 1:4)
      KK=KK+1
      CALL WriteLine(C5, Line1, Line2, KK, I, J, 1, &
        Delta(1,1,1), 2, 4, Delta(2 : 1 : -1, 1, 1 : 4))
      WRITE(11,'(A)')Line1(1:I) // C5 // Line2(1 : J)
    END DO
    CLOSE(12, STATUS = 'DELETE')
    CLOSE(13, STATUS = 'DELETE')
  END DO
  Close(11)
  IF(EstimateCI95==1)THEN
    WRITE(10,'(//A,I6/9A)') &
      '(3)R 95%Confidence interval by bootstrapping over loci with #BootsSamples=',&
      BootsSamples,'Pair,Ind1,Ind2,TrioEst,WEst,LLEst,LREst,REst,QGEst,MEst'
    J=Len_Trim(FName)
    TempFileName=Trim(FName(1:J-10)) // "RelatednessCI95.Txt"
    OPEN(11,File=Trim(TempFileName))
    KK=0
    DO M = 0, NumThreads - 1
      WRITE(ThreadString, '(I6)') M
      OPEN(12, File = 'RCI95_' // Trim(AdjustL(ThreadString)))
      DO
        READ(12, *, IOSTAT = K) I, J, (RCI95(1 : 2, L, 1), L = 1, 7)
        IF(K /= 0) EXIT
        KK = KK + 1
        CALL WriteLine(C5, Line1, Line2, KK, I, J, 1, &
          Delta(1,1,1), 2, 7, RCI95(1:2, 1 : 7, 1))
        WRITE(10,'(A)')Line1(1:I) // Line2(1 : J)
        WRITE(11,'(A)')Line1(1:I) // C5 // Line2(1 : J)
      END DO
      CLOSE(12, STATUS = 'DELETE')
    END DO
    Close(11)
    DO L = 1, 2
      J = Len_Trim(FName)
      TempFileName = Trim(FName(1:J-10)) // "Delta" // ACHAR(48+9-L) // "CI95.Txt"
      OPEN(11, File = Trim(TempFileName))
      WRITE(10,'(//A,I1,A)')'Delta',9-L,' 95%Confidence interval by bootstrapping over loci'
      WRITE(10,'(9A)')'Pair,Ind1,Ind2,TrioEst,WEst,LREst,MEst'
      KK=0
      WRITE(C1, '(I1)') L
      DO M = 0, NumThreads - 1
        WRITE(ThreadString, '(I6)') M
        OPEN(12, File = 'DCI95_' // Trim(AdjustL(ThreadString)) // '_' // C1)
        DO
          READ(12, *, IOSTAT = K) I, J, (DCI95(1 : 2, L, JJ, 1), JJ = 1, 4)
          IF(K /= 0) EXIT
          KK = KK + 1
          CALL WriteLine(C5, Line1, Line2, KK, I, J, 1, &
            Delta(1,1,1), 2, 4, DCI95(1 : 2, L, 1 : 4, 1))
          WRITE(10,'(A)')Line1(1:I) // Line2(1 : J)
          WRITE(11,'(A)')Line1(1:I) // C5 // Line2(1 : J)
        END DO
        CLOSE(12, STATUS = 'DELETE')
      END DO
      Close(11)
    END DO
  END IF
  WRITE(10,'(/A/A/A/A/A/A/A/A/A)')'NOTE:','  TrioEst: New Trio estmator',&
    "  WEst   : Wang's Estimator","  LLEst  : Li&Lynch's Estimator",&
    "  LREst  : Lynch&Ritland's Estimator","  REst   : Ritland's Estimator",&
    "  QGEst  : Queller&Goodnight's Estimator", &
    "  MEst   : Milligan's likelihood estimator"
  WRITE(10,'(//A)')'Inbreeding coefficient:'
  J = Len_Trim(FName)
  TempFileName = Trim(FName(1:J-10)) // "InbreedingEstimates.Txt"
  OPEN(11, File = Trim(TempFileName))
  TempFileName = Trim(FName(1:J-10)) // "InbreedingCI95.Txt"
  OPEN(12, File = Trim(TempFileName))
  IF(Constrain == 1)THEN
    WRITE(10,'(A)')'Individual    LH-estmator    LR-estmator'
    DO I = 1, SmplSize
      WRITE(10,'(I10,2F15.4)')I, FEst(2 : 3, I)
      WRITE(11,'(A,2F15.4)')ADJUSTR(IndivID(I)), FEst(2 : 3, I)
    END DO
    IF(EstimateCI95 == 1)THEN
      WRITE(10,'(//A/A)')'CI95 of F from bootstrapping over loci', &
        'Individual         LH-estmator         LR-estmator'
      DO I = 1, SmplSize
        WRITE(10,'(I10,4F10.4)')I, (FCI95(1 : 2, J, I), J = 2, 3)
        WRITE(12,'(I10,4F10.4)')I, (FCI95(1 : 2, J, I), J = 2, 3)
      END DO
    END IF
    WRITE(10,'(/A/A/A/A)')'NOTE:','  LH: Li&Horvits estmator','  LR: Lynch&Ritland estmator'
  ELSE
    WRITE(10,'(A)')'Individual    L3-estmator    LH-estmator    LR-estmator     L2-estmator'
    DO I=1,SmplSize
      WRITE(10,'(I10,4F15.4)')I,FEst(1:4,I)
      WRITE(11,'(A,4F15.4)')ADJUSTR(IndivID(I)),FEst(2:3,I),FEst(1:4:3,I)
    END DO
    IF(EstimateCI95==1)THEN
      WRITE(10,'(//A/3A)')'CI95 of F from bootstrapping over loci', &
        'Individual         L3-estmator         LH-estmator',&
        '         LR-estmator         L2-estmator'
      DO I=1,SmplSize
        WRITE(10,'(I10,8F10.4)')I,(FCI95(1:2,J,I),J=1,4)
        WRITE(12,'(I10,8F10.4)')I,(FCI95(1:2,J,I),J=2,3),(FCI95(1:2,J,I),J=1,4,3)
      END DO
    END IF
    WRITE(10,'(/A/A/A/A/A/A)')'NOTE:','  L3: Trio likelihood estmator',&
      '  LH: Li&Horvits estmator','  LR: Lynch&Ritland estmator',&
      '  L2: Dyad likelihood estmator'
  END IF
  Close(11)
  Close(12)
  WRITE(10,'(//)')
  WRITE(10,*)'Program started:  Date(d/m/y)=',date(7:8),'/',date(5:6),&
    '/',date(1:4),', Time(s/m/hr)=',time(5:6),'/',time(3:4),'/',time(1:2)
  CALL DATE_AND_TIME(DATE,TIME)
  WRITE(10,*)'Program finished: Date(d/m/y)=',date(7:8),'/',date(5:6),&
    '/',date(1:4),', Time(s/m/hr)=',time(5:6),'/',time(3:4),'/',time(1:2)
  WRITE(*,'(2A)')'Output to file: ',FName
  CLOSE(10)
END SUBROUTINE OUTPUT
!
Subroutine WriteLine(C5, Line1, Line2, K, I, J, N, R, N1, N2, R2)
  USE variables
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: I, J
  INTEGER, INTENT(IN) :: K, N, N1, N2
  REAL, INTENT(IN) :: R(N), R2(N1, N2)
  CHARACTER(Len = 5), INTENT(OUT) :: C5
  CHARACTER(Len = 100), INTENT(OUT) :: Line1, Line2
  CHARACTER(Len = 20) :: C20
  INTEGER :: L, M
  WRITE(C20, '(I12)') K
  C5 = IndivID(I)(1:2)//IndivID(J)(1:2) // ','
  Line1 = Trim(AdjustL(C20)) // ',' // Trim(ADJUSTL(IndivID(I))) // &
    ',' // Trim(ADJUSTL(IndivID(J))) // ','
  I = Len_Trim(Line1)
  Line2 =''
  IF(N > 1) THEN
    DO M = 1, N
      WRITE(C20, '(F11.4)') R(M)
      Line2 = Trim(Line2) // Trim(AdjustL(C20)) // Merge(',', ' ', M < N)
    END DO
  ELSE
    DO M = 1, N2
      DO L = 1, N1
        WRITE(C20, '(F11.4)') R2(L, M)
        Line2 = Trim(Line2) // Trim(AdjustL(C20)) // Merge(',', ' ', M < N2 .OR. L < N1)
      END DO
    END DO
  END IF
  J = Len_Trim(Line2)
END Subroutine WriteLine
!
SUBROUTINE WANG(I1,I2,WEst,LLEst,WDelta)
  USE VARIABLES
  USE WangVariables
  IMPLICIT NONE
  INTEGER:: J, K, L, I1, I2, MaxNA
  INTEGER(I1B) :: Flag(NumLoci)
  INTEGER(I1B), Pointer :: G1Ptr(:, :), G2Ptr(:, :)
  REAL(DP):: Est(4), A(3,2), DIVISION, S(3), U, UU, LL(5), X
  REAL:: WEst,LLEst,WDelta(2)
  G1Ptr => IndivGType(:, :, I1)
  G2Ptr => IndivGType(:, :, I2)
  A(1 : 3, 1 : 2) = 0._DP
  U = 0._DP
  MaxNA = 0
  S(1 : 3) = 0._DP
  LL(:) = 0._DP
  Flag(:) = 0_I1B
  DO J = 1, NumLoci
    L = LocusSelectPtr(J)
    IF(Flag(L) == 1_I1B) CYCLE
    Flag(L) = 1_I1B
    IF(G1Ptr(1,L ) < 1_I1B .OR. G2Ptr(1, L) < 1_I1B) CYCLE
    IF(Weight(L) <= Small) CYCLE
    IF(ObsNumAllele(L) > MaxNA) MaxNA = ObsNumAllele(L)
    A(1,1)=A(1,1) + Fac(1, L) * LocusRepeat(L)
    A(1,2)=A(1,2) + Fac(2, L) * LocusRepeat(L)
    A(2,1)=A(2,1) + Fac(3, L) * LocusRepeat(L)
    A(2,2)=A(2,2) + Fac(4, L) * LocusRepeat(L)
    A(3,1)=A(3,1) + Fac(5, L) * LocusRepeat(L)
    A(3,2)=A(3,2) + Fac(6, L) * LocusRepeat(L)
    X = Weight(L) * LocusRepeat(L)
    U = U + X
    LL(1) = LL(1) + LocusRepeat(L)
    LL(2) = LL(2) + LocusRepeat(L) / Weight(L)
    IF(ALL(G1Ptr(1:2,L)==G2Ptr(1:2,L)))THEN
      S(3) = S(3) + X
      LL(3) = LL(3) + LocusRepeat(L)
    ELSE
      K = Count(G1Ptr(1, L) == G2Ptr(1 : 2, L)) + &
          Count(G1Ptr(2, L) == G2Ptr(1 : 2, L))
      IF(K == 2)THEN
        S(2) = S(2) + X
        LL(4) = LL(4) + LocusRepeat(L)
      ELSE IF(K == 1)THEN
        S(1) = S(1) + X
        LL(5) = LL(5) + LocusRepeat(L)
      END IF
    END IF
  END DO
  LL(3) = LL(3) + LL(4) * .75_DP + LL(5) * .5_DP
  UU = 1._DP / U
  A = A * UU
  S(1 : 3) = S(1 : 3) * UU
  IF(MaxNA == 2) A(1, 1 : 2) = 1.E-50_DP
  DIVISION=(A(1,2)**2*A(2,1)*(-1+A(3,1))*(-1+A(2,1)+A(3,1))    &
    -2*A(1,1)*A(1,2)*A(2,1)*(-1+A(3,1))*(A(2,2)+A(3,2))+A(1,1) &
    *(A(2,2)**2*(-1+A(3,1))*(-1+A(1,1)+A(3,1))-2*A(2,1)*A(2,2) &
    *(-1+A(3,1))*A(3,2)+A(2,1)*(A(1,1)+A(2,1))*A(3,2)**2))
  Est(1)=(A(1,2)*A(2,1)*(-1+A(3,1))*(-1+A(2,1)+A(3,1))*S(1)+A(1,1)  &
    **2*(A(2,2)*(-1+A(3,1))*S(2)+A(2,1)*A(3,2)*(-1+S(3)))+A(1,1)*(  &
    -(A(2,1)*A(2,2))+A(2,1)*A(2,2)*A(3,1)-A(2,1)**2*A(3,2)+A(2,1)   &
    *A(2,2)*S(1)-A(2,1)*A(2,2)*A(3,1)*S(1)+A(2,1)*A(3,2)*S(1)-A(2,1)&
    *A(3,1)*A(3,2)*S(1)+A(2,2)*S(2)-2*A(2,2)*A(3,1)*S(2)+A(2,2)     &
    *A(3,1)**2*S(2)+A(2,1)*A(3,2)*S(2)-A(2,1)*A(3,1)*A(3,2)*S(2)    &
    +A(2,1)*A(2,2)*S(3)-A(2,1)*A(2,2)*A(3,1)*S(3)+A(2,1)**2*A(3,2)  &
    *S(3)-A(1,2)*A(2,1)*(-1+A(3,1))*(-1+S(2)+S(3))))/DIVISION
  Est(2)=(A(1,2)**2*A(2,1)*(-1+A(2,1)+A(3,1))*(A(3,1)-S(3))+A(1,1)  &
    *(A(2,1)*A(3,2)**2*(A(1,1)+A(2,1)-S(1)-S(2))+A(2,2)**2*(-1      &
    +A(1,1)+A(3,1))*(A(3,1)-S(3))+A(2,2)*A(3,2)*(A(2,1)-2*A(2,1)    &
    *A(3,1)-A(2,1)*S(1)-S(2)+A(1,1)*S(2)+A(3,1)*S(2)+A(2,1)*S(3)))  &
    +A(1,2)*A(2,1)*((-1+A(2,1)+A(3,1))*A(3,2)*S(1)+A(1,1)*(-2*A(2,2)&
    *A(3,1)+A(3,2)-2*A(3,1)*A(3,2)-A(3,2)*S(2)+2*A(2,2)*S(3)+A(3,2) &
    *S(3))))/DIVISION
  WDelta=REAL(Est(1:2))
  WEst=REAL(0.5_DP*Est(1)+Est(2))
  LLEst=REAL((LL(3)-LL(2))/(LL(1)-LL(2)))
END SUBROUTINE WANG
!
SUBROUTINE LynchRitland(I1,I2,LREst,LRDelta)
  USE VARIABLES
  IMPLICIT NONE
  INTEGER:: I,J,L,I1,I2,KAB,KAC,KAD,KBC,KBD,NDX(2,2),IG(2,2),II(2)
  INTEGER(I1B), Pointer :: G1Ptr(:, :), G2Ptr(:, :)
  INTEGER(I1B) :: IFlag(NumLoci)
  DOUBLE PRECISION P(2),WDD,WRR,DDE,RRE,S(3),X(4), Y, Z
  REAL:: LREst,LRDelta(2)
  DATA NDX/1,2,2,1/
  II(1:2)=(/I1,I2/)
  DO I = 1, 2
    G1Ptr => IndivGType(1:2, :, II(NDX(1,I)))
    G2Ptr => IndivGType(1:2, :, II(NDX(2,I)))
    S(2:3)=0.D0
    WDD=0.D0
    WRR=0.D0
    IFlag(:) = 0_I1B
    DO J = 1, NumLoci
      L = LocusSelectPtr(J)
      IF(IFlag(L) == 1_I1B) CYCLE
      IFlag(L) = 1_I1B
      IF(ObsNumAllele(L)<2) CYCLE
      IG(1,1:2) = G1Ptr(1:2, L)
      IG(2,1:2) = G2Ptr(1:2, L)
      IF(ANY(IG(1 : 2, 1) < 1)) CYCLE
      P(1)=ObsAlleleFre(IG(1,1), L)
      P(2)=ObsAlleleFre(IG(1,2), L)
      IF(ObsNumAllele(L)==2 .AND. P(1)>0.499d0 .AND. P(1)<0.501d0) &
        P(1 : 2)=(/0.501D0, 0.499D0/)
      KBC = Count(IG(1 : 1, 2) == IG(2, 1))
      KBD = Count(IG(1 : 1, 2) == IG(2, 2))
      KAC = Count(IG(1 : 1, 1) == IG(2, 1))
      KAD = Count(IG(1 : 1, 1) == IG(2, 2))
      KAB = Count(IG(1 : 1, 1) == IG(1, 2))
      X(1) = (P(1) + P(1)) * P(2)
      X(2) = X(1) + X(1)
      Y = (1 + KAB) * (P(1) + P(2)) - X(2)
      S(1) = Y / X(1)
      WRR = WRR + S(1) * LocusRepeat(L)
      Z = P(1) * (KBC + KBD) + P(2) * (KAC + KAD)
      S(2) = S(2) + LocusRepeat(L) * S(1) * (Z - X(2)) / Y
      Y = (1 + KAB) * (1.0 - P(1) - P(2)) + X(1)
      S(1)= Y / X(1)
      WDD = WDD + S(1) * LocusRepeat(L)
      S(3) = S(3) + LocusRepeat(L) * S(1) * (X(1)- Z + KAC * KBD + KAD * KBC) / Y
    END DO
    IF(I == 1)THEN
      DDE=S(3)/WDD*0.5d0
      RRE=S(2)/WRR*0.5D0
    ELSE
      DDE=DDE+S(3)/WDD*0.5D0
      RRE=RRE+S(2)/WRR*0.5d0
    END IF
  END DO
  LREst=REAL(RRE)
  LRDelta(2)=REAL(DDE)
  LRDelta(1)=REAL(2*(RRE-DDE))
END SUBROUTINE LynchRitland
!
SUBROUTINE RITLAND(I1,I2,RitEst,QGEst)
  USE VARIABLES
  IMPLICIT NONE
  INTEGER:: I, J, L, I1, I2, LOCI, IG(2,2), NumLoci2, N(2)
  INTEGER(I1B), Pointer :: G1Ptr(:, :), G2Ptr(:, :)
  INTEGER(I1B) :: IFlag(NumLoci)
  REAL(DP),PARAMETER:: RANGE=0.001_DP
  REAL(DP),PARAMETER:: MinP = .5_DP-RANGE, MaxP = .5_DP+RANGE
  REAL(DP):: X, P(MaxNumAllele), Z(2)
  REAL:: RitEst, QGEst
  G1Ptr => IndivGType(1:2, :, I1)
  G2Ptr => IndivGType(1:2, :, I2)
  NumLoci2=0
  RitEst=0._DP
  QGEst=0._DP
  Z(1) = 0._DP
  N(:) = 0
  IFlag(:) = 0_I1B
  DO L = 1, NumLoci
    LOCI = LocusSelectPtr(L)
    IF(IFlag(LOCI) == 1_I1B) CYCLE
    IFlag(LOCI) = 1_I1B
    IF(ObsNumAllele(LOCI) < 2) CYCLE
    IG(1:2, 1) = G1Ptr(1:2, Loci)
    IG(1:2, 2) = G2Ptr(1:2, Loci)
    IF(ANY(IG(1, 1:2) < 1)) CYCLE
    FORALL(J = 1 : ObsNumAllele(Loci)) P(J) = ObsAlleleFre(J, Loci)
    NumLoci2 = NumLoci2 + LocusRepeat(Loci)
    X = 0._DP
    DO I = 1, ObsNumAllele(LOCI)
      IF(ALL(IG(1:2, 1)==I))THEN
        IF(ALL(IG(1:2, 2)==I))THEN
          X = X + 1._DP / P(I)
        ELSE IF(.NOT. ALL(IG(1:2, 2) /= I)) THEN
          X = X + 0.5_DP / P(I)
        END IF
      ELSE IF(.NOT. ALL(IG(1:2, 1) /= I)) THEN
        IF(ALL(IG(1:2, 2)==I))THEN
          X = X + 0.5_DP / P(I)
        ELSE IF(.NOT. ALL(IG(1:2, 2)/=I))THEN
          X = X + 0.25_DP / P(I)
        END IF
      END IF
    END DO
    X=2._DP * (X - 1._DP) / (ObsNumAllele(LOCI)-1._DP) * LocusRepeat(Loci)
    RitEst=RitEst+REAL(X)
    N(1) = N(1) + (Count(IG(1,1)==IG(1 : 2, 2)) + &
                   Count(IG(2,1)==IG(1 : 2, 2))) * LocusRepeat(Loci)
    N(2) = N(2) + (2 + Count(IG(1,1)==IG(2 : 2, 1)) + &
                       Count(IG(1,2)==IG(2 : 2, 2))) * LocusRepeat(Loci)
    Z(1)=Z(1)-(P(IG(1, 1))+P(IG(2,1))+P(IG(1,2))+P(IG(2,2))) * LocusRepeat(Loci)
  END DO
  Z(2) = Z(1) + N(2)
  Z(1) = Z(1) + N(1)
  IF(NumLoci2 > 1) RitEst=RitEst/REAL(NumLoci2)
  IF(ABS(Z(2))>0.01_DP) QGEst=REAL(Z(1)/Z(2))
END SUBROUTINE RITLAND
!
SUBROUTINE InbreedingMT(I, F)
  USE VARIABLES
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: I
  INTEGER :: K, L, M(2)
  INTEGER(I1B), Pointer :: G1(:), G2(:)
  REAL(DP) :: X(2), Y(2), Z
  REAL:: F(2:3)
  G1 => IndivGType(1, :, I)
  G2 => IndivGType(2, :, I)
  F(2 : 3) = 0.0
  DO K = 1, NumLoci
    IF(G1(LocusSelectPtr(K)) < 1_I1B .OR. &
      LocusRepeat(LocusSelectPtr(K)) == 0_I2B) CYCLE
    F(2) = -1.0
    EXIT
  END DO
  IF(F(2) > -0.5) RETURN
  X = 0._DP
  Y = 0._DP
  Z = 0._DP
  M(:) = 0
  IF(CI95Q) THEN
    DO K = 1, NumLoci
      L = LocusSelectPtr(K)
      IF(LocusRepeat(L) == 0_I2B) CYCLE
      IF(G1(L) > 0_I1B) THEN
        IF(G1(L) == G2(L)) THEN
          X(2) = X(2) + RecAlleleFre(G1(L), L) * LocusRepeat(L)
          M(1) = M(1) + LocusRepeat(L)
        ELSE
          Z = Z + (RecAlleleFre(G1(L), L) + RecAlleleFre(G2(L), L)) * &
            LocusRepeat(L)
          M(2) = M(2) + LocusRepeat(L)
        END IF
        X(1) = X(1) + ObsNumAllele(L) * LocusRepeat(L)
      END IF
      LocusRepeat(L) = 0_I2B
    END DO
  ELSE
    DO K = 1, NumLoci
      L = LocusSelectPtr(K)
      IF(G1(L) < 1_I1B) CYCLE
      IF(G1(L) == G2(L)) THEN
        X(2) = X(2) + RecAlleleFre(G1(L), L)
        M(1) = M(1) + 1
      ELSE
        Z = Z + RecAlleleFre(G1(L), L) + RecAlleleFre(G2(L), L)
        M(2) = M(2) + 1
      END IF
      X(1) = X(1) + ObsNumAllele(L)
    END DO
  END IF
  Y(:) = X(2) - M(1)
  X(1) = X(1) - Sum(M(1:2))
  X(2) = X(2) - Sum(M(1:2))
  Y(2) = Y(2) - M(2)
  Y(1) = Y(1) + Z * 0.5_DP - M(2)
  F(2) = REAL(X(2) / X(1))
  F(3) = REAL(Y(2) / Y(1))
END SUBROUTINE InbreedingMT
!
FUNCTION IDyad(I1,I2)
  USE VARIABLES
  IMPLICIT NONE
  INTEGER:: I1,I2,IDyad
  IDyad=I2-(I1*(1 + I1))/2 + (I1-1)*SmplSize
END FUNCTION IDyad
!
Function CalculateRfromIBD(J)
  USE VARIABLES
  IMPLICIT NONE
  INTEGER:: J,KK(0:11,8)
  REAL(DP):: CalculateRfromIBD
  DATA KK(0:,:)/4,1,2,8,20,7*0,7,9,10,12,15,33,36,41,4*0, &
    7,3,5,11,23,26,32,53,4*0,8,16,21,29,37,42,47,54,60,3*0, &
    7,4,6,13,24,27,34,56,4*0,8,17,22,30,35,43,48,55,61,3*0, &
    5,7,14,18,38,46,6*0,11,19,25,28,31,39,40,49,50,51,59,63/
  CalculateRfromIBD=Sum(PCOM(KK(1:KK(0,J),J)))
END Function CalculateRfromIBD
!
SUBROUTINE TrioEst(I1,I2,KK,NumTrios1,NumTrios2,IPointEst,TrioREst, &
                   TrioFEst,TrioDelta)
  USE VARIABLES
  USE IBDVariables
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: I1, I2, KK, NumTrios1, NumTrios2, IPointEst
  REAL, INTENT(OUT):: TrioREst, TrioFEst(2), TrioDelta(2)
  INTEGER:: I3,J,J1,K,L,M,iter, Nest,indx(NumTrios)
  REAL(DP):: Y,X(9),fret,CalculateRfromIBD,U(4)
  REAL:: TLEst(3,NumTrios)
  INTEGER(I2B) :: IBD_Config2
  EXTERNAL Compute_TrioLogL, IBD_Config2, CalculateRfromIBD
  K = KK
  NDIM=66
  ALLOCATE(pcom(0:NDIM))
  IF(Constrain==0)TrioFEst=1.E30
  PCOM=0._DP
  X=0._DP
  TrioREst=1.E30
  DO J = 1, NumLoci
    L = LocusSelectPtr(J)
    G(1 : 2, J) = IndivGType(1 : 2, L, I1)
    G(3 : 4, J) = IndivGType(1 : 2, L, I2)
    IF(ANY(G(1 : 3 : 2, J) < 1_I1B))THEN
      IBDConfig2(J) = 0
    ELSE
      A(1:4,J)=ObsAlleleFre(G(1 : 4, J), L)
      IBDConfig2(J)=IBD_Config2(G(1 : 4, J))
    END IF
  END DO
  Nest=NumTrios
  DO I3=NumTrios1,NumTrios2
    DO J=1,NumLoci
      IF(IBDConfig2(J)==0)CYCLE
      L=LocusSelectPtr(J)
      G(5 : 6, J) = GUnrelated(1 : 2, L, I3)
      A(5 : 6, J) = ObsAlleleFre(G(5 : 6, J), L)
    END DO
    CALL Compute_TrioLogL(iter,fret,K)
    DO J=1,7,2
      X(J)=CalculateRfromIBD(J)
    END DO
    X(8)=CalculateRfromIBD(8)
    TLEst(3,I3)=REAL(X(1)+X(1)+Sum(X(3:7:2))+X(8)*.5_DP)
    TLEst(1:2,I3)=REAL(X(8:7:-1))
    IF(Constrain/=1)THEN
      DO J=2,6,2
        X(J)=CalculateRfromIBD(J)
      END DO
      Y=Sum(X(1:4))
      IF(Y<TrioFest(1))TrioFest(1)=REAL(Y)
      Y=Sum(X(1:2))+Sum(X(5:6))
      IF(Y<TrioFest(2))TrioFest(2)=REAL(Y)
    END IF
    IF(IPointEst==1)THEN
      IF(TLEst(3,I3)<TrioREst) THEN
        TrioREst=TLEst(3,I3)
        IF(EstimateCI95==1)TrioPointer = I3
      END IF
      IF(TLEst(3,I3)<0.0001_DP)THEN
        CALL TestTLest(M,I3)
        IF(M==0)THEN
          Nest=I3
          EXIT
        END IF
      END IF
    ELSE
      TrioREst=TLEST(3,I3)
      TrioDelta(1:2)=TLEst(1:2,I3)
    END IF
  END DO
  IF(IPointEst==1)THEN
    IF(Nest<NumTrios)THEN
      TrioREst=TLEST(3,Nest)
      TrioDelta(1:2)=TLEst(1:2,Nest)
    ELSE
      I3=NINT(.01*Nest)
      CALL INDEXX(Nest,TLEst(3,:),indx)
      U=0._DP
      DO J1=1,NumTrios
        CALL TestTLest(M,indx(J1))
        IF(M==0)THEN
          U(4)=U(4)+1._DP
          DO J=1,3
            U(J)=U(J)+Merge(1._DP/TLEST(J,Indx(J1)), 1.E30_DP, TLEST(J,Indx(J1))>1.E-8_DP)
          END DO
          IF(U(4)>Small.AND.J1>=I3)EXIT
        END IF
      END DO
      IF(U(4)>Small)THEN
        TrioREst=REAL(U(4)/U(3))
        DO J=1,2
          TrioDelta(J)=REAL(U(4)/U(J))
        END DO
      ELSE
        TrioREst=TLEST(3,Indx(1))
        DO J=1,2
          TrioDelta(J)=TLEST(J,Indx(1))
        END DO
      END IF
    END IF
  END IF
  IF(Allocated(pcom)) DEALLOCATE(pcom)
!
CONTAINS
  Subroutine TestTLest(M,IndxIndiv)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: M
    INTEGER:: I,J,L,Iter,IndxIndiv
    REAL(DP):: fret,V
    DO J=1,NumLoci
      IF(IBDConfig2(J)==0)CYCLE
      L=LocusSelectPtr(J)
      G(5 : 6, J) = GUnrelated(1 : 2, L, IndxIndiv)
      A(5 : 6, J) = ObsAlleleFre(G(5 : 6, J), L)
    END DO
    M=0
    DO I=1,10
      CALL Compute_TrioLogL(iter,fret,K)
      DO J=1,7,2
        X(J)=CalculateRfromIBD(J)
      END DO
      V=X(1)+X(1)+Sum(X(3:7:2))+CalculateRfromIBD(8)*.5_DP
      IF(ABS(TLEST(3,IndxIndiv)-V)>0.02)THEN
        M=1
        EXIT
      END IF
    END DO
  END Subroutine TestTLest
END SUBROUTINE TrioEst
!
SUBROUTINE Compute_TrioLogL(iter,fret,KDyad)
  USE VARIABLES
  USE IBDVariables
  IMPLICIT NONE
  REAL(DP),PARAMETER:: FTOL=1.0E-6_DP
  REAL(DP):: func3,fret
  INTEGER(I2B):: IBD_CONFIG3
  INTEGER:: iter,I,J,KDyad
  EXTERNAL:: FUNC3,IBD_CONFIG3,Powell,Coefficient
  DO J = 1, NumLoci
    IF(IBDConfig2(J) == 0)THEN
      GConfig(J) = 0_I2B
    ELSE
      GConfig(J) = IBD_CONFIG3(G(1 : 6, J), IBDConfig2(J))
    END IF
    IF(IBDConfig2(J) > 0 .AND. delta(2,KDyad,4)<0.05 .AND. &
        ObsNumAllele(LocusSelectPtr(J)) > 1)THEN   !Delta(2,:,:)=Delta7
      IF(ANY(G(5, J) == G(1 : 4, J)) .AND. ANY(G(6, J) == G(1 : 4, J)))THEN
        I = G(6, J)
        G(6, J)= -9_I1B
        GConfig(J) = IBD_CONFIG3(G(1 : 6, J), IBDConfig2(J))
        G(6, J) = Int(I, I1B)
      END IF
    END IF
  END DO
  NumMarkVarPtr => NumMarkVar
  MarkVarPtr => MarkVar
  Call Coefficient
  Call Powell(FTOL, iter, fret, FUNC3)
END SUBROUTINE Compute_TrioLogL
!
FUNCTION FUNC3(X)
  USE VARIABLES
  USE IBDVariables
  IMPLICIT NONE
  REAL(DP):: FUNC3, X, Y
  INTEGER:: J, L
  IF(X < 0._DP .OR. X > PCOM(0))THEN
    FUNC3 = 1.E30_DP
    RETURN
  ELSE IF(X < 1.E-6_DP) THEN
    FUNC3 = FUNC0
    RETURN
  END IF
  J = NINT(X * MUL)
  IF(J /= IBDValue)THEN
    IBDValue = J
    FUNC3 = LogLValue
    Y = 1._DP
    DO L = 1, NumLoci1
      Y = Y * (Prob_Locus2(L) + coeXij2(L) * X)
      IF(Y >= 1.E-280_DP) CYCLE
      IF(Y >= Small) THEN
        FUNC3 = FUNC3 + 644.72382603833279152503761_DP
        Y = Y * 1.E280_DP
      ELSE
        FUNC3 = 1.E100_DP
        FUNC3_Old = FUNC3
        RETURN
      END IF
    END DO
    FUNC3 = FUNC3 - Log(Y)
    FUNC3_Old = FUNC3
  ELSE
    FUNC3=FUNC3_Old
  END IF
END FUNCTION FUNC3
!
SUBROUTINE Milligan(I1,I2,MEst,FEst4,MDelta)
  USE VARIABLES
  USE IBDVariables
  IMPLICIT NONE
  REAL(DP),PARAMETER:: FTOL=1.0E-6_DP
  INTEGER:: J,L,I1,I2,Indiv(3),iter,NCount
  INTEGER(I2B):: IBD_Config2
  INTEGER(I1B), Pointer :: G1Ptr(:, :), G2Ptr(:, :)
  REAL(DP):: fret,P(9),BestFunc, FUNC3
  REAL:: MEst,FEst4(2),MDelta(2)
  EXTERNAL IBD_Config2, Coefficient2, FUNC3
  COMMON /IndivPair/Indiv
  G1Ptr => IndivGType(1:2, :, I1)
  G2Ptr => IndivGType(1:2, :, I2)
  NDIM=9
  ALLOCATE(pcom(0:NDIM))
  Indiv(1:2)=(/I1,I2/)
  DO J = 1, NumLoci
    L = LocusSelectPtr(J)
    G(1 : 2, J) = G1Ptr(1 : 2, L)
    G(3 : 4, J) = G2Ptr(1 : 2, L)
    GConfig(J) = Merge(IBD_Config2(G(1 : 4, J)), 0_I2B, ALL(G(1 : 3 : 2, J) > 0))
  END DO
  NumMarkVarPtr => NumMarkVar2
  MarkVarPtr => MarkVar2
  call Coefficient2
  NCount=0
  BestFunc=1.E100_DP
  DO J=1,100
    call Powell(FTOL,iter,fret,FUNC3)
    IF(Fret+1.E-4_DP<BestFunc)THEN
      NCount=0
      BestFunc=Fret
      FORALL(L = 1 : NumMarkVar2) P(MarkVar2(L)) = PCOM(MarkVar2(L))
    END IF
    IF(ABS(BestFunc-Fret)<1.E-2_DP)  NCount=NCount+1
    IF(NCount==3)EXIT
  END DO
  IF(Constrain/=1)THEN
    FEst4(1)=REAL(Sum(P(1:4)))
    FEst4(2)=REAL(Sum(P(1:2))+Sum(P(5:6)))
    MEst=REAL(P(1)+P(1)+Sum(P(3:7:2))+P(8)*.5_DP)
  ELSE
    MEst=REAL(p(7)+p(8)*.5)
  END IF
  MDelta(1:2)=REAL(p(8:7:-1))
  DEALLOCATE(pcom)
END SUBROUTINE Milligan
!
SUBROUTINE Powell(ftol, iter, fret, fun)
  USE Variables
  USE IBDVariables
  IMPLICIT NONE
  INTEGER:: iter, ITMAX, i, j, K, II, JJ, IV, JV, &
    ORDER(NumMarkVarPtr), Jdx(NumMarkVarPtr)
  Logical :: IBD(NumMarkVarPtr)
  REAL(DP):: fret, ftol, fun, tiny1, fp, xmin, brent, X, Y
  REAL(DP):: C(NumLoci)
  PARAMETER (ITMAX = 200, tiny1 = 1.e-25_DP)
  LOGICAL:: LogLInitial
  REAL:: RAN1, RNumber
  EXTERNAL fun, brent, RAN1
  DO I = 1, NumMarkVarPtr
    IBD(I) = ANY(coe(1 : NumLoci, MarkVarPtr(I)) > Small)
  END DO
  Fret = 0._DP
  DO J = 1, NumMarkVarPtr
    IF(IBD(J))THEN
      RNumber = RAN1() + .5
      PCOM(MarkVarPtr(J)) = RNumber
      Fret = Fret + RNumber
    ELSE
      PCOM(MarkVarPtr(J)) = 0._DP
    END IF
  END DO
  Fret = 1._DP/Fret
  DO J = 1, NumMarkVarPtr
    PCOM(MarkVarPtr(J)) = PCOM(MarkVarPtr(J)) * Fret
  END DO
  DO J = 1, NumLoci
    Prob_Locus(J) = 0._DP
    IF(GConfig(J) /= 0_I2B) THEN
      DO I = 1, NumMarkVarPtr
        Prob_Locus(J) = Prob_Locus(J) + coe(J, MarkVarPtr(I)) * PCOM(MarkVarPtr(I))
      END DO
    END IF
  END DO
  LogLInitial = .TRUE.
  Jdx(1 : NumMarkVarPtr) = 0
  K = NumMarkVarPtr - 1
  DO J = 1, NumMarkVarPtr - 1
    DO
      RNumber = RAN1()
      I = NINT(RNumber * K) + 1
      IF(Jdx(I) == 0)EXIT
    END DO
    ORDER(J) = I
    Jdx(I) = J
  END DO
  DO J = 1, NumMarkVarPtr
    IF(Jdx(J) /= 0) CYCLE
    ORDER(NumMarkVarPtr) = J
    EXIT
  END DO
  fret = 1.E100_DP
  DO iter = 1, ITMAX
    fp = fret
    DO II = 1, NumMarkVarPtr
      I = ORDER(II)
      IF(.NOT. IBD(I)) CYCLE
      IV = MarkVarPtr(I)
      DO K = 1, NumLoci
        Prob_Locus(K) = Prob_Locus(K) - coe(K, IV) * PCOM(IV)
      END DO
      DO JJ = II + 1, NumMarkVarPtr
        J = ORDER(JJ)
        IF(.NOT. IBD(J)) CYCLE
        JV = MarkVarPtr(J)
        PCOM(0) = PCOM(IV) + PCOM(JV)
        IF(PCOM(0) < 1.E-4_DP) CYCLE
        IBDValue = -1
        NumLoci1 = 0
        LogLValue = 0._DP
        IF(LogLInitial)THEN
          LogLInitial = .FALSE.
          DO K = 1, NumLoci
            IF(GConfig(K) == 0_I2B) CYCLE
            C(K) = coe(K, JV) * pcom(0)
            Prob_Locus(K) = Prob_Locus(K) - coe(K, JV) * PCOM(JV) + C(K)
            NumLoci1 = NumLoci1 + 1
            CoeXij2(NumLoci1) = coe(K, IV) - coe(K, JV)
            Prob_Locus2(NumLoci1) = Prob_Locus(K)
          END DO
        ELSE
          X = 1._DP
          DO K = 1, NumLoci
            IF(GConfig(K) < 1_I2B) CYCLE
            Y = coe(K, IV) - coe(K, JV)
            C(K) = coe(K, JV) * pcom(0)
            Prob_Locus(K) = Prob_Locus(K) - coe(K, JV) * PCOM(JV) + C(K)
            IF(ABS(Y) > 1.E-20_DP)THEN
              NumLoci1 = NumLoci1 + 1
              CoeXij2(NumLoci1) = Y
              Prob_Locus2(NumLoci1) = Prob_Locus(K)
            ELSE
              X = X * Prob_Locus(K)
              IF(X < Small)THEN
                LogLValue = 1.E100_DP
                X = 1._DP
                EXIT
              ELSE IF(X < 1.E-280_DP)THEN
                LogLValue = LogLValue - Log(X)
                X = 1._DP
              END IF
            END IF
          END DO
          LogLValue = LogLValue - Log(X)
        END IF
        Y = pcom(JV)
        IF(NumLoci1 > 0) THEN
          K = 0
          IF(PCOM(JV) < 1.E-6_DP) THEN
            X = FUN(PCOM(IV) - 1.E-5_DP)
            IF(X > fret) K = 1
          ELSE IF(PCOM(IV) < 1.E-6_DP) THEN
            X = FUN(PCOM(IV) + 1.E-5_DP)
            IF(X > fret) K = 1
          END IF
          IF(K == 0) THEN
            FUNC0 = LogLValue
            X = 1._DP
            DO K = 1, NumLoci1
              X = X * Prob_Locus2(K)
              IF(X >= 1.E-280_DP) CYCLE
              IF(X >= Small) THEN
                FUNC0 = FUNC0 + 644.72382603833279152503761_DP
                X = X * 1.E280_DP
              ELSE
                FUNC0 = 1.E100_DP
                X = -1._DP
                EXIT
              END IF
            END DO
            IF(X > Small) FUNC0 = FUNC0 - Log(X)
            fret = Brent(fun, ftol, xmin, pcom(IV), fret)
            pcom(IV) = xmin
            pcom(JV) = pcom(0) - xmin
          END IF
        END IF
        DO K = 1, NumLoci
          IF(GConfig(K) > 0_I2B) &
          Prob_Locus(K) = Prob_Locus(K) + coe(K, JV) * pcom(JV) - C(K)
        END DO
      END DO
      FORALL(K = 1 : NumLoci, GConfig(K) /= 0_I2B) &
        Prob_Locus(K) = Prob_Locus(K) +  coe(K, IV) * PCOM(IV)
    END DO
    fret = 0._DP
    X = 1._DP
    DO K = 1, NumLoci
      IF(GConfig(K) == 0_I2B) CYCLE
      IF(ABS(Prob_Locus(K)) <= Small)THEN
        fret = 1.E100_DP
        EXIT
      ELSE
        X = X * Prob_Locus(K)
        IF(X > 1.E-280_DP) CYCLE
        fret = fret + 644.72382603833279152503761_DP
        X = X * 1.E280_DP
      END IF
    END DO
    fret = fret - Log(X)
    if(2._DP * (fp - fret) <= ftol * (abs(fp) + abs(fret)) + tiny1) EXIT
  END DO
END SUBROUTINE Powell
!
FUNCTION brent(f,tol,xmin,v, f0)
  USE Variables
  IMPLICIT NONE
  INTEGER ITMAX,iter
  REAL(DP), INTENT(IN) :: f0
  REAL(DP) brent,tol,xmin,f,CGOLD,ZEPS,a,b,d,e,etemp, &
           fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,tiny2
  PARAMETER (ITMAX=100,CGOLD=.3819660_DP,ZEPS=1.0e-10_DP,tiny2=1.E-30_DP)
  EXTERNAL f
  a=0._DP
  b=PCOM(0)
  w=v
  x=v
  e=0._DP
  IF(f0 < 1.E99_DP) THEN
    fx = f0
  ELSE
    fx=f(x)
  END IF
  fv=fx
  fw=fx
  DO iter=1,ITMAX
    xm=0.5_DP*(a+b)
    tol1=tol*abs(x)+ZEPS
    tol2=2._DP*tol1
    if(abs(x-xm)<=(tol2-.5_DP*(b-a))) goto 3
    if(abs(e)>tol1) then
      r=(x-w)*(fx-fv)
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=2._DP*(q-r)
      if(q>0._DP) p=-p
      q=abs(q)
      etemp=e
      e=d
      if(abs(p)>=abs(.5_DP*q*etemp).or.p<=q*(a-x).or.p>=q*(b-x)) goto 1
      d=p/q
      u=x+d
      if(u-a<tol2 .or. b-u<tol2) d=sign(tol1,xm-x)
      goto 2
    endif
 1  if(x>=xm) then
      e=a-x
    else
      e=b-x
    endif
    d=CGOLD*e
 2  if(abs(d)>=tol1) then
      u=x+d
    else
      u=x+sign(tol1,d)
    endif
    fu=f(u)
    if(fu<=fx) then
      if(u>=x) then
        a=x
      else
        b=x
      endif
      v=w
      fv=fw
      w=x
      fw=fx
      x=u
      fx=fu
    else
      if(u<x) then
        a=u
      else
        b=u
      endif
      if(fu<=fw .or. ABS(w-x)<tiny2) then
        v=w
        fv=fw
        w=u
        fw=fu
      else if(fu<=fv .or. ABS(v-x)<tiny2 .or. ABS(v-w)<tiny2) then
        v=u
        fv=fu
      endif
    endif
  END DO
3 xmin=x
  brent=fx
END FUNCTION brent
!
      SUBROUTINE indexx(n,arr,indx)
      USE Variables, ONLY: ErrorMessage
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)THEN
          ErrorMessage(1) = 'Program stops because NSTACK too small in indexx'
          call StopOnError(1, 'indexx')
        end if
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
!
SUBROUTINE sort(n,arr)
      USE Variables, ONLY: ErrorMessage
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=l-1
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)THEN
          ErrorMessage(1) = 'Program stops because NSTACK too small in sort'
          call StopOnError(1, 'sort')
        end if
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
END SUBROUTINE sort
!
FUNCTION IBD_Config2(G)
  USE Variables, ONLY : I1B, I2B
  IMPLICIT NONE
  INTEGER(I2B) :: IBD_Config2
  INTEGER(I1B) :: G(4)
  IF(G(1)==G(2))THEN
    IF(G(3)==G(4))THEN
      IBD_Config2=Merge(1_I2B, 2_I2B, G(1)==G(3))
    ELSE
      IBD_Config2=Merge(3_I2B, 4_I2B, ANY(G(1)==G(3:4)))
    END IF
  ELSE
    IF(G(3)==G(4))THEN
      IBD_Config2=Merge(5_I2B, 6_I2B, ANY(G(3)==G(1:2)))
    ELSE
      IF(G(1)==G(3))THEN
        IBD_Config2=Merge(7_I2B, 8_I2B, G(2)==G(4))
      ELSE IF(G(1)==G(4) .OR. ANY(G(2)==G(3:4)))THEN
        IBD_Config2=8_I2B
      ELSE
        IBD_Config2=9_I2B
      END IF
    END IF
  END IF
END FUNCTION IBD_Config2
!
FUNCTION IBD_Config3(G,IBD12)
  USE Variables, ONLY : I1B, I2B
  IMPLICIT NONE
  INTEGER(I2B) :: IBD_Config3
  INTEGER(I1B) :: G(6)
  INTEGER :: IBD12
  SELECT CASE (IBD12)
  CASE (1)
    IF(G(5)==G(6))THEN
      IBD_Config3=Merge(1_I2B, 8_I2B, G(5)==G(1))
    ELSE
      IBD_Config3=Merge(2_I2B, 20_I2B, ANY(G(1)==G(5:6)))
    END IF
  CASE (2)
    IF(G(5)==G(6))THEN
      IF(G(5)==G(1))THEN
        IBD_Config3=9_I2B
      ELSE IF(G(5)==G(3))THEN
        IBD_Config3=10_I2B
      ELSE
        IBD_Config3=15_I2B
      END IF
    ELSE
      IF(ANY(G(1)==G(5:6)))THEN
        IBD_Config3=Merge(12_I2B, 33_I2B, ANY(G(3)==G(5:6)))
      ELSE IF(ANY(G(3)==G(5:6)))THEN
        IBD_Config3=36_I2B
      ELSE
        IBD_Config3=41_I2B
      END IF
    END IF
  CASE (3)
    IF(G(5)==G(6))THEN
      IF(G(5)==G(1))THEN
        IBD_Config3=3_I2B
      ELSE IF(ANY(G(5)==G(3:4)))THEN
        IBD_Config3=11_I2B
      ELSE
        IBD_Config3=32_I2B
      END IF
    ELSE
      IF(ALL(G(3:4)==G(5:6)))THEN
        IBD_Config3=5_I2B
      ELSE IF(ANY(G(1)==G(5:6)))THEN
        IBD_Config3=23_I2B
      ELSE IF(ANY(G(3)==G(5:6)).OR.ANY(G(4)==G(5:6)))THEN
        IBD_Config3=26_I2B
      ELSE
        IBD_Config3=53_I2B
      END IF
    END IF
  CASE (4)
    IF(G(5)==G(6))THEN
      IF(G(5)==G(1))THEN
        IBD_Config3=21_I2B
      ELSE IF(ANY(G(5)==G(3:4)))THEN
        IBD_Config3=37_I2B
      ELSE
        IBD_Config3=42_I2B
      END IF
    ELSE
      IF(ALL(G(5:6)==G(3:4)))THEN
        IBD_Config3=16_I2B
      ELSE IF(ANY(G(5)==G(3:4)))THEN
        IBD_Config3=Merge(29_I2B, 47_I2B, G(1)==G(6))
      ELSE IF(ANY(G(6)==G(3:4)))THEN
        IBD_Config3=Merge(29_I2B, 47_I2B, G(1)==G(5))
      ELSE IF(ANY(G(1)==G(5:6)))THEN
        IBD_Config3=54_I2B
      ELSE
        IBD_Config3=60_I2B
      END IF
    END IF
  CASE (5)
    IF(G(5)==G(6))THEN
      IF(G(5)==G(3))THEN
        IBD_Config3=4_I2B
      ELSE IF(ANY(G(5)==G(1:2)))THEN
        IBD_Config3=13_I2B
      ELSE
        IBD_Config3=34_I2B
      END IF
    ELSE
      IF(ALL(G(5:6)==G(1:2)))THEN
        IBD_Config3=6_I2B
      ELSE IF(ANY(G(3)==G(5:6)))THEN
        IBD_Config3=24_I2B
      ELSE IF(ANY(G(5)==G(1:2)).OR.ANY(G(6)==G(1:2)))THEN
        IBD_Config3=27_I2B
      ELSE
        IBD_Config3=56_I2B
      END IF
    END IF
  CASE (6)
    IF(G(5)==G(6))THEN
      IF(G(5)==G(3))THEN
        IBD_Config3=22_I2B
      ELSE IF(ANY(G(5)==G(1:2)))THEN
        IBD_Config3=35_I2B
      ELSE
        IBD_Config3=43_I2B
      END IF
    ELSE
      IF(ALL(G(5:6)==G(1:2)))THEN
        IBD_Config3=17_I2B
      ELSE IF(ANY(G(5)==G(1:2)))THEN
        IBD_Config3=Merge(30_I2B, 48_I2B, G(3)==G(6))
      ELSE IF(ANY(G(6)==G(1:2)))THEN
        IBD_Config3=Merge(30_I2B, 48_I2B, G(3)==G(5))
      ELSE IF(ANY(G(3)==G(5:6)))THEN
        IBD_Config3=55_I2B
      ELSE
        IBD_Config3=61_I2B
      END IF
    END IF
  CASE (7)
    IF(G(5)==G(6))THEN
      IBD_Config3=Merge(7_I2B, 18_I2B, ANY(G(5)==G(1:2)))
    ELSE
      IF(ALL(G(5:6)==G(1:2)))THEN
        IBD_Config3=14_I2B
      ELSE IF(ANY(G(5)==G(1:2)).OR.ANY(G(6)==G(1:2)))THEN
        IBD_Config3=38_I2B
      ELSE
        IBD_Config3=46_I2B
      END IF
    END IF
  CASE (8)
    IF(G(5)==G(6))THEN
      IF(ANY(G(5)==G(1:2)).AND.ANY(G(5)==G(3:4)))THEN
        IBD_Config3=25_I2B
      ELSE IF(ANY(G(5)==G(1:2)))THEN
        IBD_Config3=28_I2B
      ELSE IF(ANY(G(5)==G(3:4)))THEN
        IBD_Config3=31_I2B
      ELSE
        IBD_Config3=49_I2B
      END IF
    ELSE
      IF(ALL(G(5:6)==G(1:2)))THEN
        IBD_Config3=39_I2B
      ELSE IF(ALL(G(5:6)==G(3:4)))THEN
        IBD_Config3=40_I2B
      ELSE IF(ANY(G(5)==G(1:2)))THEN
        IF(ANY(G(5)==G(3:4)))THEN
          IBD_Config3=59_I2B
        ELSE IF(ANY(G(6)==G(3:4)))THEN
          IBD_Config3=19_I2B
        ELSE
          IBD_Config3=50_I2B
        END IF
      ELSE IF(ANY(G(6)==G(1:2)))THEN
        IF(ANY(G(6)==G(3:4)))THEN
          IBD_Config3=59_I2B
        ELSE IF(ANY(G(5)==G(3:4)))THEN
          IBD_Config3=19_I2B
        ELSE
          IBD_Config3=50_I2B
        END IF
      ELSE IF(ANY(G(5)==G(3:4)).OR.ANY(G(6)==G(3:4)))THEN
        IBD_Config3=51_I2B
      ELSE
        IBD_Config3=63_I2B
      END IF
    END IF
  CASE (9)
    IF(G(5)==G(6))THEN
      IF(ANY(G(5)==G(1:2)))THEN
        IBD_Config3=57_I2B
      ELSE IF(ANY(G(5)==G(3:4)))THEN
        IBD_Config3=58_I2B
      ELSE
        IBD_Config3=62_I2B
      END IF
    ELSE
      IF(ALL(G(5:6)==G(1:2)))THEN
        IBD_Config3=45_I2B
      ELSE IF(ALL(G(5:6)==G(3:4)))THEN
        IBD_Config3=44_I2B
      ELSE IF(ANY(G(5)==G(1:2)))THEN
        IF(ANY(G(6)==G(3:4)))THEN
          IBD_Config3=52_I2B
        ELSE
          IBD_Config3=64_I2B
        END IF
      ELSE IF(ANY(G(6)==G(1:2)))THEN
        IF(ANY(G(5)==G(3:4)))THEN
          IBD_Config3=52_I2B
        ELSE
          IBD_Config3=64_I2B
        END IF
      ELSE IF(ANY(G(5)==G(3:4)).OR.ANY(G(6)==G(3:4)))THEN
        IBD_Config3=65_I2B
      ELSE
        IBD_Config3=66_I2B
      END IF
    END IF
  END SELECT
END FUNCTION IBD_Config3
!
SUBROUTINE choldc(A,N,NP,P)
  USE Variables, ONLY: ErrorMessage
  DOUBLE PRECISION A(NP,NP),P(NP),SUM
  DO I=1,N
    DO J=I,N
      SUM=A(I,J)
      DO K=I-1,1,-1
        SUM=SUM-A(I,K)*A(J,K)
      END DO
      IF(I.EQ.J)THEN
        IF(SUM.LE.0.D0)THEN
          ErrorMessage(1) = 'Program stops because CHOLDC FAILED!'
          call StopOnError(1, 'choldc')
        end if
        P(I)=SQRT(SUM)
      ELSE
        A(J,I)=SUM/P(I)
      ENDIF
    END DO
  END DO
END SUBROUTINE choldc
!
SUBROUTINE cholsl(A,N,NP,P,B,X)
  DOUBLE PRECISION A(NP,NP),P(NP),B(NP),X(NP),SUM
  DO I=1,N
    SUM=B(I)
    DO K=I-1,1,-1
      SUM=SUM-A(I,K)*X(K)
    END DO
    X(I)=SUM/P(I)
  END DO
  DO I=N,1,-1
    SUM=X(I)
    DO K=I+1,N
      SUM=SUM-A(K,I)*X(K)
    END DO
    X(I)=SUM/P(I)
  END DO
END SUBROUTINE cholsl
!
SUBROUTINE Coefficient
  USE VARIABLES
  USE IBDVariables
  IMPLICIT NONE
  INTEGER:: L
  REAL(DP):: Q(99),QQ,e1,e2,e3,e4,e5,X(4)
  DO L = 1, NumLoci
    IF(GConfig(L) == 0_I2B) THEN
      coe(L, :) = 0._DP
      CYCLE
    END IF
    e1 = err(LocusSelectPtr(L))
    e2 = e1 / ObsNumAllele(L)
    e3 = 1-e1+e2
    e4 = -1+e1-2*e2
    e5 = -1+e1
    SELECT CASE (GConfig(L))
      CASE (1_I2B)
        Q(1:21)=(/e3**4*A(1,L)**2,e3**3*A(1,L)*(e2+A(1,L)-e1*A(1,L)), &
          e3**2*A(1,L)*(e2**2+e4*e5*A(1,L)),&
          e3*A(1,L)*(e2**3-e5*(3*e2**2-3*e2*e5+e5**2)*A(1,L)), &
          A(1,L)*(e2**2+e4*e5*A(1,L))**2,&
          e3**2*A(1,L)**2*(e2**2+e4*e5*A(1,L)), &
          e3**2*A(1,L)*(e2+A(1,L)-e1*A(1,L))**2,&
          e3**3*A(1,L)**2*(e2+A(1,L)-e1*A(1,L)),&
          e3*A(1,L)**2*(e2**3-e5*(3*e2**2-3*e2*e5+e5**2)*A(1,L)),&
          e3*A(1,L)*(e2+A(1,L)-e1*A(1,L))*(e2**2+e4*e5*A(1,L)),&
          A(1,L)*(e2**4-4*e2**3*e5*A(1,L)+e5**4*A(1,L)**2+ &
          3*e2**2*e5**2*A(1,L)*(1+A(1,L))-e2*e5**3*A(1,L)*(1+3*A(1,L))),&
          A(1,L)**2*(e2**2+e4*e5*A(1,L))**2,&
          A(1,L)*(e2+A(1,L)-e1*A(1,L))**2*(e2**2+e4*e5*A(1,L)), &
          e3**2*A(1,L)**2*(e2+A(1,L)-e1*A(1,L))**2,&
          e3*A(1,L)**2*(e2+A(1,L)-e1*A(1,L))*(e2**2+e4*e5*A(1,L)),&
          A(1,L)**2*(e2**4-4*e2**3*e5*A(1,L)+e5**4*A(1,L)**2+ &
          3*e2**2*e5**2*A(1,L)*(1+A(1,L))-e2*e5**3*A(1,L)*(1+3*A(1,L))),&
          e3*A(1,L)*(e2+A(1,L)-e1*A(1,L))**3,&
          A(1,L)**2*(e2+A(1,L)-e1*A(1,L))**2*(e2**2+e4*e5*A(1,L)), &
          e3*A(1,L)**2*(e2+A(1,L)-e1*A(1,L))**3,&
          e2**4-4*e2**3*e5*A(1,L)+6*e2**2*e5**2*A(1,L)- &
          4*e2*e5**3*A(1,L)+e5**4*A(1,L),(e2+A(1,L)-e1*A(1,L))**4/)
        coe(L,:)=(/e3**4*A(1,L),Q(1),Q(2),Q(2),Q(1),Q(1),Q(3), &
          A(1,L)*Q(20),Q(3),Q(3),Q(4),Q(1),Q(4),Q(1),Q(5),Q(6),Q(6),Q(5),Q(6),&
          A(1,L)**2*Q(20),Q(7),Q(7),Q(8),Q(8),Q(7),Q(9),Q(9),Q(10),&
          Q(8),Q(8),Q(10),Q(11),Q(6),Q(11),Q(10),Q(6),Q(10),Q(6),Q(8),&
          Q(8),Q(12),Q(13),Q(13),Q(14),Q(14),Q(12),Q(15),Q(15),Q(13),&
          Q(15),Q(15),Q(14),Q(16),Q(14),Q(14),Q(16),Q(17),Q(17),&
          Q(14),Q(18),Q(18),A(1,L)*Q(21),Q(18),Q(19),Q(19),A(1,L)**2*Q(21)/)
      CASE (2_I2B)
        QQ=Merge(A(6,L), A(5,L), G(1,L)==G(5,L))
        Q(1:14)=(/e2*e3*(2*e2**2-2*e2*e5+e5**2)*QQ*A(1,L), &
          2*e2**2*e3**2*QQ*A(1,L),2*e2*e3*QQ*A(1,L)*(e2**2+e4*e5*A(1,L)),&
          -(e4*(e2-e1*e2+e2**2+e5**2)*QQ*A(1,L)*(e2+A(1,L)-e1*A(1,L))),&
          -(e4*QQ*A(1,L)*(e2**3-e5*(3*e2**2-3*e2*e5+e5**2)*A(1,L))), &
          -(e2*e3*e4*QQ*A(1,L)*(e2+A(1,L)-e1*A(1,L))),&
          (2*e2**2-2*e2*e5+e5**2)*QQ*A(1,L)*(e2**2+e4*e5*A(1,L)), &
          2*QQ*A(1,L)*(e2**2+e4*e5*A(1,L))**2,&
          2*e2*e3*QQ*A(1,L)*(e2+A(1,L)-e1*A(1,L))**2, &
          -(e4*QQ*A(1,L)*(e2+A(1,L)-e1*A(1,L))*(e2**2+e4*e5*A(1,L))),&
          2*QQ*A(1,L)*(e2**4-4*e2**3*e5*A(1,L)+e5**4*A(1,L)**2+ &
          3*e2**2*e5**2*A(1,L)*(1+A(1,L))-e2*e5**3*A(1,L)*(1+3*A(1,L))),&
          (2*e2**2-2*e2*e5+e5**2)*QQ*A(1,L)*(e2+A(1,L)-e1*A(1,L))**2,&
          2*QQ*A(1,L)*(e2+A(1,L)-e1*A(1,L))**2*(e2**2+e4*e5*A(1,L)), &
          -(e4*QQ*A(1,L)*(e2+A(1,L)-e1*A(1,L))**3)/)
        coe(L,:)=(/0._DP,(2*e2**4-4*e2**3*e5+6*e2**2*e5**2- &
          4*e2*e5**3+e5**4)*QQ*A(1,L),0._DP,0._DP,Q(1), &
          Q(1),0._DP,0._DP,0._DP,0._DP,0._DP,Q(2),0._DP,Q(2),0._DP,&
          Q(3),Q(3),0._DP,Q(3),2*QQ*A(1,L)*(e2**4-4*e2**3*e5*A(1,L)+ &
          6*e2**2*e5**2*A(1,L)-4*e2*e5**3*A(1,L)+e5**4*A(1,L)),0._DP,0._DP,&
          Q(4),Q(4),0._DP,Q(5),Q(5),0._DP,Q(6),Q(6),0._DP,0._DP,Q(7), &
          0._DP,0._DP,Q(7),0._DP,Q(7),Q(6),Q(6),Q(8),0._DP,0._DP,Q(9), &
          Q(9),Q(8),Q(10),Q(10),0._DP,Q(10),Q(10),&
          Q(9),Q(11),Q(12),Q(12),Q(11),0._DP,0._DP,Q(12),Q(13),&
          Q(13),0._DP,Q(13),Q(14),Q(14),2*QQ*A(1,L)*(e2+A(1,L)-e1*A(1,L))**4/)
      CASE (3_I2B)
        QQ=Merge(A(4,L), A(3,L), G(1,L)==G(3,L))
        Q(1:13)=(/2*e2*e3**3*A(1,L)**2,2*e2*A(1,L)*(e2**2+ &
          e4*e5*A(1,L))*(e2-e5*(QQ+A(1,L))),&
          2*e2*e3*A(1,L)**2*(e2**2+e4*e5*A(1,L)),&
          2*e2*e3**2*A(1,L)**2*(e2-e5*(QQ+A(1,L))),&
          e2*e3*A(1,L)**2*(2*e2**2+e5**2*(QQ+2*A(1,L))-e2*e5*(1+QQ+3*A(1,L))),&
          e3**2*A(1,L)**2*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(1,L))),&
          2*e2*e3**2*A(1,L)**2*(e2+A(1,L)-e1*A(1,L)),&
          2*e2*A(1,L)**2*(e2**2+e4*e5*A(1,L))*(e2-e5*(QQ+A(1,L))),&
          2*e2*e3*A(1,L)**2*(e2+A(1,L)-e1*A(1,L))**2,&
          2*e3**2*(e2+QQ-e1*QQ)*A(1,L)**2*(e2+A(1,L)-e1*A(1,L)),&
          e3*A(1,L)**2*(e2+A(1,L)-e1*A(1,L))*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(1,L))), &
          e2+A(1,L)-e1*A(1,L),e2+QQ-e1*QQ/)
        coe(L,1:34)=(/2*e2*e3**3*A(1,L), Q(1), &
          e3**2*A(1,L)*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(1,L))), &
          2*e2*e3**2*A(1,L)*Q(12), Q(1), Q(1), &
          e2*e3*A(1,L)*(2*e2**2+e5**2*(QQ+2*A(1,L))-e2*e5*(1+QQ+3*A(1,L))),&
          2*e2*A(1,L)*(e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(QQ+3*A(1,L))),&
          2*e2*e3**2*A(1,L)*(e2-e5*(QQ+A(1,L))), 2*e2*e3*A(1,L)*(e2**2+e4*e5*A(1,L)),&
          e2*A(1,L)*(2*e2**3-2*e5**3*A(1,L)-e2**2*e5*(1+QQ+5*A(1,L))+ &
          e2*e5**2*(QQ+6*A(1,L))), Q(1),&
          2*e2*e3*A(1,L)*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L))), &
          Q(1), Q(2), Q(3), Q(4), Q(2), Q(5),&
          2*e2*A(1,L)**2*(e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)- &
          e2**2*e5*(QQ+3*A(1,L))), 2*e3**2*A(1,L)*Q(12)*Q(13),&
          2*e2*e3*A(1,L)*Q(12)**2, Q(6), Q(7), &
          e3*A(1,L)*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(1,L)))*Q(12),&
          e2*A(1,L)**2*(2*e2**3-2*e5**3*A(1,L)-e2**2*e5*(1+QQ+ &
          5*A(1,L))+e2*e5**2*(QQ+6*A(1,L))),&
          2*e2*e3*A(1,L)**2*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L))),&
          e3*A(1,L)*(2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+&
          3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L))), Q(6), Q(7),&
          e2*A(1,L)*(2*e2**2+e5**2*(QQ+2*A(1,L))-e2*e5*(1+QQ+3*A(1,L)))*Q(12),&
          A(1,L)*(2*e2**4+e5**4*QQ*A(1,L)+2*e2**2*e5**2*A(1,L)*(2+ &
          2*QQ+A(1,L))-e2*e5**3*A(1,L)*(1+3*QQ+A(1,L))-2*e2**3*e5*(QQ+3*A(1,L))),&
          Q(4), 2*e2*A(1,L)*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L)))*Q(12)/)
          coe(L,35:66)=(/2*e2*e3*A(1,L)*(e2-e5*(QQ+A(1,L)))*Q(12),Q(3),&
          A(1,L)*(e2**2+e4*e5*A(1,L))*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(1,L))),&
          Q(5),Q(6),Q(7),Q(8),2*A(1,L)*(e2**2+e4*e5*A(1,L))*Q(12)*Q(13),&
          2*e2*A(1,L)*(e2-e5*(QQ+A(1,L)))*Q(12)**2,Q(9),Q(10),Q(8),&
          A(1,L)**2*(e2**2+e4*e5*A(1,L))*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(1,L))),&
          2*e2*e3*A(1,L)**2*(e2-e5*(QQ+A(1,L)))*Q(12),&
          A(1,L)*(2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+ &
          3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L)))*Q(12),&
          e3*A(1,L)**2*(2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+ &
          3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L))),&
          e2*A(1,L)**2*(2*e2**2+e5**2*(QQ+2*A(1,L))- &
          e2*e5*(1+QQ+3*A(1,L)))*Q(12),Q(11),&
          A(1,L)**2*(2*e2**4+e5**4*QQ*A(1,L)+2*e2**2*e5**2*A(1,L)*(2+2*QQ+A(1,L))-&
          e2*e5**3*A(1,L)*(1+3*QQ+A(1,L))-2*e2**3*e5*(QQ+3*A(1,L))),Q(10),Q(9),&
          2*e2*A(1,L)**2*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+ &
          2*A(1,L)))*Q(12),2*e3*A(1,L)*Q(12)**2*Q(13),&
          A(1,L)*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(1,L)))* &
          Q(12)**2,Q(11),2*A(1,L)**2*(e2**2+e4*e5*A(1,L))*Q(12)*Q(13),&
          2*e2*A(1,L)**2*(e2-e5*(QQ+A(1,L)))*Q(12)**2,2*A(1,L)*Q(12)**3*Q(13),&
          A(1,L)**2*(2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+ &
          3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L)))*Q(12),&
          2*e3*A(1,L)**2*Q(12)**2*Q(13),A(1,L)**2*(2*e2**2+ &
          e5**2*QQ-e2*e5*(1+QQ+A(1,L)))*Q(12)**2,2*A(1,L)**2*Q(12)**3*Q(13)/)
      CASE (4_I2B)
        QQ=Merge(A(2,L), A(1,L), G(1,L)==G(3,L))
        Q(1:13)=(/2*e2*e3**3*A(3,L)**2,2*e2*A(3,L)*(e2**2+ &
          e4*e5*A(3,L))*(e2-e5*(QQ+A(3,L))),&
          2*e2*e3**2*A(3,L)**2*(e2-e5*(QQ+A(3,L))), &
          2*e2*e3*A(3,L)**2*(e2**2+e4*e5*A(3,L)),&
          e2*e3*A(3,L)**2*(2*e2**2+e5**2*(QQ+2*A(3,L))-e2*e5*(1+QQ+3*A(3,L))),&
          2*e2*e3**2*A(3,L)**2*(e2+A(3,L)-e1*A(3,L)), &
          e3**2*A(3,L)**2*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(3,L))),&
          2*e2*A(3,L)**2*(e2**2+e4*e5*A(3,L))*(e2-e5*(QQ+A(3,L))),&
          2*e3**2*(e2+QQ-e1*QQ)*A(3,L)**2*(e2+A(3,L)-e1*A(3,L)), &
          2*e2*e3*A(3,L)**2*(e2+A(3,L)-e1*A(3,L))**2,&
          e3*A(3,L)**2*(e2+A(3,L)-e1*A(3,L))*(2*e2**2+e5**2*QQ- &
          e2*e5*(1+QQ+A(3,L))),e2+A(3,L)-e1*A(3,L),e2+QQ-e1*QQ/)
        coe(L,1:42)=(/2*e2*e3**3*A(3,L),Q(1),2*e2*e3**2*A(3,L)*Q(12), &
          e3**2*A(3,L)*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(3,L))),Q(1),&
          Q(1),e2*e3*A(3,L)*(2*e2**2+e5**2*(QQ+2*A(3,L))-e2*e5*(1+QQ+3*A(3,L))),&
          2*e2*A(3,L)*(e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(QQ+3*A(3,L))),&
          2*e2*e3*A(3,L)*(e2**2+e4*e5*A(3,L)),2*e2*e3**2*A(3,L)*(e2-e5*(QQ+A(3,L))),&
          2*e2*e3*A(3,L)*(e2**2+e5**2*A(3,L)-e2*e5*(QQ+2*A(3,L))),Q(1),&
          e2*A(3,L)*(2*e2**3-2*e5**3*A(3,L)-e2**2*e5*(1+QQ+5*A(3,L))+ &
          e2*e5**2*(QQ+6*A(3,L))),Q(1),Q(2),Q(3),Q(4),&
          Q(2),Q(5),2*e2*A(3,L)**2*(e2**3+3*e2*e5**2*A(3,L)- &
          e5**3*A(3,L)-e2**2*e5*(QQ+3*A(3,L))),2*e2*e3*A(3,L)*Q(12)**2,&
          2*e3**2*A(3,L)*Q(12)*Q(13),Q(6),Q(7), &
          e3*A(3,L)*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(3,L)))*Q(12),&
          2*e2*e3*A(3,L)**2*(e2**2+e5**2*A(3,L)-e2*e5*(QQ+2*A(3,L))),&
          e2*A(3,L)**2*(2*e2**3-2*e5**3*A(3,L)-e2**2*e5*(1+QQ+ &
          5*A(3,L))+e2*e5**2*(QQ+6*A(3,L))),&
          e2*A(3,L)*(2*e2**2+e5**2*(QQ+2*A(3,L))-e2*e5*(1+QQ+ &
          3*A(3,L)))*Q(12),Q(6),Q(7),&
          e3*A(3,L)*(2*e2**3-e5**3*QQ*A(3,L)+e2*e5**2*A(3,L)*(1+ &
          3*QQ+A(3,L))-2*e2**2*e5*(QQ+2*A(3,L))),&
          2*e2*A(3,L)*(e2**2+e5**2*A(3,L)-e2*e5*(QQ+2*A(3,L)))*Q(12),Q(4),&
          A(3,L)*(2*e2**4+e5**4*QQ*A(3,L)+2*e2**2*e5**2*A(3,L)*(2+ &
          2*QQ+A(3,L))-e2*e5**3*A(3,L)*(1+3*QQ+A(3,L))-&
          2*e2**3*e5*(QQ+3*A(3,L))),A(3,L)*(e2**2+e4*e5*A(3,L))* &
          (2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(3,L))),Q(3),&
          2*e2*e3*A(3,L)*(e2-e5*(QQ+A(3,L)))*Q(12),Q(5),Q(6),Q(7),Q(8), &
          2*e2*A(3,L)*(e2-e5*(QQ+A(3,L)))*Q(12)**2/)
          coe(L,43:66)=(/2*A(3,L)*(e2**2+e4*e5*A(3,L))*Q(12)*Q(13),Q(9),Q(10),Q(8),&
          2*e2*e3*A(3,L)**2*(e2-e5*(QQ+A(3,L)))*Q(12),&
          A(3,L)**2*(e2**2+e4*e5*A(3,L))*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(3,L))),&
          A(3,L)*(2*e2**3-e5**3*QQ*A(3,L)+e2*e5**2*A(3,L)*(1+ &
          3*QQ+A(3,L))-2*e2**2*e5*(QQ+2*A(3,L)))*Q(12),&
          e2*A(3,L)**2*(2*e2**2+e5**2*(QQ+2*A(3,L))-e2*e5*(1+QQ+3*A(3,L)))*Q(12),&
          e3*A(3,L)**2*(2*e2**3-e5**3*QQ*A(3,L)+e2*e5**2*A(3,L)* &
          (1+3*QQ+A(3,L))-2*e2**2*e5*(QQ+2*A(3,L))),Q(11),&
          2*e2*A(3,L)**2*(e2**2+e5**2*A(3,L)-e2*e5*(QQ+2*A(3,L)))*Q(12),Q(10),Q(9),&
          A(3,L)**2*(2*e2**4+e5**4*QQ*A(3,L)+2*e2**2*e5**2*A(3,L)*(2+2*QQ+A(3,L))-&
             e2*e5**3*A(3,L)*(1+3*QQ+A(3,L))-2*e2**3*e5*(QQ+3*A(3,L))),&
          A(3,L)*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(3,L)))*Q(12)**2, &
          2*e3*A(3,L)*Q(12)**2*Q(13),Q(11),&
          2*e2*A(3,L)**2*(e2-e5*(QQ+A(3,L)))*Q(12)**2, &
          2*A(3,L)**2*(e2**2+e4*e5*A(3,L))*Q(12)*Q(13),2*A(3,L)*Q(12)**3*Q(13),&
          A(3,L)**2*(2*e2**3-e5**3*QQ*A(3,L)+e2*e5**2*A(3,L)* &
          (1+3*QQ+A(3,L))-2*e2**2*e5*(QQ+2*A(3,L)))*Q(12),&
          A(3,L)**2*(2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(3,L)))*Q(12)**2, &
          2*e3*A(3,L)**2*Q(12)**2*Q(13),2*A(3,L)**2*Q(12)**3*Q(13)/)
      CASE (5_I2B)
        QQ=Merge(A(4,L), A(3,L), G(1,L)==G(3,L))
        Q(1:5)=(/2*e2*e3*(2*e2**2-2*e2*e5+e5**2)*QQ*A(1,L), &
          -2*e2*e3*e4*QQ*A(1,L)*(e2+A(1,L)-e1*A(1,L)),&
          4*e2*QQ*A(1,L)*(e2**2+e4*e5*A(1,L))*(e2-e5*(QQ+A(1,L))), &
          e2+A(1,L)-e1*A(1,L),e2+QQ-e1*QQ/)
        coe(L,1:37)=(/0._DP,Q(1),0._DP,0._DP,(2*e2**2-2*e2*e5+e5**2)**2*QQ*A(1,L),&
          4*e2**2*e3**2*QQ*A(1,L),&
          0._DP,0._DP,0._DP,0._DP,0._DP,Q(1),0._DP,Q(1),0._DP,&
          2*(2*e2**2-2*e2*e5+e5**2)*QQ*A(1,L)*(e2**2+e4*e5*A(1,L)), &
          4*e2**2*e3*QQ*A(1,L)*(e2-e5*(QQ+A(1,L))),0._DP,&
          QQ*A(1,L)*(4*e2**4-4*e2*e5**3*A(1,L)+e5**4*A(1,L)- &
          2*e2**3*e5*(2+QQ+3*A(1,L))+e2**2*e5**2*(1+2*QQ+8*A(1,L))),&
          4*e2*QQ*A(1,L)*(e2**3+3*e2*e5**2*A(1,L)- &
          e5**3*A(1,L)-e2**2*e5*(QQ+3*A(1,L))),0._DP,0._DP, &
          -(e4*QQ*A(1,L)*(2*e2**3-e5**3*QQ+e2*e5**2*(1+QQ+ &
          A(1,L))-e2**2*e5*(2+QQ+A(1,L)))),Q(2),0._DP, &
          -(e4*QQ*A(1,L)*(2*e2**3+4*e2*e5**2*A(1,L)- &
          e5**3*A(1,L)-e2**2*e5*(QQ+5*A(1,L)))), &
          -2*e2*e4*QQ*A(1,L)*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L))),0._DP, &
          -(e4*QQ*A(1,L)*(2*e2**3-e5**3*A(1,L)+e2*e5**2*(1+ &
          QQ+A(1,L))-e2**2*e5*(2+QQ+A(1,L)))),Q(2),0._DP,0._DP, &
          2*e2*(2*e2**2-2*e2*e5+e5**2)*QQ*A(1,L)*(e2-e5*(QQ+A(1,L))), &
          0._DP,0._DP,4*e2*e3*QQ*A(1,L)*(e2**2+e4*e5*A(1,L)),0._DP/)
         coe(L,38:66)=(/e2*QQ*A(1,L)*(4*e2**3-e5**3*(QQ+3*A(1,L))-2*e2**2*e5*(2+ &
          QQ+3*A(1,L))+e2*e5**2*(1+2*QQ+8*A(1,L))), &
          -(e2*e3*e4*QQ*A(1,L)*(2*e2-e5*(QQ+A(1,L)))), &
          -(e4*(2*e2**2-2*e2*e5+e5**2)*QQ*A(1,L)*Q(4)),Q(3),0._DP,0._DP, &
          2*(2*e2**2-2*e2*e5+e5**2)*QQ*A(1,L)*Q(4)**2, &
          4*e2*e3*QQ*A(1,L)*Q(4)*Q(5),Q(3), &
          -(e4*QQ*A(1,L)*(e2**2+e4*e5*A(1,L))*(2*e2-e5*(QQ+A(1,L)))), &
          -2*e2*e4*QQ*A(1,L)*(e2-e5*(QQ+A(1,L)))*Q(4),0._DP, &
          -(e4*QQ*A(1,L)*(2*e2**3-e5**3*QQ*A(1,L)+ &
          e2*e5**2*A(1,L)*(1+3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L)))), &
          -(e4*QQ*A(1,L)*(2*e2**2+e5**2*A(1,L)-e2*e5*(QQ+3*A(1,L)))*Q(4)), &
          QQ*A(1,L)*(4*e2**3-e5**3*A(1,L)-2*e2**2*e5*(2+ &
          QQ+A(1,L))+e2*e5**2*(1+2*QQ+2*A(1,L)))*Q(4), &
          2*QQ*A(1,L)*(2*e2**4+e5**4*QQ*A(1,L)+2*e2**2*e5**2* &
          A(1,L)*(2+2*QQ+A(1,L))- &
          e2*e5**3*A(1,L)*(1+3*QQ+A(1,L))-2*e2**3*e5*(QQ+3*A(1,L))), &
          2*(2*e2**2-2*e2*e5+e5**2)*QQ*A(1,L)*Q(4)*Q(5), &
          4*e2*e3*QQ*A(1,L)*Q(4)**2,4*e2*QQ*A(1,L)*(e2**2+ &
          e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L)))*Q(4),0._DP,0._DP, &
          QQ*A(1,L)*(4*e2**3-e5**3*QQ-2*e2**2*e5*(2+QQ+A(1,L))+ &
          e2*e5**2*(1+2*QQ+2*A(1,L)))*Q(4), &
          4*QQ*A(1,L)*(e2**2+e4*e5*A(1,L))*Q(4)*Q(5), &
          4*e2*QQ*A(1,L)*(e2-e5*(QQ+A(1,L)))*Q(4)**2,0._DP, &
          2*QQ*A(1,L)*(2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+ &
          3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L)))*Q(4), &
          -2*e4*QQ*A(1,L)*Q(4)**2*Q(5),-(e4*QQ*A(1,L)*(2*e2-e5*(QQ+ &
          A(1,L)))*Q(4)**2),4*QQ*A(1,L)*Q(4)**3*Q(5)/)
      CASE (6_I2B)
        QQ=Merge(A(2,L), A(1,L), G(1,L)==G(3,L))
        Q(1:5)=(/2*e2*e3*(2*e2**2-2*e2*e5+e5**2)*QQ*A(3,L), &
          -2*e2*e3*e4*QQ*A(3,L)*(e2+A(3,L)-e1*A(3,L)), &
          4*e2*QQ*A(3,L)*(e2**2+e4*e5*A(3,L))*(e2-e5*(QQ+A(3,L))), &
          e2+A(3,L)-e1*A(3,L),e2+QQ-e1*QQ/)
        coe(L,1:35)=(/0._DP,Q(1),0._DP,0._DP,4*e2**2*e3**2*QQ*A(3,L),&
          (2*e2**2-2*e2*e5+e5**2)**2*QQ*A(3,L),0._DP,0._DP,0._DP,&
          0._DP,0._DP,Q(1),0._DP,Q(1),0._DP, &
          4*e2**2*e3*QQ*A(3,L)*(e2-e5*(QQ+A(3,L))), &
          2*(2*e2**2-2*e2*e5+e5**2)*QQ*A(3,L)*(e2**2+e4*e5*A(3,L)),0._DP, &
          QQ*A(3,L)*(4*e2**4-4*e2*e5**3*A(3,L)+ &
          e5**4*A(3,L)-2*e2**3*e5*(2+QQ+3*A(3,L))+ &
          e2**2*e5**2*(1+2*QQ+8*A(3,L))),4*e2*QQ*A(3,L)* &
          (e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(QQ+3*A(3,L))), &
          0._DP,0._DP,Q(2), &
          -(e4*QQ*A(3,L)*(2*e2**3-e5**3*QQ+e2*e5**2*(1+ &
          QQ+A(3,L))-e2**2*e5*(2+QQ+A(3,L)))),0._DP, &
          -2*e2*e4*QQ*A(3,L)*(e2**2+e5**2*A(3,L)-e2*e5*(QQ+2*A(3,L))), &
          -(e4*QQ*A(3,L)*(2*e2**3+4*e2*e5**2*A(3,L)- &
          e5**3*A(3,L)-e2**2*e5*(QQ+5*A(3,L)))),0._DP,Q(2), &
          -(e4*QQ*A(3,L)*(2*e2**3-e5**3*A(3,L)+e2*e5**2*(1+ &
          QQ+A(3,L))-e2**2*e5*(2+QQ+A(3,L)))),0._DP,0._DP, &
          4*e2*e3*QQ*A(3,L)*(e2**2+e4*e5*A(3,L)),0._DP,0._DP/)
         coe(L,36:66)=(/2*e2*(2*e2**2-2*e2*e5+e5**2)*QQ*A(3,L)*(e2-e5*(QQ+A(3,L))),0._DP, &
          e2*QQ*A(3,L)*(4*e2**3-e5**3*(QQ+3*A(3,L))- &
          2*e2**2*e5*(2+QQ+3*A(3,L))+e2*e5**2*(1+2*QQ+8*A(3,L))), &
          -(e4*(2*e2**2-2*e2*e5+e5**2)*QQ*A(3,L)*Q(4)), &
          -(e2*e3*e4*QQ*A(3,L)*(2*e2-e5*(QQ+A(3,L)))),Q(3),0._DP,0._DP, &
          4*e2*e3*QQ*A(3,L)*Q(4)*Q(5),2*(2*e2**2-2*e2*e5+e5**2)*QQ*A(3,L)*Q(4)**2, &
          Q(3),-2*e2*e4*QQ*A(3,L)*(e2-e5*(QQ+A(3,L)))*Q(4), &
          -(e4*QQ*A(3,L)*(e2**2+e4*e5*A(3,L))*(2*e2-e5*(QQ+A(3,L)))),0._DP, &
          -(e4*QQ*A(3,L)*(2*e2**2+e5**2*A(3,L)-e2*e5*(QQ+3*A(3,L)))*Q(4)), &
          -(e4*QQ*A(3,L)*(2*e2**3-e5**3*QQ*A(3,L)+ &
          e2*e5**2*A(3,L)*(1+3*QQ+A(3,L))-2*e2**2*e5*(QQ+2*A(3,L)))), &
          QQ*A(3,L)*(4*e2**3-e5**3*A(3,L)-2*e2**2*e5*(2+QQ+A(3,L))+ &
          e2*e5**2*(1+2*QQ+2*A(3,L)))*Q(4), &
          4*e2*QQ*A(3,L)*(e2**2+e5**2*A(3,L)-e2*e5*(QQ+2*A(3,L)))*Q(4),&
          4*e2*e3*QQ*A(3,L)*Q(4)**2, &
          2*(2*e2**2-2*e2*e5+e5**2)*QQ*A(3,L)*Q(4)*Q(5), &
          2*QQ*A(3,L)*(2*e2**4+e5**4*QQ*A(3,L)+2*e2**2*e5**2*A(3,L)*(2+2*QQ+A(3,L))- &
             e2*e5**3*A(3,L)*(1+3*QQ+A(3,L))-2*e2**3*e5*(QQ+3*A(3,L))),0._DP,0._DP, &
          QQ*A(3,L)*(4*e2**3-e5**3*QQ-2*e2**2*e5*(2+QQ+A(3,L))+ &
          e2*e5**2*(1+2*QQ+2*A(3,L)))*Q(4), &
          4*e2*QQ*A(3,L)*(e2-e5*(QQ+A(3,L)))*Q(4)**2, &
          4*QQ*A(3,L)*(e2**2+e4*e5*A(3,L))*Q(4)*Q(5),0._DP, &
          2*QQ*A(3,L)*(2*e2**3-e5**3*QQ*A(3,L)+e2*e5**2*A(3,L)* &
          (1+3*QQ+A(3,L))-2*e2**2*e5*(QQ+2*A(3,L)))*Q(4), &
          -(e4*QQ*A(3,L)*(2*e2-e5*(QQ+A(3,L)))*Q(4)**2), &
          -2*e4*QQ*A(3,L)*Q(4)**2*Q(5),4*QQ*A(3,L)*Q(4)**3*Q(5)/)
      CASE (7_I2B)
        QQ=Merge(A(2,L), A(1,L), G(1,L)==G(5,L))
        Q(1:14)=(/e3, e2 + QQ - e1*QQ, e2 + A(5,L) - e1*A(5,L), e2 - e5*(QQ + &
          A(5,L)), e2**2 - 2*e2*e5*(QQ + A(5,L)) + e5**2*(QQ + A(5,L)), 2*e2**2 + &
          e5**2*QQ - e2*e5*(1 + QQ + A(5,L)), -2*e2**2 - e5**2*QQ + e2*e5*(1 + QQ + &
          A(5,L)), 2*e2**3 - e5**3*QQ + e2*e5**2*(3*QQ + 2*A(5,L)) - e2**2*e5*(1 + 3*QQ &
          + 3*A(5,L)), -4*e2**4 + 4*e2*e5**3*QQ - e5**4*QQ + 4*e2**3*e5*(1 + QQ + &
          A(5,L)) - e2**2*e5**2*(1 + 7*QQ + 3*A(5,L)), 4*e2**4 - 4*e2*e5**3*QQ + &
          e5**4*QQ - 4*e2**3*e5*(1 + QQ + A(5,L)) + e2**2*e5**2*(1 + 7*QQ + 3*A(5,L)), &
          2*e2**3 - 2*e5**3*QQ*A(5,L) - 4*e2**2*e5*(QQ + A(5,L)) + e2*e5**2*(QQ + QQ**2 &
          + A(5,L) + 4*QQ*A(5,L) + A(5,L)**2), 2*e2**4 - 4*e2*e5**3*QQ*A(5,L) + &
          e5**4*QQ*A(5,L) - 4*e2**3*e5*(QQ + A(5,L)) + e2**2*e5**2*(QQ + QQ**2 + A(5,L) &
          + 6*QQ*A(5,L) + A(5,L)**2), 4*e2**4 - 8*e2**3*e5*(QQ + A(5,L)) + &
          e5**4*QQ*A(5,L)*(QQ + A(5,L)) - 4*e2*e5**3*QQ*A(5,L)*(1 + QQ + A(5,L)) + &
          e2**2*e5**2*(QQ + 3*QQ**2 + A(5,L) + 14*QQ*A(5,L) + 3*A(5,L)**2), 4*e2**4 + &
          e5**4*QQ*A(5,L) - 2*e2**3*e5*(1 + 3*QQ + 3*A(5,L)) - e2*e5**3*QQ*(1 + QQ + &
          5*A(5,L)) + e2**2*e5**2*(QQ**2 + A(5,L)*(3 + A(5,L)) + QQ*(5 + 6*A(5,L)))/)
        coe(L,:)=(/4*e2**2*Q(1)**2, 4*e2**2*A(5,L)*Q(1)**2, 2*e2*Q(1)*Q(6), &
          2*e2*Q(1)*Q(6), 4*e2**2*A(5,L)*Q(1)**2, 4*e2**2*A(5,L)*Q(1)**2, Q(10), &
          4*e2**2*Q(5), 4*e2**2*Q(1)*Q(4), 4*e2**2*Q(1)*Q(4), 2*e2*Q(8), &
          4*e2**2*A(5,L)*Q(1)**2, 2*e2*Q(8), 4*e2**2*A(5,L)*Q(1)**2, 4*e2**2*Q(4)**2, &
          4*e2**2*A(5,L)*Q(1)*Q(4), 4*e2**2*A(5,L)*Q(1)*Q(4), 2*Q(12), A(5,L)*Q(10), &
          4*e2**2*A(5,L)*Q(5), 4*e2*Q(1)*Q(2)*Q(3), 4*e2*Q(1)*Q(2)*Q(3), &
          2*e2*A(5,L)*Q(1)*Q(6), 2*e2*A(5,L)*Q(1)*Q(6), Q(6)**2, 2*e2*A(5,L)*Q(8), &
          2*e2*A(5,L)*Q(8), Q(14), 2*e2*A(5,L)*Q(1)*Q(6), 2*e2*A(5,L)*Q(1)*Q(6), Q(14), &
          2*e2*Q(11), 4*e2**2*A(5,L)*Q(1)*Q(4), 2*e2*Q(11), 2*e2*Q(4)*Q(6), &
          4*e2**2*A(5,L)*Q(1)*Q(4), 2*e2*Q(4)*Q(6), A(5,L)*Q(10), &
          2*e2*A(5,L)*Q(1)*Q(6), 2*e2*A(5,L)*Q(1)*Q(6), 4*e2**2*A(5,L)*Q(4)**2, &
          4*e2*Q(2)*Q(3)*Q(4), 4*e2*Q(2)*Q(3)*Q(4), 4*e2*A(5,L)*Q(1)*Q(2)*Q(3), &
          4*e2*A(5,L)*Q(1)*Q(2)*Q(3), 2*A(5,L)*Q(12), 2*e2*A(5,L)*Q(4)*Q(6), &
          2*e2*A(5,L)*Q(4)*Q(6), Q(13), A(5,L)*Q(14), A(5,L)*Q(14), A(5,L)*Q(6)**2, &
          2*e2*A(5,L)*Q(11), 4*e2*A(5,L)*Q(1)*Q(2)*Q(3), 4*e2*A(5,L)*Q(1)*Q(2)*Q(3), &
          2*e2*A(5,L)*Q(11), 2*Q(2)*Q(3)*Q(6), 2*Q(2)*Q(3)*Q(6), A(5,L)*Q(6)**2, &
          4*e2*A(5,L)*Q(2)*Q(3)*Q(4), 4*e2*A(5,L)*Q(2)*Q(3)*Q(4), 4*Q(2)**2*Q(3)**2, &
          A(5,L)*Q(13), 2*A(5,L)*Q(2)*Q(3)*Q(6), 2*A(5,L)*Q(2)*Q(3)*Q(6), &
          4*A(5,L)*Q(2)**2*Q(3)**2/)
      CASE (8_I2B)
        Q(1:20)=(/e2**4*A(5,L),e2**3*(e2+A(1,L)-e1*A(1,L)), &
          e2**2*(e2**2+e4*e5*A(1,L)), &
          e2*(e2**3-e5*(3*e2**2-3*e2*e5+e5**2)*A(1,L)), &
          (e2**2+e4*e5*A(1,L))**2,e2**2*(e2**2+e4*e5*A(1,L))*A(5,L), &
          e2**2*(e2+A(1,L)-e1*A(1,L))**2,e2**3*(e2+A(1,L)-e1*A(1,L))*A(5,L), &
          e2*(e2**3-e5*(3*e2**2-3*e2*e5+e5**2)*A(1,L))*A(5,L), &
          e2*(e2+A(1,L)-e1*A(1,L))*(e2**2+e4*e5*A(1,L)), &
          e2**4-4*e2**3*e5*A(1,L)+e5**4*A(1,L)**2+ &
          3*e2**2*e5**2*A(1,L)*(1+A(1,L))-e2*e5**3*A(1,L)*(1+3*A(1,L)), &
          (e2**2+e4*e5*A(1,L))**2*A(5,L), &
          (e2+A(1,L)-e1*A(1,L))**2*(e2**2+e4*e5*A(1,L)), &
          e2**2*(e2+A(1,L)-e1*A(1,L))**2*A(5,L), &
          e2*(e2+A(1,L)-e1*A(1,L))*(e2**2+e4*e5*A(1,L))*A(5,L), &
          (e2**4-4*e2**3*e5*A(1,L)+e5**4*A(1,L)**2+3*e2**2*e5**2*A(1,L)* &
          (1+A(1,L))-e2*e5**3*A(1,L)*(1+3*A(1,L)))*A(5,L), &
          e2*(e2+A(1,L)-e1*A(1,L))**3, &
          (e2+A(1,L)-e1*A(1,L))**2*(e2**2+e4*e5*A(1,L))*A(5,L), &
          e2*(e2+A(1,L)-e1*A(1,L))**3*A(5,L), &
          e2**4-4*e2**3*e5*A(1,L)+6*e2**2*e5**2*A(1,L)-4*e2*e5**3*A(1,L)+ &
          e5**4*A(1,L)/)
        coe(L,:)=(/e2**4,Q(1),Q(2),Q(2),Q(1),Q(1),Q(3),Q(20),Q(3), &
          Q(3),Q(4),Q(1),Q(4),Q(1),Q(5),Q(6),Q(6),Q(5),Q(6),A(5,L)*Q(20), &
          Q(7),Q(7),Q(8),Q(8),Q(7),Q(9),Q(9),Q(10),Q(8),Q(8),Q(10), &
          Q(11),Q(6),Q(11),Q(10),Q(6),Q(10),Q(6),Q(8),Q(8),Q(12),Q(13), &
          Q(13),Q(14),Q(14),Q(12),Q(15),Q(15),Q(13),Q(15),Q(15),Q(14), &
          Q(16),Q(14),Q(14),Q(16),Q(17),Q(17),Q(14),Q(18),Q(18), &
          (e2+A(1,L)-e1*A(1,L))**4,Q(18),Q(19),Q(19), &
          (e2+A(1,L)-e1*A(1,L))**4*A(5,L)/)
      CASE (9_I2B)
        Q(1:11)=(/e2**2*e3**2*A(1,L),e2**2*A(1,L)*(e2**2+e4*e5*A(1,L)), &
          e3**2*A(1,L)*(e2**2+e4*e5*A(3,L)), &
          e2**2*e3*A(1,L)*(e2-e5*(A(1,L)+A(3,L))), &
          e2*e3**2*A(1,L)*(e2+A(3,L)-e1*A(3,L)), &
          e2**2*e3*A(1,L)*(e2+A(1,L)-e1*A(1,L)), &
          e2**2*A(1,L)*(e2+A(1,L)-e1*A(1,L))**2, &
          e3**2*A(1,L)*(e2+A(3,L)-e1*A(3,L))**2, &
          e2*e3*A(1,L)*(e2+A(1,L)-e1*A(1,L))*(e2+A(3,L)-e1*A(3,L)), &
          e2+A(3,L)-e1*A(3,L),e2+A(1,L)-e1*A(1,L)/)
        coe(L,:)=(/e2**2*e3**2,Q(1),e2*e3**2*Q(10),e2**2*e3*Q(11), &
          Q(1),Q(1),e2**2*e3*(e2-e5*(A(1,L)+A(3,L))), &
          e2**2*(e2**2-2*e2*e5*(A(1,L)+A(3,L))+e5**2*(A(1,L)+A(3,L))), &
          e3**2*(e2**2+e4*e5*A(3,L)),e2**2*(e2**2+e4*e5*A(1,L)), &
          e2**2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L))),Q(1), &
          e2*e3*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L))),Q(1), &
          (e2**2+e4*e5*A(1,L))*(e2**2+e4*e5*A(3,L)),Q(2),Q(3), &
          e2**2*(e2-e5*(A(1,L)+A(3,L)))**2,Q(4), &
          e2**2*A(1,L)*(e2**2-2*e2*e5*(A(1,L)+A(3,L))+e5**2*(A(1,L)+A(3,L))), &
          e3**2*Q(10)**2,e2**2*Q(11)**2,Q(5),Q(6),e2*e3*Q(10)*Q(11), &
          e2**2*A(1,L)*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L))), &
          e2*e3*A(1,L)*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L))), &
          e2*e3*(e2-e5*(A(1,L)+A(3,L)))*Q(10),Q(5),Q(6), &
          e2**2*(e2-e5*(A(1,L)+A(3,L)))*Q(11), &
          e2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)))*Q(10),Q(3), &
          e2*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L)))*Q(11), &
          e3*(e2**2+e4*e5*A(3,L))*Q(11),Q(2), &
          e2*(e2**2+e4*e5*A(1,L))*Q(10),Q(4),Q(5),Q(6), &
          A(1,L)*(e2**2+e4*e5*A(1,L))*(e2**2+e4*e5*A(3,L)), &
          (e2**2+e4*e5*A(1,L))*Q(10)**2,(e2**2+e4*e5*A(3,L))*Q(11)**2, &
          Q(7),Q(8),e2**2*A(1,L)*(e2-e5*(A(1,L)+A(3,L)))**2, &
          e2*A(1,L)*(e2**2+e4*e5*A(1,L))*Q(10), &
          e3*A(1,L)*(e2**2+e4*e5*A(3,L))*Q(11), &
          e2*(e2-e5*(A(1,L)+A(3,L)))*Q(10)*Q(11), &
          e2*e3*A(1,L)*(e2-e5*(A(1,L)+A(3,L)))*Q(10), &
          e2**2*A(1,L)*(e2-e5*(A(1,L)+A(3,L)))*Q(11),Q(9), &
          e2*A(1,L)*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)))*Q(10), &
          Q(8),Q(7),e2*A(1,L)*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L)))*Q(11), &
          e3*Q(10)**2*Q(11),e2*Q(10)*Q(11)**2,Q(9), &
          A(1,L)*(e2**2+e4*e5*A(1,L))*Q(10)**2, &
          A(1,L)*(e2**2+e4*e5*A(3,L))*Q(11)**2,Q(10)**2*Q(11)**2, &
          e2*A(1,L)*(e2-e5*(A(1,L)+A(3,L)))*Q(10)*Q(11), &
          e3*A(1,L)*Q(10)**2*Q(11),e2*A(1,L)*Q(10)*Q(11)**2, &
          A(1,L)*Q(10)**2*Q(11)**2/)
      CASE (10_I2B)
        Q(1:11)=(/e2**2*e3**2*A(3,L),e3**2*(e2**2+e4*e5*A(1,L))*A(3,L), &
          e2**2*A(3,L)*(e2**2+e4*e5*A(3,L)), &
          e2**2*e3*A(3,L)*(e2-e5*(A(1,L)+A(3,L))), &
          e2**2*e3*A(3,L)*(e2+A(3,L)-e1*A(3,L)), &
          e2*e3**2*(e2+A(1,L)-e1*A(1,L))*A(3,L), &
          e3**2*(e2+A(1,L)-e1*A(1,L))**2*A(3,L), &
          e2**2*A(3,L)*(e2+A(3,L)-e1*A(3,L))**2, &
          e2*e3*(e2+A(1,L)-e1*A(1,L))*A(3,L)*(e2+A(3,L)-e1*A(3,L)), &
          e2+A(3,L)-e1*A(3,L),e2+A(1,L)-e1*A(1,L)/)
        coe(L,:)=(/e2**2*e3**2,Q(1),e2**2*e3*Q(10),e2*e3**2*Q(11),Q(1), &
          Q(1),e2**2*e3*(e2-e5*(A(1,L)+A(3,L))), &
          e2**2*(e2**2-2*e2*e5*(A(1,L)+A(3,L))+e5**2*(A(1,L)+A(3,L))), &
          e2**2*(e2**2+e4*e5*A(3,L)),e3**2*(e2**2+e4*e5*A(1,L)), &
          e2*e3*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L))),Q(1), &
          e2**2*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L))),Q(1), &
          (e2**2+e4*e5*A(1,L))*(e2**2+e4*e5*A(3,L)),Q(2),Q(3), &
          e2**2*(e2-e5*(A(1,L)+A(3,L)))**2,Q(4), &
          e2**2*A(3,L)*(e2**2-2*e2*e5*(A(1,L)+A(3,L))+e5**2*(A(1,L)+A(3,L))), &
          e2**2*Q(10)**2,e3**2*Q(11)**2,Q(5),Q(6),e2*e3*Q(10)*Q(11), &
          e2*e3*A(3,L)*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L))), &
          e2**2*A(3,L)*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L))), &
          e2**2*(e2-e5*(A(1,L)+A(3,L)))*Q(10),Q(5),Q(6), &
          e2*e3*(e2-e5*(A(1,L)+A(3,L)))*Q(11), &
          e2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)))*Q(10),Q(3), &
          e2*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L)))*Q(11), &
          e2*(e2**2+e4*e5*A(3,L))*Q(11),Q(2), &
          e3*(e2**2+e4*e5*A(1,L))*Q(10),Q(4),Q(5),Q(6), &
          (e2**2+e4*e5*A(1,L))*A(3,L)*(e2**2+e4*e5*A(3,L)), &
          (e2**2+e4*e5*A(1,L))*Q(10)**2,(e2**2+e4*e5*A(3,L))*Q(11)**2, &
          Q(7),Q(8),e2**2*A(3,L)*(e2-e5*(A(1,L)+A(3,L)))**2, &
          e3*(e2**2+e4*e5*A(1,L))*A(3,L)*Q(10), &
          e2*A(3,L)*(e2**2+e4*e5*A(3,L))*Q(11), &
          e2*(e2-e5*(A(1,L)+A(3,L)))*Q(10)*Q(11), &
          e2**2*A(3,L)*(e2-e5*(A(1,L)+A(3,L)))*Q(10), &
          e2*e3*A(3,L)*(e2-e5*(A(1,L)+A(3,L)))*Q(11),Q(9), &
          e2*A(3,L)*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)))*Q(10), &
          Q(8),Q(7),e2*A(3,L)*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L)))*Q(11),&
          e2*Q(10)**2*Q(11),e3*Q(10)*Q(11)**2,Q(9), &
          (e2**2+e4*e5*A(1,L))*A(3,L)*Q(10)**2, &
          A(3,L)*(e2**2+e4*e5*A(3,L))*Q(11)**2,Q(10)**2*Q(11)**2, &
          e2*A(3,L)*(e2-e5*(A(1,L)+A(3,L)))*Q(10)*Q(11), &
          e2*A(3,L)*Q(10)**2*Q(11),e3*A(3,L)*Q(10)*Q(11)**2, &
          A(3,L)*Q(10)**2*Q(11)**2/)
      CASE (11_I2B)
        Q(1:13)=(/2*e2**3*e3*A(5,L), &
          2*e2*(e2**2+e4*e5*A(1,L))*(e2-e5*(A(1,L)+A(5,L))), &
          2*e2*e3*(e2**2+e4*e5*A(1,L))*A(5,L), &
          2*e2**3*A(5,L)*(e2-e5*(A(1,L)+A(5,L))), &
          e2*A(5,L)*(2*e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)- &
          e2**2*e5*(1+3*A(1,L)+A(5,L))), &
          e2**2*A(5,L)*(2*e2**2+e5**2*A(1,L)-e2*e5*(1+A(1,L)+A(5,L))),&
          2*e2**2*e3*(e2+A(1,L)-e1*A(1,L))*A(5,L), &
          2*e2*(e2**2+e4*e5*A(1,L))*A(5,L)*(e2-e5*(A(1,L)+A(5,L))), &
          2*e2*e3*(e2+A(1,L)-e1*A(1,L))**2*A(5,L), &
          2*e2**2*(e2+A(1,L)-e1*A(1,L))*A(5,L)*(e2+A(5,L)-e1*A(5,L)), &
          e2*(e2+A(1,L)-e1*A(1,L))*A(5,L)*(2*e2**2+e5**2*A(1,L)- &
          e2*e5*(1+A(1,L)+A(5,L))),e2+A(5,L)-e1*A(5,L), &
          e2+A(1,L)-e1*A(1,L)/)
        coe(L,:)=(/2*e2**3*e3,Q(1),e2**2*(2*e2**2+e5**2*A(1,L)- &
          e2*e5*(1+A(1,L)+A(5,L))),2*e2**2*e3*Q(13),Q(1),Q(1), &
          e2*(2*e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(1+3*A(1,L)+A(5,L))), &
          2*e2*(e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(3*A(1,L)+A(5,L))), &
          2*e2**3*(e2-e5*(A(1,L)+A(5,L))),2*e2*e3*(e2**2+e4*e5*A(1,L)), &
          2*e2**4+7*e2**2*e5**2*A(1,L)-4*e2*e5**3*A(1,L)+e5**4*A(1,L)- &
          e2**3*e5*(1+5*A(1,L)+A(5,L)),Q(1), &
          2*e2**2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(5,L))),Q(1),Q(2),Q(3),Q(4), &
          Q(2),Q(5),2*e2*A(5,L)*(e2**3+3*e2*e5**2*A(1,L)- &
          e5**3*A(1,L)-e2**2*e5*(3*A(1,L)+A(5,L))),2*e2**2*Q(12)*Q(13), &
          2*e2*e3*Q(13)**2,Q(6),Q(7),e2*(2*e2**2+ &
          e5**2*A(1,L)-e2*e5*(1+A(1,L)+A(5,L)))*Q(13), &
          A(5,L)*(2*e2**4+7*e2**2*e5**2*A(1,L)-4*e2*e5**3*A(1,L)+ &
          e5**4*A(1,L)-e2**3*e5*(1+5*A(1,L)+A(5,L))), &
          2*e2**2*A(5,L)*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(5,L))), &
          e2*(2*e2**3-e5**3*A(1,L)*A(5,L)-2*e2**2*e5*(2*A(1,L)+ &
          A(5,L))+e2*e5**2*A(1,L)*(1+A(1,L)+3*A(5,L))),Q(6), &
          Q(7),(2*e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(1+ &
          3*A(1,L)+A(5,L)))*Q(13), &
          2*e2**4+e5**4*A(1,L)*A(5,L)-2*e2**3*e5*(3*A(1,L)+A(5,L))+ &
          2*e2**2*e5**2*A(1,L)*(2+A(1,L)+2*A(5,L))- &
          e2*e5**3*A(1,L)*(1+A(1,L)+3*A(5,L)),Q(4), &
          2*e2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(5,L)))*Q(13), &
          2*e2**2*(e2-e5*(A(1,L)+A(5,L)))*Q(13),Q(3), &
          (e2**2+e4*e5*A(1,L))*(2*e2**2+e5**2*A(1,L)-e2*e5*(1+A(1,L)+A(5,L))), &
          Q(5),Q(6),Q(7),Q(8),2*(e2**2+e4*e5*A(1,L))*Q(12)*Q(13), &
          2*e2*(e2-e5*(A(1,L)+A(5,L)))*Q(13)**2,Q(9),Q(10),Q(8), &
          (e2**2+e4*e5*A(1,L))*A(5,L)*(2*e2**2+e5**2*A(1,L)-e2*e5*(1+A(1,L)+A(5,L))), &
          2*e2**2*A(5,L)*(e2-e5*(A(1,L)+A(5,L)))*Q(13), &
          (2*e2**3-e5**3*A(1,L)*A(5,L)-2*e2**2*e5*(2*A(1,L)+ &
          A(5,L))+e2*e5**2*A(1,L)*(1+A(1,L)+3*A(5,L)))*Q(13), &
          e2*A(5,L)*(2*e2**3-e5**3*A(1,L)*A(5,L)-2*e2**2*e5*(2*A(1,L)+ &
          A(5,L))+e2*e5**2*A(1,L)*(1+A(1,L)+3*A(5,L))), &
          A(5,L)*(2*e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)- &
          e2**2*e5*(1+3*A(1,L)+A(5,L)))*Q(13),Q(11), &
          A(5,L)*(2*e2**4+e5**4*A(1,L)*A(5,L)-2*e2**3*e5*(3*A(1,L)+A(5,L))+ &
          2*e2**2*e5**2*A(1,L)*(2+A(1,L)+2*A(5,L))- &
          e2*e5**3*A(1,L)*(1+A(1,L)+3*A(5,L))),Q(10),Q(9), &
          2*e2*A(5,L)*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(5,L)))*Q(13), &
          2*e2*Q(12)*Q(13)**2, &
          (2*e2**2+e5**2*A(1,L)-e2*e5*(1+A(1,L)+A(5,L)))*Q(13)**2,Q(11), &
          2*(e2**2+e4*e5*A(1,L))*A(5,L)*Q(12)*Q(13), &
          2*e2*A(5,L)*(e2-e5*(A(1,L)+A(5,L)))*Q(13)**2,2*Q(12)*Q(13)**3, &
          A(5,L)*(2*e2**3-e5**3*A(1,L)*A(5,L)-2*e2**2*e5*(2*A(1,L)+ &
          A(5,L))+e2*e5**2*A(1,L)*(1+A(1,L)+3*A(5,L)))*Q(13), &
          2*e2*A(5,L)*Q(12)*Q(13)**2,A(5,L)*(2*e2**2+e5**2*A(1,L)- &
          e2*e5*(1+A(1,L)+A(5,L)))*Q(13)**2,2*A(5,L)*Q(12)*Q(13)**3/)
      CASE (12_I2B)
        Q(1:6)=(/2*e2**2*e3**2,e2*e3*(2*e2**2-2*e2*e5+e5**2), &
          -(e2*e3*e4*(e2+A(3,L)-e1*A(3,L))),-(e2*e3*e4*(e2+A(1,L)-e1*A(1,L))),&
          e2+A(3,L)-e1*A(3,L),e2+A(1,L)-e1*A(1,L)/)
        coe(L,:)=(/0._DP,Q(1),0._DP,0._DP,Q(2),Q(2),0._DP,0._DP,0._DP, &
          0._DP,0._DP,2*e2**4-4*e2**3*e5+6*e2**2*e5**2-4*e2*e5**3+e5**4,&
          0._DP,Q(1),0._DP,2*e2*e3*(e2**2+e4*e5*A(1,L)), &
          2*e2*e3*(e2**2+e4*e5*A(3,L)),0._DP, &
          e2*(2*e2**2-2*e2*e5+e5**2)*(e2-e5*(A(1,L)+A(3,L))), &
          2*e2**2*(e2**2-2*e2*e5*(A(1,L)+A(3,L))+e5**2*(A(1,L)+A(3,L))), &
          0._DP,0._DP,Q(3),Q(4),0._DP, &
          -(e2*e4*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)))), &
          -(e2*e4*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L)))), &
          0._DP,-(e4*(e2-e1*e2+e2**2+e5**2)*Q(5)), &
          -(e4*(e2-e1*e2+e2**2+e5**2)*Q(6)),0._DP,0._DP, &
          (2*e2**2-2*e2*e5+e5**2)*(e2**2+e4*e5*A(3,L)),0._DP,0._DP, &
          (2*e2**2-2*e2*e5+e5**2)*(e2**2+e4*e5*A(1,L)),0._DP, &
          2*e2**2*e3*(e2-e5*(A(1,L)+A(3,L))),Q(3),Q(4), &
          2*(e2**2+e4*e5*A(1,L))*(e2**2+e4*e5*A(3,L)),0._DP,0._DP, &
          2*e2*e3*Q(6)**2,2*e2*e3*Q(5)**2, &
          2*e2**2*(e2-e5*(A(1,L)+A(3,L)))**2,-(e4*(e2**2+e4*e5*A(1,L))*Q(5)),&
          -(e4*(e2**2+e4*e5*A(3,L))*Q(6)),0._DP, &
          -(e2*e4*(e2-e5*(A(1,L)+A(3,L)))*Q(5)), &
          -(e2*e4*(e2-e5*(A(1,L)+A(3,L)))*Q(6)), &
          (2*e2**2-2*e2*e5+e5**2)*Q(5)*Q(6), &
          2*e2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)))*Q(5), &
          (2*e2**2-2*e2*e5+e5**2)*Q(5)**2,(2*e2**2-2*e2*e5+e5**2)*Q(6)**2, &
          2*e2*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L)))*Q(6),&
          0._DP,0._DP,2*e2*e3*Q(5)*Q(6),2*(e2**2+e4*e5*A(1,L))*Q(5)**2, &
          2*(e2**2+e4*e5*A(3,L))*Q(6)**2,0._DP, &
          2*e2*(e2-e5*(A(1,L)+A(3,L)))*Q(5)*Q(6),-(e4*Q(5)**2*Q(6)), &
          -(e4*Q(5)*Q(6)**2),2*Q(5)**2*Q(6)**2/)
      CASE (13_I2B)
        Q(1:13)=(/2*e2**3*e3*A(5,L),2*e2*(e2**2+e4*e5*A(3,L))*(e2-e5*(A(3,L)+A(5,L))), &
        2*e2**3*A(5,L)*(e2-e5*(A(3,L)+A(5,L))),2*e2*e3*(e2**2+e4*e5*A(3,L))*A(5,L), &
        e2*A(5,L)*(2*e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(1+3*A(3,L)+A(5,L))), &
        2*e2**2*e3*(e2+A(3,L)-e1*A(3,L))*A(5,L),e2**2*A(5,L)*(2*e2**2+e5**2*A(3,L)-e2*e5*(1+A(3,L)+A(5,L))), &
        2*e2*(e2**2+e4*e5*A(3,L))*A(5,L)*(e2-e5*(A(3,L)+A(5,L))), &
        2*e2**2*(e2+A(3,L)-e1*A(3,L))*A(5,L)*(e2+A(5,L)-e1*A(5,L)),2*e2*e3*(e2+A(3,L)-e1*A(3,L))**2*A(5,L), &
        e2*(e2+A(3,L)-e1*A(3,L))*A(5,L)*(2*e2**2+e5**2*A(3,L)-e2*e5*(1+A(3,L)+A(5,L))),e2+A(3,L)-e1*A(3,L), &
        e2+A(5,L)-e1*A(5,L)/)
        coe(L,:)=(/2*e2**3*e3,Q(1),2*e2**2*e3*Q(12), &
        e2**2*(2*e2**2+e5**2*A(3,L)-e2*e5*(1+A(3,L)+A(5,L))),Q(1),Q(1), &
        e2*(2*e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(1+3*A(3,L)+A(5,L))), &
        2*e2*(e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(3*A(3,L)+A(5,L))),2*e2*e3*(e2**2+e4*e5*A(3,L)), &
        2*e2**3*(e2-e5*(A(3,L)+A(5,L))),2*e2**2*(e2**2+e5**2*A(3,L)-e2*e5*(2*A(3,L)+A(5,L))),Q(1), &
        2*e2**4+7*e2**2*e5**2*A(3,L)-4*e2*e5**3*A(3,L)+e5**4*A(3,L)-e2**3*e5*(1+5*A(3,L)+A(5,L)), &
        Q(1),Q(2),Q(3),Q(4),Q(2),Q(5), &
        2*e2*A(5,L)*(e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(3*A(3,L)+A(5,L))),2*e2*e3*Q(12)**2, &
        2*e2**2*Q(12)*Q(13),Q(6),Q(7),e2*(2*e2**2+e5**2*A(3,L)-e2*e5*(1+A(3,L)+A(5,L)))*Q(12), &
        2*e2**2*A(5,L)*(e2**2+e5**2*A(3,L)-e2*e5*(2*A(3,L)+A(5,L))), &
        A(5,L)*(2*e2**4+7*e2**2*e5**2*A(3,L)-4*e2*e5**3*A(3,L)+e5**4*A(3,L)-e2**3*e5*(1+5*A(3,L)+A(5,L))), &
        (2*e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(1+3*A(3,L)+A(5,L)))*Q(12),Q(6),Q(7), &
        e2*(2*e2**3-e5**3*A(3,L)*A(5,L)-2*e2**2*e5*(2*A(3,L)+A(5,L))+e2*e5**2*A(3,L)*(1+A(3,L)+3*A(5,L))), &
        2*e2*(e2**2+e5**2*A(3,L)-e2*e5*(2*A(3,L)+A(5,L)))*Q(12),Q(4), &
        2*e2**4+e5**4*A(3,L)*A(5,L)-2*e2**3*e5*(3*A(3,L)+A(5,L))+2*e2**2*e5**2*A(3,L)*(2+A(3,L)+2*A(5,L))- &
         e2*e5**3*A(3,L)*(1+A(3,L)+3*A(5,L)),(e2**2+e4*e5*A(3,L))*(2*e2**2+e5**2*A(3,L)-e2*e5*(1+A(3,L)+A(5,L))), &
        Q(3),2*e2**2*(e2-e5*(A(3,L)+A(5,L)))*Q(12),Q(5),Q(6),Q(7),Q(8),2*e2*(e2-e5*(A(3,L)+A(5,L)))*Q(12)**2, &
        2*(e2**2+e4*e5*A(3,L))*Q(12)*Q(13),Q(9),Q(10),Q(8),2*e2**2*A(5,L)*(e2-e5*(A(3,L)+A(5,L)))*Q(12), &
        (e2**2+e4*e5*A(3,L))*A(5,L)*(2*e2**2+e5**2*A(3,L)-e2*e5*(1+A(3,L)+A(5,L))), &
        (2*e2**3-e5**3*A(3,L)*A(5,L)-2*e2**2*e5*(2*A(3,L)+A(5,L))+e2*e5**2*A(3,L)*(1+A(3,L)+3*A(5,L)))*Q(12), &
        A(5,L)*(2*e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(1+3*A(3,L)+A(5,L)))*Q(12), &
        e2*A(5,L)*(2*e2**3-e5**3*A(3,L)*A(5,L)-2*e2**2*e5*(2*A(3,L)+A(5,L))+e2*e5**2*A(3,L)*(1+A(3,L)+3*A(5,L))), &
        Q(11),2*e2*A(5,L)*(e2**2+e5**2*A(3,L)-e2*e5*(2*A(3,L)+A(5,L)))*Q(12),Q(10),Q(9), &
        A(5,L)*(2*e2**4+e5**4*A(3,L)*A(5,L)-2*e2**3*e5*(3*A(3,L)+A(5,L))+ &
           2*e2**2*e5**2*A(3,L)*(2+A(3,L)+2*A(5,L))-e2*e5**3*A(3,L)*(1+A(3,L)+3*A(5,L))), &
        (2*e2**2+e5**2*A(3,L)-e2*e5*(1+A(3,L)+A(5,L)))*Q(12)**2,2*e2*Q(12)**2*Q(13),Q(11), &
        2*e2*A(5,L)*(e2-e5*(A(3,L)+A(5,L)))*Q(12)**2,2*(e2**2+e4*e5*A(3,L))*A(5,L)*Q(12)*Q(13),2*Q(12)**3*Q(13), &
        A(5,L)*(2*e2**3-e5**3*A(3,L)*A(5,L)-2*e2**2*e5*(2*A(3,L)+A(5,L))+e2*e5**2*A(3,L)*(1+A(3,L)+3*A(5,L)))* &
         Q(12),A(5,L)*(2*e2**2+e5**2*A(3,L)-e2*e5*(1+A(3,L)+A(5,L)))*Q(12)**2,2*e2*A(5,L)*Q(12)**2*Q(13), &
        2*A(5,L)*Q(12)**3*Q(13)/)
      CASE (14_I2B)
        Q(1:14)=(/8*e2**2*e3**2,4*e2*e3*(2*e2**2-2*e2*e5+e5**2), &
        4*e2*(2*e2**2-2*e2*e5+e5**2)*(e2-e5*(A(1,L)+A(2,L))), &
        -2*e2*e3*e4*(2*e2-e5*(A(1,L)+A(2,L))), &
        -2*e2*e4*(2*e2**2-3*e2*e5*(A(1,L)+A(2,L))+e5**2*(A(1,L)+A(2,L))), &
        8*e2**2*e3*(e2-e5*(A(1,L)+A(2,L))),-(e4*(2*e2**2-2*e2*e5+e5**2)*(2*e2-e5*(A(1,L)+A(2,L)))), &
        4*(2*e2**2-2*e2*e5+e5**2)*(e2+A(1,L)-e1*A(1,L))*(e2+A(2,L)-e1*A(2,L)), &
        -2*e2*e4*(e2-e5*(A(1,L)+A(2,L)))*(2*e2-e5*(A(1,L)+A(2,L))), &
        -(e4*(4*e2**3-2*e5**3*A(1,L)*A(2,L)-6*e2**2*e5*(A(1,L)+A(2,L))+ &
             e2*e5**2*(A(1,L)+A(1,L)**2+A(2,L)+6*A(1,L)*A(2,L)+A(2,L)**2))), &
        4*e2*(2*e2**3-2*e5**3*A(1,L)*A(2,L)-4*e2**2*e5*(A(1,L)+A(2,L))+ &
           e2*e5**2*(A(1,L)+A(1,L)**2+A(2,L)+4*A(1,L)*A(2,L)+A(2,L)**2)), &
        8*e2*e3*(e2+A(1,L)-e1*A(1,L))*(e2+A(2,L)-e1*A(2,L)), &
        8*e2*(e2+A(1,L)-e1*A(1,L))*(e2+A(2,L)-e1*A(2,L))*(e2-e5*(A(1,L)+A(2,L))), &
        -2*e4*(e2+A(1,L)-e1*A(1,L))*(e2+A(2,L)-e1*A(2,L))*(2*e2-e5*(A(1,L)+A(2,L)))/)
        coe(L,:)=(/0._DP,Q(1),0._DP,0._DP,Q(2),Q(2),0._DP,0._DP, &
          0._DP,0._DP,0._DP,Q(1),0._DP,2*(2*e2**2-2*e2*e5+e5**2)**2,&
          0._DP,Q(3),Q(3),0._DP, &
        2*e2*(4*e2**3-2*e5**3*(A(1,L)+A(2,L))-4*e2**2*e5*(1+A(1,L)+A(2,L))+e2*e5**2*(1+5*A(1,L)+5*A(2,L))), &
        8*e2**2*(e2**2-2*e2*e5*(A(1,L)+A(2,L))+e5**2*(A(1,L)+A(2,L))),0._DP,0._DP,&
        Q(4),Q(4),0._DP,Q(5),Q(5),0._DP,Q(4),Q(4),0._DP,0._DP,Q(6),0._DP,0._DP,Q(6),0._DP, &
        8*e2**4-4*e2*e5**3*(A(1,L)+A(2,L))+e5**4*(A(1,L)+A(2,L))-8*e2**3*e5*(1+A(1,L)+A(2,L))+ &
         2*e2**2*e5**2*(1+5*A(1,L)+5*A(2,L)),Q(7),Q(7), &
         8*e2**2*(e2-e5*(A(1,L)+A(2,L)))**2,0._DP,0._DP,Q(8),Q(8), &
        4*(2*e2**4-4*e2*e5**3*A(1,L)*A(2,L)+e5**4*A(1,L)*A(2,L)-4*e2**3*e5*(A(1,L)+A(2,L))+ &
           e2**2*e5**2*(A(1,L)+A(1,L)**2+A(2,L)+6*A(1,L)*A(2,L)+A(2,L)**2)),Q(9),Q(9),0._DP,Q(10),Q(10), &
        2*(2*e2**2+e5**2*A(1,L)-e2*e5*(1+A(1,L)+A(2,L)))*(2*e2**2+e5**2*A(2,L)-e2*e5*(1+A(1,L)+A(2,L))),Q(11), &
        Q(12),Q(12),Q(11),0._DP,0._DP,8*e2**4-8*e2**3*e5*(1+A(1,L)+A(2,L))- &
         2*e2*e5**3*(A(1,L)+A(2,L))*(1+A(1,L)+A(2,L))+e5**4*(A(1,L)**2+A(2,L)**2)+ &
         2*e2**2*e5**2*(1+A(1,L)**2+2*A(1,L)*(2+A(2,L))+A(2,L)*(4+A(2,L))),Q(13),Q(13),0._DP, &
        2*(4*e2**4-8*e2**3*e5*(A(1,L)+A(2,L))+e5**4*A(1,L)*A(2,L)*(A(1,L)+A(2,L))- &
           4*e2*e5**3*A(1,L)*A(2,L)*(1+A(1,L)+A(2,L))+ &
           e2**2*e5**2*(A(1,L)+3*A(1,L)**2+A(2,L)+14*A(1,L)*A(2,L)+3*A(2,L)**2)),Q(14),Q(14), &
        8*(e2+A(1,L)-e1*A(1,L))**2*(e2+A(2,L)-e1*A(2,L))**2/)
      CASE (15_I2B)
        Q(1:11)=(/e2**4*A(5,L),e2**2*(e2**2+e4*e5*A(1,L))*A(5,L),e2**2*(e2**2+e4*e5*A(3,L))*A(5,L), &
        e2**3*(e2-e5*(A(1,L)+A(3,L)))*A(5,L), &
        e2**3*(e2+A(3,L)-e1*A(3,L))*A(5,L),e2**3*(e2+A(1,L)-e1*A(1,L))*A(5,L), &
        e2**2*(e2+A(1,L)-e1*A(1,L))**2*A(5,L),e2**2*(e2+A(3,L)-e1*A(3,L))**2*A(5,L), &
        e2**2*(e2+A(1,L)-e1*A(1,L))*(e2+A(3,L)-e1*A(3,L))*A(5,L),e2+A(3,L)-e1*A(3,L),e2+A(1,L)-e1*A(1,L)/)
        coe(L,:)=(/e2**4,Q(1),e2**3*Q(10),e2**3*Q(11),Q(1),Q(1),e2**3*(e2-e5*(A(1,L)+A(3,L))), &
        e2**2*(e2**2-2*e2*e5*(A(1,L)+A(3,L))+e5**2*(A(1,L)+A(3,L))),e2**2*(e2**2+e4*e5*A(3,L)), &
        e2**2*(e2**2+e4*e5*A(1,L)),e2**2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L))),Q(1), &
        e2**2*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L))),Q(1),&
        (e2**2+e4*e5*A(1,L))*(e2**2+e4*e5*A(3,L)),Q(2),Q(3), &
        e2**2*(e2-e5*(A(1,L)+A(3,L)))**2,Q(4), &
        e2**2*(e2**2-2*e2*e5*(A(1,L)+A(3,L))+e5**2*(A(1,L)+A(3,L)))*A(5,L), &
        e2**2*Q(10)**2,e2**2*Q(11)**2,Q(5),Q(6),e2**2*Q(10)*Q(11), &
        e2**2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)))*A(5,L), &
        e2**2*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L)))*A(5,L), &
        e2**2*(e2-e5*(A(1,L)+A(3,L)))*Q(10),Q(5),Q(6), &
        e2**2*(e2-e5*(A(1,L)+A(3,L)))*Q(11),e2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)))*Q(10),Q(3), &
        e2*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L)))*Q(11),e2*(e2**2+e4*e5*A(3,L))*Q(11),Q(2), &
        e2*(e2**2+e4*e5*A(1,L))*Q(10),Q(4),Q(5),Q(6),(e2**2+e4*e5*A(1,L))*(e2**2+e4*e5*A(3,L))*A(5,L), &
        (e2**2+e4*e5*A(1,L))*Q(10)**2,(e2**2+e4*e5*A(3,L))*Q(11)**2, &
        Q(7),Q(8),e2**2*(e2-e5*(A(1,L)+A(3,L)))**2*A(5,L), &
        e2*(e2**2+e4*e5*A(1,L))*A(5,L)*Q(10),e2*(e2**2+e4*e5*A(3,L))*A(5,L)*Q(11), &
        e2*(e2-e5*(A(1,L)+A(3,L)))*Q(10)*Q(11),e2**2*(e2-e5*(A(1,L)+A(3,L)))*A(5,L)*Q(10), &
        e2**2*(e2-e5*(A(1,L)+A(3,L)))*A(5,L)*Q(11),Q(9), &
        e2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)))*A(5,L)*Q(10), &
        Q(8),Q(7),e2*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L)))*A(5,L)*Q(11), &
        e2*Q(10)**2*Q(11),e2*Q(10)*Q(11)**2,Q(9), &
        (e2**2+e4*e5*A(1,L))*A(5,L)*Q(10)**2,(e2**2+e4*e5*A(3,L))*A(5,L)*Q(11)**2,Q(10)**2*Q(11)**2, &
        e2*(e2-e5*(A(1,L)+A(3,L)))*A(5,L)*Q(10)*Q(11),e2*A(5,L)*Q(10)**2*Q(11),e2*A(5,L)*Q(10)*Q(11)**2, &
        A(5,L)*Q(10)**2*Q(11)**2/)
      CASE (16_I2B)
        Q(1:10)=(/4*e2**3*e3,2*e2**2*(2*e2**2-2*e2*e5+e5**2),4*e2**3*(e2-e5*(A(3,L)+A(4,L))), &
        -(e2**2*e4*(2*e2-e5*(2*A(1,L)+A(3,L)+A(4,L)))),-(e2**2*e4*(2*e2-e5*(A(3,L)+A(4,L)))), &
        4*e2**2*e3*(e2+A(1,L)-e1*A(1,L)),4*e2**2*(e2+A(3,L)-e1*A(3,L))*(e2+A(4,L)-e1*A(4,L)), &
        -(e2*e4*(e2+A(1,L)-e1*A(1,L))*(2*e2-e5*(A(3,L)+A(4,L)))),e2+A(3,L)-e1*A(3,L),e2+A(4,L)-e1*A(4,L)/)
        coe(L,:)=(/0._DP,Q(1),0._DP,0._DP,Q(2),Q(1),0._DP,0._DP,0._DP,0._DP,0._DP,&
         Q(1),0._DP,Q(2),0._DP,2*(2*e2**2-2*e2*e5+e5**2)*(e2**2+e4*e5*A(1,L)),Q(3),0._DP,Q(4), &
        4*e2**2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)+A(4,L))),0._DP,0._DP,Q(5),Q(6),0._DP, &
        -(e2*e4*(2*e2**2+2*e5**2*A(1,L)-e2*e5*(4*A(1,L)+A(3,L)+A(4,L)))),4*e2**3*(e2-e5*(A(1,L)+A(3,L)+A(4,L))), &
        0._DP,Q(5),Q(6),0._DP,0._DP,Q(3),0._DP,0._DP,4*e2*e3*(e2**2+e4*e5*A(1,L)),0._DP,Q(4),Q(5), &
        2*e2*(2*e2**2-2*e2*e5+e5**2)*(e2+A(1,L)-e1*A(1,L)), &
        4*e2*(e2**2+e4*e5*A(1,L))*(e2-e5*(A(3,L)+A(4,L))),0._DP,0._DP, &
        2*(2*e2**2-2*e2*e5+e5**2)*(e2+A(1,L)-e1*A(1,L))**2,Q(7), &
        4*e2**2*(e2-e5*(A(1,L)+A(3,L)))*(e2-e5*(A(1,L)+A(4,L))), &
        -(e4*(e2**2+e4*e5*A(1,L))*(2*e2-e5*(A(3,L)+A(4,L)))), &
        4*e2**2*(e2+A(1,L)-e1*A(1,L))*(e2-e5*(A(3,L)+A(4,L))),0._DP, &
        2*e2**2*(2*e2**2-2*e2*e5*(A(1,L)+A(3,L)+A(4,L))+e5**2*(2*A(3,L)*A(4,L)+A(1,L)*(A(3,L)+A(4,L)))), &
        -(e2*e4*(e2+A(1,L)-e1*A(1,L))*(2*e2-e5*(2*A(1,L)+A(3,L)+A(4,L)))),Q(8), &
        2*e2*(2*e2**3-e5**3*A(1,L)*(A(3,L)+A(4,L))-2*e2**2*e5*(2*A(1,L)+A(3,L)+A(4,L))+ &
           2*e2*e5**2*(A(3,L)*A(4,L)+A(1,L)*(1+A(3,L)+A(4,L)))),Q(7),4*e2*e3*(e2+A(1,L)-e1*A(1,L))**2, &
        4*e2**2*(e2+A(1,L)-e1*A(1,L))*(e2-e5*(A(1,L)+A(3,L)+A(4,L))),&
        0._DP,0._DP,Q(8),4*(e2**2+e4*e5*A(1,L))*Q(9)*Q(10), &
        4*e2*(e2+A(1,L)-e1*A(1,L))**2*(e2-e5*(A(3,L)+A(4,L))),0._DP, &
        2*e2*(e2+A(1,L)-e1*A(1,L))*(2*e2**2-2*e2*e5*(A(1,L)+A(3,L)+A(4,L))+ &
           e5**2*(2*A(3,L)*A(4,L)+A(1,L)*(A(3,L)+A(4,L)))),4*e2*(e2+A(1,L)-e1*A(1,L))*Q(9)*Q(10), &
        -(e4*(e2+A(1,L)-e1*A(1,L))**2*(2*e2-e5*(A(3,L)+A(4,L)))),4*(e2+A(1,L)-e1*A(1,L))**2*Q(9)*Q(10)/)
      CASE (17_I2B)
        Q(1:11)=(/4*e2**3*e3,2*e2**2*(2*e2**2-2*e2*e5+e5**2),4*e2**3*(e2-e5*(A(1,L)+A(2,L))), &
        -(e2**2*e4*(2*e2-e5*(A(1,L)+A(2,L)+2*A(3,L)))),4*e2**2*e3*(e2+A(3,L)-e1*A(3,L)), &
        -(e2**2*e4*(2*e2-e5*(A(1,L)+A(2,L)))),4*e2**2*(e2+A(1,L)-e1*A(1,L))*(e2+A(2,L)-e1*A(2,L)), &
        -(e2*e4*(2*e2-e5*(A(1,L)+A(2,L)))*(e2+A(3,L)-e1*A(3,L))),e2+A(3,L)-e1*A(3,L),e2+A(2,L)-e1*A(2,L), &
        e2+A(1,L)-e1*A(1,L)/)
        coe(L,:)=(/0._DP,Q(1),0._DP,0._DP,Q(1),Q(2),0._DP,0._DP,0._DP,&
        0._DP,0._DP,Q(1),0._DP,Q(2),0._DP,Q(3),2*(2*e2**2-2*e2*e5+e5**2)*(e2**2+e4*e5*A(3,L)),0._DP,Q(4),&
        4*e2**2*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+A(2,L)+2*A(3,L))),0._DP,0._DP,Q(5),Q(6),0._DP, &
        4*e2**3*(e2-e5*(A(1,L)+A(2,L)+A(3,L))), &
        -(e2*e4*(2*e2**2+2*e5**2*A(3,L)-e2*e5*(A(1,L)+A(2,L)+4*A(3,L)))),0._DP,Q(5),Q(6),0._DP,0._DP, &
        4*e2*e3*(e2**2+e4*e5*A(3,L)),0._DP,0._DP,Q(3),0._DP,Q(4),2*e2*(2*e2**2-2*e2*e5+e5**2)*Q(9),Q(6), &
        4*e2*(e2-e5*(A(1,L)+A(2,L)))*(e2**2+e4*e5*A(3,L)),0._DP,0._DP,Q(7),2*(2*e2**2-2*e2*e5+e5**2)*Q(9)**2, &
        4*e2**2*(e2-e5*(A(1,L)+A(3,L)))*(e2-e5*(A(2,L)+A(3,L))),4*e2**2*(e2-e5*(A(1,L)+A(2,L)))*Q(9), &
        -(e4*(2*e2-e5*(A(1,L)+A(2,L)))*(e2**2+e4*e5*A(3,L))),0._DP,-(e2*e4*(2*e2-e5*(A(1,L)+A(2,L)+2*A(3,L)))*Q(9)), &
        2*e2**2*(2*e2**2-2*e2*e5*(A(1,L)+A(2,L)+A(3,L))+e5**2*(2*A(1,L)*A(2,L)+(A(1,L)+A(2,L))*A(3,L))),Q(8), &
        4*e2**2*(e2-e5*(A(1,L)+A(2,L)+A(3,L)))*Q(9),4*e2*e3*Q(9)**2,Q(7), &
        2*e2*(e4*e5*(2*e2-e5*(A(1,L)+A(2,L)))*A(3,L)+2*e2*Q(10)*Q(11)),0._DP,0._DP,Q(8), &
        4*e2*(e2-e5*(A(1,L)+A(2,L)))*Q(9)**2,4*(e2**2+e4*e5*A(3,L))*Q(10)*Q(11),0._DP, &
        2*e2*(2*e2**2-2*e2*e5*(A(1,L)+A(2,L)+A(3,L))+e5**2*(2*A(1,L)*A(2,L)+(A(1,L)+A(2,L))*A(3,L)))*Q(9), &
        -(e4*(2*e2-e5*(A(1,L)+A(2,L)))*Q(9)**2),4*e2*Q(9)*Q(10)*Q(11),4*Q(9)**2*Q(10)*Q(11)/)
      CASE (18_I2B)
        Q(1:11)=(/e2 + A(1,L) - e1*A(1,L), e2 + A(2,L) - e1*A(2,L), e2 - e5*(A(1,L) &
          + A(2,L)), 2*e2 - e5*(A(1,L) + A(2,L)), 4*e2**2 - 4*e2*e5*(A(1,L) + A(2,L)) + &
          e5**2*(A(1,L) + A(2,L)), 2*e2**2 - 3*e2*e5*(A(1,L) + A(2,L)) + e5**2*(A(1,L) &
          + A(2,L)), e2**2 - 2*e2*e5*(A(1,L) + A(2,L)) + e5**2*(A(1,L) + A(2,L)), &
          2*e2**3 - 2*e5**3*A(1,L)*A(2,L) - 4*e2**2*e5*(A(1,L) + A(2,L)) + &
          e2*e5**2*(A(1,L) + A(1,L)**2 + A(2,L) + 4*A(1,L)*A(2,L) + A(2,L)**2), 4*e2**3 &
          - 2*e5**3*A(1,L)*A(2,L) - 6*e2**2*e5*(A(1,L) + A(2,L)) + e2*e5**2*(A(1,L) + &
          A(1,L)**2 + A(2,L) + 6*A(1,L)*A(2,L) + A(2,L)**2), 2*e2**4 - &
          4*e2*e5**3*A(1,L)*A(2,L) + e5**4*A(1,L)*A(2,L) - 4*e2**3*e5*(A(1,L) + A(2,L)) &
          + e2**2*e5**2*(A(1,L) + A(1,L)**2 + A(2,L) + 6*A(1,L)*A(2,L) + A(2,L)**2), &
          4*e2**4 - 8*e2**3*e5*(A(1,L) + A(2,L)) + e5**4*A(1,L)*A(2,L)*(A(1,L) + &
          A(2,L)) - 4*e2*e5**3*A(1,L)*A(2,L)*(1 + A(1,L) + A(2,L)) + &
          e2**2*e5**2*(A(1,L) + 3*A(1,L)**2 + A(2,L) + 14*A(1,L)*A(2,L) + &
          3*A(2,L)**2)/)
        coe(L,:)=(/4*e2**4, 4*e2**4*A(5,L), 2*e2**3*Q(4), 2*e2**3*Q(4), &
          4*e2**4*A(5,L), 4*e2**4*A(5,L), e2**2*Q(5), 4*e2**2*Q(7), 4*e2**3*Q(3), &
          4*e2**3*Q(3), 2*e2**2*Q(6), 4*e2**4*A(5,L), 2*e2**2*Q(6), 4*e2**4*A(5,L), &
          4*e2**2*Q(3)**2, 4*e2**3*A(5,L)*Q(3), 4*e2**3*A(5,L)*Q(3), 2*Q(10), &
          e2**2*A(5,L)*Q(5), 4*e2**2*A(5,L)*Q(7), 4*e2**2*Q(1)*Q(2), 4*e2**2*Q(1)*Q(2), &
          2*e2**3*A(5,L)*Q(4), 2*e2**3*A(5,L)*Q(4), e2**2*(-2*e2 + e5*(A(1,L) + &
          A(2,L)))**2, 2*e2**2*A(5,L)*Q(6), 2*e2**2*A(5,L)*Q(6), e2*Q(9), &
          2*e2**3*A(5,L)*Q(4), 2*e2**3*A(5,L)*Q(4), e2*Q(9), 2*e2*Q(8), &
          4*e2**3*A(5,L)*Q(3), 2*e2*Q(8), 2*e2**2*Q(3)*Q(4), 4*e2**3*A(5,L)*Q(3), &
          2*e2**2*Q(3)*Q(4), e2**2*A(5,L)*Q(5), 2*e2**3*A(5,L)*Q(4), &
          2*e2**3*A(5,L)*Q(4), 4*e2**2*A(5,L)*Q(3)**2, 4*e2*Q(1)*Q(2)*Q(3), &
          4*e2*Q(1)*Q(2)*Q(3), 4*e2**2*A(5,L)*Q(1)*Q(2), 4*e2**2*A(5,L)*Q(1)*Q(2), &
          2*A(5,L)*Q(10), 2*e2**2*A(5,L)*Q(3)*Q(4), 2*e2**2*A(5,L)*Q(3)*Q(4), Q(11), &
          e2*A(5,L)*Q(9), e2*A(5,L)*Q(9), e2**2*(-2*e2 + e5*(A(1,L) + &
          A(2,L)))**2*A(5,L), 2*e2*A(5,L)*Q(8), 4*e2**2*A(5,L)*Q(1)*Q(2), &
          4*e2**2*A(5,L)*Q(1)*Q(2), 2*e2*A(5,L)*Q(8), 2*e2*Q(1)*Q(2)*Q(4), &
          2*e2*Q(1)*Q(2)*Q(4), e2**2*(-2*e2 + e5*(A(1,L) + A(2,L)))**2*A(5,L), &
          4*e2*A(5,L)*Q(1)*Q(2)*Q(3), 4*e2*A(5,L)*Q(1)*Q(2)*Q(3), 4*Q(1)**2*Q(2)**2, &
          A(5,L)*Q(11), 2*e2*A(5,L)*Q(1)*Q(2)*Q(4), 2*e2*A(5,L)*Q(1)*Q(2)*Q(4), &
          4*A(5,L)*Q(1)**2*Q(2)**2/)
      CASE (19_I2B)
        IF(G(1,L)==G(3,L))THEN
          X(1)=A(4,L)
          X(2)=A(1,L)
          X(3)=A(2,L)
        ELSE IF(G(1,L)==G(4,L))THEN
          X(1)=A(3,L)
          X(2)=A(1,L)
          X(3)=A(2,L)
        ELSE IF(G(1,L)==G(5,L))THEN
          X(1)=A(6,L)
          X(2)=A(2,L)
          X(3)=A(1,L)
        ELSE
          X(1)=A(5,L)
          X(2)=A(2,L)
          X(3)=A(1,L)
        END IF
        Q(1:8)=(/2*e2**2*e4**2,-4*e2**2*e4*(e2-e5*(x(2)+x(3))),-4*e2**2*e4*(e2-e5*(x(1)+x(2))), &
        -4*e2*e4*(e2+x(2)-e1*x(2))*(e2+x(3)-e1*x(3)),-4*e2*e4*(e2+x(1)-e1*x(1))*(e2+x(2)-e1*x(2)), &
        e2+x(1)-e1*x(1),e2+x(2)-e1*x(2),e2+x(3)-e1*x(3)/)
        coe(L,:)=(/0._DP,8*e2**3*e3,0._DP,0._DP,Q(1),Q(1),0._DP,0._DP,0._DP,0._DP,0._DP, &
        4*e2**2*(2*e2**2-2*e2*e5+e5**2),0._DP,Q(1),0._DP,Q(2),Q(3),0._DP, &
        8*e2**4-4*e2*e5**3*x(2)+e5**4*x(2)-4*e2**3*e5*(1+x(1)+2*x(2)+x(3))+ &
         e2**2*e5**2*(1+x(1)+8*x(2)+x(3)),8*e2**2*(e2**2+e5**2*x(2)-e2*e5*(x(1)+2*x(2)+x(3))),0._DP,0._DP, &
        2*e2**2*(4*e2**2+e5**2*(x(1)+2*x(2))-e2*e5*(3+2*x(1)+2*x(2))), &
        2*e2**2*(4*e2**2+e5**2*(2*x(2)+x(3))-e2*e5*(3+2*x(2)+2*x(3))),0._DP, &
        2*e2*(4*e2**3-e5**3*x(2)+e2*e5**2*(4*x(2)+x(3))-e2**2*e5*(1+2*x(1)+6*x(2)+4*x(3))), &
        2*e2*(4*e2**3-e5**3*x(2)+e2*e5**2*(x(1)+4*x(2))-e2**2*e5*(1+4*x(1)+6*x(2)+2*x(3))),0._DP, &
        2*e2*(4*e2**3-e5**3*x(2)+e2*e5**2*(1+x(1)+2*x(2))-e2**2*e5*(3+2*x(1)+2*x(2))), &
        2*e2*(4*e2**3-e5**3*x(2)+e2*e5**2*(1+2*x(2)+x(3))-e2**2*e5*(3+2*x(2)+2*x(3))), &
        0._DP,0._DP,Q(3),0._DP,0._DP,Q(2),0._DP, &
        e2*(8*e2**3-2*e5**3*x(2)-4*e2**2*e5*(1+x(1)+2*x(2)+x(3))+e2*e5**2*(x(1)+8*x(2)+x(3))), &
        -(e2*e4*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(1)+2*x(2)))), &
        -(e2*e4*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(2)+2*x(3)))), &
        8*e2**2*(e2-e5*(x(1)+x(2)))*(e2-e5*(x(2)+x(3))),0._DP,0._DP,Q(4),Q(5), &
        4*e2*(2*e2**3-e5**3*x(2)*(x(1)+x(3))-2*e2**2*e5*(x(1)+2*x(2)+x(3))+ &
           e2*e5**2*(x(1)*(3*x(2)+x(3))+x(2)*(1+x(2)+3*x(3)))), &
        2*e2*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(1)+2*x(2)))*(e2-e5*(x(2)+x(3))), &
        2*e2*(e2-e5*(x(1)+x(2)))*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(2)+2*x(3))),0._DP, &
        8*e2**4+e5**4*x(1)*x(2)-e2*e5**3*x(2)*(1+5*x(1)+x(2))-2*e2**3*e5*(1+4*x(1)+6*x(2)+2*x(3))+ &
         2*e2**2*e5**2*(x(2)*(3+x(2)+x(3))+x(1)*(1+5*x(2)+x(3))), &
        8*e2**4+e5**4*x(2)*x(3)-2*e2**3*e5*(1+2*x(1)+6*x(2)+4*x(3))-e2*e5**3*x(2)*(1+x(2)+5*x(3))+ &
         2*e2**2*e5**2*(x(2)*(3+x(1)+x(2))+(1+x(1)+5*x(2))*x(3)), &
        8*e2**4+e5**4*x(2)**2-4*e2**3*e5*(1+x(1)+2*x(2)+x(3))-e2*e5**3*x(2)*(2+x(1)+2*x(2)+x(3))+ &
         e2**2*e5**2*(1+x(1)+2*x(2)*(3+x(1)+x(2))+(1+2*x(1)+2*x(2))*x(3)), &
        4*e2*(2*e2**3-e5**3*x(1)*x(2)-2*e2**2*e5*(x(1)+2*x(2)+x(3))+ &
           e2*e5**2*(x(2)*(1+x(2)+x(3))+x(1)*(3*x(2)+x(3)))),Q(5),Q(4), &
        4*e2*(2*e2**3-e5**3*x(2)*x(3)-2*e2**2*e5*(x(1)+2*x(2)+x(3))+ &
           e2*e5**2*(x(1)*(x(2)+x(3))+x(2)*(1+x(2)+3*x(3)))),0._DP,0._DP, &
        e2*(8*e2**3-e5**3*x(2)*(x(1)+2*x(2)+x(3))-4*e2**2*e5*(1+x(1)+2*x(2)+x(3))+ &
           e2*e5**2*(x(1)+2*x(2)*(3+x(1)+x(2))+(1+2*x(1)+2*x(2))*x(3))),8*e2*Q(6)*Q(7)*(e2-e5*(x(2)+x(3))), &
        8*e2*Q(7)*Q(8)*(e2-e5*(x(1)+x(2))),0._DP,2* &
         (4*e2**4+e5**4*x(1)*x(2)*x(3)-4*e2**3*e5*(x(1)+2*x(2)+x(3))- &
           e2*e5**3*x(2)*(x(1)+2*x(1)*x(2)+(1+4*x(1)+2*x(2))*x(3))+ &
           e2**2*e5**2*(3*x(2)**2+3*x(1)*x(3)+x(2)*(1+7*x(1)+7*x(3)))), &
        2*Q(6)*Q(7)*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(2)+2*x(3))), &
        2*Q(7)*Q(8)*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(1)+2*x(2))),8*Q(6)*Q(7)**2*Q(8)/)
      CASE (20_I2B)
        Q(1:10)=(/2*e2**4,2*e2**2*(e2**2+e4*e5*A(1,L)),2*e2**3*(e2+A(1,L)-e1*A(1,L)), &
        2*e2*(e2**3-e5*(3*e2**2-3*e2*e5+e5**2)*A(1,L)),2*(e2**2+e4*e5*A(1,L))**2, &
        2*e2**2*(e2+A(1,L)-e1*A(1,L))**2,2*e2*(e2+A(1,L)-e1*A(1,L))*(e2**2+e4*e5*A(1,L)), &
        2*(e2**4-4*e2**3*e5*A(1,L)+e5**4*A(1,L)**2+3*e2**2*e5**2*A(1,L)*(1+A(1,L))-e2*e5**3*A(1,L)*(1+3*A(1,L))), &
        2*(e2+A(1,L)-e1*A(1,L))**2*(e2**2+e4*e5*A(1,L)),2*e2*(e2+A(1,L)-e1*A(1,L))**3/)
        coe(L,:)=(/0._DP,Q(1),0._DP,0._DP,Q(1),Q(1),0._DP,0._DP,0._DP,0._DP,&
         0._DP,Q(1),0._DP,Q(1),0._DP,Q(2),Q(2),0._DP,Q(2), &
        2*(e2**4-4*e2**3*e5*A(1,L)+6*e2**2*e5**2*A(1,L)-4*e2*e5**3*A(1,L)+ &
        e5**4*A(1,L)),0._DP,0._DP,Q(3),Q(3),0._DP,Q(4),Q(4),0._DP, &
        Q(3),Q(3),0._DP,0._DP,Q(2),0._DP,0._DP,Q(2),0._DP,Q(2),Q(3),Q(3),Q(5),0._DP,0._DP,Q(6),&
        Q(6),Q(5),Q(7),Q(7),0._DP,Q(7),Q(7),Q(6),Q(8),Q(6),Q(6),Q(8), &
        0._DP,0._DP,Q(6),Q(9),Q(9),0._DP,Q(9),Q(10),Q(10),2*(e2+A(1,L)-e1*A(1,L))**4/)
      CASE (21_I2B)
        Q(1:12)=(/2*e2**2*e3**2*A(1,L),2*e2**2*A(1,L)*(e2**2+e4*e5*A(1,L)),&
        2*e2*e3**2*A(1,L)*(e2-e5*(A(3,L)+A(4,L))), &
        e2**2*e3*A(1,L)*(2*e2-e5*(2*A(1,L)+A(3,L)+A(4,L))),e2*e3**2*A(1,L)*(2*e2-e5*(A(3,L)+A(4,L))), &
        2*e2**2*e3*A(1,L)*(e2+A(1,L)-e1*A(1,L)),2*e2**2*A(1,L)*(e2+A(1,L)-e1*A(1,L))**2, &
        2*e3**2*A(1,L)*(e2+A(3,L)-e1*A(3,L))*(e2+A(4,L)-e1*A(4,L)), &
        e2*e3*A(1,L)*(e2+A(1,L)-e1*A(1,L))*(2*e2-e5*(A(3,L)+A(4,L))), &
        e2+A(1,L)-e1*A(1,L),e2+A(3,L)-e1*A(3,L),e2+A(4,L)-e1*A(4,L)/)
        coe(L,:)=(/2*e2**2*e3**2,Q(1),e2*e3**2*(2*e2-e5*(A(3,L)+A(4,L))),2*e2**2*e3*Q(10),Q(1),Q(1), &
        e2**2*e3*(2*e2-e5*(2*A(1,L)+A(3,L)+A(4,L))), &
        2*e2**2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)+A(4,L))),2*e2*e3**2*(e2-e5*(A(3,L)+A(4,L))), &
        2*e2**2*(e2**2+e4*e5*A(1,L)),e2**2*(2*e2**2+2*e5**2*A(1,L)-e2*e5*(4*A(1,L)+A(3,L)+A(4,L))),Q(1), &
        2*e2**2*e3*(e2-e5*(A(1,L)+A(3,L)+A(4,L))),Q(1),2*e2*(e2**2+e4*e5*A(1,L))*(e2-e5*(A(3,L)+A(4,L))),Q(2), &
        Q(3),2*e2**2*(e2-e5*(A(1,L)+A(3,L)))*(e2-e5*(A(1,L)+A(4,L))),Q(4), &
        2*e2**2*A(1,L)*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)+A(4,L))), &
        2*e3**2*Q(11)*Q(12),2*e2**2*Q(10)**2,Q(5), &
        Q(6),e2*e3*(2*e2-e5*(A(3,L)+A(4,L)))*Q(10), &
        e2**2*A(1,L)*(2*e2**2+2*e5**2*A(1,L)-e2*e5*(4*A(1,L)+A(3,L)+A(4,L))), &
        2*e2**2*e3*A(1,L)*(e2-e5*(A(1,L)+A(3,L)+A(4,L))), &
        e2*e3*(2*e2**2-2*e2*e5*(A(1,L)+A(3,L)+A(4,L))+ &
        e5**2*(2*A(3,L)*A(4,L)+A(1,L)*(A(3,L)+A(4,L)))),Q(5),Q(6), &
        e2**2*(2*e2-e5*(2*A(1,L)+A(3,L)+A(4,L)))*Q(10), &
        e2*(2*e2**3-e5**3*A(1,L)*(A(3,L)+A(4,L))-2*e2**2*e5*(2*A(1,L)+A(3,L)+A(4,L))+ &
           2*e2*e5**2*(A(3,L)*A(4,L)+A(1,L)*(1+A(3,L)+A(4,L)))),Q(3), &
        2*e2**2*(e2-e5*(A(1,L)+A(3,L)+A(4,L)))*Q(10),2*e2*e3*(e2-e5*(A(3,L)+A(4,L)))*Q(10),Q(2), &
        e2*(e2**2+e4*e5*A(1,L))*(2*e2-e5*(A(3,L)+A(4,L))),Q(4),Q(5),Q(6), &
        2*e2*A(1,L)*(e2**2+e4*e5*A(1,L))*(e2-e5*(A(3,L)+A(4,L))),2*(e2**2+e4*e5*A(1,L))*Q(11)*Q(12), &
        2*e2*(e2-e5*(A(3,L)+A(4,L)))*Q(10)**2,Q(7),Q(8), &
        2*e2**2*A(1,L)*(e2-e5*(A(1,L)+A(3,L)))*(e2-e5*(A(1,L)+A(4,L))), &
        e2*A(1,L)*(e2**2+e4*e5*A(1,L))*(2*e2-e5*(A(3,L)+A(4,L))),2*e2*e3*A(1,L)*(e2-e5*(A(3,L)+A(4,L)))*Q(10), &
        e2*(2*e2**2-2*e2*e5*(A(1,L)+A(3,L)+A(4,L))+e5**2*(2*A(3,L)*A(4,L)+A(1,L)*(A(3,L)+A(4,L))))*Q(10), &
        e2*e3*A(1,L)*(2*e2**2-2*e2*e5*(A(1,L)+A(3,L)+A(4,L))+e5**2*(2*A(3,L)*A(4,L)+A(1,L)*(A(3,L)+A(4,L)))), &
        e2**2*A(1,L)*(2*e2-e5*(2*A(1,L)+A(3,L)+A(4,L)))*Q(10),Q(9), &
        e2*A(1,L)*(2*e2**3-e5**3*A(1,L)*(A(3,L)+A(4,L))-2*e2**2*e5*(2*A(1,L)+A(3,L)+A(4,L))+ &
           2*e2*e5**2*(A(3,L)*A(4,L)+A(1,L)*(1+A(3,L)+A(4,L)))),Q(8),Q(7), &
        2*e2**2*A(1,L)*(e2-e5*(A(1,L)+A(3,L)+A(4,L)))*Q(10),2*e3*Q(10)*Q(11)*Q(12), &
        e2*(2*e2-e5*(A(3,L)+A(4,L)))*Q(10)**2,Q(9),2*A(1,L)*(e2**2+e4*e5*A(1,L))*Q(11)*Q(12), &
        2*e2*A(1,L)*(e2-e5*(A(3,L)+A(4,L)))*Q(10)**2,2*Q(10)**2*Q(11)*Q(12), &
        e2*A(1,L)*(2*e2**2-2*e2*e5*(A(1,L)+A(3,L)+A(4,L))+ &
        e5**2*(2*A(3,L)*A(4,L)+A(1,L)*(A(3,L)+A(4,L))))*Q(10), &
        2*e3*A(1,L)*Q(10)*Q(11)*Q(12),e2*A(1,L)*(2*e2-e5*(A(3,L)+ &
        A(4,L)))*Q(10)**2,2*A(1,L)*Q(10)**2*Q(11)*Q(12)/)
      CASE (22_I2B)
        Q(1:12)=(/2*e2**2*e3**2*A(3,L),2*e2*e3**2*(e2-e5*(A(1,L)+ &
        A(2,L)))*A(3,L),2*e2**2*A(3,L)*(e2**2+e4*e5*A(3,L)), &
        e2**2*e3*A(3,L)*(2*e2-e5*(A(1,L)+A(2,L)+2*A(3,L))), &
        2*e2**2*e3*A(3,L)*(e2+A(3,L)-e1*A(3,L)), &
        e2*e3**2*(2*e2-e5*(A(1,L)+A(2,L)))*A(3,L), &
        2*e3**2*(e2+A(1,L)-e1*A(1,L))*(e2+A(2,L)-e1*A(2,L))*A(3,L), &
        2*e2**2*A(3,L)*(e2+A(3,L)-e1*A(3,L))**2, &
        e2*e3*(2*e2-e5*(A(1,L)+A(2,L)))*A(3,L)*(e2+A(3,L)-e1*A(3,L)), &
        e2+A(3,L)-e1*A(3,L),e2+A(2,L)-e1*A(2,L),e2+A(1,L)-e1*A(1,L)/)
        coe(L,:)=(/2*e2**2*e3**2,Q(1),2*e2**2*e3*Q(10),e2*e3**2*(2*e2-e5*(A(1,L)+A(2,L))),Q(1),Q(1), &
        e2**2*e3*(2*e2-e5*(A(1,L)+A(2,L)+2*A(3,L))), &
        2*e2**2*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+A(2,L)+2*A(3,L))),2*e2**2*(e2**2+e4*e5*A(3,L)), &
        2*e2*e3**2*(e2-e5*(A(1,L)+A(2,L))),2*e2**2*e3*(e2-e5*(A(1,L)+A(2,L)+A(3,L))),Q(1), &
        e2**2*(2*e2**2+2*e5**2*A(3,L)-e2*e5*(A(1,L)+A(2,L)+4*A(3,L))),Q(1), &
        2*e2*(e2-e5*(A(1,L)+A(2,L)))*(e2**2+e4*e5*A(3,L)),Q(2),Q(3), &
        2*e2**2*(e2-e5*(A(1,L)+A(3,L)))*(e2-e5*(A(2,L)+A(3,L))),Q(4), &
        2*e2**2*A(3,L)*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+A(2,L)+2*A(3,L))),&
        2*e2**2*Q(10)**2,2*e3**2*Q(11)*Q(12),Q(5), &
        Q(6),e2*e3*(2*e2-e5*(A(1,L)+A(2,L)))*Q(10),2*e2**2*e3*A(3,L)*(e2-e5*(A(1,L)+A(2,L)+A(3,L))), &
        e2**2*A(3,L)*(2*e2**2+2*e5**2*A(3,L)-e2*e5*(A(1,L)+A(2,L)+4*A(3,L))), &
        e2**2*(2*e2-e5*(A(1,L)+A(2,L)+2*A(3,L)))*Q(10),Q(5),Q(6), &
        e2*e3*(2*e2**2-2*e2*e5*(A(1,L)+A(2,L)+A(3,L))+e5**2*(2*A(1,L)*A(2,L)+(A(1,L)+A(2,L))*A(3,L))), &
        2*e2**2*(e2-e5*(A(1,L)+A(2,L)+A(3,L)))*Q(10),Q(3), &
        e2*(e4*e5*(2*e2-e5*(A(1,L)+A(2,L)))*A(3,L)+2*e2*Q(11)*Q(12)), &
        e2*(2*e2-e5*(A(1,L)+A(2,L)))*(e2**2+e4*e5*A(3,L)), &
        Q(2),2*e2*e3*(e2-e5*(A(1,L)+A(2,L)))*Q(10),Q(4),Q(5), &
        Q(6),2*e2*(e2-e5*(A(1,L)+A(2,L)))*A(3,L)*(e2**2+e4*e5*A(3,L)),&
        2*e2*(e2-e5*(A(1,L)+A(2,L)))*Q(10)**2, &
        2*(e2**2+e4*e5*A(3,L))*Q(11)*Q(12),Q(7),Q(8), &
        2*e2**2*A(3,L)*(e2-e5*(A(1,L)+A(3,L)))*(e2-e5*(A(2,L)+A(3,L))), &
        2*e2*e3*(e2-e5*(A(1,L)+A(2,L)))*A(3,L)*Q(10),&
        e2*(2*e2-e5*(A(1,L)+A(2,L)))*A(3,L)*(e2**2+e4*e5*A(3,L)), &
        e2*(2*e2**2-2*e2*e5*(A(1,L)+A(2,L)+A(3,L))+e5**2* &
        (2*A(1,L)*A(2,L)+(A(1,L)+A(2,L))*A(3,L)))*Q(10), &
        e2**2*A(3,L)*(2*e2-e5*(A(1,L)+A(2,L)+2*A(3,L)))*Q(10), &
        e2*e3*A(3,L)*(2*e2**2-2*e2*e5*(A(1,L)+A(2,L)+A(3,L))+ &
        e5**2*(2*A(1,L)*A(2,L)+(A(1,L)+A(2,L))*A(3,L))), &
        Q(9),2*e2**2*A(3,L)*(e2-e5*(A(1,L)+A(2,L)+A(3,L)))*Q(10),Q(8),Q(7), &
        e2*A(3,L)*(e4*e5*(2*e2-e5*(A(1,L)+A(2,L)))*A(3,L)+2*e2*Q(11)*Q(12)), &
        e2*(2*e2-e5*(A(1,L)+A(2,L)))*Q(10)**2, &
        2*e3*Q(10)*Q(11)*Q(12),Q(9),2*e2*(e2-e5*(A(1,L)+A(2,L)))*A(3,L)*Q(10)**2, &
        2*A(3,L)*(e2**2+e4*e5*A(3,L))*Q(11)*Q(12),2*Q(10)**2*Q(11)*Q(12), &
        e2*A(3,L)*(2*e2**2-2*e2*e5*(A(1,L)+A(2,L)+A(3,L))+ &
        e5**2*(2*A(1,L)*A(2,L)+(A(1,L)+A(2,L))*A(3,L)))*Q(10), &
        e2*(2*e2-e5*(A(1,L)+A(2,L)))*A(3,L)*Q(10)**2, &
        2*e3*A(3,L)*Q(10)*Q(11)*Q(12),2*A(3,L)*Q(10)**2*Q(11)*Q(12)/)
      CASE (23_I2B)
        QQ=Merge(A(4,L), A(3,L), G(1,L)==G(3,L))
        Q(1:6)=(/-2*e2**2*e3*e4,-2*e2*e4*(e2**2+e4*e5*A(1,L)), &
        4*e2*(e2**2+e4*e5*A(1,L))*(e2-e5*(QQ+A(1,L))), &
        -2*e2*e4*(e2+A(1,L)-e1*A(1,L))**2,e2+A(1,L)-e1*A(1,L),e2+QQ-e1*QQ/)
        coe(L,:)=(/0._DP,-2*e2*e4*(e2-e1*e2+e2**2+e5**2),0._DP,0._DP,&
        -(e2*e4*(2*e2**2-2*e2*e5+e5**2)),Q(1),0._DP,0._DP,0._DP,&
        0._DP,0._DP,Q(1),0._DP,Q(1),0._DP, &
        Q(2),4*e2**2*e3*(e2-e5*(QQ+A(1,L))),0._DP,e2* &
         (4*e2**3-e5**3*A(1,L)+2*e2*e5**2*(QQ+3*A(1,L))-e2**2*e5*(3+2*QQ+6*A(1,L))), &
        4*e2*(e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(QQ+3*A(1,L))),0._DP,0._DP, &
        4*e2**4+e5**4*QQ-e2*e5**3*(1+3*QQ+A(1,L))-e2**3*e5*(5+2*QQ+2*A(1,L))+ &
         e2**2*e5**2*(4+3*QQ+2*A(1,L)),2*e2*(2*e2**2-2*e2*e5+e5**2)*Q(5),0._DP, &
        e2*(4*e2**3-3*e5**3*A(1,L)+e2*e5**2*(QQ+10*A(1,L))-e2**2*e5*(1+2*QQ+10*A(1,L))), &
        -2*e2*e4*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L))),0._DP, &
        e2*(4*e2**3-e5**3*(QQ+A(1,L))-e2**2*e5*(5+2*QQ+2*A(1,L))+ &
        e2*e5**2*(2+3*QQ+2*A(1,L))),4*e2**2*e3*Q(5), &
        0._DP,0._DP,2*e2*(2*e2**2-2*e2*e5+e5**2)*(e2-e5*(QQ+A(1,L))),0._DP,0._DP,Q(2),0._DP, &
        e2*(4*e2**3-e5**3*(QQ+2*A(1,L))+e2*e5**2*(1+2*QQ+6*A(1,L))-e2**2*e5*(3+2*QQ+6*A(1,L))), &
        e2*e3*(4*e2**2+e5**2*QQ-e2*e5*(1+2*QQ+2*A(1,L))),&
        e2*e4**2*Q(5),Q(3),0._DP,0._DP,Q(4),4*e2*e3*Q(5)*Q(6),Q(3), &
        (e2**2+e4*e5*A(1,L))*(4*e2**2+e5**2*QQ-e2*e5*(1+2*QQ+2*A(1,L))),&
        -2*e2*e4*(e2-e5*(QQ+A(1,L)))*Q(5),0._DP, &
        -(e4*(2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L)))), &
        e2*(4*e2**2+e5**2*(QQ+3*A(1,L))-e2*e5*(1+2*QQ+6*A(1,L)))*Q(5), &
        e2*(4*e2**2+e5**2*(2*QQ+A(1,L))-e2*e5*(3+2*QQ+2*A(1,L)))*Q(5), &
        2*(2*e2**4+e5**4*QQ*A(1,L)+2*e2**2*e5**2*A(1,L)*(2+2*QQ+A(1,L))-e2*e5**3*A(1,L)*(1+3*QQ+A(1,L))- &
           2*e2**3*e5*(QQ+3*A(1,L))),2*(2*e2**2-2*e2*e5+e5**2)*Q(5)*Q(6),Q(4), &
        4*e2*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L)))*Q(5),0._DP,0._DP, &
        (4*e2**3-e5**3*QQ+e2*e5**2*(1+2*QQ+A(1,L))-e2**2*e5*(3+2*QQ+2*A(1,L)))*Q(5), &
        4*(e2**2+e4*e5*A(1,L))*Q(5)*Q(6),4*e2*(e2-e5*(QQ+A(1,L)))*Q(5)**2,0._DP, &
        2*(2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L)))*Q(5), &
        -2*e4*Q(5)**2*Q(6),(4*e2**2+e5**2*QQ-e2*e5*(1+2*QQ+2*A(1,L)))*Q(5)**2,4*Q(5)**3*Q(6)/)
      CASE (24_I2B)
        QQ=Merge(A(2,L), A(1,L), G(1,L)==G(3,L))
        Q(1:6)=(/-2*e2**2*e3*e4,-2*e2*e4*(e2**2+e4*e5*A(3,L)),4*e2*(e2**2+e4*e5*A(3,L))*(e2-e5*(QQ+A(3,L))), &
        -2*e2*e4*(e2+A(3,L)-e1*A(3,L))**2,e2+A(3,L)-e1*A(3,L),e2+QQ-e1*QQ/)
        coe(L,:)=(/0._DP,-2*e2*e4*(e2-e1*e2+e2**2+e5**2),0._DP,0._DP,Q(1),-(e2*e4*(2*e2**2-2*e2*e5+e5**2)),&
        0._DP,0._DP,0._DP,0._DP,0._DP,Q(1),0._DP,Q(1),0._DP, &
        4*e2**2*e3*(e2-e5*(QQ+A(3,L))),Q(2),0._DP,e2* &
         (4*e2**3-e5**3*A(3,L)+2*e2*e5**2*(QQ+3*A(3,L))-e2**2*e5*(3+2*QQ+6*A(3,L))), &
        4*e2*(e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(QQ+3*A(3,L))), &
        0._DP,0._DP,2*e2*(2*e2**2-2*e2*e5+e5**2)*Q(5), &
        4*e2**4+e5**4*QQ-e2*e5**3*(1+3*QQ+A(3,L))-e2**3*e5*(5+2*QQ+2*A(3,L))+ &
         e2**2*e5**2*(4+3*QQ+2*A(3,L)),0._DP,-2*e2*e4*(e2**2+e5**2*A(3,L)-e2*e5*(QQ+2*A(3,L))), &
        e2*(4*e2**3-3*e5**3*A(3,L)+e2*e5**2*(QQ+10*A(3,L))-e2**2*e5*(1+2*QQ+10*A(3,L))),0._DP,4*e2**2*e3*Q(5), &
        e2*(4*e2**3-e5**3*(QQ+A(3,L))-e2**2*e5*(5+2*QQ+2*A(3,L))+e2*e5**2*(2+3*QQ+2*A(3,L))), &
        0._DP,0._DP,Q(2),0._DP,0._DP,2*e2*(2*e2**2-2*e2*e5+e5**2)*(e2-e5*(QQ+A(3,L))),0._DP, &
        e2*(4*e2**3-e5**3*(QQ+2*A(3,L))+e2*e5**2*(1+2*QQ+6*A(3,L))-e2**2*e5*(3+2*QQ+6*A(3,L))),e2*e4**2*Q(5), &
        e2*e3*(4*e2**2+e5**2*QQ-e2*e5*(1+2*QQ+2*A(3,L))),Q(3),0._DP,0._DP,4*e2*e3*Q(5)*Q(6),Q(4),Q(3), &
        -2*e2*e4*(e2-e5*(QQ+A(3,L)))*Q(5),(e2**2+e4*e5*A(3,L))*(4*e2**2+e5**2*QQ-e2*e5*(1+2*QQ+2*A(3,L))),0._DP, &
        e2*(4*e2**2+e5**2*(QQ+3*A(3,L))-e2*e5*(1+2*QQ+6*A(3,L)))*Q(5), &
        -(e4*(2*e2**3-e5**3*QQ*A(3,L)+e2*e5**2*A(3,L)*(1+3*QQ+A(3,L))-2*e2**2*e5*(QQ+2*A(3,L)))), &
        e2*(4*e2**2+e5**2*(2*QQ+A(3,L))-e2*e5*(3+2*QQ+2*A(3,L)))*Q(5), &
        4*e2*(e2**2+e5**2*A(3,L)-e2*e5*(QQ+2*A(3,L)))*Q(5),Q(4),2*(2*e2**2-2*e2*e5+e5**2)*Q(5)*Q(6), &
        2*(2*e2**4+e5**4*QQ*A(3,L)+2*e2**2*e5**2*A(3,L)*(2+2*QQ+A(3,L))-e2*e5**3*A(3,L)*(1+3*QQ+A(3,L))- &
           2*e2**3*e5*(QQ+3*A(3,L))),0._DP,0._DP,(4*e2**3-e5**3*QQ+e2*e5**2*(1+2*QQ+A(3,L))- &
           e2**2*e5*(3+2*QQ+2*A(3,L)))*Q(5),4*e2*(e2-e5*(QQ+A(3,L)))*Q(5)**2,&
           4*(e2**2+e4*e5*A(3,L))*Q(5)*Q(6),0._DP, &
        2*(2*e2**3-e5**3*QQ*A(3,L)+e2*e5**2*A(3,L)*(1+3*QQ+A(3,L))-2*e2**2*e5*(QQ+2*A(3,L)))*Q(5), &
        (4*e2**2+e5**2*QQ-e2*e5*(1+2*QQ+2*A(3,L)))*Q(5)**2,-2*e4*Q(5)**2*Q(6),4*Q(5)**3*Q(6)/)
      CASE (25_I2B)
        IF(ALL(G(5,L)==G(1:3:2,L)))THEN
          x(1:3)=(/A(4,L),A(1,L),A(2,L)/)
        ELSE IF(ALL(G(5,L)==G(1:4:3,L)))THEN
          x(1:3)=(/A(3,L),A(1,L),A(2,L)/)
        ELSE IF(ALL(G(5,L)==G(2:3,L)))THEN
          x(1:3)=(/A(4,L),A(2,L),A(1,L)/)
        ELSE
          x(1:3)=(/A(3,L),A(2,L),A(1,L)/)
        END IF
         Q(1:21)=(/e3, e2 + x(1) - e1*x(1), -e2 + e5*x(2), e2 + x(2) - e1*x(2), e2 - &
         e5*(x(1) + x(2)), 2*e2**2 + e5**2*x(1) - e2*e5*(1 + x(1) + x(2)), -2*e2**2 - &
         e5**2*x(1) + e2*e5*(1 + x(1) + x(2)), e2 + x(3) - e1*x(3), e2 - e5*(x(2) + &
         x(3)), 2*e2**2 + e5**2*x(3) - e2*e5*(1 + x(2) + x(3)), -2*e2**2 - e5**2*x(3) &
         + e2*e5*(1 + x(2) + x(3)), e2**2 + e5**2*x(2) - e2*e5*(x(1) + 2*x(2) + x(3)), &
         2*e2**2 + e5**2*(x(1) + 2*x(2) + x(3)) - e2*e5*(1 + 2*x(1) + 3*x(2) + x(3)), &
         2*e2**2 + e5**2*(x(1) + 2*x(2) + x(3)) - e2*e5*(1 + x(1) + 3*x(2) + 2*x(3)), &
         4*e2**3 - e5**3*(x(1) + x(3)) - 2*e2**2*e5*(2 + x(1) + 2*x(2) + x(3)) + &
         e2*e5**2*(1 + 3*x(1) + 3*x(2) + 3*x(3)), 2*e2**3 - e5**3*x(1)*x(2) - &
         2*e2**2*e5*(x(1) + 2*x(2) + x(3)) + e2*e5**2*(x(2)*(1 + 3*x(1) + x(2)) + &
         (x(1) + x(2))*x(3)), 2*e2**3 - e5**3*x(2)*x(3) - 2*e2**2*e5*(x(1) + 2*x(2) + &
         x(3)) + e2*e5**2*(x(2)*(1 + x(1) + x(2)) + (x(1) + 3*x(2))*x(3)), 2*e2**3 - &
         e5**3*x(2)*(x(1) + x(3)) - 2*e2**2*e5*(x(1) + 2*x(2) + x(3)) + &
         e2*e5**2*(x(2)*(1 + 3*x(1) + x(2)) + (x(1) + 3*x(2))*x(3)), 4*e2**3 - &
         2*e2**2*e5*(1 + 2*x(1) + 3*x(2) + x(3)) - e5**3*(3*x(1)*x(2) + (x(1) + &
         x(2))*x(3)) + e2*e5**2*(2*x(3) + x(2)*(3 + x(2) + x(3)) + x(1)*(2 + 5*x(2) + &
         x(3))), 4*e2**3 - 2*e2**2*e5*(1 + x(1) + 3*x(2) + 2*x(3)) - &
         e5**3*(3*x(2)*x(3) + x(1)*(x(2) + x(3))) + e2*e5**2*(2*x(3) + x(1)*(2 + x(2) &
         + x(3)) + x(2)*(3 + x(2) + 5*x(3))), 4*e2**4 + e5**4*x(1)*x(2)*x(3) - &
         4*e2**3*e5*(x(1) + 2*x(2) + x(3)) - e2*e5**3*x(2)*(x(3) + 2*x(2)*x(3) + &
         x(1)*(1 + 2*x(2) + 4*x(3))) + e2**2*e5**2*(3*x(2)**2 + 3*x(1)*x(3) + x(2)*(1 &
         + 7*x(1) + 7*x(3)))/)
         coe(L,:)=(/4*e2**2*Q(1)**2, 4*e2**2*Q(1)**2*x(2), 2*e2*Q(1)*Q(6), &
         2*e2*Q(1)*Q(10), 4*e2**2*Q(1)**2*x(2), 4*e2**2*Q(1)**2*x(2), e2*Q(15), &
         4*e2**2*Q(12), 4*e2**2*Q(1)*Q(5), 4*e2**2*Q(1)*Q(9), 2*e2**2*Q(14), &
         4*e2**2*Q(1)**2*x(2), 2*e2**2*Q(13), 4*e2**2*Q(1)**2*x(2), 4*e2**2*Q(5)*Q(9), &
         4*e2**2*Q(1)*Q(9)*x(2), 4*e2**2*Q(1)*Q(5)*x(2), 2*e2*(2*e2**3 - &
         e5**3*x(2)*(x(1) + x(3)) - 2*e2**2*e5*(x(1) + 2*x(2) + x(3)) + &
         e2*e5**2*(x(1)*(3*x(2) + x(3)) + x(2)*(1 + x(2) + 3*x(3)))), e2*Q(15)*x(2), &
         4*e2**2*Q(12)*x(2), 4*e2*Q(1)*Q(2)*Q(4), 4*e2*Q(1)*Q(4)*Q(8), &
         2*e2*Q(1)*Q(6)*x(2), 2*e2*Q(1)*Q(10)*x(2), Q(6)*Q(10), 2*e2**2*Q(14)*x(2), &
         2*e2**2*Q(13)*x(2), e2*(4*e2**3 - 2*e2**2*e5*(1 + 2*x(1) + 3*x(2) + x(3)) - &
         e5**3*(x(1)*x(3) + x(2)*(3*x(1) + x(3))) + e2*e5**2*(x(2)**2 + 2*x(3) + &
         x(1)*(2 + x(3)) + x(2)*(3 + 5*x(1) + x(3)))), 2*e2*Q(1)*Q(6)*x(2), &
         2*e2*Q(1)*Q(10)*x(2), e2*(4*e2**3 - 2*e2**2*e5*(1 + x(1) + 3*x(2) + 2*x(3)) - &
         e5**3*(3*x(2)*x(3) + x(1)*(x(2) + x(3))) + e2*e5**2*(x(2)**2 + 2*x(3) + &
         x(1)*(2 + x(3)) + x(2)*(3 + x(1) + 5*x(3)))), 2*e2*(2*e2**3 - e5**3*x(1)*x(2) &
         - 2*e2**2*e5*(x(1) + 2*x(2) + x(3)) + e2*e5**2*(x(2)*(1 + x(2) + x(3)) + &
         x(1)*(3*x(2) + x(3)))), 4*e2**2*Q(1)*Q(5)*x(2), 2*e2*(2*e2**3 - &
         e5**3*x(2)*x(3) - 2*e2**2*e5*(x(1) + 2*x(2) + x(3)) + e2*e5**2*(x(1)*(x(2) + &
         x(3)) + x(2)*(1 + x(2) + 3*x(3)))), 2*e2*Q(5)*Q(10), 4*e2**2*Q(1)*Q(9)*x(2), &
         2*e2*Q(6)*Q(9), e2*Q(15)*x(2), 2*e2*Q(1)*Q(6)*x(2), 2*e2*Q(1)*Q(10)*x(2), &
         4*e2**2*Q(5)*Q(9)*x(2), 4*e2*Q(2)*Q(4)*Q(9), 4*e2*Q(4)*Q(5)*Q(8), &
         4*e2*Q(1)*Q(4)*Q(8)*x(2), 4*e2*Q(1)*Q(2)*Q(4)*x(2), 2*e2*x(2)*(2*e2**3 - &
         e5**3*x(2)*(x(1) + x(3)) - 2*e2**2*e5*(x(1) + 2*x(2) + x(3)) + &
         e2*e5**2*(x(1)*(3*x(2) + x(3)) + x(2)*(1 + x(2) + 3*x(3)))), &
         2*e2*Q(6)*Q(9)*x(2), 2*e2*Q(5)*Q(10)*x(2), 4*e2**4 + e5**4*x(1)*x(2)*x(3) - &
         4*e2**3*e5*(x(1) + 2*x(2) + x(3)) - e2*e5**3*x(2)*(x(1) + 2*x(1)*x(2) + (1 + &
         4*x(1) + 2*x(2))*x(3)) + e2**2*e5**2*(3*x(2)**2 + 3*x(1)*x(3) + x(2)*(1 + &
         7*x(1) + 7*x(3))), e2*x(2)*(4*e2**3 - 2*e2**2*e5*(1 + 2*x(1) + 3*x(2) + x(3)) &
         - e5**3*(x(1)*x(3) + x(2)*(3*x(1) + x(3))) + e2*e5**2*(x(2)**2 + 2*x(3) + &
         x(1)*(2 + x(3)) + x(2)*(3 + 5*x(1) + x(3)))), e2*x(2)*(4*e2**3 - &
         2*e2**2*e5*(1 + x(1) + 3*x(2) + 2*x(3)) - e5**3*(3*x(2)*x(3) + x(1)*(x(2) + &
         x(3))) + e2*e5**2*(x(2)**2 + 2*x(3) + x(1)*(2 + x(3)) + x(2)*(3 + x(1) + &
         5*x(3)))), Q(6)*Q(10)*x(2), 2*e2*x(2)*(2*e2**3 - e5**3*x(1)*x(2) - &
         2*e2**2*e5*(x(1) + 2*x(2) + x(3)) + e2*e5**2*(x(2)*(1 + x(2) + x(3)) + &
         x(1)*(3*x(2) + x(3)))), 4*e2*Q(1)*Q(2)*Q(4)*x(2), 4*e2*Q(1)*Q(4)*Q(8)*x(2), &
         2*e2*x(2)*(2*e2**3 - e5**3*x(2)*x(3) - 2*e2**2*e5*(x(1) + 2*x(2) + x(3)) + &
         e2*e5**2*(x(1)*(x(2) + x(3)) + x(2)*(1 + x(2) + 3*x(3)))), 2*Q(2)*Q(4)*Q(10), &
         2*Q(4)*Q(6)*Q(8), Q(6)*Q(10)*x(2), 4*e2*Q(2)*Q(4)*Q(9)*x(2), &
         4*e2*Q(4)*Q(5)*Q(8)*x(2), 4*Q(2)*Q(4)**2*Q(8), x(2)*(4*e2**4 + &
         e5**4*x(1)*x(2)*x(3) - 4*e2**3*e5*(x(1) + 2*x(2) + x(3)) - &
         e2*e5**3*x(2)*(x(1) + 2*x(1)*x(2) + (1 + 4*x(1) + 2*x(2))*x(3)) + &
         e2**2*e5**2*(3*x(2)**2 + 3*x(1)*x(3) + x(2)*(1 + 7*x(1) + 7*x(3)))), &
         2*Q(2)*Q(4)*Q(10)*x(2), 2*Q(4)*Q(6)*Q(8)*x(2), 4*Q(2)*Q(4)**2*Q(8)*x(2)/)
      CASE (26_I2B)
        QQ=Merge(A(5,L), A(6,L), ANY(G(5,L)==G(3:4,L)))
        Q(1:12)=(/-2*e2**3*e4,-2*e2*e4*(e2**2+e4*e5*A(1,L)),4*e2**3*(e2-e5*(QQ+A(1,L))), &
        e2*(4*e2**3+4*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(1+2*QQ+6*A(1,L))), &
        e2**2*(4*e2**2+e5**2*A(1,L)-e2*e5*(1+2*QQ+2*A(1,L))),-2*e2**2*e4*(e2+A(1,L)-e1*A(1,L)), &
        4*e2*(e2**2+e4*e5*A(1,L))*(e2-e5*(QQ+A(1,L))),-2*e2*e4*(e2+A(1,L)-e1*A(1,L))**2, &
        4*e2**2*(e2+QQ-e1*QQ)*(e2+A(1,L)-e1*A(1,L)), &
        e2*(e2+A(1,L)-e1*A(1,L))*(4*e2**2+e5**2*A(1,L)-e2*e5*(1+2*QQ+2*A(1,L))),e2+QQ-e1*QQ, &
        e2+A(1,L)-e1*A(1,L)/)
        coe(L,:)=(/0._DP,Q(1),0._DP,0._DP,Q(1),Q(1),0._DP,0._DP,0._DP,0._DP,&
         0._DP,Q(1),0._DP,Q(1),0._DP,Q(2),Q(3),0._DP,Q(4), &
        4*e2*(e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(QQ+3*A(1,L))),0._DP,0._DP,Q(5),Q(6),0._DP, &
        4*e2**4+11*e2**2*e5**2*A(1,L)-5*e2*e5**3*A(1,L)+e5**4*A(1,L)-e2**3*e5*(1+2*QQ+10*A(1,L)), &
        4*e2**2*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L))),0._DP,Q(5),Q(6),&
        0._DP,0._DP,Q(3),0._DP,0._DP,Q(2),0._DP,Q(4),Q(5),Q(6),Q(7),0._DP,0._DP,Q(8), &
        Q(9),Q(7),(e2**2+e4*e5*A(1,L))*(4*e2**2+e5**2*A(1,L)-e2*e5*(1+2*QQ+2*A(1,L))), &
        4*e2**2*(e2-e5*(QQ+A(1,L)))*Q(12),0._DP,2*e2* &
         (2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L))), &
        (4*e2**3+4*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(1+2*QQ+6*A(1,L)))*Q(12),Q(10), &
        2*(2*e2**4+e5**4*QQ*A(1,L)+2*e2**2*e5**2*A(1,L)*(2+2*QQ+A(1,L))-e2*e5**3*A(1,L)*(1+3*QQ+A(1,L))- &
           2*e2**3*e5*(QQ+3*A(1,L))),Q(9),Q(8), &
           4*e2*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L)))*Q(12),0._DP,0._DP,Q(10), &
        4*(e2**2+e4*e5*A(1,L))*Q(11)*Q(12),4*e2*(e2-e5*(QQ+A(1,L)))*Q(12)**2,0._DP, &
        2*(2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L)))*Q(12), &
        4*e2*Q(11)*Q(12)**2,(4*e2**2+e5**2*A(1,L)-e2*e5*(1+2*QQ+2*A(1,L)))*Q(12)**2,4*Q(11)*Q(12)**3/)
      CASE (27_I2B)
        QQ=Merge(A(5,L), A(6,L), ANY(G(5,L)==G(1:2,L)))
        Q(1:12)=(/-2*e2**3*e4,4*e2**3*(e2-e5*(QQ+A(3,L))),-2*e2*e4*(e2**2+e4*e5*A(3,L)), &
          e2*(4*e2**3+4*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(1+2*QQ+6*A(3,L))), &
          -2*e2**2*e4*(e2+A(3,L)-e1*A(3,L)),e2**2*(4*e2**2+e5**2*A(3,L)-e2*e5*(1+2*QQ+2*A(3,L))), &
          4*e2*(e2**2+e4*e5*A(3,L))*(e2-e5*(QQ+A(3,L))),4*e2**2*(e2+QQ-e1*QQ)*(e2+A(3,L)-e1*A(3,L)), &
          -2*e2*e4*(e2+A(3,L)-e1*A(3,L))**2,e2*(e2+A(3,L)-e1*A(3,L))* &
         (4*e2**2+e5**2*A(3,L)-e2*e5*(1+2*QQ+2*A(3,L))),e2+QQ-e1*QQ,e2+A(3,L)-e1*A(3,L)/)
        coe(L,:)=(/0._DP,Q(1),0._DP,0._DP,Q(1),Q(1),0._DP,&
          0._DP,0._DP,0._DP,0._DP,Q(1),0._DP,Q(1),0._DP,Q(2),Q(3),0._DP,Q(4), &
          4*e2*(e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(QQ+3*A(3,L))),0._DP,0._DP,Q(5),Q(6),0._DP, &
          4*e2**2*(e2**2+e5**2*A(3,L)-e2*e5*(QQ+2*A(3,L))), &
          4*e2**4+11*e2**2*e5**2*A(3,L)-5*e2*e5**3*A(3,L)+ &
          e5**4*A(3,L)-e2**3*e5*(1+2*QQ+10*A(3,L)),0._DP,Q(5),Q(6),0._DP,0._DP, &
          Q(3),0._DP,0._DP,Q(2),0._DP,Q(4),Q(5),Q(6),Q(7),0._DP,0._DP,&
          Q(8),Q(9),Q(7),4*e2**2*(e2-e5*(QQ+A(3,L)))*Q(12), &
          (e2**2+e4*e5*A(3,L))*(4*e2**2+e5**2*A(3,L)-e2*e5*(1+2*QQ+2*A(3,L))),0._DP, &
          (4*e2**3+4*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(1+2*QQ+6*A(3,L)))*Q(12), &
          2*e2*(2*e2**3-e5**3*QQ*A(3,L)+e2*e5**2*A(3,L)*(1+3*QQ+A(3,L))-2*e2**2*e5*(QQ+2*A(3,L))),Q(10), &
          4*e2*(e2**2+e5**2*A(3,L)-e2*e5*(QQ+2*A(3,L)))*Q(12),Q(9),Q(8), &
          2*(2*e2**4+e5**4*QQ*A(3,L)+2*e2**2*e5**2*A(3,L)*(2+ &
          2*QQ+A(3,L))-e2*e5**3*A(3,L)*(1+3*QQ+A(3,L))- &
          2*e2**3*e5*(QQ+3*A(3,L))),0._DP,0._DP,Q(10), &
          4*e2*(e2-e5*(QQ+A(3,L)))*Q(12)**2,4*(e2**2+e4*e5*A(3,L))*Q(11)*Q(12), &
          0._DP,2*(2*e2**3-e5**3*QQ*A(3,L)+e2*e5**2*A(3,L)*(1+ &
          3*QQ+A(3,L))-2*e2**2*e5*(QQ+2*A(3,L)))*Q(12), &
          (4*e2**2+e5**2*A(3,L)-e2*e5*(1+2*QQ+2*A(3,L)))*Q(12)**2,&
          4*e2*Q(11)*Q(12)**2,4*Q(11)*Q(12)**3/)
      CASE (28_I2B)
        X(1:2)=Merge(A(3:4,L), A(4:3:-1,L), ANY(G(3,L)==G(1:2,L)))
        Q(1:12)=(/4*e2**3*e3*A(5,L),4*e2**3*A(5,L)*(e2-e5*(A(5,L)+x(1))),&
        4*e2**2*e3*A(5,L)*(e2-e5*(x(1)+x(2))), &
        e2*A(5,L)*(4*e2**3-e5**3*x(1)-2*e2**2*e5*(1+A(5,L)+2*x(1)+x(2))+e2*e5**2*(4*x(1)+x(2))), &
        2*e2**2*e3*A(5,L)*(2*e2-e5*(x(1)+x(2))),&
        2*e2**2*A(5,L)*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1))), &
        4*e2**2*A(5,L)*(e2+A(5,L)-e1*A(5,L))*(e2+x(1)-e1*x(1)), &
        4*e2*e3*A(5,L)*(e2+x(1)-e1*x(1))*(e2+x(2)-e1*x(2)), &
        e2*A(5,L)*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1)))*(2*e2-e5*(x(1)+x(2))),&
        e2+A(5,L)-e1*A(5,L),e2+x(1)-e1*x(1),e2+x(2)-e1*x(2)/)
        coe(L,:)=(/4*e2**3*e3,Q(1),2*e2**2*e3*(2*e2-e5*(x(1)+x(2))),&
        2*e2**2*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1))), &
        Q(1),Q(1),e2*(4*e2**3-e5**3*x(1)-2*e2**2*e5*(1+A(5,L)+2*x(1)+x(2))+e2*e5**2*(4*x(1)+x(2))), &
        4*e2**2*(e2**2+e5**2*x(1)-e2*e5*(A(5,L)+2*x(1)+x(2))),4*e2**2*e3*(e2-e5*(x(1)+x(2))), &
        4*e2**3*(e2-e5*(A(5,L)+x(1))),2*e2**2*(2*e2**2+e5**2*x(1)-e2*e5*(2*A(5,L)+3*x(1)+x(2))),Q(1), &
        2*e2*(2*e2**3-e5**3*x(1)+e2*e5**2*(3*x(1)+x(2))-e2**2*e5*(1+A(5,L)+3*x(1)+2*x(2))),Q(1), &
        4*e2**2*(e2-e5*(A(5,L)+x(1)))*(e2-e5*(x(1)+x(2))),Q(2),Q(3), &
        2*e2*(2*e2**3-e5**3*x(1)*(A(5,L)+x(2))-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
           e2*e5**2*(A(5,L)*(3*x(1)+x(2))+x(1)*(1+x(1)+3*x(2)))),Q(4), &
        4*e2**2*A(5,L)*(e2**2+e5**2*x(1)-e2*e5*(A(5,L)+2*x(1)+x(2))),&
        4*e2*e3*Q(11)*Q(12),4*e2**2*Q(10)*Q(11),Q(5), &
        Q(6),e2*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1)))*(2*e2-e5*(x(1)+x(2))), &
        2*e2**2*A(5,L)*(2*e2**2+e5**2*x(1)-e2*e5*(2*A(5,L)+3*x(1)+x(2))), &
        2*e2*A(5,L)*(2*e2**3-e5**3*x(1)+e2*e5**2*(3*x(1)+x(2))-e2**2*e5*(1+A(5,L)+3*x(1)+2*x(2))), &
        4*e2**4+e5**4*x(1)*x(2)-2*e2**3*e5*(1+A(5,L)+3*x(1)+2*x(2))-e2*e5**3*x(1)*(1+x(1)+4*x(2))+ &
         e2**2*e5**2*(x(1)**2+(2+A(5,L))*x(2)+x(1)*(5+A(5,L)+5*x(2))),Q(5),Q(6), &
        e2*(4*e2**3-e5**3*A(5,L)*x(1)-2*e2**2*e5*(2*A(5,L)+3*x(1)+x(2))+ &
           e2*e5**2*(x(1)*(1+x(1)+x(2))+A(5,L)*(5*x(1)+x(2)))), &
        2*e2*(2*e2**3-e5**3*x(1)*x(2)-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
           e2*e5**2*(A(5,L)*(x(1)+x(2))+x(1)*(1+x(1)+3*x(2)))),Q(3), &
        2*e2*(2*e2**3-e5**3*A(5,L)*x(1)-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
           e2*e5**2*(x(1)*(1+x(1)+x(2))+A(5,L)*(3*x(1)+x(2)))), &
        2*e2*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1)))*(e2-e5*(x(1)+x(2))),Q(2), &
        2*e2**2*(e2-e5*(A(5,L)+x(1)))*(2*e2-e5*(x(1)+x(2))),Q(4),Q(5),Q(6), &
        4*e2**2*A(5,L)*(e2-e5*(A(5,L)+x(1)))*(e2-e5*(x(1)+x(2))),&
        4*e2*Q(11)*Q(12)*(e2-e5*(A(5,L)+x(1))), &
        4*e2*Q(10)*Q(11)*(e2-e5*(x(1)+x(2))),Q(7),Q(8), &
        2*e2*A(5,L)*(2*e2**3-e5**3*x(1)*(A(5,L)+x(2))-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
           e2*e5**2*(A(5,L)*(3*x(1)+x(2))+x(1)*(1+x(1)+3*x(2)))), &
        2*e2**2*A(5,L)*(e2-e5*(A(5,L)+x(1)))*(2*e2-e5*(x(1)+x(2))), &
        2*e2*A(5,L)*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1)))*(e2-e5*(x(1)+x(2))), &
        4*e2**4+e5**4*A(5,L)*x(1)*x(2)-4*e2**3*e5*(A(5,L)+2*x(1)+x(2))- &
         e2*e5**3*x(1)*(A(5,L)+2*A(5,L)*x(1)+(1+4*A(5,L)+2*x(1))*x(2))+ &
         e2**2*e5**2*(3*x(1)**2+3*A(5,L)*x(2)+x(1)*(1+7*A(5,L)+7*x(2))), &
        A(5,L)*(4*e2**4+e5**4*x(1)*x(2)-2*e2**3*e5*(1+ &
        A(5,L)+3*x(1)+2*x(2))-e2*e5**3*x(1)*(1+x(1)+4*x(2))+ &
           e2**2*e5**2*(x(1)**2+(2+A(5,L))*x(2)+x(1)*(5+A(5,L)+5*x(2)))), &
        e2*A(5,L)*(4*e2**3-e5**3*A(5,L)*x(1)-2*e2**2*e5*(2*A(5,L)+3*x(1)+x(2))+ &
           e2*e5**2*(x(1)*(1+x(1)+x(2))+A(5,L)*(5*x(1)+x(2)))),Q(9), &
        2*e2*A(5,L)*(2*e2**3-e5**3*x(1)*x(2)-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
           e2*e5**2*(A(5,L)*(x(1)+x(2))+x(1)*(1+x(1)+3*x(2)))),Q(8),Q(7), &
        2*e2*A(5,L)*(2*e2**3-e5**3*A(5,L)*x(1)-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
           e2*e5**2*(x(1)*(1+x(1)+x(2))+A(5,L)*(3*x(1)+x(2)))), &
        2*Q(11)*Q(12)*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1))),&
        2*e2*Q(10)*Q(11)*(2*e2-e5*(x(1)+x(2))),Q(9), &
        4*e2*A(5,L)*Q(11)*Q(12)*(e2-e5*(A(5,L)+x(1))),&
        4*e2*A(5,L)*Q(10)*Q(11)*(e2-e5*(x(1)+x(2))), &
        4*Q(10)*Q(11)**2*Q(12),A(5,L)*(4*e2**4+e5**4*A(5,L)* &
        x(1)*x(2)-4*e2**3*e5*(A(5,L)+2*x(1)+x(2))- &
           e2*e5**3*x(1)*(A(5,L)+2*A(5,L)*x(1)+(1+4*A(5,L)+2*x(1))*x(2))+ &
           e2**2*e5**2*(3*x(1)**2+3*A(5,L)*x(2)+x(1)*(1+7*A(5,L)+7*x(2)))), &
        2*A(5,L)*Q(11)*Q(12)*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1))), &
        2*e2*A(5,L)*Q(10)*Q(11)*(2*e2-e5*(x(1)+x(2))),4*A(5,L)*Q(10)*Q(11)**2*Q(12)/)
      CASE (29_I2B)
        X(1:2)=Merge(A(3:4,L), A(4:3:-1,L), ANY(G(3,L)==G(5:6,L)))
        Q(1:6)=(/-2*e2**2*e3*e4,-2*e2*e4*(e2**2+e4*e5*A(1,L)),&
          -2*e2*e4*(e2+A(1,L)-e1*A(1,L))**2,e2+A(1,L)-e1*A(1,L), &
          e2+x(1)-e1*x(1),e2+x(2)-e1*x(2)/)
        coe(L,:)=(/0._DP,Q(1),0._DP,0._DP,-(e2*e4*(2*e2**2-2*e2*e5+e5**2)),&
          Q(1),0._DP,0._DP,0._DP,0._DP,0._DP,-2*e2*e4*(e2-e1*e2+e2**2+e5**2),&
          0._DP,Q(1),0._DP,Q(2),4*e2**2*e3*(e2-e5*(x(1)+x(2))),0._DP,e2* &
         (4*e2**3-e5**3*(A(1,L)+x(2))+e2*e5**2*(1+3*A(1,L)+x(1)+2*x(2))- &
          e2**2*e5*(3+4*A(1,L)+2*x(1)+2*x(2))), &
          4*e2**2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+x(1)+x(2))),0._DP,0._DP, &
          e2*(4*e2**3-e5**3*(x(1)+x(2))-e2**2*e5*(5+2*x(1)+2*x(2))+ &
          e2*e5**2*(2+2*x(1)+3*x(2))),4*e2**2*e3*Q(4),0._DP, &
          e2*(4*e2**3-e5**3*A(1,L)+e2*e5**2*(6*A(1,L)+ &
          x(2))-e2**2*e5*(1+8*A(1,L)+2*x(1)+2*x(2))), &
          -2*e2**2*e4*(e2-e5*(A(1,L)+x(1)+x(2))),0._DP, &
          4*e2**4+e5**4*x(2)-e2**3*e5*(5+2*x(1)+2*x(2))-e2*e5**3*(1+x(1)+3*x(2))+ &
           e2**2*e5**2*(4+2*x(1)+3*x(2)),2*e2*(2*e2**2-2*e2*e5+e5**2)*Q(4),0._DP,0._DP, &
          2*e2*(2*e2**2-2*e2*e5+e5**2)*(e2-e5*(x(1)+x(2))),0._DP,0._DP,Q(2),0._DP, &
          e2**2*(4*e2**2+e5**2*(3*A(1,L)+x(1)+2*x(2))-e2*e5*(3+4*A(1,L)+2*x(1)+2*x(2))), &
          e2*e3*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(1)+2*x(2))),e2*e4**2*Q(4), &
          4*e2*(e2**2+e4*e5*A(1,L))*(e2-e5*(x(1)+x(2))),0._DP,0._DP,Q(3),4*e2*e3*Q(5)*Q(6), &
          4*e2**2*(e2-e5*(A(1,L)+x(1)))*(e2-e5*(A(1,L)+x(2))), &
          (e2**2+e4*e5*A(1,L))*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(1)+2*x(2))),&
          -2*e2*e4*Q(4)*(e2-e5*(x(1)+x(2))),0._DP, &
          -(e2*e4*(2*e2**2-2*e2*e5*(A(1,L)+x(1)+x(2))+e5**2*(2*x(1)*x(2)+A(1,L)*(x(1)+x(2))))), &
          e2*Q(4)*(4*e2**2+e5**2*(A(1,L)+x(2))-e2*e5*(1+4*A(1,L)+2*x(1)+2*x(2))), &
          Q(4)*(4*e2**3-e5**3*x(2)+e2*e5**2*(1+x(1)+2*x(2))-e2**2*e5*(3+2*x(1)+2*x(2))), &
          2*e2*(2*e2**3-e5**3*A(1,L)*(x(1)+x(2))-2*e2**2*e5*(2*A(1,L)+x(1)+x(2))+ &
          2*e2*e5**2*(x(1)*x(2)+A(1,L)*(1+x(1)+x(2)))),2*(2*e2**2-2*e2*e5+e5**2)*Q(5)*Q(6),Q(3), &
          4*e2**2*Q(4)*(e2-e5*(A(1,L)+x(1)+x(2))),0._DP,0._DP, &
          e2*Q(4)*(4*e2**2+e5**2*(x(1)+2*x(2))-e2*e5*(3+2*x(1)+2*x(2))),&
          4*(e2**2+e4*e5*A(1,L))*Q(5)*Q(6), &
          4*e2*Q(4)**2*(e2-e5*(x(1)+x(2))),0._DP,2*e2*Q(4)* &
          (2*e2**2-2*e2*e5*(A(1,L)+x(1)+x(2))+e5**2*(2*x(1)*x(2)+A(1,L)*(x(1)+x(2)))),&
          -2*e4*Q(4)*Q(5)*Q(6),Q(4)**2*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(1)+2*x(2))),&
          4*Q(4)**2*Q(5)*Q(6)/)
      CASE (30_I2B)
        X(1:2)=Merge(A(1:2,L), A(2:1:-1,L), ANY(G(1,L)==G(5:6,L)))
        Q(1:6)=(/-2*e2**2*e3*e4,-2*e2*e4*(e2**2+e4*e5*A(3,L)), &
          -2*e2*e4*(e2+A(3,L)-e1*A(3,L))**2,e2+A(3,L)-e1*A(3,L), &
          e2+x(1)-e1*x(1),e2+x(2)-e1*x(2)/)
        coe(L,:)=(/0._DP,Q(1),0._DP,0._DP,Q(1),-(e2*e4*(2*e2**2-2*e2*e5+e5**2)),&
          0._DP,0._DP,0._DP,0._DP,0._DP,-2*e2*e4*(e2-e1*e2+e2**2+e5**2),0._DP,Q(1),0._DP, &
          4*e2**2*e3*(e2-e5*(x(1)+x(2))),Q(2),0._DP,e2* &
          (4*e2**3-e5**3*(A(3,L)+x(2))+e2*e5**2*(1+3*A(3,L)+x(1)+2*x(2))- &
          e2**2*e5*(3+4*A(3,L)+2*x(1)+2*x(2))), &
          4*e2**2*(e2**2+e5**2*A(3,L)-e2*e5*(2*A(3,L)+x(1)+x(2))),0._DP,0._DP, &
          4*e2**2*e3*Q(4),e2*(4*e2**3-e5**3*(x(1)+x(2))- &
          e2**2*e5*(5+2*x(1)+2*x(2))+e2*e5**2*(2+2*x(1)+3*x(2))), &
          0._DP,-2*e2**2*e4*(e2-e5*(A(3,L)+x(1)+x(2))), &
          e2*(4*e2**3-e5**3*A(3,L)+e2*e5**2*(6*A(3,L)+x(2))- &
          e2**2*e5*(1+8*A(3,L)+2*x(1)+2*x(2))),0._DP, &
          2*e2*(2*e2**2-2*e2*e5+e5**2)*Q(4), &
          4*e2**4+e5**4*x(2)-e2**3*e5*(5+2*x(1)+2*x(2))- &
          e2*e5**3*(1+x(1)+3*x(2))+e2**2*e5**2*(4+2*x(1)+3*x(2)),0._DP,0._DP,Q(2),0._DP,0._DP, &
          2*e2*(2*e2**2-2*e2*e5+e5**2)*(e2-e5*(x(1)+x(2))),0._DP, &
          e2**2*(4*e2**2+e5**2*(3*A(3,L)+x(1)+2*x(2))-e2*e5*(3+4*A(3,L)+2*x(1)+2*x(2))),&
          e2*e4**2*Q(4),e2*e3*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(1)+2*x(2))),&
          4*e2*(e2**2+e4*e5*A(3,L))*(e2-e5*(x(1)+x(2))),0._DP,0._DP, &
          4*e2*e3*Q(5)*Q(6),Q(3),4*e2**2*(e2-e5*(A(3,L)+x(1)))*(e2-e5*(A(3,L)+x(2))), &
          -2*e2*e4*Q(4)*(e2-e5*(x(1)+x(2))),(e2**2+ &
          e4*e5*A(3,L))*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(1)+2*x(2))),0._DP, &
          e2*Q(4)*(4*e2**2+e5**2*(A(3,L)+x(2))-e2*e5*(1+4*A(3,L)+2*x(1)+2*x(2))), &
          -(e2*e4*(2*e2**2-2*e2*e5*(A(3,L)+x(1)+x(2))+e5**2*(2*x(1)*x(2)+A(3,L)*(x(1)+x(2))))), &
          Q(4)*(4*e2**3-e5**3*x(2)+e2*e5**2*(1+x(1)+2*x(2))-e2**2*e5*(3+2*x(1)+2*x(2))), &
          4*e2**2*Q(4)*(e2-e5*(A(3,L)+x(1)+x(2))),Q(3),2*(2*e2**2-2*e2*e5+e5**2)*Q(5)*Q(6), &
          2*e2*(2*e2*Q(5)*Q(6)+e4*e5*A(3,L)*(2*e2-e5*(x(1)+x(2)))),0._DP,0._DP, &
          e2*Q(4)*(4*e2**2+e5**2*(x(1)+2*x(2))-e2*e5*(3+2*x(1)+2*x(2))),&
          4*e2*Q(4)**2*(e2-e5*(x(1)+x(2))), &
          4*(e2**2+e4*e5*A(3,L))*Q(5)*Q(6),0._DP,2*e2*Q(4)* &
          (2*e2**2-2*e2*e5*(A(3,L)+x(1)+x(2))+e5**2*(2*x(1)*x(2)+A(3,L)*(x(1)+x(2)))), &
          Q(4)**2*(4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(1)+2*x(2))),&
          -2*e4*Q(4)*Q(5)*Q(6),4*Q(4)**2*Q(5)*Q(6)/)
      CASE (31_I2B)
        X(1:2)=Merge(A(1:2,L), A(2:1:-1,L), ANY(G(1,L)==G(3:4,L)))
        Q(1:12)=(/4*e2**3*e3*A(5,L),4*e2**2*e3*A(5,L)*(e2-e5*(x(1)+x(2))), &
          4*e2**3*A(5,L)*(e2-e5*(A(5,L)+x(1))),e2*A(5,L)*(4*e2**3-e5**3*x(1)- &
          2*e2**2*e5*(1+A(5,L)+2*x(1)+x(2))+e2*e5**2*(4*x(1)+x(2))), &
          2*e2**2*A(5,L)*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1))), &
          2*e2**2*e3*A(5,L)*(2*e2-e5*(x(1)+x(2))), &
          4*e2*e3*A(5,L)*(e2+x(1)-e1*x(1))*(e2+x(2)-e1*x(2)), &
          4*e2**2*A(5,L)*(e2+A(5,L)-e1*A(5,L))*(e2+x(1)-e1*x(1)), &
          e2*A(5,L)*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1)))*(2*e2-e5*(x(1)+x(2))),&
          e2+A(5,L)-e1*A(5,L),e2+x(1)-e1*x(1),e2+x(2)-e1*x(2)/)
        coe(L,:)=(/4*e2**3*e3,Q(1),2*e2**2*(2*e2**2+e5**2*x(1)- &
          e2*e5*(1+A(5,L)+x(1))),2*e2**2*e3*(2*e2-e5*(x(1)+x(2))), &
          Q(1),Q(1),e2*(4*e2**3-e5**3*x(1)-2*e2**2*e5*(1+A(5,L)+ &
          2*x(1)+x(2))+e2*e5**2*(4*x(1)+x(2))), &
          4*e2**2*(e2**2+e5**2*x(1)-e2*e5*(A(5,L)+2*x(1)+x(2))),4*e2**3*(e2-e5*(A(5,L)+x(1))), &
          4*e2**2*e3*(e2-e5*(x(1)+x(2))),2*e2*(2*e2**3-e5**3*x(1)+e2*e5**2*(3*x(1)+x(2))- &
          e2**2*e5*(1+A(5,L)+3*x(1)+2*x(2))),Q(1), &
          2*e2**2*(2*e2**2+e5**2*x(1)-e2*e5*(2*A(5,L)+3*x(1)+x(2))), &
          Q(1),4*e2**2*(e2-e5*(A(5,L)+x(1)))*(e2-e5*(x(1)+x(2))),Q(2),Q(3), &
          2*e2*(2*e2**3-e5**3*x(1)*(A(5,L)+x(2))-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
          e2*e5**2*(A(5,L)*(3*x(1)+x(2))+x(1)*(1+x(1)+3*x(2)))),Q(4), &
          4*e2**2*A(5,L)*(e2**2+e5**2*x(1)-e2*e5*(A(5,L)+2*x(1)+x(2))),&
          4*e2**2*Q(10)*Q(11),4*e2*e3*Q(11)*Q(12),Q(5),Q(6), &
          e2*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1)))*(2*e2-e5*(x(1)+x(2))), &
          2*e2*A(5,L)*(2*e2**3-e5**3*x(1)+e2*e5**2*(3*x(1)+x(2))-e2**2*e5*(1+A(5,L)+3*x(1)+2*x(2))), &
          2*e2**2*A(5,L)*(2*e2**2+e5**2*x(1)-e2*e5*(2*A(5,L)+3*x(1)+x(2))), &
          e2*(4*e2**3-e5**3*A(5,L)*x(1)-2*e2**2*e5*(2*A(5,L)+3*x(1)+x(2))+ &
          e2*e5**2*(x(1)*(1+x(1)+x(2))+A(5,L)*(5*x(1)+x(2)))),Q(5),Q(6), &
          4*e2**4+e5**4*x(1)*x(2)-2*e2**3*e5*(1+A(5,L)+3*x(1)+2*x(2))-e2*e5**3*x(1)*(1+x(1)+4*x(2))+ &
          e2**2*e5**2*(x(1)**2+(2+A(5,L))*x(2)+x(1)*(5+A(5,L)+5*x(2))), &
          2*e2*(2*e2**3-e5**3*A(5,L)*x(1)-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
          e2*e5**2*(x(1)*(1+x(1)+x(2))+A(5,L)*(3*x(1)+x(2)))),Q(3), &
          2*e2*(2*e2**3-e5**3*x(1)*x(2)-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
          e2*e5**2*(A(5,L)*(x(1)+x(2))+x(1)*(1+x(1)+3*x(2)))), &
          2*e2**2*(e2-e5*(A(5,L)+x(1)))*(2*e2-e5*(x(1)+x(2))),Q(2), &
          2*e2*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1)))*(e2-e5*(x(1)+x(2))),Q(4),Q(5),Q(6), &
          4*e2**2*A(5,L)*(e2-e5*(A(5,L)+x(1)))*(e2-e5*(x(1)+x(2))), &
          4*e2*Q(10)*Q(11)*(e2-e5*(x(1)+x(2))), &
          4*e2*Q(11)*Q(12)*(e2-e5*(A(5,L)+x(1))),Q(7),Q(8), &
          2*e2*A(5,L)*(2*e2**3-e5**3*x(1)*(A(5,L)+x(2))-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
          e2*e5**2*(A(5,L)*(3*x(1)+x(2))+x(1)*(1+x(1)+3*x(2)))), &
          2*e2*A(5,L)*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1)))*(e2-e5*(x(1)+x(2))), &
          2*e2**2*A(5,L)*(e2-e5*(A(5,L)+x(1)))*(2*e2-e5*(x(1)+x(2))), &
          4*e2**4+e5**4*A(5,L)*x(1)*x(2)-4*e2**3*e5*(A(5,L)+2*x(1)+x(2))- &
          e2*e5**3*x(1)*(A(5,L)+2*A(5,L)*x(1)+(1+4*A(5,L)+2*x(1))*x(2))+ &
          e2**2*e5**2*(3*x(1)**2+3*A(5,L)*x(2)+x(1)*(1+7*A(5,L)+7*x(2))), &
          e2*A(5,L)*(4*e2**3-e5**3*A(5,L)*x(1)-2*e2**2*e5*(2*A(5,L)+3*x(1)+x(2))+ &
          e2*e5**2*(x(1)*(1+x(1)+x(2))+A(5,L)*(5*x(1)+x(2)))), &
          A(5,L)*(4*e2**4+e5**4*x(1)*x(2)-2*e2**3*e5*(1+ &
          A(5,L)+3*x(1)+2*x(2))-e2*e5**3*x(1)*(1+x(1)+4*x(2))+ &
          e2**2*e5**2*(x(1)**2+(2+A(5,L))*x(2)+x(1)*(5+A(5,L)+5*x(2)))),Q(9), &
          2*e2*A(5,L)*(2*e2**3-e5**3*A(5,L)*x(1)-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
          e2*e5**2*(x(1)*(1+x(1)+x(2))+A(5,L)*(3*x(1)+x(2)))),Q(8),Q(7), &
          2*e2*A(5,L)*(2*e2**3-e5**3*x(1)*x(2)-2*e2**2*e5*(A(5,L)+2*x(1)+x(2))+ &
          e2*e5**2*(A(5,L)*(x(1)+x(2))+x(1)*(1+x(1)+3*x(2)))),2*e2*Q(10)*Q(11)*(2*e2-e5*(x(1)+x(2))), &
          2*Q(11)*Q(12)*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1))), &
          Q(9),4*e2*A(5,L)*Q(10)*Q(11)*(e2-e5*(x(1)+x(2))), &
          4*e2*A(5,L)*Q(11)*Q(12)*(e2-e5*(A(5,L)+x(1))),4*Q(10)*Q(11)**2*Q(12), &
          A(5,L)*(4*e2**4+e5**4*A(5,L)*x(1)*x(2)-4*e2**3*e5*(A(5,L)+2*x(1)+x(2))- &
          e2*e5**3*x(1)*(A(5,L)+2*A(5,L)*x(1)+(1+4*A(5,L)+2*x(1))*x(2))+ &
          e2**2*e5**2*(3*x(1)**2+3*A(5,L)*x(2)+x(1)*(1+7*A(5,L)+7*x(2)))), &
          2*e2*A(5,L)*Q(10)*Q(11)*(2*e2-e5*(x(1)+x(2))), &
          2*A(5,L)*Q(11)*Q(12)*(2*e2**2+e5**2*x(1)-e2*e5*(1+A(5,L)+x(1))), &
          4*A(5,L)*Q(10)*Q(11)**2*Q(12)/)
      CASE (32_I2B)
        QQ=Merge(A(4,L), A(3,L), G(1,L)==G(3,L))
        Q(1:13)=(/2*e2**4*A(5,L),2*e2*(e2**2+e4*e5*A(1,L))*(e2-e5*(QQ+A(1,L))),&
          2*e2**2*(e2**2+e4*e5*A(1,L))*A(5,L),2*e2**3*(e2-e5*(QQ+A(1,L)))*A(5,L), &
          e2**2*(2*e2**2+e5**2*A(1,L)-e2*e5*(QQ+3*A(1,L)))*A(5,L), &
          e2**3*(2*e2-e5*(QQ+A(1,L)))*A(5,L),2*e2**3*(e2+A(1,L)-e1*A(1,L))*A(5,L), &
          2*e2*(e2**2+e4*e5*A(1,L))*(e2-e5*(QQ+A(1,L)))*A(5,L),&
          2*e2**2*(e2+A(1,L)-e1*A(1,L))**2*A(5,L), &
          2*e2**2*(e2+QQ-e1*QQ)*(e2+A(1,L)-e1*A(1,L))*A(5,L), &
          e2**2*(e2+A(1,L)-e1*A(1,L))*(2*e2-e5*(QQ+A(1,L)))*A(5,L),&
          e2+QQ-e1*QQ,e2+A(1,L)-e1*A(1,L)/)
        coe(L,:)=(/2*e2**4,Q(1),e2**3*(2*e2-e5*(QQ+A(1,L))),2*e2**3*Q(13),Q(1),Q(1), &
          e2**2*(2*e2**2+e5**2*A(1,L)-e2*e5*(QQ+3*A(1,L))), &
          2*e2*(e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(QQ+3*A(1,L))),&
          2*e2**3*(e2-e5*(QQ+A(1,L))),2*e2**2*(e2**2+e4*e5*A(1,L)),&
          e2*(2*e2**3+4*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(QQ+5*A(1,L))),Q(1), &
          2*e2**2*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L))),Q(1),Q(2),Q(3),Q(4),Q(2),Q(5), &
          2*e2*(e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(QQ+3*A(1,L)))*A(5,L),&
          2*e2**2*Q(12)*Q(13),2*e2**2*Q(13)**2,Q(6),Q(7),e2**2*(2*e2-e5*(QQ+A(1,L)))*Q(13), &
          e2*(2*e2**3+4*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(QQ+5*A(1,L)))*A(5,L), &
          2*e2**2*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L)))*A(5,L), &
          e2*(2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L))),&
          Q(6),Q(7),e2*(2*e2**2+e5**2*A(1,L)-e2*e5*(QQ+3*A(1,L)))*Q(13), &
          2*e2**4+e5**4*QQ*A(1,L)+2*e2**2*e5**2*A(1,L)*(2+2*QQ+A(1,L))- &
          e2*e5**3*A(1,L)*(1+3*QQ+A(1,L))-2*e2**3*e5*(QQ+3*A(1,L)),Q(4), &
          2*e2*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L)))*Q(13), &
          2*e2**2*(e2-e5*(QQ+A(1,L)))*Q(13),Q(3), &
          e2*(e2**2+e4*e5*A(1,L))*(2*e2-e5*(QQ+A(1,L))),Q(5),Q(6),Q(7),Q(8), &
          2*(e2**2+e4*e5*A(1,L))*Q(12)*Q(13),2*e2*(e2-e5*(QQ+A(1,L)))*Q(13)**2,Q(9),Q(10),Q(8), &
          e2*(e2**2+e4*e5*A(1,L))*(2*e2-e5*(QQ+A(1,L)))*A(5,L), &
          2*e2**2*(e2-e5*(QQ+A(1,L)))*A(5,L)*Q(13), &
          (2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L)))*Q(13), &
          e2*(2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+3*QQ+A(1,L))-2*e2**2*e5*(QQ+2*A(1,L)))*A(5,L), &
          e2*(2*e2**2+e5**2*A(1,L)-e2*e5*(QQ+3*A(1,L)))*A(5,L)*Q(13),Q(11), &
          (2*e2**4+e5**4*QQ*A(1,L)+2*e2**2*e5**2*A(1,L)*(2+2*QQ+A(1,L))-e2*e5**3* &
          A(1,L)*(1+3*QQ+A(1,L))-2*e2**3*e5*(QQ+3*A(1,L)))*A(5,L),Q(10),Q(9), &
          2*e2*(e2**2+e5**2*A(1,L)-e2*e5*(QQ+2*A(1,L)))*A(5,L)*Q(13), &
          2*e2*Q(12)*Q(13)**2,e2*(2*e2-e5*(QQ+A(1,L)))*Q(13)**2, &
          Q(11),2*(e2**2+e4*e5*A(1,L))*A(5,L)*Q(12)*Q(13), &
          2*e2*(e2-e5*(QQ+A(1,L)))*A(5,L)*Q(13)**2,2*Q(12)*Q(13)**3, &
          (2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+3*QQ+A(1,L))- &
          2*e2**2*e5*(QQ+2*A(1,L)))*A(5,L)*Q(13),2*e2*A(5,L)*Q(12)*Q(13)**2, &
          e2*(2*e2-e5*(QQ+A(1,L)))*A(5,L)*Q(13)**2,2*A(5,L)*Q(12)*Q(13)**3/)
      CASE (33_I2B)
        Q(1:10)=(/e2**2*(2*e2**2-2*e2*e5+e5**2),2*e2**3*e3, &
          2*e2**2*(e2**2+e4*e5*A(1,L)),-(e2**2*e4*(e2-e5*(A(1,L)+A(3,L)))), &
          e2*(2*e2**2-2*e2*e5+e5**2)*(e2+A(3,L)-e1*A(3,L)), &
          -(e2**2*e4*(e2+A(1,L)-e1*A(1,L))),2*e2**2*(e2+A(1,L)-e1*A(1,L))**2, &
          -(e2*e4*(e2+A(1,L)-e1*A(1,L))*(e2+A(3,L)-e1*A(3,L))), &
          e2+A(3,L)-e1*A(3,L),e2+A(1,L)-e1*A(1,L)/)
        coe(L,:)=(/0._DP,Q(1),0._DP,0._DP,Q(1),Q(2),0._DP,0._DP,0._DP,0._DP,&
          0._DP,Q(1),0._DP,Q(2),0._DP,Q(3),2*e2*e3*(e2**2+e4*e5*A(3,L)),0._DP,Q(4), &
          2*e2**2*(e2**2-2*e2*e5*(A(1,L)+A(3,L))+e5**2*(A(1,L)+A(3,L))),0._DP, &
          0._DP,Q(5),Q(6),0._DP,2*e2**2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L))), &
          -(e2*e4*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L)))),0._DP,Q(5),Q(6),0._DP,0._DP, &
          (2*e2**2-2*e2*e5+e5**2)*(e2**2+e4*e5*A(3,L)),0._DP,&
          0._DP,Q(3),0._DP,Q(4),2*e2**2*e3*Q(9),Q(6), &
          2*(e2**2+e4*e5*A(1,L))*(e2**2+e4*e5*A(3,L)),0._DP,0._DP,Q(7),&
          2*e2*e3*Q(9)**2,2*e2**2*(e2-e5*(A(1,L)+A(3,L)))**2, &
          2*e2*(e2**2+e4*e5*A(1,L))*Q(9),-(e4*(e2**2+e4*e5*A(3,L))*Q(10)),&
          0._DP,-(e2*e4*(e2-e5*(A(1,L)+A(3,L)))*Q(9)), &
          2*e2**2*(e2-e5*(A(1,L)+A(3,L)))*Q(10),Q(8),&
          2*e2*(e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)))*Q(9), &
          (2*e2**2-2*e2*e5+e5**2)*Q(9)**2,Q(7), &
          2*e2*(e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L)))*Q(10),0._DP,0._DP,Q(8), &
          2*(e2**2+e4*e5*A(1,L))*Q(9)**2,2*(e2**2+e4*e5*A(3,L))*Q(10)**2,&
          0._DP,2*e2*(e2-e5*(A(1,L)+A(3,L)))*Q(9)*Q(10), &
          -(e4*Q(9)**2*Q(10)),2*e2*Q(9)*Q(10)**2,2*Q(9)**2*Q(10)**2/)
      CASE (34_I2B)
        QQ=Merge(A(2,L), A(1,L), G(1,L)==G(3,L))
        Q(1:22)=(/2*e2**4*A(5,L),2*e2*(e2**2+e4*e5*A(3,L))*(e2-e5*(QQ+A(3,L))),&
          2*e2**3*(e2-e5*(QQ+A(3,L)))*A(5,L),2*e2**2*(e2**2+e4*e5*A(3,L))*A(5,L),&
          e2**2*(2*e2**2+e5**2*A(3,L)-e2*e5*(QQ+3*A(3,L)))*A(5,L), &
          2*e2**3*(e2+A(3,L)-e1*A(3,L))*A(5,L),e2**3*(2*e2-e5*(QQ+A(3,L)))*A(5,L), &
          2*e2*(e2**2+e4*e5*A(3,L))*(e2-e5*(QQ+A(3,L)))*A(5,L), &
          2*e2**2*(e2+QQ-e1*QQ)*(e2+A(3,L)-e1*A(3,L))*A(5,L), &
          2*e2**2*(e2+A(3,L)-e1*A(3,L))**2*A(5,L), &
          e2**2*(e2+A(3,L)-e1*A(3,L))*(2*e2-e5*(QQ+A(3,L)))*A(5,L), &
          e2+A(3,L)-e1*A(3,L),e2+QQ-e1*QQ,e2-e5*(QQ+A(3,L)), &
          e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(QQ+3*A(3,L)),&
          2*e2-e5*(QQ+A(3,L)),e2**2+e4*e5*A(3,L), &
          e2**2+e5**2*A(3,L)-e2*e5*(QQ+2*A(3,L)), &
          2*e2**3-e5**3*QQ*A(3,L)+e2*e5**2*A(3,L)*(1+3*QQ+A(3,L))-2*e2**2*e5*(QQ+2*A(3,L)), &
          2*e2**4+e5**4*QQ*A(3,L)+2*e2**2*e5**2*A(3,L)*(2+2*QQ+ &
          A(3,L))-e2*e5**3*A(3,L)*(1+3*QQ+A(3,L))-2*e2**3*e5*(QQ+3*A(3,L)),&
          2*e2**2+e5**2*A(3,L)-e2*e5*(QQ+3*A(3,L)), &
          2*e2**3+4*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(QQ+5*A(3,L))/)
        coe(L,:)=(/2*e2**4,Q(1),2*e2**3*Q(12),e2**3*Q(16),Q(1),Q(1), &
          e2**2*Q(21),2*e2*Q(15),2*e2**2*Q(17),2*e2**3*Q(14), &
          2*e2**2*Q(18),Q(1),e2*Q(22),Q(1),Q(2),Q(3),Q(4),Q(2),Q(5),&
          2*e2*A(5,L)*Q(15),2*e2**2*Q(12)**2,2*e2**2*Q(12)*Q(13), &
          Q(6),Q(7),e2**2*Q(12)*Q(16),2*e2**2*A(5,L)*Q(18),e2*A(5,L)*Q(22),&
          e2*Q(12)*Q(21),Q(6),Q(7),e2*Q(19),2*e2*Q(12)*Q(18), &
          Q(4),Q(20),e2*Q(16)*Q(17),Q(3),2*e2**2*Q(12)*Q(14),Q(5),Q(6),&
          Q(7),Q(8),2*e2*Q(12)**2*Q(14),2*Q(12)*Q(13)*Q(17),Q(9), &
          Q(10),Q(8),2*e2**2*A(5,L)*Q(12)*Q(14),e2*A(5,L)*Q(16)*Q(17),&
          Q(12)*Q(19),e2*A(5,L)*Q(12)*Q(21),e2*A(5,L)*Q(19),Q(11), &
          2*e2*A(5,L)*Q(12)*Q(18),Q(10),Q(9),A(5,L)*Q(20),&
          e2*Q(12)**2*Q(16),2*e2*Q(12)**2*Q(13),Q(11), &
          2*e2*A(5,L)*Q(12)**2*Q(14),2*A(5,L)*Q(12)*Q(13)*Q(17),&
          2*Q(12)**3*Q(13),A(5,L)*Q(12)*Q(19),e2*A(5,L)*Q(12)**2*Q(16), &
          2*e2*A(5,L)*Q(12)**2*Q(13),2*A(5,L)*Q(12)**3*Q(13)/)
      CASE (35_I2B)
        QQ=Merge(A(2,L), A(1,L), G(1,L)==G(5,L))
        Q(1:15)=(/e3, e2+QQ-e1*QQ, e2+A(3,L)-e1*A(3,L), e2**2+&
          e4*e5*A(3,L), e2-e5*(QQ+A(3,L)), e2+A(5,L)-e1*A(5,L), e2-e5*(QQ+&
          A(5,L)), 2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(5,L)), e2-e5*(A(3,L)+&
          A(5,L)), e2-e5*(QQ+A(3,L)+A(5,L)), e2**2+e5**2*A(3,L)-e2*e5*(QQ+&
          2*A(3,L)+A(5,L)), 2*e2**2+e5**2*(QQ+A(3,L))-e2*e5*(1+QQ+2*A(3,L) &
          + A(5,L)), 2*e2**3-e5**3*A(3,L)+e2*e5**2*(QQ+4*A(3,L))-e2**2*e5*(1+&
          QQ+4*A(3,L)+A(5,L)), 2*e2**2-2*e2*e5*(QQ+A(3,L)+A(5,L))+&
          e5**2*(QQ*A(3,L)+(2*QQ+A(3,L))*A(5,L)), 2*e2**3-e5**3*A(3,L)*(QQ+&
          A(5,L))-2*e2**2*e5*(QQ+2*A(3,L)+A(5,L))+2*e2*e5**2*(QQ*A(5,L)+&
          A(3,L)*(1+QQ+A(5,L)))/)
        coe(L,:)=(/2*e2**3*Q(1), 2*e2**3*A(5,L)*Q(1), 2*e2**2*Q(1)*Q(3), e2**2*Q(8), &
          2*e2**3*A(5,L)*Q(1), 2*e2**3*A(5,L)*Q(1), e2**2*Q(12), 2*e2**2*Q(11), &
          2*e2*Q(1)*Q(4), 2*e2**3*Q(7), 2*e2**3*Q(10), 2*e2**3*A(5,L)*Q(1), e2*Q(13), &
          2*e2**3*A(5,L)*Q(1), 2*e2*Q(4)*Q(7), 2*e2**3*A(5,L)*Q(7), &
          2*e2*A(5,L)*Q(1)*Q(4), 2*e2**2*Q(5)*Q(9), e2**2*A(5,L)*Q(12), &
          2*e2**2*A(5,L)*Q(11), 2*e2*Q(1)*Q(3)**2, 2*e2**2*Q(2)*Q(6), &
          2*e2**2*A(5,L)*Q(1)*Q(3), e2**2*A(5,L)*Q(8), e2*Q(3)*Q(8), &
          2*e2**3*A(5,L)*Q(10), e2*A(5,L)*Q(13), e2*Q(3)*Q(12), &
          2*e2**2*A(5,L)*Q(1)*Q(3), e2**2*A(5,L)*Q(8), e2**2*(2*e2**2-2*e2*e5*(QQ+&
          A(3,L)+A(5,L))+e5**2*(A(3,L)*A(5,L)+QQ*(A(3,L)+2*A(5,L)))), &
          2*e2**2*Q(3)*Q(10), 2*e2*A(5,L)*Q(1)*Q(4), e2*Q(15), Q(4)*Q(8), &
          2*e2**3*A(5,L)*Q(7), 2*e2**2*Q(3)*Q(7), e2**2*A(5,L)*Q(12), &
          2*e2**2*A(5,L)*Q(1)*Q(3), e2**2*A(5,L)*Q(8), 2*e2*A(5,L)*Q(4)*Q(7), &
          2*e2*Q(3)**2*Q(7), 2*Q(2)*Q(4)*Q(6), 2*e2**2*A(5,L)*Q(2)*Q(6), &
          2*e2*A(5,L)*Q(1)*Q(3)**2, 2*e2**2*A(5,L)*Q(5)*Q(9), 2*e2**2*A(5,L)*Q(3)*Q(7), &
          A(5,L)*Q(4)*Q(8), e2*(2*e2**2-2*e2*e5*(QQ+A(3,L)+A(5,L))+&
          e5**2*(A(3,L)*A(5,L)+QQ*(A(3,L)+2*A(5,L))))*Q(3), e2*A(5,L)*Q(3)*Q(12), &
          e2**2*A(5,L)*(2*e2**2-2*e2*e5*(QQ+A(3,L)+A(5,L))+e5**2*(A(3,L)*A(5,L) &
          + QQ*(A(3,L)+2*A(5,L)))), e2*A(5,L)*Q(3)*Q(8), 2*e2**2*A(5,L)*Q(3)*Q(10), &
          2*e2*A(5,L)*Q(1)*Q(3)**2, 2*e2**2*A(5,L)*Q(2)*Q(6), e2*A(5,L)*Q(15), &
          Q(3)**2*Q(8), 2*e2*Q(2)*Q(3)*Q(6), e2*A(5,L)*Q(3)*Q(8), &
          2*e2*A(5,L)*Q(3)**2*Q(7), 2*A(5,L)*Q(2)*Q(4)*Q(6), 2*Q(2)*Q(3)**2*Q(6), &
          e2*A(5,L)*(2*e2**2-2*e2*e5*(QQ+A(3,L)+A(5,L))+e5**2*(A(3,L)*A(5,L)+&
          QQ*(A(3,L)+2*A(5,L))))*Q(3), A(5,L)*Q(3)**2*Q(8), &
          2*e2*A(5,L)*Q(2)*Q(3)*Q(6), 2*A(5,L)*Q(2)*Q(3)**2*Q(6)/)
      CASE (36_I2B)
      Q(1:11)=(/2*e2**2-2*e2*e5+e5**2, e3, 1-e1+2*e2, e2+A(1,L)-&
        e1*A(1,L), e2**2+e4*e5*A(1,L), e2+A(3,L)-e1*A(3,L), e2**2+&
        e4*e5*A(3,L), e2-e5*(A(1,L)+A(3,L)), e2**2-2*e2*e5*(A(1,L)+A(3,L))+&
        e5**2*(A(1,L)+A(3,L)), e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)), &
        e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+2*A(3,L))/)
      coe(L,:)=(/0._DP, e2**2*(e2**2+Q(2)**2), 0._DP, &
        0._DP, 2*e2**3*Q(2), e2**2*(e2**2+&
        Q(2)**2), 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, &
        e2**2*(e2**2+Q(2)**2), 0._DP, 2*e2**3*Q(2), 0._DP, &
        2*e2*Q(2)*Q(5), 2*e2**2*Q(7), 0._DP, e2**2*Q(3)*Q(8), 2*e2**2*Q(9), 0._DP, 0._DP, &
        e2**2*Q(3)*Q(6), e2*Q(1)*Q(4), 0._DP, e2*Q(3)*Q(10), 2*e2**2*Q(11), 0._DP, &
        e2**2*Q(3)*Q(6), e2*Q(1)*Q(4), 0._DP, 0._DP, &
        2*e2**2*Q(7), 0._DP, 0._DP, Q(1)*Q(5), 0._DP, &
        e2**2*Q(3)*Q(8), e2**2*Q(3)*Q(6), 2*e2**2*Q(2)*Q(4), 2*Q(5)*Q(7), 0._DP, 0._DP, &
        2*e2*Q(2)*Q(4)**2, 2*e2**2*Q(6)**2, 2*e2**2*Q(8)**2, Q(3)*Q(5)*Q(6), &
        2*e2*Q(4)*Q(7), 0._DP, 2*e2**2*Q(6)*Q(8), e2*Q(3)*Q(4)*Q(8), e2*Q(3)*Q(4)*Q(6), &
        2*e2*Q(6)*Q(10), 2*e2**2*Q(6)**2, Q(1)*Q(4)**2, 2*e2*Q(4)*Q(11), 0._DP, 0._DP, &
        e2*Q(3)*Q(4)*Q(6), 2*Q(5)*Q(6)**2, 2*Q(4)**2*Q(7), 0._DP, 2*e2*Q(4)*Q(6)*Q(8), &
        2*e2*Q(4)*Q(6)**2, Q(3)*Q(4)**2*Q(6), 2*Q(4)**2*Q(6)**2/)
      CASE (37_I2B)
        QQ=Merge(A(4,L), A(3,L), G(3,L)==G(5,L))
        Q(1:15)=(/e3, e2+QQ-e1*QQ, e2+A(1,L)-e1*A(1,L), e2**2+&
          e4*e5*A(1,L), e2-e5*(QQ+A(1,L)), e2+A(5,L)-e1*A(5,L), e2-e5*(QQ+&
          A(5,L)), 2*e2**2+e5**2*QQ-e2*e5*(1+QQ+A(5,L)), e2-e5*(A(1,L)+&
          A(5,L)), e2-e5*(QQ+A(1,L)+A(5,L)), e2**2+e5**2*A(1,L)-e2*e5*(QQ+&
          2*A(1,L)+A(5,L)), 2*e2**2+e5**2*(QQ+A(1,L))-e2*e5*(1+QQ+2*A(1,L) &
          + A(5,L)), 2*e2**3-e5**3*A(1,L)+e2*e5**2*(QQ+4*A(1,L))-e2**2*e5*(1+&
          QQ+4*A(1,L)+A(5,L)), 2*e2**2-2*e2*e5*(QQ+A(1,L)+A(5,L))+&
          e5**2*(QQ*A(1,L)+(2*QQ+A(1,L))*A(5,L)), 2*e2**3-e5**3*A(1,L)*(QQ+&
          A(5,L))-2*e2**2*e5*(QQ+2*A(1,L)+A(5,L))+2*e2*e5**2*(QQ*A(5,L)+&
          A(1,L)*(1+QQ+A(5,L)))/)
        coe(L,:)=(/2*e2**3*Q(1), 2*e2**3*A(5,L)*Q(1), e2**2*Q(8), 2*e2**2*Q(1)*Q(3), &
          2*e2**3*A(5,L)*Q(1), 2*e2**3*A(5,L)*Q(1), e2**2*Q(12), 2*e2**2*Q(11), &
          2*e2**3*Q(7), 2*e2*Q(1)*Q(4), e2*Q(13), 2*e2**3*A(5,L)*Q(1), 2*e2**3*Q(10), &
          2*e2**3*A(5,L)*Q(1), 2*e2*Q(4)*Q(7), 2*e2*A(5,L)*Q(1)*Q(4), &
          2*e2**3*A(5,L)*Q(7), 2*e2**2*Q(5)*Q(9), e2**2*A(5,L)*Q(12), &
          2*e2**2*A(5,L)*Q(11), 2*e2**2*Q(2)*Q(6), 2*e2*Q(1)*Q(3)**2, &
          e2**2*A(5,L)*Q(8), 2*e2**2*A(5,L)*Q(1)*Q(3), e2*Q(3)*Q(8), e2*A(5,L)*Q(13), &
          2*e2**3*A(5,L)*Q(10), e2**2*(2*e2**2-2*e2*e5*(QQ+A(1,L)+A(5,L))+&
          e5**2*(2*QQ*A(5,L)+A(1,L)*(QQ+A(5,L)))), e2**2*A(5,L)*Q(8), &
          2*e2**2*A(5,L)*Q(1)*Q(3), e2*Q(3)*Q(12), e2*Q(15), 2*e2**3*A(5,L)*Q(7), &
          2*e2**2*Q(3)*Q(10), 2*e2**2*Q(3)*Q(7), 2*e2*A(5,L)*Q(1)*Q(4), Q(4)*Q(8), &
          e2**2*A(5,L)*Q(12), e2**2*A(5,L)*Q(8), 2*e2**2*A(5,L)*Q(1)*Q(3), &
          2*e2*A(5,L)*Q(4)*Q(7), 2*Q(2)*Q(4)*Q(6), 2*e2*Q(3)**2*Q(7), &
          2*e2*A(5,L)*Q(1)*Q(3)**2, 2*e2**2*A(5,L)*Q(2)*Q(6), 2*e2**2*A(5,L)*Q(5)*Q(9), &
          A(5,L)*Q(4)*Q(8), 2*e2**2*A(5,L)*Q(3)*Q(7), e2*(2*e2**2-2*e2*e5*(QQ+&
          A(1,L)+A(5,L))+e5**2*(2*QQ*A(5,L)+A(1,L)*(QQ+A(5,L))))*Q(3), &
          e2**2*A(5,L)*(2*e2**2-2*e2*e5*(QQ+A(1,L)+A(5,L))+e5**2*(2*QQ*A(5,L)+&
          A(1,L)*(QQ+A(5,L)))), e2*A(5,L)*Q(3)*Q(12), e2*A(5,L)*Q(3)*Q(8), &
          e2*A(5,L)*Q(15), 2*e2**2*A(5,L)*Q(2)*Q(6), 2*e2*A(5,L)*Q(1)*Q(3)**2, &
          2*e2**2*A(5,L)*Q(3)*Q(10), 2*e2*Q(2)*Q(3)*Q(6), Q(3)**2*Q(8), &
          e2*A(5,L)*Q(3)*Q(8), 2*A(5,L)*Q(2)*Q(4)*Q(6), 2*e2*A(5,L)*Q(3)**2*Q(7), &
          2*Q(2)*Q(3)**2*Q(6), e2*A(5,L)*(2*e2**2-2*e2*e5*(QQ+A(1,L)+A(5,L))+&
          e5**2*(2*QQ*A(5,L)+A(1,L)*(QQ+A(5,L))))*Q(3), 2*e2*A(5,L)*Q(2)*Q(3)*Q(6), &
          A(5,L)*Q(3)**2*Q(8), 2*A(5,L)*Q(2)*Q(3)**2*Q(6)/)
      CASE (38_I2B)
        X(1:2)=Merge(A(1:2,L), A(2:1:-1,L), ANY(G(1,L)==G(5:6,L)))
        Q(1:20)=(/2*e2**2-2*e2*e5+e5**2, e3, 1-e1+2*e2, e2+X(1)-e1*X(1), &
          e2+X(2)-e1*X(2), e2-e5*(X(1)+X(2)), 2*e2-e5*(X(1)+X(2)), 2*e2**2 &
          + e5**2*X(2)-2*e2*e5*(X(1)+X(2)), e2**2-2*e2*e5*(X(1)+X(2))+&
          e5**2*(X(1)+X(2)), 2*e2**2+e5**2*X(2)-e2*e5*(1+X(1)+X(2)), 4*e2**2 &
          + e5**2*X(2)-e2*e5*(1+2*X(1)+2*X(2)), 4*e2**2+e5**2*(X(1)+2*X(2))-&
          e2*e5*(3+2*X(1)+2*X(2)), 4*e2**3-e5**3*X(2)+e2*e5**2*(1+X(1)+&
          2*X(2))-e2**2*e5*(3+2*X(1)+2*X(2)), 4*e2**3-e5**3*X(2)+&
          e2*e5**2*(3*X(1)+4*X(2))-e2**2*e5*(1+6*X(1)+6*X(2)), 8*e2**4-&
          4*e2*e5**3*X(2)+e5**4*X(2)-4*e2**3*e5*(1+2*X(1)+2*X(2))+&
          e2**2*e5**2*(1+4*X(1)+8*X(2)), 2*e2**3-2*e5**3*X(1)*X(2)-&
          4*e2**2*e5*(X(1)+X(2))+e2*e5**2*(X(1)+X(1)**2+X(2)+4*X(1)*X(2)+&
          X(2)**2), 2*e2**4-4*e2*e5**3*X(1)*X(2)+e5**4*X(1)*X(2)-4*e2**3*e5*(X(1) &
          + X(2))+e2**2*e5**2*(X(1)+X(1)**2+X(2)+6*X(1)*X(2)+X(2)**2), &
          4*e2**4-8*e2**3*e5*(X(1)+X(2))+e5**4*X(1)*X(2)*(X(1)+X(2))-&
          4*e2*e5**3*X(1)*X(2)*(1+X(1)+X(2))+e2**2*e5**2*(X(1)+3*X(1)**2+X(2) &
          + 14*X(1)*X(2)+3*X(2)**2), 8*e2**4+e5**4*X(2)**2-2*e2*e5**3*X(2)*(1+&
          X(1)+X(2))-4*e2**3*e5*(1+2*X(1)+2*X(2))+e2**2*e5**2*(1+2*X(1)**2 &
          + 2*X(2)*(3+X(2))+X(1)*(2+4*X(2))), 8*e2**4+e5**4*X(1)*X(2)-&
          e2*e5**3*X(2)*(1+7*X(1)+X(2))-2*e2**3*e5*(1+6*X(1)+6*X(2))+&
          2*e2**2*e5**2*(X(1)**2+X(2)*(3+X(2))+X(1)*(2+6*X(2)))/)
        coe(L,:)=(/0._DP, 4*e2**2*(e2**2+Q(2)**2), 0._DP, &
          0._DP, 2*e2**2*e4**2, 2*e2**2*e4**2, &
          0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 8*e2**3*Q(2), &
          0._DP, 2*e2**2*e4**2, 0._DP, 4*e2**2*Q(3)*Q(6), &
          4*e2**2*Q(3)*Q(6), 0._DP, 2*e2*Q(3)*Q(8), &
          8*e2**2*Q(9), 0._DP, 0._DP, 2*e2*Q(13), &
          2*e2*Q(13), 0._DP, 2*e2*Q(14), 2*e2*Q(14), &
          0._DP, 2*e2**2*Q(12), 2*e2**2*Q(12), 0._DP, 0._DP, &
          4*e2**2*Q(3)*Q(6), 0._DP, 0._DP, 4*e2**2*Q(3)*Q(6), 0._DP, Q(15), e2*Q(3)*Q(11), &
          e2*Q(3)*Q(11), 8*e2**2*Q(6)**2, 0._DP, 0._DP, 4*e2*Q(3)*Q(4)*Q(5), &
          4*e2*Q(3)*Q(4)*Q(5), 4*Q(17), 2*e2*Q(6)*Q(11), 2*e2*Q(6)*Q(11), 0._DP, Q(20), &
          Q(20), 2*e2*Q(7)*Q(10), 4*e2*Q(16), 4*e2*Q(3)*Q(4)*Q(5), 4*e2*Q(3)*Q(4)*Q(5), &
          4*e2*Q(16), 0._DP, 0._DP, Q(19), 8*e2*Q(4)*Q(5)*Q(6), 8*e2*Q(4)*Q(5)*Q(6), 0._DP, &
          2*Q(18), 2*Q(4)*Q(5)*Q(11), 2*Q(4)*Q(5)*Q(11), 8*Q(4)**2*Q(5)**2/)
      CASE (39_I2B)
        IF(G(1,L)==G(3,L))THEN
          x(1:3)=(/A(2,L),A(1,L),A(4,L)/)
        ELSE IF(G(2,L)==G(3,L))THEN
          x(1:3)=(/A(1,L),A(2,L),A(4,L)/)
        ELSE IF(G(1,L)==G(4,L))THEN
          x(1:3)=(/A(2,L),A(1,L),A(3,L)/)
        ELSE
          x(1:3)=(/A(1,L),A(2,L),A(3,L)/)
        END IF
        Q(1:25)=(/2*e2**2-2*e2*e5+e5**2, e3, 1-e1+2*e2, e2+X(1)-e1*X(1), &
          e2+X(2)-e1*X(2), e2-e5*(X(1)+X(2)), 2*e2-e5*(X(1)+X(2)), 4*e2**3 &
          - e5**3*X(1)-2*e2**2*e5*(2+X(1)+X(2))+e2*e5**2*(1+2*X(1)+2*X(2)), &
          4*e2**3-e5**3*X(2)-2*e2**2*e5*(2+X(1)+X(2))+e2*e5**2*(1+2*X(1)+&
          2*X(2)), e2+X(3)-e1*X(3), e2-e5*(X(2)+X(3)), e2**2+e5**2*X(2)-&
          e2*e5*(X(1)+2*X(2)+X(3)), 4*e2**2+e5**2*X(3)-e2*e5*(1+2*X(2)+&
          2*X(3)), 2*e2**2+e5**2*X(2)-e2*e5*(X(1)+3*X(2)+2*X(3)), 4*e2**2+&
          e5**2*(X(1)+3*X(2)+X(3))-e2*e5*(1+4*X(1)+6*X(2)+2*X(3)), 8*e2**3 &
          - e5**3*(X(1)+X(2)+X(3))-2*e2**2*e5*(3+2*X(1)+4*X(2)+2*X(3))+&
          e2*e5**2*(1+3*X(1)+7*X(2)+4*X(3)), 8*e2**3-e5**3*(2*X(2)+X(3))-&
          2*e2**2*e5*(3+2*X(1)+4*X(2)+2*X(3))+e2*e5**2*(1+3*X(1)+7*X(2)+&
          4*X(3)), 2*e2**3-e5**3*X(1)*X(2)-2*e2**2*e5*(X(1)+2*X(2)+X(3))+&
          e2*e5**2*(X(2)*(1+3*X(1)+X(2))+(X(1)+X(2))*X(3)), 2*e2**3-&
          e5**3*X(2)*X(3)-2*e2**2*e5*(X(1)+2*X(2)+X(3))+e2*e5**2*(X(2)*(1+&
          X(1)+X(2))+(X(1)+3*X(2))*X(3)), 2*e2**3-e5**3*X(2)*(X(1)+X(3))-&
          2*e2**2*e5*(X(1)+2*X(2)+X(3))+e2*e5**2*(X(2)*(1+3*X(1)+X(2))+&
          (X(1)+3*X(2))*X(3)), 4*e2**3-e5**3*X(2)*X(3)-2*e2**2*e5*(X(1)+3*X(2) &
          + 2*X(3))+e2*e5**2*(X(2)*(1+X(1)+X(2))+(X(1)+5*X(2))*X(3)), 8*e2**3 &
          - 2*e2**2*e5*(1+4*X(1)+6*X(2)+2*X(3))-e5**3*(4*X(1)*X(2)+(X(1)+&
          X(2))*X(3))+2*e2*e5**2*(X(3)+X(2)*(2+X(2)+X(3))+X(1)*(1+5*X(2)+&
          X(3))), 8*e2**4+e5**4*X(1)*X(3)-2*e2**3*e5*(3+2*X(1)+4*X(2)+2*X(3)) &
          - e2*e5**3*(X(1)+X(1)*X(2)+X(2)**2+X(3)+2*(X(1)+X(2))*X(3))+&
          e2**2*e5**2*(1+4*X(3)+X(1)*(3+2*X(2)+2*X(3))+X(2)*(5+2*X(2)+&
          2*X(3))), 8*e2**4+e5**4*X(2)*X(3)-2*e2**3*e5*(3+2*X(1)+4*X(2)+&
          2*X(3))-e2*e5**3*(X(2)*(1+X(1)+X(2))+X(3)+2*(X(1)+X(2))*X(3))+&
          e2**2*e5**2*(1+4*X(3)+X(1)*(3+2*X(2)+2*X(3))+X(2)*(5+2*X(2)+&
          2*X(3))), 4*e2**4+e5**4*X(1)*X(2)*X(3)-4*e2**3*e5*(X(1)+2*X(2)+X(3)) &
          - e2*e5**3*X(2)*(X(3)+2*X(2)*X(3)+X(1)*(1+2*X(2)+4*X(3)))+&
          e2**2*e5**2*(3*X(2)**2+3*X(1)*X(3)+X(2)*(1+7*X(1)+7*X(3)))/)
        coe(L,:)=(/0._DP, 4*e2**2*Q(2)*Q(3), 0._DP, 0._DP, &
          4*e2**2*Q(2)*Q(3), 2*e2*(e2**2+&
          Q(2)**2)*Q(3), 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, &
          4*e2**2*Q(2)*Q(3), 0._DP, 2*e2*(e2**2+&
          Q(2)**2)*Q(3), 0._DP, 4*e2**2*Q(3)*Q(6), 4*e2*Q(1)*Q(11), 0._DP, e2*Q(17), &
          8*e2**2*Q(12), 0._DP, 0._DP, 2*e2*Q(2)*Q(13), 2*e2*Q(8), 0._DP, 2*e2**2*Q(15), &
          2*e2*Q(3)*Q(14), 0._DP, 2*e2*Q(2)*Q(13), &
          2*e2*Q(9), 0._DP, 0._DP, 8*e2**2*Q(2)*Q(11), 0._DP, &
          0._DP, 4*e2**2*Q(3)*Q(6), 0._DP, e2*Q(16), Q(1)*Q(13), e2*e4**2*Q(7), &
          8*e2**2*Q(6)*Q(11), 0._DP, 0._DP, 4*e2*Q(3)*Q(4)*Q(5), 4*Q(1)*Q(5)*Q(10), 4*e2*Q(20), &
          2*e2*Q(6)*Q(13), 2*e2*Q(3)*Q(7)*Q(11), 0._DP, Q(3)*Q(21), e2*(8*e2**3-&
          2*e2**2*e5*(1+4*X(1)+6*X(2)+2*X(3))-e5**3*(4*X(1)*X(2)+(X(1)+&
          X(2))*X(3))+2*e2*e5**2*(X(1)+X(2)*(2+5*X(1)+X(2))+X(3)+(X(1)+&
          X(2))*X(3))), 8*e2**4+e5**4*X(2)*X(3)-2*e2**3*e5*(3+2*X(1)+4*X(2)+&
          2*X(3))-e2*e5**3*(X(2)*(1+X(1)+X(2))+X(3)+2*(X(1)+X(2))*X(3))+&
          e2**2*e5**2*(1+3*X(1)+5*X(2)+2*X(1)*X(2)+2*X(2)**2+2*(2+X(1)+&
          X(2))*X(3)), 4*e2*Q(19), 8*e2*Q(2)*Q(5)*Q(10), 4*e2*Q(3)*Q(4)*Q(5), &
          4*e2*Q(18), 0._DP, 0._DP, 8*e2**4+e5**4*X(1)*X(3)-2*e2**3*e5*(3+2*X(1)+4*X(2) &
          + 2*X(3))+e2**2*e5**2*(1+3*X(1)+5*X(2)+2*X(1)*X(2)+2*X(2)**2+2*(2 &
          + X(1)+X(2))*X(3))-e2*e5**3*(X(1)+X(3)+(X(1)+X(2))*(X(2)+&
          2*X(3))), 8*e2*Q(5)*Q(6)*Q(10), 8*e2*Q(4)*Q(5)*Q(11), 0._DP, 2*Q(25), &
          2*Q(3)*Q(5)*Q(7)*Q(10), 2*Q(4)*Q(5)*Q(13), 8*Q(4)*Q(5)**2*Q(10)/)
      CASE (40_I2B)
        IF(G(5,L)==G(1,L))THEN
          x(1:3)=(/A(2,L),A(5,L),A(6,L)/)
        ELSE IF(G(5,L)==G(2,L))THEN
          x(1:3)=(/A(1,L),A(5,L),A(6,L)/)
        ELSE IF(G(6,L)==G(1,L))THEN
          x(1:3)=(/A(2,L),A(6,L),A(5,L)/)
        ELSE
          x(1:3)=(/A(1,L),A(6,L),A(5,L)/)
        END IF
        Q(1:25)=(/2*e2**2 - 2*e2*e5 + e5**2, e3, 1 - e1 + 2*e2, e2 + X(1) - e1*X(1), &
          e2 + X(2) - e1*X(2), e2 - e5*(X(1) + X(2)), 4*e2**2 + e5**2*X(1) - e2*e5*(1 + &
          2*X(1) + 2*X(2)), e2 + X(3) - e1*X(3), e2 - e5*(X(2) + X(3)), 2*e2 - e5*(X(2) &
          + X(3)), e2**2 + e5**2*X(2) - e2*e5*(X(1) + 2*X(2) + X(3)), 2*e2**2 + &
          e5**2*X(2) - e2*e5*(2*X(1) + 3*X(2) + X(3)), 4*e2**3 - e5**3*X(2) - &
          2*e2**2*e5*(2 + X(2) + X(3)) + e2*e5**2*(1 + 2*X(2) + 2*X(3)), 4*e2**3 - &
          e5**3*X(3) - 2*e2**2*e5*(2 + X(2) + X(3)) + e2*e5**2*(1 + 2*X(2) + 2*X(3)), &
          8*e2**3 - e5**3*(X(1) + 2*X(2)) - 2*e2**2*e5*(3 + 2*X(1) + 4*X(2) + 2*X(3)) + &
          e2*e5**2*(1 + 4*X(1) + 7*X(2) + 3*X(3)), 8*e2**3 - e5**3*(X(1) + X(2) + X(3)) &
          - 2*e2**2*e5*(3 + 2*X(1) + 4*X(2) + 2*X(3)) + e2*e5**2*(1 + 4*X(1) + 7*X(2) + &
          3*X(3)), 4*e2**2 + e5**2*(X(1) + 3*X(2) + X(3)) - e2*e5*(1 + 2*X(1) + 6*X(2) &
          + 4*X(3)), 2*e2**3 - e5**3*X(1)*X(2) - 2*e2**2*e5*(X(1) + 2*X(2) + X(3)) + &
          e2*e5**2*(X(2)*(1 + 3*X(1) + X(2)) + (X(1) + X(2))*X(3)), 4*e2**3 - &
          e5**3*X(1)*X(2) - 2*e2**2*e5*(2*X(1) + 3*X(2) + X(3)) + e2*e5**2*(X(2)*(1 + &
          5*X(1) + X(2)) + (X(1) + X(2))*X(3)), 2*e2**3 - e5**3*X(2)*X(3) - &
          2*e2**2*e5*(X(1) + 2*X(2) + X(3)) + e2*e5**2*(X(2)*(1 + X(1) + X(2)) + (X(1) &
          + 3*X(2))*X(3)), 2*e2**3 - e5**3*X(2)*(X(1) + X(3)) - 2*e2**2*e5*(X(1) + &
          2*X(2) + X(3)) + e2*e5**2*(X(2)*(1 + 3*X(1) + X(2)) + (X(1) + 3*X(2))*X(3)), &
          8*e2**4 + e5**4*X(1)*X(3) - 2*e2**3*e5*(3 + 2*X(1) + 4*X(2) + 2*X(3)) - &
          e2*e5**3*(X(3) + X(2)*(X(2) + X(3)) + X(1)*(1 + 2*X(2) + 2*X(3))) + &
          e2**2*e5**2*(1 + 3*X(3) + 2*X(1)*(2 + X(2) + X(3)) + X(2)*(5 + 2*X(2) + &
          2*X(3))), 8*e2**4 + e5**4*X(1)*X(2) - 2*e2**3*e5*(3 + 2*X(1) + 4*X(2) + &
          2*X(3)) - e2*e5**3*(X(2)*(1 + X(2) + X(3)) + X(1)*(1 + 2*X(2) + 2*X(3))) + &
          e2**2*e5**2*(1 + 3*X(3) + 2*X(1)*(2 + X(2) + X(3)) + X(2)*(5 + 2*X(2) + &
          2*X(3))), 8*e2**3 - 2*e2**2*e5*(1 + 2*X(1) + 6*X(2) + 4*X(3)) - &
          e5**3*(4*X(2)*X(3) + X(1)*(X(2) + X(3))) + 2*e2*e5**2*(X(3) + X(1)*(1 + X(2) &
          + X(3)) + X(2)*(2 + X(2) + 5*X(3))), 4*e2**4 + e5**4*X(1)*X(2)*X(3) - &
          4*e2**3*e5*(X(1) + 2*X(2) + X(3)) - e2*e5**3*X(2)*(X(3) + 2*X(2)*X(3) + &
          X(1)*(1 + 2*X(2) + 4*X(3))) + e2**2*e5**2*(3*X(2)**2 + 3*X(1)*X(3) + X(2)*(1 &
          + 7*X(1) + 7*X(3)))/)
        coe(L,:)=(/0._DP, 4*e2**2*Q(2)*Q(3), 0._DP, 0._DP, 2*e2*(e2**2 + Q(2)**2)*Q(3), &
          4*e2**2*Q(2)*Q(3), 0._DP, 0._DP, 0._DP, 0._DP, &
          0._DP, 4*e2**2*Q(2)*Q(3), 0._DP, 2*e2*(e2**2 + &
          Q(2)**2)*Q(3), 0._DP, 4*e2*Q(1)*Q(6), 4*e2**2*Q(3)*Q(9), 0._DP, e2*Q(15), &
          8*e2**2*Q(11), 0._DP, 0._DP, 2*e2*Q(14), 2*e2*Q(2)*Q(7), 0._DP, 2*e2*Q(3)*Q(12), &
          2*e2**2*Q(17), 0._DP, 2*e2*Q(13), 2*e2*Q(2)*Q(7), &
          0._DP, 0._DP, 4*e2**2*Q(3)*Q(9), 0._DP, 0._DP, &
          8*e2**2*Q(2)*Q(6), 0._DP, e2*Q(16), e2*e4**2*Q(10), Q(1)*Q(7), 8*e2**2*Q(6)*Q(9), &
          0._DP, 0._DP, 4*Q(1)*Q(4)*Q(5), 4*e2*Q(3)*Q(5)*Q(8), 4*e2*Q(21), &
          2*e2*Q(3)*Q(6)*Q(10), 2*e2*Q(7)*Q(9), 0._DP, e2*Q(24), Q(3)*Q(19), Q(23), &
          4*e2*Q(20), 4*e2*Q(3)*Q(5)*Q(8), 8*e2*Q(2)*Q(4)*Q(5), 4*e2*Q(18), 0._DP, 0._DP, &
          Q(22), 8*e2*Q(5)*Q(6)*Q(8), 8*e2*Q(4)*Q(5)*Q(9), 0._DP, 2*Q(25), &
          2*Q(5)*Q(7)*Q(8), 2*Q(3)*Q(4)*Q(5)*Q(10), 8*Q(4)*Q(5)**2*Q(8)/)
      CASE (41_I2B)
        Q(1:8)=(/e2+A(1,L)-e1*A(1,L), e2**2+e4*e5*A(1,L), e2+A(3,L)-&
          e1*A(3,L), e2**2+e4*e5*A(3,L), e2-e5*(A(1,L)+A(3,L)), e2**2-&
          2*e2*e5*(A(1,L)+A(3,L))+e5**2*(A(1,L)+A(3,L)), e2**2+e5**2*A(1,L)-&
          e2*e5*(2*A(1,L)+A(3,L)), e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+&
          2*A(3,L))/)
        coe(L,:)=(/0._DP, 2*e2**4, 0._DP, 0._DP, 2*e2**4, 2*e2**4, &
          0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 2*e2**4, 0._DP, &
          2*e2**4, 0._DP, 2*e2**2*Q(2), 2*e2**2*Q(4), 0._DP, &
          2*e2**3*Q(5), 2*e2**2*Q(6), 0._DP, 0._DP, &
          2*e2**3*Q(3), 2*e2**3*Q(1), 0._DP, 2*e2**2*Q(7), &
          2*e2**2*Q(8), 0._DP, 2*e2**3*Q(3), &
          2*e2**3*Q(1), 0._DP, 0._DP, 2*e2**2*Q(4), 0._DP, &
          0._DP, 2*e2**2*Q(2), 0._DP, 2*e2**3*Q(5), &
          2*e2**3*Q(3), 2*e2**3*Q(1), 2*Q(2)*Q(4), &
          0._DP, 0._DP, 2*e2**2*Q(1)**2, &
          2*e2**2*Q(3)**2, 2*e2**2*Q(5)**2, 2*e2*Q(2)*Q(3), 2*e2*Q(1)*Q(4), 0._DP, &
          2*e2**2*Q(3)*Q(5), 2*e2**2*Q(1)*Q(5), 2*e2**2*Q(1)*Q(3), 2*e2*Q(3)*Q(7), &
          2*e2**2*Q(3)**2, 2*e2**2*Q(1)**2, 2*e2*Q(1)*Q(8), &
          0._DP, 0._DP, 2*e2**2*Q(1)*Q(3), &
          2*Q(2)*Q(3)**2, 2*Q(1)**2*Q(4), 0._DP, &
          2*e2*Q(1)*Q(3)*Q(5), 2*e2*Q(1)*Q(3)**2, &
          2*e2*Q(1)**2*Q(3), 2*Q(1)**2*Q(3)**2/)
      CASE (42_I2B)
        Q(1:14)=(/e2+A(1,L)-e1*A(1,L), e2**2+e4*e5*A(1,L), e2+A(3,L)-&
        e1*A(3,L), e2-e5*(A(1,L)+A(3,L)), e2+A(4,L)-e1*A(4,L), e2-&
        e5*(A(1,L)+A(4,L)), e2-e5*(A(3,L)+A(4,L)), 2*e2-e5*(A(3,L)+A(4,L)), &
        e2-e5*(A(1,L)+A(3,L)+A(4,L)), 2*e2-e5*(2*A(1,L)+A(3,L)+A(4,L)), &
        e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)+A(4,L)), 2*e2**2+&
        2*e5**2*A(1,L)-e2*e5*(4*A(1,L)+A(3,L)+A(4,L)), 2*e2**2-&
        2*e2*e5*(A(1,L)+A(3,L)+A(4,L))+e5**2*(2*A(3,L)*A(4,L)+A(1,L)*(A(3,L) &
        + A(4,L))), 2*e2**3-e5**3*A(1,L)*(A(3,L)+A(4,L))-2*e2**2*e5*(2*A(1,L)+&
        A(3,L)+A(4,L))+2*e2*e5**2*(A(3,L)*A(4,L)+A(1,L)*(1+A(3,L)+&
        A(4,L)))/)
        coe(L,:)=(/2*e2**4, 2*e2**4*A(5,L), e2**3*Q(8), 2*e2**3*Q(1), &
        2*e2**4*A(5,L), 2*e2**4*A(5,L), e2**3*Q(10), 2*e2**2*Q(11), 2*e2**3*Q(7), &
        2*e2**2*Q(2), e2**2*Q(12), 2*e2**4*A(5,L), 2*e2**3*Q(9), 2*e2**4*A(5,L), &
        2*e2*Q(2)*Q(7), 2*e2**2*A(5,L)*Q(2), 2*e2**3*A(5,L)*Q(7), 2*e2**2*Q(4)*Q(6), &
        e2**3*A(5,L)*Q(10), 2*e2**2*A(5,L)*Q(11), 2*e2**2*Q(3)*Q(5), 2*e2**2*Q(1)**2, &
        e2**3*A(5,L)*Q(8), 2*e2**3*A(5,L)*Q(1), e2**2*Q(1)*Q(8), e2**2*A(5,L)*Q(12), &
        2*e2**3*A(5,L)*Q(9), e2**2*Q(13), e2**3*A(5,L)*Q(8), 2*e2**3*A(5,L)*Q(1), &
        e2**2*Q(1)*Q(10), e2*Q(14), 2*e2**3*A(5,L)*Q(7), 2*e2**2*Q(1)*Q(9), &
        2*e2**2*Q(1)*Q(7), 2*e2**2*A(5,L)*Q(2), e2*Q(2)*Q(8), e2**3*A(5,L)*Q(10), &
        e2**3*A(5,L)*Q(8), 2*e2**3*A(5,L)*Q(1), 2*e2*A(5,L)*Q(2)*Q(7), &
        2*Q(2)*Q(3)*Q(5), 2*e2*Q(1)**2*Q(7), 2*e2**2*A(5,L)*Q(1)**2, &
        2*e2**2*A(5,L)*Q(3)*Q(5), 2*e2**2*A(5,L)*Q(4)*Q(6), e2*A(5,L)*Q(2)*Q(8), &
        2*e2**2*A(5,L)*Q(1)*Q(7), e2*Q(1)*Q(13), e2**2*A(5,L)*Q(13), &
        e2**2*A(5,L)*Q(1)*Q(10), e2**2*A(5,L)*Q(1)*Q(8), e2*A(5,L)*Q(14), &
        2*e2**2*A(5,L)*Q(3)*Q(5), 2*e2**2*A(5,L)*Q(1)**2, 2*e2**2*A(5,L)*Q(1)*Q(9), &
        2*e2*Q(1)*Q(3)*Q(5), e2*Q(1)**2*Q(8), e2**2*A(5,L)*Q(1)*Q(8), &
        2*A(5,L)*Q(2)*Q(3)*Q(5), 2*e2*A(5,L)*Q(1)**2*Q(7), 2*Q(1)**2*Q(3)*Q(5), &
        e2*A(5,L)*Q(1)*Q(13), 2*e2*A(5,L)*Q(1)*Q(3)*Q(5), e2*A(5,L)*Q(1)**2*Q(8), &
        2*A(5,L)*Q(1)**2*Q(3)*Q(5)/)
      CASE (43_I2B)
        Q(1:14)=(/e2+A(1,L)-e1*A(1,L), e2+A(2,L)-e1*A(2,L), e2-e5*(A(1,L) &
        + A(2,L)), 2*e2-e5*(A(1,L)+A(2,L)), e2+A(3,L)-e1*A(3,L), e2**2+&
        e4*e5*A(3,L), 2*e2*(e2+A(1,L)-e1*A(1,L))*(e2+A(2,L)-e1*A(2,L))-(1-&
        e1+2*e2)*e5*(2*e2-e5*(A(1,L)+A(2,L)))*A(3,L), e2-e5*(A(1,L)+&
        A(3,L)), e2-e5*(A(2,L)+A(3,L)), e2-e5*(A(1,L)+A(2,L)+A(3,L)), 2*e2 &
        - e5*(A(1,L)+A(2,L)+2*A(3,L)), e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+&
        A(2,L)+2*A(3,L)), 2*e2**2+2*e5**2*A(3,L)-e2*e5*(A(1,L)+A(2,L)+&
        4*A(3,L)), 2*e2**2-2*e2*e5*(A(1,L)+A(2,L)+A(3,L))+&
        e5**2*(2*A(1,L)*A(2,L)+(A(1,L)+A(2,L))*A(3,L))/)
        coe(L,:)=(/2*e2**4, 2*e2**4*A(5,L), 2*e2**3*Q(5), e2**3*Q(4), &
        2*e2**4*A(5,L), 2*e2**4*A(5,L), e2**3*Q(11), 2*e2**2*Q(12), 2*e2**2*Q(6), &
        2*e2**3*Q(3), 2*e2**3*Q(10), 2*e2**4*A(5,L), e2**2*Q(13), 2*e2**4*A(5,L), &
        2*e2*Q(3)*Q(6), 2*e2**3*A(5,L)*Q(3), 2*e2**2*A(5,L)*Q(6), 2*e2**2*Q(8)*Q(9), &
        e2**3*A(5,L)*Q(11), 2*e2**2*A(5,L)*Q(12), 2*e2**2*Q(5)**2, 2*e2**2*Q(1)*Q(2), &
        2*e2**3*A(5,L)*Q(5), e2**3*A(5,L)*Q(4), e2**2*Q(4)*Q(5), &
        2*e2**3*A(5,L)*Q(10), e2**2*A(5,L)*Q(13), e2**2*Q(5)*Q(11), &
        2*e2**3*A(5,L)*Q(5), e2**3*A(5,L)*Q(4), e2**2*Q(14), 2*e2**2*Q(5)*Q(10), &
        2*e2**2*A(5,L)*Q(6), e2*Q(7), e2*Q(4)*Q(6), 2*e2**3*A(5,L)*Q(3), &
        2*e2**2*Q(3)*Q(5), e2**3*A(5,L)*Q(11), 2*e2**3*A(5,L)*Q(5), &
        e2**3*A(5,L)*Q(4), 2*e2*A(5,L)*Q(3)*Q(6), 2*e2*Q(3)*Q(5)**2, &
        2*Q(1)*Q(2)*Q(6), 2*e2**2*A(5,L)*Q(1)*Q(2), 2*e2**2*A(5,L)*Q(5)**2, &
        2*e2**2*A(5,L)*Q(8)*Q(9), 2*e2**2*A(5,L)*Q(3)*Q(5), e2*A(5,L)*Q(4)*Q(6), &
        e2*Q(5)*Q(14), e2**2*A(5,L)*Q(5)*Q(11), e2**2*A(5,L)*Q(14), &
        e2**2*A(5,L)*Q(4)*Q(5), 2*e2**2*A(5,L)*Q(5)*Q(10), 2*e2**2*A(5,L)*Q(5)**2, &
        2*e2**2*A(5,L)*Q(1)*Q(2), e2*A(5,L)*Q(7), e2*Q(4)*Q(5)**2, &
        2*e2*Q(1)*Q(2)*Q(5), e2**2*A(5,L)*Q(4)*Q(5), 2*e2*A(5,L)*Q(3)*Q(5)**2, &
        2*A(5,L)*Q(1)*Q(2)*Q(6), 2*Q(1)*Q(2)*Q(5)**2, e2*A(5,L)*Q(5)*Q(14), &
        e2*A(5,L)*Q(4)*Q(5)**2, 2*e2*A(5,L)*Q(1)*Q(2)*Q(5), &
        2*A(5,L)*Q(1)*Q(2)*Q(5)**2/)
      CASE (44_I2B)
        Q(1:21)=(/2*e2**2-2*e2*e5+e5**2, e3, 1-e1+2*e2, e2+A(1,L)-&
        e1*A(1,L), e2+A(2,L)-e1*A(2,L), e2-e5*(A(1,L)+A(2,L)), 2*e2-&
        e5*(A(1,L)+A(2,L)), e2+A(3,L)-e1*A(3,L), e2+A(4,L)-e1*A(4,L), e2-&
        e5*(A(3,L)+A(4,L)), 2*e2-e5*(A(3,L)+A(4,L)), e2-e5*(A(1,L)+A(2,L)+&
        A(3,L)+A(4,L)), 2*e2-e5*(A(1,L)+A(2,L)+A(3,L)+A(4,L)), 2*e2-&
        e5*(2*A(1,L)+2*A(2,L)+A(3,L)+A(4,L)), 2*e2**2-2*e2*e5*(A(1,L)+&
        A(2,L)+A(3,L)+A(4,L))+e5**2*((A(1,L)+A(2,L))*A(3,L)+(A(1,L)+&
        A(2,L)+2*A(3,L))*A(4,L)), 2*e2-e5*(A(1,L)+A(2,L)+2*(A(3,L)+&
        A(4,L))), 4*e2**2+e5**2*((A(1,L)+A(2,L))*A(3,L)+(A(1,L)+A(2,L)+&
        4*A(3,L))*A(4,L))-2*e2*e5*(A(1,L)+A(2,L)+2*(A(3,L)+A(4,L))), 2*e2**2 &
        - 2*e2*e5*(A(1,L)+A(2,L)+A(3,L)+A(4,L))+e5**2*(A(2,L)*(A(3,L)+&
        A(4,L))+A(1,L)*(2*A(2,L)+A(3,L)+A(4,L))), 2*e2**2-2*e2*e5*(A(1,L)+&
        A(2,L)+A(3,L)+A(4,L))+e5**2*(2*A(3,L)*A(4,L)+A(2,L)*(A(3,L)+A(4,L)) &
        + A(1,L)*(2*A(2,L)+A(3,L)+A(4,L))), 4*e2**2-2*e2*e5*(2*A(1,L)+&
        2*A(2,L)+A(3,L)+A(4,L))+e5**2*(A(2,L)*(A(3,L)+A(4,L))+&
        A(1,L)*(4*A(2,L)+A(3,L)+A(4,L))), 4*e2**3-4*e2**2*e5*(A(1,L)+A(2,L)+&
        A(3,L)+A(4,L))+e2*e5**2*(4*A(3,L)*A(4,L)+3*A(2,L)*(A(3,L)+A(4,L))+&
        A(1,L)*(4*A(2,L)+3*(A(3,L)+A(4,L))))-2*e5**3*(A(2,L)*A(3,L)*A(4,L)+&
        A(1,L)*(A(3,L)*A(4,L)+A(2,L)*(A(3,L)+A(4,L))))/)
        coe(L,:)=(/0._DP, 8*e2**3*Q(2), 0._DP, 0._DP, 4*e2**2*(e2**2+Q(2)**2), 8*e2**3*Q(2), &
        0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 8*e2**3*Q(2), 0._DP, &
        4*e2**2*(e2**2+Q(2)**2), 0._DP, 4*e2*Q(1)*Q(6), &
        8*e2**3*Q(10), 0._DP, 2*e2**2*Q(3)*Q(13), 8*e2**3*Q(12), 0._DP, 0._DP, &
        2*e2**2*Q(3)*Q(11), 4*e2**2*Q(2)*Q(7), 0._DP, 2*e2**2*Q(3)*Q(14), 4*e2**3*Q(16), &
        0._DP, 2*e2**2*Q(3)*Q(11), 4*e2**2*Q(2)*Q(7), 0._DP, 0._DP, 8*e2**3*Q(10), 0._DP, 0._DP, &
        8*e2**2*Q(2)*Q(6), 0._DP, 2*e2**2*Q(3)*Q(13), 2*e2**2*Q(3)*Q(11), 2*e2*Q(1)*Q(7), &
        8*e2**2*Q(6)*Q(10), 0._DP, 0._DP, 4*Q(1)*Q(4)*Q(5), 8*e2**2*Q(8)*Q(9), 4*e2**2*Q(19), &
        2*e2*Q(3)*Q(6)*Q(11), 4*e2**2*Q(7)*Q(10), 0._DP, 2*e2**2*Q(17), e2*Q(3)*Q(20), &
        e2*Q(3)*Q(7)*Q(11), 4*e2**2*Q(15), 8*e2**2*Q(8)*Q(9), 8*e2*Q(2)*Q(4)*Q(5), &
        4*e2**2*Q(18), 0._DP, 0._DP, e2*Q(3)*Q(7)*Q(11), 8*e2*Q(6)*Q(8)*Q(9), &
        8*e2*Q(4)*Q(5)*Q(10), 0._DP, 2*e2*Q(21), 4*e2*Q(7)*Q(8)*Q(9), &
        2*Q(3)*Q(4)*Q(5)*Q(11), 8*Q(4)*Q(5)*Q(8)*Q(9)/)
      CASE (45_I2B)
        Q(1:21)=(/2*e2**2-2*e2*e5+e5**2, e3, 1-e1+2*e2, e2+A(1,L)-&
        e1*A(1,L), e2+A(2,L)-e1*A(2,L), e2-e5*(A(1,L)+A(2,L)), 2*e2-&
        e5*(A(1,L)+A(2,L)), e2+A(3,L)-e1*A(3,L), e2+A(4,L)-e1*A(4,L), e2-&
        e5*(A(3,L)+A(4,L)), 2*e2-e5*(A(3,L)+A(4,L)), e2-e5*(A(1,L)+A(2,L)+&
        A(3,L)+A(4,L)), 2*e2-e5*(A(1,L)+A(2,L)+A(3,L)+A(4,L)), 2*e2-&
        e5*(2*A(1,L)+2*A(2,L)+A(3,L)+A(4,L)), 2*e2**2-2*e2*e5*(A(1,L)+&
        A(2,L)+A(3,L)+A(4,L))+e5**2*((A(1,L)+A(2,L))*A(3,L)+(A(1,L)+&
        A(2,L)+2*A(3,L))*A(4,L)), 2*e2-e5*(A(1,L)+A(2,L)+2*(A(3,L)+&
        A(4,L))), 4*e2**2+e5**2*((A(1,L)+A(2,L))*A(3,L)+(A(1,L)+A(2,L)+&
        4*A(3,L))*A(4,L))-2*e2*e5*(A(1,L)+A(2,L)+2*(A(3,L)+A(4,L))), 2*e2**2 &
        - 2*e2*e5*(A(1,L)+A(2,L)+A(3,L)+A(4,L))+e5**2*(A(2,L)*(A(3,L)+&
        A(4,L))+A(1,L)*(2*A(2,L)+A(3,L)+A(4,L))), 2*e2**2-2*e2*e5*(A(1,L)+&
        A(2,L)+A(3,L)+A(4,L))+e5**2*(2*A(3,L)*A(4,L)+A(2,L)*(A(3,L)+A(4,L)) &
        + A(1,L)*(2*A(2,L)+A(3,L)+A(4,L))), 4*e2**2-2*e2*e5*(2*A(1,L)+&
        2*A(2,L)+A(3,L)+A(4,L))+e5**2*(A(2,L)*(A(3,L)+A(4,L))+&
        A(1,L)*(4*A(2,L)+A(3,L)+A(4,L))), 4*e2**3-4*e2**2*e5*(A(1,L)+A(2,L)+&
        A(3,L)+A(4,L))+e2*e5**2*(4*A(3,L)*A(4,L)+3*A(2,L)*(A(3,L)+A(4,L))+&
        A(1,L)*(4*A(2,L)+3*(A(3,L)+A(4,L))))-2*e5**3*(A(2,L)*A(3,L)*A(4,L)+&
        A(1,L)*(A(3,L)*A(4,L)+A(2,L)*(A(3,L)+A(4,L))))/)
        coe(L,:)=(/0._DP, 8*e2**3*Q(2), 0._DP, 0._DP, 8*e2**3*Q(2), 4*e2**2*(e2**2+Q(2)**2), &
        0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 8*e2**3*Q(2), 0._DP, &
        4*e2**2*(e2**2+Q(2)**2), 0._DP, 8*e2**3*Q(6), &
        4*e2*Q(1)*Q(10), 0._DP, 2*e2**2*Q(3)*Q(13), 8*e2**3*Q(12), 0._DP, 0._DP, &
        4*e2**2*Q(2)*Q(11), 2*e2**2*Q(3)*Q(7), 0._DP, 4*e2**3*Q(14), 2*e2**2*Q(3)*Q(16), &
        0._DP, 4*e2**2*Q(2)*Q(11), 2*e2**2*Q(3)*Q(7), 0._DP, 0._DP, 8*e2**2*Q(2)*Q(10), 0._DP, 0._DP, &
        8*e2**3*Q(6), 0._DP, 2*e2**2*Q(3)*Q(13), 2*e2*Q(1)*Q(11), 2*e2**2*Q(3)*Q(7), &
        8*e2**2*Q(6)*Q(10), 0._DP, 0._DP, 8*e2**2*Q(4)*Q(5), 4*Q(1)*Q(8)*Q(9), 4*e2**2*Q(19), &
        4*e2**2*Q(6)*Q(11), 2*e2*Q(3)*Q(7)*Q(10), 0._DP, e2*Q(3)*Q(17), 2*e2**2*Q(20), &
        e2*Q(3)*Q(7)*Q(11), 4*e2**2*Q(15), 8*e2*Q(2)*Q(8)*Q(9), 8*e2**2*Q(4)*Q(5), &
        4*e2**2*Q(18), 0._DP, 0._DP, e2*Q(3)*Q(7)*Q(11), 8*e2*Q(6)*Q(8)*Q(9), &
        8*e2*Q(4)*Q(5)*Q(10), 0._DP, 2*e2*Q(21), 2*Q(3)*Q(7)*Q(8)*Q(9), &
        4*e2*Q(4)*Q(5)*Q(11), 8*Q(4)*Q(5)*Q(8)*Q(9)/)
      CASE (46_I2B)
        Q(1:11)=(/e2+A(1,L)-e1*A(1,L), e2+A(2,L)-e1*A(2,L), e2-e5*(A(1,L) &
        + A(2,L)), 2*e2-e5*(A(1,L)+A(2,L)), 4*e2**2-4*e2*e5*(A(1,L)+A(2,L))+&
        e5**2*(A(1,L)+A(2,L)), 2*e2**2-3*e2*e5*(A(1,L)+A(2,L))+e5**2*(A(1,L) &
        + A(2,L)), e2**2-2*e2*e5*(A(1,L)+A(2,L))+e5**2*(A(1,L)+A(2,L)), &
        2*e2**3-2*e5**3*A(1,L)*A(2,L)-4*e2**2*e5*(A(1,L)+A(2,L))+&
        e2*e5**2*(A(1,L)+A(1,L)**2+A(2,L)+4*A(1,L)*A(2,L)+A(2,L)**2), 4*e2**3 &
        - 2*e5**3*A(1,L)*A(2,L)-6*e2**2*e5*(A(1,L)+A(2,L))+e2*e5**2*(A(1,L)+&
        A(1,L)**2+A(2,L)+6*A(1,L)*A(2,L)+A(2,L)**2), 2*e2**4-&
        4*e2*e5**3*A(1,L)*A(2,L)+e5**4*A(1,L)*A(2,L)-4*e2**3*e5*(A(1,L)+A(2,L)) &
        + e2**2*e5**2*(A(1,L)+A(1,L)**2+A(2,L)+6*A(1,L)*A(2,L)+A(2,L)**2), &
        4*e2**4-8*e2**3*e5*(A(1,L)+A(2,L))+e5**4*A(1,L)*A(2,L)*(A(1,L)+&
        A(2,L))-4*e2*e5**3*A(1,L)*A(2,L)*(1+A(1,L)+A(2,L))+&
        e2**2*e5**2*(A(1,L)+3*A(1,L)**2+A(2,L)+14*A(1,L)*A(2,L)+&
        3*A(2,L)**2)/)
        coe(L,:)=(/0._DP, 8*e2**4, 0._DP, 0._DP, 8*e2**4, 8*e2**4, &
        0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 8*e2**4, 0._DP, &
        8*e2**4, 0._DP, 8*e2**3*Q(3), 8*e2**3*Q(3), 0._DP, 2*e2**2*Q(5), 8*e2**2*Q(7), 0._DP, 0._DP, &
        4*e2**3*Q(4), 4*e2**3*Q(4), 0._DP, 4*e2**2*Q(6), 4*e2**2*Q(6), 0._DP, 4*e2**3*Q(4), &
        4*e2**3*Q(4), 0._DP, 0._DP, 8*e2**3*Q(3), 0._DP, 0._DP, 8*e2**3*Q(3), 0._DP, 2*e2**2*Q(5), &
        4*e2**3*Q(4), 4*e2**3*Q(4), 8*e2**2*Q(3)**2, 0._DP, 0._DP, 8*e2**2*Q(1)*Q(2), &
        8*e2**2*Q(1)*Q(2), 4*Q(10), 4*e2**2*Q(3)*Q(4), 4*e2**2*Q(3)*Q(4), 0._DP, &
        2*e2*Q(9), 2*e2*Q(9), 2*e2**2*(-2*e2+e5*(A(1,L)+A(2,L)))**2, 4*e2*Q(8), &
        8*e2**2*Q(1)*Q(2), 8*e2**2*Q(1)*Q(2), 4*e2*Q(8), 0._DP, 0._DP, 2*e2**2*(-2*e2+&
        e5*(A(1,L)+A(2,L)))**2, 8*e2*Q(1)*Q(2)*Q(3), 8*e2*Q(1)*Q(2)*Q(3), 0._DP, &
        2*Q(11), 4*e2*Q(1)*Q(2)*Q(4), 4*e2*Q(1)*Q(2)*Q(4), 8*Q(1)**2*Q(2)**2/)
      CASE (47_I2B)
        X(1)=Merge(A(4,L), A(3,L), ANY(G(3,L)==G(5:6,L)))
        X(2:3)=Merge(A(5:6,L), A(6:5:-1,L), ANY(G(5,L)==G(3:4,L)))
        Q(1:15)=(/1-e1+2*e2, e2+A(1,L)-e1*A(1,L), e2**2+e4*e5*A(1,L), e2+&
        x(1)-e1*x(1), e2-e5*(A(1,L)+x(1)), e2+x(2)-e1*x(2), e2-e5*(A(1,L) &
        + x(2)), e2-e5*(x(1)+x(2)), e2-e5*(A(1,L)+x(1)+x(2)), e2**2+&
        e5**2*A(1,L)-e2*e5*(2*A(1,L)+x(1)+x(2)), 4*e2**2+e5**2*x(1)-&
        e2*e5*(1+2*x(1)+2*x(2)), 4*e2**2+e5**2*(A(1,L)+x(1))-e2*e5*(1+&
        4*A(1,L)+2*x(1)+2*x(2)), 4*e2**3-e5**3*A(1,L)+e2*e5**2*(6*A(1,L)+&
        x(1))-e2**2*e5*(1+8*A(1,L)+2*x(1)+2*x(2)), 2*e2**2-2*e2*e5*(A(1,L) &
        + x(1)+x(2))+e5**2*(2*x(1)*x(2)+A(1,L)*(x(1)+x(2))), 2*e2**3-&
        e5**3*A(1,L)*(x(1)+x(2))-2*e2**2*e5*(2*A(1,L)+x(1)+x(2))+&
        2*e2*e5**2*(x(1)*x(2)+A(1,L)*(1+x(1)+x(2)))/)
        coe(L,:)=(/0._DP, 2*e2**3*Q(1), 0._DP, 0._DP, 2*e2**3*Q(1), &
        2*e2**3*Q(1), 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, &
        2*e2**3*Q(1), 0._DP, 2*e2**3*Q(1), 0._DP, 2*e2*Q(1)*Q(3), 4*e2**3*Q(8), 0._DP, &
        e2**2*Q(12), 4*e2**2*Q(10), 0._DP, 0._DP, e2**2*Q(11), 2*e2**2*Q(1)*Q(2), 0._DP, &
        e2*Q(13), 4*e2**3*Q(9), 0._DP, e2**2*Q(11), 2*e2**2*Q(1)*Q(2), 0._DP, 0._DP, &
        4*e2**3*Q(8), 0._DP, 0._DP, 2*e2*Q(1)*Q(3), 0._DP, e2**2*Q(12), e2**2*Q(11), &
        2*e2**2*Q(1)*Q(2), 4*e2*Q(3)*Q(8), 0._DP, 0._DP, 2*e2*Q(1)*Q(2)**2, &
        4*e2**2*Q(4)*Q(6), 4*e2**2*Q(5)*Q(7), Q(3)*Q(11), 4*e2**2*Q(2)*Q(8), 0._DP, &
        2*e2**2*Q(14), e2*Q(2)*Q(12), e2*Q(2)*Q(11), 2*e2*Q(15), 4*e2**2*Q(4)*Q(6), &
        2*e2*Q(1)*Q(2)**2, 4*e2**2*Q(2)*Q(9), 0._DP, 0._DP, e2*Q(2)*Q(11), 4*Q(3)*Q(4)*Q(6), &
        4*e2*Q(2)**2*Q(8), 0._DP, 2*e2*Q(2)*Q(14), 4*e2*Q(2)*Q(4)*Q(6), Q(2)**2*Q(11), &
        4*Q(2)**2*Q(4)*Q(6)/)
      CASE (48_I2B)
        X(3)=Merge(A(6,L), A(5,L), ANY(G(5,L)==G(1:2,L)))
        X(1:2)=Merge(A(1:2,L), A(2:1:-1,L), ANY(G(1,L)==G(5:6,L)))
        Q(1:15)=(/1-e1+2*e2, e2+A(3,L)-e1*A(3,L), e2**2+e4*e5*A(3,L), e2+&
        x(1)-e1*x(1), e2-e5*(A(3,L)+x(1)), e2+x(2)-e1*x(2), e2-e5*(A(3,L) &
        + x(2)), e2-e5*(x(1)+x(2)), e2-e5*(A(3,L)+x(1)+x(2)), e2**2+&
        e5**2*A(3,L)-e2*e5*(2*A(3,L)+x(1)+x(2)), 4*e2**2+e5**2*x(2)-&
        e2*e5*(1+2*x(1)+2*x(2)), 4*e2**2+e5**2*(A(3,L)+x(2))-e2*e5*(1+&
        4*A(3,L)+2*x(1)+2*x(2)), 4*e2**3-e5**3*A(3,L)+e2*e5**2*(6*A(3,L)+&
        x(2))-e2**2*e5*(1+8*A(3,L)+2*x(1)+2*x(2)), 2*e2**2-2*e2*e5*(A(3,L) &
        + x(1)+x(2))+e5**2*(2*x(1)*x(2)+A(3,L)*(x(1)+x(2))), 2*e2**3-&
        e5**3*A(3,L)*(x(1)+x(2))-2*e2**2*e5*(2*A(3,L)+x(1)+x(2))+&
        2*e2*e5**2*(x(1)*x(2)+A(3,L)*(1+x(1)+x(2)))/)
        coe(L,:)=(/0._DP, 2*e2**3*Q(1), 0._DP, 0._DP, 2*e2**3*Q(1), &
        2*e2**3*Q(1), 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, &
        2*e2**3*Q(1), 0._DP, 2*e2**3*Q(1), 0._DP, 4*e2**3*Q(8), 2*e2*Q(1)*Q(3), 0._DP, &
        e2**2*Q(12), 4*e2**2*Q(10), 0._DP, 0._DP, 2*e2**2*Q(1)*Q(2), e2**2*Q(11), 0._DP, &
        4*e2**3*Q(9), e2*Q(13), 0._DP, 2*e2**2*Q(1)*Q(2), e2**2*Q(11), 0._DP, 0._DP, &
        2*e2*Q(1)*Q(3), 0._DP, 0._DP, 4*e2**3*Q(8), 0._DP, e2**2*Q(12), 2*e2**2*Q(1)*Q(2), &
        e2**2*Q(11), 4*e2*Q(3)*Q(8), 0._DP, 0._DP, 4*e2**2*Q(4)*Q(6), 2*e2*Q(1)*Q(2)**2, &
        4*e2**2*Q(5)*Q(7), 4*e2**2*Q(2)*Q(8), Q(3)*Q(11), 0._DP, e2*Q(2)*Q(12), &
        2*e2**2*Q(14), e2*Q(2)*Q(11), 4*e2**2*Q(2)*Q(9), 2*e2*Q(1)*Q(2)**2, &
        4*e2**2*Q(4)*Q(6), 2*e2*(2*e2*Q(4)*Q(6)-e5*A(3,L)*Q(1)*(2*e2-e5*(x(1)+&
        x(2)))), 0._DP, 0._DP, e2*Q(2)*Q(11), 4*e2*Q(2)**2*Q(8), 4*Q(3)*Q(4)*Q(6), 0._DP, &
        2*e2*Q(2)*Q(14), Q(2)**2*Q(11), 4*e2*Q(2)*Q(4)*Q(6), 4*Q(2)**2*Q(4)*Q(6)/)
      CASE (49_I2B)
        X(3)=Merge(A(4,L), A(3,L), ANY(G(3,L)==G(1:2,L)))
        X(1:2)=Merge(A(1:2,L), A(2:1:-1,L), ANY(G(1,L)==G(3:4,L)))
        Q(1:17)=(/e2+x(1)-e1*x(1), e2+x(2)-e1*x(2), e2-e5*(x(1)+x(2)), &
        2*e2-e5*(x(1)+x(2)), e2+x(3)-e1*x(3), e2-e5*(x(1)+x(3)), 2*e2-&
        e5*(x(1)+x(3)), 4*e2**2+e5**2*x(1)-2*e2*e5*(2*x(1)+x(2)+x(3)), &
        e2**2+e5**2*x(1)-e2*e5*(2*x(1)+x(2)+x(3)), 2*e2**2+e5**2*x(1)-&
        e2*e5*(3*x(1)+2*x(2)+x(3)), 2*e2**2+e5**2*x(1)-e2*e5*(3*x(1)+x(2)+&
        2*x(3)), 2*e2**3-e5**3*x(1)*x(2)-2*e2**2*e5*(2*x(1)+x(2)+x(3))+&
        e2*e5**2*(x(1)*(1+x(1)+3*x(2))+(x(1)+x(2))*x(3)), 4*e2**3-&
        e5**3*x(1)*x(2)-2*e2**2*e5*(3*x(1)+2*x(2)+x(3))+e2*e5**2*(x(1)*(1+&
        x(1)+5*x(2))+(x(1)+x(2))*x(3)), 2*e2**3-e5**3*x(1)*x(3)-&
        2*e2**2*e5*(2*x(1)+x(2)+x(3))+e2*e5**2*(x(1)*(1+x(1)+x(2))+&
        (3*x(1)+x(2))*x(3)), 2*e2**3-e5**3*x(1)*(x(2)+x(3))-&
        2*e2**2*e5*(2*x(1)+x(2)+x(3))+e2*e5**2*(x(1)*(1+x(1)+3*x(2))+&
        (3*x(1)+x(2))*x(3)), 4*e2**3-e5**3*x(1)*x(3)-2*e2**2*e5*(3*x(1)+x(2) &
        + 2*x(3))+e2*e5**2*(x(1)*(1+x(1)+x(2))+(5*x(1)+x(2))*x(3)), 4*e2**4 &
        + e5**4*x(1)*x(2)*x(3)-4*e2**3*e5*(2*x(1)+x(2)+x(3))-&
        e2*e5**3*x(1)*(x(3)+2*x(1)*x(3)+x(2)*(1+2*x(1)+4*x(3)))+&
        e2**2*e5**2*(3*x(1)**2+3*x(2)*x(3)+x(1)*(1+7*x(2)+7*x(3)))/)
        coe(L,:)=(/4*e2**4, 4*e2**4*A(5,L), 2*e2**3*Q(7), 2*e2**3*Q(4), &
        4*e2**4*A(5,L), 4*e2**4*A(5,L), e2**2*Q(8), 4*e2**2*Q(9), 4*e2**3*Q(6), &
        4*e2**3*Q(3), 2*e2**2*Q(10), 4*e2**4*A(5,L), 2*e2**2*Q(11), 4*e2**4*A(5,L), &
        4*e2**2*Q(3)*Q(6), 4*e2**3*A(5,L)*Q(3), 4*e2**3*A(5,L)*Q(6), 2*e2*Q(15), &
        e2**2*A(5,L)*Q(8), 4*e2**2*A(5,L)*Q(9), 4*e2**2*Q(1)*Q(5), 4*e2**2*Q(1)*Q(2), &
        2*e2**3*A(5,L)*Q(7), 2*e2**3*A(5,L)*Q(4), e2**2*Q(4)*Q(7), &
        2*e2**2*A(5,L)*Q(10), 2*e2**2*A(5,L)*Q(11), e2*Q(16), 2*e2**3*A(5,L)*Q(7), &
        2*e2**3*A(5,L)*Q(4), e2*Q(13), 2*e2*Q(14), 4*e2**3*A(5,L)*Q(6), 2*e2*Q(12), &
        2*e2**2*Q(4)*Q(6), 4*e2**3*A(5,L)*Q(3), 2*e2**2*Q(3)*Q(7), e2**2*A(5,L)*Q(8), &
        2*e2**3*A(5,L)*Q(7), 2*e2**3*A(5,L)*Q(4), 4*e2**2*A(5,L)*Q(3)*Q(6), &
        4*e2*Q(1)*Q(3)*Q(5), 4*e2*Q(1)*Q(2)*Q(6), 4*e2**2*A(5,L)*Q(1)*Q(2), &
        4*e2**2*A(5,L)*Q(1)*Q(5), 2*e2*A(5,L)*Q(15), 2*e2**2*A(5,L)*Q(3)*Q(7), &
        2*e2**2*A(5,L)*Q(4)*Q(6), Q(17), e2*A(5,L)*Q(16), e2*A(5,L)*Q(13), &
        e2**2*A(5,L)*Q(4)*Q(7), 2*e2*A(5,L)*Q(14), 4*e2**2*A(5,L)*Q(1)*Q(5), &
        4*e2**2*A(5,L)*Q(1)*Q(2), 2*e2*A(5,L)*Q(12), 2*e2*Q(1)*Q(4)*Q(5), &
        2*e2*Q(1)*Q(2)*Q(7), e2**2*A(5,L)*Q(4)*Q(7), 4*e2*A(5,L)*Q(1)*Q(3)*Q(5), &
        4*e2*A(5,L)*Q(1)*Q(2)*Q(6), 4*Q(1)**2*Q(2)*Q(5), A(5,L)*Q(17), &
        2*e2*A(5,L)*Q(1)*Q(4)*Q(5), 2*e2*A(5,L)*Q(1)*Q(2)*Q(7), &
        4*A(5,L)*Q(1)**2*Q(2)*Q(5)/)
      CASE (50_I2B)
        X(4)=Merge(A(6,L), A(5,L), ANY(G(5,L)==G(1:2,L)))
        X(3)=Merge(A(4,L), A(3,L), ANY(G(3,L)==G(1:2,L)))
        X(1:2)=Merge(A(1:2,L), A(2:1:-1,L), ANY(G(1,L)==G(3:4,L)))
        Q(1:18)=(/1-e1+2*e2, e2+x(1)-e1*x(1), e2+x(2)-e1*x(2), e2-&
        e5*(x(1)+x(2)), 4*e2**2+e5**2*x(1)-e2*e5*(1+2*x(1)+2*x(2)), e2+&
        x(3)-e1*x(3), e2-e5*(x(1)+x(3)), 2*e2-e5*(x(1)+x(3)), e2**2+&
        e5**2*x(1)-e2*e5*(2*x(1)+x(2)+x(3)), 2*e2**2+e5**2*x(1)-&
        e2*e5*(3*x(1)+2*x(2)+x(3)), 8*e2**3-e5**3*x(1)+e2*e5**2*(5*x(1)+&
        x(3))-2*e2**2*e5*(1+4*x(1)+2*x(2)+2*x(3)), 4*e2**3-e5**3*x(1)+&
        e2*e5**2*(4*x(1)+x(3))-e2**2*e5*(1+6*x(1)+2*x(2)+4*x(3)), 2*e2**3-&
        e5**3*x(1)*x(2)-2*e2**2*e5*(2*x(1)+x(2)+x(3))+e2*e5**2*(x(1)*(1+&
        x(1)+3*x(2))+(x(1)+x(2))*x(3)), 4*e2**3-e5**3*x(1)*x(2)-&
        2*e2**2*e5*(3*x(1)+2*x(2)+x(3))+e2*e5**2*(x(1)*(1+x(1)+5*x(2))+&
        (x(1)+x(2))*x(3)), 2*e2**3-e5**3*x(1)*x(3)-2*e2**2*e5*(2*x(1)+x(2)+&
        x(3))+e2*e5**2*(x(1)*(1+x(1)+x(2))+(3*x(1)+x(2))*x(3)), 2*e2**3-&
        e5**3*x(1)*(x(2)+x(3))-2*e2**2*e5*(2*x(1)+x(2)+x(3))+&
        e2*e5**2*(x(1)*(1+x(1)+3*x(2))+(3*x(1)+x(2))*x(3)), 8*e2**4+&
        e5**4*x(1)*x(3)-2*e2**3*e5*(1+6*x(1)+2*x(2)+4*x(3))-&
        e2*e5**3*x(1)*(1+x(1)+5*x(3))+2*e2**2*e5**2*(x(1)*(3+x(1)+x(2))+&
        (1+5*x(1)+x(2))*x(3)), 4*e2**4+e5**4*x(1)*x(2)*x(3)-&
        4*e2**3*e5*(2*x(1)+x(2)+x(3))-e2*e5**3*x(1)*(x(3)+2*x(1)*x(3)+&
        x(2)*(1+2*x(1)+4*x(3)))+e2**2*e5**2*(3*x(1)**2+3*x(2)*x(3)+x(1)*(1 &
        + 7*x(2)+7*x(3)))/)
        coe(L,:)=(/0._DP, 4*e2**3*Q(1), 0._DP, 0._DP, 4*e2**3*Q(1), &
        4*e2**3*Q(1), 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, &
        4*e2**3*Q(1), 0._DP, 4*e2**3*Q(1), 0._DP, 8*e2**3*Q(4), 4*e2**2*Q(1)*Q(7), 0._DP, &
        e2*Q(11), 8*e2**2*Q(9), 0._DP, 0._DP, 2*e2**2*Q(1)*Q(8), 2*e2**2*Q(5), 0._DP, &
        4*e2**2*Q(10), 2*e2*Q(12), 0._DP, 2*e2**2*Q(1)*Q(8), 2*e2**2*Q(5), 0._DP, 0._DP, &
        4*e2**2*Q(1)*Q(7), 0._DP, 0._DP, 8*e2**3*Q(4), 0._DP, e2*Q(11), 2*e2**2*Q(1)*Q(8), &
        2*e2**2*Q(5), 8*e2**2*Q(4)*Q(7), 0._DP, 0._DP, 8*e2**2*Q(2)*Q(3), &
        4*e2*Q(1)*Q(2)*Q(6), 4*e2*Q(16), 4*e2**2*Q(4)*Q(8), 2*e2*Q(5)*Q(7), 0._DP, Q(17), &
        2*e2*Q(14), e2*Q(5)*Q(8), 4*e2*Q(15), 4*e2*Q(1)*Q(2)*Q(6), 8*e2**2*Q(2)*Q(3), &
        4*e2*Q(13), 0._DP, 0._DP, e2*Q(5)*Q(8), 8*e2*Q(2)*Q(4)*Q(6), 8*e2*Q(2)*Q(3)*Q(7), 0._DP, &
        2*Q(18), 2*Q(2)*Q(5)*Q(6), 4*e2*Q(2)*Q(3)*Q(8), 8*Q(2)**2*Q(3)*Q(6)/)
      CASE (51_I2B)
        X(1:3:2)=Merge(A(1:2,L), A(2:1:-1,L), ANY(G(1,L)==G(3:4,L)))
        X(2)=Merge(A(4,L), A(3,L), ANY(G(4,L)==G(5:6,L)))
        Q(1:18)=(/1-e1+2*e2, e2+x(1)-e1*x(1), e2+x(2)-e1*x(2), e2-&
        e5*(x(1)+x(2)), 4*e2**2+e5**2*x(1)-e2*e5*(1+2*x(1)+2*x(2)), e2+&
        x(3)-e1*x(3), e2-e5*(x(1)+x(3)), 2*e2-e5*(x(1)+x(3)), e2**2+&
        e5**2*x(1)-e2*e5*(2*x(1)+x(2)+x(3)), 2*e2**2+e5**2*x(1)-&
        e2*e5*(3*x(1)+2*x(2)+x(3)), 8*e2**3-e5**3*x(1)+e2*e5**2*(5*x(1)+&
        x(3))-2*e2**2*e5*(1+4*x(1)+2*x(2)+2*x(3)), 4*e2**3-e5**3*x(1)+&
        e2*e5**2*(4*x(1)+x(3))-e2**2*e5*(1+6*x(1)+2*x(2)+4*x(3)), 2*e2**3-&
        e5**3*x(1)*x(2)-2*e2**2*e5*(2*x(1)+x(2)+x(3))+e2*e5**2*(x(1)*(1+&
        x(1)+3*x(2))+(x(1)+x(2))*x(3)), 4*e2**3-e5**3*x(1)*x(2)-&
        2*e2**2*e5*(3*x(1)+2*x(2)+x(3))+e2*e5**2*(x(1)*(1+x(1)+5*x(2))+&
        (x(1)+x(2))*x(3)), 2*e2**3-e5**3*x(1)*x(3)-2*e2**2*e5*(2*x(1)+x(2)+&
        x(3))+e2*e5**2*(x(1)*(1+x(1)+x(2))+(3*x(1)+x(2))*x(3)), 2*e2**3-&
        e5**3*x(1)*(x(2)+x(3))-2*e2**2*e5*(2*x(1)+x(2)+x(3))+&
        e2*e5**2*(x(1)*(1+x(1)+3*x(2))+(3*x(1)+x(2))*x(3)), 8*e2**4+&
        e5**4*x(1)*x(3)-2*e2**3*e5*(1+6*x(1)+2*x(2)+4*x(3))-&
        e2*e5**3*x(1)*(1+x(1)+5*x(3))+2*e2**2*e5**2*(x(1)*(3+x(1)+x(2))+&
        (1+5*x(1)+x(2))*x(3)), 4*e2**4+e5**4*x(1)*x(2)*x(3)-&
        4*e2**3*e5*(2*x(1)+x(2)+x(3))-e2*e5**3*x(1)*(x(3)+2*x(1)*x(3)+&
        x(2)*(1+2*x(1)+4*x(3)))+e2**2*e5**2*(3*x(1)**2+3*x(2)*x(3)+x(1)*(1 &
        + 7*x(2)+7*x(3)))/)
        coe(L,:)=(/0._DP, 4*e2**3*Q(1), 0._DP, 0._DP, 4*e2**3*Q(1), &
        4*e2**3*Q(1), 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, &
        4*e2**3*Q(1), 0._DP, 4*e2**3*Q(1), 0._DP, 4*e2**2*Q(1)*Q(7), 8*e2**3*Q(4), 0._DP, &
        e2*Q(11), 8*e2**2*Q(9), 0._DP, 0._DP, 2*e2**2*Q(5), 2*e2**2*Q(1)*Q(8), 0._DP, 2*e2*Q(12), &
        4*e2**2*Q(10), 0._DP, 2*e2**2*Q(5), 2*e2**2*Q(1)*Q(8), 0._DP, 0._DP, 8*e2**3*Q(4), 0._DP, 0._DP, &
        4*e2**2*Q(1)*Q(7), 0._DP, e2*Q(11), 2*e2**2*Q(5), 2*e2**2*Q(1)*Q(8), &
        8*e2**2*Q(4)*Q(7), 0._DP, 0._DP, 4*e2*Q(1)*Q(2)*Q(6), 8*e2**2*Q(2)*Q(3), &
        4*e2*(2*e2**3-e5**3*x(1)*(x(2)+x(3))-2*e2**2*e5*(2*x(1)+x(2)+x(3)) &
        + e2*e5**2*(x(2)*(3*x(1)+x(3))+x(1)*(1+x(1)+3*x(3)))), &
        2*e2*Q(5)*Q(7), 4*e2**2*Q(4)*Q(8), 0._DP, 2*e2*(4*e2**3-e5**3*x(1)*x(2)-&
        2*e2**2*e5*(3*x(1)+2*x(2)+x(3))+e2*e5**2*(x(1)*(1+x(1)+x(3))+&
        x(2)*(5*x(1)+x(3)))), 8*e2**4+e5**4*x(1)*x(3)-2*e2**3*e5*(1+6*x(1)+&
        2*x(2)+4*x(3))-e2*e5**3*x(1)*(1+x(1)+5*x(3))+2*e2**2*e5**2*(x(3)+&
        x(2)*(x(1)+x(3))+x(1)*(3+x(1)+5*x(3))), e2*Q(5)*Q(8), 4*e2*(2*e2**3-&
        e5**3*x(1)*x(2)-2*e2**2*e5*(2*x(1)+x(2)+x(3))+e2*e5**2*(x(1)*(1+&
        x(1)+x(3))+x(2)*(3*x(1)+x(3)))), 8*e2**2*Q(2)*Q(3), &
        4*e2*Q(1)*Q(2)*Q(6), 4*e2*(2*e2**3-e5**3*x(1)*x(3)-2*e2**2*e5*(2*x(1)+&
        x(2)+x(3))+e2*e5**2*(x(2)*(x(1)+x(3))+x(1)*(1+x(1)+3*x(3)))), 0._DP, &
        0._DP, e2*Q(5)*Q(8), 8*e2*Q(2)*Q(3)*Q(7), 8*e2*Q(2)*Q(4)*Q(6), 0._DP, 2*(4*e2**4+&
        e5**4*x(1)*x(2)*x(3)-4*e2**3*e5*(2*x(1)+x(2)+x(3))-&
        e2*e5**3*x(1)*(x(2)+2*x(1)*x(2)+(1+2*x(1)+4*x(2))*x(3))+&
        e2**2*e5**2*(3*x(1)**2+3*x(2)*x(3)+x(1)*(1+7*x(2)+7*x(3)))), &
        4*e2*Q(2)*Q(3)*Q(8), 2*Q(2)*Q(5)*Q(6), 8*Q(2)**2*Q(3)*Q(6)/)
      CASE (52_I2B)
        IF(ALL(G(1:3:2,L)==G(5:6,L)).OR.ALL(G(1:3:2,L)==G(6:5:-1,L)))THEN
          X(1:2)=A(1:3:2,L)
          X(3:4)=A(2:4:2,L)
        ELSE IF(ALL(G(1:4:3,L)==G(5:6,L)).OR.ALL(G(1:4:3,L)==G(6:5:-1,L)))THEN
          X(1:2)=A(1:4:3,L)
          X(3:4)=A(2:3,L)
        ELSE IF(ALL(G(2:3,L)==G(5:6,L)).OR.ALL(G(2:3,L)==G(6:5:-1,L)))THEN
          X(1:2)=A(2:3,L)
          X(3:4)=A(1:4:3,L)
        ELSE
          X(1:2)=A(2:4:2,L)
          X(3:4)=A(1:3:2,L)
        END IF
        Q(1:28)=(/2*e2**2-2*e2*e5+e5**2, e3, 1-e1+2*e2, e2+x(1)-e1*x(1), &
        e2+x(2)-e1*x(2), e2+x(3)-e1*x(3), e2-e5*(x(1)+x(3)), 4*e2**2+&
        e5**2*x(3)-e2*e5*(1+2*x(1)+2*x(3)), 4*e2**2+e5**2*(x(1)+2*x(3))-&
        e2*e5*(3+2*x(1)+2*x(3)), 4*e2**3-e5**3*x(3)+e2*e5**2*(1+x(1)+&
        2*x(3))-e2**2*e5*(3+2*x(1)+2*x(3)), e2+x(4)-e1*x(4), e2-e5*(x(2) &
        + x(4)), e2-e5*(x(1)+x(2)+x(3)+x(4)), 4*e2**2+e5**2*x(4)-e2*e5*(1 &
        + 2*x(2)+2*x(4)), 4*e2**2+e5**2*(x(2)+2*x(4))-e2*e5*(3+2*x(2)+&
        2*x(4)), 4*e2**3-e5**3*x(4)+e2*e5**2*(1+x(2)+2*x(4))-e2**2*e5*(3+&
        2*x(2)+2*x(4)), 4*e2**2+e5**2*(x(1)+x(3)+x(4))-e2*e5*(1+4*x(1)+&
        2*x(2)+4*x(3)+2*x(4)), 8*e2**3-e5**3*(x(3)+x(4))-4*e2**2*e5*(1+&
        x(1)+x(2)+x(3)+x(4))+e2*e5**2*(1+x(1)+x(2)+3*x(3)+3*x(4)), &
        4*e2**2+e5**2*(x(2)+x(3)+x(4))-e2*e5*(1+2*x(1)+4*x(2)+2*x(3)+&
        4*x(4)), 2*e2**2-2*e2*e5*(x(1)+x(2)+x(3)+x(4))+e5**2*(x(2)*(x(1)+&
        x(3))+(x(1)+2*x(2)+x(3))*x(4)), 8*e2**2-4*e2*e5*(1+x(1)+x(2)+&
        x(3)+x(4))+e5**2*(x(1)+x(2)+3*(x(3)+x(4))), 2*e2**2-2*e2*e5*(x(1) &
        + x(2)+x(3)+x(4))+e5**2*(x(2)*x(3)+(2*x(2)+x(3))*x(4)+x(1)*(x(2) &
        + 2*x(3)+x(4))), 2*e2**2-2*e2*e5*(x(1)+x(2)+x(3)+x(4))+&
        e5**2*(x(3)*(x(2)+x(4))+x(1)*(x(2)+2*x(3)+x(4))), 8*e2**3-&
        2*e2**2*e5*(1+4*x(1)+2*x(2)+4*x(3)+2*x(4))-e5**3*(2*x(1)*x(3)+&
        (x(1)+x(3))*x(4))+2*e2*e5**2*(x(4)+x(3)*(1+x(2)+x(4))+x(1)*(1+&
        x(2)+4*x(3)+x(4))), 8*e2**4+e5**4*x(3)*x(4)-4*e2**3*e5*(1+x(1)+&
        x(2)+x(3)+x(4))-e2*e5**3*((1+x(1))*x(4)+x(3)*(1+x(2)+2*x(4)))+&
        e2**2*e5**2*(1+x(2)+3*x(3)+3*x(4)+2*x(3)*(x(2)+x(4))+x(1)*(1+&
        2*x(2)+2*x(4))), 8*e2**3-4*e2**2*e5*(1+x(1)+x(2)+x(3)+x(4))-&
        e5**3*(x(2)*x(3)+(x(1)+2*x(3))*x(4))+e2*e5**2*(x(2)+2*x(3)*(x(2)+&
        x(4))+3*(x(3)+x(4))+x(1)*(1+2*x(2)+2*x(4))), 8*e2**3-&
        2*e2**2*e5*(1+2*x(1)+4*x(2)+2*x(3)+4*x(4))-e5**3*(x(3)*x(4)+&
        x(2)*(x(3)+2*x(4)))+2*e2*e5**2*(x(3)+(1+x(1)+x(3))*x(4)+x(2)*(1+&
        x(1)+x(3)+4*x(4))), 4*e2**3-4*e2**2*e5*(x(1)+x(2)+x(3)+x(4))+&
        e2*e5**2*(3*x(2)*x(3)+4*x(2)*x(4)+3*x(3)*x(4)+x(1)*(3*x(2)+4*x(3)+&
        3*x(4)))-2*e5**3*(x(2)*x(3)*x(4)+x(1)*(x(2)*x(3)+(x(2)+&
        x(3))*x(4)))/)
        coe(L,:)=(/0._DP, 8*e2**3*Q(2), 0._DP, 0._DP, 2*e2**2*e4**2, &
        2*e2**2*e4**2, 0._DP, 0._DP, 0._DP, 0._DP, &
        0._DP, 4*e2**2*(e2**2+Q(2)**2), 0._DP, 2*e2**2*e4**2, 0._DP, 4*e2**2*Q(3)*Q(7), &
        4*e2**2*Q(3)*Q(12), 0._DP, e2*Q(18), 8*e2**3*Q(13), 0._DP, 0._DP, 2*e2**2*Q(15), &
        2*e2**2*Q(9), 0._DP, 2*e2**2*Q(17), 2*e2**2*Q(19), 0._DP, 2*e2*Q(16), 2*e2*Q(10), 0._DP, &
        0._DP, 4*e2**2*Q(3)*Q(12), 0._DP, 0._DP, 4*e2**2*Q(3)*Q(7), 0._DP, e2**2*(8*e2**2-&
        4*e2*e5*(1+x(1)+x(2)+x(3)+x(4))+e5**2*(x(1)+x(2)+3*x(3)+&
        3*x(4))), e2*Q(3)*Q(14), e2*Q(3)*Q(8), 8*e2**2*Q(7)*Q(12), 0._DP, 0._DP, &
        4*e2*Q(3)*Q(4)*Q(6), 4*e2*Q(3)*Q(5)*Q(11), 4*e2**2*(2*e2**2-2*e2*e5*(x(1)+&
        x(2)+x(3)+x(4))+e5**2*(2*x(2)*x(4)+x(3)*(x(2)+x(4))+x(1)*(x(2)+&
        2*x(3)+x(4)))), 2*e2*Q(7)*Q(14), 2*e2*Q(8)*Q(12), 0._DP, e2*(8*e2**3-&
        2*e2**2*e5*(1+2*x(1)+4*x(2)+2*x(3)+4*x(4))-e5**3*(2*x(2)*x(4)+&
        x(3)*(x(2)+x(4)))+2*e2*e5**2*(4*x(2)*x(4)+(1+x(1))*(x(2)+x(4))+&
        x(3)*(1+x(2)+x(4)))), e2*(2*(1-e1+4*e2)*Q(4)*Q(6)-&
        2*e2*e5*x(2)*(2*e2-e5*(x(1)+x(3)))-e5*Q(3)*(2*e2-e5*(x(1)+&
        x(3)))*x(4)), 8*e2**4+e5**4*x(3)*x(4)-4*e2**3*e5*(1+x(1)+x(2)+x(3) &
        + x(4))-e2*e5**3*((1+x(1))*x(4)+x(3)*(1+x(2)+2*x(4)))+&
        e2**2*e5**2*(1+x(2)+3*x(4)+x(1)*(1+2*x(2)+2*x(4))+x(3)*(3+&
        2*x(2)+2*x(4))), 4*e2**2*Q(20), 4*e2*Q(3)*Q(5)*Q(11), 4*e2*Q(3)*Q(4)*Q(6), &
        4*e2**2*Q(23), 0._DP, 0._DP, e2*(8*e2**3-4*e2**2*e5*(1+x(1)+x(2)+x(3)+x(4)) &
        - e5**3*(x(1)*x(4)+x(3)*(x(2)+2*x(4)))+e2*e5**2*(x(2)+3*x(4)+&
        x(1)*(1+2*x(2)+2*x(4))+x(3)*(3+2*x(2)+2*x(4)))), &
        8*e2*Q(5)*Q(7)*Q(11), 8*e2*Q(4)*Q(6)*Q(12), 0._DP, 2*e2*(4*e2**3-&
        4*e2**2*e5*(x(1)+x(2)+x(3)+x(4))+e2*e5**2*(4*x(2)*x(4)+3*x(3)*(x(2) &
        + x(4))+x(1)*(4*x(3)+3*(x(2)+x(4))))-2*e5**3*(x(2)*x(3)*x(4)+&
        x(1)*(x(2)*x(4)+x(3)*(x(2)+x(4))))), 2*Q(5)*Q(8)*Q(11), &
        2*Q(4)*Q(6)*Q(14), 8*Q(4)*Q(5)*Q(6)*Q(11)/)
      CASE (53_I2B)
        QQ=Merge(A(4,L), A(3,L), G(1,L)==G(3,L))
        Q(1:11)=(/e2+QQ-e1*QQ, e2+A(1,L)-e1*A(1,L), e2**2+e4*e5*A(1,L), e2 &
        - e5*(QQ+A(1,L)), 2*e2-e5*(QQ+A(1,L)), e2**2+e5**2*A(1,L)-e2*e5*(QQ &
        + 2*A(1,L)), 2*e2**3-e5**3*QQ*A(1,L)+e2*e5**2*A(1,L)*(1+3*QQ+A(1,L)) &
        - 2*e2**2*e5*(QQ+2*A(1,L)), 2*e2**2+e5**2*A(1,L)-e2*e5*(QQ+3*A(1,L)), &
        e2**3+3*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(QQ+3*A(1,L)), 2*e2**4 &
        + e5**4*QQ*A(1,L)+2*e2**2*e5**2*A(1,L)*(2+2*QQ+A(1,L))-&
        e2*e5**3*A(1,L)*(1+3*QQ+A(1,L))-2*e2**3*e5*(QQ+3*A(1,L)), 2*e2**3+&
        4*e2*e5**2*A(1,L)-e5**3*A(1,L)-e2**2*e5*(QQ+5*A(1,L))/)
        coe(L,:)=(/0._DP, 4*e2**4, 0._DP, 0._DP, 4*e2**4, 4*e2**4, &
        0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 4*e2**4, 0._DP, &
        4*e2**4, 0._DP, 4*e2**2*Q(3), 4*e2**3*Q(4), 0._DP, 2*e2**2*Q(8), 4*e2*Q(9), 0._DP, 0._DP, &
        2*e2**3*Q(5), 4*e2**3*Q(2), 0._DP, 2*e2*Q(11), 4*e2**2*Q(6), 0._DP, 2*e2**3*Q(5), &
        4*e2**3*Q(2), 0._DP, 0._DP, 4*e2**3*Q(4), 0._DP, 0._DP, 4*e2**2*Q(3), 0._DP, 2*e2**2*Q(8), &
        2*e2**3*Q(5), 4*e2**3*Q(2), 4*e2*Q(3)*Q(4), 0._DP, 0._DP, 4*e2**2*Q(2)**2, &
        4*e2**2*Q(1)*Q(2), 4*e2*Q(3)*Q(4), 2*e2*Q(3)*Q(5), 4*e2**2*Q(2)*Q(4), 0._DP, &
        2*e2*Q(7), 2*e2*Q(2)*Q(8), 2*e2**2*Q(2)*Q(5), 2*Q(10), 4*e2**2*Q(1)*Q(2), &
        4*e2**2*Q(2)**2, 4*e2*Q(2)*Q(6), 0._DP, 0._DP, 2*e2**2*Q(2)*Q(5), 4*Q(1)*Q(2)*Q(3), &
        4*e2*Q(2)**2*Q(4), 0._DP, 2*Q(2)*Q(7), 4*e2*Q(1)*Q(2)**2, 2*e2*Q(2)**2*Q(5), &
        4*Q(1)*Q(2)**3/)
      CASE (54_I2B)
        QQ=Merge(A(6,L), A(5,L), G(1,L)==G(5,L))
        Q(1:17)=(/2*e2**2-2*e2*e5+e5**2, e3, 1-e1+2*e2, e2+A(1,L)-&
        e1*A(1,L), e2**2+e4*e5*A(1,L), e2+A(3,L)-e1*A(3,L), e2-e5*(A(1,L)+&
        A(3,L)), e2+A(4,L)-e1*A(4,L), e2-e5*(A(1,L)+A(4,L)), e2-e5*(A(3,L) &
        + A(4,L)), 2*e2-e5*(A(3,L)+A(4,L)), e2-e5*(A(1,L)+A(3,L)+A(4,L)), &
        2*e2-e5*(2*A(1,L)+A(3,L)+A(4,L)), e2**2+e5**2*A(1,L)-&
        e2*e5*(2*A(1,L)+A(3,L)+A(4,L)), 2*e2**2+2*e5**2*A(1,L)-&
        e2*e5*(4*A(1,L)+A(3,L)+A(4,L)), 2*e2**2-2*e2*e5*(A(1,L)+A(3,L)+&
        A(4,L))+e5**2*(2*A(3,L)*A(4,L)+A(1,L)*(A(3,L)+A(4,L))), 2*e2**3-&
        e5**3*A(1,L)*(A(3,L)+A(4,L))-2*e2**2*e5*(2*A(1,L)+A(3,L)+A(4,L))+&
        2*e2*e5**2*(A(3,L)*A(4,L)+A(1,L)*(1+A(3,L)+A(4,L)))/)
        coe(L,:)=(/0._DP, 2*e2**2*(e2**2+Q(2)**2), 0._DP, 0._DP, 2*e2**2*(e2**2+Q(2)**2), &
        4*e2**3*Q(2), 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, &
        2*e2**2*(e2**2+Q(2)**2), 0._DP, 4*e2**3*Q(2), 0._DP, &
        4*e2**2*Q(5), 4*e2**2*Q(2)*Q(10), 0._DP, e2**2*Q(3)*Q(13), 4*e2**2*Q(14), 0._DP, 0._DP, &
        e2*Q(1)*Q(11), 2*e2**2*Q(3)*Q(4), 0._DP, 2*e2**2*Q(15), 2*e2**2*Q(3)*Q(12), 0._DP, &
        e2*Q(1)*Q(11), 2*e2**2*Q(3)*Q(4), 0._DP, 0._DP, 2*e2*Q(1)*Q(10), 0._DP, 0._DP, 4*e2**2*Q(5), &
        0._DP, e2**2*Q(3)*Q(13), 2*e2**2*Q(2)*Q(11), 2*e2**2*Q(3)*Q(4), 4*e2*Q(5)*Q(10), &
        0._DP, 0._DP, 4*e2**2*Q(4)**2, 4*e2*Q(2)*Q(6)*Q(8), 4*e2**2*Q(7)*Q(9), &
        2*e2*Q(5)*Q(11), 2*e2*Q(3)*Q(4)*Q(10), 0._DP, e2*Q(3)*Q(16), 2*e2**2*Q(4)*Q(13), &
        e2*Q(3)*Q(4)*Q(11), 2*e2*Q(17), 2*Q(1)*Q(6)*Q(8), 4*e2**2*Q(4)**2, &
        4*e2**2*Q(4)*Q(12), 0._DP, 0._DP, e2*Q(3)*Q(4)*Q(11), 4*Q(5)*Q(6)*Q(8), &
        4*e2*Q(4)**2*Q(10), 0._DP, 2*e2*Q(4)*Q(16), 2*Q(3)*Q(4)*Q(6)*Q(8), &
        2*e2*Q(4)**2*Q(11), 4*Q(4)**2*Q(6)*Q(8)/)
      CASE (55_I2B)
        QQ=Merge(A(6,L), A(5,L), G(3,L)==G(5,L))
        Q(1:17)=(/2*e2**2-2*e2*e5+e5**2, e3, 1-e1+2*e2, e2+A(1,L)-&
        e1*A(1,L), e2+A(2,L)-e1*A(2,L), e2-e5*(A(1,L)+A(2,L)), 2*e2-&
        e5*(A(1,L)+A(2,L)), e2+A(3,L)-e1*A(3,L), e2**2+e4*e5*A(3,L), 2*e2*(e2 &
        + A(1,L)-e1*A(1,L))*(e2+A(2,L)-e1*A(2,L))-(1-e1+2*e2)*e5*(2*e2-&
        e5*(A(1,L)+A(2,L)))*A(3,L), e2-e5*(A(1,L)+A(3,L)), e2-e5*(A(2,L)+&
        A(3,L)), e2-e5*(A(1,L)+A(2,L)+A(3,L)), 2*e2-e5*(A(1,L)+A(2,L)+&
        2*A(3,L)), e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+A(2,L)+2*A(3,L)), 2*e2**2 &
        + 2*e5**2*A(3,L)-e2*e5*(A(1,L)+A(2,L)+4*A(3,L)), 2*e2**2-&
        2*e2*e5*(A(1,L)+A(2,L)+A(3,L))+e5**2*(2*A(1,L)*A(2,L)+(A(1,L)+&
        A(2,L))*A(3,L))/)
        coe(L,:)=(/0._DP, 2*e2**2*(e2**2+Q(2)**2), 0._DP, 0._DP, 4*e2**3*Q(2), 2*e2**2*(e2**2 &
        + Q(2)**2), 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 2*e2**2*(e2**2+Q(2)**2), 0._DP, 4*e2**3*Q(2), 0._DP, &
        4*e2**2*Q(2)*Q(6), 4*e2**2*Q(9), 0._DP, e2**2*Q(3)*Q(14), 4*e2**2*Q(15), 0._DP, 0._DP, &
        2*e2**2*Q(3)*Q(8), e2*Q(1)*Q(7), 0._DP, 2*e2**2*Q(3)*Q(13), 2*e2**2*Q(16), 0._DP, &
        2*e2**2*Q(3)*Q(8), e2*Q(1)*Q(7), 0._DP, 0._DP, 4*e2**2*Q(9), 0._DP, 0._DP, 2*e2*Q(1)*Q(6), 0._DP, &
        e2**2*Q(3)*Q(14), 2*e2**2*Q(3)*Q(8), 2*e2**2*Q(2)*Q(7), 4*e2*Q(6)*Q(9), 0._DP, 0._DP, &
        4*e2*Q(2)*Q(4)*Q(5), 4*e2**2*Q(8)**2, 4*e2**2*Q(11)*Q(12), &
        2*e2*Q(3)*Q(6)*Q(8), 2*e2*Q(7)*Q(9), 0._DP, 2*e2**2*Q(8)*Q(14), e2*Q(3)*Q(17), &
        e2*Q(3)*Q(7)*Q(8), 4*e2**2*Q(8)*Q(13), 4*e2**2*Q(8)**2, 2*Q(1)*Q(4)*Q(5), &
        2*e2*Q(10), 0._DP, 0._DP, e2*Q(3)*Q(7)*Q(8), 4*e2*Q(6)*Q(8)**2, 4*Q(4)*Q(5)*Q(9), 0._DP, &
        2*e2*Q(8)*Q(17), 2*e2*Q(7)*Q(8)**2, 2*Q(3)*Q(4)*Q(5)*Q(8), &
        4*Q(4)*Q(5)*Q(8)**2/)
      CASE (56_I2B)
        QQ=Merge(A(2,L), A(1,L), G(3,L)==G(1,L))
        Q(1:11)=(/e2+QQ-e1*QQ, e2+A(3,L)-e1*A(3,L), e2**2+e4*e5*A(3,L), e2 &
        - e5*(QQ+A(3,L)), 2*e2-e5*(QQ+A(3,L)), e2**2+e5**2*A(3,L)-e2*e5*(QQ &
        + 2*A(3,L)), 2*e2**3-e5**3*QQ*A(3,L)+e2*e5**2*A(3,L)*(1+3*QQ+A(3,L)) &
        - 2*e2**2*e5*(QQ+2*A(3,L)), 2*e2**2+e5**2*A(3,L)-e2*e5*(QQ+3*A(3,L)), &
        e2**3+3*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(QQ+3*A(3,L)), 2*e2**4 &
        + e5**4*QQ*A(3,L)+2*e2**2*e5**2*A(3,L)*(2+2*QQ+A(3,L))-&
        e2*e5**3*A(3,L)*(1+3*QQ+A(3,L))-2*e2**3*e5*(QQ+3*A(3,L)), 2*e2**3+&
        4*e2*e5**2*A(3,L)-e5**3*A(3,L)-e2**2*e5*(QQ+5*A(3,L))/)
        coe(L,:)=(/0._DP, 4*e2**4, 0._DP, 0._DP, 4*e2**4, 4*e2**4, &
        0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 4*e2**4, 0._DP, &
        4*e2**4, 0._DP, 4*e2**3*Q(4), 4*e2**2*Q(3), 0._DP, 2*e2**2*Q(8), 4*e2*Q(9), 0._DP, 0._DP, &
        4*e2**3*Q(2), 2*e2**3*Q(5), 0._DP, 4*e2**2*Q(6), 2*e2*Q(11), 0._DP, 4*e2**3*Q(2), &
        2*e2**3*Q(5), 0._DP, 0._DP, 4*e2**2*Q(3), 0._DP, 0._DP, 4*e2**3*Q(4), 0._DP, 2*e2**2*Q(8), &
        4*e2**3*Q(2), 2*e2**3*Q(5), 4*e2*Q(3)*Q(4), 0._DP, 0._DP, 4*e2**2*Q(1)*Q(2), &
        4*e2**2*Q(2)**2, 4*e2*Q(3)*Q(4), 4*e2**2*Q(2)*Q(4), 2*e2*Q(3)*Q(5), 0._DP, &
        2*e2*Q(2)*Q(8), 2*e2*Q(7), 2*e2**2*Q(2)*Q(5), 4*e2*Q(2)*Q(6), &
        4*e2**2*Q(2)**2, 4*e2**2*Q(1)*Q(2), 2*Q(10), 0._DP, 0._DP, 2*e2**2*Q(2)*Q(5), &
        4*e2*Q(2)**2*Q(4), 4*Q(1)*Q(2)*Q(3), 0._DP, 2*Q(2)*Q(7), 2*e2*Q(2)**2*Q(5), &
        4*e2*Q(1)*Q(2)**2, 4*Q(1)*Q(2)**3/)
      CASE (57_I2B)
        QQ=Merge(A(2,L), A(1,L), G(5,L)==G(1,L))
        Q(1:19)=(/e3, e2+QQ-e1*QQ, e2+A(3,L)-e1*A(3,L), e2+A(4,L)-&
        e1*A(4,L), e2-e5*(A(3,L)+A(4,L)), 2*e2-e5*(A(3,L)+A(4,L)), e2+&
        A(5,L)-e1*A(5,L), e2-e5*(QQ+A(5,L)), 2*e2**2+e5**2*QQ-e2*e5*(1+QQ &
        + A(5,L)), e2-e5*(QQ+A(3,L)+A(4,L)+A(5,L)), 4*e2**2+e5**2*(2*QQ+&
        A(3,L)+A(4,L))-2*e2*e5*(1+QQ+A(3,L)+A(4,L)+A(5,L)), 2*e2**2+&
        e5**2*(QQ+A(3,L)+A(4,L))-e2*e5*(1+QQ+2*A(3,L)+2*A(4,L)+A(5,L)), &
        2*e2-e5*(2*QQ+A(3,L)+A(4,L)+2*A(5,L)), 2*e2**2-2*e2*e5*(QQ+A(3,L) &
        + A(4,L)+A(5,L))+e5**2*(2*A(3,L)*A(4,L)+QQ*(A(3,L)+A(4,L))+(A(3,L) &
        + A(4,L))*A(5,L)), 2*e2**2-2*e2*e5*(QQ+A(3,L)+A(4,L)+A(5,L))+&
        e5**2*(QQ*(A(3,L)+A(4,L))+(2*QQ+A(3,L)+A(4,L))*A(5,L)), 4*e2**2-&
        2*e2*e5*(2*QQ+A(3,L)+A(4,L)+2*A(5,L))+e5**2*(QQ*(A(3,L)+A(4,L))+&
        (4*QQ+A(3,L)+A(4,L))*A(5,L)), 4*e2**3-4*e2**2*e5*(QQ+A(3,L)+A(4,L) &
        + A(5,L))+e2*e5**2*(4*A(3,L)*A(4,L)+3*QQ*(A(3,L)+A(4,L))+4*QQ*A(5,L) &
        + 3*(A(3,L)+A(4,L))*A(5,L))-2*e5**3*(QQ*A(3,L)*A(4,L)+(QQ*A(3,L)+(QQ &
        + A(3,L))*A(4,L))*A(5,L)), 4*e2**3-e5**3*(2*A(3,L)*A(4,L)+QQ*(A(3,L)+&
        A(4,L)))-2*e2**2*e5*(1+QQ+2*A(3,L)+2*A(4,L)+A(5,L))+&
        e2*e5**2*(QQ*(2+A(3,L)+A(4,L))+A(4,L)*(2+A(5,L))+A(3,L)*(2+&
        4*A(4,L)+A(5,L))), 2*e2**2-2*e2*e5*(QQ+A(3,L)+A(4,L)+A(5,L))+&
        e5**2*(2*A(3,L)*A(4,L)+(A(3,L)+A(4,L))*A(5,L)+QQ*(A(3,L)+A(4,L)+&
        2*A(5,L)))/)
        coe(L,:)=(/4*e2**3*Q(1), 4*e2**3*A(5,L)*Q(1), 2*e2**2*Q(1)*Q(6), &
        2*e2**2*Q(9), 4*e2**3*A(5,L)*Q(1), 4*e2**3*A(5,L)*Q(1), e2**2*Q(11), &
        4*e2**3*Q(10), 4*e2**2*Q(1)*Q(5), 4*e2**3*Q(8), 2*e2**3*Q(13), &
        4*e2**3*A(5,L)*Q(1), 2*e2**2*Q(12), 4*e2**3*A(5,L)*Q(1), 4*e2**2*Q(5)*Q(8), &
        4*e2**3*A(5,L)*Q(8), 4*e2**2*A(5,L)*Q(1)*Q(5), 2*e2**2*Q(19), &
        e2**2*A(5,L)*Q(11), 4*e2**3*A(5,L)*Q(10), 4*e2*Q(1)*Q(3)*Q(4), &
        4*e2**2*Q(2)*Q(7), 2*e2**2*A(5,L)*Q(1)*Q(6), 2*e2**2*A(5,L)*Q(9), &
        e2*Q(6)*Q(9), 2*e2**3*A(5,L)*Q(13), 2*e2**2*A(5,L)*Q(12), e2*Q(18), &
        2*e2**2*A(5,L)*Q(1)*Q(6), 2*e2**2*A(5,L)*Q(9), e2**2*(4*e2**2-2*e2*e5*(2*QQ &
        + A(3,L)+A(4,L)+2*A(5,L))+e5**2*((A(3,L)+A(4,L))*A(5,L)+QQ*(A(3,L) &
        + A(4,L)+4*A(5,L)))), 2*e2**2*Q(14), 4*e2**2*A(5,L)*Q(1)*Q(5), &
        2*e2**2*(2*e2**2-2*e2*e5*(QQ+A(3,L)+A(4,L)+A(5,L))+e5**2*((A(3,L)+&
        A(4,L))*A(5,L)+QQ*(A(3,L)+A(4,L)+2*A(5,L)))), 2*e2*Q(5)*Q(9), &
        4*e2**3*A(5,L)*Q(8), 2*e2**2*Q(6)*Q(8), e2**2*A(5,L)*Q(11), &
        2*e2**2*A(5,L)*Q(1)*Q(6), 2*e2**2*A(5,L)*Q(9), 4*e2**2*A(5,L)*Q(5)*Q(8), &
        4*e2*Q(3)*Q(4)*Q(8), 4*e2*Q(2)*Q(5)*Q(7), 4*e2**2*A(5,L)*Q(2)*Q(7), &
        4*e2*A(5,L)*Q(1)*Q(3)*Q(4), 2*e2**2*A(5,L)*Q(19), 2*e2**2*A(5,L)*Q(6)*Q(8), &
        2*e2*A(5,L)*Q(5)*Q(9), e2*(4*e2**3-4*e2**2*e5*(QQ+A(3,L)+A(4,L)+&
        A(5,L))+e2*e5**2*(4*A(3,L)*A(4,L)+3*(A(3,L)+A(4,L))*A(5,L)+&
        QQ*(3*(A(3,L)+A(4,L))+4*A(5,L)))-2*e5**3*(A(3,L)*A(4,L)*A(5,L)+&
        QQ*(A(3,L)*A(4,L)+(A(3,L)+A(4,L))*A(5,L)))), e2*A(5,L)*Q(18), &
        e2**2*A(5,L)*(4*e2**2-2*e2*e5*(2*QQ+A(3,L)+A(4,L)+2*A(5,L))+&
        e5**2*((A(3,L)+A(4,L))*A(5,L)+QQ*(A(3,L)+A(4,L)+4*A(5,L)))), &
        e2*A(5,L)*Q(6)*Q(9), 2*e2**2*A(5,L)*Q(14), 4*e2*A(5,L)*Q(1)*Q(3)*Q(4), &
        4*e2**2*A(5,L)*Q(2)*Q(7), 2*e2**2*A(5,L)*(2*e2**2-2*e2*e5*(QQ+A(3,L)+&
        A(4,L)+A(5,L))+e5**2*((A(3,L)+A(4,L))*A(5,L)+QQ*(A(3,L)+A(4,L)+&
        2*A(5,L)))), 2*Q(3)*Q(4)*Q(9), 2*e2*Q(2)*Q(6)*Q(7), e2*A(5,L)*Q(6)*Q(9), &
        4*e2*A(5,L)*Q(3)*Q(4)*Q(8), 4*e2*A(5,L)*Q(2)*Q(5)*Q(7), &
        4*Q(2)*Q(3)*Q(4)*Q(7), e2*A(5,L)*(4*e2**3-4*e2**2*e5*(QQ+A(3,L)+A(4,L) &
        + A(5,L))+e2*e5**2*(4*A(3,L)*A(4,L)+3*(A(3,L)+A(4,L))*A(5,L)+&
        QQ*(3*(A(3,L)+A(4,L))+4*A(5,L)))-2*e5**3*(A(3,L)*A(4,L)*A(5,L)+&
        QQ*(A(3,L)*A(4,L)+(A(3,L)+A(4,L))*A(5,L)))), 2*A(5,L)*Q(3)*Q(4)*Q(9), &
        2*e2*A(5,L)*Q(2)*Q(6)*Q(7), 4*A(5,L)*Q(2)*Q(3)*Q(4)*Q(7)/)
      CASE (58_I2B)
        QQ=Merge(A(4,L), A(3,L), G(5,L)==G(3,L))
        Q(1:19)=(/e3, e2+QQ-e1*QQ, e2+A(1,L)-e1*A(1,L), e2+A(2,L)-&
          e1*A(2,L), e2-e5*(A(1,L)+A(2,L)), 2*e2-e5*(A(1,L)+A(2,L)), e2+&
          A(5,L)-e1*A(5,L), e2-e5*(QQ+A(5,L)), 2*e2**2+e5**2*QQ-e2*e5*(1+QQ &
          + A(5,L)), e2-e5*(QQ+A(1,L)+A(2,L)+A(5,L)), 4*e2**2+e5**2*(2*QQ+&
          A(1,L)+A(2,L))-2*e2*e5*(1+QQ+A(1,L)+A(2,L)+A(5,L)), 2*e2**2+&
          e5**2*(QQ+A(1,L)+A(2,L))-e2*e5*(1+QQ+2*A(1,L)+2*A(2,L)+A(5,L)), &
          2*e2-e5*(2*QQ+A(1,L)+A(2,L)+2*A(5,L)), 2*e2**2-2*e2*e5*(QQ+A(1,L) &
          + A(2,L)+A(5,L))+e5**2*(2*A(1,L)*A(2,L)+QQ*(A(1,L)+A(2,L))+(A(1,L) &
          + A(2,L))*A(5,L)), 2*e2**2-2*e2*e5*(QQ+A(1,L)+A(2,L)+A(5,L))+&
          e5**2*(QQ*(A(1,L)+A(2,L))+(2*QQ+A(1,L)+A(2,L))*A(5,L)), 4*e2**2-&
          2*e2*e5*(2*QQ+A(1,L)+A(2,L)+2*A(5,L))+e5**2*(QQ*(A(1,L)+A(2,L))+&
          (4*QQ+A(1,L)+A(2,L))*A(5,L)), 4*e2**3-4*e2**2*e5*(QQ+A(1,L)+A(2,L) &
          + A(5,L))+e2*e5**2*(4*A(1,L)*A(2,L)+3*QQ*(A(1,L)+A(2,L))+4*QQ*A(5,L) &
          + 3*(A(1,L)+A(2,L))*A(5,L))-2*e5**3*(QQ*A(1,L)*A(2,L)+(QQ*A(1,L)+(QQ &
          + A(1,L))*A(2,L))*A(5,L)), 4*e2**3-e5**3*(2*A(1,L)*A(2,L)+QQ*(A(1,L)+&
          A(2,L)))-2*e2**2*e5*(1+QQ+2*A(1,L)+2*A(2,L)+A(5,L))+&
          e2*e5**2*(QQ*(2+A(1,L)+A(2,L))+A(2,L)*(2+A(5,L))+A(1,L)*(2+&
          4*A(2,L)+A(5,L))), 2*e2**2-2*e2*e5*(QQ+A(1,L)+A(2,L)+A(5,L))+&
          e5**2*(2*A(1,L)*A(2,L)+(A(1,L)+A(2,L))*A(5,L)+QQ*(A(1,L)+A(2,L)+&
          2*A(5,L)))/)
        coe(L,:)=(/4*e2**3*Q(1), 4*e2**3*A(5,L)*Q(1), 2*e2**2*Q(9), &
        2*e2**2*Q(1)*Q(6), 4*e2**3*A(5,L)*Q(1), 4*e2**3*A(5,L)*Q(1), e2**2*Q(11), &
        4*e2**3*Q(10), 4*e2**3*Q(8), 4*e2**2*Q(1)*Q(5), 2*e2**2*Q(12), &
        4*e2**3*A(5,L)*Q(1), 2*e2**3*(2*e2-e5*(A(1,L)+A(2,L)+2*(QQ+A(5,L)))), &
        4*e2**3*A(5,L)*Q(1), 4*e2**2*Q(5)*Q(8), 4*e2**2*A(5,L)*Q(1)*Q(5), &
        4*e2**3*A(5,L)*Q(8), 2*e2**2*(2*e2**2-2*e2*e5*(QQ+A(1,L)+A(2,L)+&
        A(5,L))+e5**2*(2*QQ*A(5,L)+A(2,L)*(QQ+A(5,L))+A(1,L)*(QQ+2*A(2,L)+&
        A(5,L)))), e2**2*A(5,L)*Q(11), 4*e2**3*A(5,L)*Q(10), 4*e2**2*Q(2)*Q(7), &
        4*e2*Q(1)*Q(3)*Q(4), 2*e2**2*A(5,L)*Q(9), 2*e2**2*A(5,L)*Q(1)*Q(6), &
        e2*Q(6)*Q(9), 2*e2**2*A(5,L)*Q(12), 2*e2**3*A(5,L)*(2*e2-e5*(A(1,L)+&
        A(2,L)+2*(QQ+A(5,L)))), e2**2*(4*e2**2+e5**2*(QQ*(A(1,L)+A(2,L))+&
        (4*QQ+A(1,L)+A(2,L))*A(5,L))-2*e2*e5*(A(1,L)+A(2,L)+2*(QQ+&
        A(5,L)))), 2*e2**2*A(5,L)*Q(9), 2*e2**2*A(5,L)*Q(1)*Q(6), e2*(4*e2**3-&
        e5**3*(2*A(1,L)*A(2,L)+QQ*(A(1,L)+A(2,L)))-2*e2**2*e5*(1+QQ+&
        2*A(1,L)+2*A(2,L)+A(5,L))+e2*e5**2*(2*QQ+A(2,L)*(2+QQ+A(5,L))+&
        A(1,L)*(2+QQ+4*A(2,L)+A(5,L)))), 2*e2**2*Q(15), 4*e2**3*A(5,L)*Q(8), &
        2*e2**2*(2*e2**2-2*e2*e5*(QQ+A(1,L)+A(2,L)+A(5,L))+&
        e5**2*(A(2,L)*(QQ+A(5,L))+A(1,L)*(QQ+2*A(2,L)+A(5,L)))), &
        2*e2**2*Q(6)*Q(8), 4*e2**2*A(5,L)*Q(1)*Q(5), 2*e2*Q(5)*Q(9), &
        e2**2*A(5,L)*Q(11), 2*e2**2*A(5,L)*Q(9), 2*e2**2*A(5,L)*Q(1)*Q(6), &
        4*e2**2*A(5,L)*Q(5)*Q(8), 4*e2*Q(2)*Q(5)*Q(7), 4*e2*Q(3)*Q(4)*Q(8), &
        4*e2*A(5,L)*Q(1)*Q(3)*Q(4), 4*e2**2*A(5,L)*Q(2)*Q(7), 2*e2**2*A(5,L)*(2*e2**2 &
        - 2*e2*e5*(QQ+A(1,L)+A(2,L)+A(5,L))+e5**2*(2*QQ*A(5,L)+A(2,L)*(QQ+&
        A(5,L))+A(1,L)*(QQ+2*A(2,L)+A(5,L)))), 2*e2*A(5,L)*Q(5)*Q(9), &
        2*e2**2*A(5,L)*Q(6)*Q(8), e2*(4*e2**3-4*e2**2*e5*(QQ+A(1,L)+A(2,L)+&
        A(5,L))+e2*e5**2*(4*QQ*A(5,L)+3*A(2,L)*(QQ+A(5,L))+A(1,L)*(4*A(2,L)+&
        3*(QQ+A(5,L))))-2*e5**3*(QQ*A(2,L)*A(5,L)+A(1,L)*(QQ*A(5,L)+&
        A(2,L)*(QQ+A(5,L))))), e2**2*A(5,L)*(4*e2**2+e5**2*(QQ*(A(1,L)+A(2,L)) &
        + (4*QQ+A(1,L)+A(2,L))*A(5,L))-2*e2*e5*(A(1,L)+A(2,L)+2*(QQ+&
        A(5,L)))), e2*A(5,L)*(4*e2**3-e5**3*(2*A(1,L)*A(2,L)+QQ*(A(1,L)+&
        A(2,L)))-2*e2**2*e5*(1+QQ+2*A(1,L)+2*A(2,L)+A(5,L))+&
        e2*e5**2*(2*QQ+A(2,L)*(2+QQ+A(5,L))+A(1,L)*(2+QQ+4*A(2,L)+&
        A(5,L)))), e2*A(5,L)*Q(6)*Q(9), 2*e2**2*A(5,L)*Q(15), &
        4*e2**2*A(5,L)*Q(2)*Q(7), 4*e2*A(5,L)*Q(1)*Q(3)*Q(4), &
        2*e2**2*A(5,L)*(2*e2**2 &
        - 2*e2*e5*(QQ+A(1,L)+A(2,L)+A(5,L))+e5**2*(A(2,L)*(QQ+A(5,L))+&
        A(1,L)*(QQ+2*A(2,L)+A(5,L)))), 2*e2*Q(2)*Q(6)*Q(7), 2*Q(3)*Q(4)*Q(9), &
        e2*A(5,L)*Q(6)*Q(9), 4*e2*A(5,L)*Q(2)*Q(5)*Q(7), 4*e2*A(5,L)*Q(3)*Q(4)*Q(8), &
        4*Q(2)*Q(3)*Q(4)*Q(7), e2*A(5,L)*(4*e2**3-4*e2**2*e5*(QQ+A(1,L)+A(2,L) &
        + A(5,L))+e2*e5**2*(4*QQ*A(5,L)+3*A(2,L)*(QQ+A(5,L))+A(1,L)*(4*A(2,L) &
        + 3*(QQ+A(5,L))))-2*e5**3*(QQ*A(2,L)*A(5,L)+A(1,L)*(QQ*A(5,L)+&
        A(2,L)*(QQ+A(5,L))))), 2*e2*A(5,L)*Q(2)*Q(6)*Q(7), 2*A(5,L)*Q(3)*Q(4)*Q(9), &
        4*A(5,L)*Q(2)*Q(3)*Q(4)*Q(7)/)
      CASE (59_I2B)
        X(1)=Merge(A(5,L), A(6,L), ANY(G(5,L)==G(1:2,L)))
        X(2)=Merge(A(2,L), A(1,L), ANY(G(5:6,L)==G(1,L)))
        X(3)=Merge(A(4,L), A(3,L), ANY(G(5:6,L)==G(3,L)))
        Q(1:27)=(/2*e2**2-2*e2*e5+e5**2, e3, 1-e1+2*e2, e2+x(1)-e1*x(1), &
        e2+x(2)-e1*x(2), e2-e5*(x(1)+x(2)), 4*e2**2+e5**2*x(2)-e2*e5*(1+&
        2*x(1)+2*x(2)), 4*e2**2+e5**2*(x(1)+2*x(2))-e2*e5*(3+2*x(1)+&
        2*x(2)), 4*e2**3-e5**3*x(2)+e2*e5**2*(1+x(1)+2*x(2))-e2**2*e5*(3+&
        2*x(1)+2*x(2)), e2+x(3)-e1*x(3), e2-e5*(x(1)+x(3)), e2**2+&
        e5**2*x(1)-e2*e5*(2*x(1)+x(2)+x(3)), 4*e2**2+e5**2*x(3)-e2*e5*(1+&
        2*x(1)+2*x(3)), 4*e2**2+e5**2*(x(1)+2*x(3))-e2*e5*(3+2*x(1)+&
        2*x(3)), 4*e2**3-e5**3*x(3)+e2*e5**2*(1+x(1)+2*x(3))-e2**2*e5*(3+&
        2*x(1)+2*x(3)), 4*e2**2+e5**2*(3*x(1)+x(2)+x(3))-e2*e5*(1+6*x(1) &
        + 4*x(2)+2*x(3)), 8*e2**3-e5**3*(x(2)+x(3))-4*e2**2*e5*(1+2*x(1)+&
        x(2)+x(3))+e2*e5**2*(1+4*x(1)+3*x(2)+3*x(3)), 4*e2**2+&
        e5**2*(3*x(1)+x(2)+x(3))-e2*e5*(1+6*x(1)+2*x(2)+4*x(3)), 2*e2**3 &
        - e5**3*x(1)*x(2)-2*e2**2*e5*(2*x(1)+x(2)+x(3))+e2*e5**2*(x(1)*(1+&
        x(1)+3*x(2))+(x(1)+x(2))*x(3)), 8*e2**3-2*e2**2*e5*(1+6*x(1)+&
        4*x(2)+2*x(3))-e5**3*(4*x(1)*x(2)+(x(1)+x(2))*x(3))+&
        2*e2*e5**2*(x(2)+x(1)*(2+x(1)+5*x(2))+x(3)+(x(1)+x(2))*x(3)), &
        2*e2**3-e5**3*x(1)*x(3)-2*e2**2*e5*(2*x(1)+x(2)+x(3))+&
        e2*e5**2*(x(1)*(1+x(1)+x(2))+(3*x(1)+x(2))*x(3)), 2*e2**3-&
        e5**3*x(1)*(x(2)+x(3))-2*e2**2*e5*(2*x(1)+x(2)+x(3))+&
        e2*e5**2*(x(1)*(1+x(1)+3*x(2))+(3*x(1)+x(2))*x(3)), 8*e2**2-&
        4*e2*e5*(1+2*x(1)+x(2)+x(3))+e5**2*(4*x(1)+3*(x(2)+x(3))), &
        8*e2**3-4*e2**2*e5*(1+2*x(1)+x(2)+x(3))+e2*e5**2*(2*x(1)*(1+x(1) &
        + x(2))+2*(x(1)+x(2))*x(3)+3*(x(2)+x(3)))-e5**3*(2*x(2)*x(3)+&
        x(1)*(x(2)+x(3))), 8*e2**4+e5**4*x(2)*x(3)-4*e2**3*e5*(1+2*x(1)+&
        x(2)+x(3))+e2**2*e5**2*(1+3*x(2)+2*x(1)*(1+x(1)+x(2))+3*x(3)+&
        2*(x(1)+x(2))*x(3))-e2*e5**3*((1+x(1))*x(3)+x(2)*(1+x(1)+&
        2*x(3))), 8*e2**3-2*e2**2*e5*(1+6*x(1)+2*x(2)+4*x(3))+&
        2*e2*e5**2*(x(2)+x(1)*(2+x(1)+x(2))+x(3)+(5*x(1)+x(2))*x(3))-&
        e5**3*(x(2)*x(3)+x(1)*(x(2)+4*x(3))), 4*e2**4+e5**4*x(1)*x(2)*x(3)-&
        4*e2**3*e5*(2*x(1)+x(2)+x(3))-e2*e5**3*x(1)*(x(3)+2*x(1)*x(3)+&
        x(2)*(1+2*x(1)+4*x(3)))+e2**2*e5**2*(3*x(1)**2+3*x(2)*x(3)+x(1)*(1 &
        + 7*x(2)+7*x(3)))/)
        coe(L,:)=(/0._DP, 4*e2**2*(e2**2+Q(2)**2), 0._DP, 0._DP, 2*e2**2*e4**2, 2*e2**2*e4**2, &
        0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 8*e2**3*Q(2), 0._DP, 2*e2**2*e4**2, 0._DP, 4*e2**2*Q(3)*Q(6), &
        4*e2**2*Q(3)*Q(11), 0._DP, e2**2*Q(23), 8*e2**2*Q(12), 0._DP, 0._DP, 2*e2*Q(15), &
        2*e2*Q(9), 0._DP, 2*e2**2*Q(16), 2*e2**2*Q(18), 0._DP, 2*e2**2*Q(14), 2*e2**2*Q(8), &
        0._DP, 0._DP, 4*e2**2*Q(3)*Q(11), 0._DP, 0._DP, 4*e2**2*Q(3)*Q(6), 0._DP, e2*Q(17), &
        e2*Q(3)*Q(13), e2*Q(3)*Q(7), 8*e2**2*Q(6)*Q(11), 0._DP, 0._DP, 4*e2*Q(3)*Q(4)*Q(5), &
        4*e2*Q(3)*Q(4)*Q(10), 4*e2*(2*e2**3-e5**3*x(1)*(x(2)+x(3))-&
        2*e2**2*e5*(2*x(1)+x(2)+x(3))+e2*e5**2*(x(1)+x(1)**2+x(2)*x(3)+&
        3*x(1)*(x(2)+x(3)))), 2*e2*Q(6)*Q(13), 2*e2*Q(7)*Q(11), 0._DP, e2*(8*e2**3-&
        2*e2**2*e5*(1+6*x(1)+2*x(2)+4*x(3))-e5**3*(4*x(1)*x(3)+x(2)*(x(1)+&
        x(3)))+2*e2*e5**2*(x(1)**2+x(3)+x(2)*(1+x(1)+x(3))+x(1)*(2+&
        5*x(3)))), e2*(8*e2**3-2*e2**2*e5*(1+6*x(1)+4*x(2)+2*x(3))-&
        e5**3*(x(1)*x(3)+x(2)*(4*x(1)+x(3)))+2*e2*e5**2*(x(1)**2+x(3)+&
        x(1)*(2+x(3))+x(2)*(1+5*x(1)+x(3)))), e2*(8*e2**3-4*e2**2*e5*(1+&
        2*x(1)+x(2)+x(3))-e5**3*(2*x(2)*x(3)+x(1)*(x(2)+x(3)))+&
        e2*e5**2*(2*x(1)**2+3*x(2)+3*x(3)+2*x(2)*x(3)+2*x(1)*(1+x(2)+&
        x(3)))), 4*e2*(2*e2**3-e5**3*x(1)*x(3)-2*e2**2*e5*(2*x(1)+x(2)+x(3)) &
        + e2*e5**2*(x(2)*(x(1)+x(3))+x(1)*(1+x(1)+3*x(3)))), &
        4*e2*Q(3)*Q(4)*Q(10), 4*e2*Q(3)*Q(4)*Q(5), 4*e2*(2*e2**3-e5**3*x(1)*x(2)-&
        2*e2**2*e5*(2*x(1)+x(2)+x(3))+e2*e5**2*(x(1)*(1+x(1)+x(3))+&
        x(2)*(3*x(1)+x(3)))), 0._DP, 0._DP, 8*e2**4+e5**4*x(2)*x(3)-4*e2**3*e5*(1+&
        2*x(1)+x(2)+x(3))+e2**2*e5**2*(1+2*x(1)**2+3*x(2)+3*x(3)+&
        2*x(2)*x(3)+2*x(1)*(1+x(2)+x(3)))-e2*e5**3*((1+x(1))*x(3)+x(2)*(1 &
        + x(1)+2*x(3))), 8*e2*Q(4)*Q(6)*Q(10), 8*e2*Q(4)*Q(5)*Q(11), 0._DP, 2*(4*e2**4 &
        + e5**4*x(1)*x(2)*x(3)-4*e2**3*e5*(2*x(1)+x(2)+x(3))-&
        e2*e5**3*x(1)*(x(2)+x(3)+4*x(2)*x(3)+2*x(1)*(x(2)+x(3)))+&
        e2**2*e5**2*(x(1)+3*x(1)**2+3*x(2)*x(3)+7*x(1)*(x(2)+x(3)))), &
        2*Q(4)*Q(7)*Q(10), 2*Q(4)*Q(5)*Q(13), 8*Q(4)**2*Q(5)*Q(10)/)
      CASE (60_I2B)
        Q(1:14)=(/e2+A(1,L)-e1*A(1,L), e2**2+e4*e5*A(1,L), e2+A(3,L)-&
        e1*A(3,L), e2-e5*(A(1,L)+A(3,L)), e2+A(4,L)-e1*A(4,L), e2-&
        e5*(A(1,L)+A(4,L)), e2-e5*(A(3,L)+A(4,L)), 2*e2-e5*(A(3,L)+A(4,L)), &
        e2-e5*(A(1,L)+A(3,L)+A(4,L)), 2*e2-e5*(2*A(1,L)+A(3,L)+A(4,L)), &
        e2**2+e5**2*A(1,L)-e2*e5*(2*A(1,L)+A(3,L)+A(4,L)), 2*e2**2+&
        2*e5**2*A(1,L)-e2*e5*(4*A(1,L)+A(3,L)+A(4,L)), 2*e2**2-&
        2*e2*e5*(A(1,L)+A(3,L)+A(4,L))+e5**2*(2*A(3,L)*A(4,L)+A(1,L)*(A(3,L) &
        + A(4,L))), 2*e2**3-e5**3*A(1,L)*(A(3,L)+A(4,L))-2*e2**2*e5*(2*A(1,L)+&
        A(3,L)+A(4,L))+2*e2*e5**2*(A(3,L)*A(4,L)+A(1,L)*(1+A(3,L)+&
        A(4,L)))/)
        coe(L,:)=(/0._DP, 4*e2**4, 0._DP, 0._DP, 4*e2**4, 4*e2**4, &
        0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 4*e2**4, 0._DP, &
        4*e2**4, 0._DP, 4*e2**2*Q(2), 4*e2**3*Q(7), 0._DP, 2*e2**3*Q(10), 4*e2**2*Q(11), 0._DP, &
        0._DP, 2*e2**3*Q(8), 4*e2**3*Q(1), 0._DP, 2*e2**2*Q(12), 4*e2**3*Q(9), 0._DP, &
        2*e2**3*Q(8), 4*e2**3*Q(1), 0._DP, 0._DP, 4*e2**3*Q(7), 0._DP, 0._DP, 4*e2**2*Q(2), 0._DP, &
        2*e2**3*Q(10), 2*e2**3*Q(8), 4*e2**3*Q(1), 4*e2*Q(2)*Q(7), 0._DP, 0._DP, &
        4*e2**2*Q(1)**2, 4*e2**2*Q(3)*Q(5), 4*e2**2*Q(4)*Q(6), 2*e2*Q(2)*Q(8), &
        4*e2**2*Q(1)*Q(7), 0._DP, 2*e2**2*Q(13), 2*e2**2*Q(1)*Q(10), 2*e2**2*Q(1)*Q(8), &
        2*e2*Q(14), 4*e2**2*Q(3)*Q(5), 4*e2**2*Q(1)**2, 4*e2**2*Q(1)*Q(9), 0._DP, 0._DP, &
        2*e2**2*Q(1)*Q(8), 4*Q(2)*Q(3)*Q(5), 4*e2*Q(1)**2*Q(7), 0._DP, 2*e2*Q(1)*Q(13), &
        4*e2*Q(1)*Q(3)*Q(5), 2*e2*Q(1)**2*Q(8), 4*Q(1)**2*Q(3)*Q(5)/)
      CASE (61_I2B)
        Q(1:14)=(/e2+A(1,L)-e1*A(1,L), e2+A(2,L)-e1*A(2,L), e2-e5*(A(1,L) &
        + A(2,L)), 2*e2-e5*(A(1,L)+A(2,L)), e2+A(3,L)-e1*A(3,L), e2**2+&
        e4*e5*A(3,L), 2*e2*(e2+A(1,L)-e1*A(1,L))*(e2+A(2,L)-e1*A(2,L))-(1-&
        e1+2*e2)*e5*(2*e2-e5*(A(1,L)+A(2,L)))*A(3,L), e2-e5*(A(1,L)+&
        A(3,L)), e2-e5*(A(2,L)+A(3,L)), e2-e5*(A(1,L)+A(2,L)+A(3,L)), 2*e2 &
        - e5*(A(1,L)+A(2,L)+2*A(3,L)), e2**2+e5**2*A(3,L)-e2*e5*(A(1,L)+&
        A(2,L)+2*A(3,L)), 2*e2**2+2*e5**2*A(3,L)-e2*e5*(A(1,L)+A(2,L)+&
        4*A(3,L)), 2*e2**2-2*e2*e5*(A(1,L)+A(2,L)+A(3,L))+&
        e5**2*(2*A(1,L)*A(2,L)+(A(1,L)+A(2,L))*A(3,L))/)
        coe(L,:)=(/0._DP, 4*e2**4, 0._DP, 0._DP, 4*e2**4, 4*e2**4, &
        0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 4*e2**4, 0._DP, &
        4*e2**4, 0._DP, 4*e2**3*Q(3), 4*e2**2*Q(6), 0._DP, 2*e2**3*Q(11), 4*e2**2*Q(12), 0._DP, &
        0._DP, 4*e2**3*Q(5), 2*e2**3*Q(4), 0._DP, 4*e2**3*Q(10), 2*e2**2*Q(13), 0._DP, &
        4*e2**3*Q(5), 2*e2**3*Q(4), 0._DP, 0._DP, 4*e2**2*Q(6), 0._DP, 0._DP, 4*e2**3*Q(3), 0._DP, &
        2*e2**3*Q(11), 4*e2**3*Q(5), 2*e2**3*Q(4), 4*e2*Q(3)*Q(6), 0._DP, 0._DP, &
        4*e2**2*Q(1)*Q(2), 4*e2**2*Q(5)**2, 4*e2**2*Q(8)*Q(9), 4*e2**2*Q(3)*Q(5), &
        2*e2*Q(4)*Q(6), 0._DP, 2*e2**2*Q(5)*Q(11), 2*e2**2*Q(14), 2*e2**2*Q(4)*Q(5), &
        4*e2**2*Q(5)*Q(10), 4*e2**2*Q(5)**2, 4*e2**2*Q(1)*Q(2), 2*e2*Q(7), 0._DP, 0._DP, &
        2*e2**2*Q(4)*Q(5), 4*e2*Q(3)*Q(5)**2, 4*Q(1)*Q(2)*Q(6), 0._DP, 2*e2*Q(5)*Q(14), &
        2*e2*Q(4)*Q(5)**2, 4*e2*Q(1)*Q(2)*Q(5), 4*Q(1)*Q(2)*Q(5)**2/)
      CASE (62_I2B)
        Q(1:18)=(/e2+A(1,L)-e1*A(1,L), e2+A(2,L)-e1*A(2,L), e2-e5*(A(1,L) &
        + A(2,L)), 2*e2-e5*(A(1,L)+A(2,L)), e2+A(3,L)-e1*A(3,L), e2+A(4,L) &
        - e1*A(4,L), e2-e5*(A(3,L)+A(4,L)), 2*e2-e5*(A(3,L)+A(4,L)), e2-&
        e5*(A(1,L)+A(2,L)+A(3,L)+A(4,L)), 2*e2-e5*(A(1,L)+A(2,L)+A(3,L)+&
        A(4,L)), 2*e2-e5*(2*A(1,L)+2*A(2,L)+A(3,L)+A(4,L)), 2*e2**2-&
        2*e2*e5*(A(1,L)+A(2,L)+A(3,L)+A(4,L))+e5**2*((A(1,L)+A(2,L))*A(3,L) &
        + (A(1,L)+A(2,L)+2*A(3,L))*A(4,L)), 2*e2-e5*(A(1,L)+A(2,L)+&
        2*(A(3,L)+A(4,L))), 4*e2**2+e5**2*((A(1,L)+A(2,L))*A(3,L)+(A(1,L)+&
        A(2,L)+4*A(3,L))*A(4,L))-2*e2*e5*(A(1,L)+A(2,L)+2*(A(3,L)+A(4,L))), &
        2*e2**2-2*e2*e5*(A(1,L)+A(2,L)+A(3,L)+A(4,L))+e5**2*(A(2,L)*(A(3,L) &
        + A(4,L))+A(1,L)*(2*A(2,L)+A(3,L)+A(4,L))), 2*e2**2-2*e2*e5*(A(1,L)+&
        A(2,L)+A(3,L)+A(4,L))+e5**2*(2*A(3,L)*A(4,L)+A(2,L)*(A(3,L)+A(4,L)) &
        + A(1,L)*(2*A(2,L)+A(3,L)+A(4,L))), 4*e2**2-2*e2*e5*(2*A(1,L)+&
        2*A(2,L)+A(3,L)+A(4,L))+e5**2*(A(2,L)*(A(3,L)+A(4,L))+&
        A(1,L)*(4*A(2,L)+A(3,L)+A(4,L))), 4*e2**3-4*e2**2*e5*(A(1,L)+A(2,L)+&
        A(3,L)+A(4,L))+e2*e5**2*(4*A(3,L)*A(4,L)+3*A(2,L)*(A(3,L)+A(4,L))+&
        A(1,L)*(4*A(2,L)+3*(A(3,L)+A(4,L))))-2*e5**3*(A(2,L)*A(3,L)*A(4,L)+&
        A(1,L)*(A(3,L)*A(4,L)+A(2,L)*(A(3,L)+A(4,L))))/)
        coe(L,:)=(/4*e2**4, 4*e2**4*A(5,L), 2*e2**3*Q(8), 2*e2**3*Q(4), &
        4*e2**4*A(5,L), 4*e2**4*A(5,L), 2*e2**3*Q(10), 4*e2**3*Q(9), 4*e2**3*Q(7), &
        4*e2**3*Q(3), 2*e2**3*Q(11), 4*e2**4*A(5,L), 2*e2**3*Q(13), 4*e2**4*A(5,L), &
        4*e2**2*Q(3)*Q(7), 4*e2**3*A(5,L)*Q(3), 4*e2**3*A(5,L)*Q(7), 2*e2**2*Q(16), &
        2*e2**3*A(5,L)*Q(10), 4*e2**3*A(5,L)*Q(9), 4*e2**2*Q(5)*Q(6), &
        4*e2**2*Q(1)*Q(2), 2*e2**3*A(5,L)*Q(8), 2*e2**3*A(5,L)*Q(4), e2**2*Q(4)*Q(8), &
        2*e2**3*A(5,L)*Q(11), 2*e2**3*A(5,L)*Q(13), e2**2*Q(14), 2*e2**3*A(5,L)*Q(8), &
        2*e2**3*A(5,L)*Q(4), e2**2*Q(17), 2*e2**2*Q(12), 4*e2**3*A(5,L)*Q(7), &
        2*e2**2*Q(15), 2*e2**2*Q(4)*Q(7), 4*e2**3*A(5,L)*Q(3), 2*e2**2*Q(3)*Q(8), &
        2*e2**3*A(5,L)*Q(10), 2*e2**3*A(5,L)*Q(8), 2*e2**3*A(5,L)*Q(4), &
        4*e2**2*A(5,L)*Q(3)*Q(7), 4*e2*Q(3)*Q(5)*Q(6), 4*e2*Q(1)*Q(2)*Q(7), &
        4*e2**2*A(5,L)*Q(1)*Q(2), 4*e2**2*A(5,L)*Q(5)*Q(6), 2*e2**2*A(5,L)*Q(16), &
        2*e2**2*A(5,L)*Q(3)*Q(8), 2*e2**2*A(5,L)*Q(4)*Q(7), e2*Q(18), &
        e2**2*A(5,L)*Q(14), e2**2*A(5,L)*Q(17), e2**2*A(5,L)*Q(4)*Q(8), &
        2*e2**2*A(5,L)*Q(12), 4*e2**2*A(5,L)*Q(5)*Q(6), 4*e2**2*A(5,L)*Q(1)*Q(2), &
        2*e2**2*A(5,L)*Q(15), 2*e2*Q(4)*Q(5)*Q(6), 2*e2*Q(1)*Q(2)*Q(8), &
        e2**2*A(5,L)*Q(4)*Q(8), 4*e2*A(5,L)*Q(3)*Q(5)*Q(6), &
        4*e2*A(5,L)*Q(1)*Q(2)*Q(7), 4*Q(1)*Q(2)*Q(5)*Q(6), e2*A(5,L)*Q(18), &
        2*e2*A(5,L)*Q(4)*Q(5)*Q(6), 2*e2*A(5,L)*Q(1)*Q(2)*Q(8), &
        4*A(5,L)*Q(1)*Q(2)*Q(5)*Q(6)/)
      CASE (63_I2B)
        X(1:2)=Merge(A(1:2,L), A(2:1:-1,L), ANY(G(1,L)==G(3:4,L)))
        X(3)=Merge(A(4,L), A(3,L), ANY(G(3,L)==G(1:2,L)))
        Q(1:17)=(/e2+x(1)-e1*x(1), e2+x(2)-e1*x(2), e2-e5*(x(1)+x(2)), &
        2*e2-e5*(x(1)+x(2)), e2+x(3)-e1*x(3), e2-e5*(x(1)+x(3)), 2*e2-&
        e5*(x(1)+x(3)), 4*e2**2+e5**2*x(1)-2*e2*e5*(2*x(1)+x(2)+x(3)), &
        e2**2+e5**2*x(1)-e2*e5*(2*x(1)+x(2)+x(3)), 2*e2**2+e5**2*x(1)-&
        e2*e5*(3*x(1)+2*x(2)+x(3)), 2*e2**2+e5**2*x(1)-e2*e5*(3*x(1)+x(2)+&
        2*x(3)), 2*e2**3-e5**3*x(1)*x(2)-2*e2**2*e5*(2*x(1)+x(2)+x(3))+&
        e2*e5**2*(x(1)*(1+x(1)+3*x(2))+(x(1)+x(2))*x(3)), 4*e2**3-&
        e5**3*x(1)*x(2)-2*e2**2*e5*(3*x(1)+2*x(2)+x(3))+e2*e5**2*(x(1)*(1+&
        x(1)+5*x(2))+(x(1)+x(2))*x(3)), 2*e2**3-e5**3*x(1)*x(3)-&
        2*e2**2*e5*(2*x(1)+x(2)+x(3))+e2*e5**2*(x(1)*(1+x(1)+x(2))+&
        (3*x(1)+x(2))*x(3)), 2*e2**3-e5**3*x(1)*(x(2)+x(3))-&
        2*e2**2*e5*(2*x(1)+x(2)+x(3))+e2*e5**2*(x(1)*(1+x(1)+3*x(2))+&
        (3*x(1)+x(2))*x(3)), 4*e2**3-e5**3*x(1)*x(3)-2*e2**2*e5*(3*x(1)+x(2) &
        + 2*x(3))+e2*e5**2*(x(1)*(1+x(1)+x(2))+(5*x(1)+x(2))*x(3)), 4*e2**4 &
        + e5**4*x(1)*x(2)*x(3)-4*e2**3*e5*(2*x(1)+x(2)+x(3))-&
        e2*e5**3*x(1)*(x(3)+2*x(1)*x(3)+x(2)*(1+2*x(1)+4*x(3)))+&
        e2**2*e5**2*(3*x(1)**2+3*x(2)*x(3)+x(1)*(1+7*x(2)+7*x(3)))/)
        coe(L,:)=(/0._DP, 8*e2**4, 0._DP, 0._DP, 8*e2**4, 8*e2**4, &
        0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 8*e2**4, 0._DP, &
        8*e2**4, 0._DP, 8*e2**3*Q(3), 8*e2**3*Q(6), 0._DP, 2*e2**2*Q(8), 8*e2**2*Q(9), 0._DP, 0._DP, &
        4*e2**3*Q(7), 4*e2**3*Q(4), 0._DP, 4*e2**2*Q(10), 4*e2**2*Q(11), 0._DP, 4*e2**3*Q(7), &
        4*e2**3*Q(4), 0._DP, 0._DP, 8*e2**3*Q(6), 0._DP, 0._DP, 8*e2**3*Q(3), 0._DP, 2*e2**2*Q(8), &
        4*e2**3*Q(7), 4*e2**3*Q(4), 8*e2**2*Q(3)*Q(6), 0._DP, 0._DP, 8*e2**2*Q(1)*Q(2), &
        8*e2**2*Q(1)*Q(5), 4*e2*Q(15), 4*e2**2*Q(3)*Q(7), 4*e2**2*Q(4)*Q(6), 0._DP, &
        2*e2*Q(16), 2*e2*Q(13), 2*e2**2*Q(4)*Q(7), 4*e2*Q(14), 8*e2**2*Q(1)*Q(5), &
        8*e2**2*Q(1)*Q(2), 4*e2*Q(12), 0._DP, 0._DP, 2*e2**2*Q(4)*Q(7), 8*e2*Q(1)*Q(3)*Q(5), &
        8*e2*Q(1)*Q(2)*Q(6), 0._DP, 2*Q(17), 4*e2*Q(1)*Q(4)*Q(5), 4*e2*Q(1)*Q(2)*Q(7), &
        8*Q(1)**2*Q(2)*Q(5)/)
      CASE (64_I2B)
        X(1:2)=Merge(A(1:2,L), A(2:1:-1,L), ANY(G(1,L)==G(5:6,L)))
        Q(1:19)=(/1-e1+2*e2, e2+A(3,L)-e1*A(3,L), e2+A(4,L)-e1*A(4,L), &
        e2-e5*(A(3,L)+A(4,L)), 2*e2-e5*(A(3,L)+A(4,L)), e2+x(1)-e1*x(1), &
        e2+x(2)-e1*x(2), 2*(1-e1+4*e2)*(e2+A(3,L)-e1*A(3,L))*(e2+A(4,L) &
        - e1*A(4,L))-2*e2*e5*(2*e2-e5*(A(3,L)+A(4,L)))*x(1)-(1-e1+&
        2*e2)*e5*(2*e2-e5*(A(3,L)+A(4,L)))*x(2), e2-e5*(x(1)+x(2)), e2-&
        e5*(A(3,L)+A(4,L)+x(1)+x(2)), 4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(1) &
        + 2*x(2)), 8*e2**2+e5**2*(A(3,L)+A(4,L)+2*x(2))-2*e2*e5*(1+2*A(3,L) &
        + 2*A(4,L)+2*x(1)+2*x(2)), 4*e2**2+e5**2*(A(3,L)+A(4,L)+x(2))-&
        e2*e5*(1+4*A(3,L)+4*A(4,L)+2*x(1)+2*x(2)), 2*e2**2-2*e2*e5*(A(3,L) &
        + A(4,L)+x(1)+x(2))+e5**2*((A(3,L)+A(4,L))*x(1)+(A(3,L)+A(4,L)+&
        2*x(1))*x(2)), 2*e2-e5*(A(3,L)+A(4,L)+2*(x(1)+x(2))), 4*e2**2+&
        e5**2*((A(3,L)+A(4,L))*x(1)+(A(3,L)+A(4,L)+4*x(1))*x(2))-&
        2*e2*e5*(A(3,L)+A(4,L)+2*(x(1)+x(2))), 2*e2**2-2*e2*e5*(A(3,L)+&
        A(4,L)+x(1)+x(2))+e5**2*(A(4,L)*(x(1)+x(2))+A(3,L)*(2*A(4,L)+x(1) &
        + x(2))), 2*e2**2-2*e2*e5*(A(3,L)+A(4,L)+x(1)+x(2))+&
        e5**2*(2*x(1)*x(2)+A(4,L)*(x(1)+x(2))+A(3,L)*(2*A(4,L)+x(1)+x(2))), &
        4*e2**3-4*e2**2*e5*(A(3,L)+A(4,L)+x(1)+x(2))+e2*e5**2*(4*x(1)*x(2) &
        + 3*A(4,L)*(x(1)+x(2))+A(3,L)*(4*A(4,L)+3*(x(1)+x(2))))-&
        2*e5**3*(A(4,L)*x(1)*x(2)+A(3,L)*(x(1)*x(2)+A(4,L)*(x(1)+x(2))))/)
        coe(L,:)=(/0._DP, 4*e2**3*Q(1), 0._DP, 0._DP, 4*e2**3*Q(1), &
        4*e2**3*Q(1), 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, &
        4*e2**3*Q(1), 0._DP, 4*e2**3*Q(1), 0._DP, 8*e2**3*Q(9), 4*e2**2*Q(1)*Q(4), 0._DP, &
        e2**2*Q(12), 8*e2**3*Q(10), 0._DP, 0._DP, 2*e2**2*Q(1)*Q(5), 2*e2**2*Q(11), 0._DP, &
        4*e2**3*(2*e2-e5*(A(3,L)+A(4,L)+2*x(1)+2*x(2))), 2*e2**2*Q(13), 0._DP, &
        2*e2**2*Q(1)*Q(5), 2*e2**2*Q(11), 0._DP, 0._DP, 4*e2**2*Q(1)*Q(4), 0._DP, 0._DP, &
        8*e2**3*Q(9), 0._DP, e2**2*Q(12), 2*e2**2*Q(1)*Q(5), 2*e2**2*Q(11), &
        8*e2**2*Q(4)*Q(9), 0._DP, 0._DP, 8*e2**2*Q(6)*Q(7), 4*e2*Q(1)*Q(2)*Q(3), &
        4*e2**2*(2*e2**2-2*e2*e5*(A(3,L)+A(4,L)+x(1)+x(2))+&
        e5**2*(2*A(3,L)*A(4,L)+(A(3,L)+A(4,L))*x(2)+x(1)*(A(3,L)+A(4,L)+&
        2*x(2)))), 4*e2**2*Q(5)*Q(9), 2*e2*Q(4)*Q(11), 0._DP, e2*(8*e2**3-2*e2**2*e5*(1 &
        + 4*A(3,L)+4*A(4,L)+2*x(1)+2*x(2))-e5**3*(2*A(3,L)*A(4,L)+(A(3,L)+&
        A(4,L))*x(2))+2*e2*e5**2*(4*A(3,L)*A(4,L)+(A(3,L)+A(4,L))*(1+x(1))+&
        (1+A(3,L)+A(4,L))*x(2))), 2*e2**2*(4*e2**2-2*e2*e5*(A(3,L)+A(4,L)+&
        2*x(1)+2*x(2))+e5**2*((A(3,L)+A(4,L))*x(2)+x(1)*(A(3,L)+A(4,L)+&
        4*x(2)))), e2*Q(5)*Q(11), 4*e2**2*(2*e2**2-2*e2*e5*(A(3,L)+A(4,L)+x(1) &
        + x(2))+e5**2*(A(3,L)*(x(1)+x(2))+A(4,L)*(2*A(3,L)+x(1)+x(2)))), &
        4*e2*Q(1)*Q(2)*Q(3), 8*e2**2*Q(6)*Q(7), 4*e2**2*(2*e2**2-2*e2*e5*(A(3,L)+&
        A(4,L)+x(1)+x(2))+e5**2*((A(3,L)+A(4,L))*x(2)+x(1)*(A(3,L)+A(4,L) &
        + 2*x(2)))), 0._DP, 0._DP, e2*Q(5)*Q(11), 8*e2*Q(2)*Q(3)*Q(9), 8*e2*Q(4)*Q(6)*Q(7), &
        0._DP, 2*e2*(4*e2**3-4*e2**2*e5*(A(3,L)+A(4,L)+x(1)+x(2))+&
        e2*e5**2*(4*A(3,L)*A(4,L)+3*(A(3,L)+A(4,L))*x(2)+x(1)*(3*(A(3,L)+&
        A(4,L))+4*x(2)))-2*e5**3*(A(3,L)*A(4,L)*x(2)+x(1)*(A(3,L)*A(4,L)+&
        (A(3,L)+A(4,L))*x(2)))), 2*Q(2)*Q(3)*Q(11), 4*e2*Q(5)*Q(6)*Q(7), &
        8*Q(2)*Q(3)*Q(6)*Q(7)/)
      CASE (65_I2B)
        X(1:2)=Merge(A(3:4,L), A(4:3:-1,L), ANY(G(3,L)==G(5:6,L)))
        Q(1:19)=(/1-e1+2*e2, e2+A(1,L)-e1*A(1,L), e2+A(2,L)-e1*A(2,L), &
        e2-e5*(A(1,L)+A(2,L)), 2*e2-e5*(A(1,L)+A(2,L)), e2+x(1)-e1*x(1), &
        e2+x(2)-e1*x(2), 2*(1-e1+4*e2)*(e2+A(1,L)-e1*A(1,L))*(e2+A(2,L) &
        - e1*A(2,L))-2*e2*e5*(2*e2-e5*(A(1,L)+A(2,L)))*x(1)-(1-e1+&
        2*e2)*e5*(2*e2-e5*(A(1,L)+A(2,L)))*x(2), e2-e5*(x(1)+x(2)), e2-&
        e5*(A(1,L)+A(2,L)+x(1)+x(2)), 4*e2**2+e5**2*x(2)-e2*e5*(1+2*x(1) &
        + 2*x(2)), 8*e2**2+e5**2*(A(1,L)+A(2,L)+2*x(2))-2*e2*e5*(1+2*A(1,L) &
        + 2*A(2,L)+2*x(1)+2*x(2)), 4*e2**2+e5**2*(A(1,L)+A(2,L)+x(2))-&
        e2*e5*(1+4*A(1,L)+4*A(2,L)+2*x(1)+2*x(2)), 2*e2**2-2*e2*e5*(A(1,L) &
        + A(2,L)+x(1)+x(2))+e5**2*((A(1,L)+A(2,L))*x(1)+(A(1,L)+A(2,L)+&
        2*x(1))*x(2)), 2*e2-e5*(A(1,L)+A(2,L)+2*(x(1)+x(2))), 4*e2**2+&
        e5**2*((A(1,L)+A(2,L))*x(1)+(A(1,L)+A(2,L)+4*x(1))*x(2))-&
        2*e2*e5*(A(1,L)+A(2,L)+2*(x(1)+x(2))), 2*e2**2-2*e2*e5*(A(1,L)+&
        A(2,L)+x(1)+x(2))+e5**2*(A(2,L)*(x(1)+x(2))+A(1,L)*(2*A(2,L)+x(1) &
        + x(2))), 2*e2**2-2*e2*e5*(A(1,L)+A(2,L)+x(1)+x(2))+&
        e5**2*(2*x(1)*x(2)+A(2,L)*(x(1)+x(2))+A(1,L)*(2*A(2,L)+x(1)+x(2))), &
        4*e2**3-4*e2**2*e5*(A(1,L)+A(2,L)+x(1)+x(2))+e2*e5**2*(4*x(1)*x(2) &
        + 3*A(2,L)*(x(1)+x(2))+A(1,L)*(4*A(2,L)+3*(x(1)+x(2))))-&
        2*e5**3*(A(2,L)*x(1)*x(2)+A(1,L)*(x(1)*x(2)+A(2,L)*(x(1)+x(2))))/)
        coe(L,:)=(/0._DP, 4*e2**3*Q(1), 0._DP, 0._DP, 4*e2**3*Q(1), &
        4*e2**3*Q(1), 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, &
        4*e2**3*Q(1), 0._DP, 4*e2**3*Q(1), 0._DP, 4*e2**2*Q(1)*Q(4), 8*e2**3*Q(9), 0._DP, &
        e2**2*Q(12), 8*e2**3*Q(10), 0._DP, 0._DP, 2*e2**2*Q(11), 2*e2**2*Q(1)*Q(5), 0._DP, &
        2*e2**2*Q(13), 4*e2**3*Q(15), 0._DP, 2*e2**2*Q(11), 2*e2**2*Q(1)*Q(5), 0._DP, 0._DP, &
        8*e2**3*Q(9), 0._DP, 0._DP, 4*e2**2*Q(1)*Q(4), 0._DP, e2**2*Q(12), 2*e2**2*Q(11), &
        2*e2**2*Q(1)*Q(5), 8*e2**2*Q(4)*Q(9), 0._DP, 0._DP, 4*e2*Q(1)*Q(2)*Q(3), &
        8*e2**2*Q(6)*Q(7), 4*e2**2*Q(18), 2*e2*Q(4)*Q(11), 4*e2**2*Q(5)*Q(9), 0._DP, &
        2*e2**2*Q(16), e2*Q(8), e2*Q(5)*Q(11), 4*e2**2*Q(14), 8*e2**2*Q(6)*Q(7), &
        4*e2*Q(1)*Q(2)*Q(3), 4*e2**2*Q(17), 0._DP, 0._DP, e2*Q(5)*Q(11), 8*e2*Q(4)*Q(6)*Q(7), &
        8*e2*Q(2)*Q(3)*Q(9), 0._DP, 2*e2*Q(19), 4*e2*Q(5)*Q(6)*Q(7), 2*Q(2)*Q(3)*Q(11), &
        8*Q(2)*Q(3)*Q(6)*Q(7)/)
      CASE (66_I2B)
        Q(1:18)=(/e2+A(1,L)-e1*A(1,L), e2+A(2,L)-e1*A(2,L), e2-e5*(A(1,L) &
        + A(2,L)), 2*e2-e5*(A(1,L)+A(2,L)), e2+A(3,L)-e1*A(3,L), e2+A(4,L) &
        - e1*A(4,L), e2-e5*(A(3,L)+A(4,L)), 2*e2-e5*(A(3,L)+A(4,L)), e2-&
        e5*(A(1,L)+A(2,L)+A(3,L)+A(4,L)), 2*e2-e5*(A(1,L)+A(2,L)+A(3,L)+&
        A(4,L)), 2*e2-e5*(2*A(1,L)+2*A(2,L)+A(3,L)+A(4,L)), 2*e2**2-&
        2*e2*e5*(A(1,L)+A(2,L)+A(3,L)+A(4,L))+e5**2*((A(1,L)+A(2,L))*A(3,L) &
        + (A(1,L)+A(2,L)+2*A(3,L))*A(4,L)), 2*e2-e5*(A(1,L)+A(2,L)+&
        2*(A(3,L)+A(4,L))), 4*e2**2+e5**2*((A(1,L)+A(2,L))*A(3,L)+(A(1,L)+&
        A(2,L)+4*A(3,L))*A(4,L))-2*e2*e5*(A(1,L)+A(2,L)+2*(A(3,L)+A(4,L))), &
        2*e2**2-2*e2*e5*(A(1,L)+A(2,L)+A(3,L)+A(4,L))+e5**2*(A(2,L)*(A(3,L) &
        + A(4,L))+A(1,L)*(2*A(2,L)+A(3,L)+A(4,L))), 2*e2**2-2*e2*e5*(A(1,L)+&
        A(2,L)+A(3,L)+A(4,L))+e5**2*(2*A(3,L)*A(4,L)+A(2,L)*(A(3,L)+A(4,L)) &
        + A(1,L)*(2*A(2,L)+A(3,L)+A(4,L))), 4*e2**2-2*e2*e5*(2*A(1,L)+&
        2*A(2,L)+A(3,L)+A(4,L))+e5**2*(A(2,L)*(A(3,L)+A(4,L))+&
        A(1,L)*(4*A(2,L)+A(3,L)+A(4,L))), 4*e2**3-4*e2**2*e5*(A(1,L)+A(2,L)+&
        A(3,L)+A(4,L))+e2*e5**2*(4*A(3,L)*A(4,L)+3*A(2,L)*(A(3,L)+A(4,L))+&
        A(1,L)*(4*A(2,L)+3*(A(3,L)+A(4,L))))-2*e5**3*(A(2,L)*A(3,L)*A(4,L)+&
        A(1,L)*(A(3,L)*A(4,L)+A(2,L)*(A(3,L)+A(4,L))))/)
        coe(L,:)=(/0._DP, 8*e2**4, 0._DP, 0._DP, 8*e2**4, 8*e2**4, &
        0._DP, 0._DP, 0._DP, 0._DP, 0._DP, 8*e2**4, 0._DP, &
        8*e2**4, 0._DP, 8*e2**3*Q(3), 8*e2**3*Q(7), 0._DP, 4*e2**3*Q(10), 8*e2**3*Q(9), 0._DP, 0._DP, &
        4*e2**3*Q(8), 4*e2**3*Q(4), 0._DP, 4*e2**3*Q(11), 4*e2**3*Q(13), 0._DP, 4*e2**3*Q(8), &
        4*e2**3*Q(4), 0._DP, 0._DP, 8*e2**3*Q(7), 0._DP, 0._DP, 8*e2**3*Q(3), 0._DP, 4*e2**3*Q(10), &
        4*e2**3*Q(8), 4*e2**3*Q(4), 8*e2**2*Q(3)*Q(7), 0._DP, 0._DP, 8*e2**2*Q(1)*Q(2), &
        8*e2**2*Q(5)*Q(6), 4*e2**2*Q(16), 4*e2**2*Q(3)*Q(8), 4*e2**2*Q(4)*Q(7), 0._DP, &
        2*e2**2*Q(14), 2*e2**2*Q(17), 2*e2**2*Q(4)*Q(8), 4*e2**2*Q(12), &
        8*e2**2*Q(5)*Q(6), 8*e2**2*Q(1)*Q(2), 4*e2**2*Q(15), 0._DP, 0._DP, 2*e2**2*Q(4)*Q(8), &
        8*e2*Q(3)*Q(5)*Q(6), 8*e2*Q(1)*Q(2)*Q(7), 0._DP, 2*e2*Q(18), 4*e2*Q(4)*Q(5)*Q(6), &
        4*e2*Q(1)*Q(2)*Q(8), 8*Q(1)*Q(2)*Q(5)*Q(6)/)
    END SELECT
  END DO
END Subroutine Coefficient
!
SUBROUTINE Coefficient2
  USE VARIABLES
  USE IBDVariables
  IMPLICIT NONE
  INTEGER:: LL, L
  REAL(DP):: Y, P(4), e2, e0, Q(3)
  DO L = 1, NumLoci
    IF(GConfig(L) == 0) THEN
      Coe(L, 1 : 9) = 0._DP
      CYCLE
    END IF
    LL = LocusSelectPtr(L)
    IF(ANY(G(1 : 4, L) < 1_I1B))CYCLE
    e2 = Err(LL) / ObsNumAllele(LL)
    e0 = -1._DP + Err(LL)
    P(1:4)=ObsAlleleFre(G(1 : 4, L), LL)
    SELECT CASE (GConfig(L))
    CASE (1_I2B)
      Y = (e2 - e0*p(1))
      Coe(L, 1) = e2 ** 4 + e0*(e0 - 2*e2)*(e0 ** 2 - 2*e0*e2 + 2*e2 ** 2)* p(1)
      Coe(L, 2) = (e2 ** 2 + e0 ** 2*p(1) - 2*e0*e2*p(1)) ** 2
      Coe(L, 3) = Y*(e2 ** 3 - e0*(e0 ** 2 - 3*e0*e2 + 3*e2 ** 2)*p(1))
      Coe(L, 4) = Y ** 2*(e2 ** 2 + e0 ** 2*p(1) - 2*e0*e2*p(1))
      Coe(L, 5) = Coe(L, 3)
      Coe(L, 6) = Coe(L, 4)
      Coe(L, 7) = Coe(L, 2)
      Coe(L, 8) = Coe(L, 4)
      Coe(L, 9) = Y ** 4
    CASE (2_I2B)
      Coe(L, 1) = e2 ** 2*(e2 ** 2 + e0 ** 2*(p(1) + p(3)) - 2*e0*e2*(p(1) + p(3)))
      Coe(L, 2) = (e2 ** 2 + e0 ** 2*p(1) - 2*e0*e2*p(1))*(e2 ** 2 + &
        e0 ** 2*p(3) - 2*e0*e2*p(3))
      Coe(L, 3) = e2*(e2 - e0*p(3))*(e2 ** 2 + e0 ** 2*p(1) - e0*e2*(2*p(1) + p(3)))
      Coe(L, 4) = (e2 ** 2 + e0 ** 2*p(1) - 2*e0*e2*p(1))*(e2 - e0*p(3)) ** 2
      Coe(L, 5) = e2*(e2 - e0*p(1))*(e2 ** 2 + e0 ** 2*p(3) - e0*e2*(p(1) + 2*p(3)))
      Coe(L, 6) = (e2 - e0*p(1)) ** 2*(e2 ** 2 + e0 ** 2*p(3) - 2*e0*e2*p(3))
      Coe(L, 7) = e2 ** 2*(e2 - e0*(p(1) + p(3))) ** 2
      Coe(L, 8) = e2*(e2 - e0*p(1))*(e2 - e0*p(3))*(e2 - e0*(p(1) + p(3)))
      Coe(L, 9) = (e2 - e0*p(1)) ** 2*(e2 - e0*p(3)) ** 2
    CASE (3_I2B)
      Y=Merge(p(4), p(3), G(1, LL)==G(3, LL))
      Coe(L, 1) = 2*e2 ** 3*(e2 - e0*y) - &
        2*e0*e2*(e0 ** 2 - 3*e0*e2 + 3*e2 ** 2) * p(1)
      Coe(L, 2) = 2*e2*(e2 ** 2 + e0 ** 2*p(1) - 2*e0*e2*p(1))*(e2 - e0*(y + p(1)))
      Coe(L, 3) = 2*e2 ** 4 + e0 ** 4*y*p(1) + 2*e0 ** 2*e2 ** 2*p(1)*(2 + 2*y + p(1)) - &
        e0 ** 3*e2*p(1)*(1 + 3*y + p(1)) - 2*e0*e2 ** 3*(y + 3*p(1))
      Coe(L, 4) = 2*(e2 - e0*y)*(e2 - e0*p(1))*(e2 ** 2 + &
        e0 ** 2*p(1) - 2*e0*e2*p(1))
      Coe(L, 5) = 2*e2*(e2 - e0*p(1))*(e2 ** 2 + &
        e0 ** 2*p(1) - e0*e2*(y + 2*p(1)))
      Coe(L, 6) = 2*e2*(e2 - e0*p(1)) ** 2*(e2 - e0*(y + p(1)))
      Coe(L, 7) = 2*e2*(e2 ** 2 + e0 ** 2*p(1) - 2*e0*e2*p(1))*(e2 - e0*(y + p(1)))
      Coe(L, 8) = (e2 - e0*p(1))*(2*e2 ** 3 - &
        e0 ** 3*y*p(1) + e0 ** 2*e2*p(1)*(1 + 3*y + p(1)) - 2*e0*e2 ** 2*(y + 2*p(1)))
      Coe(L, 9) = 2*(e2 - e0*y)*(e2 - e0*p(1)) ** 3
    CASE (4_I2B)
      Coe(L, 1) = 2*e2 ** 2*(e2 ** 2 + e0 ** 2*p(1) - e0*e2*(2*p(1) + p(3) + p(4)))
      Coe(L, 2) = 2*e2*(e2 ** 2 + e0 ** 2*p(1) - &
        2*e0*e2*p(1))*(e2 - e0*(p(3) + p(4)))
      Coe(L, 3) = e2*(2*e2 ** 3 - e0 ** 3*p(1)*(p(3) + p(4)) - &
        2*e0*e2 ** 2*(2*p(1) + p(3) + p(4)) + &
        2*e0 ** 2*e2*(p(3)*p(4) + p(1)*(1 + p(3) + p(4))))
      Coe(L, 4) = 2*(e2 ** 2 + e0 ** 2*p(1) - 2*e0*e2*p(1))*(e2 - e0*p(3))*(e2 - e0*p(4))
      Coe(L, 5) = 2*e2 ** 2*(e2 - e0*p(1))*(e2 - e0*(p(1) + p(3) + p(4)))
      Coe(L, 6) = 2*e2*(e2 - e0*p(1)) ** 2*(e2 - e0*(p(3) + p(4)))
      Coe(L, 7) = 2*e2 ** 2*(e2 - e0*(p(1) + p(3)))*(e2 - e0*(p(1) + p(4)))
      Coe(L, 8) = e2*(e2 - e0*p(1))*(2*e2 ** 2 - 2*e0*e2*(p(1) + p(3) + p(4)) + &
        e0 ** 2*(2*p(3)*p(4) + p(1)*(p(3) + p(4))))
      Coe(L, 9) = 2*(e2 - e0*p(1)) ** 2*(e2 - e0*p(3))*(e2 - e0*p(4))
    CASE (5_I2B)
      Y=Merge(p(2), p(1), G(1, LL)==G(3, LL))
      Coe(L, 1) = 2*e2 ** 3*(e2 - e0*y) - &
        2*e0*e2*(e0 ** 2 - 3*e0*e2 + 3*e2 ** 2)*p(3)
      Coe(L, 2) = 2*e2*(e2 ** 2 + e0 ** 2*p(3) - 2*e0*e2*p(3))*(e2 - e0*(y + p(3)))
      Coe(L, 3) = 2*e2*(e2 - e0*p(3))*(e2 ** 2 + e0 ** 2*p(3) - &
        e0*e2*(y + 2*p(3)))
      Coe(L, 4) = 2*e2*(e2 - e0*p(3)) ** 2*(e2 - e0*(y + p(3)))
      Coe(L, 5) = 2*e2 ** 4 + e0 ** 4*y*p(3) + 2*e0 ** 2*e2 ** 2*p(3)*(2 + &
        2*y + p(3)) - e0 ** 3*e2*p(3)*(1 + 3*y + p(3)) - &
        2*e0*e2 ** 3*(y + 3*p(3))
      Coe(L, 6) = 2*(e2 - e0*y)*(e2 - e0*p(3))*(e2 ** 2 + &
        e0 ** 2*p(3) - 2*e0*e2*p(3))
      Coe(L, 7) = 2*e2*(e2 ** 2 + e0 ** 2*p(3) - &
        2*e0*e2*p(3))*(e2 - e0*(y + p(3)))
      Coe(L, 8) = (e2 - e0*p(3))*(2*e2 ** 3 - &
        e0 ** 3*y*p(3) + e0 ** 2*e2*p(3)*(1 + 3*y + p(3)) - &
        2*e0*e2 ** 2*(y + 2*p(3)))
      Coe(L, 9) = 2*(e2 - e0*y)*(e2 - e0*p(3)) ** 3
    CASE (6_I2B)
      Coe(L, 1) = 2*e2 ** 2*(e2 ** 2 + e0 ** 2*p(3) - e0*e2*(p(1) + p(2) + 2*p(3)))
      Coe(L, 2) = 2*e2*(e2 - e0*(p(1) + p(2)))*(e2 ** 2 + e0 ** 2*p(3) - 2*e0*e2*p(3))
      Coe(L, 3) = 2*e2 ** 2*(e2 - e0*p(3))*(e2 - e0*(p(1) + p(2) + p(3)))
      Coe(L, 4) = 2*e2*(e2 - e0*(p(1) + p(2)))*(e2 - e0*p(3)) ** 2
      Coe(L, 5) = e2*(2*e2*(e2 - e0*p(1))*(e2 - e0*p(2)) - &
        e0*(e0 - 2*e2)*(-2*e2 + e0*(p(1) + p(2)))*p(3))
      Coe(L, 6) = 2*(e2 - e0*p(1))*(e2 - &
        e0*p(2))*(e2 ** 2 + e0 ** 2*p(3) - 2*e0*e2*p(3))
      Coe(L, 7) = 2*e2 ** 2*(e2 - e0*(p(1) + p(3)))*(e2 - e0*(p(2) + p(3)))
      Coe(L, 8) = e2*(e2 - e0*p(3))*(2*(e2 - e0*p(1))*(e2 - e0*p(2)) + &
        e0*(-2*e2 + e0*(p(1) + p(2)))*p(3))
      Coe(L, 9) = 2*(e2 - e0*p(1))*(e2 - e0*p(2))*(e2 - e0*p(3)) ** 2
   CASE (7_I2B)
     Coe(L, 1) = 4*e2 ** 2*(e2 ** 2 + e0 ** 2*(p(1) + p(2)) - 2*e0*e2*(p(1) + p(2)))
     Coe(L, 2) = 4*e2 ** 2*(e2 - e0*(p(1) + p(2))) ** 2
     Coe(L, 3) = 2*e2*(2*e2 ** 3 - 2*e0 ** 3*p(1)*p(2) - &
       4*e0*e2 ** 2*(p(1) + p(2)) + e0 ** 2*e2*(p(1) + p(1) ** 2 + &
       p(2) + 4*p(1)*p(2) + p(2) ** 2))
     Coe(L, 4) = 4*e2*(e2 - e0*p(1))*(e2 - e0*p(2))*(e2 - e0*(p(1) + p(2)))
     Coe(L, 5) = 2*e2*(2*e2 ** 3 - 2*e0 ** 3*p(1)*p(2) - &
       4*e0*e2 ** 2*(p(1) + p(2)) + e0 ** 2*e2*(p(1) + p(1) ** 2 + &
       p(2) + 4*p(1)*p(2) + p(2) ** 2))
     Coe(L, 6) = Coe(L, 4)
     Coe(L, 7) = 2*(2*e2 ** 4 + e0 ** 4*p(1)*p(2) - &
       4*e0 ** 3*e2*p(1)*p(2) - 4*e0*e2 ** 3*(p(1) + p(2)) + &
       e0 ** 2*e2 ** 2*(p(1) + p(1) ** 2 + p(2) + 6*p(1)*p(2) + p(2) ** 2))
     Coe(L, 8) = 4*e2 ** 4 - 8*e0*e2 ** 3*(p(1) + p(2)) + &
       e0 ** 4*p(1)*p(2)*(p(1) + p(2)) - &
       4*e0 ** 3*e2*p(1)*p(2)*(1 + p(1) + p(2)) + &
       e0 ** 2*e2 ** 2*(p(1) + 3*p(1) ** 2 + p(2) + 14*p(1)*p(2) + 3*p(2) ** 2)
     Coe(L, 9) = 4*(e2 - e0*p(1)) ** 2*(e2 - e0*p(2)) ** 2
    CASE (8_I2B)
      IF(G(1, LL)==G(3, LL))THEN
        Q(1:2)=P(1:2)
        Q(3)=p(4)
      ELSE IF(G(1, LL)==G(4, LL))THEN
        Q(1 : 3)=p(1 : 3)
      ELSE IF(G(2, LL)==G(3, LL))THEN
        Q(1)=p(2)
        Q(2)=p(1)
        Q(3)=p(4)
      ELSE IF(G(2, LL)==G(4, LL))THEN
        Q(1)=p(2)
        Q(2)=p(1)
        Q(3)=p(3)
      END IF
      Coe(L, 1) = 4*e2 ** 2*(e2 ** 2 + e0 ** 2*q(1) - e0*e2*(2*q(1) + q(2) + q(3)))
      Coe(L, 2) = 4*e2 ** 2*(e2 - e0*(q(1) + q(2)))*(e2 - e0*(q(1) + q(3)))
      Coe(L, 3) = 2*e2*(2*e2 ** 3 - e0 ** 3*q(1)*q(3) - 2*e0*e2 ** 2*(2*q(1) + &
        q(2) + q(3)) + e0 ** 2*e2*(q(1)*(1 + q(1) + q(2)) + (3*q(1) + q(2))*q(3)))
      Coe(L, 4) = 4*e2*(e2 - e0*q(1))*(e2 - e0*(q(1) + q(2)))*(e2 - e0*q(3))
      Coe(L, 5) = 2*e2*(2*e2 ** 3 - e0 ** 3*q(1)*q(2) - 2*e0*e2 ** 2*(2*q(1) + &
        q(2) + q(3)) + e0 ** 2*e2*(q(1)*(1 + q(1) + 3*q(2)) + (q(1) + q(2))*q(3)))
      Coe(L, 6) = 4*e2*(e2 - e0*q(1))*(e2 - e0*q(2))*(e2 - e0*(q(1) + q(3)))
      Coe(L, 7) = 2*e2*(2*e2 ** 3 - e0 ** 3*q(1)*(q(2) + q(3)) - &
        2*e0*e2 ** 2*(2*q(1) + q(2) + q(3)) + &
        e0 ** 2*e2*(q(1)*(1 + q(1) + 3*q(2)) + (3*q(1) + q(2))*q(3)))
      Coe(L, 8) = 4*e2 ** 4 + e0 ** 4*q(1)*q(2)*q(3) - 4*e0*e2 ** 3*(2*q(1) + &
        q(2) + q(3)) - e0 ** 3*e2*q(1)*(q(3) + 2*q(1)*q(3) + q(2)*(1 + &
        2*q(1) + 4*q(3))) + e0 ** 2*e2 ** 2*(3*q(1) ** 2 + &
        3*q(2)*q(3) + q(1)*(1 + 7*q(2) + 7*q(3)))
      Coe(L, 9) = 4*(e2 - e0*q(1)) ** 2*(e2 - e0*q(2))*(e2 - e0*q(3))
    CASE (9)
      Coe(L, 1) = 4*e2 ** 3*(e2 - e0*(p(1) + p(2) + p(3) + p(4)))
      Coe(L, 2) = 4*e2 ** 2*(e2 - e0*(p(1) + p(2)))*(e2 - e0*(p(3) + p(4)))
      Coe(L, 3) = 2*e2 ** 2*(2*e2 ** 2 - 2*e0*e2*(p(1) + p(2) + p(3) + p(4)) + &
        e0 ** 2*((p(1) + p(2))*p(3) + (p(1) + p(2) + 2*p(3))*p(4)))
      Coe(L, 4) = 4*e2*(e2 - e0*(p(1) + p(2)))*(e2 - e0*p(3))*(e2 - e0*p(4))
      Coe(L, 5) = 2*e2 ** 2*(2*e2 ** 2 - 2*e0*e2*(p(1) + p(2) + p(3) + p(4)) + &
        e0 ** 2*(p(2)*(p(3) + p(4)) + p(1)*(2*p(2) + p(3) + p(4))))
      Coe(L, 6) = 4*e2*(e2 - e0*p(1))*(e2 - e0*p(2))*(e2 - e0*(p(3) + p(4)))
      Coe(L, 7) = 2*e2 ** 2*(2*e2 ** 2 - 2*e0*e2*(p(1) + p(2) + p(3) + p(4)) + &
        e0 ** 2*(2*p(3)*p(4) + p(2)*(p(3) + p(4)) + &
        p(1)*(2*p(2) + p(3) + p(4))))
      Coe(L, 8) = e2*(4*e2 ** 3 - 4*e0*e2 ** 2*(p(1) + p(2) + p(3) + p(4)) + &
        e0 ** 2*e2*(4*p(3)*p(4) + 3*p(2)*(p(3) + p(4)) + &
        p(1)*(4*p(2) + 3*(p(3) + p(4)))) - 2*e0 ** 3*(p(2)*p(3)*p(4) + &
        p(1)*(p(3)*p(4) + p(2)*(p(3) + p(4)))))
      Coe(L, 9) = 4*(e2 - e0*p(1))*(e2 - e0*p(2))*(e2 - e0*p(3))*(e2 - e0*p(4))
    END SELECT
  END DO
END SUBROUTINE Coefficient2
!
FUNCTION Ran1()
  Use Variables, Only : IDUM
  IMPLICIT NONE
  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
  REAL :: Ran1
  INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
  REAL, SAVE :: am
  INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
  if (idum <= 0 .or. iy < 0) then
      am=nearest(1.0,-1.0)/Real(IM)
      iy=ior(ieor(888889999,abs(idum)),1)
      ix=ieor(777755555,abs(idum))
      idum=abs(idum)+1
  end if
  ix=ieor(ix,ishft(ix,13))
  ix=ieor(ix,ishft(ix,-17))
  ix=ieor(ix,ishft(ix,5))
  k=iy/IQ
  iy=IA*(iy-k*IQ)-IR*k
  if (iy < 0) iy=iy+IM
  Ran1=am*ior(iand(IM,ieor(ix,iy)),1)
END FUNCTION Ran1
