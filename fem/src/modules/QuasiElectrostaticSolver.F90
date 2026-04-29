!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  This solver implements a scalar electroquasistatic potential equation,
! *
! *    div(sigma grad(Potential) + eps grad(d Potential / dt)) = -S,
! *
! *  using the weak form
! *
! *    (grad w, eps grad(d Potential / dt))
! *  + (grad w, sigma grad(Potential))
! *  = (w, S) + <w, Electric Flux>.
! *
! *  It is intended for transient weighting-potential calculations in
! *  conductive dielectrics.
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Initialize the scalar electroquasistatic potential solver.
!------------------------------------------------------------------------------
SUBROUTINE QuasiElectrostaticSolver_Init(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
!------------------------------------------------------------------------------
  Params => GetSolverParams()

  CALL ListAddNewString(Params, 'Variable', 'Potential')
  CALL ListAddNewInteger(Params, 'Time Derivative Order', 1)
!------------------------------------------------------------------------------
END SUBROUTINE QuasiElectrostaticSolver_Init
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solve the scalar electroquasistatic equation for electric potential.
!------------------------------------------------------------------------------
SUBROUTINE QuasiElectrostaticSolver(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params
  REAL(KIND=dp) :: Norm, RelativeChange
  INTEGER :: iter, NonlinearIter
  LOGICAL :: Found
  CHARACTER(*), PARAMETER :: Caller = 'QuasiElectrostaticSolver'
!------------------------------------------------------------------------------
  Params => GetSolverParams()

  NonlinearIter = GetInteger(Params, 'Nonlinear System Max Iterations', Found)
  IF (.NOT. Found) NonlinearIter = 1

  DO iter = 1, NonlinearIter
    CALL DefaultInitialize()
    CALL BulkAssembly()
    CALL DefaultFinishBulkAssembly()
    CALL BoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    Norm = DefaultSolve()
    RelativeChange = Solver % Variable % NonlinChange

    WRITE(Message, *) 'Result Norm      : ', Norm
    CALL Info(Caller, Message, Level=4)
    WRITE(Message, *) 'Relative Change  : ', RelativeChange
    CALL Info(Caller, Message, Level=4)

    IF (Solver % Variable % NonlinConverged == 1) EXIT
  END DO

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
    INTEGER :: tLocal, nLocal, ndLocal
!------------------------------------------------------------------------------
!$omp parallel do private(Element,nLocal,ndLocal)
    DO tLocal = 1, GetNOFActive()
      Element => GetActiveElement(tLocal)
      nLocal = GetElementNOFNodes(Element)
      ndLocal = GetElementNOFDOFs(Element)
      CALL LocalMatrix(Element, nLocal, ndLocal)
    END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE BoundaryAssembly()
!------------------------------------------------------------------------------
    INTEGER :: tLocal, nLocal, ndLocal
!------------------------------------------------------------------------------
!$omp parallel do private(Element,nLocal,ndLocal)
    DO tLocal = 1, GetNOFBoundaryElements()
      Element => GetBoundaryElement(tLocal)
      IF (.NOT. ActiveBoundaryElement(Element)) CYCLE
      nLocal = GetElementNOFNodes(Element)
      ndLocal = GetElementNOFDOFs(Element)
      CALL BoundaryLocalMatrix(Element, nLocal, ndLocal)
    END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE BoundaryAssembly
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Element, n, nd)
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Element
    INTEGER :: n, nd
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: Mass(nd, nd), Stiff(nd, nd), Force(nd)
    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd, 3), DetJ, Weight
    REAL(KIND=dp) :: Sigma(3, 3), Eps(3, 3), Source
    REAL(KIND=dp) :: SigmaN(3, 3, n), EpsN(3, 3, n), SourceN(n)
    LOGICAL :: IsScalar, Stat, FoundSource
    INTEGER :: i, j, p, q, t
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: Material, BodyForce
    TYPE(Nodes_t), SAVE :: Nodes
!$omp threadprivate(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes, Element)
    Material => GetMaterial(Element)
    BodyForce => GetBodyForce(Element)

    CALL GetConductivity(Material, SigmaN, IsScalar, Element)
    CALL GetPermittivity(Material, EpsN, IsScalar, Element)

    SourceN = 0.0_dp
    IF (ASSOCIATED(BodyForce)) THEN
      SourceN = GetReal(BodyForce, 'Current Source', FoundSource, Element)
      IF (.NOT. FoundSource) SourceN = 0.0_dp
    END IF

    Mass = 0.0_dp
    Stiff = 0.0_dp
    Force = 0.0_dp

    IP = GaussPoints(Element)
    DO t = 1, IP % n
      Stat = ElementInfo(Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
                         DetJ, Basis, dBasisdx)
      IF (.NOT. Stat) CYCLE

      Sigma = 0.0_dp
      Eps = 0.0_dp
      DO i = 1, 3
        DO j = 1, 3
          Sigma(i, j) = SUM(SigmaN(i, j, 1:n) * Basis(1:n))
          Eps(i, j) = SUM(EpsN(i, j, 1:n) * Basis(1:n))
        END DO
      END DO

      Source = SUM(SourceN(1:n) * Basis(1:n))
      Weight = IP % s(t) * DetJ

      DO p = 1, nd
        DO q = 1, nd
          Mass(p, q) = Mass(p, q) + Weight * &
              SUM(MATMUL(Eps, dBasisdx(q, :)) * dBasisdx(p, :))
          Stiff(p, q) = Stiff(p, q) + Weight * &
              SUM(MATMUL(Sigma, dBasisdx(q, :)) * dBasisdx(p, :))
        END DO
        Force(p) = Force(p) + Weight * Source * Basis(p)
      END DO
    END DO

    IF (Transient) THEN
      CALL Default1stOrderTime(Mass, Stiff, Force, UElement=Element)
    END IF
    CALL DefaultUpdateEquations(Stiff, Force, UElement=Element)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE BoundaryLocalMatrix(Element, n, nd)
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Element
    INTEGER :: n, nd
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: Mass(nd, nd), Stiff(nd, nd), Force(nd)
    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd, 3), DetJ, Weight
    REAL(KIND=dp) :: Flux, FluxN(n)
    LOGICAL :: FoundFlux, Stat
    INTEGER :: p, t
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Nodes_t), SAVE :: Nodes
!$omp threadprivate(Nodes)
!------------------------------------------------------------------------------
    BC => GetBC(Element)
    IF (.NOT. ASSOCIATED(BC)) RETURN

    FluxN = GetReal(BC, 'Electric Flux', FoundFlux, Element)
    IF (.NOT. FoundFlux) RETURN

    CALL GetElementNodes(Nodes, Element)

    Mass = 0.0_dp
    Stiff = 0.0_dp
    Force = 0.0_dp

    IP = GaussPoints(Element)
    DO t = 1, IP % n
      Stat = ElementInfo(Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
                         DetJ, Basis, dBasisdx)
      IF (.NOT. Stat) CYCLE

      Flux = SUM(FluxN(1:n) * Basis(1:n))
      Weight = IP % s(t) * DetJ

      DO p = 1, nd
        Force(p) = Force(p) + Weight * Flux * Basis(p)
      END DO
    END DO

    IF (Transient) THEN
      CALL Default1stOrderTime(Mass, Stiff, Force, UElement=Element)
    END IF
    CALL DefaultUpdateEquations(Stiff, Force, UElement=Element)
!------------------------------------------------------------------------------
  END SUBROUTINE BoundaryLocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE GetConductivity(Material, Tensor, IsScalar, Element)
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: Tensor(:, :, :)
    LOGICAL :: IsScalar
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    LOGICAL :: Found
!------------------------------------------------------------------------------
    CALL ReadTensor(Material, Tensor, IsScalar, 'Electric Conductivity', &
                    Element, Found)
    IF (.NOT. Found) THEN
      CALL ReadTensor(Material, Tensor, IsScalar, 'Electrical Conductivity', &
                      Element, Found)
    END IF
    IF (.NOT. Found) Tensor = 0.0_dp
!------------------------------------------------------------------------------
  END SUBROUTINE GetConductivity
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE GetPermittivity(Material, Tensor, IsScalar, Element)
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: Tensor(:, :, :)
    LOGICAL :: IsScalar
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Eps0
    LOGICAL :: Found, FoundEps0
!------------------------------------------------------------------------------
    CALL ReadTensor(Material, Tensor, IsScalar, 'Relative Permittivity', &
                    Element, Found)
    IF (Found) THEN
      Eps0 = ListGetConstReal(Model % Constants, 'Permittivity Of Vacuum', &
                              FoundEps0)
      IF (.NOT. FoundEps0) Eps0 = 8.854187817e-12_dp
      Tensor = Eps0 * Tensor
      RETURN
    END IF

    CALL ReadTensor(Material, Tensor, IsScalar, 'Electric Permittivity', &
                    Element, Found)
    IF (.NOT. Found) THEN
      CALL Fatal(Caller, 'Missing material property: Relative Permittivity ' // &
                         'or Electric Permittivity')
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetPermittivity
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReadTensor(Material, Tensor, IsScalar, Name, Element, Found)
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: Tensor(:, :, :)
    LOGICAL :: IsScalar, Found
    CHARACTER(LEN=*) :: Name
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    INTEGER :: i, j, n
    REAL(KIND=dp), POINTER :: Work(:, :, :) => NULL()
!$omp threadprivate(Work)
!------------------------------------------------------------------------------
    Tensor = 0.0_dp
    IsScalar = .TRUE.

    n = Element % TYPE % NumberOfNodes
    CALL ListGetRealArray(Material, Name, Work, n, Element % NodeIndexes, Found)
    IF (.NOT. Found) RETURN

    IsScalar = SIZE(Work, 1) == 1 .AND. SIZE(Work, 2) == 1
    IF (IsScalar) THEN
      DO i = 1, SIZE(Tensor, 1)
        Tensor(i, i, 1:n) = Work(1, 1, 1:n)
      END DO
    ELSE IF (SIZE(Work, 1) == 1) THEN
      DO i = 1, MIN(SIZE(Tensor, 1), SIZE(Work, 2))
        Tensor(i, i, 1:n) = Work(1, i, 1:n)
      END DO
    ELSE IF (SIZE(Work, 2) == 1) THEN
      DO i = 1, MIN(SIZE(Tensor, 1), SIZE(Work, 1))
        Tensor(i, i, 1:n) = Work(i, 1, 1:n)
      END DO
    ELSE
      DO i = 1, MIN(SIZE(Tensor, 1), SIZE(Work, 1))
        DO j = 1, MIN(SIZE(Tensor, 2), SIZE(Work, 2))
          Tensor(i, j, 1:n) = Work(i, j, 1:n)
        END DO
      END DO
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ReadTensor
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE QuasiElectrostaticSolver
!------------------------------------------------------------------------------
