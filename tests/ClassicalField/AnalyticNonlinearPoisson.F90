!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a nonlinear Poisson equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example ClassicalField/Poisson/AnalyticNonlinearPoisson/src/NonlinearPoissonExample.F90
!! Example program to solve a nonlinear Poisson equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/Poisson/AnalyticNonlinearPoisson/history.html
!<

!> Main program
PROGRAM NonlinearPoissonExample

  USE OpenCMISS

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(OC_RP), PARAMETER :: HEIGHT=0.5_OC_RP
  REAL(OC_RP), PARAMETER :: WIDTH=0.5_OC_RP
  REAL(OC_RP), PARAMETER :: LENGTH=1.0_OC_RP

  INTEGER(OC_Intg), PARAMETER :: ContextUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: CoordinateSystemUserNumber=2
  INTEGER(OC_Intg), PARAMETER :: RegionUserNumber=3
  INTEGER(OC_Intg), PARAMETER :: BasisUserNumber=4
  INTEGER(OC_Intg), PARAMETER :: GeneratedMeshUserNumber=5
  INTEGER(OC_Intg), PARAMETER :: MeshUserNumber=6
  INTEGER(OC_Intg), PARAMETER :: DecompositionUserNumber=7
  INTEGER(OC_Intg), PARAMETER :: DecomposerUserNumber=8
  INTEGER(OC_Intg), PARAMETER :: GeometricFieldUserNumber=9
  INTEGER(OC_Intg), PARAMETER :: DependentFieldUserNumber=10
  INTEGER(OC_Intg), PARAMETER :: MaterialsFieldUserNumber=11
  INTEGER(OC_Intg), PARAMETER :: AnalyticFieldUserNumber=12
  INTEGER(OC_Intg), PARAMETER :: EquationsSetUserNumber=13
  INTEGER(OC_Intg), PARAMETER :: ProblemUserNumber=14
  INTEGER(OC_Intg), PARAMETER :: EquationsSetFieldUserNumber=15

  !Program variables

  INTEGER(OC_Intg) :: NUMBER_DIMENSIONS,INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI
  INTEGER(OC_Intg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(OC_Intg) :: component_idx
  INTEGER(OC_Intg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT

  LOGICAL :: EXPORT_FIELD

  !CMISS variables

  TYPE(OC_BasisType) :: Basis
  TYPE(OC_BoundaryConditionsType) :: BoundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: ComputationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(OC_DecompositionType) :: Decomposition
  TYPE(OC_DecomposerType) :: Decomposer
  TYPE(OC_EquationsType) :: Equations
  TYPE(OC_EquationsSetType) :: EquationsSet
  TYPE(OC_FieldType) :: EquationsSetField,GeometricField,DependentField,MaterialsField,AnalyticField
  TYPE(OC_FieldsType) :: Fields
  TYPE(OC_GeneratedMeshType) :: GeneratedMesh
  TYPE(OC_MeshType) :: Mesh
  TYPE(OC_ProblemType) :: Problem
  TYPE(OC_RegionType) :: Region,WorldRegion
  TYPE(OC_SolverType) :: Solver,LinearSolver
  TYPE(OC_SolverEquationsType) :: SolverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables

  INTEGER(OC_Intg) :: DecompositionIndex,EquationsSetIndex
  INTEGER(OC_Intg) :: Err
  INTEGER(OC_Intg) :: NumberOfComputationNodes,ComputationNodeNumber

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG

  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Get input arguments
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS >= 4) THEN
    !If we have enough arguments then use the first four for setting up the problem. The subsequent arguments may be used to
    !pass flags to, say, PETSc.
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HandleError("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_X_ELEMENTS
    IF(NUMBER_GLOBAL_X_ELEMENTS<=0) CALL HandleError("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HandleError("Error for command argument 2.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Y_ELEMENTS
    IF(NUMBER_GLOBAL_Y_ELEMENTS<0) CALL HandleError("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HandleError("Error for command argument 3.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Z_ELEMENTS
    IF(NUMBER_GLOBAL_Z_ELEMENTS<0) CALL HandleError("Invalid number of Z elements.")
    CALL GET_COMMAND_ARGUMENT(4,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HandleError("Error for command argument 4.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION_TYPE
    IF(INTERPOLATION_TYPE<=0) CALL HandleError("Invalid Interpolation specification.")
    IF(NUMBER_GLOBAL_Z_ELEMENTS>0) THEN
      NUMBER_DIMENSIONS=3
    ELSEIF(NUMBER_GLOBAL_Y_ELEMENTS>0) THEN
      NUMBER_DIMENSIONS=2
    ELSE
      NUMBER_DIMENSIONS=1
    ENDIF
  ELSE
    !If there are not enough arguments default the problem specification
    NUMBER_DIMENSIONS=2
    NUMBER_GLOBAL_X_ELEMENTS=5
    NUMBER_GLOBAL_Y_ELEMENTS=5
    NUMBER_GLOBAL_Z_ELEMENTS=0
    INTERPOLATION_TYPE=1
  ENDIF

  !Intialise OpenCMISS
  CALL OC_Initialise(Err)  
  !Trap all errors
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,Err)
  !Output to a file
  CALL OC_OutputSetOn("NonlinearPoisson",Err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(ContextUserNumber,context,Err)  
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  !Get the computation nodes information
  CALL OC_ComputationEnvironment_Initialise(ComputationEnvironment,Err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationNodeNumber,err)

  !Start the creation of a new RC coordinate system
  CALL OC_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL OC_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,context,CoordinateSystem,Err)
  !Set the coordinate system number of dimensions
  CALL OC_CoordinateSystem_DimensionSet(CoordinateSystem,NUMBER_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL OC_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL OC_Region_Initialise(Region,Err)
  CALL OC_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL OC_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL OC_Region_CreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL OC_Basis_Initialise(Basis,Err)
  CALL OC_Basis_CreateStart(BasisUserNumber,context,Basis,Err)
  CALL OC_Basis_NumberOfXiSet(Basis,NUMBER_DIMENSIONS,Err)
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1,2,3,4)
    CALL OC_Basis_TypeSet(Basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL OC_Basis_TypeSet(Basis,OC_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HandleError("Invalid interpolation type.")
  END SELECT
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1)
    NUMBER_OF_GAUSS_XI=2
  CASE(2)
    NUMBER_OF_GAUSS_XI=3
  CASE(3,4)
    NUMBER_OF_GAUSS_XI=4
  CASE DEFAULT
    NUMBER_OF_GAUSS_XI=0
  END SELECT
  IF(NUMBER_DIMENSIONS==1) THEN
    CALL OC_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSEIF(NUMBER_DIMENSIONS==2) THEN
    CALL OC_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSE
    CALL OC_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL OC_Basis_CreateFinish(Basis,Err)

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL OC_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL OC_GeneratedMesh_TypeSet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL OC_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)
  !Define the mesh on the region
  IF(NUMBER_DIMENSIONS==1) THEN
    CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH],Err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS],Err)
  ELSEIF(NUMBER_DIMENSIONS==2) THEN
    CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(Mesh,Err)
  CALL OC_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL OC_Decomposition_Initialise(Decomposition,Err)
  CALL OC_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Finish the decomposition
  CALL OC_Decomposition_CreateFinish(Decomposition,Err)

  !Decompose
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(decomposerUserNumber,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !Start to create a default (geometric) field on the region
  CALL OC_Field_Initialise(GeometricField,Err)
  CALL OC_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  DO component_idx=1,NUMBER_DIMENSIONS
    CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
  ENDDO
  !Finish creating the field
  CALL OC_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create the equations_set
  CALL OC_EquationsSet_Initialise(EquationsSet,Err)
  CALL OC_Field_Initialise(EquationsSetField,Err)
  CALL OC_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & OC_EQUATIONS_SET_POISSON_EQUATION_TYPE,OC_EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE], &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL OC_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL OC_Field_Initialise(DependentField,Err)
  CALL OC_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Initialise the field to zero
  CALL OC_Field_ComponentValuesInitialise(DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,0.0_OC_RP, &
    & Err)

  !Create the equations set material field variables
  CALL OC_Field_Initialise(MaterialsField,Err)
  CALL OC_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set analytic field variables
  CALL OC_Field_Initialise(AnalyticField,Err)
  IF(NUMBER_DIMENSIONS==2) THEN
    CALL OC_EquationsSet_AnalyticCreateStart(EquationsSet,OC_EQUATIONS_SET_EXPONENTIAL_POISSON_EQUATION_TWO_DIM_1, &
      & AnalyticFieldUserNumber,AnalyticField,Err)
  ELSE
    WRITE(*,'(A)') "One and three dimensions are not implemented."
    STOP
  ENDIF
  !Finish the equations set analytic field variables
  CALL OC_EquationsSet_AnalyticCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL OC_Equations_Initialise(Equations,Err)
  CALL OC_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL OC_Equations_SparsityTypeSet(Equations,OC_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_NO_OUTPUT,Err)
  !Finish the equations set equations
  CALL OC_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Start the creation of a problem.
  CALL OC_Problem_Initialise(Problem,Err)
  CALL OC_Problem_CreateStart(ProblemUserNumber,context,[OC_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & OC_PROBLEM_POISSON_EQUATION_TYPE,OC_PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL OC_Problem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL OC_Problem_ControlLoopCreateFinish(Problem,Err)

  !Start the creation of the problem solvers
  CALL OC_Solver_Initialise(Solver,Err)
  CALL OC_Solver_Initialise(LinearSolver,Err)
  CALL OC_Problem_SolversCreateStart(Problem,Err)
  CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
  !Set the solver output
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_NO_OUTPUT,Err)
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_TIMING_OUTPUT,Err)
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_SOLVER_OUTPUT,Err)
  CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_MATRIX_OUTPUT,Err)
  !Set the Jacobian type
  CALL OC_Solver_NewtonJacobianCalculationTypeSet(Solver,OC_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  !CALL OC_Solver_NewtonJacobianCalculationTypeSet(Solver,OC_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  !Get the associated linear solver
  CALL OC_Solver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL OC_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,500,Err)
  !Finish the creation of the problem solver
  CALL OC_Problem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL OC_Solver_Initialise(Solver,Err)
  CALL OC_SolverEquations_Initialise(SolverEquations,Err)
  CALL OC_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL OC_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL OC_SolverEquations_SparsityTypeSet(SolverEquations,OC_SOLVER_SPARSE_MATRICES,Err)
  !CALL OC_SolverEquations_SparsityTypeSet(SolverEquations,OC_SOLVER_FULL_MATRICES,Err)
  !Add in the equations set
  CALL OC_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Set up the boundary conditions as per the analytic solution
  CALL OC_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  CALL OC_SolverEquations_BoundaryConditionsAnalytic(SolverEquations,Err)
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL OC_Problem_Solve(Problem,Err)

  !Output Analytic analysis
  Call OC_AnalyticAnalysis_Output(DependentField,"",Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL OC_Fields_Initialise(Fields,Err)
    CALL OC_Fields_Create(Region,Fields,Err)
    CALL OC_Fields_NodesExport(Fields,"NonlinearPoisson","FORTRAN",Err)
    CALL OC_Fields_ElementsExport(Fields,"NonlinearPoisson","FORTRAN",Err)
    CALL OC_Fields_Finalise(Fields,Err)
  ENDIF

  !Destroy the context
  CALL OC_Context_Destroy(context,Err)
  !Finialise OpenCMISS
  CALL OC_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HandleError(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HandleError

END PROGRAM NonlinearPoissonExample
