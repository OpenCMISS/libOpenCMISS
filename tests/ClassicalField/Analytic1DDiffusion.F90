!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a 1D diffusion equation and compare it to the analytic solution using OpenCMISS calls.
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

!> \example ClassicalField/Diffusion/Analytic1DDiffusion/src/Analytic1DDiffusionExample.F90
!! Example program to solve a 1D diffusion equation using OpenCMISS calls.
!!
!! \htmlinclude ClassicalField/Diffusion/Analytic1DDiffusion/history.html
!<

!> Main program
PROGRAM Analytic1DDiffusionExample

  USE OpenCMISS

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(OC_RP), PARAMETER ::    PI=3.141592653589793238462643383279502884197_OC_RP

  INTEGER(OC_Intg), PARAMETER :: NUMBER_GLOBAL_X_ELEMENTS=6
  REAL(OC_RP), PARAMETER :: LENGTH=3.0_OC_RP
  REAL(OC_RP), PARAMETER :: END_TIME=0.1_OC_RP
  REAL(OC_RP), PARAMETER :: TIME_STEP=0.01_OC_RP
  REAL(OC_RP), PARAMETER :: A=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: B=PI/2.0_OC_RP
  REAL(OC_RP), PARAMETER :: C=0.0_OC_RP
  REAL(OC_RP), PARAMETER :: K=1.0_OC_RP
  
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
  INTEGER(OC_Intg), PARAMETER :: EquationsSetUserNumber=12
  INTEGER(OC_Intg), PARAMETER :: EquationsSetFieldUserNumber=13
  INTEGER(OC_Intg), PARAMETER :: ProblemUserNumber=14
  INTEGER(OC_Intg), PARAMETER :: AnalyticFieldUserNumber=15
  !Program types
  
  !Program variables
  
  !CMISS variables

  TYPE(OC_BasisType) :: Basis
  TYPE(OC_BoundaryConditionsType) :: BoundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: CoordinateSystem
  TYPE(OC_DecompositionType) :: Decomposition
  TYPE(OC_DecomposerType) :: Decomposer
  TYPE(OC_EquationsType) :: Equations
  TYPE(OC_EquationsSetType) :: EquationsSet
  TYPE(OC_FieldType) :: GeometricField,DependentField,EquationsSetField,MaterialsField,AnalyticField
  TYPE(OC_FieldsType) :: Fields
  TYPE(OC_GeneratedMeshType) :: GeneratedMesh  
  TYPE(OC_MeshType) :: Mesh
  TYPE(OC_ProblemType) :: Problem
  TYPE(OC_ControlLoopType) :: ControlLoop
  TYPE(OC_RegionType) :: Region,WorldRegion
  TYPE(OC_SolverType) :: Solver, LinearSolver
  TYPE(OC_SolverEquationsType) :: SolverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables

  INTEGER(OC_Intg) :: numberOfComputationNodes, computationNodeNumber
  INTEGER(OC_Intg) :: decompositionIndex, EquationsSetIndex
  INTEGER(OC_Intg) :: err
  
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise OpenCMISS
  CALL OC_Initialise(Err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  CALL OC_OutputSetOn("Diffusion1DAnalytic",err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(ContextUserNumber,context,Err)

  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)
    
  !Get the computation nodes information
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationNodeNumber,err)
  
  !Start the creation of a new RC coordinate system
  CALL OC_CoordinateSystem_Initialise(CoordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,context,CoordinateSystem,err)
  !Set the coordinate system to be 1D
  CALL OC_CoordinateSystem_DimensionSet(CoordinateSystem,1,err)
  !Finish the creation of the coordinate system
  CALL OC_CoordinateSystem_CreateFinish(CoordinateSystem,err)
  
  !Start the creation of the region
  CALL OC_Region_Initialise(Region,err)
  CALL OC_Region_CreateStart(RegionUserNumber,WorldRegion,Region,err)
  !Label the Region
  CALL OC_Region_LabelSet(Region,"Region",err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL OC_Region_CoordinateSystemSet(Region,CoordinateSystem,err)
  !Finish the creation of the region
  CALL OC_Region_CreateFinish(Region,err)
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL OC_Basis_Initialise(Basis,err)
  CALL OC_Basis_CreateStart(BasisUserNumber,context,Basis,err)
  !Set the basis to be a linear Lagrange basis
  CALL OC_Basis_NumberOfXiSet(Basis,1,err)
  !Finish the creation of the basis
  CALL OC_Basis_CreateFinish(BASIS,err)
  
  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(GeneratedMesh,err)
  CALL OC_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,err)
  !Set up a regular mesh
  CALL OC_GeneratedMesh_TypeSet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL OC_GeneratedMesh_BasisSet(GeneratedMesh,Basis,err)   
  !Define the mesh on the region
  CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH],err)
  CALL OC_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS],err)
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(Mesh,err)
  CALL OC_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,err)
  
  !Create a decomposition
  CALL OC_Decomposition_Initialise(Decomposition,err)
  CALL OC_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,err)
  !Finish the decomposition
  CALL OC_Decomposition_CreateFinish(Decomposition,err)
  
  !Decompose
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(decomposerUserNumber,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !Start to create a default (geometric) field on the region
  CALL OC_Field_Initialise(GeometricField,err)
  CALL OC_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(GeometricField,Decomposition,err)
  !Set the domain to be used by the field components.
  CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,1,1,err)
  !Finish creating the field
  CALL OC_Field_CreateFinish(GeometricField,err)
       
  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,err)
  
  !Create the equations_set
  CALL OC_EquationsSet_Initialise(EquationsSet,err)
  CALL OC_Field_Initialise(EquationsSetField,err)
  CALL OC_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & OC_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,OC_EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,err)
  !Finish creating the equations set
  CALL OC_EquationsSet_CreateFinish(EquationsSet,err)

  !Create the equations set dependent field variables
  CALL OC_Field_Initialise(DependentField,err)
  CALL OC_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_DependentCreateFinish(EquationsSet,err)

  !Create the equations set material field variables
  CALL OC_Field_Initialise(MaterialsField,err)
  CALL OC_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_MaterialsCreateFinish(EquationsSet,err)
  !Set the conductivity
  CALL OC_Field_ComponentValuesInitialise(MaterialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,K,err)
 
  !Create the equations set analytic field variables
  CALL OC_Field_Initialise(AnalyticField,err)
  CALL OC_EquationsSet_AnalyticCreateStart(EquationsSet,OC_EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1, & 
    & AnalyticFieldUserNumber,AnalyticField,err)
  !Finish the equations set analytic field variables
  CALL OC_EquationsSet_AnalyticCreateFinish(EquationsSet,err)
  !Set the analytic field parameters
  !Set the multiplicative constant. 
  CALL OC_Field_ComponentValuesInitialise(AnalyticField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,A,err)
  !Set the phase. 
  CALL OC_Field_ComponentValuesInitialise(AnalyticField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,B,err)
  !Set the offset. 
  CALL OC_Field_ComponentValuesInitialise(AnalyticField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3,C,err)
  !Set the length to be the length of the mesh
  CALL OC_Field_ComponentValuesInitialise(AnalyticField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,4,LENGTH,err)
  
  !Create the equations set equations
  CALL OC_Equations_Initialise(Equations,err)
  CALL OC_EquationsSet_EquationsCreateStart(EquationsSet,Equations,err)
  !Set the equations matrices sparsity type
  CALL OC_Equations_SparsityTypeSet(Equations,OC_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  !CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_NO_OUTPUT,err)
  !CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_TIMING_OUTPUT,err)
  !CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_MATRIX_OUTPUT,err)
  CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL OC_EquationsSet_EquationsCreateFinish(EquationsSet,err)

  !Create the problem
  CALL OC_Problem_Initialise(Problem,err)
  CALL OC_Problem_CreateStart(ProblemUserNumber,context,[OC_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & OC_PROBLEM_DIFFUSION_EQUATION_TYPE,OC_PROBLEM_LINEAR_DIFFUSION_SUBTYPE],Problem,err)
  !Finish the creation of a problem.
  CALL OC_Problem_CreateFinish(Problem,err)

  !Create the problem control
  CALL OC_Problem_ControlLoopCreateStart(Problem,err)
  CALL OC_ControlLoop_Initialise(ControlLoop,err)
  !Get the control loop
  CALL OC_Problem_ControlLoopGet(Problem,OC_CONTROL_LOOP_NODE,ControlLoop,err)
  !Set the times
  CALL OC_ControlLoop_TimesSet(ControlLoop,0.0_OC_RP,END_TIME,TIME_STEP,err)
  !Finish creating the problem control loop
  CALL OC_Problem_ControlLoopCreateFinish(Problem,err)

  !Start the creation of the problem solvers
  CALL OC_Solver_Initialise(Solver,err)
  CALL OC_Solver_Initialise(LinearSolver,err)
  CALL OC_Problem_SolversCreateStart(Problem,err)
  CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,err)
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_NO_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_PROGRESS_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_TIMING_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_SOLVER_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_MATRIX_OUTPUT,err)
  CALL OC_Solver_DynamicLinearSolverGet(Solver,LinearSolver,err)
  CALL OC_Solver_OutputTypeSet(LinearSolver,OC_SOLVER_PROGRESS_OUTPUT,err)
  !Finish the creation of the problem solver
  CALL OC_Problem_SolversCreateFinish(Problem,err)

  !Create the problem solver equations
  CALL OC_Solver_Initialise(Solver,err)
  CALL OC_SolverEquations_Initialise(SolverEquations,err)
  CALL OC_Problem_SolverEquationsCreateStart(Problem,err)
  !Get the solve equations
  CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,err)
  CALL OC_Solver_SolverEquationsGet(Solver,SolverEquations,err)
  !Set the solver equations sparsity
  CALL OC_SolverEquations_SparsityTypeSet(SolverEquations,OC_SOLVER_SPARSE_MATRICES,err)
  !CALL OC_SolverEquations_SparsityTypeSet(SolverEquations,OC_SOLVER_FULL_MATRICES,err)  
  !Add in the equations set
  CALL OC_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateFinish(Problem,err)

  !Create the equations set boundary conditions
  CALL OC_BoundaryConditions_Initialise(BoundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsAnalytic(SolverEquations,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,err)

  !Solve the problem
  CALL OC_Problem_Solve(Problem,err)

  !Output Analytic analysis
  !CALL OC_EquationsSet_AnalyticTimeSet(EquationsSet,END_TIME,err)
  !CALL OC_EquationsSet_AnalyticEvaluate(EquationsSet,err)
  CALL OC_AnalyticAnalysis_Output(DependentField,"Diffusion1DAnalytic",err)

  !Output fields
  CALL OC_Fields_Initialise(Fields,err)
  CALL OC_Fields_Create(Region,Fields,err)
  CALL OC_Fields_NodesExport(Fields,"Diffusion1DAnalytic","FORTRAN",err)
  CALL OC_Fields_ElementsExport(Fields,"Diffusion1DAnalytic","FORTRAN",err)
  CALL OC_Fields_Finalise(Fields,err)

  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finalise and quit
  CALL OC_Finalise(err)
  WRITE(*,'(A)') "Program successfully completed."

  STOP
  
END PROGRAM Analytic1DDiffusionExample
