!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a finite elasticity equation using OpenCMISS calls.
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
!> Contributor(s): Kumar Mithraratne, Adam Reeve
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

!> \example FiniteElasticity/SimpleShear/src/FortranExample.F90
!! Example program to solve a finite elasticity equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/SimpleShear/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/SimpleShear/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM SimpleShearExample

  USE OpenCMISS
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  REAL(OCRP) :: VALUE

  !Test program parameters

  REAL(OCRP), PARAMETER :: HEIGHT=1.0_OCRP
  REAL(OCRP), PARAMETER :: WIDTH=1.0_OCRP
  REAL(OCRP), PARAMETER :: LENGTH=1.0_OCRP
!  INTEGER(OCIntg), PARAMETER :: InterpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  INTEGER(OCIntg), PARAMETER :: InterpolationType=OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  INTEGER(OCIntg), PARAMETER :: PressureInterpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
!  LOGICAL, PARAMETER :: UsePressureBasis=.TRUE.
  LOGICAL, PARAMETER :: UsePressureBasis=.FALSE.
  INTEGER(OCIntg), PARAMETER :: NumberOfGaussXi=3

  INTEGER(OCIntg), PARAMETER :: ContextUserNumber=1
  INTEGER(OCIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(OCIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(OCIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(OCIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(OCIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(OCIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(OCIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(OCIntg), PARAMETER :: DecomposerUserNumber=1
  INTEGER(OCIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(OCIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(OCIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(OCIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(OCIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(OCIntg), PARAMETER :: EquationsSetFieldUserNumber=5
  INTEGER(OCIntg), PARAMETER :: ProblemUserNumber=1

  !Program types

  !Program variables

  INTEGER(OCIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(OCIntg) :: decompositionIndex,EquationsSetIndex
  INTEGER(OCIntg) :: NumberOfComputationNodes,NumberOfDomains,ComputationNodeNumber
  INTEGER(OCIntg) :: NodeNumber,NodeDomain,node_idx
  INTEGER(OCIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(OCIntg),ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(OCIntg),ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(OCIntg),ALLOCATABLE :: BackSurfaceNodes(:)
  INTEGER(OCIntg),ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(OCIntg),ALLOCATABLE :: TopSurfaceNodes(:)
  INTEGER(OCIntg) :: LeftNormalXi,RightNormalXi,FrontNormalXi,BackNormalXi,BottomNormalXi,TopNormalXi

  INTEGER(OCIntg), PARAMETER :: NUMBER_OF_COMPONENTS = 3 !nearly incompressible

  !CMISS variables
  TYPE(OC_BasisType) :: Basis, PressureBasis
  TYPE(OC_BoundaryConditionsType) :: BoundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: ComputationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: CoordinateSystem
  TYPE(OC_DecompositionType) :: Decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: Equations
  TYPE(OC_EquationsSetType) :: EquationsSet
  TYPE(OC_FieldType) :: GeometricField,FibreField,MaterialField,DependentField,EquationsSetField
  TYPE(OC_FieldsType) :: Fields
  TYPE(OC_GeneratedMeshType) :: GeneratedMesh
  TYPE(OC_MeshType) :: Mesh
  TYPE(OC_ProblemType) :: Problem
  TYPE(OC_RegionType) :: Region,WorldRegion
  TYPE(OC_SolverType) :: Solver,LinearSolver
  TYPE(OC_SolverEquationsType) :: SolverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables
  INTEGER(OCIntg) :: Err

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
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,Err)
  !Set all diganostic levels on for testing
  !CALL OC_DiagnosticsSetOn(OC_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_RESIDUAL_EVALUATE"],Err)
  CALL OC_OutputSetOn("SimpleShear",Err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(ContextUserNumber,context,Err)  
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)
  
  !Get the number of computation nodes and this computation node number
  CALL OC_ComputationEnvironment_Initialise(ComputationEnvironment,Err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationNodeNumber,err)

  NumberGlobalXElements=2
  NumberGlobalYElements=2
  NumberGlobalZElements=2

  !Create a 3D rectangular cartesian coordinate system
  CALL OC_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL OC_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,context,CoordinateSystem,Err)
  CALL OC_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL OC_Region_Initialise(Region,Err)
  CALL OC_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL OC_Region_LabelSet(Region,"Region",Err)
  CALL OC_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL OC_Region_CreateFinish(Region,Err)

  !Define geometric basis
  CALL OC_Basis_Initialise(Basis,Err)
  CALL OC_Basis_CreateStart(BasisUserNumber,context,Basis,Err)
  SELECT CASE(InterpolationType)
  CASE(1,2,3,4)
    CALL OC_Basis_TypeSet(Basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL OC_Basis_TypeSet(Basis,OC_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  IF(NumberGlobalZElements==0) THEN
    CALL OC_Basis_NumberOfXiSet(Basis,2,Err)
    CALL OC_Basis_InterpolationXiSet(Basis,[InterpolationType,InterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(Basis,[NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ELSE
    CALL OC_Basis_NumberOfXiSet(Basis,3,Err)
    CALL OC_Basis_InterpolationXiSet(Basis,[InterpolationType,InterpolationType,InterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(Basis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ENDIF
  CALL OC_Basis_CreateFinish(Basis,Err)

  !Define pressure basis
  IF(UsePressureBasis) THEN
    CALL OC_Basis_Initialise(PressureBasis,Err)
    CALL OC_Basis_CreateStart(PressureBasisUserNumber,context,PressureBasis,Err)
    SELECT CASE(PressureInterpolationType)
    CASE(1,2,3,4)
      CALL OC_Basis_TypeSet(PressureBasis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CASE(7,8,9)
      CALL OC_Basis_TypeSet(PressureBasis,OC_BASIS_SIMPLEX_TYPE,Err)
    END SELECT
    IF(NumberGlobalZElements==0) THEN
      CALL OC_Basis_NumberOfXiSet(PressureBasis,2,Err)
      CALL OC_Basis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType],Err)
      IF(NumberOfGaussXi>0) THEN
        CALL OC_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
      ENDIF
    ELSE
      CALL OC_Basis_NumberOfXiSet(PressureBasis,3,Err)
      CALL OC_Basis_InterpolationXiSet(PressureBasis, &
        & [PressureInterpolationType,PressureInterpolationType,PressureInterpolationType],Err)
      IF(NumberOfGaussXi>0) THEN
        CALL OC_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
      ENDIF
    ENDIF
    CALL OC_Basis_CreateFinish(PressureBasis,Err)
  ENDIF

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL OC_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL OC_GeneratedMesh_TypeSet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  IF(UsePressureBasis) THEN
    CALL OC_GeneratedMesh_BasisSet(GeneratedMesh,[Basis,PressureBasis],Err)
  ELSE
    CALL OC_GeneratedMesh_BasisSet(GeneratedMesh,[Basis],Err)
  ENDIF
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
      & NumberGlobalZElements],Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(Mesh,Err)
  CALL OC_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL OC_Decomposition_Initialise(Decomposition,Err)
  CALL OC_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL OC_Decomposition_CreateFinish(Decomposition,Err)

  !Decompose
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(decomposerUserNumber,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (default is geometry)
  CALL OC_Field_Initialise(GeometricField,Err)
  CALL OC_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL OC_Field_DecompositionSet(GeometricField,Decomposition,Err)
  CALL OC_Field_VariableLabelSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL OC_Field_ScalingTypeSet(GeometricField,OC_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL OC_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create a fibre field and attach it to the geometric field
  CALL OC_Field_Initialise(FibreField,Err)
  CALL OC_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL OC_Field_TypeSet(FibreField,OC_FIELD_FIBRE_TYPE,Err)
  CALL OC_Field_DecompositionSet(FibreField,Decomposition,Err)
  CALL OC_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL OC_Field_VariableLabelSet(FibreField,OC_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL OC_Field_CreateFinish(FibreField,Err)

  !Create the dependent field
  CALL OC_Field_Initialise(DependentField,Err)
  CALL OC_Field_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL OC_Field_TypeSet(DependentField,OC_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL OC_Field_DecompositionSet(DependentField,Decomposition,Err)
  CALL OC_Field_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL OC_Field_DependentTypeSet(DependentField,OC_FIELD_DEPENDENT_TYPE,Err)
  CALL OC_Field_NumberOfVariablesSet(DependentField,2,Err)
  CALL OC_Field_VariableLabelSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL OC_Field_NumberOfComponentsSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,Err)
  IF(UsePressureBasis) THEN
    !Set the pressure to be nodally based and use the second mesh component if required
    CALL OC_Field_ComponentInterpolationSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,4,OC_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL OC_Field_ComponentInterpolationSet(DependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & OC_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL OC_Field_ComponentMeshComponentSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,4,2,Err)
    CALL OC_Field_ComponentMeshComponentSet(DependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
  END IF
  CALL OC_Field_CreateFinish(DependentField,Err)

  !Create the material field
  CALL OC_Field_Initialise(MaterialField,Err)
  CALL OC_Field_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL OC_Field_TypeSet(MaterialField,OC_FIELD_MATERIAL_TYPE,Err)
  CALL OC_Field_DecompositionSet(MaterialField,Decomposition,Err)
  CALL OC_Field_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL OC_Field_NumberOfVariablesSet(MaterialField,1,Err)
  CALL OC_Field_VariableLabelSet(MaterialField,OC_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL OC_Field_NumberOfComponentsSet(MaterialField,OC_FIELD_U_VARIABLE_TYPE,3,Err)
  CALL OC_Field_CreateFinish(MaterialField,Err)

  !Create the equations_set
  CALL OC_Field_Initialise(EquationsSetField,Err)
  CALL OC_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[OC_EQUATIONS_SET_ELASTICITY_CLASS, &
    & OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,OC_EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE], &
!    & OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,OC_EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &
!    & OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,OC_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  CALL OC_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field
  CALL OC_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL OC_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field 
  CALL OC_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL OC_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 0.5 and 0.0 respectively. Third value is kappa (bulk modulus ???)
  CALL OC_Field_ComponentValuesInitialise(MaterialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,0.5_OCRP,Err)
  CALL OC_Field_ComponentValuesInitialise(MaterialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,0.0_OCRP,Err)
  CALL OC_Field_ComponentValuesInitialise(MaterialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 3,10000.0_OCRP,Err)

  !Create the equations set equations
  CALL OC_Equations_Initialise(Equations,Err)
  CALL OC_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL OC_Equations_SparsityTypeSet(Equations,OC_EQUATIONS_SPARSE_MATRICES,Err)
  CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_NO_OUTPUT,Err)
  CALL OC_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs 
  CALL OC_Field_ParametersToFieldParametersComponentCopy(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,Err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,Err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3,Err)

  !Define the problem
  CALL OC_Problem_Initialise(Problem,Err)
  CALL OC_Problem_CreateStart(ProblemUserNumber,context,[OC_PROBLEM_ELASTICITY_CLASS,OC_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & OC_PROBLEM_NO_SUBTYPE],Problem,Err)
  CALL OC_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(Problem,Err)
  CALL OC_Problem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL OC_Solver_Initialise(Solver,Err)
  CALL OC_Solver_Initialise(LinearSolver,Err)
  CALL OC_Problem_SolversCreateStart(Problem,Err)
  CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_PROGRESS_OUTPUT,Err)
  CALL OC_Solver_NewtonJacobianCalculationTypeSet(Solver,OC_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL OC_Solver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL OC_Solver_LinearTypeSet(LinearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL OC_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL OC_Solver_Initialise(Solver,Err)
  CALL OC_SolverEquations_Initialise(SolverEquations,Err)
  CALL OC_Problem_SolverEquationsCreateStart(Problem,Err)
  CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL OC_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL OC_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL OC_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL OC_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL OC_GeneratedMesh_SurfaceGet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_TOP_SURFACE,TopSurfaceNodes,TopNormalXi,Err)
  CALL OC_GeneratedMesh_SurfaceGet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi, &
    & Err)
  CALL OC_GeneratedMesh_SurfaceGet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL OC_GeneratedMesh_SurfaceGet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL OC_GeneratedMesh_SurfaceGet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)
  CALL OC_GeneratedMesh_SurfaceGet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_BACK_SURFACE,BackSurfaceNodes,BackNormalXi,Err)

  ! set z=WIDTH nodes to 10% x-displacement, no displacement in y- and z-direction
  DO node_idx=1,SIZE(TopSurfaceNodes,1)
    NodeNumber=TopSurfaceNodes(node_idx)
    CALL OC_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationNodeNumber) THEN
      ! x-direction
      CALL OC_Field_ParameterSetGetNode(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,&
        & VALUE,Err)
      CALL OC_BoundaryConditions_SetNode(BoundaryConditions,DependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & OC_BOUNDARY_CONDITION_FIXED,VALUE+0.1_OCRP*WIDTH,Err)
      ! y-direction
      CALL OC_Field_ParameterSetGetNode(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,2,&
        & VALUE,Err)
      CALL OC_BoundaryConditions_SetNode(BoundaryConditions,DependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & OC_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      ! z-direction
      CALL OC_Field_ParameterSetGetNode(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,3,&
        & VALUE,Err)
      CALL OC_BoundaryConditions_SetNode(BoundaryConditions,DependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & OC_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  ENDDO

  ! set z=0 nodes to no displacement in x-, y- and z-direction
  DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    NodeNumber=BottomSurfaceNodes(node_idx)
    CALL OC_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationNodeNumber) THEN
      ! x-direction
      CALL OC_Field_ParameterSetGetNode(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,&
        & VALUE,Err)
      CALL OC_BoundaryConditions_SetNode(BoundaryConditions,DependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & OC_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      ! y-direction
      CALL OC_Field_ParameterSetGetNode(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,2,&
        & VALUE,Err)
      CALL OC_BoundaryConditions_SetNode(BoundaryConditions,DependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & OC_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      ! z-direction
      CALL OC_Field_ParameterSetGetNode(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,3,&
        & VALUE,Err)
      CALL OC_BoundaryConditions_SetNode(BoundaryConditions,DependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & OC_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  ENDDO

  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL OC_Problem_Solve(Problem,Err)

  !Output solution
  CALL OC_Fields_Initialise(Fields,Err)
  CALL OC_Fields_Create(Region,Fields,Err)
  CALL OC_Fields_NodesExport(Fields,"SimpleShear","FORTRAN",Err)
  CALL OC_Fields_ElementsExport(Fields,"SimpleShear","FORTRAN",Err)
  CALL OC_Fields_Finalise(Fields,Err)

  CALL OC_Context_Destroy(context,Err)
  CALL OC_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM SimpleShearExample

