!> \file
!> \author Adam Reeve
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
!> Contributor(s): Adam Reeve
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

!> \example FiniteElasticity/Cantilever/src/CantileverExample.F90
!! Example program to solve a finite elasticity equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM CantileverExample

  USE OpenCMISS

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters
  REAL(OC_RP), PARAMETER :: Width=60.0_OC_RP
  REAL(OC_RP), PARAMETER :: Length=40.0_OC_RP
  REAL(OC_RP), PARAMETER :: Height=40.0_OC_RP
  INTEGER(OC_Intg) :: DisplacementInterpolationType
  INTEGER(OC_Intg) :: PressureInterpolationType
  INTEGER(OC_Intg) :: PressureMeshComponent
  INTEGER(OC_Intg) :: NumberOfGaussXi
  INTEGER(OC_Intg) :: ScalingType
  REAL(OC_RP), PARAMETER :: Density=9.0E-4_OC_RP !in g mm^-3
  REAL(OC_RP), PARAMETER :: Gravity(3)=[0.0_OC_RP,0.0_OC_RP,-9.8_OC_RP] !in m s^-2
  INTEGER(OC_Intg), PARAMETER :: NumberOfLoadIncrements=2

  INTEGER(OC_Intg), PARAMETER :: ContextUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: RegionUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: DisplacementBasisUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(OC_Intg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: MeshUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: DecomposerUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(OC_Intg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(OC_Intg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(OC_Intg), PARAMETER :: FieldSourceUserNumber=5
  INTEGER(OC_Intg), PARAMETER :: EquationsSetFieldUserNumber=6
  INTEGER(OC_Intg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: ProblemUserNumber=1

  !Program variables
  INTEGER(OC_Intg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(OC_Intg) :: decompositionIndex,EquationsSetIndex
  INTEGER(OC_Intg) :: NumberOfComputationNodes,NumberOfDomains,ComputationNodeNumber
  INTEGER(OC_Intg) :: NodeNumber,NodeDomain,node_idx,component_idx,deriv_idx
  INTEGER(OC_Intg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(OC_Intg) :: LeftNormalXi
  INTEGER(OC_Intg) :: NumberOfArguments,ArgumentLength,ArgStatus
  CHARACTER(LEN=255) :: CommandArgument

  !CMISS variables
  TYPE(OC_BasisType) :: DisplacementBasis,PressureBasis
  TYPE(OC_BoundaryConditionsType) :: BoundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: ComputationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_ControlLoopType) :: ControlLoop
  TYPE(OC_CoordinateSystemType) :: CoordinateSystem
  TYPE(OC_MeshType) :: Mesh
  TYPE(OC_DecompositionType) :: Decomposition
  TYPE(OC_DecomposerType) :: Decomposer
  TYPE(OC_EquationsType) :: Equations
  TYPE(OC_EquationsSetType) :: EquationsSet
  TYPE(OC_FieldType) :: GeometricField,FibreField,MaterialField,DependentField,SourceField,EquationsSetField
  TYPE(OC_FieldsType) :: Fields
  TYPE(OC_GeneratedMeshType) :: GeneratedMesh
  TYPE(OC_ProblemType) :: Problem
  TYPE(OC_RegionType) :: Region,WorldRegion
  TYPE(OC_SolverType) :: Solver,LinearSolver
  TYPE(OC_SolverEquationsType) :: SolverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  INTEGER(OC_Intg) :: Err

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

  !Intialise OpenCMISS
  CALL OC_Initialise(Err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,Err)
  CALL OC_OutputSetOn("Cantilever",Err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(ContextUserNumber,context,Err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  !Read in arguments and overwrite default values
  !Usage: CantileverExample [Displacement Interpolation Type] [X elements] [Y elements] [Z elements] [Scaling Type]
  !Defaults:
  DisplacementInterpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  NumberGlobalXElements=3
  NumberGlobalYElements=2
  NumberGlobalZElements=2
  ScalingType=OC_FIELD_ARITHMETIC_MEAN_SCALING

  NumberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(NumberOfArguments >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HandleError("Error for command argument 1.")
    READ(CommandArgument(1:ArgumentLength),*) DisplacementInterpolationType
  ENDIF
  IF(NumberOfArguments >= 2) THEN
    CALL GET_COMMAND_ARGUMENT(2,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HandleError("Error for command argument 2.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalXElements
    IF(NumberGlobalXElements<1) CALL HandleError("Invalid number of X elements.")
  ENDIF
  IF(NumberOfArguments >= 3) THEN
    CALL GET_COMMAND_ARGUMENT(3,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HandleError("Error for command argument 3.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalYElements
    IF(NumberGlobalYElements<1) CALL HandleError("Invalid number of Y elements.")
  ENDIF
  IF(NumberOfArguments >= 4) THEN
    CALL GET_COMMAND_ARGUMENT(4,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HandleError("Error for command argument 4.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalZElements
    IF(NumberGlobalZElements<1) CALL HandleError("Invalid number of Z elements.")
  ENDIF
  IF(DisplacementInterpolationType==OC_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
    IF(NumberOfArguments >= 5) THEN
      CALL GET_COMMAND_ARGUMENT(5,CommandArgument,ArgumentLength,ArgStatus)
      IF(ArgStatus>0) CALL HandleError("Error for command argument 5.")
      READ(CommandArgument(1:ArgumentLength),*) ScalingType
      IF(ScalingType<0.OR.ScalingType>5) CALL HandleError("Invalid scaling type.")
    ENDIF
  ELSE
    ScalingType=OC_FIELD_NO_SCALING
  ENDIF
  SELECT CASE(DisplacementInterpolationType)
  CASE(OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=2
    PressureMeshComponent=1
  CASE(OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=3
    PressureMeshComponent=2
    PressureInterpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  CASE(OC_BASIS_CUBIC_LAGRANGE_INTERPOLATION,OC_BASIS_CUBIC_HERMITE_INTERPOLATION)
    NumberOfGaussXi=4
    PressureMeshComponent=2
    !Should generally use quadratic interpolation but use linear to match CMISS example 5e
    PressureInterpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  CASE DEFAULT
    NumberOfGaussXi=0
    PressureMeshComponent=1
  END SELECT
  WRITE(*,'("Interpolation: ", i3)') DisplacementInterpolationType
  WRITE(*,'("Elements: ", 3 i3)') NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  WRITE(*,'("Scaling type: ", i3)') ScalingType

  !Get the number of computation nodes and this computation node number
  CALL OC_ComputationEnvironment_Initialise(ComputationEnvironment,Err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationNodeNumber,err)

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

  !Define basis function for displacement
  CALL OC_Basis_Initialise(DisplacementBasis,Err)
  CALL OC_Basis_CreateStart(DisplacementBasisUserNumber,context,DisplacementBasis,Err)
  SELECT CASE(DisplacementInterpolationType)
  CASE(1,2,3,4)
    CALL OC_Basis_TypeSet(DisplacementBasis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL OC_Basis_TypeSet(DisplacementBasis,OC_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  CALL OC_Basis_NumberOfXiSet(DisplacementBasis,3,Err)
  CALL OC_Basis_InterpolationXiSet(DisplacementBasis,[DisplacementInterpolationType,DisplacementInterpolationType, &
      & DisplacementInterpolationType],Err)
  IF(NumberOfGaussXi>0) THEN
    CALL OC_Basis_QuadratureNumberOfGaussXiSet(DisplacementBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL OC_Basis_CreateFinish(DisplacementBasis,Err)

  IF(PressureMeshComponent/=1) THEN
    !Basis for pressure
    CALL OC_Basis_Initialise(PressureBasis,Err)
    CALL OC_Basis_CreateStart(PressureBasisUserNumber,context,PressureBasis,Err)
    SELECT CASE(PressureInterpolationType)
    CASE(1,2,3,4)
      CALL OC_Basis_TypeSet(PressureBasis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CASE(7,8,9)
      CALL OC_Basis_TypeSet(PressureBasis,OC_BASIS_SIMPLEX_TYPE,Err)
    END SELECT
    CALL OC_Basis_NumberOfXiSet(PressureBasis,3,Err)
    CALL OC_Basis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType, &
        & PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
    CALL OC_Basis_CreateFinish(PressureBasis,Err)
  ENDIF

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL OC_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL OC_GeneratedMesh_TypeSet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  IF(PressureMeshComponent==1) THEN
    CALL OC_GeneratedMesh_BasisSet(GeneratedMesh,[DisplacementBasis],Err)
  ELSE
    CALL OC_GeneratedMesh_BasisSet(GeneratedMesh,[DisplacementBasis,PressureBasis],Err)
  ENDIF
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,[Width,Height],Err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,[Width,Length,Height],Err)
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
  
  !Create a field to put the geometry (defualt is geometry)
  CALL OC_Field_Initialise(GeometricField,Err)
  CALL OC_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL OC_Field_DecompositionSet(GeometricField,Decomposition,Err)
  CALL OC_Field_VariableLabelSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL OC_Field_ScalingTypeSet(GeometricField,ScalingType,Err)
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
  CALL OC_Field_ScalingTypeSet(FibreField,ScalingType,Err)
  CALL OC_Field_CreateFinish(FibreField,Err)

  !Create the equations_set
  CALL OC_Field_Initialise(EquationsSetField,Err)
  CALL OC_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[OC_EQUATIONS_SET_ELASTICITY_CLASS, &
    & OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,OC_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  CALL OC_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL OC_Field_Initialise(DependentField,Err)
  CALL OC_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL OC_Field_VariableLabelSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  DO component_idx=1,3
    CALL OC_Field_ComponentMeshComponentSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
    CALL OC_Field_ComponentMeshComponentSet(DependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,component_idx,1,Err)
  ENDDO
  CALL OC_Field_ComponentMeshComponentSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,4,PressureMeshComponent,Err)
  CALL OC_Field_ComponentMeshComponentSet(DependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,4,PressureMeshComponent,Err)
  IF(PressureMeshComponent==1) THEN
    CALL OC_Field_ComponentInterpolationSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,4, &
      & OC_FIELD_ELEMENT_BASED_INTERPOLATION, &
      & Err)
    CALL OC_Field_ComponentInterpolationSet(DependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & OC_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ELSE
    CALL OC_Field_ComponentInterpolationSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,4,OC_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL OC_Field_ComponentInterpolationSet(DependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & OC_FIELD_NODE_BASED_INTERPOLATION,Err)
  ENDIF
  CALL OC_Field_ScalingTypeSet(DependentField,ScalingType,Err)
  CALL OC_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL OC_Field_Initialise(MaterialField,Err)
  CALL OC_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL OC_Field_VariableLabelSet(MaterialField,OC_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL OC_Field_VariableLabelSet(MaterialField,OC_FIELD_V_VARIABLE_TYPE,"Density",Err)
  CALL OC_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively
  CALL OC_Field_ComponentValuesInitialise(MaterialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2.0_OC_RP,Err)
  CALL OC_Field_ComponentValuesInitialise(MaterialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,6.0_OC_RP,Err)
  CALL OC_Field_ComponentValuesInitialise(MaterialField,OC_FIELD_V_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,Density,Err)

  !Create the source field with the gravity vector
  CALL OC_Field_Initialise(SourceField,Err)
  CALL OC_EquationsSet_SourceCreateStart(EquationsSet,FieldSourceUserNumber,SourceField,Err)
  CALL OC_Field_ScalingTypeSet(SourceField,ScalingType,Err)
  CALL OC_EquationsSet_SourceCreateFinish(EquationsSet,Err)
  DO component_idx=1,3
    CALL OC_Field_ComponentValuesInitialise(SourceField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
        & component_idx,Gravity(component_idx),Err)
  ENDDO

  !Create the equations set equations
  CALL OC_Equations_Initialise(Equations,Err)
  CALL OC_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL OC_Equations_SparsityTypeSet(Equations,OC_EQUATIONS_SPARSE_MATRICES,Err)
  CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_NO_OUTPUT,Err)
  CALL OC_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL OC_Field_ParametersToFieldParametersComponentCopy(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,Err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,Err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3,Err)
  CALL OC_Field_ComponentValuesInitialise(DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,4, &
    & -14.0_OC_RP,Err)
  CALL OC_Field_ParameterSetUpdateStart(DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,Err)
  CALL OC_Field_ParameterSetUpdateFinish(DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,Err)

  !Define the problem
  CALL OC_Problem_Initialise(Problem,Err)
  CALL OC_Problem_CreateStart(ProblemUserNumber,context,[OC_PROBLEM_ELASTICITY_CLASS,OC_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & OC_PROBLEM_NO_SUBTYPE],Problem,Err)
  CALL OC_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(Problem,Err)
  CALL OC_ControlLoop_Initialise(ControlLoop,Err)
  CALL OC_Problem_ControlLoopGet(Problem,OC_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL OC_ControlLoop_MaximumIterationsSet(ControlLoop,NumberOfLoadIncrements,Err)
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
  CALL OC_SolverEquations_SparsityTypeSet(SolverEquations,OC_SOLVER_SPARSE_MATRICES,Err)
  CALL OC_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL OC_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL OC_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL OC_GeneratedMesh_SurfaceGet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)

  !Fix x=0 nodes in x, y and z
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL OC_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationNodeNumber) THEN
      DO component_idx=1,3
        CALL OC_BoundaryConditions_AddNode(BoundaryConditions,DependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & component_idx,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,Err)
        IF(DisplacementInterpolationType==OC_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
          DO deriv_idx=3,8
            CALL OC_BoundaryConditions_AddNode(BoundaryConditions,DependentField,OC_FIELD_U_VARIABLE_TYPE,1,deriv_idx, &
              & NodeNumber, &
              & component_idx,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,Err)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL OC_Problem_Solve(Problem,Err)

  !Output solution
  CALL OC_Fields_Initialise(Fields,Err)
  CALL OC_Fields_Create(Region,Fields,Err)
  CALL OC_Fields_NodesExport(Fields,"Cantilever","FORTRAN",Err)
  CALL OC_Fields_ElementsExport(Fields,"Cantilever","FORTRAN",Err)
  CALL OC_Fields_Finalise(Fields,Err)

  CALL OC_Context_Destroy(context,Err)
  CALL OC_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HandleError(ERROR_STRING)
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP
  END SUBROUTINE HandleError

END PROGRAM CantileverExample


