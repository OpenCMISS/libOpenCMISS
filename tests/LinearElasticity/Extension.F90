!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve an Analytic Laplace equation using OpenCMISS calls.
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

!> \example ClassicalField/Laplace/ANALYTIC_LINEAR_ELASTICITY/src/ANALYTIC_LINEAR_ELASTICITYExample.F90
!! Example illustrating the use of OpenCMISS to solve the Laplace problem and check with its Analytic Solution.
!! 
!! \htmlinclude ClassicalField/Laplace/ANALYTIC_LINEAR_ELASTICITY/history.html
!< 

!> Main program
PROGRAM AnalyticLinearElasticityExample
#ifndef NOMPIMOD
  USE MPI
#endif
  USE OpenCMISS
  USE TEST_FRAMEWORK_ROUTINES

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  REAL(OCRP), PARAMETER :: ORIGIN(3)=[0.0_OCRP,0.0_OCRP,0.0_OCRP]
  REAL(OCRP), PARAMETER :: LENGTH=20.0_OCRP
  REAL(OCRP), PARAMETER :: WIDTH=20.0_OCRP
  REAL(OCRP), PARAMETER :: HEIGHT=5.0_OCRP

  INTEGER(OCIntg), PARAMETER :: NumberOfDomains=1

  INTEGER(OCIntg), PARAMETER :: ContextUserNumber=1
  INTEGER(OCIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(OCIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(OCIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(OCIntg), PARAMETER :: GeneratedMeshUserNumber = 1
  INTEGER(OCIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(OCIntg), PARAMETER :: DecomposerUserNumber=1

  INTEGER(OCIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(OCIntg), PARAMETER :: FieldGeometryNumberOfVariables=1

  INTEGER(OCIntg), PARAMETER :: FieldDependentUserNumber=2
  INTEGER(OCIntg), PARAMETER :: FieldDependentNumberOfVariables=2

  INTEGER(OCIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(OCIntg), PARAMETER :: FieldMaterialNumberOfVariables=1

  INTEGER(OCIntg), PARAMETER :: FieldAnalyticUserNumber=4

  INTEGER(OCIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(OCIntg), PARAMETER :: EquationsSetFieldUserNumber=5
  INTEGER(OCIntg), PARAMETER :: ProblemUserNumber=1

  !Program types

  TYPE(OC_ContextType) :: context
  TYPE(OC_RegionType) :: WorldRegion

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
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(ContextUserNumber,context,Err)  
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  !CALL OC_DiagnosticsSetOn(OC_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_FINITE_ELEMENT_CALCULATE"],Err)

  CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(1,0,0,"LinearLagrange")
  CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(1,1,0,"BiLinearLagrange")
  CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(1,1,1,"TriLinearLagrange")
  !CALL ANALYTIC_LINEAR_ELASTICITY_TESTCASE_QUADRATIC_LAGRANGE_EXPORT(1,0,0,"QuadraticLagrange")

  CALL OC_Context_Destroy(context,Err)
  CALL OC_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

CONTAINS


  !
  !================================================================================================================================
  !  
    !>Check if the convergence of linear langrange interpolation is expected.
  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT(NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements,OutputFilename)

    !Argument variables
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalXElements !<initial number of elements per axis
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalYElements !<final number of elements per axis
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalZElements !<increment interval number of elements per axis
    CHARACTER(LEN=*), INTENT(IN) :: OutputFilename !<The Error condition string
    !Local Variables
    TYPE(OC_FieldType) :: DependentField

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements, &
      & OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION,DependentField)

    CALL OC_AnalyticAnalysis_Output(DependentField,OutputFilename,Err)
    
    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
      & GeneratedMeshUserNumber,ProblemUserNumber)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_LINEAR_LAGRANGE_EXPORT

  !
  !================================================================================================================================
  !  
    !>Check if the convergence of quadratic langrange interpolation is expected.
  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_QUADRATIC_LAGRANGE_EXPORT(NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements,OutputFilename)

    !Argument variables
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalXElements !<initial number of elements per axis
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalYElements !<final number of elements per axis
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalZElements !<increment interval number of elements per axis
    CHARACTER(LEN=*), INTENT(IN) :: OutputFilename !<The Error condition string
    !Local Variables
    TYPE(OC_FieldType) :: DependentField

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements, &
      & OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,DependentField)

    CALL OC_AnalyticAnalysis_Output(DependentField,OutputFilename,Err)
    
    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
      & GeneratedMeshUserNumber,ProblemUserNumber)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_QUADRATIC_LAGRANGE_EXPORT

  !
  !================================================================================================================================
  !  
    !>Check if the convergence of cubic langrange interpolation is expected.
  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_CUBIC_LAGRANGE_EXPORT(NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements,OutputFilename)

    !Argument variables
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalXElements !<initial number of elements per axis
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalYElements !<final number of elements per axis
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalZElements !<increment interval number of elements per axis
    CHARACTER(LEN=*), INTENT(IN) :: OutputFilename !<The Error condition string
    !Local Variables
    TYPE(OC_FieldType) :: DependentField

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements, &
      & OC_BASIS_CUBIC_LAGRANGE_INTERPOLATION,DependentField)

    CALL OC_AnalyticAnalysis_Output(DependentField,OutputFilename,Err)
    
    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
      & GeneratedMeshUserNumber,ProblemUserNumber)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_TESTCASE_CUBIC_LAGRANGE_EXPORT

  !
  !================================================================================================================================
  !   
    
  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC(NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements, &
    & InterpolationSpecifications,DependentField)
    !Argument variables 
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalXElements !<number of elements on x axis
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalYElements !<number of elements on y axis
    INTEGER(OCIntg), INTENT(IN) :: NumberGlobalZElements !<number of elements on z axis
    INTEGER(OCIntg), INTENT(IN) :: InterpolationSpecifications !<the interpolation specifications
    TYPE(OC_FieldType) :: DependentField

    !Program variables
    REAL(OCRP) :: MeshDimensions(3),MaterialParameters(6)
    INTEGER(OCIntg) :: AnalyticFunction,Interpolation(3),NumberOfGaussPoints(3),EquationSetSubtype
    INTEGER(OCIntg) :: FieldGeometryNumberOfComponents,FieldDependentNumberOfComponents,NumberOfElements(3)
    INTEGER(OCIntg) :: MPI_IERROR
    INTEGER(OCIntg) :: decompositionIndex,EquationsSetIndex,FieldComponentIndex,FieldMaterialNumberOfComponents,NumberOfXi
    INTEGER(OCIntg) :: NumberOfComputationNodes,ComputationNodeNumber

    !CMISS variables

    TYPE(OC_BasisType) :: Basis
    TYPE(OC_BoundaryConditionsType) :: BoundaryConditions
    TYPE(OC_ComputationEnvironmentType) :: ComputationEnvironment
    TYPE(OC_CoordinateSystemType) :: CoordinateSystem
    TYPE(OC_GeneratedMeshType) :: GeneratedMesh
    TYPE(OC_DecompositionType) :: Decomposition
    TYPE(OC_DecomposerType) :: Decomposer
    TYPE(OC_EquationsType) :: Equations
    TYPE(OC_EquationsSetType) :: EquationsSet
    TYPE(OC_FieldType) :: AnalyticField,EquationsSetField,GeometricField,MaterialField
    TYPE(OC_MeshType) :: Mesh
    TYPE(OC_ProblemType) :: Problem
    TYPE(OC_RegionType) :: Region
    TYPE(OC_SolverType) :: Solver
    TYPE(OC_SolverEquationsType) :: SolverEquations
    TYPE(OC_WorkGroupType) :: worldWorkGroup

    IF((NumberGlobalYElements == 0) .AND. (NumberGlobalZElements == 0)) THEN
      NumberOfXi = 1
      EquationSetSubtype = OC_EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE
      AnalyticFunction=OC_EQUATIONS_SET_LINEAR_ELASTICITY_ONE_DIM_1
      !Prescribe material properties Area,E1
      FieldMaterialNumberOfComponents = 2 !Young's Modulus & Poisson's Ratio
      MaterialParameters = [WIDTH*HEIGHT,10.0E3_OCRP,0.0_OCRP,0.0_OCRP,0.0_OCRP,0.0_OCRP]
    ELSEIF (NumberGlobalZElements == 0) THEN
      NumberOfXi = 2
      EquationSetSubtype = OC_EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE
      AnalyticFunction=OC_EQUATIONS_SET_LINEAR_ELASTICITY_TWO_DIM_1
      !Prescribe material properties h,E1,v12
      FieldMaterialNumberOfComponents = 3 !Young's Modulus & Poisson's Ratio
      MaterialParameters = [HEIGHT,10.0E3_OCRP,0.3_OCRP,0.0_OCRP,0.0_OCRP,0.0_OCRP]
    ELSE
      NumberOfXi = 3
      EquationSetSubtype = OC_EQUATIONS_SET_THREE_DIMENSIONAL_ORTHOTROPIC_SUBTYPE
      AnalyticFunction=OC_EQUATIONS_SET_LINEAR_ELASTICITY_THREE_DIM_1
      !Prescribe material properties E1,E2,E3 & v13,v23,v12
      FieldMaterialNumberOfComponents = 6 !Young's Modulus & Poisson's Ratio
      MaterialParameters = [10.0E3_OCRP,10.0E3_OCRP,10.0E3_OCRP,0.3_OCRP,0.3_OCRP,0.3_OCRP]
    ENDIF
    Interpolation = [InterpolationSpecifications,InterpolationSpecifications,InterpolationSpecifications]
    NumberOfElements = [NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements]
    MeshDimensions = [LENGTH,WIDTH,HEIGHT]
    NumberOfGaussPoints = [4,4,4]
    FieldGeometryNumberOfComponents=NumberOfXi
    FieldDependentNumberOfComponents=NumberOfXi

    !Get the number of computation nodes and this computation node number
    CALL OC_ComputationEnvironment_Initialise(ComputationEnvironment,Err)
    CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
    CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
    CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
    CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationNodes,err)
    CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationNodeNumber,err)

    !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computation nodes
    CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

    !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
    CALL OC_CoordinateSystem_Initialise(CoordinateSystem,Err)
    CALL OC_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,context,CoordinateSystem,Err)
    CALL OC_CoordinateSystem_TypeSet(CoordinateSystem,OC_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
    CALL OC_CoordinateSystem_DimensionSet(CoordinateSystem,NumberOfXi,Err)
    CALL OC_CoordinateSystem_OriginSet(CoordinateSystem,ORIGIN,Err)
    CALL OC_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

    !Create a region and assign the CS to the region
    CALL OC_Region_Initialise(Region,Err)
    CALL OC_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
    CALL OC_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
    CALL OC_Region_CreateFinish(Region,Err)

    CALL OC_Basis_Initialise(Basis,Err)
    CALL OC_Basis_CreateStart(BasisUserNumber,context,Basis,Err)
    CALL OC_Basis_TypeSet(Basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CALL OC_Basis_NumberOfXiSet(Basis,NumberOfXi,Err)
    CALL OC_Basis_InterpolationXiSet(Basis,Interpolation(1:NumberOfXi),Err)
    CALL OC_Basis_QuadratureNumberOfGaussXiSet(Basis,NumberOfGaussPoints(1:NumberOfXi),Err)
    CALL OC_Basis_CreateFinish(Basis,Err)

    !Start the creation of a generated Mesh in the Region
    CALL OC_GeneratedMesh_Initialise(GeneratedMesh,Err)
    CALL OC_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    CALL OC_GeneratedMesh_TypeSet(GeneratedMesh,1,Err)
    CALL OC_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)

    !Define the Mesh on the Region
    CALL OC_GeneratedMesh_OriginSet(GeneratedMesh,ORIGIN(1:NumberOfXi),Err)
    CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,MeshDimensions(1:NumberOfXi),Err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,NumberOfElements(1:NumberOfXi),Err)
    CALL OC_Mesh_Initialise(Mesh,Err)
    CALL OC_GeneratedMesh_CreateFinish(GeneratedMesh,1,Mesh,Err)

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
    CALL OC_Field_TypeSet(GeometricField,OC_FIELD_GEOMETRIC_TYPE,Err)  
    CALL OC_Field_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
    CALL OC_Field_NumberOfComponentsSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
    DO FieldComponentIndex=1,FieldGeometryNumberOfComponents
      CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL OC_Field_CreateFinish(GeometricField,Err)

    !Update the geometric field parameters
    CALL OC_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

    !Create a dependent field with two variables and three components
    CALL OC_Field_Initialise(DependentField,Err)
    CALL OC_Field_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
    CALL OC_Field_TypeSet(DependentField,OC_FIELD_GENERAL_TYPE,Err)  
    CALL OC_Field_DecompositionSet(DependentField,Decomposition,Err)
    CALL OC_Field_GeometricFieldSet(DependentField,GeometricField,Err) 
    CALL OC_Field_DependentTypeSet(DependentField,OC_FIELD_DEPENDENT_TYPE,Err) 
    CALL OC_Field_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)
    CALL OC_Field_NumberOfComponentsSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
    CALL OC_Field_NumberOfComponentsSet(DependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentNumberOfComponents,Err)
    DO FieldComponentIndex=1,FieldDependentNumberOfComponents
      CALL OC_Field_ComponentMeshComponentSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,FieldComponentIndex,1,Err)
      CALL OC_Field_ComponentMeshComponentSet(DependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL OC_Field_CreateFinish(DependentField,Err)

    !Create a material field and attach it to the geometric field
    CALL OC_Field_Initialise(MaterialField,Err)
    CALL OC_Field_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
    CALL OC_Field_TypeSet(MaterialField,OC_FIELD_MATERIAL_TYPE,Err)
    CALL OC_Field_DecompositionSet(MaterialField,Decomposition,Err)
    CALL OC_Field_GeometricFieldSet(MaterialField,GeometricField,Err)
    CALL OC_Field_NumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
    CALL OC_Field_NumberOfComponentsSet(MaterialField,OC_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)
    DO FieldComponentIndex=1,FieldMaterialNumberOfComponents
      CALL OC_Field_ComponentMeshComponentSet(MaterialField,OC_FIELD_U_VARIABLE_TYPE,FieldComponentIndex,1,Err)
    ENDDO !FieldComponentIndex
    CALL OC_Field_CreateFinish(MaterialField,Err)

    !Set isotropic elasticity material parameters - Young's Modulus & Poisson's Ratio
    DO FieldComponentIndex=1,FieldMaterialNumberOfComponents
      CALL OC_Field_ComponentValuesInitialise(MaterialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
        & FieldComponentIndex, &
        & MaterialParameters(FieldComponentIndex),Err)
    ENDDO !FieldComponentIndex

    !Create a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
    CALL OC_EquationsSet_Initialise(EquationsSet,Err)
    CALL OC_Field_Initialise(EquationsSetField,Err)
    CALL OC_EquationsSet_CreateStart(EquationSetUserNumber,Region,GeometricField,[OC_EQUATIONS_SET_ELASTICITY_CLASS, &
      & OC_EQUATIONS_SET_LINEAR_ELASTICITY_TYPE,EquationSetSubtype],EquationsSetFieldUserNumber,EquationsSetField, &
      & EquationsSet,Err)
    
    CALL OC_EquationsSet_CreateFinish(EquationsSet,Err)

    CALL OC_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err) 
    CALL OC_EquationsSet_DependentCreateFinish(EquationsSet,Err)

    CALL OC_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)  
    CALL OC_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

    !Create the Equations set analtyic field variables
    CALL OC_Field_Initialise(AnalyticField,Err)
    CALL OC_EquationsSet_AnalyticCreateStart(EquationsSet,AnalyticFunction,FieldAnalyticUserNumber,AnalyticField,Err)
    CALL OC_EquationsSet_AnalyticCreateFinish(EquationsSet,Err)

    !Create the equations set equations
    CALL OC_Equations_Initialise(Equations,Err)
    CALL OC_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
    CALL OC_Equations_SparsityTypeSet(EQUATIONS,OC_EQUATIONS_SPARSE_MATRICES,Err)
                                                !OC_EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the equations.
                                                !OC_EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the equations. 
    CALL OC_Equations_OutputTypeSet(EQUATIONS,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
                                              !OC_EQUATIONS_NO_OUTPUT !<No output from the equations.
                                              !OC_EQUATIONS_TIMING_OUTPUT !<Timing information output.
                                              !OC_EQUATIONS_MATRIX_OUTPUT !<All below and equation matrices output.
                                              !OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT !<All below and Element matrices output.
    CALL OC_EquationsSet_EquationsCreateFinish(EquationsSet,Err)
    
    !Define the problem
    CALL OC_Problem_Initialise(Problem,Err)
    CALL OC_Problem_CreateStart(ProblemUserNumber,context,[OC_PROBLEM_ELASTICITY_CLASS,OC_PROBLEM_LINEAR_ELASTICITY_TYPE, &
      & OC_PROBLEM_NO_SUBTYPE],Problem,Err)
    CALL OC_Problem_CreateFinish(Problem,Err)

    !Create the problem control loop
    CALL OC_Problem_ControlLoopCreateStart(Problem,Err)
    CALL OC_Problem_ControlLoopCreateFinish(Problem,Err)

    !Start the creation of the Problem Solvers
    !Create the problem Solvers
    CALL OC_Solver_Initialise(Solver,Err)
    CALL OC_Problem_SolversCreateStart(Problem,Err)
    CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
    CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_MATRIX_OUTPUT,Err)
                                        !OC_SOLVER_NO_OUTPUT !<No output from the Solver routines. \see OPENCMISS_SolverOutputTypes,OPENCMISS
                                        !OC_SOLVER_PROGRESS_OUTPUT !<Progress output from Solver routines.
                                        !OC_SOLVER_TIMING_OUTPUT !<Timing output from the Solver routines plus below.
                                        !OC_SOLVER_SOLVER_OUTPUT !<Solver specific output from the Solver routines plus below.
                                        !OC_SOLVER_MATRIX_OUTPUT !<Solver matrices output from the Solver routines plus below.
    CALL OC_Solver_LibraryTypeSet(Solver,OC_SOLVER_PETSC_LIBRARY,Err)
                                          !OC_SOLVER_OC_LIBRARY     !<CMISS (internal) Solver library.
                                          !OC_SOLVER_PETSC_LIBRARY     !<PETSc Solver library.
                                          !OC_SOLVER_MUMPS_LIBRARY     !<MUMPS Solver library.
                                          !OC_SOLVER_SUPERLU_LIBRARY   !<SuperLU Solver library.
                                          !OC_SOLVER_SPOOLES_LIBRARY !<SPOOLES Solver library.
                                          !OC_SOLVER_UMFPACK_LIBRARY   !<UMFPACK Solver library.
                                          !OC_SOLVER_LUSOL_LIBRARY     !<LUSOL Solver library.
                                          !OC_SOLVER_ESSL_LIBRARY      !<ESSL Solver library.
                                          !OC_SOLVER_LAPACK_LIBRARY    !<LAPACK Solver library.
                                          !OC_SOLVER_TAO_LIBRARY       !<TAO Solver library.
                                          !OC_SOLVER_HYPRE_LIBRARY     !<Hypre Solver library.
    CALL OC_Solver_LinearTypeSet(Solver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
                                        !OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE    !<Direct linear Solver type.
                                        !OC_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE !<Iterative linear Solver type.
    CALL OC_Problem_SolversCreateFinish(Problem,Err)

    !Create the problem Solver equations
    CALL OC_Solver_Initialise(Solver,Err)
    CALL OC_SolverEquations_Initialise(SolverEquations,Err)
    CALL OC_Problem_SolverEquationsCreateStart(Problem,Err)   
    CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
    CALL OC_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
    CALL OC_SolverEquations_SparsityTypeSet(SolverEquations,OC_SOLVER_SPARSE_MATRICES,Err)
                                                            !OC_SOLVER_SPARSE_MATRICES !<Use sparse Solver matrices.
                                                            !OC_SOLVER_FULL_MATRICES !<Use fully populated Solver matrices.
    CALL OC_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
    CALL OC_Problem_SolverEquationsCreateFinish(Problem,Err)

    !Prescribe boundary conditions
    CALL OC_BoundaryConditions_Initialise(BoundaryConditions,Err)
    CALL OC_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
    CALL OC_SolverEquations_BoundaryConditionsAnalytic(SolverEquations,Err)
    CALL OC_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

    !=SOLVE Problem==================================================================================================================
    !Solve the Problem
    CALL OC_Problem_Solve(Problem,Err)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC

  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
    & GeneratedMeshUserNumber,ProblemUserNumber)

    !Argument variables
    INTEGER(OCIntg), INTENT(IN) :: CoordinateSystemUserNumber
    INTEGER(OCIntg), INTENT(IN) :: RegionUserNumber
    INTEGER(OCIntg), INTENT(IN) :: BasisUserNumber
    INTEGER(OCIntg), INTENT(IN) :: GeneratedMeshUserNumber
    INTEGER(OCIntg), INTENT(IN) :: ProblemUserNumber

    INTEGER(OCIntg) :: contextUserNumber

    CALL OC_Context_UserNumberGet(context,contextUserNumber,err)
    CALL OC_Problem_Destroy(contextUserNumber,ProblemUserNumber,Err)
    CALL OC_GeneratedMesh_Destroy(contextUserNumber,RegionUserNumber,GeneratedMeshUserNumber,Err)
    CALL OC_Basis_Destroy(contextUserNumber,BasisUserNumber,Err)
    CALL OC_Region_Destroy(contextUserNumber,RegionUserNumber,Err)
    CALL OC_CoordinateSystem_Destroy(contextUserNumber,CoordinateSystemUserNumber,Err)

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN

END PROGRAM AnalyticLinearElasticityExample 

