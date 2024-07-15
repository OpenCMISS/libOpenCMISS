!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve an Analytic Helmholtz equation using OpenCMISS calls.
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

!> \example ClassicalField/Helmholtz/AnalyticHelmholtz/src/AnalyticHelmholtzExample.F90
!! Example illustrating the use of OpenCMISS to solve the Helmholtz problem and check with its Analytic Solution.
!! 
!! \htmlinclude ClassicalField/Helmholtz/AnalyticHelmholtz/history.html
!< 

!> Main program
PROGRAM ANALYTICHELMHOLTZEXAMPLE
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

  REAL(OC_RP), PARAMETER :: ORIGIN(2)=[-3.141592653579_OC_RP/2, -3.141592653579_OC_RP/2]
  REAL(OC_RP), PARAMETER :: HEIGHT=2.0_OC_RP
  REAL(OC_RP), PARAMETER :: WIDTH=2.0_OC_RP
  REAL(OC_RP), PARAMETER :: LENGTH=2.0_OC_RP
  REAL(OC_RP), PARAMETER :: k=1.0_OC_RP

  !Program types

  !Program variables

  TYPE(OC_ContextType) :: context
  TYPE(OC_RegionType) :: worldRegion
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables
  INTEGER(OC_Intg) :: Err

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
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(1_OC_Intg,context,Err)
  
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  CALL ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(2,10,2)
  CALL ANALYTICHELMHOLTZ_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(2,10,2)
  CALL ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_EXPORT(2,6,0)

  !Destroy the context
  CALL OC_Context_Destroy(context,Err)
  !Finalise OpenCMISS
  CALL OC_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

CONTAINS

  !
  !================================================================================================================================
  !  
    !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_EXPORT(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS)

    !Argument variables
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<number of elements in x direction
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<number of elements in y direction
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<number of elements in z direction
    !Local Variables
    TYPE(OC_FieldType) :: FIELD

    CALL ANALYTICHELMHOLTZ_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,1, &
      & FIELD)

    CALL OC_AnalyticAnalysis_Output(FIELD,"AnalyticHelmholtzBilinear",Err)
    
    CALL ANALYTICHELMHOLTZ_GENERIC_CLEAN(1,1,1,1,1)

  END SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_EXPORT
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(OC_RP) :: VALUE
    REAL(OC_RP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    
    CALL ANALYTICHELMHOLTZ_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,1,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)

    CALL TEST_FRAMEWORK_ASSERT_EQUALS(2.0_OC_RP,VALUE,0.5_OC_RP,ERR)
    
    WRITE(*,'(A)') "Analytic Helmholtz Example Testcase1 - bilinear lagrange is successfully completed."
    
  END SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(OC_RP) :: VALUE
    REAL(OC_RP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)

    CALL ANALYTICHELMHOLTZ_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,3,X_VALUES,Y_VALUES)
    
   CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)
   CALL TEST_FRAMEWORK_ASSERT_EQUALS(4.0_OC_RP,VALUE,1.0_OC_RP,Err)
   IF (Err/=0) THEN
     WRITE(*,'(A,F3.5)') "Analytic Helmholtz Example Testcase2 - bicubic Hermite failure: Convergence should be around 4.0" &
       & //", but it was ", VALUE
   ENDIF
   WRITE(*,'(A)') "Analytic Helmholtz Example Testcase2 - bicubic Hermite is successfully completed."

  END SUBROUTINE ANALYTICHELMHOLTZ_TESTCASE_BICUBIC_HERMITE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
    & NUMBER_OF_ELEMENTS_XI_INTERVAL,INTERPOLATION_SPECIFICATIONS,X_VALUES,Y_VALUES)
  
    !Argument variables 
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: INTERPOLATION_SPECIFICATIONS !<interpolation specifications
    REAL(OC_RP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    !Local Variables
    REAL(OC_RP) :: VALUE
    
    INTEGER(OC_Intg) :: i
    TYPE(OC_FieldType) :: FIELD
    
    ALLOCATE(X_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)
    ALLOCATE(Y_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)

    DO i = NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL
      
      CALL ANALYTICHELMHOLTZ_GENERIC(i,i,0,INTERPOLATION_SPECIFICATIONS,FIELD)
      CALL OC_AnalyticAnalysis_AbsoluteErrorGetNode(FIELD,1,1,1,(i+1)**2/2+1,1,VALUE,Err)

      Y_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(VALUE)
      X_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(HEIGHT/i)
      CALL ANALYTICHELMHOLTZ_GENERIC_CLEAN(1,1,1,1,1)
   
    ENDDO
  END SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CONVERGENCE
  
  
  !
  !================================================================================================================================
  !   
    
  SUBROUTINE ANALYTICHELMHOLTZ_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & INTERPOLATION_SPECIFICATIONS,DEPENDENT_FIELD)
    !Argument variables
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<number of elements on x axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<number of elements on y axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<number of elements on z axis
    INTEGER(OC_Intg), INTENT(IN) :: INTERPOLATION_SPECIFICATIONS !<the interpolation specifications
    TYPE(OC_FieldType) :: DEPENDENT_FIELD
    !Local Variables
    INTEGER(OC_Intg) :: NUMBER_OF_DOMAINS
    INTEGER(OC_Intg) :: MPI_IERROR

    INTEGER(OC_Intg) :: AnalyticFunction
    INTEGER(OC_Intg) :: decompositionIndex,EquationsSetIndex

    INTEGER(OC_Intg), PARAMETER :: CoordinateSystemUserNumber=1
    INTEGER(OC_Intg), PARAMETER :: RegionUserNumber=1
    INTEGER(OC_Intg), PARAMETER :: BasisUserNumber=1
    INTEGER(OC_Intg), PARAMETER :: GeneratedMeshUserNumber=1
    INTEGER(OC_Intg), PARAMETER :: MeshUserNumber=1
    INTEGER(OC_Intg), PARAMETER :: DecompositionUserNumber=1
    INTEGER(OC_Intg), PARAMETER :: DecomposerUserNumber=1
    INTEGER(OC_Intg), PARAMETER :: GeometricFieldUserNumber=1
    INTEGER(OC_Intg), PARAMETER :: EquationsSetFieldUserNumber=2
    INTEGER(OC_Intg), PARAMETER :: MaterialsFieldUserNumber=3
    INTEGER(OC_Intg), PARAMETER :: DependentFieldUserNumber=4
    INTEGER(OC_Intg), PARAMETER :: AnalyticFieldUserNumber=5
    INTEGER(OC_Intg), PARAMETER :: EquationsSetUserNumber=1
    INTEGER(OC_Intg), PARAMETER :: ProblemUserNumber=1

    TYPE(OC_BasisType) :: Basis
    TYPE(OC_BoundaryConditionsType) :: BoundaryConditions
    TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
    TYPE(OC_CoordinateSystemType) :: CoordinateSystem
    TYPE(OC_DecompositionType) :: Decomposition
    TYPE(OC_DecomposerType) :: Decomposer
    TYPE(OC_EquationsType) :: Equations
    TYPE(OC_EquationsSetType) :: EquationsSet
    TYPE(OC_FieldType) :: AnalyticField,GeometricField,EquationsSetField,MaterialsField
    TYPE(OC_GeneratedMeshType) :: GeneratedMesh
    TYPE(OC_MeshType) :: Mesh
    TYPE(OC_ProblemType) :: Problem
    TYPE(OC_RegionType) :: Region
    TYPE(OC_SolverType) :: Solver
    TYPE(OC_SolverEquationsType) :: SolverEquations
    TYPE(OC_WorkGroupType) :: worldWorkGroup
   
    !Get the computation nodes information
    CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
    CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
    
    CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
    CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
    
    !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computation nodes
    CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(INTERPOLATION_SPECIFICATIONS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

    !Start the creation of a new RC coordinate system
    CALL OC_CoordinateSystem_Initialise(CoordinateSystem,Err)
    CALL OC_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,context,CoordinateSystem,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the coordinate system to be 2D
      CALL OC_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
    ELSE
      !Set the coordinate system to be 3D
      CALL OC_CoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
    ENDIF
    !Finish the creation of the coordinate system
    CALL OC_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

    !Start the creation of the region
    CALL OC_Region_Initialise(Region,Err)
    CALL OC_Region_CreateStart(RegionUserNumber,worldRegion,Region,Err)
    !Set the regions coordinate system to the 2D RC coordinate system that we have created
    CALL OC_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
    !Finish the creation of the region
    CALL OC_Region_CreateFinish(Region,Err)
  
    !Start the creation of a basis (default is trilinear lagrange)
    CALL OC_Basis_Initialise(Basis,Err)
    CALL OC_Basis_CreateStart(BasisUserNumber,context,Basis,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the basis to be a bilinear basis
      CALL OC_Basis_NumberOfXiSet(Basis,2,Err)
      CALL OC_Basis_InterpolationXiSet(Basis,[INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS],Err)
    ELSE
      !Set the basis to be a trilinear basis
      CALL OC_Basis_NumberOfXiSet(Basis,3,Err)
      CALL OC_Basis_InterpolationXiSet(Basis,[INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS, &
        & INTERPOLATION_SPECIFICATIONS],Err)
    ENDIF
    !Finish the creation of the basis
    CALL OC_Basis_CreateFinish(Basis,Err)

    !Start the creation of a generated mesh in the region
    CALL OC_GeneratedMesh_Initialise(GeneratedMesh,Err)
    CALL OC_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    !Set up a regular mesh
    CALL OC_GeneratedMesh_TypeSet(GeneratedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
    CALL OC_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)
    !Define the mesh on the region
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
      CALL OC_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS], &
        & Err)
      CALL OC_GeneratedMesh_OriginSet(GeneratedMesh,ORIGIN,Err)
    ELSE
      CALL OC_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
      CALL OC_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS, &
        & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS],Err)
    ENDIF
    !Finish the creation of a generated mesh in the region
    CALL OC_Mesh_Initialise(Mesh,Err)
    CALL OC_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
    
    !Create a decomposition
    CALL OC_Decomposition_Initialise(Decomposition,Err)
    CALL OC_Decomposition_CreateStart(1,Mesh,Decomposition,Err)
    CALL OC_Decomposition_CreateFinish(Decomposition,Err)

    !Decompose
    CALL OC_Decomposer_Initialise(decomposer,err)
    CALL OC_Decomposer_CreateStart(1,region,worldWorkGroup,decomposer,err)
    !Add in the decomposition
    CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
    !Finish the decomposer
    CALL OC_Decomposer_CreateFinish(decomposer,err)
  
    !Start to create a default (geometric) field on the region
    CALL OC_Field_Initialise(GeometricField,Err)
    CALL OC_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
    !Set the decomposition to use
    CALL OC_Field_DecompositionSet(GeometricField,Decomposition,Err)
    !Set the domain to be used by the field components
    !NB these are needed now as the default mesh component number is 1
    CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,1,1,Err)
    CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,2,1,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,3,1,Err)
    ENDIF
    !Finish creating the field
    CALL OC_Field_CreateFinish(GeometricField,Err)

    !Update the geometric field parameters
    CALL OC_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

    !Create the equations set
    CALL OC_EquationsSet_Initialise(EquationsSet,Err)
    CALL OC_Field_Initialise(EquationsSetField,Err)
    !Set the equations set to be a standard Helmholtz problem
    CALL OC_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
      & OC_EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE,OC_EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE],EquationsSetFieldUserNumber, &
      & EquationsSetField,EquationsSet,Err)
    !Finish creating the equations set
    CALL OC_EquationsSet_CreateFinish(EquationsSet,Err)
  
    !Create the equations set dependent field variables
    CALL OC_Field_Initialise(DEPENDENT_FIELD,Err)
    CALL OC_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DEPENDENT_FIELD,Err)
    !Finish the equations set dependent field variables
    CALL OC_EquationsSet_DependentCreateFinish(EquationsSet,Err)

    !Create the equations set material field variables
    CALL OC_Field_Initialise(MaterialsField,Err)
    CALL OC_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
    CALL OC_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)
    !Set wave number, k
    CALL OC_Field_ComponentValuesInitialise(MaterialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,k,Err)

    !Create the equations set analytic field variables
    AnalyticFunction=OC_EQUATIONS_SET_HELMHOLTZ_EQUATION_TWO_DIM_1
    CALL OC_Field_Initialise(AnalyticField,Err)
    CALL OC_EquationsSet_AnalyticCreateStart(EquationsSet,AnalyticFunction,AnalyticFieldUserNumber,AnalyticField,Err)
    !Finish the equations set analtyic field variables
    CALL OC_EquationsSet_AnalyticCreateFinish(EquationsSet,Err)

    !Create the equations set equations
    CALL OC_Equations_Initialise(Equations,Err)
    CALL OC_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
    !Set the equations matrices sparsity type
    CALL OC_Equations_SparsityTypeSet(Equations,OC_EQUATIONS_SPARSE_MATRICES,Err)
    CALL OC_EquationsSet_EquationsCreateFinish(EquationsSet,Err)
  
    !Create the problem
    CALL OC_Problem_Initialise(Problem,Err)
    CALL OC_Problem_CreateStart(ProblemUserNumber,context,[OC_PROBLEM_CLASSICAL_FIELD_CLASS, &
      & OC_PROBLEM_HELMHOLTZ_EQUATION_TYPE,OC_PROBLEM_STANDARD_HELMHOLTZ_SUBTYPE],Problem,Err)
    !Finish creating the problem
    CALL OC_Problem_CreateFinish(Problem,Err)

    !Create the problem control loop
    CALL OC_Problem_ControlLoopCreateStart(Problem,Err)
    !Finish creating the problem control
    CALL OC_Problem_ControlLoopCreateFinish(Problem,Err)

    !Start the creation of the problem solvers
    CALL OC_Solver_Initialise(Solver,Err)
    CALL OC_Problem_SolversCreateStart(Problem,Err)
    CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
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

  END SUBROUTINE ANALYTICHELMHOLTZ_GENERIC

  SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber,GeneratedMeshUserNumber, &
    & ProblemUserNumber)

    !Argument variables    
    INTEGER(OC_Intg), INTENT(IN) :: CoordinateSystemUserNumber
    INTEGER(OC_Intg), INTENT(IN) :: RegionUserNumber
    INTEGER(OC_Intg), INTENT(IN) :: BasisUserNumber
    INTEGER(OC_Intg), INTENT(IN) :: GeneratedMeshUserNumber
    INTEGER(OC_Intg), INTENT(IN) :: ProblemUserNumber

    INTEGER(OC_Intg) :: contextUserNumber

    CALL OC_Context_UserNumberGet(context,contextUserNumber,err)

    CALL OC_Problem_Destroy(contextUserNumber,ProblemUserNumber,Err)
    CALL OC_GeneratedMesh_Destroy(contextUserNumber,RegionUserNumber,GeneratedMeshUserNumber,Err)
    CALL OC_Basis_Destroy(contextUserNumber,BasisUserNumber,Err)
    CALL OC_Region_Destroy(contextUserNumber,RegionUserNumber,Err)
    CALL OC_CoordinateSystem_Destroy(contextUserNumber,CoordinateSystemUserNumber,Err)

  END SUBROUTINE ANALYTICHELMHOLTZ_GENERIC_CLEAN



END PROGRAM ANALYTICHELMHOLTZEXAMPLE
