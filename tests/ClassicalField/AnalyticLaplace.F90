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

!> \example ClassicalField/Laplace/AnalyticLaplace/src/AnalyticLaplaceExample.F90
!! Example illustrating the use of OpenCMISS to solve the Laplace problem and check with its Analytic Solution.
!! 
!! \htmlinclude ClassicalField/Laplace/AnalyticLaplace/history.html
!< 

!> Main program
PROGRAM AnalyticLaplaceExample
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

  REAL(OC_RP), PARAMETER :: ORIGIN(2)=[-3.141592653579_OC_RP/2.0_OC_RP, -3.141592653579_OC_RP/2.0_OC_RP]
  REAL(OC_RP), PARAMETER :: HEIGHT=2.0_OC_RP
  REAL(OC_RP), PARAMETER :: WIDTH=2.0_OC_RP
  REAL(OC_RP), PARAMETER :: LENGTH=2.0_OC_RP

  !Program types

  !Program variables

  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_RegionType) :: worldRegion
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  INTEGER(OC_Intg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS,INTERPOLATION
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT

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
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(1_OC_Intg,context,Err)

  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  CALL OC_Context_RandomSeedsSet(context,9999,err)
    
  !CALL OC_DiagnosticsSetOn(OC_ALL_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",[""],Err)
  
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION

    SELECT CASE(INTERPOLATION)
    CASE(1)
      CALL ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(2,6,2)
      CALL ANALYTICLAPLACE_TESTCASE_LINEAR_LAGRANGE_EXPORT(2,2,0)
    CASE(4)
      CALL ANALYTICLAPLACE_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(2,10,2)
      CALL ANALYTICLAPLACE_TESTCASE_CUBIC_HERMITE_EXPORT(2,2,0)
    CASE(7)
      CALL ANALYTICLAPLACE_TESTCASE_BILINEAR_SIMPLEX_CONVERGENCE(2,6,2)
      CALL ANALYTICLAPLACE_TESTCASE_LINEAR_SIMPLEX_EXPORT(2,2,0)
    CASE DEFAULT
      CALL HANDLE_ERROR("Invalid interpolation specified.")
    END SELECT
  ELSE
    !Run all tests
    CALL ANALYTICLAPLACE_TESTCASE_BILINEAR_SIMPLEX_CONVERGENCE(2,6,2)
    CALL ANALYTICLAPLACE_TESTCASE_LINEAR_SIMPLEX_EXPORT(2,2,0)
    CALL ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(2,6,2)
    CALL ANALYTICLAPLACE_TESTCASE_LINEAR_LAGRANGE_EXPORT(2,2,0)
    CALL ANALYTICLAPLACE_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(2,10,2)
    CALL ANALYTICLAPLACE_TESTCASE_CUBIC_HERMITE_EXPORT(2,2,0)
  ENDIF

  CALL OC_Context_Destroy(context,Err)
  CALL OC_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

CONTAINS

  !
  !================================================================================================================================
  !  

  !>Export analytic analyis
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_CUBIC_HERMITE_EXPORT(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS)

    !Argument variables
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<initial number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<final number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<increment interval number of elements per axis
    !Local Variables
    INTEGER(OC_Intg) :: contextUserNumber
    TYPE(OC_FieldType) :: FIELD

    CALL OC_Context_UserNumberGet(context,contextUserNumber,err)
    CALL OC_Field_Initialise(FIELD,Err)
    CALL ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,4, &
      & FIELD)

    CALL OC_AnalyticAnalysis_Output(FIELD,"AnalyticLaplaceCubicHermite",Err)
    
    CALL ANALYTICLAPLACE_GENERIC_CLEAN(contextUserNumber,1,1,1,1,1)

  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_CUBIC_HERMITE_EXPORT
  
  !
  !================================================================================================================================
  !
  
  !>Export analytic analyis
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_LINEAR_LAGRANGE_EXPORT(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS)

    !Argument variables
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<initial number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<final number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<increment interval number of elements per axis
    !Local Variables
    INTEGER(OC_Intg) :: contextUserNumber
    TYPE(OC_FieldType) :: FIELD

    CALL OC_Context_UserNumberGet(context,contextUserNumber,err)
    CALL OC_Field_Initialise(FIELD,Err)
    CALL ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,1, &
      & FIELD)

    CALL OC_AnalyticAnalysis_Output(FIELD,"AnalyticLaplaceLinearLagrange",Err)
    
    CALL ANALYTICLAPLACE_GENERIC_CLEAN(contextUserNumber,1,1,1,1,1)

  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_LINEAR_LAGRANGE_EXPORT
  
  !
  !================================================================================================================================
  !
  
  !>Export analytic analyis
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_LINEAR_SIMPLEX_EXPORT(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS)

    !Argument variables
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<initial number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<final number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<increment interval number of elements per axis
    !Local Variables
    INTEGER(OC_Intg) :: contextUserNumber
    TYPE(OC_FieldType) :: FIELD

    CALL OC_Context_UserNumberGet(context,contextUserNumber,err)
    CALL OC_Field_Initialise(FIELD,Err)
    CALL ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,7, &
      & FIELD)

    CALL OC_AnalyticAnalysis_Output(FIELD,"AnalyticLaplaceLinearSimplex",Err)
    
    CALL ANALYTICLAPLACE_GENERIC_CLEAN(contextUserNumber,1,1,1,1,1)

  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_LINEAR_SIMPLEX_EXPORT
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear Simplex interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_SIMPLEX_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(OC_RP) :: VALUE
    REAL(OC_RP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    
    CALL ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,7,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)

    CALL TEST_FRAMEWORK_ASSERT_EQUALS(2.0_OC_RP,VALUE,0.5_OC_RP,ERR)
    
    WRITE(*,'(A)') "Analytic Laplace Example Testcase 1 - bilinear Simplex has successfully completed."
    
  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_SIMPLEX_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(OC_RP) :: VALUE
    REAL(OC_RP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    
    CALL ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,1,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)

    CALL TEST_FRAMEWORK_ASSERT_EQUALS(2.0_OC_RP,VALUE,0.5_OC_RP,ERR)
    
    WRITE(*,'(A)') "Analytic Laplace Example Testcase 2 - bilinear Lagrange has successfully completed."
    
  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(OC_RP) :: VALUE
    REAL(OC_RP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)

    !Note INTERPOLATION_SPECIFICATIONS of 4 is Cubic Hermite
    CALL ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,4,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)
    !This test is superconvergent so look for a slope of 5 rather than 4. Should really test >= 4
    CALL TEST_FRAMEWORK_ASSERT_EQUALS(5.0_OC_RP,VALUE,1.0_OC_RP,Err)
    IF (Err/=0) THEN
      WRITE(*,'(A,F6.3)') "Analytic Laplace Example Testcase 3 - bicubic Hermite failure: Convergence should be around 4.0" &
        & //", but it was ", VALUE
    ENDIF
    WRITE(*,'(A)') "Analytic Laplace Example Testcase 3 - bicubic Hermite has successfully completed."
    
  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_BICUBIC_HERMITE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of the specified interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
    & NUMBER_OF_ELEMENTS_XI_INTERVAL,INTERPOLATION_SPECIFICATIONS,X_VALUES,Y_VALUES)
  
    !Argument variables 
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    INTEGER(OC_Intg), INTENT(IN) :: INTERPOLATION_SPECIFICATIONS !<interpolation specifications
    REAL(OC_RP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    !Local Variables
    REAL(OC_RP) :: VALUE
    
    INTEGER(OC_Intg) :: i,contextUserNumber
    TYPE(OC_FieldType) :: FIELD


    CALL OC_Context_UserNumberGet(context,contextUserNumber,err)
    ALLOCATE(X_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)
    ALLOCATE(Y_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)

    DO i = NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL

      CALL OC_Field_Initialise(FIELD,Err)
      CALL ANALYTICLAPLACE_GENERIC(i,i,0,INTERPOLATION_SPECIFICATIONS,FIELD)
      CALL OC_AnalyticAnalysis_AbsoluteErrorGetNode(FIELD,1,1,1,(i+1)**2/2+1,1,VALUE,Err)

      Y_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(VALUE)
      X_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(HEIGHT/i)
      CALL ANALYTICLAPLACE_GENERIC_CLEAN(contextUserNumber,1,1,1,1,1)
   
    ENDDO
  END SUBROUTINE ANALYTICLAPLACE_GENERIC_CONVERGENCE
  
  
  !
  !================================================================================================================================
  !   
    
  SUBROUTINE ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & INTERPOLATION_SPECIFICATIONS,DEPENDENT_FIELD)
    !Argument variables 
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<number of elements on x axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<number of elements on y axis
    INTEGER(OC_Intg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<number of elements on z axis
    INTEGER(OC_Intg), INTENT(IN) :: INTERPOLATION_SPECIFICATIONS !<the interpolation specifications
    TYPE(OC_FieldType) :: DEPENDENT_FIELD
    !Local Variables
    INTEGER(OC_Intg) :: MPI_IERROR
    INTEGER(OC_Intg) :: ANALYTIC_FUNCTION
    INTEGER(OC_Intg) :: decompositionIndex, EquationsSetIndex

    TYPE(OC_BasisType) :: Basis
    TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
    TYPE(OC_ContextType) :: context
    TYPE(OC_CoordinateSystemType) :: CoordinateSystem
    TYPE(OC_BoundaryConditionsType) :: BoundaryConditions
    TYPE(OC_GeneratedMeshType) :: GENERATED_MESH
    TYPE(OC_MeshType) :: MESH
    TYPE(OC_DecompositionType) :: DECOMPOSITION
    TYPE(OC_DecomposerType) :: decomposer
    TYPE(OC_EquationsType) :: EQUATIONS
    TYPE(OC_EquationsSetType) :: EQUATIONS_SET
    TYPE(OC_FieldType) :: ANALYTIC_FIELD,GEOMETRIC_FIELD,EquationsSetField
    TYPE(OC_ProblemType) :: PROBLEM
    TYPE(OC_RegionType) :: REGION
    TYPE(OC_SolverType) :: SOLVER
    TYPE(OC_SolverEquationsType) :: SOLVER_EQUATIONS
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
    CALL OC_CoordinateSystem_CreateStart(1,context,CoordinateSystem,Err)
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
    CALL OC_Region_Initialise(REGION,Err)
    CALL OC_Region_CreateStart(1,worldRegion,REGION,Err)
    !Set the regions coordinate system to the 2D RC coordinate system that we have created
    CALL OC_Region_CoordinateSystemSet(REGION,CoordinateSystem,Err)
    !Finish the creation of the region
    CALL OC_Region_CreateFinish(REGION,Err)
  
    !Start the creation of a basis (default is trilinear lagrange)
    CALL OC_Basis_Initialise(Basis,Err)
    CALL OC_Basis_CreateStart(1,context,Basis,Err)
    SELECT CASE(INTERPOLATION_SPECIFICATIONS)
    CASE(OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION,OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
      & OC_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
      & OC_BASIS_CUBIC_HERMITE_INTERPOLATION)
      !Do nothing
    CASE(OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION,OC_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
      & OC_BASIS_CUBIC_SIMPLEX_INTERPOLATION)
      CALL OC_Basis_TypeSet(Basis,OC_BASIS_SIMPLEX_TYPE,Err)
    CASE DEFAULT
      WRITE(*,'(A)') "Invalid interpolation specification."
      STOP
    END SELECT
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
    !Set the number of Gauss points
    IF(INTERPOLATION_SPECIFICATIONS==OC_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(Basis,[3,3],Err)
    ENDIF
    !Finish the creation of the basis
    CALL OC_Basis_CreateFinish(Basis,Err)

    !Start the creation of a generated mesh in the region
    CALL OC_GeneratedMesh_Initialise(GENERATED_MESH,Err)
    CALL OC_GeneratedMesh_CreateStart(1,REGION,GENERATED_MESH,Err)
    !Set up a regular 100x100 mesh
    CALL OC_GeneratedMesh_TypeSet(GENERATED_MESH,1,Err)
    CALL OC_GeneratedMesh_BasisSet(GENERATED_MESH,Basis,Err)
    !Define the mesh on the region
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      CALL OC_GeneratedMesh_ExtentSet(GENERATED_MESH,[WIDTH,HEIGHT],Err)
      CALL OC_GeneratedMesh_NumberOfElementsSet(GENERATED_MESH,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS], &
        & Err)
      CALL OC_GeneratedMesh_OriginSet(GENERATED_MESH,ORIGIN,Err)
    ELSE
      CALL OC_GeneratedMesh_ExtentSet(GENERATED_MESH,[WIDTH,HEIGHT,LENGTH],Err)
      CALL OC_GeneratedMesh_NumberOfElementsSet(GENERATED_MESH,[NUMBER_GLOBAL_X_ELEMENTS, &
        & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS],Err)
    ENDIF
    !Finish the creation of a generated mesh in the region
    CALL OC_Mesh_Initialise(MESH,Err)
    CALL OC_GeneratedMesh_CreateFinish(GENERATED_MESH,1,MESH,Err)
    
    !Create a decomposition
    CALL OC_Decomposition_Initialise(DECOMPOSITION,Err)
    CALL OC_Decomposition_CreateStart(1,MESH,DECOMPOSITION,Err)
    CALL OC_Decomposition_CreateFinish(DECOMPOSITION,Err)

    !Decompose
    CALL OC_Decomposer_Initialise(decomposer,err)
    CALL OC_Decomposer_CreateStart(1,region,worldWorkGroup,decomposer,err)
    !Add in the decomposition
    CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
    !Finish the decomposer
    CALL OC_Decomposer_CreateFinish(decomposer,err)
  
    !Start to create a default (geometric) field on the region
    CALL OC_Field_Initialise(GEOMETRIC_FIELD,Err)
    CALL OC_Field_CreateStart(1,REGION,GEOMETRIC_FIELD,Err)
    !Set the decomposition to use
    CALL OC_Field_DecompositionSet(GEOMETRIC_FIELD,DECOMPOSITION,Err)
    !Set the domain to be used by the field components
    !NB these are needed now as the default mesh component number is 1
    CALL OC_Field_ComponentMeshComponentSet(GEOMETRIC_FIELD,OC_FIELD_U_VARIABLE_TYPE,1,1,Err)
    CALL OC_Field_ComponentMeshComponentSet(GEOMETRIC_FIELD,OC_FIELD_U_VARIABLE_TYPE,2,1,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      CALL OC_Field_ComponentMeshComponentSet(GEOMETRIC_FIELD,OC_FIELD_U_VARIABLE_TYPE,3,1,Err)
    ENDIF
    IF(INTERPOLATION_SPECIFICATIONS==OC_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
      CALL OC_Field_ScalingTypeSet(GEOMETRIC_FIELD,OC_FIELD_ARITHMETIC_MEAN_SCALING,Err)
    ENDIF
    !Finish creating the field
    CALL OC_Field_CreateFinish(GEOMETRIC_FIELD,Err)

    !Update the geometric field parameters
    CALL OC_GeneratedMesh_GeometricParametersCalculate(GENERATED_MESH,GEOMETRIC_FIELD,Err)

    !Create the equations_set
    CALL OC_EquationsSet_Initialise(EQUATIONS_SET,Err)
    CALL OC_Field_Initialise(EquationsSetField,Err)
    !Set the equations set to be a standard Laplace problem
    CALL OC_EquationsSet_CreateStart(1,REGION,GEOMETRIC_FIELD,[OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
      & OC_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,OC_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE],8, &
      & EquationsSetField,EQUATIONS_SET,Err)
    
    !Finish creating the equations set
    CALL OC_EquationsSet_CreateFinish(EQUATIONS_SET,Err)
  
    !Create the equations set dependent field variables
    CALL OC_Field_Initialise(DEPENDENT_FIELD,Err)
    CALL OC_EquationsSet_DependentCreateStart(EQUATIONS_SET,2,DEPENDENT_FIELD,Err)
    !Finish the equations set dependent field variables
    CALL OC_EquationsSet_DependentCreateFinish(EQUATIONS_SET,Err)

    !Create the equations set analytic field variables
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      ANALYTIC_FUNCTION=OC_EQUATIONS_SET_STANDARD_LAPLACE_EQUATION_THREE_DIM_2
    ELSE
      ANALYTIC_FUNCTION=OC_EQUATIONS_SET_STANDARD_LAPLACE_EQUATION_TWO_DIM_2
    ENDIF
    CALL OC_Field_Initialise(ANALYTIC_FIELD,Err)
    CALL OC_EquationsSet_AnalyticCreateStart(EQUATIONS_SET,ANALYTIC_FUNCTION,3,ANALYTIC_FIELD,Err)
    !Finish the equations set analtyic field variables
    CALL OC_EquationsSet_AnalyticCreateFinish(EQUATIONS_SET,Err)

    !Create the equations set equations
    CALL OC_Equations_Initialise(EQUATIONS,Err)
    CALL OC_EquationsSet_EquationsCreateStart(EQUATIONS_SET,EQUATIONS,Err)
    !Set the equations matrices sparsity type
    CALL OC_Equations_SparsityTypeSet(EQUATIONS,OC_EQUATIONS_SPARSE_MATRICES,Err)
    CALL OC_EquationsSet_EquationsCreateFinish(EQUATIONS_SET,Err)

    !Create the problem
    CALL OC_Problem_Initialise(PROBLEM,Err)
    CALL OC_Problem_CreateStart(1,context,[OC_PROBLEM_CLASSICAL_FIELD_CLASS,OC_PROBLEM_LAPLACE_EQUATION_TYPE, &
      & OC_PROBLEM_STANDARD_LAPLACE_SUBTYPE],PROBLEM,Err)
    !Finish creating the problem
    CALL OC_Problem_CreateFinish(PROBLEM,Err)

    !Create the problem control loop
    CALL OC_Problem_ControlLoopCreateStart(PROBLEM,Err)
    !Finish creating the problem control
    CALL OC_Problem_ControlLoopCreateFinish(PROBLEM,Err)

    !Start the creation of the problem solvers
    CALL OC_Solver_Initialise(Solver,Err)
    CALL OC_Problem_SolversCreateStart(Problem,Err)
    CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
    !Set solver to direct type
    CALL OC_Solver_LinearTypeSet(Solver,OC_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL OC_Solver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_OC_RP,Err)
    CALL OC_Solver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_OC_RP,Err)
    !CALL OC_Solver_LinearTypeSet(Solver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    !CALL OC_Solver_LibraryTypeSet(Solver,OC_SOLVER_MUMPS_LIBRARY,Err)
    !Finish the creation of the problem solver
    CALL OC_Problem_SolversCreateFinish(Problem,Err)

    !Start the creation of the problem solver equations
    CALL OC_Solver_Initialise(Solver,Err)
    CALL OC_SolverEquations_Initialise(Solver_Equations,Err)
    CALL OC_Problem_SolverEquationsCreateStart(Problem,Err)
    !Get the solve equations
    CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
    CALL OC_Solver_SolverEquationsGet(Solver,Solver_Equations,Err)
    !Set the solver equations sparsity
    CALL OC_SolverEquations_SparsityTypeSet(Solver_Equations,OC_SOLVER_SPARSE_MATRICES,Err)
    !CALL OC_SolverEquations_SparsityTypeSet(SolverEquations,OC_SOLVER_FULL_MATRICES,Err)
    !Add in the equations set
    CALL OC_SolverEquations_EquationsSetAdd(Solver_Equations,Equations_Set,EquationsSetIndex,Err)
    !Finish the creation of the problem solver equations
    CALL OC_Problem_SolverEquationsCreateFinish(Problem,Err)

    !Set up the boundary conditions as per the analytic solution
    CALL OC_BoundaryConditions_Initialise(BoundaryConditions,Err)
    CALL OC_SolverEquations_BoundaryConditionsCreateStart(Solver_Equations,BoundaryConditions,Err)
    CALL OC_SolverEquations_BoundaryConditionsAnalytic(Solver_Equations,Err)
    CALL OC_SolverEquations_BoundaryConditionsCreateFinish(Solver_Equations,Err)

    !Solve the problem
    CALL OC_Problem_Solve(PROBLEM,Err)

  END SUBROUTINE ANALYTICLAPLACE_GENERIC

  SUBROUTINE ANALYTICLAPLACE_GENERIC_CLEAN(contextUserNumber,CoordinateSystemUserNumber,RegionUserNumber, &
    & BasisUserNumber,GeneratedMeshUserNumber,ProblemUserNumber)

    !Argument variables
    INTEGER(OC_Intg), INTENT(IN) :: contextUserNumber
    INTEGER(OC_Intg), INTENT(IN) :: CoordinateSystemUserNumber
    INTEGER(OC_Intg), INTENT(IN) :: RegionUserNumber
    INTEGER(OC_Intg), INTENT(IN) :: BasisUserNumber
    INTEGER(OC_Intg), INTENT(IN) :: GeneratedMeshUserNumber
    INTEGER(OC_Intg), INTENT(IN) :: ProblemUserNumber

    CALL OC_Problem_Destroy(contextUserNumber,ProblemUserNumber,Err)
    CALL OC_GeneratedMesh_Destroy(contextUserNumber,RegionUserNumber,GeneratedMeshUserNumber,Err)
    CALL OC_Basis_Destroy(contextUserNumber,BasisUserNumber,Err)
    CALL OC_Region_Destroy(contextUserNumber,RegionUserNumber,Err)
    CALL OC_CoordinateSystem_Destroy(contextUserNumber,CoordinateSystemUserNumber,Err)

  END SUBROUTINE ANALYTICLAPLACE_GENERIC_CLEAN

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

END PROGRAM AnalyticLaplaceExample 
