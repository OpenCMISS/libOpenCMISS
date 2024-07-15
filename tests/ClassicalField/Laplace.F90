!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a Laplace equation using OpenCMISS calls.
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

!> \example ClassicalField/Laplace/Laplace/Fortran/src/LaplaceExample.F90
!! Example program to solve a Laplace equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/Laplace/Laplace/history.html
!!
!<

!> Main program
PROGRAM LaplaceExample

  USE OpenCMISS
  
#ifdef WITH_MPI
#ifdef WITH_F08_MPI
  USE MPI_F08
#elif WITH_F90_MPI 
  USE MPI
#endif  
#endif


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef WITH_MPI  
#ifdef WITH_F77_MPI
#include "mpif.h"
#endif
#endif  


  !Test program parameters

  REAL(OC_RP), PARAMETER :: HEIGHT=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: WIDTH=1.0_OC_RP
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
  INTEGER(OC_Intg), PARAMETER :: EquationsSetFieldUserNumber=10
  INTEGER(OC_Intg), PARAMETER :: DependentFieldUserNumber=11
  INTEGER(OC_Intg), PARAMETER :: EquationsSetUserNumber=12
  INTEGER(OC_Intg), PARAMETER :: ProblemUserNumber=13
 
  !Program types
  
  !Program variables

  INTEGER(OC_Intg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  INTEGER(OC_Intg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT,Filename

  !OpenCMISS variables

  TYPE(OC_BasisType) :: Basis
  TYPE(OC_BoundaryConditionsType) :: BoundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: ComputationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: CoordinateSystem
  TYPE(OC_DecompositionType) :: Decomposition
  TYPE(OC_DecomposerType) :: Decomposer
  TYPE(OC_EquationsType) :: Equations
  TYPE(OC_EquationsSetType) :: EquationsSet
  TYPE(OC_FieldType) :: GeometricField,EquationsSetField,DependentField
  TYPE(OC_FieldsType) :: Fields
  TYPE(OC_GeneratedMeshType) :: GeneratedMesh  
  TYPE(OC_MeshType) :: Mesh
  TYPE(OC_NodesType) :: Nodes
  TYPE(OC_ProblemType) :: Problem
  TYPE(OC_RegionType) :: Region,WorldRegion
  TYPE(OC_SolverType) :: Solver
  TYPE(OC_SolverEquationsType) :: SolverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic OpenOC variables
  
  INTEGER(OC_Intg) :: NumberOfComputationNodes,ComputationNodeNumber
  INTEGER(OC_Intg) :: decompositionIndex,EquationsSetIndex
  INTEGER(OC_Intg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(OC_Intg) :: FirstNodeDomain,LastNodeDomain
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
    IF(NUMBER_GLOBAL_Y_ELEMENTS<=0) CALL HandleError("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HandleError("Error for command argument 3.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Z_ELEMENTS
    IF(NUMBER_GLOBAL_Z_ELEMENTS<0) CALL HandleError("Invalid number of Z elements.")
    CALL GET_COMMAND_ARGUMENT(4,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HandleError("Error for command argument 4.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION_TYPE
    IF(INTERPOLATION_TYPE<=0) CALL HandleError("Invalid Interpolation specification.")
  ELSE
    !If there are not enough arguments default the problem specification 
    NUMBER_GLOBAL_X_ELEMENTS=1
    NUMBER_GLOBAL_Y_ELEMENTS=3
    NUMBER_GLOBAL_Z_ELEMENTS=1
!    INTERPOLATION_TYPE=1
    
    INTERPOLATION_TYPE=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
!    INTERPOLATION_TYPE=OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
!    INTERPOLATION_TYPE=OC_BASIS_CUBIC_LAGRANGE_INTERPOLATION    
    
  ENDIF
  
  !Intialise OpenCMISS
  CALL OC_Initialise(Err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,Err)
  CALL OC_DiagnosticsSetOn(OC_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"],Err)
  WRITE(Filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "Laplace",NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS,INTERPOLATION_TYPE  
  CALL OC_OutputSetOn(Filename,Err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(ContextUserNumber,context,Err)
  
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)
  
  CALL OC_Context_RandomSeedsSet(context,9999,Err)
  
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
  CALL OC_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL OC_Region_LabelSet(Region,"LaplaceRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL OC_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL OC_Region_CreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL OC_Basis_Initialise(Basis,Err)
  CALL OC_Basis_CreateStart(BasisUserNumber,context,Basis,Err)
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
    NUMBER_OF_GAUSS_XI=0 !Don't set number of Gauss points for tri/tet
  END SELECT
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bi-interpolation basis
    CALL OC_Basis_NumberOfXiSet(Basis,2,Err)
    CALL OC_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSE
    !Set the basis to be a tri-interpolation basis
    CALL OC_Basis_NumberOfXiSet(Basis,3,Err)
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
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
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
  
  !Destory the mesh now that we have decomposed it
  !CALL OC_Mesh_Destroy(Mesh,Err)
 
  !Start to create a default (geometric) field on the region
  CALL OC_Field_Initialise(GeometricField,Err)
  CALL OC_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL OC_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
  
  !Create the Standard Laplace Equations set
  CALL OC_EquationsSet_Initialise(EquationsSet,Err)
  CALL OC_Field_Initialise(EquationsSetField,Err)
  CALL OC_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & OC_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,OC_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL OC_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL OC_Field_Initialise(DependentField,Err)
  CALL OC_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Set the DOFs to be contiguous across components
  CALL OC_Field_DOFOrderTypeSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  CALL OC_Field_DOFOrderTypeSet(DependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,OC_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Initialise the field with an initial guess
  CALL OC_Field_ComponentValuesInitialise(DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,0.5_OC_RP, &
    & Err)

  !Create the equations set equations
  CALL OC_Equations_Initialise(Equations,Err)
  CALL OC_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL OC_Equations_SparsityTypeSet(Equations,OC_EQUATIONS_SPARSE_MATRICES,Err)
  !CALL OC_Equations_SparsityTypeSet(Equations,OC_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_NO_OUTPUT,Err)
  !CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL OC_EquationsSet_EquationsCreateFinish(EquationsSet,Err)
  
  !Start the creation of a problem.
  CALL OC_Problem_Initialise(Problem,Err)
  CALL OC_Problem_CreateStart(ProblemUserNumber,context,[OC_PROBLEM_CLASSICAL_FIELD_CLASS,OC_PROBLEM_LAPLACE_EQUATION_TYPE, &
    & OC_PROBLEM_STANDARD_LAPLACE_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL OC_Problem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL OC_Problem_ControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL OC_Solver_Initialise(Solver,Err)
  CALL OC_Problem_SolversCreateStart(Problem,Err)
  CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_NO_OUTPUT,Err)
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_TIMING_OUTPUT,Err)
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_SOLVER_OUTPUT,Err)
  !CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_MATRIX_OUTPUT,Err)
  
  CALL OC_Solver_LinearTypeSet(Solver,OC_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
  CALL OC_Solver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_OC_RP,Err)
  CALL OC_Solver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_OC_RP,Err)

  !CALL OC_Solver_LinearTypeSet(Solver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)  
  !CALL OC_Solver_LinearTypeSet(Solver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  !CALL OC_Solver_LibraryTypeSet(Solver,OC_SOLVER_MUMPS_LIBRARY,Err)
  !CALL OC_Solver_LibraryTypeSet(Solver,OC_SOLVER_LAPACK_LIBRARY,Err)
  !CALL OC_Solver_LibraryTypeSet(Solver,OC_SOLVER_SUPERLU_LIBRARY,Err)
  !CALL OC_Solver_LibraryTypeSet(Solver,OC_SOLVER_PASTIX_LIBRARY,Err)
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

  !Start the creation of the equations set boundary conditions
  CALL OC_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  CALL OC_Nodes_Initialise(Nodes,Err)
  CALL OC_Region_NodesGet(Region,Nodes,Err)
  CALL OC_Nodes_NumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL OC_Decomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL OC_Decomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationNodeNumber) THEN
    CALL OC_BoundaryConditions_SetNode(BoundaryConditions,DependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
      & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationNodeNumber) THEN
    CALL OC_BoundaryConditions_SetNode(BoundaryConditions,DependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
      & OC_BOUNDARY_CONDITION_FIXED,1.0_OC_RP,Err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL OC_Problem_Solve(Problem,Err)

  !Export results
  CALL OC_Fields_Initialise(Fields,Err)
  CALL OC_Fields_Create(Region,Fields,Err)
  CALL OC_Fields_NodesExport(Fields,"Laplace","FORTRAN",Err)
  CALL OC_Fields_ElementsExport(Fields,"Laplace","FORTRAN",Err)
  CALL OC_Fields_Finalise(Fields,Err)
  
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
    
END PROGRAM LaplaceExample
