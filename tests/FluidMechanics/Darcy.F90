!> \file
!> \author Christian Michler
!> \brief This is an example program to solve an analytic Darcy equation using OpenCMISS calls.
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

!> \example FluidMechanics/Darcy/Analytic/src/AnalyticExample.F90
!! Example program to solve an analytic Darcy equation using OpenCMISS calls.
!!
!! \htmlinclude FluidMechanics/Darcy/Analytic/history.html
!<

! ! 
! !  This example considers an analytic Darcy problem.
! ! 

!> Main program

PROGRAM DarcyAnalyticExample

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OpenCMISS
  USE FLUID_MECHANICS_IO_ROUTINES
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWINCMISS
#endif

  !
  !================================================================================================================================
  !

  !PROGRAM VARIABLES AND TYPES

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  INTEGER(OC_Intg), PARAMETER :: ContextUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: CoordinateSystemUserNumber=2
  INTEGER(OC_Intg), PARAMETER :: RegionUserNumber=3
  INTEGER(OC_Intg), PARAMETER :: MeshUserNumber=4
  INTEGER(OC_Intg), PARAMETER :: DecompositionUserNumber=5
  INTEGER(OC_Intg), PARAMETER :: DecomposerUserNumber=6
  INTEGER(OC_Intg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(OC_Intg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(OC_Intg), PARAMETER :: DependentFieldUserNumberDarcy=9
  INTEGER(OC_Intg), PARAMETER :: MaterialsFieldUserNumberDarcy=10
  INTEGER(OC_Intg), PARAMETER :: AnalyticFieldUserNumberDarcy=11
  INTEGER(OC_Intg), PARAMETER :: EquationsSetUserNumberDarcy=12
  INTEGER(OC_Intg), PARAMETER :: ProblemUserNumber=13

  INTEGER(OC_Intg), PARAMETER :: DomainUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: SolverDarcyUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: MaterialsFieldUserNumberDarcyPorosity=1
  INTEGER(OC_Intg), PARAMETER :: MaterialsFieldUserNumberDarcyPermOverVis=2

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(OC_Intg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(OC_Intg) :: BASIS_TYPE
  INTEGER(OC_Intg) :: BASIS_NUMBER_GEOMETRY
  INTEGER(OC_Intg) :: BASIS_NUMBER_VELOCITY
  INTEGER(OC_Intg) :: BASIS_NUMBER_PRESSURE
  INTEGER(OC_Intg) :: BASIS_XI_GAUSS_GEOMETRY
  INTEGER(OC_Intg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(OC_Intg) :: BASIS_XI_GAUSS_PRESSURE
  INTEGER(OC_Intg) :: BASIS_XI_INTERPOLATION_GEOMETRY
  INTEGER(OC_Intg) :: BASIS_XI_INTERPOLATION_VELOCITY
  INTEGER(OC_Intg) :: BASIS_XI_INTERPOLATION_PRESSURE
  INTEGER(OC_Intg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(OC_Intg) :: MESH_COMPONENT_NUMBER_GEOMETRY
  INTEGER(OC_Intg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(OC_Intg) :: MESH_COMPONENT_NUMBER_PRESSURE
  INTEGER(OC_Intg) :: NUMBER_OF_NODES_GEOMETRY
  INTEGER(OC_Intg) :: NUMBER_OF_NODES_VELOCITY
  INTEGER(OC_Intg) :: NUMBER_OF_NODES_PRESSURE
  INTEGER(OC_Intg) :: NUMBER_OF_ELEMENT_NODES_GEOMETRY
  INTEGER(OC_Intg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY
  INTEGER(OC_Intg) :: NUMBER_OF_ELEMENT_NODES_PRESSURE
  INTEGER(OC_Intg) :: TOTAL_NUMBER_OF_NODES
  INTEGER(OC_Intg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(OC_Intg) :: MAXIMUM_ITERATIONS
  INTEGER(OC_Intg) :: RESTART_VALUE
!   INTEGER(OC_Intg) :: MPI_IERROR

  INTEGER(OC_Intg) :: EQUATIONS_DARCY_OUTPUT
  INTEGER(OC_Intg) :: COMPONENT_NUMBER
  INTEGER(OC_Intg) :: NODE_NUMBER
  INTEGER(OC_Intg) :: ELEMENT_NUMBER

  INTEGER(OC_Intg) :: LINEAR_SOLVER_DARCY_OUTPUT_TYPE

  INTEGER(OC_Intg) :: ANALYTICAL_TYPE
  INTEGER(OC_Intg) :: INPUT_TYPE

  REAL(OC_RP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2, DOMAIN_LENGTH
  REAL(OC_RP) :: GEOMETRY_TOLERANCE

  REAL(OC_RP) :: INITIAL_FIELD_DARCY(3)
  REAL(OC_RP) :: DIVERGENCE_TOLERANCE
  REAL(OC_RP) :: RELATIVE_TOLERANCE
  REAL(OC_RP) :: ABSOLUTE_TOLERANCE
  REAL(OC_RP) :: LINESEARCH_ALPHA
  REAL(OC_RP) :: VALUE
  REAL(OC_RP) :: POROSITY_PARAM_DARCY, PERM_OVER_VIS_PARAM_DARCY

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_DARCY_DIRECT_FLAG

  CHARACTER *15 BUFFER
  CHARACTER *15 OUTPUT_STRING

  !CMISS variables

  !Context
  TYPE(OC_ContextType) :: context
  !Computation environment
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  !Regions  
  TYPE(OC_RegionType) :: Region
  TYPE(OC_RegionType) :: WorldRegion
  !Coordinate systems
  TYPE(OC_CoordinateSystemType) :: CoordinateSystem
  !Basis
  TYPE(OC_BasisType) :: BasisGeometry
  TYPE(OC_BasisType) :: BasisVelocity
  TYPE(OC_BasisType) :: BasisPressure
  !Nodes
  TYPE(OC_NodesType) :: Nodes
  !Elements
  TYPE(OC_MeshElementsType) :: MeshElementsGeometry
  TYPE(OC_MeshElementsType) :: MeshElementsVelocity
  TYPE(OC_MeshElementsType) :: MeshElementsPressure
  !Meshes
  TYPE(OC_MeshType) :: Mesh
  !Decompositions
  TYPE(OC_DecompositionType) :: Decomposition
  TYPE(OC_DecomposerType) :: Decomposer
  !Fields
  TYPE(OC_FieldsType) :: Fields
  !Field types
  TYPE(OC_FieldType) :: GeometricField
  TYPE(OC_FieldType) :: EquationsSetField
  TYPE(OC_FieldType) :: DependentFieldDarcy
  TYPE(OC_FieldType) :: MaterialsFieldDarcy
  TYPE(OC_FieldType) :: AnalyticFieldDarcy
  !Equations sets
  TYPE(OC_EquationsSetType) :: EquationsSetDarcy
  !Equations
  TYPE(OC_EquationsType) :: EquationsDarcy
  !Problems
  TYPE(OC_ProblemType) :: Problem
  !Control loops
  TYPE(OC_ControlLoopType) :: ControlLoop
  !Solvers
  TYPE(OC_SolverType) :: LinearSolverDarcy
  !Solver equations
  TYPE(OC_SolverEquationsType) :: SolverEquationsDarcy
  TYPE(OC_BoundaryConditionsType) :: BoundaryConditions
  !Work group
  TYPE(OC_WorkGroupType) :: worldWorkGroup

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(OC_Intg) :: decompositionIndex,EquationsSetIndex
  INTEGER(OC_Intg) :: computationNodeNumber,numberOfComputationNodes
  INTEGER(OC_Intg) :: Err


  INTEGER(OC_Intg) :: DIAG_LEVEL_LIST(5)
  CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(1) !,TIMING_ROUTINE_LIST(1)

  
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

  !
  !================================================================================================================================
  !

  !PROBLEM CONTROL PANEL

  !INITIALISE OPENCMISS
  CALL OC_Initialise(Err)  
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,Err)
 
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
  
  !  
  !================================================================================================================================
  !

  !Import cmHeart mesh information
  CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)  
  BASIS_NUMBER_GEOMETRY=CM%ID_M
  BASIS_NUMBER_VELOCITY=CM%ID_V
  BASIS_NUMBER_PRESSURE=CM%ID_P
  NUMBER_OF_DIMENSIONS=CM%D
  BASIS_TYPE=CM%IT_T
  BASIS_XI_INTERPOLATION_GEOMETRY=CM%IT_M
  BASIS_XI_INTERPOLATION_VELOCITY=CM%IT_V
  BASIS_XI_INTERPOLATION_PRESSURE=CM%IT_P
  NUMBER_OF_NODES_GEOMETRY=CM%N_M
  NUMBER_OF_NODES_VELOCITY=CM%N_V
  NUMBER_OF_NODES_PRESSURE=CM%N_P
  TOTAL_NUMBER_OF_NODES=CM%N_T
  TOTAL_NUMBER_OF_ELEMENTS=CM%E_T
  NUMBER_OF_ELEMENT_NODES_GEOMETRY=CM%EN_M
  NUMBER_OF_ELEMENT_NODES_VELOCITY=CM%EN_V
  NUMBER_OF_ELEMENT_NODES_PRESSURE=CM%EN_P
  !Set domain dimensions
  DOMAIN_X1 = -5.0_OC_RP
  DOMAIN_X2 =  5.0_OC_RP
  DOMAIN_Y1 = -5.0_OC_RP
  DOMAIN_Y2 =  5.0_OC_RP
  DOMAIN_Z1 = -5.0_OC_RP
  DOMAIN_Z2 =  5.0_OC_RP
  !Set domain length
  DOMAIN_LENGTH = 10.0_OC_RP
  !Set geometric tolerance
  GEOMETRY_TOLERANCE = 1.0E-12_OC_RP
  !Set initial values
  INITIAL_FIELD_DARCY(1)=0.0_OC_RP
  INITIAL_FIELD_DARCY(2)=0.0_OC_RP
  INITIAL_FIELD_DARCY(3)=0.0_OC_RP
  !Set material parameters
!   POROSITY_PARAM_DARCY=0.3_OC_RP  !??? Also try with porosity unequal 1.0
  POROSITY_PARAM_DARCY=1.0_OC_RP
  PERM_OVER_VIS_PARAM_DARCY=1.0_OC_RP  !The value of 1.0 is also hard-coded as PERM_OVER_VIS_PARAM in DARCY_EQUATION_ANALYTIC_CALCULATE
  !Set number of Gauss points (Mind that also material field may be interpolated)
  BASIS_XI_GAUSS_GEOMETRY=3 !4
  BASIS_XI_GAUSS_VELOCITY=3 !4
  BASIS_XI_GAUSS_PRESSURE=3 !4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=OC_SOLVER_SOLVER_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DARCY_OUTPUT=OC_EQUATIONS_NO_OUTPUT
  !Set solver parameters
  LINEAR_SOLVER_DARCY_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-14_OC_RP !default: 1.0E-05_OC_RP
  ABSOLUTE_TOLERANCE=1.0E-14_OC_RP !default: 1.0E-10_OC_RP
  DIVERGENCE_TOLERANCE=1.0E5_OC_RP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000_OC_Intg !default: 100000
  RESTART_VALUE=30_OC_Intg !default: 30
  LINESEARCH_ALPHA=1.0_OC_RP

  !
  !================================================================================================================================
  !

  ! Available analytical cases:
  ! 1=OC_EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1
  ! 2=OC_EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2
  ! 3=OC_EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3
  ! 4=OC_EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1
  ! 5=OC_EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2
  ! 6=OC_EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3


  WRITE(*,*)'1=POLYNOM, 2=EXP, 3=COS/SIN:'
!   READ(*,*) 

  IF(COMMAND_ARGUMENT_COUNT()==2) THEN
    CALL GET_COMMAND_ARGUMENT(1,BUFFER)
    READ(BUFFER,*) INPUT_TYPE
    CALL GET_COMMAND_ARGUMENT(2,BUFFER)
    READ(BUFFER,*) OUTPUT_STRING
  ELSE
    !TODO more detailed error message
    WRITE(*,*)'INPUT ERROR!!!'
  ENDIF


  IF(INPUT_TYPE==1.AND.NUMBER_OF_DIMENSIONS==2) ANALYTICAL_TYPE=OC_EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_1
  IF(INPUT_TYPE==2.AND.NUMBER_OF_DIMENSIONS==2) ANALYTICAL_TYPE=OC_EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_2
  IF(INPUT_TYPE==3.AND.NUMBER_OF_DIMENSIONS==2) ANALYTICAL_TYPE=OC_EQUATIONS_SET_DARCY_EQUATION_TWO_DIM_3
  IF(INPUT_TYPE==1.AND.NUMBER_OF_DIMENSIONS==3) ANALYTICAL_TYPE=OC_EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_1
  IF(INPUT_TYPE==2.AND.NUMBER_OF_DIMENSIONS==3) ANALYTICAL_TYPE=OC_EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_2
  IF(INPUT_TYPE==3.AND.NUMBER_OF_DIMENSIONS==3) ANALYTICAL_TYPE=OC_EQUATIONS_SET_DARCY_EQUATION_THREE_DIM_3

  !
  !================================================================================================================================
  !

  !Set diagnostics

  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5

  DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_FINITE_ELEMENT_CALCULATE"

  !OC_ALL_DIAG_TYPE/OC_IN_DIAG_TYPE/OC_FROM_DIAG_TYPE
!   CALL OC_DiagnosticsSetOn(OC_IN_DIAG_TYPE,DIAG_LEVEL_LIST,"DarcyDiagnostics",DIAG_ROUTINE_LIST,Err)

  !OC_ALL_TIMING_TYPE/OC_IN_TIMING_TYPE/OC_FROM_TIMING_TYPE
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL OC_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL OC_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,context,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL OC_CoordinateSystem_DimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL OC_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION

  !Start the creation of a new region
  CALL OC_Region_Initialise(Region,Err)
  CALL OC_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL OC_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL OC_Region_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !BASES

  !Start the creation of new bases
  MESH_NUMBER_OF_COMPONENTS=1
  CALL OC_Basis_Initialise(BasisGeometry,Err)
  CALL OC_Basis_CreateStart(BASIS_NUMBER_GEOMETRY,context,BasisGeometry,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL OC_Basis_TypeSet(BasisGeometry,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL OC_Basis_NumberOfXiSet(BasisGeometry,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS==2) THEN
    CALL OC_Basis_InterpolationXiSet(BasisGeometry,[BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY],Err)
    CALL OC_Basis_QuadratureNumberOfGaussXiSet(BasisGeometry,[BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY],Err)
  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
    CALL OC_Basis_InterpolationXiSet(BasisGeometry,[BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, & 
      & BASIS_XI_INTERPOLATION_GEOMETRY],Err)                         
    CALL OC_Basis_QuadratureNumberOfGaussXiSet(BasisGeometry,[BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
      & BASIS_XI_GAUSS_GEOMETRY],Err)
  ENDIF
  !Finish the creation of the basis
  CALL OC_Basis_CreateFinish(BasisGeometry,Err)
  !
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisVelocity=BasisGeometry
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new velocity basis
    CALL OC_Basis_Initialise(BasisVelocity,Err)
    !Start the creation of a basis
    CALL OC_Basis_CreateStart(BASIS_NUMBER_VELOCITY,context,BasisVelocity,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL OC_Basis_TypeSet(BasisVelocity,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL OC_Basis_NumberOfXiSet(BasisVelocity,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL OC_Basis_InterpolationXiSet(BasisVelocity,[BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY],Err)
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(BasisVelocity,[BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY],Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL OC_Basis_InterpolationXiSet(BasisVelocity,[BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY, & 
        & BASIS_XI_INTERPOLATION_VELOCITY],Err)                         
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(BasisVelocity,[BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY, & 
        & BASIS_XI_GAUSS_VELOCITY],Err)
    ENDIF
    !Finish the creation of the basis
    CALL OC_Basis_CreateFinish(BasisVelocity,Err)
  ENDIF
  !
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisPressure=BasisGeometry
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_VELOCITY) THEN
    BasisPressure=BasisVelocity
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new pressure basis
    CALL OC_Basis_Initialise(BasisPressure,Err)
    !Start the creation of a basis
    CALL OC_Basis_CreateStart(BASIS_NUMBER_PRESSURE,context,BasisPressure,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL OC_Basis_TypeSet(BasisPressure,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL OC_Basis_NumberOfXiSet(BasisPressure,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL OC_Basis_InterpolationXiSet(BasisPressure,[BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE],Err)
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(BasisPressure,[BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE],Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL OC_Basis_InterpolationXiSet(BasisPressure,[BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE, & 
        & BASIS_XI_INTERPOLATION_PRESSURE],Err)                         
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(BasisPressure,[BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE, & 
        & BASIS_XI_GAUSS_PRESSURE],Err)
    ENDIF
    !Finish the creation of the basis
    CALL OC_Basis_CreateFinish(BasisPressure,Err)
  ENDIF

  !
  !================================================================================================================================
  !

  !MESH

  !Start the creation of mesh nodes
  CALL OC_Nodes_Initialise(Nodes,Err)
  CALL OC_Nodes_CreateStart(Region,TOTAL_NUMBER_OF_NODES,Nodes,Err)
  CALL OC_Nodes_CreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL OC_Mesh_Initialise(Mesh,Err)
  CALL OC_Mesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL OC_Mesh_NumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL OC_Mesh_NumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)
  !Specify spatial mesh component
  CALL OC_MeshElements_Initialise(MeshElementsGeometry,Err)
  CALL OC_MeshElements_Initialise(MeshElementsVelocity,Err)
  CALL OC_MeshElements_Initialise(MeshElementsPressure,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_PRESSURE=1
  CALL OC_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL OC_MeshElements_NodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
  ENDDO
  CALL OC_MeshElements_CreateFinish(MeshElementsGeometry,Err)
  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsVelocity=MeshElementsGeometry
  ELSE
    MESH_COMPONENT_NUMBER_VELOCITY=MESH_COMPONENT_NUMBER_GEOMETRY+1
    CALL OC_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_VELOCITY,BasisVelocity,MeshElementsVelocity,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL OC_MeshElements_NodesSet(MeshElementsVelocity,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_VELOCITY),Err)
    ENDDO
    CALL OC_MeshElements_CreateFinish(MeshElementsVelocity,Err)
  ENDIF
  !Specify pressure mesh component
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsPressure=MeshElementsGeometry
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_GEOMETRY
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_VELOCITY) THEN
    MeshElementsPressure=MeshElementsVelocity
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_VELOCITY
  ELSE
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_VELOCITY+1
    CALL OC_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_PRESSURE,BasisPressure,MeshElementsPressure,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL OC_MeshElements_NodesSet(MeshElementsPressure,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_PRESSURE),Err)
    ENDDO
    CALL OC_MeshElements_CreateFinish(MeshElementsPressure,Err)
  ENDIF
  !Finish the creation of the mesh
  CALL OC_Mesh_CreateFinish(Mesh,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

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
  !Set the field type
  CALL OC_Field_TypeSet(GeometricField,OC_FIELD_GEOMETRIC_TYPE,Err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL OC_Field_ScalingTypeSet(GeometricField,OC_FIELD_NO_SCALING,Err)
  !Set the mesh component to be used by the field components.

  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
  ENDDO

  !Finish creating the field
  CALL OC_Field_CreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL OC_Field_ParameterSetUpdateNode(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1, & 
        & OC_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL OC_Field_ParameterSetUpdateStart(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,Err)
  CALL OC_Field_ParameterSetUpdateFinish(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for ALE Darcy
  CALL OC_EquationsSet_Initialise(EquationsSetDarcy,Err)
  CALL OC_Field_Initialise(EquationsSetField,Err)
  CALL OC_EquationsSet_CreateStart(EquationsSetUserNumberDarcy,Region,GeometricField,[OC_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & OC_EQUATIONS_SET_DARCY_EQUATION_TYPE,OC_EQUATIONS_SET_STANDARD_DARCY_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSetDarcy,Err)
  !Finish creating the equations set
  CALL OC_EquationsSet_CreateFinish(EquationsSetDarcy,Err)

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for ALE Darcy
  CALL OC_Field_Initialise(DependentFieldDarcy,Err)
  CALL OC_EquationsSet_DependentCreateStart(EquationsSetDarcy,DependentFieldUserNumberDarcy, & 
    & DependentFieldDarcy,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL OC_Field_ComponentMeshComponentSet(DependentFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
    CALL OC_Field_ComponentMeshComponentSet(DependentFieldDarcy,OC_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
  ENDDO
  COMPONENT_NUMBER=NUMBER_OF_DIMENSIONS+1
    CALL OC_Field_ComponentMeshComponentSet(DependentFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
    CALL OC_Field_ComponentMeshComponentSet(DependentFieldDarcy,OC_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_DependentCreateFinish(EquationsSetDarcy,Err)

  !Initialise dependent field (velocity components)
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL OC_Field_ComponentValuesInitialise(DependentFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, & 
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for ALE Darcy
  CALL OC_Field_Initialise(MaterialsFieldDarcy,Err)
  CALL OC_EquationsSet_MaterialsCreateStart(EquationsSetDarcy,MaterialsFieldUserNumberDarcy, & 
    & MaterialsFieldDarcy,Err)
  !Finish the equations set materials field variables
  CALL OC_EquationsSet_MaterialsCreateFinish(EquationsSetDarcy,Err)
  CALL OC_Field_ComponentValuesInitialise(MaterialsFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
  CALL OC_Field_ComponentValuesInitialise(MaterialsFieldDarcy,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)

  !
  !================================================================================================================================
  !

  !ANALYTIC FIELDS

  !Create the equations set analytic field variables for static Darcy
  CALL OC_Field_Initialise(AnalyticFieldDarcy,Err)
  !--- No distinction wrt. number dimension necessary ???
  IF(NUMBER_OF_DIMENSIONS==2) THEN  
    CALL OC_EquationsSet_AnalyticCreateStart(EquationsSetDarcy,ANALYTICAL_TYPE,AnalyticFieldUserNumberDarcy, &
      & AnalyticFieldDarcy,Err)
  ELSE
    CALL OC_EquationsSet_AnalyticCreateStart(EquationsSetDarcy,ANALYTICAL_TYPE,AnalyticFieldUserNumberDarcy, &
      & AnalyticFieldDarcy,Err)
  ENDIF
  !Finish the equations set analytic field variables
  CALL OC_EquationsSet_AnalyticCreateFinish(EquationsSetDarcy,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS

  !Create the equations set equations
  CALL OC_Equations_Initialise(EquationsDarcy,Err)
  CALL OC_EquationsSet_EquationsCreateStart(EquationsSetDarcy,EquationsDarcy,Err)
  !Set the equations matrices sparsity type
  CALL OC_Equations_SparsityTypeSet(EquationsDarcy,OC_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL OC_Equations_OutputTypeSet(EquationsDarcy,EQUATIONS_DARCY_OUTPUT,Err)
  !Finish the equations set equations
  CALL OC_EquationsSet_EquationsCreateFinish(EquationsSetDarcy,Err)

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL OC_Problem_Initialise(Problem,Err)
  CALL OC_ControlLoop_Initialise(ControlLoop,Err)
  CALL OC_Problem_CreateStart(ProblemUserNumber,context,[OC_PROBLEM_FLUID_MECHANICS_CLASS,OC_PROBLEM_DARCY_EQUATION_TYPE, &
    & OC_PROBLEM_STANDARD_DARCY_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL OC_Problem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL OC_Problem_ControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL OC_Solver_Initialise(LinearSolverDarcy,Err)
  CALL OC_Problem_SolversCreateStart(Problem,Err)
  !Get the Darcy solver
  CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,SolverDarcyUserNumber,LinearSolverDarcy,Err)
  !Set the output type
  CALL OC_Solver_OutputTypeSet(LinearSolverDarcy,LINEAR_SOLVER_DARCY_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_DARCY_DIRECT_FLAG) THEN
    CALL OC_Solver_LinearTypeSet(LinearSolverDarcy,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL OC_Solver_LibraryTypeSet(LinearSolverDarcy,OC_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL OC_Solver_LinearTypeSet(LinearSolverDarcy,OC_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL OC_Solver_LinearIterativeMaximumIterationsSet(LinearSolverDarcy,MAXIMUM_ITERATIONS,Err)
    CALL OC_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverDarcy,DIVERGENCE_TOLERANCE,Err)
    CALL OC_Solver_LinearIterativeRelativeToleranceSet(LinearSolverDarcy,RELATIVE_TOLERANCE,Err)
    CALL OC_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverDarcy,ABSOLUTE_TOLERANCE,Err)
    CALL OC_Solver_LinearIterativeGMRESRestartSet(LinearSolverDarcy,RESTART_VALUE,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL OC_Problem_SolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL OC_Solver_Initialise(LinearSolverDarcy,Err)
  CALL OC_SolverEquations_Initialise(SolverEquationsDarcy,Err)

  CALL OC_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the Darcy solver equations
  CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,SolverDarcyUserNumber,LinearSolverDarcy,Err)
  CALL OC_Solver_SolverEquationsGet(LinearSolverDarcy,SolverEquationsDarcy,Err)
  !Set the solver equations sparsity
  CALL OC_SolverEquations_SparsityTypeSet(SolverEquationsDarcy,OC_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL OC_SolverEquations_EquationsSetAdd(SolverEquationsDarcy,EquationsSetDarcy,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Set up the boundary conditions as per the analytic solution
  CALL OC_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditions,Err)
  CALL OC_SolverEquations_BoundaryConditionsAnalytic(SolverEquationsDarcy,Err)
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsDarcy,Err)

  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL OC_Problem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  !OUTPUT

  !Output Analytic Analysis
  CALL OC_AnalyticAnalysis_Output(DependentFieldDarcy,OUTPUT_STRING,Err)


  EXPORT_FIELD_IO=.FALSE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL OC_Fields_Initialise(Fields,Err)
    CALL OC_Fields_Create(Region,Fields,Err)
    CALL OC_Fields_NodesExport(Fields,"Darcy","FORTRAN",Err)
    CALL OC_Fields_ElementsExport(Fields,"Darcy","FORTRAN",Err)
    CALL OC_Fields_Finalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF

  !Destroy the context
  CALL OC_Context_Destroy(context,Err)
  !Finialise OpenCMISS
  CALL OC_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM DarcyAnalyticExample
