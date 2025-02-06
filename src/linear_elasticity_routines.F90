!> \file
!> \author Chris Bradley
!> \brief This module handles all linear elasticity routines.
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
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s): Chris Bradley
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

!>This module handles all linear elasticity routines.
MODULE LinearElasticityRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE CoordinateSystemAccessRoutines
  USE DataProjectionAccessRoutines
  USE DecompositionAccessRoutines
  USE DistributedMatrixVector
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE MatrixVector
  USE ProblemAccessRoutines
  USE ProfilingRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE Timer
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC LinearElasticity_BoundaryConditionsAnalyticCalculate

  PUBLIC LinearElasticity_EquationsSetDerivedVariableCalculate
  
  PUBLIC LinearElasticity_EquationsSetSetup
  
  PUBLIC LinearElasticity_EquationsSetSolutionMethodSet
  
  PUBLIC LinearElasticity_EquationsSetSpecificationSet
  
  PUBLIC LinearElasticity_FiniteElementCalculate

  PUBLIC LinearElasticity_ProblemSpecificationSet
  
  PUBLIC LinearElasticity_ProblemSetup

CONTAINS

  !
  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE LinearElasticity_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the analytic boundary conditions for
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<The boundary conditions to set the analytic boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_VARIABLES=99
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_COMPONENTS=99
    INTEGER(INTG) :: analyticFunctionType,componentIdx,derivativeIdx,dimensionIdx,equationsSetSubtype,esSpecification(3), &
      & globalDerivativeIndex,localDOFIdx,nodeIdx,numberOfAxisNodes(MAX_NUMBER_OF_COMPONENTS,MAX_NUMBER_OF_VARIABLES), &
      & numberOfComponents,numberOfDimensions,numberOfNodes,numberOfNodeDerivatives,numberOfVariables,numberofVersions, &
      & variableIdx,variableType,versionIdx
    REAL(DP) :: analyticDisplacement,analyticValue,E,F,H,Iyy,L,maxExtent(MAX_NUMBER_OF_COMPONENTS),normal(3), &
      & position(3,MAXIMUM_GLOBAL_DERIV_NUMBER),tangents(3,2),W,x(3)
    REAL(DP), PARAMETER :: GEOMETRIC_TOLERANCE=1.0E-6_DP
    REAL(DP), POINTER :: analyticParameters(:)
    LOGICAL :: onAxis,setBC
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField
    TYPE(FieldVariableType), POINTER :: analyticVariable,dependentVariable,geometricVariable
    TYPE(VARYING_STRING) :: localError,localWarning
    
    ENTERS("LinearElasticity_BoundaryConditionsAnalyticCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    equationsSetSubType=esSpecification(3)
    CALL EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*999)
    CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)

    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(geometricVariable)
    CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
    NULLIFY(analyticField)
    CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
    NULLIFY(analyticVariable)
    NULLIFY(analyticParameters)
    IF(ASSOCIATED(analyticField)) THEN
      CALL Field_VariableGet(analyticField,FIELD_U_VARIABLE_TYPE,analyticVariable,err,error,*999)
      CALL FieldVariable_ParameterSetDataGet(analyticVariable,FIELD_VALUES_SET_TYPE,analyticParameters,err,error,*999)
    ENDIF
    !
    ! Identify the dimensions of the domain in order to work out how many nodes, elements etc. 
    !
    IF(numberOfVariables>MAX_NUMBER_OF_VARIABLES) THEN
      localError="The number of variables of "//TRIM(NumberToVString(numberOfVariables,"*",err,error))// &
        & " is greater than the maximum number of variables. Increase MAX_NUMBER_OF_VARIABLES."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    numberOfAxisNodes = 0
    maxExtent = 0.0_DP
    DO variableIdx=1,numberOfVariables      
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
      IF(numberOfComponents>MAX_NUMBER_OF_COMPONENTS) THEN
        localError="The number of components of "//TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
          & " is greater than the maximum number of components. Increase MAX_NUMBER_OF_COMPONENTS."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DO componentIdx=1,numberOfComponents
        CALL FieldVariable_ComponentInterpolationCheck(dependentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
          & err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        !Loop over the local nodes excluding the ghosts.
        CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
        DO nodeIdx=1,numberOfNodes
          !Find the position of the node (interpolated since the geometric field might not have these nodes)
          CALL FieldVariable_PositionNormalTangentsCalculateNode(dependentVariable,componentIdx,nodeIdx, &
            & position,normal,tangents,err,error,*999)
          x(1:numberOfDimensions)=position(1:numberOfDimensions,NO_PART_DERIV)
          onAxis=.TRUE.
          DO dimensionIdx=1,numberOfDimensions
            IF(dimensionIdx/=componentIdx) THEN
              onAxis=onAxis.AND.(ABS(x(dimensionIdx))<=GEOMETRIC_TOLERANCE)
            ENDIF
          ENDDO !dimensionIdx
          IF(onAxis) THEN
            numberOfAxisNodes(componentIdx,variableIdx)=numberOfAxisNodes(componentIdx,variableIdx)+1
            IF(x(componentIdx)>maxExtent(componentIdx)) maxExtent(componentIdx)=x(componentIdx)
          ENDIF
        ENDDO !nodeIdx
      ENDDO !componentIdx
    ENDDO !variableIdx
    !Santity check
    SELECT CASE(equationsSetSubType)
    CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
      SELECT CASE(analyticFunctionType)
      CASE DEFAULT
        !Do nothing
      END SELECT
    CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE)     
      SELECT CASE(analyticFunctionType)
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_CANTILEVER_END_LOAD)
        !Check we have nodes
        DO variableIdx=1,numberOfVariables
          DO dimensionIdx=1,numberOfDimensions
            IF(numberOfAxisNodes(dimensionIdx,variableIdx)==0) THEN
              localError="The number of nodes along the axis of direction "//TRIM(NumberToVString(dimensionIdx,"*",err,error))// &
                & " of variable "//TRIM(NumberToVString(variableIdx,"*",err,error))//" of "// &
                & TRIM(NumberToVString(numberOfAxisNodes(dimensionIdx,variableIdx),"*",err,error))// &
                & " is invalid. There must be > 0 axis nodes."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !dimnsionIdx
        ENDDO !variableIdx
       !Check dimensions
        L=analyticParameters(1)
        W=analyticParameters(2)
        H=analyticParameters(3)        
        IF(ABS(maxExtent(1)-L)>GEOMETRIC_TOLERANCE) THEN
          localWarning="The maximum extent in the x direction of "//TRIM(NumberToVString(maxExtent(1),"*",err,error))// &
            & " is different from the analytic field L parameter of "//TRIM(NumberToVString(L,"*",err,error))//"."
          CALL FlagWarning(localWarning,err,error,*999)
        ENDIF
        IF(ABS(maxExtent(2)-W)>GEOMETRIC_TOLERANCE) THEN
          localWarning="The maximum extent in the y direction of "//TRIM(NumberToVString(maxExtent(2),"*",err,error))// &
            & " is different from the analytic field W parameter of "//TRIM(NumberToVString(W,"*",err,error))//"."
          CALL FlagWarning(localWarning,err,error,*999)
        ENDIF
        IF(ABS(maxExtent(3)-H)>GEOMETRIC_TOLERANCE) THEN
          localWarning="The maximum extent in the z direction of "//TRIM(NumberToVString(maxExtent(3),"*",err,error))// &
            & " is different from the analytic field H parameter of "//TRIM(NumberToVString(H,"*",err,error))//"."
          CALL FlagWarning(localWarning,err,error,*999)
        ENDIF
      CASE DEFAULT
        !Do nothing
      END SELECT
    CASE DEFAULT
      !Do nothing
    END SELECT
    !Now loop over and set boundary conditions and analytic values    
    DO variableIdx=1,numberOfVariables !U & T
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        CALL FieldVariable_ComponentInterpolationCheck(dependentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
          & err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        !Loop over the local nodes excluding the ghosts.
        CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
        DO nodeIdx=1,numberOfNodes
          !Find the position of the node
          CALL FieldVariable_PositionNormalTangentsCalculateNode(dependentVariable,componentIdx,nodeIdx, &
            & position,normal,tangents,err,error,*999)
          x(1:numberOfDimensions)=position(1:numberOfDimensions,NO_PART_DERIV)
          !Loop over the derivatives
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          DO derivativeIdx=1,numberOfNodeDerivatives
            CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
            CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
            DO versionIdx=1,numberOfVersions             
              setBC = .FALSE.            
              SELECT CASE(equationsSetSubType)
              CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
                & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
                SELECT CASE(analyticFunctionType)
                CASE DEFAULT
                  localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
                    & " is invalid for a 2D plane stress or plane strain linear elasticity equations set."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE)     
                SELECT CASE(analyticFunctionType)
                CASE(EQUATIONS_SET_LINEAR_ELASTICITY_CANTILEVER_END_LOAD)
                  !y(x)=2Fx^2(3L-x)/(EWH^3)
                  L=analyticParameters(1)
                  W=analyticParameters(2)
                  H=analyticParameters(3)
                  E=analyticParameters(4)
                  F=analyticParameters(5)
                  Iyy=W*H*H*H/12.0_DP
                  analyticDisplacement=F*x(1)*x(1)*(3.0_DP*L-x(1))/(6.0_DP*E*Iyy)                
                  IF(ABS(x(1))<=GEOMETRIC_TOLERANCE) THEN
                    !Build in end
                    SELECT CASE(variableType)
                    CASE(FIELD_U_VARIABLE_TYPE)
                      analyticValue=0.0_DP
                      setBC=.TRUE.
                    CASE(FIELD_T_VARIABLE_TYPE)
                      IF(ABS(x(3))<=GEOMETRIC_TOLERANCE) THEN
                        !Bottom built-in edge. Reaction force
                        IF(componentIdx==3) THEN
                          analyticValue=-F/REAL(numberOfAxisNodes(3,2),DP)                          
                        ELSE
                          analyticValue=0.0_DP                          
                        ENDIF
                      ELSE
                        analyticValue=0.0_DP
                      ENDIF
                    CASE DEFAULT
                      !Do nothing
                    END SELECT
                  ELSE IF(ABS(x(1)-L)<=GEOMETRIC_TOLERANCE) THEN
                    !Free end
                    IF(ABS(x(3)-H)<=GEOMETRIC_TOLERANCE) THEN
                      !Upper end edge
                      SELECT CASE(variableType)
                      CASE(FIELD_U_VARIABLE_TYPE)
                        IF(componentIdx==3) THEN
                          analyticValue=analyticDisplacement
                        ELSE
                          analyticValue=0.0_DP
                        ENDIF
                      CASE(FIELD_T_VARIABLE_TYPE)
                        IF(componentIdx==3) THEN
                          analyticValue=F/REAL(numberOfAxisNodes(3,2),DP)
                          setBC=.TRUE.
                        ELSE
                          analyticValue=0.0_DP
                        ENDIF
                      CASE DEFAULT
                        !Do nothing
                      END SELECT
                    ELSE
                      IF(variableType==FIELD_U_VARIABLE_TYPE.AND.componentIdx==3) THEN
                        analyticValue=analyticDisplacement
                      ELSE
                        analyticValue=0.0_DP
                      ENDIF
                    ENDIF
                  ELSE
                    IF(variableType==FIELD_U_VARIABLE_TYPE.AND.componentIdx==3) THEN
                      analyticValue=analyticDisplacement
                    ELSE
                      analyticValue=0.0_DP
                    ENDIF
                  ENDIF
                CASE DEFAULT
                  !Do nothing
                END SELECT
              CASE DEFAULT
                !Do nothing
              END SELECT
              CALL FieldVariable_LocalNodeDOFGet(dependentVariable,versionIdx,derivativeIdx,nodeIdx,componentIdx,localDOFIdx, &
                & err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
                & analyticValue,err,error,*999)
              IF(setBC) THEN
                !Default to version 1 of each node derivative
                CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx, &
                  & err,error,*999)
                IF(variableType == FIELD_T_VARIABLE_TYPE) THEN
                  CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx, &
                    & BOUNDARY_CONDITION_NEUMANN_POINT,analyticValue,err,error,*999)
                ELSE                  
                  CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx, &
                    & BOUNDARY_CONDITION_FIXED,analyticValue,err,error,*999)
                ENDIF
              ENDIF
            ENDDO !versionIdx
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    IF(ASSOCIATED(analyticField)) & 
      & CALL FieldVariable_ParameterSetDataRestore(analyticVariable,FIELD_VALUES_SET_TYPE,analyticParameters,err,error,*999)
    
    EXITS("LinearElasticity_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("LinearElasticity_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("LinearElasticity_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE LinearElasticity_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !

  !>Calculated an output field for a linear elasticity equations set.
  SUBROUTINE LinearElasticity_EquationsSetDerivedVariableCalculate(equationsSet,derivedType,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate the output for
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived field type to calculate. \see EquationsSetRoutines_DerivedTypes.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(FieldVariableType), POINTER :: derivedVariable

    ENTERS("LinearElasticity_EquationsSetDerivedVariableCalculate",err,error,*999)
    
    NULLIFY(derivedVariable)
    CALL EquationsSet_DerivedTypeVariableGet(equationsSet,derivedType,derivedVariable,err,error,*999)
    CALL LinearElasticity_StressStrainCalculate(equationsSet,derivedType,derivedVariable,err,error,*999)
   
    EXITS("LinearElasticity_EquationsSetDerivedVariableCalculate")
    RETURN
999 ERRORS("LinearElasticity_EquationsSetDerivedVariableCalculate",err,error)
    EXITS("LinearElasticity_EquationsSetDerivedVariableCalculate")
    RETURN 1
    
  END SUBROUTINE LinearElasticity_EquationsSetDerivedVariableCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a linear elasticity finite element equations set.
  SUBROUTINE LinearElasticity_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_ELEMENT_PARAMETERS=64
    INTEGER(INTG) :: colsComponentIdx,colsVariableType,columnElementDOFIdx,columnElementParameterIdx, &
      & esSpecification(3),gaussPointIdx,numberOfColumnElementParameters(3), &
      & numberOfColsComponents,numberOfDimensions,numberOfGauss,numberOfRowsComponents,numberOfRowElementParameters(3), &
      & numberOfXi,rowsComponentIdx,rowElementParameterIdx,rowElementDOFIdx,rowsVariableType,scalingType, &
      & totalNumberOfColumnElementParameters,totalNumberOfRowElementParameters,voigtIdx
    !INTEGER(INTG) :: offDiagonalComponents(3),offDiagonalDependentVariable(2,2,3),diagonalSubMatrixLocation(3), &
    !  & offDiagonalSubMatLocation(2,3)
    REAL(DP) :: BtCB,BtCMatrix(6,MAX_NUMBER_OF_ELEMENT_PARAMETERS*3),columnBMatrix(6,MAX_NUMBER_OF_ELEMENT_PARAMETERS*3), &
      & gaussWeight,jacobian,jacobianGaussWeight,C(6,6),rowBMatrix(6,MAX_NUMBER_OF_ELEMENT_PARAMETERS*3),sourceParam,thickness
    LOGICAL :: update,updateMatrix,updateRHS,updateSource
    TYPE(BasisType), POINTER :: dependentBasis,geometricBasis
    TYPE(BasisPtrType) :: columnBases(3),rowBases(3)
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,fibreField,geometricField,materialsField,sourceField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,fibreInterpParameters,geometricInterpParameters, &
      & materialsInterpParameters,rowsInterpParameters,sourceInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint,fibreInterpPoint,materialsInterpPoint,sourceInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,geometricVariable,rowsVariable
    TYPE(QuadratureSchemeType), POINTER :: dependentQuadratureScheme,geometricQuadratureScheme
    TYPE(QuadratureSchemePtrType) :: columnQuadratureSchemes(3),rowQuadratureSchemes(3)
    TYPE(VARYING_STRING) :: localError
    !TYPE DPHI_DX_COMP_TYPE !A type to store dPhidX for each mesh component
    !REAL(DP) :: dPhidX(64,3)
    !END TYPE DPHI_DX_COMP_TYPE
    !TYPE(DPHI_DX_COMP_TYPE) :: dPhidXComponent(3)

    ENTERS("LinearElasticity_FiniteElementCalculate",err,error,*999)
    
!!Have a look at XPES40.f in the old CMISS code.
!!Q - CPB: Need to think about anisotropic materials with fibre fields.
!!Q - CPB: why store this dPhidX(columnElementParameterIdx,xiIdx) as opposed to just using it directly? A - to minimize operations - Otherwise it would be calculated many more times than
!!         necessary within the loops below 
!!Q - TPBG: Need to be able to use different Quadrature schemes with different bases? A - No use highest quadrature scheme for all directions
!!TODO:: Check whether quadrature scheme being used is suffient to interpolate highest order basis function    
!!Q - TPBG: Need to be able to use different Interpolation for Geometric & Dependent field?

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE, &
      & EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE, &
      & EQUATIONS_SET_THREE_DIMENSIONAL_ORTHOTROPIC_SUBTYPE)
      !OK
    CASE(EQUATIONS_SET_PLATE_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_SHELL_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a linear elasticity equation type of a elasticty equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(linearMapping)
    CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,1,sourceMapping,err,error,*999)
    ENDIF
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    NULLIFY(equationsMatrix)
    CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,equationsMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(equationsMatrix,updateMatrix,err,error,*999)
    NULLIFY(sourceVectors)
    NULLIFY(sourceVector)
    updateSource=.FALSE.
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
      CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
      CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateSource,err,error,*999)
    ENDIF
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF

    update=(updateMatrix.OR.updateSource.OR.updateRHS)

    IF(update) THEN
      
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
      
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
      NULLIFY(geometricDomain)
      CALL Decomposition_DomainGet(geometricDecomposition,0,geometricDomain,err,error,*999)
      NULLIFY(geometricDomainTopology)
      CALL Domain_DomainTopologyGet(geometricDomain,geometricDomainTopology,err,error,*999)
      NULLIFY(geometricDomainElements)
      CALL DomainTopology_DomainElementsGet(geometricDomainTopology,geometricDomainElements,err,error,*999)
      NULLIFY(geometricBasis)
      CALL DomainElements_ElementBasisGet(geometricDomainElements,elementNumber,geometricBasis,err,error,*999)
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfXi,err,error,*999)
      NULLIFY(geometricQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(geometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,geometricQuadratureScheme,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpParameters,err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters,err,error,*999)

      NULLIFY(fibreField)
      NULLIFY(fibreInterpParameters)
      NULLIFY(fibreInterpPoint)
      CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
      IF(ASSOCIATED(fibreField)) THEN
        CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & fibreInterpParameters,err,error,*999)
        CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & fibreInterpPoint,err,error,*999)
      ENDIF

      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(dependentDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
      NULLIFY(dependentDomainTopology)
      CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
      NULLIFY(dependentDomainElements)
      CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
      NULLIFY(dependentBasis)
      CALL DomainElements_ElementBasisGet(dependentDomainElements,elementNumber,dependentBasis,err,error,*999)
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)

      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(materialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
            
      NULLIFY(sourceField)
      NULLIFY(sourceInterpParameters)
      NULLIFY(sourceInterpPoint)
      CALL EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*999)
      IF(ASSOCIATED(sourceField)) THEN
        CALL EquationsInterpolation_SourceParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & sourceInterpParameters,err,error,*999)
        CALL EquationsInterpolation_SourcePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & sourceInterpPoint,err,error,*999)
      ENDIF
      
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)      
     
      !Create an array of bases with each rows component 
      DO rowsComponentIdx=1,numberOfRowsComponents
        NULLIFY(rowDomain)
        CALL FieldVariable_ComponentDomainGet(rowsVariable,rowsComponentIdx,rowDomain,err,error,*999)
        NULLIFY(rowDomainTopology)
        CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
        NULLIFY(rowDomainElements)
        CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
        NULLIFY(rowBases(rowsComponentIdx)%ptr)
        CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBases(rowsComponentIdx)%ptr,err,error,*999)
        CALL Basis_NumberOfElementParametersGet(rowBases(rowsComponentIdx)%ptr,numberOfRowElementParameters(rowsComponentIdx), &
          & err,error,*999)
        NULLIFY(rowQuadratureSchemes(rowsComponentIdx)%ptr)
        CALL Basis_QuadratureSchemeGet(rowBases(rowsComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
          & rowQuadratureSchemes(rowsComponentIdx)%ptr,err,error,*999)
      ENDDO !rowsComponentIdx
      totalNumberOfRowElementParameters = SUM(numberOfRowElementParameters(1:numberOfRowsComponents))
      IF(totalNumberOfRowElementParameters>MAX_NUMBER_OF_ELEMENT_PARAMETERS) THEN
        localError="The total/sum of the row element parameters of "// &
          & TRIM(NumberToVString(totalNumberOfRowElementParameters,"*",err,error))// &
          & " is greater than the maximum number of element parameters of "// &
          & TRIM(NumberToVString(MAX_NUMBER_OF_ELEMENT_PARAMETERS,"*",err,error))// &
          & ". Increase MAX_NUMBER_OF_ELEMENT_PARAMETERS."
        CALL FlagError(localError,err,error,*999)
      ENDIF

      IF(updateMatrix) THEN
        !Create an array of bases with each component 
        DO colsComponentIdx=1,numberOfColsComponents
          NULLIFY(columnDomain)
          CALL FieldVariable_ComponentDomainGet(colsVariable,colsComponentIdx,columnDomain,err,error,*999)
          NULLIFY(columnDomainTopology)
          CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
          NULLIFY(columnDomainElements)
          CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
          NULLIFY(columnBases(colsComponentIdx)%ptr)
          CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBases(colsComponentIdx)%ptr,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(columnBases(colsComponentIdx)%ptr, &
            & numberOfColumnElementParameters(colsComponentIdx),err,error,*999)
          NULLIFY(columnQuadratureSchemes(colsComponentIdx)%ptr)
          CALL Basis_QuadratureSchemeGet(columnBases(colsComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
            & columnQuadratureSchemes(colsComponentIdx)%ptr,err,error,*999)
        ENDDO !colsComponentIdx
        totalNumberOfColumnElementParameters = SUM(numberOfColumnElementParameters(1:numberOfColsComponents))
        IF(totalNumberOfColumnElementParameters>MAX_NUMBER_OF_ELEMENT_PARAMETERS) THEN
          localError="The total/sum of the column element parameters of "// &
            & TRIM(NumberToVString(totalNumberOfColumnElementParameters,"*",err,error))// &
            & " is greater than the maximum number of element parameters of "// &
            & TRIM(NumberToVString(MAX_NUMBER_OF_ELEMENT_PARAMETERS,"*",err,error))// &
            & ". Increase MAX_NUMBER_OF_ELEMENT_PARAMETERS."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF

      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE, &
        & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
        & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE, &
        & EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE, &
        & EQUATIONS_SET_THREE_DIMENSIONAL_ORTHOTROPIC_SUBTYPE)
        !
        !ONE, TWO & THREE DIMENSIONAL LINEAR ELASTICITY
        !
        !Loop over gauss points & integrate upper triangular portion of Stiffness matrix
        DO gaussPointIdx=1,numberOfGauss !Gauss point index
          
          !Interpolate geometric, fibre and material fields at gauss points & calculate geometric field metrics
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
            & err,error,*999)
          CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
          IF(ASSOCIATED(fibreField)) &
            & CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
            & err,error,*999)
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
            & err,error,*999)
          IF(esSpecification(3)==EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE.OR. &
            & esSpecification(3)==EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE) THEN
            thickness=materialsInterpPoint%values(3,NO_PART_DERIV)
          ENDIF
          IF(ASSOCIATED(sourceField)) THEN
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,sourceInterpPoint, &
              & err,error,*999)
            sourceParam=sourceInterpPoint%values(1,NO_PART_DERIV)
          ENDIF

          !Calculate jacobianGaussWeight.
          CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
          CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
          jacobianGaussWeight=jacobian*gaussWeight

          IF(updateMatrix) THEN
            !Calculate B matrices
            CALL LinearElasticity_StrainMatrixCalculateGauss(numberOfDimensions,numberOfXi, &
              & numberOfRowElementParameters(1:numberOfDimensions),geometricInterpPointMetrics,gaussPointIdx, &
              & rowQuadratureSchemes(1:numberOfDimensions),rowBMatrix,err,error,*999)
            CALL LinearElasticity_StrainMatrixCalculateGauss(numberOfDimensions,numberOfXi, &
              & numberOfColumnElementParameters(1:numberOfDimensions),geometricInterpPointMetrics,gaussPointIdx, &
              & columnQuadratureSchemes(1:numberOfDimensions),columnBMatrix,err,error,*999)
          
            !Create Linear Elasticity Tensor C
            CALL LinearElasticity_ElasticityTensor(esSpecification(3),numberOfDimensions,materialsInterpPoint,C,err,error,*999)
            !Compute B^T.C
            DO rowElementDOFIdx=1,totalNumberOfRowElementParameters
              DO voigtIdx=1,NUMBER_OF_VOIGT(numberOfDimensions)
                BtCMatrix(voigtIdx,rowElementDOFIdx)= &
                  & DOT_PRODUCT(rowBMatrix(1:NUMBER_OF_VOIGT(numberOfDimensions),rowElementDOFIdx), &
                  & C(1:NUMBER_OF_VOIGT(numberOfDimensions),voigtIdx))
              ENDDO !voigtIdx
            ENDDO !rowElementDOFIdx
          ENDIF

!!TODO: EXPLOIT STIFFNESS MATRIX SYMMETRY
          
          rowElementDOFIdx=0
          DO rowsComponentIdx=1,numberOfRowsComponents
            DO rowElementParameterIdx=1,numberOfRowElementParameters(rowsComponentIdx)
              rowElementDOFIdx=rowElementDOFIdx+1
              IF(updateMatrix) THEN
                columnElementDOFIdx=0
                DO colsComponentIdx=1,numberOfColsComponents
                  DO columnElementParameterIdx=1,numberOfColumnElementParameters(colsComponentIdx)
                    columnElementDOFIdx=columnElementDOFIdx+1
                    !Compute (B^T.C).B
                    BtCB=DOT_PRODUCT(BtCMatrix(1:NUMBER_OF_VOIGT(numberOfDimensions),rowElementDOFIdx), &
                      & columnBMatrix(1:NUMBER_OF_VOIGT(numberOfDimensions),columnElementDOFIdx))
                    equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+BtCB*jacobianGaussWeight
                  ENDDO !columnElementParameterIdx
                ENDDO !colsComponentIdx
              ENDIF !updateMatrix
              IF(updateRHS) THEN
                rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
              ENDIF !updateRHS
            ENDDO !rowElementParameterIdx
          ENDDO !rowsComponentIdx                     
          IF(esSpecification(3)==EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE.OR. &
            & esSpecification(3)==EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE) THEN
            equationsMatrix%elementMatrix%matrix=equationsMatrix%elementMatrix%matrix*thickness
          ENDIF
        ENDDO !gaussPointIdx
        
        !Scale factor adjustment
        CALL Field_ScalingTypeGet(dependentField,scalingType,err,error,*999)
        IF(scalingType/=FIELD_NO_SCALING) THEN
          NULLIFY(rowsInterpParameters)
          CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,rowsVariableType,rowsInterpParameters, &
            & err,error,*999)
          NULLIFY(colsInterpParameters)
          CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,colsVariableType,colsInterpParameters, &
            & err,error,*999)
          CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,rowsInterpParameters,err,error,*999)
          CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,colsInterpParameters,err,error,*999)
          !Loop over element rows
          rowElementDOFIdx=0          
          DO rowsComponentIdx=1,numberOfRowsComponents
            DO rowElementParameterIdx=1,numberOfRowElementParameters(rowsComponentIdx)
              rowElementDOFIdx=rowElementDOFIdx+1                    
              IF(updateMatrix) THEN
                !Loop over element columns
                columnElementDOFIdx=0
                DO colsComponentIdx=1,numberOfColsComponents
                  DO columnElementParameterIdx=1,numberOfColumnElementParameters(colsComponentIdx)
                    columnElementDOFIdx=columnElementDOFIdx+1
                    equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                      & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowsComponentIdx)* &
                      & colsInterpParameters%scaleFactors(columnElementParameterIdx,colsComponentIdx)                  
                  ENDDO !columnElementParameterIdx
                ENDDO !colsComponentIdx
              ENDIF !update matrix
              IF(updateRHS) THEN
                rhsVector%elementVector%vector(rowElementDOFIdx)=rhsVector%elementVector%vector(rowElementDOFIdx)* &
                  & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowsComponentIdx)
              ENDIF
            ENDDO !rowElementParameterIdx
          ENDDO !rowsComponentIdx
        ENDIF !scaling
      
!       !Loop over gauss points & integrate upper triangular portion of Stiffness matrix
      
!       DO gaussPointIdx=1,numberOfGauss !Gauss point index

!           !Parameters for number of off diagonal stress/strain terms for a given number of xi directions and order of
!           !calculation for shear terms
!           !These parameters do not change for 1D,2D,3D Linear Elasticity
!           offDiagonalComponents = [0,1,3]
!           offDiagonalDependentVariable(1,1,:) = [1,1,2]
!           offDiagonalDependentVariable(1,2,:) = [2,3,3]
!           offDiagonalDependentVariable(2,1,:) = offDiagonalDependentVariable(1,2,:)
!           offDiagonalDependentVariable(2,2,:) = offDiagonalDependentVariable(1,1,:)
!           !
!           diagonalSubMatrixLocation(:) = [0,numberDependentElementParameters(1),SUM(numberDependentElementParameters(1:2))]
!           offDiagonalSubMatLocation(1,:) = [0,0,numberDependentElementParameters(1)]
!           offDiagonalSubMatLocation(2,:) = [numberDependentElementParameters(1),diagonalSubMatrixLocation(3), &
!             & diagonalSubMatrixLocation(3)]

!           !Interpolate geometric, fibre and material fields at gauss points & calculate geometric field metrics
!           CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
!             & err,error,*999)
! !!TODO:: Add option to only evaluate required metrics
!           CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
!           IF(ASSOCIATED(fibreField)) &
!             & CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
!             & err,error,*999)
!           CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
!             & err,error,*999)
!           IF(ASSOCIATED(sourceField)) THEN
!             CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,sourceInterpPoint, &
!               & err,error,*999)
!             sourceParam=sourceInterpPoint%values(1,NO_PART_DERIV)
!           ENDIF

!           !Calculate jacobianGaussWeight.
!           CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
!           CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
!           jacobianGaussWeight=jacobian*gaussWeight
          
!           DO colsComponentIdx=1,numberOfColsComponents
!             dPhidXComponent(colsComponentIdx)%dPhidX = 0.0_DP
!             DO columnElementParameterIdx=1,numberDependentElementParameters(colsComponentIdx)
!               DO columnXiIdx=1,numberOfXi
!                 CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureSchemes(colsComponentIdx)%ptr, &
!                   & columnElementParameterIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx, &
!                   & colsdPhidXi,err,error,*999)
!                 DO rowXiIdx=1,numberOfXi
!                   !!TODO: second index in dXidX should be component not xi?
!                   dPhidXComponent(colsComponentIdx)%dPhidX(columnElementParameterIdx,rowXiIdx) = &
!                     & dPhidXComponent(columnXiIdx)%dPhidX(columnElementParameterIdx,rowXiIdx)+ &
!                     & geometricInterpPointMetrics%dXidX(rowXiIdx,columnXiIdx)*colsdPhidXi
!                 ENDDO !rowXiIdx
!               ENDDO !columnXiIdx
!             ENDDO !columnElementParameterIdx
!           ENDDO !xiIdx
!           !TODO:: what about fibres?
!           !Create Linear Elasticity Tensor C
!           CALL LinearElasticity_ElasticityTensor(esSpecification(3),numberOfDimensions,materialsInterpPoint,C,err,error,*999)
!           !Store Elasticity Tensor diagonal & off diagonal stress coefficients
!           jacobianGaussWeightDiagC(3,:) = [jacobianGaussWeight*C(4,4),jacobianGaussWeight*C(5,5),jacobianGaussWeight*C(3,3)]
!           jacobianGaussWeightDiagC(2,:) = [jacobianGaussWeight*C(6,6),jacobianGaussWeight*C(2,2),jacobianGaussWeightDiagC(3,2)]
!           jacobianGaussWeightDiagC(1,:) = [jacobianGaussWeight*C(1,1),jacobianGaussWeightDiagC(2,1),jacobianGaussWeightDiagC(3,1)]
!           jacobianGaussWeightOffDiagC(1,:) = [jacobianGaussWeight*C(1,2),jacobianGaussWeight*C(1,3),jacobianGaussWeight*C(2,3)]
!           jacobianGaussWeightOffDiagC(2,:) = [jacobianGaussWeight*C(6,6),jacobianGaussWeight*C(4,4),jacobianGaussWeight*C(5,5)]
!           !Construct Element Matrix Diagonal Terms
!           DO xiIdx=1,numberOfXi
!             DO columnElementParameterIdx=1,numberDependentElementParameters(xiIdx)
!               DO rowElementParameterIdx=columnElementParameterIdx,numberDependentElementParameters(xiIdx)
!                 equationsMatrix%elementMatrix%matrix(diagonalSubMatrixLocation(xiIdx)+ &
!                   & columnElementParameterIdx,diagonalSubMatrixLocation(xiIdx)+rowElementParameterIdx) = &
!                   & equationsMatrix%elementMatrix%matrix(diagonalSubMatrixLocation(xiIdx)+ &
!                   & columnElementParameterIdx,diagonalSubMatrixLocation(xiIdx)+rowElementParameterIdx) + &
!                   & DOT_PRODUCT(dPhidXComponent(xiIdx)%dPhidX(columnElementParameterIdx,1:numberOfXi)* &
!                   & dPhidXComponent(xiIdx)%dPhidX(rowElementParameterIdx,1:numberOfXi), &
!                   & jacobianGaussWeightDiagC(xiIdx,1:numberOfXi))
!               ENDDO !rowElementParameterIdx
!             ENDDO !columnElementParameterIdx
!           ENDDO !xiIdx
!           !Construct Element Matrix Off-Diagonal Terms
!           DO xiIdx=1,offDiagonalComponents(numberOfXi)
!             DO columnElementParameterIdx=1,numberDependentElementParameters(offDiagonalDependentVariable(1,1,xiIdx))
!               DO rowElementParameterIdx=1,numberDependentElementParameters(offDiagonalDependentVariable(1,2,xiIdx))
!                 equationsMatrix%elementMatrix%matrix(offDiagonalSubMatLocation(1,xiIdx)+ &
!                   & columnElementParameterIdx,offDiagonalSubMatLocation(2,xiIdx)+rowElementParameterIdx) = &
!                   & equationsMatrix%elementMatrix%matrix(offDiagonalSubMatLocation(1,xiIdx)+ &
!                   & columnElementParameterIdx,offDiagonalSubMatLocation(2,xiIdx)+rowElementParameterIdx)+ &
!                   & DOT_PRODUCT(dPhidXComponent(offDiagonalDependentVariable(1,1,xiIdx))% &
!                   & dPhidX(columnElementParameterIdx,offDiagonalDependentVariable(1,:,xiIdx))* &
!                   & dPhidXComponent(offDiagonalDependentVariable(1,2,xiIdx))% &
!                   & dPhidX(rowElementParameterIdx,offDiagonalDependentVariable(2,:,xiIdx)), &
!                   & jacobianGaussWeightOffDiagC(:,xiIdx))
!               ENDDO !rowElementParameterIdx
!             ENDDO !columnElementParameterIdx
!           ENDDO !xiIdx
          
!           !Below is the full form of constructing the off-Diagonal terms. This will be documented in the linear elasticity
!           !equation set page on doxygen for clarity
                        
!           !Expanding the DOT_PRODUCT terms
          
!           ! offDiagonalDependentVariable(1,:) = [1,1,2]
!           ! offDiagonalDependentVariable(2,:) = [2,3,3]
!           ! DO xiIdx=1,offDiagonalComponents(numberOfXi)
!           !   DO columnElementParameterIdx=1,numberDependentElementParameters(offDiagonalDependentVariable(1,1,xiIdx))
!           !     DO rowElementParameterIdx=1,numberDependentElementParameters(offDiagonalDependentVariable(1,2,xiIdx))
!           !       equatiocolumnElementParameterIdxMatrix%elementMatrix%matrix(offDiagonalSubMatLocation(1,xiIdx)+ &
!           !         & columnElementParameterIdx,offDiagonalSubMatLocation(2,xiIdx)+rowElementParameterIdx) = &
!           !         & equationsMatrix%elementMatrix%matrix(offDiagonalSubMatLocation(1,xiIdx)+ &
!           !         & columnElementParameterIdx,offDiagonalSubMatLocation(2,xiIdx)+rowElementParameterIdx)+ &
!           !         & jacobianGaussWeightOffDiagC(1,xiIdx)*dPhidXComponent(offDiagonalDependentVariable(1,xiIdx))% &
!           !         & dPhidX(columnElementParameterIdx,offDiagonalDependentVariable(1,xiIdx))* &
!           !         & dPhidXComponent(offDiagonalDependentVariable(2,xiIdx))% &
!           !         & dPhidX(rowElementParameterIdx,offDiagonalDependentVariable(2,xiIdx))+ &
!           !         & jacobianGaussWeightOffDiagC(2,xiIdx)*dPhidXComponent(offDiagonalDependentVariable(1,xiIdx))% &
!           !         & dPhidX(columnElementParameterIdx,offDiagonalDependentVariable(2,xiIdx))* &
!           !         & dPhidXComponent(offDiagonalDependentVariable(2,xiIdx))% &
!           !         & dPhidX(rowElementParameterIdx,offDiagonalDependentVariable(1,xiIdx))
!           !     ENDDO !rowElementParameterIdx
!           !   ENDDO !columnElementParameterIdx
!           ! ENDDO !xiIdx
          
!           ! !Expanding the xiIdx loop above
          
!           ! DO columnElementParameterIdx=1,numberDependentElementParameters(1)
!           !   DO rowElementParameterIdx=1,numberDependentElementParameters(2)
!           !     equationsMatrix%elementMatrix%matrix(columnElementParameterIdx,rowElementParameterIdx+ &
!           !       & numberDependentElementParameters(1))= &
!           !       & equationsMatrix%elementMatrix%matrix(columnElementParameterIdx,rowElementParameterIdx+ &
!           !       & numberDependentElementParameters(1))+ &
!           !       & jacobianGaussWeight*C(1,2)*dPhidXComponent(1)%dPhidX(columnElementParameterIdx,1)* &
!           !       & dPhidXComponent(2)%dPhidX(rowElementParameterIdx,2)+ &
!           !       & jacobianGaussWeight*C(6,6)*dPhidXComponent(1)%dPhidX(columnElementParameterIdx,2)* &
!           !       & dPhidXComponent(2)%dPhidX(rowElementParameterIdx,1)
!           !   ENDDO !columnElementParameterIdx
!           ! ENDDO !rowElementParameterIdx
!           ! DO columnElementParameterIdx=1,numberDependentElementParameters(1)
!           !   DO rowElementParameterIdx=1,numberDependentElementParameters(3)
!           !     equationsMatrix%elementMatrix%matrix(columnElementParameterIdx,rowElementParameterIdx+ &
!           !       & numberDependentElementParameters(1)+numberDependentElementParameters(2))  &
!           !       & equationsMatrix%elementMatrix%matrix(columnElementParameterIdx,rowElementParameterIdx+ &
!           !       & numberDependentElementParameters(1)+numberDependentElementParameters(2)) + &
!           !       & jacobianGaussWeight*C(1,3)*dPhidXComponent(1)%dPhidX(columnElementParameterIdx,1)* &
!           !       & dPhidXComponent(3)%dPhidX(rowElementParameterIdx,3)+ &
!           !       & jacobianGaussWeight*C(4,4)*dPhidXComponent(1)%dPhidX(columnElementParameterIdx,3)* &
!           !       & dPhidXComponent(3)%dPhidX(rowElementParameterIdx,1)
!           !   ENDDO !columnElementParameterIdx
!           ! ENDDO !rowElementParameterIdx
!           ! DO columnElementParameterIdx=1,numberDependentElementParameters(2)
!           !   DO rowElementParameterIdx=1,numberDependentElementParameters(3)
!           !     equationsMatrix%elementMatrix%matrix(columnElementParameterIdx+ &
!           !       & numberDependentElementParameters(1),rowElementParameterIdx+numberDependentElementParameters(1)+ &
!           !       & numberDependentElementParameters(2))= &
!           !       & equationsMatrix%elementMatrix%matrix(columnElementParameterIdx+ &
!           !       & numberDependentElementParameters(1),rowElementParameterIdx+numberDependentElementParameters(1)+ &
!           !       & numberDependentElementParameters(2))+ &
!           !       & jacobianGaussWeight*C(2,3)*dPhidXComponent(2)%dPhidX(columnElementParameterIdx,2)* &
!           !       & dPhidXComponent(3)%dPhidX(rowElementParameterIdx,3)+ &
!           !       & jacobianGaussWeight*C(5,5)*dPhidXComponent(2)%dPhidX(columnElementParameterIdx,3)* &
!           !       & dPhidXComponent(3)%dPhidX(rowElementParameterIdx,2)
!           !   ENDDO !columnElementParameterIdx
!           ! ENDDO !rowElementParameterIdx
          
!         ENDDO !gaussPointIdx
        
!         !!If Plane Stress/Strain problem multiply equation matrix by thickness
!         !IF(esSpecification(3) == EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE .OR. &
!         !  & esSpecification(3) == EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE .OR. & 
!         !  & esSpecification(3) == EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE) THEN
!         !  DO rowElementDOFIdx=1,totalNumberDependentElementParameters
!         !    DO columnElementDOFIdx=rowElementDOFIdx,totalNumberDependentElementParameters
! !!TODO::Bring 2D plane stress/strain element thickness in through a field - element constant when it can be exported by field i/o.
! !!      Currently brought in through material field (Temporary)
              
!         !      equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
!         !        & equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
!         !        & materialsInterpPoint%values(1,NO_PART_DERIV)
!         !    ENDDO !columnElementDOFIdx
!         !  ENDDO !rowElementDOFIdx
!         !ENDIF
      
! !!TODO:: Is this RHS Vector update required? find out/check - RHS not used - BC are prescribed during assembling
! !!       eg update RHS only when BC change - stiffness matrix should be the same
!         IF(updateRHS) THEN
!           rhsVector%elementVector%vector=0.0_DP
!         ENDIF
        
!         !Scale factor adjustment, Application of scale factors is symmetric
!         CALL Field_ScalingTypeGet(dependentField,scalingType,err,error,*999)
!         IF(scalingType/=FIELD_NO_SCALING) THEN
!           NULLIFY(rowsInterpParameters)
!           CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,rowsVariableType,rowsInterpParameters, &
!           & err,error,*999)
!           NULLIFY(colsInterpParameters)
!           CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,colsVariableType,colsInterpParameters, &
!             & err,error,*999)
!           CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,rowsInterpParameters,err,error,*999)
!           CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,colsInterpParameters,err,error,*999)
!           DO xiIdx=1,numberOfXi
!             rowsSF(diagonalSubMatrixLocation(xiIdx)+1:SUM(numberDependentElementParameters(1:xiIdx)))= &
!               & rowsInterpParameters%scaleFactors(:,xiIdx)
!             colsSF(diagonalSubMatrixLocation(xiIdx)+1:SUM(numberDependentElementParameters(1:xiIdx)))= &
!               & colsInterpParameters%scaleFactors(:,xiIdx)
!           ENDDO !xiIdx
!           DO rowElementDOFIdx=1,totalNumberDependentElementParameters
!             IF(updateMatrix) THEN
!               DO columnElementDOFIdx=rowElementDOFIdx,totalNumberDependentElementParameters
!                 equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
!                 & equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
!                 & rowsSF(rowElementDOFIdx)*colsSF(columnElementDOFIdx)
!               ENDDO !columnElementDOFIdx
!             ENDIF
! !!TODO:: Check if RHS update required for Linear Elasticity ie is the RHS the force terms but they are set during assembling and not here?
!             IF(updateRHS) THEN
!               rhsVector%elementVector%vector(rowElementDOFIdx)= &
!                 & rhsVector%elementVector%vector(rowElementDOFIdx)*rowsSF(rowElementDOFIdx)
!             ENDIF
!           ENDDO !rowElementDOFIdx
          
!           IF(updateMatrix) THEN
!             !Transpose upper triangular portion of Stiffness matrix to give lower triangular portion.
!             !Has to be done after scale factors are applied
! !!TODO:: Use symmetric linear equation solver or alternatively traspose to give full matrix when asemmbling or when
! !!       creating solver matrices        
! !!TODO:: Better to use SIZE(equationsMatrix%elementMatrix%matrix,1) as apposed to totalNumberDependentElementParameters?
! !!       Is the size re-calculated at end of every loop?
!             DO rowElementDOFIdx=2,totalNumberDependentElementParameters
!               DO columnElementDOFIdx=1,rowElementDOFIdx-1
!                 equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
!                   &  equationsMatrix%elementMatrix%matrix(columnElementDOFIdx,rowElementDOFIdx)
!               ENDDO !columnElementDOFIdx
!             ENDDO !rowElementDOFIdx
!           ENDIF
!         ENDIF !scaling
        
      CASE(EQUATIONS_SET_PLATE_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_SHELL_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a linear elasticity equation type of a elasticty equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
    ENDIF !update

    EXITS("LinearElasticity_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("LinearElasticity_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE LinearElasticity_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the linear elasticity tensor
  SUBROUTINE LinearElasticity_ElasticityTensor(equationsSetSubtype,numberOfDimensions,materialsInterpolatedPoint, &
    & elasticityTensor,err,error,*)

    !Argument variables    
    INTEGER(INTG), INTENT(IN) :: equationsSetSubtype !<The subtype of the particular equation set being used
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions for the equations set
    TYPE(FieldInterpolatedPointType), POINTER :: materialsInterpolatedPoint
    REAL(DP), INTENT(OUT) :: elasticityTensor(:,:) !<The Linear Elasticity Tensor C
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: delta,E,E1,E2,E3,v,v13,v23,v12,v31,v32,v21
    REAL(DP) :: G12,G13,G23
    TYPE(VARYING_STRING) :: localError

    ENTERS("LinearElasticity_ElasticityTensor",err,error,*999)
    
    elasticityTensor(1:NUMBER_OF_VOIGT(numberOfDimensions),1:NUMBER_OF_VOIGT(numberOfDimensions))=0.0_DP
    SELECT CASE(equationsSetSubtype)
    CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE)
      !Isotropic 3D linear elasticity tensor, 2 independent elastic constants
      E = materialsInterpolatedPoint%values(1,1)
      v = materialsInterpolatedPoint%values(2,1)
      delta = E/((1.0_DP+v)*(1.0_DP-2.0_DP*v))
      elasticityTensor(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,1)) = (1.0_DP-v)*delta
      elasticityTensor(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,2)) = v*delta
      elasticityTensor(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,3)) = v*delta 
      elasticityTensor(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,1)) = v*delta
      elasticityTensor(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,2)) = (1.0_DP-v)*delta
      elasticityTensor(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,3)) = v*delta
      elasticityTensor(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,1)) = v*delta
      elasticityTensor(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,2)) = v*delta
      elasticityTensor(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,3)) = (1.0_DP-v)*delta
      elasticityTensor(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,3)) = (1.0_DP-2.0_DP*v)*delta/2.0_DP
      elasticityTensor(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,3)) = (1.0_DP-2.0_DP*v)*delta/2.0_DP
      elasticityTensor(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,2)) = (1.0_DP-2.0_DP*v)*delta/2.0_DP
    CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ORTHOTROPIC_SUBTYPE)
      !General orthotropic 3D linear elasticity tensor, 9 independent elastic constants
      E1 = materialsInterpolatedPoint%values(1,1)
      E2 = materialsInterpolatedPoint%values(2,1)
      E3 = materialsInterpolatedPoint%values(3,1)
      v23 = materialsInterpolatedPoint%values(4,1)
      v13 = materialsInterpolatedPoint%values(5,1)
      v12 = materialsInterpolatedPoint%values(6,1)
      G23 = materialsInterpolatedPoint%values(7,1)
      G13 = materialsInterpolatedPoint%values(8,1)
      G12 = materialsInterpolatedPoint%values(9,1)
      v32 = v23*E3/E2
      v31 = v13*E3/E1
      v21 = v12*E2/E1
      delta = 1.0_DP-v12*v21-v23*v32-v13*v31-v12*v23*v31-v21*v32*v13
      elasticityTensor(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,1)) = E1*(1.0_DP-v23*v32)/delta
      elasticityTensor(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,2)) = E1*(v21+v23*v31)/delta
      elasticityTensor(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,3)) = E1*(v31+v21*v32)/delta 
      elasticityTensor(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,1)) = E2*(v12+v13*v32)/delta
      elasticityTensor(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,2)) = E2*(1.0_DP-v13*v31)/delta
      elasticityTensor(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,3)) = E2*(v32+v12*v31)/delta
      elasticityTensor(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,1)) = E3*(v13+v12*v23)/delta
      elasticityTensor(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,2)) = E3*(v23+v13*v21)/delta
      elasticityTensor(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,3)) = E3*(1.0_DP-v12*v21)/delta
      elasticityTensor(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,3)) = G23
      elasticityTensor(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,3)) = G13 
      elasticityTensor(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,2)) = G12
    CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE)
      !Plane stress 2D isotropic elasticity tensor, 2 independent elastic constants
      E = materialsInterpolatedPoint%values(1,1)
      v = materialsInterpolatedPoint%values(2,1)
      delta = E/(1.0_DP-v*v)
      elasticityTensor(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,1)) = delta
      elasticityTensor(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,2)) = v*delta
      elasticityTensor(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,1)) = v*delta
      elasticityTensor(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,2)) = delta
      elasticityTensor(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,2)) = (1.0_DP-v)*delta/2.0_DP
    CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
      !Plane strain 2D isotropic linear elasticity tensor, 2 independent elastic constants
      E = materialsInterpolatedPoint%values(1,1)
      v = materialsInterpolatedPoint%values(2,1)
      delta = E/((1.0_DP+v)*(1.0_DP-2.0_DP*v))
      elasticityTensor(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,1)) = (1.0_DP-v)*delta
      elasticityTensor(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,2)) = v*delta
      elasticityTensor(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,1)) = v*delta
      elasticityTensor(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,2)) = (1.0_DP-v)*delta
      elasticityTensor(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,2)) = (1.0_DP-2.0_DP*v)*delta/2.0_DP
    CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
      !1D linear elasticity tensor
      E = materialsInterpolatedPoint%values(2,1)
      elasticityTensor(TENSOR_TO_VOIGT1(1,1),TENSOR_TO_VOIGT1(1,1)) = E
    CASE(EQUATIONS_SET_PLATE_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_SHELL_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(equationsSetSubtype,"*",err,error))// &
            & " is not valid for a linear elasticity equation type of a elasticty equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      SELECT CASE(equationsSetSubtype)
      CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE)
       CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Elasticity tensor: 3D Isotropic",err,error,*999)
      CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ORTHOTROPIC_SUBTYPE)
       CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Elasticity tensor: 3D Orthotropic",err,error,*999)
      CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE)
       CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Elasticity tensor: 2D plane stress",err,error,*999)
      CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
       CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Elasticity tensor: 3D plane strain",err,error,*999)
      CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
       CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Elasticity tensor: 1D",err,error,*999)
      CASE(EQUATIONS_SET_PLATE_SUBTYPE)
       CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Elasticity tensor: plate",err,error,*999)
      CASE(EQUATIONS_SET_SHELL_SUBTYPE)
       CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Elasticity tensor: shell",err,error,*999)
      CASE DEFAULT
       CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Elasticity tensor: unknown",err,error,*999)
      END SELECT
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_VOIGT(numberOfDimensions),1,1, &
        & NUMBER_OF_VOIGT(numberOfDimensions),8,8,elasticityTensor(1:NUMBER_OF_VOIGT(numberOfDimensions), &
        & 1:NUMBER_OF_VOIGT(numberOfDimensions)),WRITE_STRING_MATRIX_NAME_AND_INDICES, &
        & '("  C','(",I1,",:)',':",8(X,E13.6))','9X,8(X,E13.6))',err,error,*999)
    ENDIF

    EXITS("LinearElasticity_ElasticityTensor")
    RETURN
999 ERRORSEXITS("LinearElasticity_ElasticityTensor",err,error)
    RETURN 1
    
  END SUBROUTINE LinearElasticity_ElasticityTensor

  !
  !================================================================================================================================
  !

  !>Sets up the Linear elasticity equation type of an elasticity equations set class.
  SUBROUTINE LinearElasticity_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup a linear elasticity equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,derivedIdx,esSpecification(3),geometricComponentNumber,geometricMeshComponent, &
      & geometricScalingType,numberOfAnalyticComponents,numberOfComponents,numberOfDimensions,numberOfTensorComponents, &
      & solutionMethod,sparsityType,variableIdx,variableType
    INTEGER(INTG), ALLOCATABLE :: variableTypes(:)
    REAL(DP) :: E
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsSetDerivedType), POINTER :: equationsDerived
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: analyticField,geometricField,materialsField
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("LinearElasticity_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE, &
      & EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE, &
      & EQUATIONS_SET_THREE_DIMENSIONAL_ORTHOTROPIC_SUBTYPE)
      !OK
    CASE(EQUATIONS_SET_PLATE_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_SHELL_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a linear elasticity equation type of an elasticity equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)

    SELECT CASE(equationsSetSetup%setupType)
    CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
      !-----------------------------------------------------------------
      ! I n i t i a l   s e t u p
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        !Default to FEM solution
        CALL LinearElasticity_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !-----------------------------------------------------------------
      ! G e o m e t r i c   f i e l d
      !-----------------------------------------------------------------
      !Do nothing???
!!TODO: Check 1, 2, 3 D for the various subtypes
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          !Create the auto created dependent field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
          NULLIFY(geometricDecomposition)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,2,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_T_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_T_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_T_VARIABLE_TYPE,FIELD_DP_TYPE, &
            & err,error,*999)
          numberOfComponents=numberOfDimensions
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & numberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_T_VARIABLE_TYPE, &
            & numberOfComponents,err,error,*999)
          !Default to the geometric interpolation setup
          DO componentIdx=1,numberOfDimensions
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
              & err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_T_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_T_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            !Default the scaling to the geometric field scaling
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsSet%dependent%dependentField,geometricScalingType,err,error,*999)
          CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_T_VARIABLE_TYPE], &
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_T_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_T_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          numberOfComponents=numberOfDimensions
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_T_VARIABLE_TYPE,numberOfComponents, &
            & err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfDimensions
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_T_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
          CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity equation"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
      !-----------------------------------------------------------------
      ! M a t e r i a l s   f i e l d 
      !-----------------------------------------------------------------
      NULLIFY(equationsMaterials)
      CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
        numberOfComponents=1
      CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
        & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
        numberOfComponents=3
      CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE)
        numberOfComponents=2
      CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ORTHOTROPIC_SUBTYPE)
        numberOfComponents=9
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a linear elasticity equation type of an elasticity equation set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Create the auto created materials field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          NULLIFY(geometricDecomposition)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,numberOfComponents, &
            & err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
          DO componentIdx=1,numberOfComponents
            !Default to to the first geometric component with constant interpolation
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          ENDDO !componentIdx
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values for the materials field
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
            !1 independent elastic constant
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,1,30.0E6_DP,err,error,*999) !Young's Modulus
          CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
            & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
            !2 independent elastic constants+thickness
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,1,30.0E6_DP,err,error,*999) !Young's Modulus
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,2,0.25_DP,err,error,*999) !Poisson's Ratio
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,3,1.0_DP,err,error,*999) !thickness           
          CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE)
            !2 independent elastic constants
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,1,30.0E6_DP,err,error,*999) !Young's Modulus
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,2,0.25_DP,err,error,*999) !Poisson's Ratio
          CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ORTHOTROPIC_SUBTYPE)
            !9 independent elastic constants
            DO componentIdx=1,3
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & componentIdx,30.0E6_DP,err,error,*999) !Youngs's modulus
            ENDDO !componentIdx
            DO componentIdx=4,6
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & componentIdx,0.25_DP,err,error,*999) !Poisson's ratio
            ENDDO !componentIdx
            DO componentIdx=7,9
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & componentIdx,30.0E6_DP,err,error,*999) !Shear modulus
            ENDDO !componentIdx
          CASE DEFAULT
            localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is not valid for a linear elasticity equation type of an elasticity equation set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      !-----------------------------------------------------------------
      ! S o u r c e   f i e l d 
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        !Do nothing
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      !-----------------------------------------------------------------
      ! A n a l y t i c   f i e l d
      !-----------------------------------------------------------------
      numberOfAnalyticComponents=0
      NULLIFY(equationsAnalytic)
      CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
          & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
          SELECT CASE(equationsAnalytic%analyticFunctionType)
          CASE DEFAULT
            !Do nothing
          END SELECT
        CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE)
          SELECT CASE(equationsSetSetup%analyticFunctionType)
          CASE(EQUATIONS_SET_LINEAR_ELASTICITY_CANTILEVER_END_LOAD)
          !Check that we are in 3D
          IF(numberOfDimensions/=3) THEN
            localError="The number of geometric dimensions of "// &
              & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid. The analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " requires that there be 3 geometric dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Check the materials values are constant
          NULLIFY(materialsField)
          CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
          !E
          CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,1,FIELD_CONSTANT_INTERPOLATION, &
            & err,error,*999)
          !nu
          CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,2,FIELD_CONSTANT_INTERPOLATION, &
            & err,error,*999)
          !Set the number of analytic components
          numberOfAnalyticComponents=5 !L,H,W,E, and F
          !Set analtyic function type
          equationsAnalytic%analyticFunctionType=EQUATIONS_SET_LINEAR_ELASTICITY_CANTILEVER_END_LOAD            
          CASE DEFAULT
            !Do nothing
          END SELECT
        CASE DEFAULT
          !Do nothing
        END SELECT
        !Set temporal nature
        CALL EquationsSet_AnalyticIsTemporalSet(equationsSet,.FALSE.,err,error,*999)
        !Create analytic field if required
        IF(numberOfAnalyticComponents>1) THEN
          IF(equationsAnalytic%analyticFieldAutoCreated) THEN
            !Create the auto created analytic field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsAnalytic%analyticField,err,error,*999)
            CALL Field_LabelSet(equationsAnalytic%analyticField,"Analytic Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsAnalytic%analyticField,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsAnalytic%analyticField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsAnalytic%analyticField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsAnalytic%analyticField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsAnalytic%analyticField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL Field_VariableLabelSet(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,"Analytic",err,error,*999)
            IF(numberOfAnalyticComponents==1) THEN
              CALL Field_DimensionSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
            ELSE
              CALL Field_DimensionSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
            ENDIF
            CALL Field_DataTypeSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            !Set the number of analytic components
            CALL Field_NumberOfComponentsSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
              & numberOfAnalyticComponents,err,error,*999)
            !Default the analytic components to the 1st geometric interpolation setup with constant interpolation
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            DO componentIdx=1,numberOfAnalyticComponents
              CALL Field_ComponentMeshComponentSet(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentInterpolationSet(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsAnalytic%analyticField,geometricScalingType,err,error,*999)
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            IF(numberOfAnalyticComponents==1) THEN
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            ELSE
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            ENDIF
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfAnalyticComponents, &
              & err,error,*999)
          ENDIF
        ENDIF 
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        NULLIFY(analyticField)
        CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
        IF(ASSOCIATED(analyticField)) THEN
          IF(equationsAnalytic%analyticFieldAutoCreated) THEN
            CALL Field_CreateFinish(analyticField,err,error,*999)
            !Set the default values for the analytic field
            SELECT CASE(esSpecification(3))
            CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
              & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
              SELECT CASE(equationsAnalytic%analyticFunctionType)
              CASE DEFAULT
                !Do nothing
              END SELECT
            CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE)
              SELECT CASE(equationsAnalytic%analyticFunctionType)
              CASE(EQUATIONS_SET_LINEAR_ELASTICITY_CANTILEVER_END_LOAD)
                !Set number of analytic field components
                numberOfAnalyticComponents=5
                !L
                CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,10.0_DP,err,error,*999)
                !W
                CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,2,2.0_DP,err,error,*999)
                !H
                CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,3,2.0_DP,err,error,*999)
                !E
                !Default the Youngs modulus values to the same as the materials field
                NULLIFY(materialsField)
                CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
                CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,E, &
                  & err,error,*999)
                CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,4,E,err,error,*999)
                !F
                CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,5,-100.0_DP,err,error,*999)
               CASE DEFAULT
                !Do nothing
              END SELECT
            CASE DEFAULT
              !Do nothing
            END SELECT
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Linear Elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s 
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
        !Create the equations
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          !Finish the equations creation
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_CreateFinish(equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          !Create the equations mapping.
          NULLIFY(vectorMapping)
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
          CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
            & err,error,*999)
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_T_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMappingVector_RHSCoefficientSet(vectorMapping,-1.0_DP,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          NULLIFY(vectorMatrices)
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
          SELECT CASE(sparsityType)
          CASE(EQUATIONS_MATRICES_FULL_MATRICES)
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
          CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
            CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
          CASE DEFAULT
            localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_DERIVED_TYPE)
      !-----------------------------------------------------------------
      ! D e r i v e d   f i e l d
      !-----------------------------------------------------------------
      ! We want to be able to set which derived variables are calculated before finishing the derived
      ! field, so don't create field variables or check the provided field until the finish action.
      NULLIFY(equationsDerived)
      CALL EquationsSet_DerivedGet(equationsSet,equationsDerived,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsDerived%derivedFieldAutoCreated) THEN
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsDerived%derivedField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsDerived%derivedField,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_LabelSet(equationsDerived%derivedField,"Derived Field",err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsDerived%derivedField,FIELD_DEPENDENT_TYPE,err,error,*999)
          NULLIFY(geometricDecomposition)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsDerived%derivedField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsDerived%derivedField,geometricField,err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        ALLOCATE(variableTypes(equationsDerived%numberOfVariables),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate derived field variable types.",err,error,*999)
        variableIdx=0
        DO derivedIdx=1,EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES
          IF(equationsDerived%variableTypes(derivedIdx)/=0) THEN
            variableIdx=variableIdx+1
            variableTypes(variableIdx)=equationsDerived%variableTypes(derivedIdx)
          END IF
        ENDDO !derivedIdx
        numberOfTensorComponents=NUMBER_OF_VOIGT(numberOfDimensions)
        IF(equationsDerived%derivedFieldAutoCreated) THEN
          CALL Field_NumberOfVariablesSetAndLock(equationsDerived%derivedField,equationsDerived%numberOfVariables,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsDerived%derivedField,variableTypes,err,error,*999)
          DO derivedIdx=1,EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES
            variableType=equationsDerived%variableTypes(derivedIdx)
            IF(variableType/=0) THEN
              CALL Field_DataTypeSetAndLock(equationsDerived%derivedField,variableType,FIELD_DP_TYPE,err,error,*999)
              SELECT CASE(derivedIdx)
              CASE(EQUATIONS_SET_SMALL_STRAIN_TENSOR)
                CALL Field_DimensionSetAndLock(equationsDerived%derivedField,variableType,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_VariableLabelSet(equationsDerived%derivedField,variableType,"Strain",err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(equationsDerived%derivedField,variableType,numberOfTensorComponents, &
                  & err,error,*999)
              CASE(EQUATIONS_SET_CAUCHY_STRESS_TENSOR)
                CALL Field_DimensionSetAndLock(equationsDerived%derivedField,variableType,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_VariableLabelSet(equationsDerived%derivedField,variableType,"Stress",err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(equationsDerived%derivedField,variableType,numberOfTensorComponents, &
                  & err,error,*999)
              CASE(EQUATIONS_SET_ELASTIC_WORK)
                CALL Field_DimensionSetAndLock(equationsDerived%derivedField,variableType,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_VariableLabelSet(equationsDerived%derivedField,variableType,"Work",err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(equationsDerived%derivedField,variableType,1,err,error,*999)
              CASE DEFAULT
                CALL FlagError("The specified derived field type of "//TRIM(NumberToVString(derivedIdx,"*",err,error))// &
                  & " is not supported for a linear elasticity equations set type.",err,error,*999)
              END SELECT
            ENDIF
          ENDDO !derivedIdx
          !Finish creating the derived field
          CALL Field_CreateFinish(equationsDerived%derivedField,err,error,*999)
        ELSE
          !Check the user specified derived field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)            
          DO derivedIdx=1,EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES
            variableType=equationsDerived%variableTypes(derivedIdx)
            IF(variableType/=0) THEN
              CALL Field_DataTypeCheck(equationsSetSetup%field,variableType,FIELD_DP_TYPE,err,error,*999)
              SELECT CASE(derivedIdx)
              CASE(EQUATIONS_SET_SMALL_STRAIN_TENSOR)
                CALL Field_DimensionCheck(equationsDerived%derivedField,variableType,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsDerived%derivedField,variableType,numberOfTensorComponents, &
                  & err,error,*999)
              CASE(EQUATIONS_SET_CAUCHY_STRESS_TENSOR)
                CALL Field_DimensionCheck(equationsDerived%derivedField,variableType,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsDerived%derivedField,variableType,numberOfTensorComponents, &
                  & err,error,*999)
              CASE(EQUATIONS_SET_ELASTIC_WORK)
                CALL Field_DimensionCheck(equationsDerived%derivedField,variableType,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsDerived%derivedField,variableType,1,err,error,*999)
              CASE DEFAULT
                CALL FlagError("The specified derived field type of "//TRIM(NumberToVString(derivedIdx,"*",err,error))// &
                  & " is not supported for a linear elasticity equations set type.",err,error,*999)
              END SELECT
            ENDIF
          ENDDO !derivedIdx
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a linear elasticity equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("LinearElasticity_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("LinearElasticity_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE LinearElasticity_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a linear elasticity equation type of an elasticity equations set class.
  SUBROUTINE LinearElasticity_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError

    ENTERS("LinearElasticity_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE, &
      & EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE, & 
      & EQUATIONS_SET_THREE_DIMENSIONAL_ORTHOTROPIC_SUBTYPE)
      SELECT CASE(solutionMethod)
      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
        equationsSet%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
      CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The specified solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_PLATE_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_SHELL_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a linear elasticity equation type of an elasticity equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("LinearElasticity_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("LinearElasticity_EquationsSetSolutionMethodSet",err,error)
    EXITS("LinearElasticity_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE LinearElasticity_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a linear elasticity equation type of an elasticity equations set class.
  SUBROUTINE LinearElasticity_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("LinearElasticity_EquationsSetSpecificationSet",err,error,*999)

    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    
    SELECT CASE(specification(3))
    CASE(EQUATIONS_SET_THREE_DIMENSIONAL_ISOTROPIC_SUBTYPE, &
      & EQUATIONS_SET_THREE_DIMENSIONAL_ORTHOTROPIC_SUBTYPE, & 
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE, &
      & EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
      !ok
    CASE(EQUATIONS_SET_PLATE_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_SHELL_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
        & " is not valid for a linear elasticity equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    CALL EquationsSet_SpecificationSet(equationsSet,3,[EQUATIONS_SET_ELASTICITY_CLASS, &
      & EQUATIONS_SET_LINEAR_ELASTICITY_TYPE,specification(3)],err,error,*999)
 
    EXITS("LinearElasticity_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("LinearElasticity_EquationsSetSpecificationSet",err,error)
    EXITS("LinearElasticity_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE LinearElasticity_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the linear elasticity problem.
  SUBROUTINE LinearElasticity_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem set to setup a linear elasticity equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("LinearElasticity_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_NO_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a linear elasticity type of an elasticity problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    SELECT CASE(problemSetup%setupType)
    CASE(PROBLEM_SETUP_INITIAL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Do nothing????
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Do nothing????
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a simple control loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the control loops
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)            
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVERS_TYPE)
      !Get the control loop
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoop)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Start the solvers creation
        NULLIFY(solvers)
        CALL Solvers_CreateStart(controlLoop,solvers,err,error,*999)
        CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
        !Set the solver to be a linear solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
        !Set solver defaults
        CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solvers
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      !Get the control loop
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoop)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      !Get the solver
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      NULLIFY(solver)
      CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Create the solver equations
        NULLIFY(solverEquations)
        CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
        CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
        CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solver equations
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        !Finish the solver equations creation
        CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a linear elasticity problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("LinearElasticity_ProblemSetup")
    RETURN
999 ERRORSEXITS("LinearElasticity_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE LinearElasticity_ProblemSetup

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a linear elasticity type problem.
  SUBROUTINE LinearElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specifiation to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("LinearElasticity_ProblemSpecificationSet",err,error,*999)

    IF(SIZE(problemSpecification,1)<3) THEN
      localError="The size of the specified problem specification array of "// &
        & TRIM(NumberToVString(SIZE(problemSpecification,1),"*",err,error))// &
        & " is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    problemSubtype=problemSpecification(3)
    
    SELECT CASE(problemSubtype)
    CASE(PROBLEM_NO_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a linear elasticity problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Set full specification
    CALL Problem_SpecificationSet(problem,3,[PROBLEM_ELASTICITY_CLASS, PROBLEM_LINEAR_ELASTICITY_TYPE, problemSubtype], &
      & err,error,*999)

    EXITS("LinearElasticity_ProblemSpecificationSet")
    RETURN
999 ERRORS("LinearElasticity_ProblemSpecificationSet",err,error)
    EXITS("LinearElasticity_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE LinearElasticity_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculated the strain matrix, B, for a linear elasticity equations set at a Gauss point.
  SUBROUTINE LinearElasticity_StrainMatrixCalculateGauss(numberOfDimensions,numberOfXi,numberOfElementParameters, &
    & interpPointMetrics,gaussPointIdx,quadratureSchemes,BMatrix,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: numberOfXi !<The number of xi dimensions
    INTEGER(INTG), INTENT(IN) :: numberOfElementParameters(:) !<numberOfElementParameters(componentIdx). The number of interpolation element parameters for componentIdx'th component.
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: interpPointMetrics !<A pointer to the interpolated point metrics evaluated at the gauss point.
    INTEGER(INTG), INTENT(IN) :: gaussPointIdx !<The gauss point number to calculated the strain matrix at.
    TYPE(QuadratureSchemePtrType), INTENT(IN) :: quadratureSchemes(:) !<quadratureSchemes(componentIdx)%ptr. A pointer to the quadrature scheme for componentIdx'th component.
    REAL(DP), INTENT(INOUT) :: BMatrix(:,:) !<BMatrix(voigtIdx,elementDOFIdx). On return, the strain matrix at the Gauss point.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_ELEMENT_PARAMETERS=64
    INTEGER(INTG) :: componentIdx1,componentIdx2,elementDOFIdx,elementParameterIdx,otherDirection,otherDirectionIdx, &
      & totalNumberOfElementParameters,xiIdx
    REAL(DP) :: dPhidX(3,3*MAX_NUMBER_OF_ELEMENT_PARAMETERS),dPhidXi(3,3*MAX_NUMBER_OF_ELEMENT_PARAMETERS)

    ENTERS("LinearElasticity_StrainMatrixCalculateGauss",err,error,*999)

#ifdef PRECHECKS
    IF(numberOfDimensions<1.OR.numberOfDimensions>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
        & " is invalid. The number of dimensions must be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(numberOfXi<1.OR.numberOfXi>numberOfDimensions) THEN
      localError="The specified number of xi of "//TRIM(NumberToVString(numberOfXi,"*",err,error))// &
        & " is invalid. The number of xi must be >= 1 and <= the number of dimensions of "// &
        & TRIM(NumberToVString(numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(numberOfDimensions>SIZE(numberOfElementParameters,1)) THEN
      localError="The size of the supplied number of element parameters array of "// &
        & TRIM(NumberToVString(SIZE(numberOfElementParameters,1),"*",err,error))// &
        & " is invalid. The size must be >= to the number of dimensions of "// &
        & TRIM(NumberToVString(numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(interpPointMetrics)) CALL FlagError("Interpolated point metrics is not associated.",err,error,*999)
    IF(numberOfDimensions>SIZE(quadratureSchemes,1)) THEN
      localError="The size of the supplied quadrature schemes array of "// &
        & TRIM(NumberToVString(SIZE(quadratureSchemes,1),"*",err,error))// &
        & " is invalid. The size must be >= to the number of dimensions of "// &
        & TRIM(NumberToVString(numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(NUMBER_OF_VOIGT(numberOfDimensions)>SIZE(BMatrix,1)) THEN
      localError="The size of the first dimension of the specified BMatrix array of "// &
        & TRIM(NumberToVString(SIZE(BMatrix,1),"*",err,error))// &
        & " is invalid. The size must be >= the number of Voigt terms of "// &
        & TRIM(NumberToVString(NUMBER_OF_VOIGT(numberOfDimensions),"*",err,error))//" for "// &
        & TRIM(NumberToVString(numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    totalNumberOfElementParameters=SUM(numberOfElementParameters(1:numberOfDimensions))
    IF(totalNumberOfElementParameters>SIZE(BMatrix,2)) THEN
      localError="The size of the second dimension of the specified BMatrix array of "// &
        & TRIM(NumberToVString(SIZE(BMatrix,2),"*",err,error))// &
        & "is invalid. The size must be >= the total number of element parameters of "// &
        & TRIM(NumberToVString(totalNumberOfElementParameters,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    maxNumberOfElementParameters=MAX(numberOfElementParameters(1:numberOfDimensions))
    IF(maxNumberOfElementParameters>MAX_NUMBER_OF_ELEMENT_PARAMETERS) THEN
      localError="The maximum number of element parameters across dimensions of "// &
        & TRIM(NumberToVString(maxNumberOfElementParameters,"*",err,error))// &
        & " is greater than the MAX_NUMBER_OF_ELEMENT_PARAMETERS of "// &
        & TRIM(NumberToVString(MAX_NUMBER_OF_ELEMENT_PARAMETERS,"*",err,error))// &
        & ". Increase MAX_NUMBER_OF_ELEMENT_PARAMETERS."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    !Loop over the dimensions
    elementDOFIdx=0
    DO componentIdx1=1,numberOfDimensions
      !Loop over element parameters
      DO elementParameterIdx=1,numberOfElementParameters(componentIdx1)
        elementDOFIdx=elementDOFIdx+1
        DO xiIdx=1,numberOfXi
          CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureSchemes(componentIdx1)%ptr,elementParameterIdx, &
            & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,dPhidXi(xiIdx,elementDOFIdx), &
            & err,error,*999)
        ENDDO !xiIdx
        DO componentIdx2=1,numberOfDimensions
          dPhidX(componentIdx2,elementDOFIdx)= &
            & DOT_PRODUCT(dPhidXi(1:numberOfXi,elementDOFIdx),interpPointMetrics%dXidX(1:numberOfXi,componentIdx2))
        ENDDO !componentIdx2
      ENDDO !elementParameterIdx
    ENDDO !componentIdx1
    !Loop over the dimensions
    elementDOFIdx=0
    !Loop over element parameters
    DO componentIdx1=1,numberOfDimensions
      DO elementParameterIdx=1,numberOfElementParameters(componentIdx1)
        elementDOFIdx=elementDOFIdx+1
        BMatrix(1:NUMBER_OF_VOIGT(numberOfDimensions),elementDOFIdx)=0.0_DP
        BMatrix(TENSOR_TO_VOIGT(componentIdx1,componentIdx1,numberOfDimensions),elementDOFIdx)= &
          & dPhidX(componentIdx1,elementDOFIdx)
        DO otherDirectionIdx=1,NUMBER_OTHER_DIRECTIONS(numberOfDimensions)
          otherDirection=OTHER_DIRECTIONS(componentIdx1,otherDirectionIdx,numberOfDimensions)
          BMatrix(TENSOR_TO_VOIGT(componentIdx1,otherDirection,numberOfDimensions),elementDOFIdx)=&
            & dPhidX(otherDirection,elementDOFIdx)
        ENDDO !otherDirectionIdx
      ENDDO !elementParameterIdx
    ENDDO !componentIdx1

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Strain matrix:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Gauss point = ",gaussPointIdx,err,error,*999)
      totalNumberOfElementParameters=SUM(numberOfElementParameters(1:numberOfDimensions))
      CALL WriteStringMatrixTranspose(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_VOIGT(numberOfDimensions),1,1, &
        & totalNumberOfElementParameters,8,8,BMatrix(1:NUMBER_OF_VOIGT(numberOfDimensions), &
        & 1:totalNumberOfElementParameters),WRITE_STRING_MATRIX_NAME_AND_INDICES, &
        & '("  B^t','(:,",I3,")',':",8(X,E13.6))','(13X,8(X,E13.6))',err,error,*999)
    ENDIF
   
    EXITS("LinearElasticity_StrainMatrixCalculateGauss")
    RETURN
999 ERRORS("LinearElasticity_StrainMatrixCalculateGauss",err,error)
    EXITS("LinearElasticity_StrainMatrixCalculateGauss")
    RETURN 1
    
  END SUBROUTINE LinearElasticity_StrainMatrixCalculateGauss

  !
  !================================================================================================================================
  !

  !>Calculates the strain field for a linear elasticity linear element equations set.
  SUBROUTINE LinearElasticity_StressStrainCalculate(equationsSet,derivedType,fieldVariable,err,error,*)

    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate stress/strain for.
    INTEGER(INTG), INTENT(IN) :: derivedType !<The type of derived field to calculate.     
    TYPE(FieldVariableType), POINTER, INTENT(INOUT) :: fieldVariable !<The field variable to store the stress/strain in.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dataPointIdx,dataPointNumber,dependentVariableType,elementIdx,elementNumber, &
      & esSpecification(3),fieldInterpolation,fieldVarType,finishIdx,gaussPointIdx,numberOfComponents, &
      & numberOfDataPoints,numberOfDependentComponents,numberOfDimensions,numberOfElementXi,numberOfGauss,numberOfTimes, &
      & numberOfXi,outputType,partIdx,startIdx
    REAL(DP) :: xi(3),values(3,3)
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime4(1),userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1)
    TYPE(BasisType), POINTER :: basis
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(DataProjectionType), POINTER :: dataProjection
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: interpolation
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: field,dependentField,geometricField,materialsField
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
    TYPE(FieldInterpolationParametersType), POINTER :: geometricInterpParameters,dependentInterpParameters, &
      & materialsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint,dependentInterpPoint,materialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("LinearElasticity_StressStrainCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    !Get the coordinate system
    NULLIFY(coordinateSystem)
    CALL EquationsSet_CoordinateSystemGet(equationsSet,coordinateSystem,err,error,*999)
    CALL CoordinateSystem_DimensionGet(coordinateSystem,numberOfDimensions,err,error,*999)
    !Check the provided strain field variable has appropriate components and interpolation
    IF(derivedType == EQUATIONS_SET_ELASTIC_WORK) THEN
      numberOfComponents=1
    ELSE
      numberOfComponents=NUMBER_OF_VOIGT(numberOfDimensions)
    ENDIF
    NULLIFY(field)
    CALL FieldVariable_FieldGet(fieldVariable,field,err,error,*999)
    CALL FieldVariable_VariableTypeGet(fieldVariable,fieldVarType,err,error,*999)
    CALL FieldVariable_NumberOfComponentsCheck(fieldVariable,numberOfComponents,err,error,*999)
    CALL FieldVariable_ComponentInterpolationGet(fieldVariable,1,fieldInterpolation,err,error,*999)
    !Check the interpolation type
    SELECT CASE(fieldInterpolation)
    CASE(FIELD_CONSTANT_INTERPOLATION)
      CALL FlagError("Can not calculate stress or strain for a field with constant interpolation.",err,error,*999)
    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
      !OK
    CASE(FIELD_NODE_BASED_INTERPOLATION)
      CALL FlagError("Stress/strain calculation is not implemented for a field with node based interpolation.",err,error,*999)
    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
      CALL FlagError("Stress/strain calculation is not implemented for a field with grid point based interpolation.", &
        & err,error,*999)
    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
      !OK
    CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
      !OK
    CASE DEFAULT
      localError="The field interpolation type for component 1 of "//TRIM(NumberToVString(fieldInterpolation,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Check that all the components have the same interpolation type
    DO componentIdx=2,numberOfComponents
      CALL FieldVariable_ComponentInterpolationCheck(fieldVariable,componentIdx,fieldInterpolation,err,error,*999)
    ENDDO !componentIdx

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(linearMapping)
    CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
    NULLIFY(dependentVariable)
    CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,dependentVariable,err,error,*999)
    CALL FieldVariable_VariableTypeGet(dependentVariable,dependentVariableType,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)      
   
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*999)

    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    
    !Grab interpolation points
    NULLIFY(interpolation)
    CALL Equations_InterpolationGet(equations,interpolation,err,error,*999)    
    NULLIFY(geometricInterpParameters)
    CALL EquationsInterpolation_GeometricParametersGet(interpolation,FIELD_U_VARIABLE_TYPE,geometricInterpParameters, &
      & err,error,*999)
    NULLIFY(geometricInterpPoint)
    CALL EquationsInterpolation_GeometricPointGet(interpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint,err,error,*999)
    NULLIFY(geometricInterpPointMetrics)
    CALL EquationsInterpolation_GeometricPointMetricsGet(interpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPointMetrics, &
      & err,error,*999)
    NULLIFY(dependentInterpParameters)
    CALL EquationsInterpolation_DependentParametersGet(interpolation,dependentVariableType,dependentInterpParameters, &
      & err,error,*999)
    NULLIFY(dependentInterpPoint)
    CALL EquationsInterpolation_DependentPointGet(interpolation,dependentVariableType,dependentInterpPoint,err,error,*999)
    NULLIFY(materialsInterpParameters)
    NULLIFY(materialsInterpPoint)
    IF(ASSOCIATED(materialsField)) THEN
      CALL EquationsInterpolation_MaterialsParametersGet(interpolation,FIELD_U_VARIABLE_TYPE,materialsInterpParameters, &
        & err,error,*999)
      CALL EquationsInterpolation_MaterialsPointGet(interpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint,err,error,*999)
    ENDIF
 
    numberOfTimes=0
    elementUserElapsed=0.0_SP
    elementSystemElapsed=0.0_SP
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)

    !Loop over the two parts: 1 - boundary and ghost elements, 2 - internal
    DO partIdx=1,2          
      
      IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
      ENDIF
      
      IF(partIdx==1) THEN
        CALL DomainMapping_BoundaryStartGet(elementsMapping,startIdx,err,error,*999)
        CALL DomainMapping_GhostFinishGet(elementsMapping,finishIdx,err,error,*999)
      ELSE
        CALL DomainMapping_InternalStartGet(elementsMapping,startIdx,err,error,*999)
        CALL DomainMapping_InternalFinishGet(elementsMapping,finishIdx,err,error,*999)
      ENDIF
      
      !Loop over (1) the boundary and ghost elements, (2) the internal elements
      DO elementIdx=startIdx,finishIdx
        
        numberOfTimes=numberOfTimes+1
        CALL DomainMapping_NumberGet(elementsMapping,elementIdx,elementNumber,err,error,*999)
        
        IF(diagnostics1) THEN
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",elementNumber,err,error,*999)
        ENDIF

        NULLIFY(basis)
        CALL DomainElements_ElementBasisGet(domainElements,elementNumber,basis,err,error,*999)
        CALL Basis_NumberOfXiGet(basis,numberOfXi,err,error,*999)
      
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
 
        SELECT CASE(fieldInterpolation)
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
          !Interpolate dependent, geometric, fibre etc. fields
          xi=[0.5_DP,0.5_DP,0.5_DP]
          CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),geometricInterpPoint,err,error,*999)
          CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
          CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),dependentInterpPoint,err,error,*999)
          CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),materialsInterpPoint,err,error,*999)
       
          CALL LinearElasticity_StressStrainPoint(equationsSet,derivedType,numberOfDimensions,numberOfXi,gaussPointIdx, &
            & elementNumber,geometricInterpPoint,geometricInterpPointMetrics,dependentInterpPoint, &
            & materialsInterpPoint,values,err,error,*999)

          IF(derivedType == EQUATIONS_SET_ELASTIC_WORK) THEN
            CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
              & 1,values(1,1),err,error,*999)
          ELSE
            !We only want to store the independent components 
            SELECT CASE(numberOfDimensions)
            CASE(1)
              ! 1 dimensional problem
              CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                & TENSOR_TO_VOIGT1(1,1),values(1,1),err,error,*999)
            CASE(2)
              ! 2 dimensional problem
              CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                & TENSOR_TO_VOIGT2(1,1),values(1,1),err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                & TENSOR_TO_VOIGT2(1,2),values(1,2),err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                & TENSOR_TO_VOIGT2(2,2),values(2,2),err,error,*999)
            CASE(3)
              ! 3 dimensional problem
              CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                & TENSOR_TO_VOIGT3(1,1),values(1,1),err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                & TENSOR_TO_VOIGT3(1,2),values(1,2),err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                & TENSOR_TO_VOIGT3(1,3),values(1,3),err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                & TENSOR_TO_VOIGT3(2,2),values(2,2),err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                & TENSOR_TO_VOIGT3(2,3),values(2,3),err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                & TENSOR_TO_VOIGT3(3,3),values(3,3),err,error,*999)
            CASE DEFAULT
              localError="The number of dimensions of "//TRIM(NumberToVString(numberofDimensions,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
          
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)            
            
          NULLIFY(quadratureScheme)               
          CALL Basis_QuadratureSchemeGet(basis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
          CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme,numberOfGauss,err,error,*999)
          
          !Loop over gauss points        
          DO gaussPointIdx=1,numberOfGauss
            
            IF(diagnostics1) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Gauss point number = ",gaussPointIdx,err,error,*999)
            ENDIF
            
            !Interpolate dependent, geometric, fibre etc. fields
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
              & err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
              & err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
              & err,error,*999)
            
            CALL LinearElasticity_StressStrainPoint(equationsSet,derivedType,numberOfDimensions,numberOfXi,gaussPointIdx, &
              & elementNumber,geometricInterpPoint,geometricInterpPointMetrics,dependentInterpPoint, &
              & materialsInterpPoint,values,err,error,*999)
            
            IF(derivedType == EQUATIONS_SET_ELASTIC_WORK) THEN
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,1,values(1,1),err,error,*999)
            ELSE
              !We only want to store the independent components 
              SELECT CASE(numberOfDimensions)
              CASE(1)
                ! 1 dimensional problem
                CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                  & elementNumber,TENSOR_TO_VOIGT1(1,1),values(1,1),err,error,*999)
              CASE(2)
                ! 2 dimensional problem
                CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                  & elementNumber,TENSOR_TO_VOIGT2(1,1),values(1,1),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                  & elementNumber,TENSOR_TO_VOIGT2(1,2),values(1,2),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                  & elementNumber,TENSOR_TO_VOIGT2(2,2),values(2,2),err,error,*999)
              CASE(3)
                ! 3 dimensional problem
                CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                  & elementNumber,TENSOR_TO_VOIGT3(1,1),values(1,1),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                  & elementNumber,TENSOR_TO_VOIGT3(1,2),values(1,2),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                  & elementNumber,TENSOR_TO_VOIGT3(1,3),values(1,3),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                  & elementNumber,TENSOR_TO_VOIGT3(2,2),values(2,2),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                  & elementNumber,TENSOR_TO_VOIGT3(2,3),values(2,3),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                  & elementNumber,TENSOR_TO_VOIGT3(3,3),values(3,3),err,error,*999)
              CASE DEFAULT
                localError="The number of dimensions of "//TRIM(NumberToVString(numberofDimensions,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          ENDDO !gaussPointIdx
          
        CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
          
          NULLIFY(dataProjection)
          CALL Field_DataProjectionGet(field,dataProjection,err,error,*999)
          NULLIFY(dataPoints)
          CALL DecompositionTopology_DecompositionDataPointsGet(decompositionTopology,dataPoints,err,error,*999)
          CALL DecompositionDataPoints_ElementNumberOfDataPointsGet(dataPoints,elementNumber,numberOfDataPoints,err,error,*999)
 
          DO dataPointIdx=1,numberOfDataPoints

            CALL DecompositionDataPoints_ElementDataGlobalNumberGet(dataPoints,dataPointIdx,elementNumber,dataPointNumber, &
              & err,error,*999)
            CALL DataProjection_ResultElementXiGet(dataProjection,dataPointNumber,numberOfElementXi,xi,err,error,*999)
            CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),geometricInterpPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
            CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),dependentInterpPoint,err,error,*999)
            CALL Field_InterpolateXi(FIRST_PART_DERIV,xi(1:numberOfXi),materialsInterpPoint,err,error,*999)
            
            CALL LinearElasticity_StressStrainPoint(equationsSet,derivedType,numberOfDimensions,numberOfXi,gaussPointIdx, &
              & elementNumber,geometricInterpPoint,geometricInterpPointMetrics,dependentInterpPoint, &
              & materialsInterpPoint,values,err,error,*999)
            
            IF(derivedType == EQUATIONS_SET_ELASTIC_WORK) THEN
              CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                & 1,values(1,1),err,error,*999)
            ELSE
              !We only want to store the independent components 
              SELECT CASE(numberOfDimensions)
              CASE(1)
                ! 1 dimensional problem
                CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                  & TENSOR_TO_VOIGT1(1,1),values(1,1),err,error,*999)
              CASE(2)
                ! 2 dimensional problem
                CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                  & TENSOR_TO_VOIGT2(1,1),values(1,1),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                  & TENSOR_TO_VOIGT2(1,2),values(1,2),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                  & TENSOR_TO_VOIGT2(2,2),values(2,2),err,error,*999)
              CASE(3)
                ! 3 dimensional problem
                CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                  & TENSOR_TO_VOIGT3(1,1),values(1,1),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                  & TENSOR_TO_VOIGT3(1,2),values(1,2),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                  & TENSOR_TO_VOIGT3(1,3),values(1,3),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                  & TENSOR_TO_VOIGT3(2,2),values(2,2),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                  & TENSOR_TO_VOIGT3(2,3),values(2,3),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalElement(fieldVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
                  & TENSOR_TO_VOIGT3(3,3),values(3,3),err,error,*999)
              CASE DEFAULT
                localError="The number of dimensions of "//TRIM(NumberToVString(numberofDimensions,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
                      
          ENDDO !dataPointIdx
          
        CASE DEFAULT
          localError="The field interpolation type for component 1 of "// &
            & TRIM(NumberToVString(fieldInterpolation,"*",err,error))//" is invalid or not implemented."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        
      ENDDO !elementIdx
      
      !Output timing information if required
      IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
        userElapsed=userTime2(1)-userTime1(1)
        systemElapsed=systemTime2(1)-systemTime1(1)
        elementUserElapsed=elementUserElapsed+userElapsed
        elementSystemElapsed=elementSystemElapsed+systemElapsed
        IF(partIdx==1) THEN
          CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
          CALL Profiling_TimingsOutput(1,"Boundary+ghost elements calculation",userElapsed,systemElapsed,err,error,*999)
        ELSE
          CALL Profiling_TimingsOutput(1,"Internal elements calculation",userElapsed,systemElapsed,err,error,*999)
          IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element calculation", &
            & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
        ENDIF
      ENDIF !outputType>=EQUATIONS_TIMING_OUTPUT
      
      IF(partIdx==1) THEN
        IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
          CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
          CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
        ENDIF
        !Start to update the field
        CALL FieldVariable_ParameterSetUpdateStart(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      ELSE
        !Finish to update the field
        CALL FieldVariable_ParameterSetUpdateFinish(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        !Output timing information if required
        IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
          CALL CPUTimer(USER_CPU,userTime4,err,error,*999)
          CALL CPUTimer(SYSTEM_CPU,systemTime4,err,error,*999)
          userElapsed=userTime4(1)-userTime3(1)
          systemElapsed=systemTime4(1)-systemTime3(1)
          CALL Profiling_TimingsOutput(1,"Parameters update transfer",userElapsed,systemElapsed,err,error,*999)
        ENDIF !equations%outputType>=EQUATIONS_TIMING_OUTPUT
      ENDIF
      
    ENDDO !partIdx
        
    EXITS("LinearElasticity_StressStrainCalculate")
    RETURN
999 ERRORSEXITS("LinearElasticity_StressStrainCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE LinearElasticity_StressStrainCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates stress and strain at a point. \TODO merge this with interpolate xi below.
  SUBROUTINE LinearElasticity_StressStrainPoint(equationsSet,evaluateType,numberOfDimensions,numberOfXi,pointNumber, &
    & elementNumber,geometricInterpPoint,geometricInterpPointMetrics,dependentInterpPoint,materialsInterpPoint,values, &
    & err,error,*)
    ! Argument variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate the tensor for
    INTEGER(INTG), INTENT(IN) :: evaluateType !<The type of tensor to evaluate.
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: numberOfXi !<The number of xi directions
    INTEGER(INTG), INTENT(IN) :: pointNumber !<The point number to evaluate the tensor for
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The user element number to evaluate the tensor for
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint !<The geometric interpolated point
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics !<The geometric interpolated point metrics
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint !<The dependent interpolated point
    TYPE(FieldInterpolatedPointType), POINTER :: materialsInterpPoint !<The materials interpolated point
    REAL(DP), INTENT(OUT) :: values(:,:) !<On exit, the interpolated tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: coordinateIdx,dimensionIdx,esSpecification(3),partialDerivativeIndex,voigtIdx1,voigtIdx2,xiIdx
    REAL(DP) :: C(6,6),dudx(3,3),eV(6),sigmaV(6),W

    ENTERS("LinearElasticity_StressStrainPoint",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    !Calculate du/dx
    DO coordinateIdx=1,numberOfDimensions
      DO dimensionIdx=1,numberOfDimensions
        dudx(coordinateIdx,dimensionIdx)=0.0_DP
        DO xiIdx=1,numberOfXi
          partialDerivativeIndex=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx)
          dudx(coordinateIdx,dimensionIdx)=dudx(coordinateIdx,dimensionIdx)+ &
            & dependentInterpPoint%values(coordinateIdx,partialDerivativeIndex)* &
            & geometricInterpPointMetrics%dXidX(xiIdx,dimensionIdx)
        ENDDO !xiIdx
      ENDDO !dimensionIdx
    ENDDO !coordinateIdx
    
    !Calculate eV, the small strain tensor in Voigt form
    eV=0.0_DP
    DO coordinateIdx=1,numberOfDimensions
      !Calculate direct strains, e_ii = du_i/dx^i
      eV(TENSOR_TO_VOIGT(coordinateIdx,coordinateIdx,numberOfDimensions))=dudx(coordinateIdx,coordinateIdx)
      !Calculate shear strains, e_ij = du_j/dx^i + du_i/dx^j, i /= j
      DO dimensionIdx=1,numberOfDimensions
        IF(dimensionIdx/=coordinateIdx) THEN
          eV(TENSOR_TO_VOIGT(coordinateIdx,dimensionIdx,numberOfDimensions))= &
            & eV(TENSOR_TO_VOIGT(coordinateIdx,dimensionIdx,numberOfDimensions))+ &
            & dudx(coordinateIdx,dimensionIdx)
        ENDIF
      ENDDO !dimensionIdx
    ENDDO !coordinateIdx
    
    IF(evaluateType /= EQUATIONS_SET_SMALL_STRAIN_TENSOR) THEN
      !Calculate c
      CALL LinearElasticity_ElasticityTensor(esSpecification(3),numberOfDimensions,materialsInterpPoint,C,err,error,*999)
      !Calculate sigma from sigma=c::e
      sigmaV=0.0_DP
      DO voigtIdx1=1,NUMBER_OF_VOIGT(numberOfDimensions)
        sigmaV(voigtIdx1)=0.0_DP
        DO voigtIdx2=1,NUMBER_OF_VOIGT(numberOfDimensions)
          sigmaV(voigtIdx1)=sigmaV(voigtIdx1)+C(voigtIdx1,voigtIdx2)*eV(voigtIdx2)
        ENDDO !voightIdx2
      ENDDO !voightIdx1
      IF(evaluateType == EQUATIONS_SET_ELASTIC_WORK) THEN
        !Calculate W=sigma::e
        W=DOT_PRODUCT(sigmaV(1:NUMBER_OF_VOIGT(numberOfDimensions)),eV(1:NUMBER_OF_VOIGT(numberOfDimensions)))
      ENDIF
    ENDIF
    
    SELECT CASE(evaluateType)
    CASE(EQUATIONS_SET_SMALL_STRAIN_TENSOR)
      CALL VoigtToTensor(numberOfDimensions,[TENSOR_COVARIANT_INDEX],eV,values,err,error,*999)
    CASE(EQUATIONS_SET_CAUCHY_STRESS_TENSOR)
      CALL VoigtToTensor(numberOfDimensions,[TENSOR_CONTRAVARIANT_INDEX],sigmaV,values,err,error,*999)
    CASE(EQUATIONS_SET_ELASTIC_WORK)
      values(1,1)=W
    CASE DEFAULT
      CALL FlagError("The tensor evalaute type of "//TRIM(NumberToVString(evaluateType,"*",err,error))//" is invalid "// &
        & "for linear elasticity equation sets.",err,error,*999)
    END SELECT
    
    EXITS("LinearElasticity_StressStrainPoint")
    RETURN
999 ERRORS("LinearElasticity_StressStrainPoint",err,error)
    EXITS("LinearElasticity_StressStrainPoint")
    RETURN 1
    
  END SUBROUTINE LinearElasticity_StressStrainPoint


END MODULE LinearElasticityRoutines
