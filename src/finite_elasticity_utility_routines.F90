!> \file
!> \author Chris Bradley
!> \brief This module handles utility routines for finite elasticity to allow access from other modules.
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

!>This module handles utility routines for finite elasticity to allow access from other modules.
MODULE FiniteElasticityUtilityRoutines

  USE BaseRoutines
  USE Constants
  USE CoordinateSystemRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE Strings
  USE Types
  
#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FiniteElasticity_DeformationGradientTensorCalculate

  PUBLIC FiniteElasticity_PolarDecompositionCalculate

  PUBLIC FiniteElasticity_StructuralTensorComponentCalculate

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the deformation gradient tensor at a given interpolated point
  SUBROUTINE FiniteElasticity_DeformationGradientTensorCalculate(dependentInterpPointMetrics,geometricInterpPointMetrics,&
    & dXdNu,F,J,FNu,JNu,err,error,*)

    !Argument variables
    TYPE(FieldInterpolatedPointMetricsType), POINTER, INTENT(IN) :: dependentInterpPointMetrics !<The interpolated point metrics of the deformed/spatial geometry
    TYPE(FieldInterpolatedPointMetricsType), POINTER, INTENT(IN) :: geometricInterpPointMetrics !<The interpolated point metrics of the undeformed/reference geometry
    REAL(DP), INTENT(IN) :: dXdNu(:,:) !<dXdNu(XIdx,NuIdx). The transformation matrix between X and Nu coordinates (identity if there are no fibres). 
    REAL(DP), INTENT(OUT) :: F(:,:) !<F(zIdx,XIdx). On return, the deformation gradient tensor F with respect to X (material geometric) coordinates.
    REAL(DP), INTENT(OUT) :: J !<On return, the Jacobian of the deformation i.e., determinant F
    REAL(DP), INTENT(OUT) :: FNu(:,:) !<F(zIdx,NuIdx). On return, the deformation gradient tensor F with respect to Nu (material fibre) coordinates
    REAL(DP), INTENT(OUT) :: JNu !<On return, the Jacobian of the deformation with respect to material fibre coordinates i.e., determinant FNu
    INTEGER(INTG), INTENT(OUT) :: err   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfXDimensions,numberOfXiDimensions,numberOfZDimensions
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticity_DeformationGradientTensorCalculate",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dependentInterpPointMetrics)) &
      & CALL FlagError("Dependent interpolated point metrics is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(geometricInterpPointMetrics)) &
      & CALL FlagError("Geometric interpolated point metrics is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(dependentInterpPointMetrics%dXdXi)) &
      & CALL FlagError("Dependent interpolated point metrics dXdXi is not allocated.",err,error,*999)
    IF(.NOT.ALLOCATED(geometricInterpPointMetrics%dXdXi)) &
      & CALL FlagError("Geometric interpolated point metrics dXdXi is not allocated.",err,error,*999)
    IF(.NOT.ALLOCATED(geometricInterpPointMetrics%dXidX)) &
      & CALL FlagError("Geometric interpolated point metrics dXidX is not allocated.",err,error,*999)
    IF(SIZE(dXdNu,1)<geometricInterpPointMetrics%numberOfXDimensions) THEN
      localError="The size of the first index of dXdNu of "//TRIM(NumberToVString(SIZE(dXdNu,1),"*",err,error))// &
        & " is too small. The size should be >= "// &
        & TRIM(NumberToVString(geometricInterpPointMetrics%numberOfXDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dXdNu,2)<geometricInterpPointMetrics%numberOfXDimensions) THEN
      localError="The size of the second index of dXdNu of "//TRIM(NumberToVString(SIZE(dXdNu,2),"*",err,error))// &
        & " is too small. The size should be >= "// &
        & TRIM(NumberToVString(geometricInterpPointMetrics%numberOfXDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(F,1)<dependentInterpPointMetrics%numberOfXDimensions) THEN
      localError="The size of the first index of F of "//TRIM(NumberToVString(SIZE(F,1),"*",err,error))// &
        & " is too small. The size should be >= "// &
        & TRIM(NumberToVString(dependentInterpPointMetrics%numberOfXDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(F,2)<geometricInterpPointMetrics%numberOfXDimensions) THEN
      localError="The size of the second index of F of "//TRIM(NumberToVString(SIZE(F,2),"*",err,error))// &
        & " is too small. The size should be >= "// &
        & TRIM(NumberToVString(geometricInterpPointMetrics%numberOfXDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(FNu,1)<dependentInterpPointMetrics%numberOfXDimensions) THEN
      localError="The size of the first index of FNu of "//TRIM(NumberToVString(SIZE(FNu,1),"*",err,error))// &
        & " is too small. The size should be >= "// &
        & TRIM(NumberToVString(dependentInterpPointMetrics%numberOfXDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(FNu,2)<geometricInterpPointMetrics%numberOfXDimensions) THEN
      localError="The size of the second index of F of "//TRIM(NumberToVString(SIZE(FNu,2),"*",err,error))// &
        & " is too small. The size should be >= "// &
        & TRIM(NumberToVString(geometricInterpPointMetrics%numberOfXDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    CALL FieldInterpolatedPointMetrics_NumberOfXDimensionsGet(geometricInterpPointMetrics,numberOfXDimensions,err,error,*999)
    CALL FieldInterpolatedPointMetrics_NumberOfXiDimensionsGet(geometricInterpPointMetrics,numberOfXiDimensions,err,error,*999)
    CALL FieldInterpolatedPointMetrics_NumberOfXDimensionsGet(dependentInterpPointMetrics,numberOfZDimensions,err,error,*999)

    CALL IdentityMatrix(F,err,error,*999)
    CALL IdentityMatrix(FNu,err,error,*999)
    
    !F = dz/dX = dz/dXi * dXi/dX (deformation gradient tensor wrt material coordinates, F)
    CALL MatrixProduct(dependentInterpPointMetrics%dXdXi(1:numberOfZDimensions,1:numberOfXiDimensions), &
      & geometricInterpPointMetrics%dXidX(1:numberOfXiDimensions,1:numberOfXDimensions), &
      & F(1:numberOfZDimensions,1:numberOfXDimensions),err,error,*999)
    CALL Determinant(F(1:numberOfZDimensions,1:numberOfXDimensions),J,err,error,*999)
    !FNu = dz/dNu = F * dX/dNu * = dz/dX * dX/dNu  (deformation gradient tensor wrt material fibre coordinates, FNu)
    CALL MatrixProduct(F(1:numberOfZDimensions,1:numberOfXDimensions),dXdNu(1:numberOfXDimensions,1:numberOfXDimensions), &
       & FNu(1:numberOfZDimensions,1:numberOfXDimensions),err,error,*999)    
    CALL Determinant(FNu(1:numberOfZDimensions,1:numberOfXDimensions),JNu,err,error,*999)

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Calculated deformation gradient tensor:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of z dimensions  = ",numberOfZDimensions,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of X dimensions  = ",numberOfXDimensions,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi dimensions = ",numberOfXiDimensions,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Deformation gradient tensor wrt X coordinates:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfZDimensions,1,1,numberOfXiDimensions, &
        & numberOfXDimensions,numberOfXDimensions,F,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
        & '("    F','(",I1,",:)','   :",3(X,E13.6))','(19X,3(X,E13.6))',err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Determinant F = ",J,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Deformation gradient tensor wrt Nu coordinates:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfZDimensions,1,1,numberOfXDimensions, &
        & numberOfXDimensions,numberOfXDimensions,F,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
        & '("    FNu','(",I1,",:)',' :",3(X,E13.6))','(19X,3(X,E13.6))',err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Determinant FNu = ",JNu,err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_DeformationGradientTensorCalculate")
    RETURN
999 ERRORS("FiniteElasticity_DeformationGradientTensorCalculate",err,error)
    EXITS("FiniteElasticity_DeformationGradientTensorCalculate")
    RETURN 1

  END SUBROUTINE FiniteElasticity_DeformationGradientTensorCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the polar decomposition of a deformation gradient tensor at a given interpolated point. 
  SUBROUTINE FiniteElasticity_PolarDecompositionCalculate(F,R,U,v,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: F(:,:) !<F(xIdx,XIdx). The deformation gradient tensor F to calculate the polar decomposition for
    REAL(DP), INTENT(OUT) :: R(:,:) !<R(xIdx,XIdx). On return, the rotation tensor R of the polar decomposition for F = R.U = v.R
    REAL(DP), INTENT(OUT) :: U(:,:) !<U(XIdx,XIdx). On return, the right stretch tensor U of the polar decomposition for F = R.U = v.R
    REAL(DP), INTENT(OUT) :: v(:,:) !<v(xIdx,xIdx). On return, the left stretch tensor v of the polar decomposition for F = R.U = v.R
    INTEGER(INTG), INTENT(OUT) :: err   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dimensionIdx,numberOfXDimensions,numberOfXiDimensions,numberOfZDimensions
    REAL(DP) :: C(SIZE(F,2),SIZE(F,2)),eValues(SIZE(F,2))
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticity_PolarDecompositionCalculate",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(SIZE(F,1)/=SIZE(F,2)) THEN
      localError="The first dimension of the specified F tensor of "//TRIM(NumberTOVString(SIZE(F,1),"*",err,error))// &
        & " does not match the second dimension of the specified F tensor of "//TRIM(NumberTOVString(SIZE(F,2),"*",err,error))// &
        & ". The F tensor must be square."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(SIZE(R,1)/=SIZE(F,1)) THEN
      localError="The first dimension of the specified R of "//TRIM(NumberToVString(SIZE(R,1),"*",err,error))// &
        & " does not match the first dimension of the specified F of "//TRIM(NumberToVString(SIZE(F,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(R,2)/=SIZE(F,2)) THEN
      localError="The second dimension of the specified R of "//TRIM(NumberToVString(SIZE(R,2),"*",err,error))// &
        & " does not match the second dimension of the specified F of "//TRIM(NumberToVString(SIZE(F,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(U,1)/=SIZE(F,2)) THEN
      localError="The first dimension of the specified U of "//TRIM(NumberToVString(SIZE(U,1),"*",err,error))// &
        & " does not match the second dimension of the specified F of "//TRIM(NumberToVString(SIZE(F,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(U,2)/=SIZE(F,2)) THEN
      localError="The second dimension of the specified U of "//TRIM(NumberToVString(SIZE(U,2),"*",err,error))// &
        & " does not match the second dimension of the specified F of "//TRIM(NumberToVString(SIZE(F,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(v,1)/=SIZE(F,1)) THEN
      localError="The first dimension of the specified v of "//TRIM(NumberToVString(SIZE(v,1),"*",err,error))// &
        & " does not match the first dimension of the specified F of "//TRIM(NumberToVString(SIZE(F,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(v,2)/=SIZE(F,1)) THEN
      localError="The second dimension of the specified v of "//TRIM(NumberToVString(SIZE(v,2),"*",err,error))// &
        & " does not match the first dimension of the specified F of "//TRIM(NumberToVString(SIZE(F,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    numberOfZDimensions = SIZE(F,1)
    numberOfXDimensions = SIZE(F,2)

    U = 0.0_DP
    
    SELECT CASE(numberOfXDimensions)
    CASE(1)
    CASE(2)
    CASE(3)

      !Compute C=F^T.F&
             
      CALL MatrixTransposeProduct(F(1:numberOfZDimensions,1:numberOfXDimensions),F(1:numberOfZDimensions,1:numberOfXDimensions), &
        & C,err,error,*999)
      !Compute eigenvalues of C
      CALL Eigenvalue(C(1:numberOfXDimensions,1:numberOfXDimensions),eValues(1:numberOfXDimensions),err,error,*999)
      !Compute U
      DO dimensionIdx=1,numberOfXDimensions
        U(dimensionIdx,dimensionIdx) = SQRT(eValues(dimensionIdx))
      ENDDO !dimensionIdx

    CASE DEFAULT
      localError="The number of dimensions of "//TRIM(NumberToVString(numberOfXDimensions,"*",err,error))// &
        & " is invalid. The number of dimensions needs to be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    IF(diagnostics1) THEN
    ENDIF

    EXITS("FiniteElasticity_PolarDecompositionCalculate")
    RETURN
999 ERRORS("FiniteElasticity_PolarDecompositionCalculate",err,error)
    EXITS("FiniteElasticity_PolarDecompositionCalculate")
    RETURN 1

  END SUBROUTINE FiniteElasticity_PolarDecompositionCalculate

  !
  !================================================================================================================================
  !

  !>Calculates a component of a structural tensor from the tensor product of two vectors
  SUBROUTINE FiniteElasticity_StructuralTensorComponentCalculate(vector1,vector2,structuralTensor,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: vector1(:) !<vector1(componentIdx). The first vector to calculate the structural tensor for
    REAL(DP), INTENT(IN) :: vector2(:) !<vector2(componentIdx). The second vector to calculate the structural tensor for
    REAL(DP), INTENT(OUT) :: structuralTensor(:,:) !<structuralTensor(componentIdx1,componentIdx2). On return, the calculated structural tensor. 
    INTEGER(INTG), INTENT(OUT) :: err   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfDimensions
    REAL(DP) :: norm1,norm2,norm
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticity_StructuralTensorComponentCalculate",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(SIZE(vector1,1)<1.OR.SIZE(vector1,1)>3) THEN
      localError="The dimension of vector 1 of "//TRIM(NumberToVString(SIZE(vector1,1),"*",err,error))// &
        & " is invalid. The dimension should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(vector1,1)/=SIZE(vector2,1)) THEN
      localError="The size of vector 1 of "//TRIM(NumberToVString(SIZE(vector1,1),"*",err,error))// &
        & " does not match the size of vector 2 of "//TRIM(NumberToVString(SIZE(vector2,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(structuralTensor,1)/=SIZE(vector1,1)) THEN
      localError="The size of the first dimension of the specified structural tensor of "// &
        & TRIM(NumberToVString(SIZE(structuralTensor,1),"*",err,error))// &
        & " does not match the size of the specified vector 1 of "// &
        & TRIM(NumberToVString(SIZE(vector1,1),"*",err,error))//". The sizes should be equal."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(structuralTensor,2)/=SIZE(vector2,1)) THEN
      localError="The size of the second dimension of the specified structural tensor of "// &
        & TRIM(NumberToVString(SIZE(structuralTensor,2),"*",err,error))// &
        & " does not match the size of the specified vector 2 of "// &
        & TRIM(NumberToVString(SIZE(vector2,1),"*",err,error))//". The sizes should be equal."
      CALL FlagError(localError,err,error,*999)
    ENDIF 
#endif

    numberOfDimensions=SIZE(vector1,1)
    SELECT CASE(numberOfDimensions)
    CASE(1)
      structuralTensor(1,1)=1.0_DP
    CASE(2)
      norm1=SQRT(vector1(1)**2+vector1(2)**2)
      norm2=SQRT(vector2(1)**2+vector2(2)**2)
      IF(ABS(norm1)<ZERO_TOLERANCE) THEN
        localError="The norm of vector 1 of "//TRIM(NumberToVString(norm1,"*",err,error))//" is invalid. The norm must be >= 0.0"
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(ABS(norm2)<ZERO_TOLERANCE) THEN
        localError="The norm of vector 2 of "//TRIM(NumberToVString(norm2,"*",err,error))//" is invalid. The norm must be >= 0.0"
        CALL FlagError(localError,err,error,*999)
      ENDIF
      norm=norm1*norm2
      structuralTensor(1,1)=vector1(1)*vector2(1)/norm
      structuralTensor(2,1)=vector1(2)*vector2(1)/norm
      structuralTensor(1,2)=vector1(1)*vector2(2)/norm
      structuralTensor(2,2)=vector1(2)*vector2(2)/norm
    CASE(3)
      norm1=SQRT(vector1(1)**2+vector1(2)**2+vector1(3)**2)
      norm2=SQRT(vector2(1)**2+vector2(2)**2+vector2(3)**2)
      IF(ABS(norm1)<ZERO_TOLERANCE) THEN
        localError="The norm of vector 1 of "//TRIM(NumberToVString(norm1,"*",err,error))//" is invalid. The norm must be >= 0.0"
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(ABS(norm2)<ZERO_TOLERANCE) THEN
        localError="The norm of vector 2 of "//TRIM(NumberToVString(norm2,"*",err,error))//" is invalid. The norm must be >= 0.0"
        CALL FlagError(localError,err,error,*999)
      ENDIF
      norm=norm1*norm2
      structuralTensor(1,1)=vector1(1)*vector2(1)/norm
      structuralTensor(2,1)=vector1(2)*vector2(1)/norm
      structuralTensor(3,1)=vector1(3)*vector2(1)/norm
      structuralTensor(1,2)=vector1(1)*vector2(2)/norm
      structuralTensor(2,2)=vector1(2)*vector2(2)/norm
      structuralTensor(3,2)=vector1(3)*vector2(2)/norm
      structuralTensor(1,3)=vector1(1)*vector2(3)/norm
      structuralTensor(2,3)=vector1(2)*vector2(3)/norm
      structuralTensor(3,3)=vector1(3)*vector2(3)/norm
    CASE DEFAULT
      localError="The number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Structural tensor:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of dimensions  = ",numberOfDimensions,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Vector 1:",err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDimensions,numberOfDimensions,numberOfDimensions, &
        & vector1,'("    v1(:)  :",3(X,E13.6))','(12X,3(X,E13.6))',err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Vector 2:",err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDimensions,numberOfDimensions,numberOfDimensions, &
        & vector2,'("    v2(:)  :",3(X,E13.6))','(12X,3(X,E13.6))',err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Structural tensor:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDimensions,1,1,numberOfDimensions, &
        & numberOfDimensions,numberOfDimensions,structuralTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
        & '("    S','(",I1,",:)',' :",3(X,E13.6))','(12X,3(X,E13.6))',err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_StructuralTensorComponentCalculate")
    RETURN
999 ERRORS("FiniteElasticity_StructuralTensorComponentCalculate",err,error)
    EXITS("FiniteElasticity_StructuralTensorComponentCalculate")
    RETURN 1

  END SUBROUTINE FiniteElasticity_StructuralTensorComponentCalculate

  !
  !================================================================================================================================
  !

  
END MODULE FiniteElasticityUtilityRoutines
