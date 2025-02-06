!> \file
!> \brief Fortran interface definitions for CellML model routines (in C++)
!>
!> The Original Code is OpenCMISS
!>
!> Please see LICENSE for license information
!>
!> Contributor(s): Chris Bradley, David Nickerson
!>
MODULE OpenCMISSCellMLModel


  USE BaseRoutines
  USE Constants
  USE ISO_C_BINDING, ONLY : C_INT, C_CHAR, C_PTR, C_NULL_PTR, C_ASSOCIATED
  USE ISO_VARYING_STRING
  USE Kinds
  USE OpenCMISSFortranC
  
#include "macros.h"
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !THESE NEED TO BE KEPT UP-TO-DATE WITH THOSE IN opencmiss_cellml_model.hpp
  
  INTEGER(C_INT), PARAMETER :: CELLML_MODEL_NO_ERROR = 0
  INTEGER(C_INT), PARAMETER :: CELLML_MODEL_COULD_NOT_CREATE_MODEL_ERROR = -9995
  INTEGER(C_INT), PARAMETER :: CELLML_MODEL_INVALID_URI_ERROR = -9996
  INTEGER(C_INT), PARAMETER :: CELLML_MODEL_PTR_NOT_NULL_ERROR = -9997
  INTEGER(C_INT), PARAMETER :: CELLML_MODEL_PTR_IS_NULL_ERROR = -9998
  INTEGER(C_INT), PARAMETER :: CELLML_MODEL_ERROR = -9999
  
  !Module types

  !Module variables

  !Module interface
  
  INTERFACE

    !>Fortran binding to the C++ CellML create model routine
    INTEGER(C_INT) FUNCTION CellMLModel_Create(uri, model) &
      & BIND(C, NAME='CellMLModel_CreateF')
      
      USE ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_PTR
      
      !Argument variables
      CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: uri(*) !<The uri string to the CellML model.
      TYPE(C_PTR), INTENT(INOUT) :: model !<On return the newly creted CellMLModel object.
      
    END FUNCTION CellMLModel_Create
    
    !>Fortran binding to the C++ CellML destroy model routine
    INTEGER(C_INT) FUNCTION CellMLModel_Destroy(model) &
      & BIND(C, NAME='CellMLModel_DestroyF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT
      
      !Argument variables
      TYPE(C_PTR), INTENT(INOUT) :: model !<A pointer to the CellML model to destroy
      
    END FUNCTION CellMLModel_Destroy

    !>Fortran binding to the C++ CellML model get initial value routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetInitialValue(model, variableName, initialValue) &
      & BIND(C,NAME='CellMLModel_GetInitialValueF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE, C_CHAR
      
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellMLModel to get the initial value from.
      CHARACTER(KIND=C_CHAR) :: variableName(*) !<The name of the variable to get the initial value for
      REAL(C_DOUBLE) :: initialValue !<On return, the initial value of the variable.
      
    END FUNCTION CellMLModel_GetInitialValue

    !>Fortran binding to the C++ CellML model get initial value by index routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetInitialValueByIndex(model, cellMLVariableType, variableIndex, initialValue) &
      & BIND(C,NAME='CellMLModel_GetInitialValueByIndexF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE
      
      !Argument variables
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellMLModel to get the initial value for
      INTEGER(C_INT), INTENT(IN), VALUE :: cellmlVariableType !<The type of CellML variable to get the initial value for
      INTEGER(C_INT), INTENT(IN), VALUE :: variableIndex !<The index of the CellML variable to get the initial value for
      REAL(C_DOUBLE), INTENT(OUT) :: initialValue !<On return, the initial value of the specified CellML variable
      
    END FUNCTION CellMLModel_GetInitialValueByIndex
    
    !>Fortran binding to the C++ CellML model get CellML variable type routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetVariableType(model, variableName, cellMLVariableType) &
      & BIND(C,NAME='CellMLModel_GetVariableTypeF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_CHAR
      
      !Argument variables
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to get the CellML variable type for
      CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: variableName(*) !<The variable name to get the CellML variable type for
      INTEGER(C_INT), INTENT(OUT) :: cellMLVariableType !<On return, the CellML variable type
      
    END FUNCTION CellMLModel_GetVariableType

    !>Fortran binding to the C++ CellML get CellML variable index routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetVariableIndex(model, variableName, variableIndex) &
      & BIND(C,NAME='CellMLModel_GetVariableIndexF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_CHAR
      
      !Argument variables
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to get the CellML variable index for
      CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: variableName(*) !<The variable name to the CellML variable index for
      INTEGER(C_INT), INTENT(OUT) :: variableIndex !<On return, the CellML variable index
      
    END FUNCTION CellMLModel_GetVariableIndex

    !>Fortran binding to the C++ CellML set variable as known routine
    INTEGER(C_INT) FUNCTION CellMLModel_SetVariableAsKnown(model, variableName) &
      & BIND(C,NAME='CellMLModel_SetVariableAsKnownF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_CHAR
      
      !Argument variables
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to set a variable as known for
      CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: variableName(*) !<The variable name to set the variable as known for
      
    END FUNCTION CellMLModel_SetVariableAsKnown
     
    !>Fortran binding to the C++ CellML set variable as known routine
    INTEGER(C_INT) FUNCTION CellMLModel_SetVariableAsWanted(model, variableName) &
      & BIND(C,NAME='CellMLModel_SetVariableAsWantedF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_CHAR
      
      !Argument variables
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to set a variable as wanted for
      CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: variableName(*) !<The variable name to set the variable as wanted for
      
    END FUNCTION CellMLModel_SetVariableAsWanted
     
    !>Fortran binding to the C++ CellML instatiate routine
    INTEGER(C_INT) FUNCTION CellMLModel_Instatiate(model) &
      & BIND(C,NAME='CellMLModel_InstantiateF')
      
      USE :: ISO_C_BINDING, ONLY: C_PTR, C_INT
      
      !Argument variables
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to instantiate
      
    END FUNCTION CellMLModel_Instatiate

    
    !>Fortran binding to the C++ CellML model get the number of algebraics routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetNumberOfAlgebraics(model, numberOfAlgebraic) &
      & BIND(C,NAME='CellMLModel_GetNumberOfAlgebraicsF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT
      
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to get the number of algebraics for
      INTEGER(C_INT), INTENT(OUT) :: numberOfAlgebraic !<On return, the number of algebraic variables.
      
    END FUNCTION CellMLModel_GetNumberOfAlgebraics
     
    !>Fortran binding to the C++ CellML model get the number of computed constants routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetNumberOfComputedConstants(model, numberOfComputedConstants) &
      & BIND(C,NAME='CellMLModel_GetNumberOfComputedConstantsF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT
      
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to get the number of computed constants for
      INTEGER(C_INT), INTENT(OUT) :: numberOfComputedConstants !<On return, the number of computed constants for the CellML model.
      
    END FUNCTION CellMLModel_GetNumberOfComputedConstants
    
    !>Fortran binding to the C++ CellML model get the number of constants routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetNumberOfConstants(model, numberOfConstants) &
      & BIND(C,NAME='CellMLModel_GetNumberOfConstantsF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT
      
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to get the number of constants for
      INTEGER(C_INT), INTENT(OUT) :: numberOfConstants !<On return, the number of constants for the CellML model.
      
    END FUNCTION CellMLModel_GetNumberOfConstants
     
    !>Fortran binding to the C++ CellML model get the number of rates routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetNumberOfRates(model, numberOfRates) &
      & BIND(C,NAME='CellMLModel_GetNumberOfRatesF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT
      
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to get the number of rates for
      INTEGER(C_INT), INTENT(OUT) :: numberOfRates !<On return, the number of rates 
      
    END FUNCTION CellMLModel_GetNumberOfRates
         
    !>Fortran binding to the C++ CellML model get the number of states routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetNumberOfStates(model, numberOfStates) &
      & BIND(C,NAME='CellMLModel_GetNumberOfStatesF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT
      
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to get the number of states for
      INTEGER(C_INT), INTENT(OUT) :: numberOfStates !<On return, the number of states
      
    END FUNCTION CellMLModel_GetNumberOfStates
    
    SUBROUTINE CellMLModel_CallComputedConstantsRoutine(model, variables) &
      & BIND(C,NAME='CellMLModel_CallComputedConstantsRoutineF')
      
      USE ISO_C_BINDING, ONLY: C_PTR,C_DOUBLE
      
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to call the computed constants routine for
      REAL(C_DOUBLE), INTENT(INOUT) :: variables(*) !<The array of variables, updated on exit
      
    END SUBROUTINE CellMLModel_CallComputedConstantsRoutine
    
    SUBROUTINE CellMLModel_CallRatesRoutine(model, variableOfIntegration, states, rates, wanted, known) &
      & BIND(C,NAME='CellMLModel_CallRatesRoutineF')
      
      USE ISO_C_BINDING, ONLY: C_PTR,C_DOUBLE
      
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to call the rates routine for
      REAL(C_DOUBLE), INTENT(IN), VALUE :: variableOfIntegration !<The value of the variable of integration
      REAL(C_DOUBLE), INTENT(INOUT) :: states(*) !<The array of state variables, updated on exit
      REAL(C_DOUBLE), INTENT(OUT)  :: rates(*) !<On return, the array of rates
      REAL(C_DOUBLE), INTENT(OUT) :: wanted(*) !<On return, the array of wanted 
      REAL(C_DOUBLE), INTENT(IN) :: known(*) !<The array of known variables
      
    END SUBROUTINE CellMLModel_CallRatesRoutine
    
   SUBROUTINE CellMLModel_CallVariablesRoutine(model, variableOfIntegration, states, rates, variables) &
      & BIND(C,NAME='CellMLModel_CallVariablesRoutineF')
      
      USE ISO_C_BINDING, ONLY: C_PTR,C_DOUBLE
      
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to call the variables routine for
      REAL(C_DOUBLE), INTENT(IN), VALUE :: variableOfIntegration !<The value of the variable of integration
      REAL(C_DOUBLE), INTENT(INOUT) :: states(*) !<The array of state variables, updated on exit
      REAL(C_DOUBLE), INTENT(INOUT)  :: rates(*) !<The array of rates.
      REAL(C_DOUBLE), INTENT(INOUT) :: variables(*) !<The array of variables, update on exit.
      
    END SUBROUTINE CellMLModel_CallVariablesRoutine
    
    !>Fortran binding to the C++ CellML model get last error routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetLastError(model, maxErrorStringLength, errorString) &
      & BIND(C,NAME='CellMLModel_GetLastErrorF')
      
      USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_CHAR
      
      TYPE(C_PTR), INTENT(IN), VALUE :: model !<A pointer to the CellML model to get the last error for
      INTEGER(C_INT), INTENT(IN), VALUE :: maxErrorStringLength !<The maximum size of the error string
      CHARACTER(LEN=1,KIND=C_CHAR), INTENT(INOUT) :: errorString(*) !<On return, the last error string
      
    END FUNCTION CellMLModel_GetLastError
    
  END INTERFACE

  PUBLIC CELLML_MODEL_NO_ERROR,CELLML_MODEL_COULD_NOT_CREATE_MODEL_ERROR,CELLML_MODEL_INVALID_URI_ERROR, &
    & CELLML_MODEL_PTR_NOT_NULL_ERROR,CELLML_MODEL_PTR_IS_NULL_ERROR,CELLML_MODEL_ERROR

  PUBLIC CellMLModel_Create,CellMLModel_Destroy

  PUBLIC CellMLModel_GetInitialValue,CellMLModel_GetInitialValueByIndex

  PUBLIC CellMLModel_GetVariableIndex,CellMLModel_GetVariableType

  PUBLIC CellMLModel_GetNumberOfAlgebraics,CellMLModel_GetNumberOfComputedConstants, &
    & CellMLModel_GetNumberOfConstants,CellMLModel_GetNumberOfRates,CellMLModel_GetNumberOfStates

  PUBLIC CellMLModel_Instatiate

  PUBLIC CellMLModel_SetVariableAsKnown,CellMLModel_SetVariableAsWanted

  PUBLIC CellMLModel_CallComputedConstantsRoutine,CellMLModel_CallRatesRoutine,CellMLModel_CallVariablesRoutine

  PUBLIC CellMLModel_GetLastError

  PUBLIC CellMLModel_CheckError
  
CONTAINS

  !>Check the CellML model error code and handle any errors.
  SUBROUTINE CellMLModel_CheckError(model,cellMLErrorCode,err,error,*)

    USE ISO_C_BINDING, ONLY : C_INT, C_CHAR, C_PTR
    
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: model !<The CellML model to check the error code for.
    INTEGER(INTG), INTENT(IN) :: cellMLErrorCode !<The CellML error code to check.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    CHARACTER(LEN=1,KIND=C_CHAR) :: cError(MAXSTRLEN)
    CHARACTER(LEN=MAXSTRLEN) :: dummyError
    
    ENTERS("CellMLModel_CheckError",err,error,*999)
    
    SELECT CASE(cellMLErrorCode)
    CASE(CELLML_MODEL_NO_ERROR)
      !No error, do nothing
    CASE(CELLML_MODEL_COULD_NOT_CREATE_MODEL_ERROR)
      CALL FlagError("Could not create new CellML model.",err,error,*999)
    CASE(CELLML_MODEL_INVALID_URI_ERROR)
      CALL FlagError("The URI to the CellML model is invalid.",err,error,*999)
    CASE(CELLML_MODEL_PTR_NOT_NULL_ERROR)
      CALL FlagError("The CellML model pointer for a new CellML model is not null.",err,error,*999)
    CASE(CELLML_MODEL_PTR_IS_NULL_ERROR)
      CALL FlagError("CellML model pointer is null.",err,error,*999)
    CASE DEFAULT
      err = CellMLModel_GetLastError(model,MAXSTRLEN,cError)
      IF(err /= CELLML_MODEL_NO_ERROR) THEN
        CALL OpenCMISSC2FString(cError,dummyError)
        CALL FlagError(dummyError,err,error,*999)
      ENDIF
    END SELECT
    
    EXITS("CellMLModel_CheckError")
    RETURN
999 ERRORSEXITS("CellMLModel_CheckError",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModel_CheckError
  
END MODULE OpenCMISSCellMLModel
