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

  USE, INTRINSIC :: ISO_C_BINDING
  
#include "macros.h"
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters
  
  !Module types

  !Module variables

  !Module interface
  
  INTERFACE

    !>Fortran binding to the C++ CellML create model routine
    FUNCTION CellMLModel_Create(uri) &
      & BIND(C, NAME='cellml_model_create_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR

      !Argument variables
      CHARACTER, INTENT(IN) :: uri(*) !<The uri string to the CellML model
      !Function return variable
      TYPE(C_PTR) :: CellML_CreateModel !<On return, a pointer to the created CellML model
      
    END FUNCTION CellMLModel_Create

    !>Fortran binding to the C++ CellML destroy model routine
    SUBROUTINE CellMLModel_Destroy(model) &
      & BIND(C, NAME='cellml_model_destroy_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
      
      !Argument variables
      TYPE(C_PTR), INTENT(INOUT) :: model !<A pointer to the CellML model to destroy
      
    END SUBROUTINE CellMLModel_Destroy

    !>Fortran binding to the C++ CellML model get initial value routine
    INTEGER (C_INT) FUNCTION CellMLModel_GetInitialValue(model,variableName,initialValue) &
      & BIND(C,NAME='cellml_model_get_initial_value_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT,C_DOUBLE
      
      TYPE (C_PTR), VALUE :: model
      CHARACTER, DIMENSION(*) :: variableName
      REAL (C_DOUBLE) :: initialValue
      
    END FUNCTION CellMLModel_GetInitialValue

    !>Fortran binding to the C++ CellML model get initial value by index routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetInitialValueByIndex(model,cellMLVariableType,variableIndex,initialValue) &
      & BIND(C,NAME='cellml_model_get_initial_value_by_index_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT,C_DOUBLE
      
      !Argument variables
      TYPE(C_PTR), VALUE :: model !<A pointer to the CellML model to get the initial value for
      INTEGER(C_INT), INTENT(IN) :: cellmlVariableType !<The type of CellML variable to get the initial value for
      INTEGER(C_INT), INTENT(IN) :: variableIndex !<The index of the CellML variable to get the initial value for
      REAL(C_DOUBLE), INTENT(OUT) :: initialValue !<On return, the initial value of the specified CellML variable
      
    END FUNCTION CellMLModel_GetInitialValueByIndex
    
    !>Fortran binding to the C++ CellML model get CellML variable type routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetVariableType(model,variableName,cellMLVariableType) &
      & BIND(C,NAME='cellml_model_get_variable_type_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT
      
      !Argument variables
      TYPE(C_PTR), VALUE :: model !<A pointer to the CellML model to get the CellML variable type for
      CHARACTER, INTENT(IN) :: variableName(*) !<The variable name to get the CellML variable type for
      INTEGER(C_INT), INTENT(OUT) :: cellMLVariableType !<On return, the CellML variable type
      
    END FUNCTION CellMLModel_GetVariableType

    !>Fortran binding to the C++ CellML get CellML variable index routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetVariableIndex(model,variableName,variableIndex) &
      & BIND(C,NAME='cellml_model_get_variable_index_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT
      
      !Argument variables
      TYPE (C_PTR), VALUE :: model !<A pointer to the CellML model to get the CellML variable index for
      CHARACTER, INTENT(IN) :: variableName(*) !<The variable name to the CellML variable index for
      INTEGER(C_INT), INTENT(OUT) :: variableIndex !<On return, the CellML variable index
      
    END FUNCTION CellMLModel_GetVariableIndex

    !>Fortran binding to the C++ CellML set variable as known routine
    INTEGER(C_INT) FUNCTION CellMLModel_SetVariableAsKnown(model,variableName) &
      & BIND(C,NAME='cellml_model_set_variable_as_known_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT
      
      !Argument variables
      TYPE(C_PTR), VALUE :: model !<A pointer to the CellML model to set a variable as known for
      CHARACTER, INTENT(IN) :: variableName(*) !<The variable name to set the variable as known for
      
    END FUNCTION CellMLModel_SetVariableAsKnown
     
    !>Fortran binding to the C++ CellML set variable as known routine
    INTEGER(C_INT) FUNCTION CellMLModel_SetVariableAsWanted(model,variableName) &
      & BIND(C,NAME='cellml_model_set_variable_as_wanted_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT
      
      !Argument variables
      TYPE(C_PTR), VALUE :: model !<A pointer to the CellML model to set a variable as wanted for
      CHARACTER, INTENT(IN) :: variableName(*) !<The variable name to set the variable as wanted for
      
    END FUNCTION CellMLModel_SetVariableAsWanted
     
    !>Fortran binding to the C++ CellML save temporary files flag routine
    SUBROUTINE CellMLModel_SetSaveTemporaryFiles(model,saveTempState) &
      & BIND(C,NAME='cellml_model_set_save_temp_files_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT
      
      !Argument variables
      TYPE(C_PTR), VALUE :: model !<A pointer to the CellML model to set the save temporary files flag for
      INTEGER(C_INT), INTENT(IN), VALUE :: saveTempState !<The flag for saving temporary flags.
      
    END SUBROUTINE CellMLModel_SetSaveTemporaryFiles
     
    !>Fortran binding to the C++ CellML get temporary files flag routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetSaveTemporaryFiles(model) &
      & BIND(C,NAME='cellml_model_get_save_temp_files_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT
      
      !Argument variables
      TYPE(C_PTR), VALUE :: model !<A pointer to the CellML model to get the save temporary files flag for.
      
    END FUNCTION CellMLModel_GetSaveTemporaryFiles

    !>Fortran binding to the C++ CellML instatiate routine
    INTEGER(C_INT) FUNCTION CellMLModel_Instatiate(model) &
      & BIND(C,NAME='cellml_model_instantiate_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT
      
      !Argument variables
      TYPE(C_PTR), VALUE :: model !<A pointer to the CellML model to instata      !Argument variables
      
    END FUNCTION CellMLModel_Instatiate

    !>Fortran binding to the C++ CellML model get the number of constants routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetNumberConstants(model) &
      & BIND(C,NAME='cellml_model_get_number_constants_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT
      
      TYPE(C_PTR), VALUE :: model !<A pointer to the CellML model to get the number of constants for
      
    END FUNCTION CellMLModel_GetNumberConstants
     
    !>Fortran binding to the C++ CellML model get the number of rates routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetNumberRates(model) &
      & BIND(C,NAME='cellml_model_get_number_rates_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT
      
      TYPE(C_PTR), VALUE :: model !<A pointer to the CellML model to get the number of rates for
      
    END FUNCTION CellMLModel_GetNumberRates

    
    !>Fortran binding to the C++ CellML model get the number of algebraic routine
    INTEGER(C_INT) FUNCTION CellMLModel_GetNumberAlgebraic(model) &
      & BIND(C,NAME='cellml_model_get_number_algebraic_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_INT
      
      TYPE(C_PTR), VALUE :: model !<A pointer to the CellML model to get the number of algebraic for
      
    END FUNCTION CellMLModel_GetNumberAlgebraic
     
    SUBROUTINE CellMLModel_CallRHSRoutine(model,variableOfIntegration,states,rates,wanted,known) &
      & BIND(C,NAME='cellml_model_call_rhs_routine_f')
      
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_DOUBLE
      
      TYPE(C_PTR), VALUE :: model !<A pointer to the CellML model to call the RHS routine for
      REAL(C_DOUBLE), INTENT(IN), VALUE :: variableOfIntegration !<The value of the variable of integration
      REAL(C_DOUBLE), INTENT(INOUT) :: states(*) !<The array of state variables, updated on exit
      REAL(C_DOUBLE), INTENT(OUT)  :: rates(*) !<On return, the array of rates
      REAL(C_DOUBLE), INTENT(OUT) :: wanted(*) !<On return, the array of wanted 
      REAL(C_DOUBLE), INTENT(IN) :: known(*) !<The array of known variables
      
    END SUBROUTINE CellMLModel_CallRHSRoutine
    
  END INTERFACE
  
CONTAINS
  
END MODULE OpenCMISSCellMLModel
