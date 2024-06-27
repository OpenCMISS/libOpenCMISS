MODULE OC

  PRIVATE

  !>An example oc_ type.
  TYPE oc_ExampleType
  END TYPE oc_ExampleType

  INTERFACE oc_Example_SomeInterface
    MODULE PROCEDURE oc_Example_SomeInterfaceNumber
    MODULE PROCEDURE oc_Example_SomeInterfaceObj
  END INTERFACE !oc_Example_SomeInterface

  INTERFACE oc_Example_CreateStart
    MODULE PROCEDURE oc_Example_CreateStartObj
    MODULE PROCEDURE oc_Example_CreateStartNumber
  END INTERFACE !oc_Example_CreateStart

  INTERFACE oc_ArrayRoutine
    MODULE PROCEDURE oc_ArrayRoutine0
    MODULE PROCEDURE oc_ArrayRoutine1
  END INTERFACE !oc_ArrayRoutine

  INTERFACE oc_StringRoutine
    MODULE PROCEDURE oc_StringRoutineCObj
    MODULE PROCEDURE oc_StringRoutineVSObj
    MODULE PROCEDURE oc_StringRoutineCNumber
    MODULE PROCEDURE oc_StringRoutineVSNumber
  END INTERFACE !oc_StringRoutine

  !> \addtogroup OpenCMISS_ExampleEnum OpenCMISS::ExampleEnum
  !> \brief Example of an enum
  !>@{
  INTEGER(INTG), PARAMETER :: OC_ENUM_ONE = FIRST_VALUE !<Description of first enum value
  INTEGER(INTG), PARAMETER :: OC_ENUM_TWO = SECOND_VALUE !<Description of second enum value
  INTEGER(INTG), PARAMETER :: OC_ENUM_THREE = THIRD_VALUE !<Description of third enum value
  !>@}

  INTEGER(INTG), PARAMETER :: UNGROUPED_CONSTANT = 1 !<Description

  INTEGER(INTG), PARAMETER :: NON_PUBLIC_CONSTANT = 1

  PUBLIC oc_ExampleType, oc_Example_SomeInterface, oc_Example_CreateStart, &
    & oc_Example_Initialise, oc_StringRoutine, oc_ArrayRoutine, OC_ENUM_ONE, &
    & OC_ENUM_TWO, OC_ENUM_THREE, UNGROUPED_CONSTANT

CONTAINS

  !>Doxygen comment describing subroutine
  SUBROUTINE oc_Example_SomeInterfaceObj(Example, InputString, OutputString, InputArray, &
      & OutputArray, InputArray2D, OutputArray2D, ArrayWithSize, InputReal, OutputReal, Err)
    TYPE(oc_ExampleType), INTENT(INOUT) :: Example !<Comment for Example
    CHARACTER(LEN=*), INTENT(IN) :: InputString !<Comment for InputString
    CHARACTER(LEN=*), INTENT(OUT) :: OutputString !<Comment for OutputString
    INTEGER(INTG), INTENT(IN) :: InputArray(:) !<Comment for InputArray
    INTEGER(INTG), INTENT(OUT) :: OutputArray(:) !<Comment for OutputArray
    INTEGER(INTG), INTENT(IN) :: InputArray2D(:,:) !<Comment for InputArray2D
    INTEGER(INTG), INTENT(OUT) :: OutputArray2D(:,:) !<Comment for OutputArray2D
    INTEGER(INTG), INTENT(IN) :: ArrayWithSize(2) !<Comment for ArrayWithSize
    REAL(DP), INTENT(IN) :: InputReal !<Comment for InputReal
    REAL(DP), INTENT(OUT) :: OutputReal !<Comment for OutputReal
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_Example_SomeInterfaceObj

  SUBROUTINE oc_Example_SomeInterfaceNumber(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_Example_SomeInterfaceNumber

  SUBROUTINE oc_StringRoutineCObj(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_StringRoutineCObj

  SUBROUTINE oc_StringRoutineVSObj(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_StringRoutineVSObj

  SUBROUTINE oc_StringRoutineCNumber(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_StringRoutineCNumber

  SUBROUTINE oc_StringRoutineVSNumber(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_StringRoutineVSNumber

  SUBROUTINE oc_ArrayRoutine0(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_ArrayRoutine0

  SUBROUTINE oc_ArrayRoutine1(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_ArrayRoutine1

  SUBROUTINE oc_Example_Initialise(Example, Err)
    TYPE(oc_ExampleType), INTENT(OUT) :: Example
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_Example_Initialise

  SUBROUTINE oc_Example_CreateStartObj(UserNumber, Example, Err)
    INTEGER(INTG), INTENT(IN) :: UserNumber
    TYPE(oc_ExampleType), INTENT(INOUT) :: Example
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_Example_CreateStartObj

  SUBROUTINE oc_Example_CreateStartNumber(UserNumber, Err)
    INTEGER(INTG), INTENT(IN) :: UserNumber
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_Example_CreateStartNumber

  SUBROUTINE oc_NonPublicRoutine(Err)
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code.
  END SUBROUTINE oc_NonPublicRoutine

END MODULE OC
