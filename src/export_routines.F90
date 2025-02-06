!> \file
!> \author Chris Bradley
!> \brief This module handles all export related routines.
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

!> This module handles all export related routines.
MODULE ExportRoutines

  USE BaseRoutines
  USE ExportAccessRoutines
  USE FieldAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types
  USE VTKRoutines
  
#include "macros.h"  
  
  IMPLICIT NONE
  
  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  !>Sets/changes the base filename for an export
  INTERFACE Export_BaseFilenameSet
    MODULE PROCEDURE Export_BaseFilenameSetC
    MODULE PROCEDURE Export_BaseFilenameSetVS
  END INTERFACE Export_BaseFilenameSet

  PUBLIC Export_BaseFilenameSet
  
  PUBLIC Export_CreateFinish,Export_CreateStart

  PUBLIC Export_Finalise,Export_Initialise

  PUBLIC Exports_Finalise,Exports_Initialise

  PUBLIC ExportVariable_Finalise,ExportVariable_Initialise

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the base filename for an export for character filenames. \see OpenCMISS::OC_Export_BaseFilenameSet
  SUBROUTINE Export_BaseFilenameSetC(export,baseFilename,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER :: export !<A pointer to the export to set/change the base filename for
    CHARACTER(LEN=*), INTENT(IN) :: baseFilename !<The base filename to set/change
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Export_BaseFilenameSetC",err,error,*999)

    CALL Export_AssertNotFinished(export,err,error,*999)
    
    export%basefilename=baseFilename(1:LEN_TRIM(baseFilename))

    EXITS("Export_BaseFilenameSetC")
    RETURN
999 ERRORSEXITS("Export_BaseFilenameSetC",err,error)
    RETURN 1

  END SUBROUTINE Export_BaseFilenameSetC

  !
  !================================================================================================================================
  !

  !>Sets/changes the base filename for an export for varying string filenames. \see OpenCMISS::OC_Export_BaseFilenameSet
  SUBROUTINE Export_BaseFilenameSetVS(export,baseFilename,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER :: export !<A pointer to the export to set/change the base filename for
    TYPE(VARYING_STRING), INTENT(IN) :: baseFilename !<The base filename to set/change
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Export_BaseFilenameSetVS",err,error,*999)

    CALL Export_AssertNotFinished(export,err,error,*999)
    
    export%baseFilename = baseFilename
 
    EXITS("Export_BaseFilenameSetVS")
    RETURN
999 ERRORSEXITS("Export_BaseFilenameSetVS",err,error)
    RETURN 1

  END SUBROUTINE Export_BaseFilenameSetVS

  !
  !================================================================================================================================
  !

  !>Finishes an export 
  SUBROUTINE Export_CreateFinish(export,err,error,*)
    
    !Argument variables
    TYPE(ExportType), POINTER, INTENT(INOUT) :: export !<The export to finish
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Export_CreateFinish",err,error,*999)

    CALL Export_AssertNotFinished(export,err,error,*999)

    !Finish the different export options    
    SELECT CASE(export%exportFormat)
    CASE(EXPORT_EXFILE_FORMAT)
      CALL FlagError("Not implemented.",err,error,*999)
      !CALL Exfile_CreateFinish(export%exfileExport,err,error,*999)
    CASE(EXPORT_VTK_FORMAT)
      CALL VTKExport_CreateFinish(export%vtkExport,err,error,*999)
    CASE DEFAULT
      localError="The export format of "//TRIM(NumberToVString(export%exportFormat,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Finish the export
    export%exportFinished = .TRUE.
    
    EXITS("Export_CreateFinish")
    RETURN
999 ERRORSEXITS("Export_CreateFinish",err,error)
    RETURN 1

  END SUBROUTINE Export_CreateFinish

  !
  !================================================================================================================================
  !
 
  !>Finishes an export 
  SUBROUTINE Export_CreateStart(exportUserNumber,exports,export,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: exportUserNumber !<The user number of the export to create
    TYPE(ExportsType), POINTER :: exports !<The exports to create the export for
    TYPE(ExportType), POINTER, INTENT(INOUT) :: export !<On return, a pointer of the created export. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,exportIdx
    TYPE(ExportType), POINTER :: newExport
    TYPE(ExportPtrType), ALLOCATABLE :: newExports(:)
    TYPE(VARYING_STRING) :: dummyError,localError

    NULLIFY(newExport)
    
    ENTERS("Export_CreateStart",err,error,*997)
    
    IF(.NOT.ASSOCIATED(exports)) CALL FlagError("Exports is not associated.",err,error,*997)
    IF(ASSOCIATED(export)) CALL FlagError("Export is already associated.",err,error,*997)

    CALL Export_UserNumberFind(exportUserNumber,exports,export,err,error,*997)
    IF(ASSOCIATED(export)) THEN
      localError="Export number "//TRIM(NumberToVString(exportUserNumber,"*",err,error))//" has already been created."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !Initialise the export
    CALL Export_Initialise(newExport,err,error,*999)
    !Set the default values
    newExport%userNumber = exportUserNumber
    newExport%globalNumber = exports%numberOfExports + 1
    newExport%exports => exports
    !Default to VTK
    newExport%exportFormat = EXPORT_VTK_FORMAT
    CALL Export_VTKInitialise(newExport,err,error,*999)
    !Add new export into list of exports
    ALLOCATE(newExports(exports%numberOfExports+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new exports.",err,error,*999)
    DO exportIdx=1,exports%numberOfExports
      newExports(exportIdx)%ptr=>exports%exports(exportIdx)%ptr
    ENDDO !exportIdx
    newExports(exports%numberOfExports+1)%ptr=>newExport
    CALL MOVE_ALLOC(newExports,exports%exports)
    exports%numberOfExports=exports%numberOfExports + 1
    !Return the pointer to the new export
    export=>newExport
    
    EXITS("Export_CreateStart")
    RETURN
999 CALL Export_Finalise(newExport,dummyErr,dummyError,*998)
998 IF(ALLOCATED(newExports)) DEALLOCATE(newExports)
    NULLIFY(export)
997 ERRORSEXITS("Export_CreateStart",err,error)
    RETURN 1

  END SUBROUTINE Export_CreateStart

  !
  !================================================================================================================================
  !

  !>Sets/changes the export format for an export
  SUBROUTINE Export_ExportFormatSet(export,exportFormat,err,error,*)
    
    !Argument variables
    TYPE(ExportType), POINTER, INTENT(INOUT) :: export !<The export to set the export format for
    INTEGER(INTG), INTENT(IN) :: exportFormat !<The export format to set. \see ExportRoutines_ExportFormatTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: previousExportFormat
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Export_ExportFormatSet",err,error,*999)

    CALL Export_AssertNotFinished(export,err,error,*999)

    previousExportFormat = export%exportFormat
    IF(exportFormat /= previousExportFormat) THEN
      !Initialise the new export format
      SELECT CASE(exportFormat)
      CASE(EXPORT_EXFILE_FORMAT)
        CALL Export_ExfileInitialise(export,err,error,*999)
      CASE(EXPORT_VTK_FORMAT)
        CALL Export_VTKInitialise(export,err,error,*999)
      CASE DEFAULT
        localError="The specified export format of "//TRIM(NumberToVString(exportFormat,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finalise the old export format
      SELECT CASE(previousExportFormat)
      CASE(EXPORT_EXFILE_FORMAT)
        CALL Export_ExfileFinalise(export%exfileExport,err,error,*999)
      CASE(EXPORT_VTK_FORMAT)
        CALL Export_VTKFinalise(export%vtkExport,err,error,*999)
      CASE DEFAULT
        localError="The export format of "//TRIM(NumberToVString(previousExportFormat,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Export_ExportFormatSet")
    RETURN
999 ERRORSEXITS("Export_ExportFormatSet",err,error)
    RETURN 1

  END SUBROUTINE Export_ExportFormatSet

  !
  !================================================================================================================================
  !
 
  !>Finishes an export 
  SUBROUTINE Export_ExportVariableAdd(export,fieldVariable,parameterSetType,startComponent,endComponent,exportName, &
    & exportVariableIndex,err,error,*)    
    !Argument variables
    TYPE(ExportType), POINTER :: export !<The export to add the export variable for
    TYPE(FieldVariableType), POINTER, INTENT(IN) :: fieldVariable !<A pointer to the field variable to add as an export variable
    INTEGER(INTG), INTENT(IN) :: parameterSetType !<The parameter set type of the field variable to add as an export variable
    INTEGER(INTG), INTENT(IN) :: startComponent !<The start component number of the variable to export.
    INTEGER(INTG), INTENT(IN) :: endComponent !<The end component number of the variable to export.
    TYPE(VARYING_STRING), INTENT(IN) :: exportName !<The name of the variable in export files.
    INTEGER(INTG), INTENT(OUT) :: exportVariableIndex !<On exit, the index of the added export variable.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: checkParameterSetType,dummyErr,exportVariableIdx
    LOGICAL :: variableOK
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(ExportVariableType), POINTER :: checkExportVariable,newExportVariable
    TYPE(ExportVariablePtrType), ALLOCATABLE :: newExportVariables(:)
    TYPE(FieldParameterSetType), POINTER :: parameterSet
    TYPE(FieldVariableType), POINTER :: checkFieldVariable
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: dummyError,localError

    NULLIFY(newExportVariable)
    
    ENTERS("Export_ExportVariableAdd",err,error,*997)

    CALL Export_AssertNotFinished(export,err,error,*999)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*997)
    
    !Check if the parameter set exists
    NULLIFY(parameterSet)
    CALL FieldVariable_ParameterSetGet(fieldVariable,parameterSetType,parameterSet,err,error,*999)
    !Check the components
    CALL FieldVariable_ComponentNumberCheck(fieldVariable,startComponent,err,error,*999)
    CALL FieldVariable_ComponentNumberCheck(fieldVariable,endComponent,err,error,*999)
    IF(endComponent<startComponent) THEN
      localError="The end component number of "//TRIM(NumberToVString(endComponent,"*",err,error))// &
        & " must be >= than the start component number of "//TRIM(NumberToVString(startComponent,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    !Check that the field variable to add is compatible with the other field variables already added.
    NULLIFY(region)
    CALL FieldVariable_RegionGet(fieldVariable,region,err,error,*999)
    NULLIFY(decomposition)
    CALL FieldVariable_DecompositionGet(fieldVariable,decomposition,err,error,*999)
    variableOK = .TRUE.
    DO exportVariableIdx = 1,export%numberOfExportVariables
      NULLIFY(checkExportVariable)
      CALL Export_ExportVariableIndexGet(export,exportVariableIdx,checkExportVariable,err,error,*999)      
      NULLIFY(checkFieldVariable)
      CALL ExportVariable_FieldVariableGet(checkExportVariable,checkFieldVariable,err,error,*999)
      !Check if we already have that parameter set in the export variable list
      IF(ASSOCIATED(fieldVariable,checkFieldVariable)) THEN
        CALL ExportVariable_ParameterSetTypeGet(checkExportVariable,checkParameterSetType,err,error,*999)
        IF(parameterSetType==checkParameterSetType) THEN
          !Parameter set already added.
          localError="The parameter set type of "//TRIM(NumberToVString(parameterSetType,"*",err,error))// &
            & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
          IF(ASSOCIATED(fieldVariable%field)) localError=localError//" of field number "// &
            & TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
          localError=localError//" has already been added to export number "// &
            & TRIM(NumberToVString(export%userNumber,"*",err,error))//" at index number "// &
            & TRIM(NumberToVString(exportVariableIdx,"*",err,error))//"."
          variableOK = .FALSE.
          EXIT
        ENDIF
      ENDIF
      !Check the field variable is in the same region as the others
      CALL FieldVariable_FieldVariableAssertSameRegion(fieldVariable,checkFieldVariable,err,error,*999)
      !Check the field variable has the same decomposition as the others
      CALL FieldVariable_FieldVariableAssertSameDecomposition(fieldVariable,checkFieldVariable,err,error,*999)
      !Check the field variable is OK for the different exports
      SELECT CASE(export%exportFormat)
      CASE(EXPORT_EXFILE_FORMAT)
        CALL FlagError("Not implemented.",err,error,*999)
        !CALL ExfileExport_CheckVariable(export%exfileExport,fieldVariable,startComponent,endComponent,err,error,*999)
      CASE(EXPORT_VTK_FORMAT)
        CALL VTKExport_CheckVariable(export%vtkExport,fieldVariable,startComponent,endComponent,err,error,*999)
      CASE DEFAULT
        localError="The export format of "//TRIM(NumberToVString(export%exportFormat,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !exportVariableIdx
    IF(.NOT.variableOK) CALL FlagError(localError,err,error,*999)

    !Field variable is OK. Add it to the export variables list.
    
    !Initialise the export
    CALL ExportVariable_Initialise(newExportVariable,err,error,*999)
    !Set the default values
    newExportVariable%exportIndex = export%numberOfExportVariables + 1
    newExportVariable%export => export
    newExportVariable%fieldVariable => fieldVariable
    newExportVariable%parameterSetType = parameterSetType
    newExportVariable%startComponent = startComponent
    newExportVariable%endComponent = endComponent
    newExportVariable%exportName = exportName
    !Add new export variable into list of export variables
    ALLOCATE(newExportVariables(export%numberOfExportVariables+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new export variables.",err,error,*999)
    DO exportVariableIdx=1,export%numberOfExportVariables
      newExportVariables(exportVariableIdx)%ptr=>export%exportVariables(exportVariableIdx)%ptr
    ENDDO !exportVariableIdx
    newExportVariables(export%numberOfExportVariables+1)%ptr=>newExportVariable
    CALL MOVE_ALLOC(newExportVariables,export%exportVariables)
    export%numberOfExportVariables=export%numberOfExportVariables + 1
    !Return the export variable index of the new export
    exportVariableIndex = export%numberOfExportVariables
    
    EXITS("Export_ExportVariableAdd")
    RETURN
999 CALL ExportVariable_Finalise(newExportVariable,dummyErr,dummyError,*998)
998 IF(ALLOCATED(newExportVariables)) DEALLOCATE(newExportVariables)
997 ERRORSEXITS("Export_ExportVariableAdd",err,error)
    RETURN 1

  END SUBROUTINE Export_ExportVariableAdd

  !
  !================================================================================================================================
  !

  !>Finalises an export and deallocates all memory
  SUBROUTINE Export_Finalise(export,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER, INTENT(INOUT) :: export !<The export to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: exportVariableIdx
 
    ENTERS("Export_Finalise",err,error,*999)

    IF(ASSOCIATED(export)) THEN
      DO exportVariableIdx = 1,SIZE(export%exportVariables,1)
        CALL ExportVariable_Finalise(export%exportVariables(exportVariableIdx)%ptr,err,error,*999)
      ENDDO !exportVariableIdx
      DEALLOCATE(export%exportVariables)
      CALL Export_ExfileFinalise(export%exfileExport,err,error,*999)
      CALL Export_VTKFinalise(export%vtkExport,err,error,*999)
      export%baseFilename=""
      DEALLOCATE(export)
    ENDIF
 
    EXITS("Export_Finalise")
    RETURN
999 ERRORSEXITS("Export_Finalise",err,error)
    RETURN 1

  END SUBROUTINE Export_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises an export 
  SUBROUTINE Export_Initialise(export,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER, INTENT(INOUT) :: export !<The export to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Export_Initialise",err,error,*999)

    IF(ASSOCIATED(export)) CALL FlagError("Export is already associated.",err,error,*999)
    
    ALLOCATE(export,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate export.",err,error,*999)
    export%exportFinished=.FALSE.
    export%exportFormat=0
    NULLIFY(export%exfileExport)
    NULLIFY(export%vtkExport)
    NULLIFY(export%controlLoop)
    export%baseFilename=""
    export%numberOfExportVariables=0
  
    EXITS("Export_Initialise")
    RETURN
999 ERRORSEXITS("Export_Initialise",err,error)
    RETURN 1

  END SUBROUTINE Export_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises an exfile export and deallocates all memory
  SUBROUTINE Export_ExfileFinalise(exfileExport,err,error,*)

    !Argument variables
    TYPE(ExfileExportType), POINTER, INTENT(INOUT) :: exfileExport !<The exfile export to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Export_ExfileFinalise",err,error,*999)

    IF(ASSOCIATED(exfileExport)) THEN
      DEALLOCATE(exfileExport)
    ENDIF
 
    EXITS("Export_ExfileFinalise")
    RETURN
999 ERRORSEXITS("Export_ExfileFinalise",err,error)
    RETURN 1

  END SUBROUTINE Export_ExfileFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an exfile export 
  SUBROUTINE Export_ExfileInitialise(export,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER, INTENT(INOUT) :: export !<The export to initialise the exfile export for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Export_ExfileInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(export)) CALL FlagError("Export is not associated.",err,error,*999)
    IF(ASSOCIATED(export%exfileExport)) CALL FlagError("Export exfile export is already associated.",err,error,*999)
    
    ALLOCATE(export%exfileExport,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate exfile export.",err,error,*999)
    export%exportFormat = EXPORT_EXFILE_FORMAT
    export%exfileExport%export=>export
     
    EXITS("Export_ExfileInitialise")
    RETURN
999 ERRORSEXITS("Export_ExfileInitialise",err,error)
    RETURN 1

  END SUBROUTINE Export_ExfileInitialise

  !
  !================================================================================================================================
  !

  !>Finalises an vtk export and deallocates all memory
  SUBROUTINE Export_VTKFinalise(vtkExport,err,error,*)

    !Argument variables
    TYPE(VTKExportType), POINTER, INTENT(INOUT) :: vtkExport !<The VTK export to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Export_VTKFinalise",err,error,*999)

    IF(ASSOCIATED(vtkExport)) THEN
      DEALLOCATE(vtkExport)
    ENDIF
 
    EXITS("Export_VTKFinalise")
    RETURN
999 ERRORSEXITS("Export_VTKFinalise",err,error)
    RETURN 1

  END SUBROUTINE Export_VTKFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an VTK export 
  SUBROUTINE Export_VTKInitialise(export,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER, INTENT(INOUT) :: export !<The export to initialise the VTK export for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Export_VTKInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(export)) CALL FlagError("Export is not associated.",err,error,*999)
    IF(ASSOCIATED(export%vtkExport)) CALL FlagError("Export VTK export is already associated.",err,error,*999)

    !Allocat and initialise a new VTK export
    ALLOCATE(export%vtkExport,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate VTK export.",err,error,*999)
    export%exportFormat = EXPORT_VTK_FORMAT
    export%vtkExport%export=>export
 
    EXITS("Export_VTKInitialise")
    RETURN
999 ERRORSEXITS("Export_VTKInitialise",err,error)
    RETURN 1

  END SUBROUTINE Export_VTKInitialise

  !
  !================================================================================================================================
  !

  !>Finalises an exports and deallocates all memory
  SUBROUTINE Exports_Finalise(exports,err,error,*)

    !Argument variables
    TYPE(ExportsType), POINTER, INTENT(INOUT) :: exports !<The exports to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: exportIdx
 
    ENTERS("Exports_Finalise",err,error,*999)

    IF(ASSOCIATED(exports)) THEN
      DO exportIdx = 1,SIZE(exports%exports,1)
        CALL Export_Finalise(exports%exports(exportIdx)%ptr,err,error,*999)
      ENDDO !exportIdx
      DEALLOCATE(exports%exports)
      DEALLOCATE(exports)
    ENDIF
 
    EXITS("Exports_Finalise")
    RETURN
999 ERRORSEXITS("Exports_Finalise",err,error)
    RETURN 1

  END SUBROUTINE Exports_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises an exports 
  SUBROUTINE Exports_Initialise(exports,err,error,*)

    !Argument variables
    TYPE(ExportsType), POINTER, INTENT(INOUT) :: exports !<The exports to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Exports_Initialise",err,error,*999)

    IF(ASSOCIATED(exports)) CALL FlagError("Exports is already associated.",err,error,*999)
    
    ALLOCATE(exports,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate exports.",err,error,*999)
    exports%numberOfExports=0
  
    EXITS("Exports_Initialise")
    RETURN
999 ERRORSEXITS("Exports_Initialise",err,error)
    RETURN 1

  END SUBROUTINE Exports_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises an export variable and deallocates all memory
  SUBROUTINE ExportVariable_Finalise(exportVariable,err,error,*)

    !Argument variables
    TYPE(ExportVariableType), POINTER, INTENT(INOUT) :: exportVariable !<The export variable to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ExportVariable_Finalise",err,error,*999)

    IF(ASSOCIATED(exportVariable)) THEN
      exportVariable%exportName=""
      DEALLOCATE(exportVariable)
    ENDIF
 
    EXITS("ExportVariable_Finalise")
    RETURN
999 ERRORSEXITS("ExportVariable_Finalise",err,error)
    RETURN 1

  END SUBROUTINE ExportVariable_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises an export variable
  SUBROUTINE ExportVariable_Initialise(exportVariable,err,error,*)

    !Argument variables
    TYPE(ExportVariableType), POINTER, INTENT(INOUT) :: exportVariable !<The export variable to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ExportVariable_Initialise",err,error,*999)

    IF(ASSOCIATED(exportVariable)) CALL FlagError("Export variable is already associated.",err,error,*999)
    
    ALLOCATE(exportVariable,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate export variable.",err,error,*999)
    NULLIFY(exportVariable%export)
    exportVariable%exportIndex=0
    NULLIFY(exportVariable%fieldVariable)
    exportVariable%parameterSetType=0
    NULLIFY(exportVariable%parameterSet)
    exportVariable%startComponent=0
    exportVariable%endComponent=0
    exportVariable%exportName=""
 
    EXITS("ExportVariable_Initialise")
    RETURN
999 ERRORSEXITS("ExportVariable_Initialise",err,error)
    RETURN 1

  END SUBROUTINE ExportVariable_Initialise

  !
  !================================================================================================================================
  !

END MODULE ExportRoutines
