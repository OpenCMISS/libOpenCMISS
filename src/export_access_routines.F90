!> \file
!> \author Chris Bradley
!> \brief This module handles all export access related routines.
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

!> This module handles all export access related routines.
MODULE ExportAccessRoutines

  USE BaseRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup ExportRoutines_ExportFormatTypes Export::ExportFormatTypes
  !> \brief Export format types
  !> \see ExportRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: EXPORT_EXFILE_FORMAT = 1 !<Exfile (cmgui) format type \see ExportRoutines_ExportFormatTypes,ExportRoutines
  INTEGER(INTG), PARAMETER :: EXPORT_VTK_FORMAT = 2 !<VTK format type \see ExportRoutines_ExportFormatTypes,ExportRoutines
  !>@}

  !Module types

  !Module variables

  !Interfaces

  !>Gets the label for a field.
  INTERFACE Export_BasefilenameGet
    MODULE PROCEDURE Export_BaseFilenameGetC
    MODULE PROCEDURE Export_BaseFilenameGetVS
  END INTERFACE Export_BasefilenameGet
 
  PUBLIC EXPORT_EXFILE_FORMAT,EXPORT_VTK_FORMAT

  PUBLIC Export_AssertIsFinished,Export_AssertNotFinished

  PUBLIC Export_BaseFilenameGet

  PUBLIC Export_ExfileExportGet

  PUBLIC Export_ExportFormatGet

  PUBLIC Export_ExportVariableIndexGet
  
  PUBLIC Export_UserNumberFind

  PUBLIC Export_VTKExportGet

  PUBLIC ExportVariable_FieldVariableGet

  PUBLIC ExportVariable_ParameterSetTypeGet

  PUBLIC VTKExport_DomainGet

  PUBLIC VTKExport_ExportGet

CONTAINS
  
  !
  !=================================================================================================================================
  !

  !>Assert that a export has been finished
  SUBROUTINE Export_AssertIsFinished(export,err,error,*)

    !Argument Variables
    TYPE(ExportType), POINTER, INTENT(IN) :: export !<The export to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Export_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(export)) CALL FlagError("Export is not associated.",err,error,*999)
#endif    

    IF(.NOT.export%exportFinished) THEN
      localError="Export number "//TRIM(NumberToVString(export%userNumber,"*",err,error))//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Export_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Export_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Export_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a export has not been finished
  SUBROUTINE Export_AssertNotFinished(export,err,error,*)

    !Argument Variables
    TYPE(ExportType), POINTER, INTENT(IN) :: export !<The export to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Export_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(export)) CALL FlagError("Export is not associated.",err,error,*999)
#endif    

    IF(export%exportFinished) THEN
      localError="Export number "//TRIM(NumberToVString(export%userNumber,"*",err,error))//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Export_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Export_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Export_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Gets the export base filename for an export for character filenames
  SUBROUTINE Export_BaseFilenameGetC(export,baseFilename,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER, INTENT(IN) :: export !<The export to get the export base filename for.
    CHARACTER(LEN=*), INTENT(OUT) :: baseFilename !<On return, the base filename for the export.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cLength,vsLength
 
    ENTERS("Export_BaseFilenameGetC",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(export)) CALL FlagError("Export is not associated.",err,error,*999)
#endif    

    !Get the export base filename
    cLength=LEN(baseFilename)
    vsLength=LEN_TRIM(export%baseFilename)
    IF(cLength>vsLength) THEN
      baseFilename = CHAR(LEN_TRIM(export%baseFilename))
    ELSE
      baseFilename = CHAR(export%baseFilename,cLength)
    ENDIF
    
    EXITS("Export_BaseFilenameGetC")
    RETURN
999 ERRORSEXITS("Export_BaseFilenameGetC",err,error)
    RETURN 1
    
  END SUBROUTINE Export_BaseFilenameGetC

  !
  !================================================================================================================================
  !

  !>Gets the export base filename for an export for varying string filenames
  SUBROUTINE Export_BaseFilenameGetVS(export,baseFilename,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER, INTENT(IN) :: export !<The export to get the export base filename for.
    TYPE(VARYING_STRING), INTENT(OUT) :: baseFilename !<On return, the base filename for the export.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Export_BaseFilenameGetVS",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(export)) CALL FlagError("Export is not associated.",err,error,*999)
#endif    

    !Get the export base filename
    baseFilename = export%baseFilename
    
    EXITS("Export_BaseFilenameGetVS")
    RETURN
999 ERRORSEXITS("Export_BaseFilenameGetVS",err,error)
    RETURN 1
    
  END SUBROUTINE Export_BaseFilenameGetVS

  !
  !================================================================================================================================
  !

  !>Gets the exfile export information from an export
  SUBROUTINE Export_ExfileExportGet(export,exfileExport,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER, INTENT(IN) :: export !<The export to get the exfile export for.
    TYPE(ExfileExportType), POINTER, INTENT(INOUT) :: exfileExport !<On exit, a pointer to the exfile export for the export. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Export_ExfileExportGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(exfileExport)) CALL FlagError("Exfile export is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(export)) CALL FlagError("Export is not associated.",err,error,*999)
#endif    

    !Get the export exfile export
    exfileExport=>export%exfileExport

#ifdef WITH_POSTCHECKS    
    !Check exfile export is associated.
    IF(.NOT.ASSOCIATED(exfileExport)) THEN
      localError="Export exfile export is not associated for export number "// &
        & TRIM(NumberToVString(export%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Export_ExfileExportGet")
    RETURN
999 NULLIFY(exfileExport)    
998 ERRORSEXITS("Export_ExfileExportGet",err,error)
    RETURN 1
    
  END SUBROUTINE Export_ExfileExportGet

  !
  !================================================================================================================================
  !

  !>Gets the export format from an export
  SUBROUTINE Export_ExportFormatGet(export,exportFormat,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER, INTENT(IN) :: export !<The export to get the export format for.
    INTEGER(INTG), INTENT(OUT) :: exportFormat !<On return, the export format for the export. \see ExportRoutines_ExportFormatTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Export_ExportFormatGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(export)) CALL FlagError("Export is not associated.",err,error,*999)
#endif    

    !Get the export format
    exportFormat=export%exportFormat
    
    EXITS("Export_ExportFormatGet")
    RETURN
999 ERRORSEXITS("Export_ExportFormatGet",err,error)
    RETURN 1
    
  END SUBROUTINE Export_ExportFormatGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a export variable of a specified index
  SUBROUTINE Export_ExportVariableIndexGet(export,variableIdx,exportVariable,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER :: export !<A pointer to the export to get the export variable for.
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The index of the export variable to set.
    TYPE(ExportVariableType), POINTER :: exportVariable !<On exit, a pointer to the export variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("Export_ExportVariableIndexGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(exportVariable)) CALL FlagError("Export variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(export)) CALL FlagError("Export is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(export%exportVariables)) CALL FlagError("Export variables is not allocated.",err,error,*999)
    IF(variableIdx<0.OR.variableIdx>SIZE(export%exportVariables,1)) THEN
      localError="The specified export variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid. The export variable index must be between 1 and "// &
        & TRIM(NumberToVString(SIZE(export%exportVariables,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    exportVariable=>export%exportVariables(variableIdx)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(exportVariable)) THEN
      localError="The export variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is not associated for export number "//TRIM(NumberToVString(export%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("Export_ExportVariableIndexGet")
    RETURN
#ifdef WITH_CHECKS    
999 NULLIFY(exportVariable)
#endif    
998 ERRORSEXITS("Export_ExportVariableIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE Export_ExportVariableIndexGet

  !
  !================================================================================================================================
  !
  
  !>Finds and returns a pointer to the export identified by user number in the given list of exports. If no export with that number exits mesh is left nullified.   
  SUBROUTINE Export_UserNumberFind(userNumber,exports,export,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the export to find
    TYPE(ExportsType), POINTER :: exports !<The list of exports containing the export.
    TYPE(ExportType), POINTER :: export !<On return, a pointer to the export of the specified user number. In no export with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: exportIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Export_UserNumberFind",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(exports)) CALL FlagError("Exports is not associated.",err,error,*999)
    IF(ASSOCIATED(export)) CALL FlagError("Export is already associated.",err,error,*999)
#endif    
    
    !Get the export from the user number
    NULLIFY(export)
    IF(ALLOCATED(exports%exports)) THEN
      DO exportIdx=1,exports%numberOfExports
        IF(ASSOCIATED(exports%exports(exportIdx)%ptr)) THEN
          IF(exports%exports(exportIdx)%ptr%userNumber==userNumber) THEN
            export=>exports%exports(exportIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The export pointer in exports is not associated for export index "// &
            & TRIM(NumberToVString(exportIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)          
        ENDIF
      ENDDO !exportIdx      
    ENDIF
    
    EXITS("Export_UserNumberFind")
    RETURN
999 ERRORSEXITS("Export_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE Export_UserNumberFind

  !
  !================================================================================================================================
  !

  !>Gets the VTK export information from an export
  SUBROUTINE Export_VTKExportGet(export,vtkExport,err,error,*)

    !Argument variables
    TYPE(ExportType), POINTER, INTENT(IN) :: export !<The export to get the VTK export for.
    TYPE(VTKExportType), POINTER, INTENT(INOUT) :: vtkExport !<On exit, a pointer to the VTK export for the export. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Export_VTKExportGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(vtkExport)) CALL FlagError("VTK export is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(export)) CALL FlagError("Export is not associated.",err,error,*999)
#endif    

    !Get the export VTK export
    vtkExport=>export%vtkExport

#ifdef WITH_POSTCHECKS    
    !Check vtk export is associated.
    IF(.NOT.ASSOCIATED(vtkExport)) THEN
      localError="Export VTK export is not associated for export number "// &
        & TRIM(NumberToVString(export%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Export_VTKExportGet")
    RETURN
999 NULLIFY(vtkExport)    
998 ERRORSEXITS("Export_VTKExportGet",err,error)
    RETURN 1
    
  END SUBROUTINE Export_VTKExportGet

  !
  !================================================================================================================================
  !

  !>Gets the field variable from an export variable
  SUBROUTINE ExportVariable_FieldVariableGet(exportVariable,fieldVariable,err,error,*)

    !Argument variables
    TYPE(ExportVariableType), POINTER, INTENT(IN) :: exportVariable !<The export variable to get the field variable for.
    TYPE(FieldVariableType), POINTER, INTENT(INOUT) :: fieldVariable !<On exit, a pointer to the field variable for the export variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("ExportVariable_FieldVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(exportVariable)) CALL FlagError("Export variable is not associated.",err,error,*999)
#endif    

    !Get the export variable field variable
    fieldVariable=>exportVariable%fieldVariable

#ifdef WITH_POSTCHECKS    
    !Check field variable is associated.
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      localError="Export variable field variable is not associated for export variable index "// &
        & TRIM(NumberToVString(exportVariable%exportIndex,"*",err,error))
      IF(ASSOCIATED(exportVariable%export)) localError=localError//" of export number "// &
        & TRIM(NumberToVString(exportVariable%export%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("ExportVariable_FieldVariableGet")
    RETURN
999 NULLIFY(fieldVariable)    
998 ERRORSEXITS("ExportVariable_FieldVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE ExportVariable_FieldVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the parameter set type from an export variable
  SUBROUTINE ExportVariable_ParameterSetTypeGet(exportVariable,parameterSetType,err,error,*)

    !Argument variables
    TYPE(ExportVariableType), POINTER, INTENT(IN) :: exportVariable !<The export variable to get the parameter set type for.
    INTEGER(INTG), INTENT(OUT) :: parameterSetType !<On exit, the parameter set type for the export variable.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ExportVariable_ParameterSetTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(exportVariable)) CALL FlagError("Export variable is not associated.",err,error,*999)
#endif    

    !Get the export variable parameter set type
    parameterSetType=exportVariable%parameterSetType
    
    EXITS("ExportVariable_ParameterSetTypeGet")
    RETURN
999 ERRORSEXITS("ExportVariable_ParameterSetTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE ExportVariable_ParameterSetTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the domain information from a VTK export
  SUBROUTINE VTKExport_DomainGet(vtkExport,domain,err,error,*)

    !Argument variables
    TYPE(VTKExportType), POINTER, INTENT(IN) :: vtkExport !<The VTK export to get the domain for.
    TYPE(DomainType), POINTER, INTENT(INOUT) :: domain !<On exit, a pointer to the domain for the VTK export. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("VTKExport_DomainGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(domain)) CALL FlagError("Domain is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vtkExport)) CALL FlagError("VTK export is not associated.",err,error,*999)
#endif    

    !Get the VTK export domain
    domain=>vtkExport%domain

#ifdef WITH_POSTCHECKS    
    !Check vtk domain is associated.
    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("VTK export domain is not associated.",err,error,*999)
#endif    
    
    EXITS("VTKExport_DomainGet")
    RETURN
999 NULLIFY(domain)    
998 ERRORSEXITS("VTKExport_DomainGet",err,error)
    RETURN 1
    
  END SUBROUTINE VTKExport_DomainGet

  !
  !================================================================================================================================
  !

  !>Gets the export information from a VTK export
  SUBROUTINE VTKExport_ExportGet(vtkExport,export,err,error,*)

    !Argument variables
    TYPE(VTKExportType), POINTER, INTENT(IN) :: vtkExport !<The VTK export to get the export for.
    TYPE(ExportType), POINTER, INTENT(INOUT) :: export !<On exit, a pointer to the export for the VTK export. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("VTKExport_ExportGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(export)) CALL FlagError("Export is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vtkExport)) CALL FlagError("VTK export is not associated.",err,error,*999)
#endif    

    !Get the VTK export export
    export=>vtkExport%export

#ifdef WITH_POSTCHECKS    
    !Check vtk export is associated.
    IF(.NOT.ASSOCIATED(export)) CALL FlagError("VTK export export is not associated.",err,error,*999)
#endif    
    
    EXITS("VTKExport_ExportGet")
    RETURN
999 NULLIFY(export)    
998 ERRORSEXITS("VTKExport_ExportGet",err,error)
    RETURN 1
    
  END SUBROUTINE VTKExport_ExportGet

  !
  !================================================================================================================================
  !
  
 END MODULE ExportAccessRoutines
