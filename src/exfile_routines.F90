!> \author Chris Bradley
!> \brief This module handles all VTK IO related routines.
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

!> This module handles all Exfile IO related routines.
MODULE ExfileRoutines
  
  USE BaseRoutines
  USE Constants
  USE ExportAccessRoutines
  USE FieldAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  
  
  IMPLICIT NONE
  
  PRIVATE

  !Module parameters
  
  !Module types 

  !Module variables

  !Interfaces

  PUBLIC Exfile_Export

  PUBLIC ExfileExport_CheckVariable

  PUBLIC ExfileExport_CreateFinish
    
CONTAINS

  !
  !================================================================================================================================
  !

  !>Exports a Exfile file
  SUBROUTINE Exfile_Export(exfileExport,err,error,*)

    !Argument variables
    TYPE(ExfileExportType), POINTER, INTENT(IN) :: exfileExport !<The exfile export information object
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ExportType), POINTER :: export
        
    ENTERS("Exfile_Export",err,error,*999)

    IF(.NOT.ASSOCIATED(exfileExport)) CALL FlagError("Exfile export is not associated.",err,error,*999)

    NULLIFY(export)
    CALL ExfileExport_ExportGet(exfileExport,export,err,error,*999)
       
    EXITS("Exfile_Export")
    RETURN
999 ERRORSEXITS("Exfile_Export",err,error)
    RETURN 1

  END SUBROUTINE Exfile_Export
 
  !
  !================================================================================================================================
  !

 !>Checks that a field variable can be handled by a exfile export
  SUBROUTINE ExfileExport_CheckVariable(exfileExport,fieldVariable,startComponent,endComponent,err,error,*)

    !Argument variables
    TYPE(ExfileExportType), POINTER, INTENT(IN) :: exfileExport !<The exfile export information object to check the variable for
    TYPE(FieldVariableType), POINTER, INTENT(IN) :: fieldVariable !<The field variable to check the export for
    INTEGER(INTG), INTENT(IN) :: startComponent !<The start component of the field variable to check
    INTEGER(INTG), INTENT(IN) :: endComponent !<The end component of the field variable to check
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
       
    ENTERS("ExfileExport_CheckVariable",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(exfileExport)) CALL FlagError("Exfile export is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
#endif    

      
    EXITS("ExfileExport_CheckVariable")
    RETURN
999 ERRORSEXITS("ExfileExport_CheckVariable",err,error)
    RETURN 1

  END SUBROUTINE ExfileExport_CheckVariable
  
  !
  !================================================================================================================================
  !

  !>Finishes the creation of a Exfile export
  SUBROUTINE ExfileExport_CreateFinish(exfileExport,err,error,*)

    !Argument variables
    TYPE(ExfileExportType), POINTER :: exfileExport !<The exfile export to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ExfileExport_CreateFinish",err,error,*999)

 
    EXITS("ExfileExport_CreateFinish")
    RETURN
999 ERRORSEXITS("ExfileExport_CreateFinish",err,error)
    RETURN 1

  END SUBROUTINE ExfileExport_CreateFinish

  !
  !================================================================================================================================
  !

END MODULE ExfileRoutines
