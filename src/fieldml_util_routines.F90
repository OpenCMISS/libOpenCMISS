!> \file
!> \author Caton Little
!> \brief This module handles non-IO FieldML logic.
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

!> Utility routines for FieldML

MODULE FIELDML_UTIL_ROUTINES

  USE BaseRoutines
  USE BasisRoutines
  USE CONSTANTS
  USE CoordinateSystemRoutines
  USE FieldRoutines
  USE FIELDML_API
  USE FIELDML_TYPES
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE RegionRoutines
  USE STRINGS
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Interfaces
  
  INTERFACE FIELDML_UTIL_CHECK_FIELDML_ERROR
    MODULE PROCEDURE FieldMLUtil_CheckFieldMLSessionErrorVS
    MODULE PROCEDURE FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC
  END INTERFACE FIELDML_UTIL_CHECK_FIELDML_ERROR

  PUBLIC :: FIELDML_UTIL_CHECK_FIELDML_ERROR, FIELDML_IO_INITIALISE, FIELDML_IO_FINALISE

CONTAINS

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldMLUtil_CheckFieldMLSessionErrorVS( ERROR_DESCRIPTION, FML_HANDLE, ERR, ERROR, * )
    TYPE(VARYING_STRING), INTENT(IN) :: ERROR_DESCRIPTION
    INTEGER(INTG), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    INTEGER(INTG) :: FML_ERR
    
    ENTERS( "FieldMLUtil_CheckFieldMLSessionErrorVS", ERR, ERROR, *999 )
    FML_ERR = Fieldml_GetLastError( FML_HANDLE )

    IF( FML_ERR == FML_ERR_NO_ERROR ) THEN
      EXITS( "FieldMLUtil_CheckFieldMLSessionErrorVS" )
      RETURN
    ENDIF
    
    CALL FlagError( ERROR_DESCRIPTION // " (error number " // TRIM(NumberToVString(FML_ERR,"*",ERR,ERROR)) // ")", &
      & ERR, ERROR, *999 )

999 ERRORSEXITS( "FieldMLUtil_CheckFieldMLSessionErrorVS", ERR, ERROR )
    RETURN 1
    
  END SUBROUTINE FieldMLUtil_CheckFieldMLSessionErrorVS
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC( ERROR_DESCRIPTION, FML_HANDLE, ERR, ERROR, * )
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_DESCRIPTION
    INTEGER(INTG), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    INTEGER(INTG) :: FML_ERR
    
    ENTERS( "FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC", ERR, ERROR, *999 )
    FML_ERR = Fieldml_GetLastError( FML_HANDLE )

    IF( FML_ERR == FML_ERR_NO_ERROR ) THEN
      EXITS( "FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC" )
      RETURN
    ENDIF
    
    CALL FlagError( ERROR_DESCRIPTION // " (error number " // TRIM(NumberToVString(FML_ERR,"*",ERR,ERROR)) // ")", &
      & ERR, ERROR, *999 )

999 ERRORSEXITS( "FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC", ERR, ERROR )
    RETURN 1
    
  END SUBROUTINE FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC
  
  !
  !================================================================================================================================
  !
  
  !<Initialise up the given FieldML parsing state.
  SUBROUTINE FIELDML_IO_INITIALISE( FIELDML_INFO, IS_OUT, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), POINTER :: FIELDML_INFO !<The FieldML parsing state to initialise.
    LOGICAL, INTENT(IN) :: IS_OUT !< True if the state is being used for output, false otherwise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS( "FIELDML_IO_INITIALISE", ERR, ERROR, *998 )

    IF(ASSOCIATED(FIELDML_INFO)) THEN
      CALL FlagError("FieldML info is already associated",ERR,ERROR,*998)
    ENDIF

    ALLOCATE(FIELDML_INFO, STAT = ERR)
    IF(ERR /= 0) THEN
      CALL FlagError("Could not allocate FieldML info",ERR,ERROR,*998)
    ENDIF

    FIELDML_INFO%IS_OUT = IS_OUT
    FIELDML_INFO%FML_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%NODES_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%MESH_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%ELEMENTS_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%XI_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%NODE_DOFS_HANDLE = FML_INVALID_HANDLE
    !fieldmlInfo%elementDofsHandle = FML_INVALID_HANDLE
    !fieldmlInfo%constantDofsHandle = FML_INVALID_HANDLE
    NULLIFY( FIELDML_INFO%COMPONENT_HANDLES )
    NULLIFY( FIELDML_INFO%BASIS_HANDLES )
    NULLIFY( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES )
    NULLIFY( FIELDML_INFO%BASIS_LAYOUT_HANDLES )
    
    CALL LIST_CREATE_START( FIELDML_INFO%COMPONENT_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%COMPONENT_HANDLES, LIST_INTG_TYPE, ERR, ERROR, *999 )
    CALL LIST_MUTABLE_SET( FIELDML_INFO%COMPONENT_HANDLES, .TRUE., ERR, ERROR, *999 )
    CALL LIST_CREATE_FINISH( FIELDML_INFO%COMPONENT_HANDLES, ERR, ERROR, *999 )
    
    CALL LIST_CREATE_START( FIELDML_INFO%BASIS_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%BASIS_HANDLES, LIST_INTG_TYPE, ERR, ERROR, *999 )
    CALL LIST_CREATE_FINISH( FIELDML_INFO%BASIS_HANDLES, ERR, ERROR, *999 )
    
    CALL LIST_CREATE_START( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, LIST_INTG_TYPE, ERR, ERROR, *999 )
    CALL LIST_CREATE_FINISH( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, ERR, ERROR, *999 )
    
    CALL LIST_CREATE_START( FIELDML_INFO%BASIS_LAYOUT_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%BASIS_LAYOUT_HANDLES, LIST_INTG_TYPE, ERR, ERROR, *999 )
    CALL LIST_CREATE_FINISH( FIELDML_INFO%BASIS_LAYOUT_HANDLES, ERR, ERROR, *999 )

    EXITS( "FIELDML_IO_INITIALISE" )
    RETURN

999 CALL FIELDML_IO_FINALISE(FIELDML_INFO, DUMMY_ERR, DUMMY_ERROR, *998)
998 ERRORSEXITS( "FIELDML_IO_INITIALISE", ERR, ERROR )
    RETURN 1
    
  END SUBROUTINE FIELDML_IO_INITIALISE

  !
  !================================================================================================================================
  !

  !<Clean up the given FieldML parsing state.
  SUBROUTINE FIELDML_IO_FINALISE( FIELDML_INFO, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), POINTER :: FIELDML_INFO !<The FieldML parsing state to clean up.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    !Locals
    INTEGER(INTG) :: FML_ERR

    ENTERS( "FIELDML_IO_FINALISE", ERR, ERROR, *999 )

    IF(ASSOCIATED(FIELDML_INFO)) THEN
      FML_ERR = Fieldml_Destroy( FIELDML_INFO%FML_HANDLE )
      CALL LIST_DESTROY( FIELDML_INFO%COMPONENT_HANDLES, ERR, ERROR, *999 )
      CALL LIST_DESTROY( FIELDML_INFO%BASIS_HANDLES, ERR, ERROR, *999 )
      CALL LIST_DESTROY( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, ERR, ERROR, *999 )
      CALL LIST_DESTROY( FIELDML_INFO%BASIS_LAYOUT_HANDLES, ERR, ERROR, *999 )
      DEALLOCATE(FIELDML_INFO)
    ENDIF

    EXITS( "FIELDML_IO_FINALISE" )
    RETURN
999 ERRORSEXITS( "FIELDML_IO_FINALISE", ERR, ERROR )
    RETURN 1
    
  END SUBROUTINE FIELDML_IO_FINALISE

  !
  !================================================================================================================================
  !

END MODULE FIELDML_UTIL_ROUTINES
