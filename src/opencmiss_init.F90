!> \file
!> \author Chris Bradley
!> \brief The top level OpenCMISS Iron module.
!>
!> \mainpage OpenCMISS Documentation
!>
!> An open source interactive computer program for Continuum Mechanics, Image analysis, Signal processing and System
!> Identification. Target usage: Bioengineering application of finite element analysis, boundary element and collocation
!> techniques.
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
!>

!> \defgroup OpenCMISS_Init OpenCMISS::OpenCMISSInit
!> The top level cmiss module.
MODULE OpenCMISSInit

  USE ISO_C_BINDING
  
  USE BaseRoutines
  USE Constants
  USE ContextRoutines
  USE ComputationAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
#ifdef WITH_MPI  
#ifdef WITH_F08_MPI
  USE MPI_F08
#elif WITH_F90_MPI
  USE MPI
#endif
  USE OpenCMISSMPI
#endif
#ifdef WITH_PETSC  
  USE OpenCMISSPETSc
#endif  
  USE Strings
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

#ifdef WITH_MPI  
#ifdef WITH_F77_MPI
#include "mpif.h"
#endif
#endif  

  !Module parameters

  !> \addtogroup OC_ErrorHandlingModes OpenCMISS::Constants::ErrorHandlingModes
  !> \brief Error handling mode parameters
  !> \see OpenCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: OC_RETURN_ERROR_CODE = 0 !<Just return the error code \see OC_ErrorHandlingModes,OpenCMISS
  INTEGER(INTG), PARAMETER :: OC_OUTPUT_ERROR = 1 !<Output the error traceback and return the error code \see OC_ErrorHandlingModes,OpenCMISS
  INTEGER(INTG), PARAMETER :: OC_TRAP_ERROR = 2 !<Trap the error by outputing the error traceback and stopping the program \see OC_ErrorHandlingModes,OpenCMISS
  !>@}
  
  !Module types

  !Module variables

  LOGICAL, SAVE :: ocFirstInit = .FALSE. !<ocFirstInit will be .TRUE. if oc has ever been initialised, .FALSE. if not.

  INTEGER(INTG), SAVE :: oc_ErrorHandlingMode !<The current error handling mode for OpenCMISS \see OC_ErrorHandlingModes

  LOGICAL, SAVE :: openCMISSMPIInitialised=.FALSE. !<Is .TRUE. if OpenCMISS has initialised MPI

  !Interfaces

  INTERFACE

    SUBROUTINE OC_InitFatalHandler() BIND(C,NAME="OC_InitFatalHandler")
    END SUBROUTINE OC_InitFatalHandler

    SUBROUTINE OC_ResetFatalHandler() BIND(C,NAME="OC_ResetFatalHandler")
    END SUBROUTINE OC_ResetFatalHandler
    
    SUBROUTINE OC_SetFatalHandler() BIND(C,NAME="OC_SetFatalHandler")
    END SUBROUTINE OC_SetFatalHandler

  END INTERFACE

  PUBLIC OC_RETURN_ERROR_CODE,OC_OUTPUT_ERROR,OC_TRAP_ERROR

  PUBLIC OC_ErrorHandlingModeGet_,OC_ErrorHandlingModeSet_
  
  PUBLIC OC_HandleError
  
  PUBLIC OC_Finalise_,OC_Initialise_

CONTAINS

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Returns the error handling mode for OpenCMISS \see OpenCMISS::OC_ErrorHandlingModeGet
  SUBROUTINE OC_ErrorHandlingModeGet_(errorHandlingMode,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: errorHandlingMode !<On return, the error handling mode. \see OC_ErrorHandlingModes,OpenCMISS
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables

    ENTERS("OC_ErrorHandlingModeGet_",err,error,*999)

    errorHandlingMode=OC_ErrorHandlingMode
    
    EXITS("OC_ErrorHandlingModeGet_")
    RETURN
999 ERRORSEXITS("OC_ErrorHandlingModeGet_",err,error)
    RETURN 1
    
  END SUBROUTINE OC_ErrorHandlingModeGet_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Sets the error handling mode for cmiss \see OpenCMISS::OC_ErrorHandlingModeSet
  SUBROUTINE OC_ErrorHandlingModeSet_(errorHandlingMode,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: errorHandlingMode !<The error handling mode to set. \see OC_ErrorHandlingModes,OpenCMISS
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("OC_ErrorHandlingModeSet_",err,error,*999)

    SELECT CASE(errorHandlingMode)
    CASE(OC_RETURN_ERROR_CODE)
      OC_ErrorHandlingMode=OC_RETURN_ERROR_CODE
    CASE(OC_OUTPUT_ERROR)
      OC_ErrorHandlingMode=OC_OUTPUT_ERROR
    CASE(OC_TRAP_ERROR)
      OC_ErrorHandlingMode=OC_TRAP_ERROR
    CASE DEFAULT
      localError="The supplied error handling mode of "//TRIM(NumberToVString(errorHandlingMode,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("OC_ErrorHandlingModeSet_")
    RETURN
999 ERRORSEXITS("OC_ErrorHandlingModeSet_",err,error)
    RETURN 1
    
  END SUBROUTINE OC_ErrorHandlingModeSet_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.
  
  !>Finalises OpenCMISS. \see OpenCMISS::OC_Finalise
  SUBROUTINE OC_Finalise_(err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables
    INTEGER(INTG) :: mpiIError
    LOGICAL :: mpiFinalised

    IF(ocFirstInit) THEN
      !Reset the signal handler
      CALL OC_ResetFatalHandler()
      !Finalise contexts
      CALL Contexts_Finalise(err,error,*999)
#ifdef WITH_PETSC      
      !Finalise PETSc
      !Call this after MPI_COMM_FREE as PETSc routines are called when some MPI comm attributes are freed.
      !CALL Petsc_LogView(PETSC_COMM_WORLD,"OpenCMISSTest.petsc",err,error,*999)
      CALL Petsc_Finalise(err,error,*999)
#endif
#ifdef WITH_MPI      
      !Finalise MPI
      IF(openCMISSMPIInitialised) THEN
        !Check if MPI has been finalised
        CALL MPI_FINALIZED(mpiFinalised,mpiIError)
        CALL MPI_ErrorCheck("MPI_FINALIZED",mpiIError,err,error,*999)
        IF(.NOT.mpiFinalised) THEN
          CALL MPI_FINALIZE(mpiIError)
          CALL MPI_ErrorCheck("MPI_FINALIZE",mpiIError,err,error,*999)
        ENDIF
      ENDIF
#endif      
      !Finalise the base routines
      CALL BaseRoutines_Finalise(err,error,*999)
      !Reset first init
      ocFirstInit=.FALSE.
    ENDIF
     
    RETURN
999 RETURN 1
    
  END SUBROUTINE OC_Finalise_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Initialises OpenCMISS. \see OpenCMISS::OC_Initialise
  SUBROUTINE OC_Initialise_(majorVersion,minorVersion,patchVersion,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: majorVersion !<OpenCMISS major version number
    INTEGER(INTG), INTENT(IN) :: minorVersion !<OpenCMISS minor version number
    INTEGER(INTG), INTENT(IN) :: patchVersion !<OpenCMISS patch version number
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mpiIError
    LOGICAL :: mpiInitialised    
    TYPE(VARYING_STRING) :: versionString

    IF(.NOT.ocFirstInit) THEN
      !Initialise error mode
      OC_ErrorHandlingMode = OC_OUTPUT_ERROR !Default for now, maybe make OC_RETURN_ERROR_CODE the default
      !Initialise the base routines
      CALL BaseRoutines_Initialise(err,error,*999)
#ifdef WITH_MPI      
      !Initialise MPI
      openCMISSMPIInitialised=.FALSE.
      CALL MPI_INITIALIZED(mpiInitialised,mpiIError)
      CALL MPI_ErrorCheck("MPI_INITIALIZED",mpiIError,err,error,*999)
      IF(.NOT.mpiInitialised) THEN
        !Initialise the MPI environment
        CALL MPI_INIT(mpiIError)
        CALL MPI_ErrorCheck("MPI_INIT",mpiIError,err,error,*999)
        openCMISSMPIInitialised=.TRUE.
      ENDIF
#endif      
#ifdef WITH_PETSC      
      !Initialise PETSc
      CALL Petsc_Initialise(PETSC_NULL_CHARACTER,err,error,*999)
#endif      
      !Initialise contexts
      CALL Contexts_Initialise(err,error,*999)
      !Setup signal handling
      CALL OC_InitFatalHandler()
      CALL OC_SetFatalHandler()

      !Write out the OpenCMISS version???
      IF(.FALSE.) THEN
        versionString="OpenCMISS version "//TRIM(NumberToVString(majorVersion,"*",err,error))
        versionString=versionString//"."
        versionString=versionString//TRIM(NumberToVString(minorVersion,"*",err,error))
        versionString=versionString//"."
        versionString=versionString//TRIM(NumberToVString(patchVersion,"*",err,error))
       !WRITE(*,'(A)') CHAR(versionString)
        versionString=""
      ENDIF

      !Set first initalised
      ocFirstInit = .TRUE.
    ENDIF
    
    RETURN
999 RETURN 1
    
  END SUBROUTINE OC_Initialise_

  !
  !================================================================================================================================
  !

  !>Handle an error condition
  SUBROUTINE OC_HandleError(err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mpiError
    
    SELECT CASE(OC_ErrorHandlingMode)
    CASE(OC_RETURN_ERROR_CODE)
      !Do nothing
    CASE(OC_OUTPUT_ERROR)
      CALL WriteError(err,error,*999)
    CASE(OC_TRAP_ERROR)
      CALL WriteError(err,error,*999)
#ifdef WITH_MPI      
      CALL MPI_ABORT(MPI_COMM_WORLD,err,mpiError)
#endif      
      STOP
    CASE DEFAULT
      !Do nothing
    END SELECT

    RETURN
999 RETURN

  END SUBROUTINE OC_HandleError
  
  !
  !================================================================================================================================
  !

END MODULE OpenCMISSInit
