!> \file
!> \author Chris Bradley
!> \brief This module handles all Helmholtz equations routines.
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

!>This module handles all Helmholtz equations routines.
MODULE HelmholtzEquationsRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BoundaryConditionsRoutines
  USE Constants
  USE ControlLoopRoutines
  USE DistributedMatrixVector
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsMappingRoutines
  USE EquationsMatricesRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE NodeRoutines
  USE ProblemAccessRoutines
  USE Strings
  USE SolverRoutines
  USE Timer
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE,HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP, &
    & HELMHOLTZ_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET,HelmholtzEquation_EquationsSetSpecificationSet, &
    & HELMHOLTZ_EQUATION_PROBLEM_SUBTYPE_SET,HELMHOLTZ_EQUATION_PROBLEM_SETUP
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Helmholtz equation finite element equations set.
  SUBROUTINE HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EquationsType), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_NO_SOURCE_HELMHOLTZ_SUBTYPE)


          
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a Helmholtz equation type of a classical field equations set class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets up the Helmholtz equation type of a classical field equations set class.
  SUBROUTINE HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Helmholtz equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_NO_SOURCE_HELMHOLTZ_SUBTYPE)
        CALL HELMHOLTZ_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)        
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Helmholtz equation type of a classical field equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Helmholtz equation type of an classical field equations set class.
  SUBROUTINE HELMHOLTZ_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HELMHOLTZ_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_NO_SUBTYPE)        
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Helmholtz equation type of an classical field equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HELMHOLTZ_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 ERRORSEXITS("HELMHOLTZ_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HELMHOLTZ_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Helmholtz equation type of a classical field equations set class.
  SUBROUTINE HelmholtzEquation_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("HelmholtzEquation_EquationsSetSpecificationSet",err,error,*999)

    IF(SIZE(specification,1)<3) THEN
      CALL FlagError("Equations set specification must have at least three entries for a Helmholtz type equations set.", &
        & err,error,*999)
    ENDIF
    
    subtype=specification(3)
    
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_NO_SOURCE_HELMHOLTZ_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
        & " is not valid for a Helmholtz type of a classical field equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Set full specification
    CALL EquationsSet_SpecificationSet(equationsSet,3,[EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
      & EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE,subtype],err,error,*999)
    
    EXITS("HelmholtzEquation_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("HelmholtzEquation_EquationsSetSpecificationSet",err,error)
    EXITS("HelmholtzEquation_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE HelmholtzEquation_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the linear Helmholtz equation.
  SUBROUTINE HELMHOLTZ_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HELMHOLTZ_EQUATION_EQUATION_SET_LINEAR_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_HELMHOLTZ_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%setupType)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",ERR,ERROR))// &
            & " is invalid for a standard Helmholtz equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a linear Helmholtz equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HELMHOLTZ_EQUATION_EQUATIONS_SET_LINEAR_SETUP")
    RETURN
999 ERRORSEXITS("HELMHOLTZ_EQUATION_EQUATIONS_SET_LINEAR_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HELMHOLTZ_EQUATION_EQUATIONS_SET_LINEAR_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the Helmholtz solution.
  SUBROUTINE HELMHOLTZ_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)
    
    !Argument variables
    TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the solutions set to setup a Helmholtz equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HELMHOLTZ_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_NO_SOURCE_HELMHOLTZ_SUBTYPE)
        CALL HELMHOLTZ_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NumberToVString(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Helmholtz equation type of a classical field problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HELMHOLTZ_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("HELMHOLTZ_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HELMHOLTZ_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a Helmholtz equation type .
  SUBROUTINE HELMHOLTZ_EQUATION_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HELMHOLTZ_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_NO_SOURCE_HELMHOLTZ_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_HELMHOLTZ_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NO_SOURCE_HELMHOLTZ_SUBTYPE     
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NumberToVString(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Helmholtz equation type of a classical field problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("HELMHOLTZ_EQUATION_PROBLEM_SUBTYPE_SET")
    RETURN
999 ERRORSEXITS("HELMHOLTZ_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HELMHOLTZ_EQUATION_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the linear Helmholtz equations solution.
  SUBROUTINE HELMHOLTZ_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HELMHOLTZ_EQUATION_PROBLEM_LINEAR_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%SUBTYPE==PROBLEM_NO_SOURCE_HELMHOLTZ_SUBTYPE) THEN
        SELECT CASE(PROBLEM_SETUP%setupType)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",ERR,ERROR))// &
            & " is invalid for a linear Helmholtz equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a linear Helmholtz equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("HELMHOLTZ_EQUATION_PROBLEM_LINEAR_SETUP")
    RETURN
999 ERRORSEXITS("HELMHOLTZ_EQUATION_PROBLEM_LINEAR_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HELMHOLTZ_EQUATION_PROBLEM_LINEAR_SETUP

  !
  !================================================================================================================================
  !
 
END MODULE HelmholtzEquationsRoutines
