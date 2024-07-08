/*
 * \file
 * \author Chris Bradley
 * \brief This is an example program to solve Laplace's equation using OpenCMISS calls from C.
 *
 * \section LICENSE
 *
 * Version: MPL 1.1/GPL 2.0/LGPL 2.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * The Original Code is OpenCMISS
 *
 * The Initial Developer of the Original Code is University of Auckland,
 * Auckland, New Zealand and University of Oxford, Oxford, United
 * Kingdom. Portions created by the University of Auckland and University
 * of Oxford are Copyright (C) 2007 by the University of Auckland and
 * the University of Oxford. All Rights Reserved.
 *
 * Contributor(s):
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU General Public License Version 2 or later (the "GPL"), or
 * the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
 * in which case the provisions of the GPL or the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of either the GPL or the LGPL, and not to allow others to
 * use your version of this file under the terms of the MPL, indicate your
 * decision by deleting the provisions above and replace them with the notice
 * and other provisions required by the GPL or the LGPL. If you do not delete
 * the provisions above, a recipient may use your version of this file under
 * the terms of any one of the MPL, the GPL or the LGPL.
 *
 */
#include <stdlib.h>
#include <stdio.h>

#include "opencmiss/opencmiss.h"

#define STRING_SIZE 20

#define HEIGHT 1.0
#define WIDTH 2.0
#define LENGTH 3.0

#define NUMBER_GLOBAL_X_ELEMENTS 5
#define NUMBER_GLOBAL_Y_ELEMENTS 5
#define NUMBER_GLOBAL_Z_ELEMENTS 5

#define CONTEXT_USER_NUMBER 1
#define COORDINATE_SYSTEM_USER_NUMBER 2
#define REGION_USER_NUMBER 3
#define BASIS_USER_NUMBER 4
#define GENERATED_MESH_USER_NUMBER 5
#define MESH_USER_NUMBER 6
#define DECOMPOSITION_USER_NUMBER 7
#define DECOMPOSER_USER_NUMBER 8
#define GEOMETRIC_FIELD_USER_NUMBER 9
#define DEPENDENT_FIELD_USER_NUMBER 10
#define EQUATIONS_SET_FIELD_USER_NUMBER 11
#define EQUATIONS_SET_USER_NUMBER 12
#define PROBLEM_USER_NUMBER 13

#define MAX_COORDINATES 3

int main()
{
  OC_BasisType basis = (OC_BasisType)NULL;
  OC_BoundaryConditionsType boundaryConditions=(OC_BoundaryConditionsType)NULL;
  OC_ComputationEnvironmentType computationEnvironment = (OC_ComputationEnvironmentType)NULL;
  OC_ContextType context = (OC_ContextType)NULL;
  OC_CoordinateSystemType coordinateSystem=(OC_CoordinateSystemType)NULL;
  OC_DecompositionType decomposition=(OC_DecompositionType)NULL;
  OC_DecomposerType decomposer=(OC_DecomposerType)NULL;
  OC_EquationsType equations=(OC_EquationsType)NULL;
  OC_EquationsSetType equationsSet=(OC_EquationsSetType)NULL;
  OC_FieldType geometricField=(OC_FieldType)NULL,dependentField=(OC_FieldType)NULL,equationsSetField=(OC_FieldType)NULL;
  OC_GeneratedMeshType generatedMesh=(OC_GeneratedMeshType)NULL;
  OC_MeshType mesh=(OC_MeshType)NULL;
  OC_ProblemType problem=(OC_ProblemType)NULL;
  OC_RegionType region=(OC_RegionType)NULL,worldRegion=(OC_RegionType)NULL;
  OC_SolverType solver=(OC_SolverType)NULL;
  OC_SolverEquationsType solverEquations=(OC_SolverEquationsType)NULL;
  OC_WorkGroupType worldWorkGroup=(OC_WorkGroupType)NULL;

  int numberOfComputationNodes,computationNodeNumber;
  int decompositionIndex,equationsSetIndex;
  int firstNodeNumber,lastNodeNumber;
  int firstNodeDomain,lastNodeDomain;

  int numberXiElements[MAX_COORDINATES];
  int controlLoopIdentifier[1];
  double meshExtent[MAX_COORDINATES];

  int equationsSetSpecification[3];
  int problemSpecification[3];

  int err;

  controlLoopIdentifier[0]=OC_CONTROL_LOOP_NODE;

  err = OC_Initialise();
  OPENCMISS_CHECK_ERROR(err,"Initialising OpenCMISS");
  err = OC_Context_Initialise(&context);
  OPENCMISS_CHECK_ERROR(err,"Initialising context");
  err = OC_Context_Create(CONTEXT_USER_NUMBER,context);
  OPENCMISS_CHECK_ERROR(err,"Creating context");
  err = OC_Region_Initialise(&worldRegion);
  OPENCMISS_CHECK_ERROR(err,"Initialising world region");
  err = OC_Context_WorldRegionGet(context,worldRegion);
  OPENCMISS_CHECK_ERROR(err,"Get world region");
  err = OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR);

  err = OC_ComputationEnvironment_Initialise(&computationEnvironment);
  err = OC_Context_ComputationEnvironmentGet(context,computationEnvironment);

  err = OC_WorkGroup_Initialise(&worldWorkGroup);
  err = OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup);
  err = OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,&numberOfComputationNodes);
  err = OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,&computationNodeNumber);

  /* Start the creation of a new RC coordinate system */
  err = OC_CoordinateSystem_Initialise(&coordinateSystem);
  err = OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem);
  if(NUMBER_GLOBAL_Z_ELEMENTS==0)
    {
      /* Set the coordinate system to be 2D */
      err = OC_CoordinateSystem_DimensionSet(coordinateSystem,2);
    }
  else
    {
      /* Set the coordinate system to be 3D */
      err = OC_CoordinateSystem_DimensionSet(coordinateSystem,3);
    }
  /* Finish the creation of the coordinate system */
  err = OC_CoordinateSystem_CreateFinish(coordinateSystem);

  /* Start the creation of the region */
  err = OC_Region_Initialise(&region);
  err = OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region);
  /* Set the regions coordinate system to the 2D RC coordinate system that we have created */
  err = OC_Region_CoordinateSystemSet(region,coordinateSystem);
  /* Finish the creation of the region */
  err = OC_Region_CreateFinish(region);

  /* Start the creation of a basis (default is trilinear lagrange) */
  err = OC_Basis_Initialise(&basis);
  err = OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis);
  if(NUMBER_GLOBAL_Z_ELEMENTS==0)
    {
      /* Set the basis to be a bilinear Lagrange basis */
      err = OC_Basis_NumberOfXiSet(basis,2);
    }
  else
    {
      /* Set the basis to be a trilinear Lagrange basis */
      err = OC_Basis_NumberOfXiSet(basis,3);
    }
  /* Finish the creation of the basis */
  err = OC_Basis_CreateFinish(basis);

  /* Start the creation of a generated mesh in the region */
  err = OC_GeneratedMesh_Initialise(&generatedMesh);
  err = OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh);
  /* Set up a regular x*y*z mesh */
  err = OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE);
  /* Set the default basis */
  err = OC_GeneratedMesh_BasisSet(generatedMesh,1,&basis);
  OPENCMISS_CHECK_ERROR(err,"Setting mesh basis");
  /* Define the mesh on the region */
  meshExtent[0]=WIDTH;
  meshExtent[1]=HEIGHT;
  numberXiElements[0]=NUMBER_GLOBAL_X_ELEMENTS;
  numberXiElements[1]=NUMBER_GLOBAL_Y_ELEMENTS;
  if(NUMBER_GLOBAL_Z_ELEMENTS!=0)
    {
      meshExtent[2]=LENGTH;
      numberXiElements[2]=NUMBER_GLOBAL_Z_ELEMENTS;
    }
  err = OC_GeneratedMesh_ExtentSet(generatedMesh,MAX_COORDINATES,meshExtent);
  err = OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,MAX_COORDINATES,numberXiElements);
  /* Finish the creation of a generated mesh in the region */
  err = OC_Mesh_Initialise(&mesh);
  /* Finish the creation of a generated mesh in the region */
  err = OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh);

  /* Create a decomposition */
  err = OC_Decomposition_Initialise(&decomposition);
  err = OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition);
  /* Finish the decomposition */
  err = OC_Decomposition_CreateFinish(decomposition);

  /* Create a decomposer */
  err = OC_Decomposer_Initialise(&decomposer);
  err = OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer);
  /* Add in the decomposition */
  err = OC_Decomposer_DecompositionAdd(decomposer,decomposition,&decompositionIndex);
  /* Finish the decomposer */
  err = OC_Decomposer_CreateFinish(decomposer);
  
  /* Start to create a default (geometric) field on the region */
  err = OC_Field_Initialise(&geometricField);
  err = OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField);
  /* Set the decomposition to use */
  err = OC_Field_DecompositionSet(geometricField,decomposition);
  /* Set the domain to be used by the field components. */
  err = OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,1,1);
  err = OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,2,1);
  if(NUMBER_GLOBAL_Z_ELEMENTS!=0)
    {
      err = OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,3,1);
    }
  /* Finish creating the field */
  err = OC_Field_CreateFinish(geometricField);

  /* Update the geometric field parameters */
  err = OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField);

  /* Create the equations_set */
  err = OC_EquationsSet_Initialise(&equationsSet);
  err = OC_Field_Initialise(&equationsSetField);
  equationsSetSpecification[0] = OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS;
  equationsSetSpecification[1] = OC_EQUATIONS_SET_LAPLACE_EQUATION_TYPE;
  equationsSetSpecification[2] = OC_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE;
  err = OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField, \
    3,equationsSetSpecification,EQUATIONS_SET_FIELD_USER_NUMBER, \
    equationsSetField,equationsSet);
  OPENCMISS_CHECK_ERROR(err,"Creating equations set");
  /* Finish creating the equations set */
  err = OC_EquationsSet_CreateFinish(equationsSet);

  /* Create the equations set dependent field variables */
  err = OC_Field_Initialise(&dependentField);
  err = OC_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField);
  /* Finish the equations set dependent field variables */
  err = OC_EquationsSet_DependentCreateFinish(equationsSet);

  /* Create the equations set equations */
  err = OC_Equations_Initialise(&equations);
  err = OC_EquationsSet_EquationsCreateStart(equationsSet,equations);
  /* Set the equations matrices sparsity type */
  err = OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES);
  /* Set the equations set output */
  /* err = OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT); */
  err = OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_TIMING_OUTPUT);
  /* err = OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_MATRIX_OUTPUT); */
  /* err = OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT); */
  /* Finish the equations set equations */
  err = OC_EquationsSet_EquationsCreateFinish(equationsSet);

  /* Start the creation of a problem, setting the problem to be a standard Laplace problem. */
  err = OC_Problem_Initialise(&problem);
  problemSpecification[0] = OC_PROBLEM_CLASSICAL_FIELD_CLASS;
  problemSpecification[1] = OC_PROBLEM_LAPLACE_EQUATION_TYPE;
  problemSpecification[2] = OC_PROBLEM_STANDARD_LAPLACE_SUBTYPE;
  err = OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,3,problemSpecification,problem);
  /* Finish the creation of a problem. */
  err = OC_Problem_CreateFinish(problem);

  /* Start the creation of the problem control loop */
  err = OC_Problem_ControlLoopCreateStart(problem);
  /* Finish creating the problem control loop */
  err = OC_Problem_ControlLoopCreateFinish(problem);

  /* Start the creation of the problem solvers */
  err = OC_Solver_Initialise(&solver);
  err = OC_Problem_SolversCreateStart(problem);
  err = OC_Problem_SolverGet(problem,1,controlLoopIdentifier,1,solver);
  /* err = OC_Solver_OutputTypeSet(solver,OC_SOLVER_NO_OUTPUT); */
  /* err = OC_Solver_OutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT); */
  /* err = OC_Solver_OutputTypeSet(solver,OC_SOLVER_TIMING_OUTPUT); */
  /* err = OC_Solver_OutputTypeSet(solver,OC_SOLVER_SOLVER_OUTPUT); */
  err = OC_Solver_OutputTypeSet(solver,OC_SOLVER_MATRIX_OUTPUT);
  err = OC_Solver_LinearTypeSet(solver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE);
  err = OC_Solver_LibraryTypeSet(solver,OC_SOLVER_MUMPS_LIBRARY);
  /* Finish the creation of the problem solver */
  err = OC_Problem_SolversCreateFinish(problem);

  /* Start the creation of the problem solver equations */
  solver=(OC_SolverType)NULL;
  err = OC_Solver_Initialise(&solver);
  err = OC_SolverEquations_Initialise(&solverEquations);
  err = OC_Problem_SolverEquationsCreateStart(problem);
  /* Get the solve equations */
  err = OC_Problem_SolverGet(problem,1,controlLoopIdentifier,1,solver);
  err = OC_Solver_SolverEquationsGet(solver,solverEquations);
  /* Set the solver equations sparsity */
  err = OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES);
  /* err = OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_FULL_MATRICES);  */
  /* Add in the equations set */
  err = OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,&equationsSetIndex);
  /* Finish the creation of the problem solver equations */
  err = OC_Problem_SolverEquationsCreateFinish(problem);

  /* Start the creation of the equations set boundary conditions */
  err = OC_BoundaryConditions_Initialise(&boundaryConditions);
  err = OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions);
  /* Set the first node to 0.0 and the last node to 1.0 */
  firstNodeNumber=1;
  if(NUMBER_GLOBAL_Z_ELEMENTS==0)
    {
      lastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1);
    }
  else
    {
      lastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1);
    }
  err = OC_Decomposition_NodeDomainGet(decomposition,firstNodeNumber,1,&firstNodeDomain);
  err = OC_Decomposition_NodeDomainGet(decomposition,lastNodeNumber,1,&lastNodeDomain);
  if(firstNodeDomain==computationNodeNumber)
    {
      err = OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,firstNodeNumber,1, \
        OC_BOUNDARY_CONDITION_FIXED,0.0);
    }
  if(lastNodeDomain==computationNodeNumber)
    {
      err = OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,lastNodeNumber,1, \
        OC_BOUNDARY_CONDITION_FIXED,1.0);
    }
  /* Finish the creation of the equations set boundary conditions */
  err = OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations);

  /* Solve the problem */
  err = OC_Problem_Solve(problem);

  /* Destroy the context */
  err = OC_Context_Destroy(context);
  /* Finalise OpenCMISS */
  err = OC_Finalise();

  return err;
}
