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
  oc_BasisType basis = (oc_BasisType)NULL;
  oc_BoundaryConditionsType boundaryConditions=(oc_BoundaryConditionsType)NULL;
  oc_ComputationEnvironmentType computationEnvironment = (oc_ComputationEnvironmentType)NULL;
  oc_ContextType context = (oc_ContextType)NULL;
  oc_CoordinateSystemType coordinateSystem=(oc_CoordinateSystemType)NULL;
  oc_DecompositionType decomposition=(oc_DecompositionType)NULL;
  oc_DecomposerType decomposer=(oc_DecomposerType)NULL;
  oc_EquationsType equations=(oc_EquationsType)NULL;
  oc_EquationsSetType equationsSet=(oc_EquationsSetType)NULL;
  oc_FieldType geometricField=(oc_FieldType)NULL,dependentField=(oc_FieldType)NULL,equationsSetField=(oc_FieldType)NULL;
  oc_GeneratedMeshType generatedMesh=(oc_GeneratedMeshType)NULL;
  oc_MeshType mesh=(oc_MeshType)NULL;
  oc_ProblemType problem=(oc_ProblemType)NULL;
  oc_RegionType region=(oc_RegionType)NULL,worldRegion=(oc_RegionType)NULL;
  oc_SolverType solver=(oc_SolverType)NULL;
  oc_SolverEquationsType solverEquations=(oc_SolverEquationsType)NULL;
  oc_WorkGroupType worldWorkGroup=(oc_WorkGroupType)NULL;

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

  err = oc_Initialise();
  OPENCMISS_CHECK_ERROR(err,"Initialising OpenCMISS");
  err = oc_Context_Initialise(&context);
  OPENCMISS_CHECK_ERROR(err,"Initialising context");
  err = oc_Context_Create(CONTEXT_USER_NUMBER,context);
  OPENCMISS_CHECK_ERROR(err,"Creating context");
  err = oc_Region_Initialise(&worldRegion);
  OPENCMISS_CHECK_ERROR(err,"Initialising world region");
  err = oc_Context_WorldRegionGet(context,worldRegion);
  OPENCMISS_CHECK_ERROR(err,"Get world region");
  err = oc_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR);

  err = oc_ComputationEnvironment_Initialise(&computationEnvironment);
  err = oc_Context_ComputationEnvironmentGet(context,computationEnvironment);

  err = oc_WorkGroup_Initialise(&worldWorkGroup);
  err = oc_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup);
  err = oc_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,&numberOfComputationNodes);
  err = oc_WorkGroup_GroupNodeNumberGet(worldWorkGroup,&computationNodeNumber);

  /* Start the creation of a new RC coordinate system */
  err = oc_CoordinateSystem_Initialise(&coordinateSystem);
  err = oc_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem);
  if(NUMBER_GLOBAL_Z_ELEMENTS==0)
    {
      /* Set the coordinate system to be 2D */
      err = oc_CoordinateSystem_DimensionSet(coordinateSystem,2);
    }
  else
    {
      /* Set the coordinate system to be 3D */
      err = oc_CoordinateSystem_DimensionSet(coordinateSystem,3);
    }
  /* Finish the creation of the coordinate system */
  err = oc_CoordinateSystem_CreateFinish(coordinateSystem);

  /* Start the creation of the region */
  err = oc_Region_Initialise(&region);
  err = oc_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region);
  /* Set the regions coordinate system to the 2D RC coordinate system that we have created */
  err = oc_Region_CoordinateSystemSet(region,coordinateSystem);
  /* Finish the creation of the region */
  err = oc_Region_CreateFinish(region);

  /* Start the creation of a basis (default is trilinear lagrange) */
  err = oc_Basis_Initialise(&basis);
  err = oc_Basis_CreateStart(BASIS_USER_NUMBER,context,basis);
  if(NUMBER_GLOBAL_Z_ELEMENTS==0)
    {
      /* Set the basis to be a bilinear Lagrange basis */
      err = oc_Basis_NumberOfXiSet(basis,2);
    }
  else
    {
      /* Set the basis to be a trilinear Lagrange basis */
      err = oc_Basis_NumberOfXiSet(basis,3);
    }
  /* Finish the creation of the basis */
  err = oc_Basis_CreateFinish(basis);

  /* Start the creation of a generated mesh in the region */
  err = oc_GeneratedMesh_Initialise(&generatedMesh);
  err = oc_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh);
  /* Set up a regular x*y*z mesh */
  err = oc_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE);
  /* Set the default basis */
  err = oc_GeneratedMesh_BasisSet(generatedMesh,1,&basis);
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
  err = oc_GeneratedMesh_ExtentSet(generatedMesh,MAX_COORDINATES,meshExtent);
  err = oc_GeneratedMesh_NumberOfElementsSet(generatedMesh,MAX_COORDINATES,numberXiElements);
  /* Finish the creation of a generated mesh in the region */
  err = oc_Mesh_Initialise(&mesh);
  /* Finish the creation of a generated mesh in the region */
  err = oc_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh);

  /* Create a decomposition */
  err = oc_Decomposition_Initialise(&decomposition);
  err = oc_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition);
  /* Finish the decomposition */
  err = oc_Decomposition_CreateFinish(decomposition);

  /* Create a decomposer */
  err = oc_Decomposer_Initialise(&decomposer);
  err = oc_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer);
  /* Add in the decomposition */
  err = oc_Decomposer_DecompositionAdd(decomposer,decomposition,&decompositionIndex);
  /* Finish the decomposer */
  err = oc_Decomposer_CreateFinish(decomposer);
  
  /* Start to create a default (geometric) field on the region */
  err = oc_Field_Initialise(&geometricField);
  err = oc_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField);
  /* Set the decomposition to use */
  err = oc_Field_DecompositionSet(geometricField,decomposition);
  /* Set the domain to be used by the field components. */
  err = oc_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,1,1);
  err = oc_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,2,1);
  if(NUMBER_GLOBAL_Z_ELEMENTS!=0)
    {
      err = oc_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,3,1);
    }
  /* Finish creating the field */
  err = oc_Field_CreateFinish(geometricField);

  /* Update the geometric field parameters */
  err = oc_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField);

  /* Create the equations_set */
  err = oc_EquationsSet_Initialise(&equationsSet);
  err = oc_Field_Initialise(&equationsSetField);
  equationsSetSpecification[0] = OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS;
  equationsSetSpecification[1] = OC_EQUATIONS_SET_LAPLACE_EQUATION_TYPE;
  equationsSetSpecification[2] = OC_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE;
  err = oc_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField, \
    3,equationsSetSpecification,EQUATIONS_SET_FIELD_USER_NUMBER, \
    equationsSetField,equationsSet);
  OPENCMISS_CHECK_ERROR(err,"Creating equations set");
  /* Finish creating the equations set */
  err = oc_EquationsSet_CreateFinish(equationsSet);

  /* Create the equations set dependent field variables */
  err = oc_Field_Initialise(&dependentField);
  err = oc_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField);
  /* Finish the equations set dependent field variables */
  err = oc_EquationsSet_DependentCreateFinish(equationsSet);

  /* Create the equations set equations */
  err = oc_Equations_Initialise(&equations);
  err = oc_EquationsSet_EquationsCreateStart(equationsSet,equations);
  /* Set the equations matrices sparsity type */
  err = oc_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES);
  /* Set the equations set output */
  /* err = oc_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT); */
  err = oc_Equations_OutputTypeSet(equations,OC_EQUATIONS_TIMING_OUTPUT);
  /* err = oc_Equations_OutputTypeSet(equations,OC_EQUATIONS_MATRIX_OUTPUT); */
  /* err = oc_Equations_OutputTypeSet(equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT); */
  /* Finish the equations set equations */
  err = oc_EquationsSet_EquationsCreateFinish(equationsSet);

  /* Start the creation of a problem, setting the problem to be a standard Laplace problem. */
  err = oc_Problem_Initialise(&problem);
  problemSpecification[0] = OC_PROBLEM_CLASSICAL_FIELD_CLASS;
  problemSpecification[1] = OC_PROBLEM_LAPLACE_EQUATION_TYPE;
  problemSpecification[2] = OC_PROBLEM_STANDARD_LAPLACE_SUBTYPE;
  err = oc_Problem_CreateStart(PROBLEM_USER_NUMBER,context,3,problemSpecification,problem);
  /* Finish the creation of a problem. */
  err = oc_Problem_CreateFinish(problem);

  /* Start the creation of the problem control loop */
  err = oc_Problem_ControlLoopCreateStart(problem);
  /* Finish creating the problem control loop */
  err = oc_Problem_ControlLoopCreateFinish(problem);

  /* Start the creation of the problem solvers */
  err = oc_Solver_Initialise(&solver);
  err = oc_Problem_SolversCreateStart(problem);
  err = oc_Problem_SolverGet(problem,1,controlLoopIdentifier,1,solver);
  /* err = oc_Solver_OutputTypeSet(solver,OC_SOLVER_NO_OUTPUT); */
  /* err = oc_Solver_OutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT); */
  /* err = oc_Solver_OutputTypeSet(solver,OC_SOLVER_TIMING_OUTPUT); */
  /* err = oc_Solver_OutputTypeSet(solver,OC_SOLVER_SOLVER_OUTPUT); */
  err = oc_Solver_OutputTypeSet(solver,OC_SOLVER_MATRIX_OUTPUT);
  err = oc_Solver_LinearTypeSet(solver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE);
  err = oc_Solver_LibraryTypeSet(solver,OC_SOLVER_MUMPS_LIBRARY);
  /* Finish the creation of the problem solver */
  err = oc_Problem_SolversCreateFinish(problem);

  /* Start the creation of the problem solver equations */
  solver=(oc_SolverType)NULL;
  err = oc_Solver_Initialise(&solver);
  err = oc_SolverEquations_Initialise(&solverEquations);
  err = oc_Problem_SolverEquationsCreateStart(problem);
  /* Get the solve equations */
  err = oc_Problem_SolverGet(problem,1,controlLoopIdentifier,1,solver);
  err = oc_Solver_SolverEquationsGet(solver,solverEquations);
  /* Set the solver equations sparsity */
  err = oc_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES);
  /* err = oc_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_FULL_MATRICES);  */
  /* Add in the equations set */
  err = oc_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,&equationsSetIndex);
  /* Finish the creation of the problem solver equations */
  err = oc_Problem_SolverEquationsCreateFinish(problem);

  /* Start the creation of the equations set boundary conditions */
  err = oc_BoundaryConditions_Initialise(&boundaryConditions);
  err = oc_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions);
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
  err = oc_Decomposition_NodeDomainGet(decomposition,firstNodeNumber,1,&firstNodeDomain);
  err = oc_Decomposition_NodeDomainGet(decomposition,lastNodeNumber,1,&lastNodeDomain);
  if(firstNodeDomain==computationNodeNumber)
    {
      err = oc_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,firstNodeNumber,1, \
        OC_BOUNDARY_CONDITION_FIXED,0.0);
    }
  if(lastNodeDomain==computationNodeNumber)
    {
      err = oc_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,lastNodeNumber,1, \
        OC_BOUNDARY_CONDITION_FIXED,1.0);
    }
  /* Finish the creation of the equations set boundary conditions */
  err = oc_SolverEquations_BoundaryConditionsCreateFinish(solverEquations);

  /* Solve the problem */
  err = oc_Problem_Solve(problem);

  /* Destroy the context */
  err = oc_Context_Destroy(context);
  /* Finalise OpenCMISS */
  err = oc_Finalise();

  return err;
}
