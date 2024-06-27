#!/usr/bin/env python

#> \file
#> \author Chris Bradley
#> \brief This is an example script to solve a finite elasticity equation using OpenCMISS calls in python.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): 
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>

#> \example FiniteElasticity/Cantilever/src/CantileverExample.py
## Example script to solve a finite elasticity equation using OpenCMISS calls in python.
## \par Latest Builds:
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-intel'>Linux Intel Build</a>
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-gnu'>Linux GNU Build</a>
#<

#> Main script
# Add Python bindings directory to PATH
import sys, os

# Intialise OpenCMISS
from opencmiss.opencmiss import opencmiss

# Set problem parameters

width = 60.0
length = 40.0
height = 40.0
density=9.0E-4 #in g mm^-3
gravity=[0.0,0.0,-9.81] #in m s^-2

UsePressureBasis = False
NumberOfGaussXi = 2

coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
pressureBasisUserNumber = 2
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
materialFieldUserNumber = 3
dependentFieldUserNumber = 4
sourceFieldUserNumber = 5
equationsSetFieldUserNumber = 6
equationsSetUserNumber = 1
problemUserNumber = 1

# Set all diganostic levels on for testing
#opencmiss.DiagnosticsSetOn(opencmiss.DiagnosticTypes.All,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

numberOfLoadIncrements = 2
numberGlobalXElements = 1
numberGlobalYElements = 1
numberGlobalZElements = 1
InterpolationType = 1
if(numberGlobalZElements==0):
    numberOfXi = 2
else:
    numberOfXi = 3

# Get the number of computational nodes and this computational node number
numberOfComputationalNodes = opencmiss.ComputationalNumberOfNodesGet()
computationalNodeNumber = opencmiss.ComputationalNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = opencmiss.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = opencmiss.Region()
region.CreateStart(regionUserNumber,opencmiss.WorldRegion)
region.LabelSet("Region")
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Define basis
basis = opencmiss.Basis()
basis.CreateStart(basisUserNumber)
if InterpolationType in (1,2,3,4):
    basis.type = opencmiss.BasisTypes.LAGRANGE_HERMITE_TP
elif InterpolationType in (7,8,9):
    basis.type = opencmiss.BasisTypes.SIMPLEX
basis.numberOfXi = numberOfXi
basis.interpolationXi = [opencmiss.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
if(NumberOfGaussXi>0):
    basis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
basis.CreateFinish()

if(UsePressureBasis):
    # Define pressure basis
    pressureBasis = opencmiss.Basis()
    pressureBasis.CreateStart(pressureBasisUserNumber)
    if InterpolationType in (1,2,3,4):
        pressureBasis.type = opencmiss.BasisTypes.LAGRANGE_HERMITE_TP
    elif InterpolationType in (7,8,9):
        pressureBasis.type = opencmiss.BasisTypes.SIMPLEX
    pressureBasis.numberOfXi = numberOfXi
    pressureBasis.interpolationXi = [opencmiss.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
    if(NumberOfGaussXi>0):
        pressureBasis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
    pressureBasis.CreateFinish()

# Start the creation of a generated mesh in the region
generatedMesh = opencmiss.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = opencmiss.GeneratedMeshTypes.REGULAR
if(UsePressureBasis):
    generatedMesh.basis = [basis,pressureBasis]
else:
    generatedMesh.basis = [basis]
if(numberGlobalZElements==0):
    generatedMesh.extent = [width,height]
    generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements]
else:
    generatedMesh.extent = [width,length,height]
    generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements,numberGlobalZElements]
# Finish the creation of a generated mesh in the region
mesh = opencmiss.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)

# Create a decomposition for the mesh
decomposition = opencmiss.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = opencmiss.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = opencmiss.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.DecompositionSet(decomposition)
geometricField.TypeSet(opencmiss.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(opencmiss.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(opencmiss.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(opencmiss.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(opencmiss.FieldVariableTypes.U,3,1)
if InterpolationType == 4:
    geometricField.fieldScalingType = opencmiss.FieldScalingTypes.ARITHMETIC_MEAN
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create a fibre field and attach it to the geometric field
fibreField = opencmiss.Field()
fibreField.CreateStart(fibreFieldUserNumber,region)
fibreField.TypeSet(opencmiss.FieldTypes.FIBRE)
fibreField.DecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(opencmiss.FieldVariableTypes.U,"Fibre")
if InterpolationType == 4:
    fibreField.fieldScalingType = opencmiss.FieldScalingTypes.ARITHMETIC_MEAN
fibreField.CreateFinish()

# Create the equations_set
equationsSetField = opencmiss.Field()
equationsSet = opencmiss.EquationsSet()
equationsSetSpecification = [opencmiss.EquationsSetClasses.ELASTICITY,
    opencmiss.EquationsSetTypes.FINITE_ELASTICITY,
    opencmiss.EquationsSetSubtypes.MOONEY_RIVLIN]
equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
    equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create the dependent field
dependentField = opencmiss.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(opencmiss.FieldVariableTypes.U,"Dependent")
dependentField.ComponentInterpolationSet(opencmiss.FieldVariableTypes.U,4,opencmiss.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.ComponentInterpolationSet(opencmiss.FieldVariableTypes.DELUDELN,4,opencmiss.FieldInterpolationTypes.ELEMENT_BASED)
if(UsePressureBasis):
    # Set the pressure to be nodally based and use the second mesh component
    if InterpolationType == 4:
        dependentField.ComponentInterpolationSet(opencmiss.FieldVariableTypes.U,4,opencmiss.FieldInterpolationTypes.NODE_BASED)
        dependentField.ComponentInterpolationSet(opencmiss.FieldVariableTypes.DELUDELN,4,opencmiss.FieldInterpolationTypes.NODE_BASED)
    dependentField.ComponentMeshComponentSet(opencmiss.FieldVariableTypes.U,4,2)
    dependentField.ComponentMeshComponentSet(opencmiss.FieldVariableTypes.DELUDELN,4,2)
if InterpolationType == 4:
    dependentField.fieldScalingType = opencmiss.FieldScalingTypes.ARITHMETIC_MEAN
equationsSet.DependentCreateFinish()


# Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
opencmiss.Field.ParametersToFieldParametersComponentCopy(
    geometricField,opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,1,
    dependentField,opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,1)
opencmiss.Field.ParametersToFieldParametersComponentCopy(
    geometricField,opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,2,
    dependentField,opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,2)
opencmiss.Field.ParametersToFieldParametersComponentCopy(
    geometricField,opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,3,
    dependentField,opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,3)
opencmiss.Field.ComponentValuesInitialiseDP(
    dependentField,opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,4,-8.0)

# Create the material field
materialField = opencmiss.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(opencmiss.FieldVariableTypes.U,"Material")
materialField.VariableLabelSet(opencmiss.FieldVariableTypes.V,"Density")
equationsSet.MaterialsCreateFinish()

# Set Mooney-Rivlin constants c10 and c01 respectively.
materialField.ComponentValuesInitialiseDP(
    opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,1,2.0)
materialField.ComponentValuesInitialiseDP(
    opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,2,2.0)
materialField.ComponentValuesInitialiseDP(
    opencmiss.FieldVariableTypes.V,opencmiss.FieldParameterSetTypes.VALUES,1,density)

#Create the source field with the gravity vector
sourceField = opencmiss.Field()
equationsSet.SourceCreateStart(sourceFieldUserNumber,sourceField)
if InterpolationType == 4:
    sourceField.fieldScalingType = opencmiss.FieldScalingTypes.ARITHMETIC_MEAN
else:
    sourceField.fieldScalingType = opencmiss.FieldScalingTypes.UNIT
equationsSet.SourceCreateFinish()

#Set the gravity vector component values
sourceField.ComponentValuesInitialiseDP(
    opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,1,gravity[0])
sourceField.ComponentValuesInitialiseDP(
    opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,2,gravity[1])
sourceField.ComponentValuesInitialiseDP(
    opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,3,gravity[2])

# Create equations
equations = opencmiss.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = opencmiss.EquationsSparsityTypes.SPARSE
equations.outputType = opencmiss.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Define the problem
problem = opencmiss.Problem()
problemSpecification = [opencmiss.ProblemClasses.ELASTICITY,
        opencmiss.ProblemTypes.FINITE_ELASTICITY,
        opencmiss.ProblemSubtypes.NONE]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create the problem control loop
problem.ControlLoopCreateStart()
controlLoop = opencmiss.ControlLoop()
problem.ControlLoopGet([opencmiss.ControlLoopIdentifiers.NODE],controlLoop)
controlLoop.MaximumIterationsSet(numberOfLoadIncrements)
problem.ControlLoopCreateFinish()

# Create problem solver
nonLinearSolver = opencmiss.Solver()
linearSolver = opencmiss.Solver()
problem.SolversCreateStart()
problem.SolverGet([opencmiss.ControlLoopIdentifiers.NODE],1,nonLinearSolver)
nonLinearSolver.outputType = opencmiss.SolverOutputTypes.PROGRESS
nonLinearSolver.NewtonJacobianCalculationTypeSet(opencmiss.JacobianCalculationTypes.EQUATIONS)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.linearType = opencmiss.LinearSolverTypes.DIRECT
#linearSolver.libraryType = opencmiss.SolverLibraries.LAPACK
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = opencmiss.Solver()
solverEquations = opencmiss.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([opencmiss.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = opencmiss.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = opencmiss.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
# Set x=0 nodes to no x displacment
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,1,1,opencmiss.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,3,1,opencmiss.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,5,1,opencmiss.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,7,1,opencmiss.BoundaryConditionsTypes.FIXED,0.0)

# Set y=0 nodes to no y displacement
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,1,2,opencmiss.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,3,2,opencmiss.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,5,2,opencmiss.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,7,2,opencmiss.BoundaryConditionsTypes.FIXED,0.0)

# Set z=0 nodes to no y displacement
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,1,3,opencmiss.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,3,3,opencmiss.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,5,3,opencmiss.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,7,3,opencmiss.BoundaryConditionsTypes.FIXED,0.0)
solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Export results
fields = opencmiss.Fields()
fields.CreateRegion(region)
fields.NodesExport("Cantilever","FORTRAN")
fields.ElementsExport("Cantilever","FORTRAN")
fields.Finalise()

