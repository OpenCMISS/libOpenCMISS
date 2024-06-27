#!/usr/bin/env python

#> \file
#> \author Chris Bradley
#> \brief This is an example script to solve a Laplace problem using OpenCMISS calls in python.
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
#> Contributor(s): Adam Reeve
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

#> \example ClassicalField/Laplace/LaplacePy/LaplaceExample.py
## Example script to solve a Laplace problem using OpenCMISS calls in python.
## \par Latest Builds:
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/Laplace/LaplacePy/build-intel'>Linux Intel Build</a>
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/Laplace/LaplacePy/build-gnu'>Linux GNU Build</a>
#<


# Add Python bindings directory to PATH
import sys, os

# Intialise OpenCMISS
from opencmiss.opencmiss import opencmiss

# Set problem parameters
height = 1.0
width = 2.0
length = 3.0

(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1,12)

numberGlobalXElements = 5
numberGlobalYElements = 5
numberGlobalZElements = 5

opencmiss.DiagnosticsSetOn(opencmiss.DiagnosticTypes.IN,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

# Get the computational nodes information
numberOfComputationalNodes = opencmiss.ComputationalNumberOfNodesGet()
computationalNodeNumber = opencmiss.ComputationalNodeNumberGet()

# Creation a RC coordinate system
coordinateSystem = opencmiss.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = opencmiss.Region()
region.CreateStart(regionUserNumber,opencmiss.WorldRegion)
region.label = "LaplaceRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-linear lagrange basis
basis = opencmiss.Basis()
basis.CreateStart(basisUserNumber)
basis.type = opencmiss.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 3
basis.interpolationXi = [opencmiss.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
basis.quadratureNumberOfGaussXi = [2]*3
basis.CreateFinish()

# Create a generated mesh
generatedMesh = opencmiss.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = opencmiss.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [width,height,length]
generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements,numberGlobalZElements]

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
geometricField.decomposition = decomposition
geometricField.ComponentMeshComponentSet(opencmiss.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(opencmiss.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(opencmiss.FieldVariableTypes.U,3,1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create standard Laplace equations set
equationsSetField = opencmiss.Field()
equationsSet = opencmiss.EquationsSet()
equationsSetSpecification = [opencmiss.EquationsSetClasses.CLASSICAL_FIELD,
        opencmiss.EquationsSetTypes.LAPLACE_EQUATION,
        opencmiss.EquationsSetSubtypes.STANDARD_LAPLACE]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = opencmiss.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.DOFOrderTypeSet(opencmiss.FieldVariableTypes.U,opencmiss.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(opencmiss.FieldVariableTypes.DELUDELN,opencmiss.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(opencmiss.FieldVariableTypes.U,opencmiss.FieldParameterSetTypes.VALUES,1,0.5)

# Create equations
equations = opencmiss.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = opencmiss.EquationsSparsityTypes.SPARSE
equations.outputType = opencmiss.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Create Laplace problem
problem = opencmiss.Problem()
problemSpecification = [opencmiss.ProblemClasses.CLASSICAL_FIELD,
        opencmiss.ProblemTypes.LAPLACE_EQUATION,
        opencmiss.ProblemSubtypes.STANDARD_LAPLACE]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = opencmiss.Solver()
problem.SolversCreateStart()
problem.SolverGet([opencmiss.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = opencmiss.SolverOutputTypes.SOLVER
solver.linearType = opencmiss.LinearSolverTypes.ITERATIVE
solver.linearIterativeAbsoluteTolerance = 1.0E-12
solver.linearIterativeRelativeTolerance = 1.0E-12
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

# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = opencmiss.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
firstNodeNumber=1
nodes = opencmiss.Nodes()
region.NodesGet(nodes)
lastNodeNumber = nodes.numberOfNodes
firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber,1)
lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber,1)
if firstNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,firstNodeNumber,1,opencmiss.BoundaryConditionsTypes.FIXED,0.0)
if lastNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,opencmiss.FieldVariableTypes.U,1,1,lastNodeNumber,1,opencmiss.BoundaryConditionsTypes.FIXED,1.0)
solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Export results
baseName = "laplace"
dataFormat = "PLAIN_TEXT"
fml = opencmiss.FieldMLIO()
fml.OutputCreate(mesh, "", baseName, dataFormat)
fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField,
    opencmiss.FieldVariableTypes.U, opencmiss.FieldParameterSetTypes.VALUES)
fml.OutputAddFieldNoType(baseName+".phi", dataFormat, dependentField,
    opencmiss.FieldVariableTypes.U, opencmiss.FieldParameterSetTypes.VALUES)
fml.OutputWrite("LaplaceExample.xml")
fml.Finalise()

opencmiss.Finalise()
