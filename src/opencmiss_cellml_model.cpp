#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#include <unistd.h>
#include <dlfcn.h>

#include <libcellml>

#include "opencmiss_cellml_model.hpp"

/* CellML variable types */
enum CellMLVariableType {
  CellMLVariableOfIntegrationType = 1,
  CellMLStateType = 2,
  CellMLConstantType = 3,
  CellMLComputedConstantType = 4,
  CellMLAlgebraicType = 5,
  CellMLExternalType = 6
};

/* OpenCMISS variable types */
enum OpenCMISSVariableType {
  OpenCMISSStateType = 1,
  OpenCMISSParameterType = 2,
  OpenCMISSIntermediateType = 3
};

/*
 * Prototype local methods
 */

CellMLModel::CellMLModel()
{
  lastErrorCode = CELLML_MODEL_NO_ERROR;
  lastErrorString = "";
  model = NULL;
  nBound = -1;
  nRates = -1;
  nAlgebraic = -1;
  nConstants = -1;
  mNumberOfWantedVariables = 0;
  mNumberOfKnownVariables = 0;
  mNumberOfIndependentVariables = 0;
  mStateCounter = 0;
  mIntermediateCounter = 0;
  mParameterCounter = 0;
}

CellMLModel::CellMLModel(const char* url) :
  modelURL(url)
{
  lastErrorCode = CELLML_MODEL_NO_ERROR;
  lastErrorString = "";
  modelFile = "";
  parser = NULL;
  model = NULL;
  validator = NULL;
  isValid = false;
  generator = NULL;
  generatorProfile = NULL;
  compileCommand = "gcc -fPIC -g -shared -x c -o";
  
  nBound = -1;
  nRates = -1;
  nAlgebraic = -1;
  nConstants = -1;
  mNumberOfWantedVariables = 0;
  mNumberOfKnownVariables = 0;
  mNumberOfIndependentVariables = 0;
  mStateCounter = 0;
  mIntermediateCounter = 0;
  mParameterCounter = 0;
  
  // Read in the CellML file
  modelFile = url;
  if(std::filesystem::exists(modelFile))
    {
      std::ifstream inFile(modelFile);
      std::stringstream inFileContents;
      inFileContents << inFile.rdbuf();
      // NOT SURE ALL OF THESE LIBCELLML CREATIONS ETC. SHOULD BE IN THE CONSTRUCTOR?
      // Parse the CellML file. NOTE: not strict mode.
      parser = libcellml::Parser::create(false);
      //parser = libcellml::Parser::create();
      model = parser->parseModel(inFileContents.str());
      // Check for any parser issues
      int numberOfParserErrors = CellMLModel::dumpIssues((char*)"Parser",parser);
      // Resolve any imports
      importer = libcellml::Importer::create(false);
      importer->resolveImports(model,"");
      // Check for any importer issues
      int numberofImporterErrors = CellMLModel::dumpIssues((char*)"Importer",importer);
      // Flatten the model
      model = importer->flattenModel(model);
      // Validate the CellML mode
      validator = libcellml::Validator::create();
      validator->validateModel(model);
      // Check for any validation issues
      int numberOfValidatorErrors = CellMLModel::dumpIssues((char*)"Validator",validator);
      // Analyse the CellML mode
      analyser = libcellml::Analyser::create();
      analyser->analyseModel(model);
      // Check for any analyser issues
      numberOfAnalyserErrors = CellMLModel::dumpIssues((char*)"Analyser",analyser);
      isValid = ( (numberOfParserErrors == 0) ||
		  (numberOfImporterErrors == 0) ||
		  (numberOfValidatorErrors == 0) ||
		  (numberOfAnalyserErrors == 0) );
    }
  else
    {
      CellMLModel::setLastError(CELLML_MODEL_ERROR,"The CellML file: " + modelFile.string() + ", does not exist!");
    }    
}

CellMLModel::~CellMLModel()
{
  // free up libCellML objects?
  //if(parser) free parser;
  //if(model) free model;
  //if(importer) free importer;
  //if(validator) free validator;
  //if(analyser) free analyser;
  //if(generator) free generator;
  //if(generatorProfile) free generatorProfile;
  if(interfaceCodeFileExists) unlink(interfaceCodeFileName.c_str());
  if(implementationCodeFileExists) unlink(implementationCodeFileName.c_str());
  if(dsoFileExists) unlink(dsoFileName.c_str());
  if(tmpDirectoryExists) rmdir(tmpDirectoryName.c_str());
}

std::string CellMLModel::getCompileCommand()
{
  return compileCommand;
}

int CellMLModel::getInitialValue(const char* name, double* value)
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  std::string variableName = name;
  
  auto variableIndex = variableNameToIndexMap.find(variableName);
  if(variableIndex != variableNameToIndexMap.end())
    {
      int index = variableNameToIndexMap[variableName];
      auto variableValueIndex = variableIndexToInitialValueMap.find(index);
      if(variableValueIndex != variableIndexToInitialValueMap.end())
	{
	  *value = variableIndexToInitialValueMap[index];
	}
      else
	{
	  *value = 0.0;
	}
    }
  else
    {
      errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Could not find the variable name: " + variableName + "!");	  
    }
  
  return errorCode;
}

int CellMLModel::getInitialValueByIndex(const int openCMISSVariableType, const int openCMISSVariableIndex, double* value)
{
  int errorCode = CELLML_MODEL_NO_ERROR;

  auto variableIndex = openCMISSToVariableIndexMap.find(std::make_pair(openCMISSVariableType,openCMISSVariableIndex));
  if(variableIndex != openCMISSToVariableIndexMap.end())
    {
      int index = openCMISSToVariableIndexMap[std::make_pair(openCMISSVariableType,openCMISSVariableIndex)];
      auto variableValueIndex = variableIndexToInitialValueMap.find(index);
      if(variableValueIndex != variableIndexToInitialValueMap.end())
	{
	  *value = variableIndexToInitialValueMap[index];
	}
      else
	{
	  *value = 0.0;
	}
    }
  else
    {
      errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Could not find the variable of type: " + std::to_string(openCMISSVariableType) + " and index: " + std::to_string(openCMISSVariableIndex) + "!");	  
    }
    
  return errorCode;
}

int CellMLModel::getNumberOfAlgebraics(
				       int* numOfAlgebraics
				       )
{
  int errorCode = CELLML_MODEL_NO_ERROR;

  if(isInstantiated)
    {
      *numOfAlgebraics = numberOfAlgebraics;
    }
  else
    {
      errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Cannot get the number of algebraic variables. The model has not been instantiated!");
    }

  return errorCode;
}

int CellMLModel::getNumberOfComputedConstants(
					      int* numOfComputedConstants
					      )
{
  int errorCode = CELLML_MODEL_NO_ERROR;

  if(isInstantiated)
    {
      *numOfComputedConstants = numberOfComputedConstants;
    }
  else
    {
      errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Cannot get the number of computed constant variables. The model has not been instantiated!");
    }

  return errorCode;
}

int CellMLModel::getNumberOfConstants(
				      int* numOfConstants
				      )
{
  int errorCode = CELLML_MODEL_NO_ERROR;

  if(isInstantiated)
    {
      *numOfConstants = numberOfConstants;
    }
  else
    {
      errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Cannot get the number of constant variables. The model has not been instantiated!");
    }

  return errorCode;
}

int CellMLModel::getNumberOfStates(
				   int* numOfStates
				   )
{
  int errorCode = CELLML_MODEL_NO_ERROR;

  if(isInstantiated)
    {
      *numOfStates = numberOfStates;
    }
  else
    {
      errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Cannot get the number of state variables. The model has not been instantiated!");
    }

  return errorCode;
}

int CellMLModel::getLastError(
			      const int maxErrorStringLength,
			      char* errorString
			      )
{
  if(!(lastErrorCode == CELLML_MODEL_NO_ERROR)) std::strncpy(errorString,lastErrorString.c_str(),maxErrorStringLength);
  
  return lastErrorCode;
}

bool CellMLModel::instantiated()
{
  return isInstantiated;
}

int CellMLModel::getVariableType(const char* name, int* openCMISSVariableType)
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  std::string variableName = name;

  auto variableIndex = variableNameToIndexMap.find(variableName);
  if(variableIndex != variableNameToIndexMap.end())
    {
      int index = variableNameToIndexMap[variableName];
      auto variableValueIndex = variableIndexToOpenCMISSTypeMap.find(index);
      if(variableValueIndex != variableIndexToOpenCMISSTypeMap.end())
	{
	  *openCMISSVariableType = variableIndexToOpenCMISSTypeMap[index];
	}
      else
	{
	  errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Could not find the variable index: " + std::to_string(index) + " corresponding to the variable name: " + variableName + "!");	  
	}
    }
  else
    {
      errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Could not find the variable name: " + variableName + "!");	  
    }
  
  return errorCode;
}

int CellMLModel::getVariableIndex(
				  const char* name,
				  int* openCMISSVariableIndex
				  )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  std::string variableName = name;

  auto variableIndex = variableNameToIndexMap.find(variableName);
  if(variableIndex != variableNameToIndexMap.end())
    {
      int index = variableNameToIndexMap[variableName];
      auto variableValueIndex = variableIndexToOpenCMISSIndexMap.find(index);
      if(variableValueIndex != variableIndexToOpenCMISSIndexMap.end())
	{
	  *openCMISSVariableIndex = variableIndexToOpenCMISSIndexMap[index];
	}
      else
	{
	  errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Could not find the variable index: " + std::to_string(index) + " corresponding to the variable name: " + variableName + "!");	  
	}
    }
  else
    {
      errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Could not find the variable name: " + variableName + "!");	  
    }
  
  return errorCode;
}

void CellMLModel::setCompileCommand(
				    const std::string& command
				    )
{
  compileCommand = command;
}

int CellMLModel::setVariableAsKnown(const char* name)
{
  int code;
  return code;
}

int CellMLModel::setVariableAsWanted(const char* name)
{
  return 0;
}

int CellMLModel::instantiate()
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  int systemCode = 0;
  
  //libcellml::PrinterPtr printer = libcellml::Printer::create();
  //std::cout << printer->printModel(model,true);
  
  //Compute the variables
  if(isValid)
    {
      analyserModel = analyser->model();
      int variableCounter = 0;
      numberOfStateVariables = analyserModel->stateCount();
      std::cout << "Number of analyser state variables : " << numberOfStateVariables << std::endl;
      // Loop over the State variables
      size_t variableIndex;
      double initialValue;
      numberOfStates = 0;
      for(variableIndex = 0; variableIndex < numberOfStateVariables; variableIndex++)
	{
	  libcellml::AnalyserVariablePtr analyserVariable = analyserModel->state(variableIndex);
	  libcellml::VariablePtr variable = analyserVariable->variable();
	  libcellml::ComponentPtr component = std::dynamic_pointer_cast<libcellml::Component>(variable->parent());
	  std::string fullVariableName = component->name() + "/" + variable->name();
	  std::cout << "State : " << variableIndex << std::endl;
	  std::cout << "   Name : " << fullVariableName << std::endl;
	  std::cout << "   Type : " << analyserVariable->typeAsString(analyserVariable->type()) << std::endl;
	  std::cout << "   Initial value : " << variable->initialValue() << std::endl;
	  variableNameToIndexMap.insert(std::make_pair(fullVariableName,variableCounter));
	  variableIndexToNameMap.insert(std::make_pair(variableCounter,fullVariableName));
	  variableIndexToCellMLTypeMap.insert(std::make_pair(variableCounter,CellMLStateType));
	  variableIndexToCellMLIndexMap.insert(std::make_pair(variableCounter,numberOfStates));
	  cellMLToVariableIndexMap.insert(std::make_pair(std::make_pair(CellMLStateType,numberOfStates),variableCounter));
	  variableIndexToOpenCMISSTypeMap.insert(std::make_pair(variableCounter,OpenCMISSStateType));
	  variableIndexToOpenCMISSIndexMap.insert(std::make_pair(variableCounter,numberOfStates));
	  openCMISSToVariableIndexMap.insert(std::make_pair(std::make_pair(OpenCMISSStateType,numberOfStates),variableCounter));
	  initialValue = std::stod(variable->initialValue());
	  variableIndexToInitialValueMap.insert(std::make_pair(variableCounter,initialValue));
	  numberOfStates++;
	  variableCounter++;
	}	  
      numberOfCellMLVariables = analyserModel->variableCount();
      // Loop over the CellML variables
      std::cout << "Number of analyser CellML variables : " << numberOfCellMLVariables << std::endl;
      totalNumberOfVariables = numberOfStateVariables + numberOfCellMLVariables;
      numberOfConstants = 0;
      numberOfComputedConstants = 0;
      numberOfAlgebraics = 0;
      numberOfExternals = 0;
      numberOfParameters = 0;
      numberOfIntermediates = 0;
      for(variableIndex = 0; variableIndex < numberOfCellMLVariables; variableIndex++)
	{
	  libcellml::AnalyserVariablePtr analyserVariable = analyserModel->variable(variableIndex);
	  libcellml::VariablePtr variable = analyserVariable->variable();
	  libcellml::ComponentPtr component = std::dynamic_pointer_cast<libcellml::Component>(variable->parent());
	  std::string fullVariableName = component->name() + "/" + variable->name();
	  libcellml::AnalyserVariable::Type analyserVariableType = analyserVariable->type();
	  std::cout << "Variable : " << variableIndex << std::endl;
	  std::cout << "   Name : " << fullVariableName << std::endl;
	  std::cout << "   Type : " << analyserVariable->typeAsString(analyserVariable->type()) << std::endl;
	  std::cout << "   Initial value : " << variable->initialValue() << std::endl;
	  variableNameToIndexMap.insert(std::make_pair(fullVariableName,variableCounter));
	  variableIndexToNameMap.insert(std::make_pair(variableCounter,fullVariableName));
	  if(analyserVariableType == libcellml::AnalyserVariable::Type::VARIABLE_OF_INTEGRATION)
	    {
	      // Do nothing
	    }		
	  else if(analyserVariableType == libcellml::AnalyserVariable::Type::CONSTANT)
	    {
	      variableIndexToCellMLTypeMap.insert(std::make_pair(variableCounter,CellMLConstantType));
	      variableIndexToCellMLIndexMap.insert(std::make_pair(variableCounter,numberOfConstants));
	      cellMLToVariableIndexMap.insert(std::make_pair(std::make_pair(CellMLConstantType,numberOfConstants),variableCounter));
	      variableIndexToOpenCMISSTypeMap.insert(std::make_pair(variableCounter,OpenCMISSParameterType));
	      variableIndexToOpenCMISSIndexMap.insert(std::make_pair(variableCounter,numberOfParameters));
	      openCMISSToVariableIndexMap.insert(std::make_pair(std::make_pair(OpenCMISSParameterType,numberOfParameters),variableCounter));
	      initialValue = std::stod(variable->initialValue());
	      variableIndexToInitialValueMap.insert(std::make_pair(variableCounter,initialValue));
	      numberOfConstants++;
	      numberOfParameters++;
	    }		
	  else if(analyserVariableType == libcellml::AnalyserVariable::Type::COMPUTED_CONSTANT)
	    {
	      variableIndexToCellMLTypeMap.insert(std::make_pair(variableCounter,CellMLComputedConstantType));
	      variableIndexToCellMLIndexMap.insert(std::make_pair(variableCounter,numberOfComputedConstants));
	      cellMLToVariableIndexMap.insert(std::make_pair(std::make_pair(CellMLComputedConstantType,numberOfComputedConstants),variableCounter));
	      variableIndexToOpenCMISSTypeMap.insert(std::make_pair(variableCounter,OpenCMISSIntermediateType));
	      variableIndexToOpenCMISSIndexMap.insert(std::make_pair(variableCounter,numberOfIntermediates));
	      openCMISSToVariableIndexMap.insert(std::make_pair(std::make_pair(OpenCMISSIntermediateType,numberOfIntermediates),variableCounter));
	      numberOfComputedConstants++;
	      numberOfIntermediates++;
	    }
	  else if(analyserVariableType == libcellml::AnalyserVariable::Type::ALGEBRAIC)
	    {
	      variableIndexToCellMLTypeMap.insert(std::make_pair(variableCounter,CellMLAlgebraicType));
	      variableIndexToCellMLIndexMap.insert(std::make_pair(variableCounter,numberOfAlgebraics));
	      cellMLToVariableIndexMap.insert(std::make_pair(std::make_pair(CellMLAlgebraicType,numberOfAlgebraics),variableCounter));
	      variableIndexToOpenCMISSTypeMap.insert(std::make_pair(variableCounter,OpenCMISSIntermediateType));
	      variableIndexToOpenCMISSIndexMap.insert(std::make_pair(variableCounter,numberOfIntermediates));
	      openCMISSToVariableIndexMap.insert(std::make_pair(std::make_pair(OpenCMISSIntermediateType,numberOfIntermediates),variableCounter));
	      numberOfAlgebraics++;
	      numberOfIntermediates++;
	    }
	  else if(analyserVariableType == libcellml::AnalyserVariable::Type::EXTERNAL)
	    {
	      // Do nothing???
	      numberOfExternals++;
	      //numberOfIntermediate++;
	    }
	  else
	    {
	      std::cerr << "Variable " << variableIndex << " has type " << analyserVariable->typeAsString(analyserVariableType) << std::endl;
	    }
	  variableCounter++;
	}
      if(errorCode == CELLML_MODEL_NO_ERROR)
	{
	  //Generate the code
	  generator = libcellml::Generator::create();
	  generatorProfile = libcellml::GeneratorProfile::create(libcellml::GeneratorProfile::Profile::C);
	  interfaceCodeFileName = model->name() + ".h";
	  implementationCodeFileName = model->name() + ".c";
	  generatorProfile->setInterfaceFileNameString(interfaceCodeFileName);
	  generator->setProfile(generatorProfile);
	  generator->setModel(analyser->model());
	  interfaceCodeString = generator->interfaceCode();
	  implementationCodeString = generator->implementationCode();
	  
	  if((interfaceCodeString.length() > 1) && (implementationCodeString.length() > 1))
	    {
	      /* We have code, so dump it out to files in a temporary
		 directory so we can have the compiled object nice and handy to
		 delete */
	      char tmpNameTemplate[64] = "tmp.cellml2code.XXXXXX";
	      if(mkdtemp(tmpNameTemplate))
		{
		  tmpDirectoryName = tmpNameTemplate;
		  tmpDirectoryExists = true;
		  std::string codeFileName = tmpDirectoryName + "/" + interfaceCodeFileName;
		  std::ofstream outFile(codeFileName);
		  outFile << interfaceCodeString;
		  outFile.close();
		  interfaceCodeFileExists = true;
		  codeFileName = tmpDirectoryName + "/" + implementationCodeFileName;
		  outFile.open(codeFileName);
		  outFile << implementationCodeString;
		  outFile.close();
		  implementationCodeFileExists = true;
		  char* dso = (char*)malloc(codeFileName.length()+10);
		  /* need to make the dso name right so that it'll load */
		  sprintf(dso,"%s%s.so",(tmpNameTemplate[0]=='/'?"":"./"),codeFileName.c_str());
		  dsoFileName = dso;
		  free(dso);
		  /* compile the code into a shared object */
		  char* compCommand = (char*)malloc(compileCommand.length()+dsoFileName.length()+codeFileName.length()+3);
		  sprintf(compCommand,"%s %s %s",compileCommand.c_str(),dsoFileName.c_str(),codeFileName.c_str());
		  systemCode = system(compCommand);
		  if(systemCode == 0)
		    {
		      dsoFileExists = true;
		      // now load the DSO back into memory and get the required functions
		      dsoHandle = dlopen(dsoFileName.c_str(),RTLD_LOCAL|RTLD_LAZY);
		      if(dsoHandle)
			{
			  /* find the required methods */
			  if(numberOfComputedConstants > 0)
			    {
			      computedConstantsRoutine = (void (*)(double*)) dlsym(dsoHandle,"computeComputedConstants");
			      if(computedConstantsRoutine == NULL)
				{
				  errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Error getting method: computeComputedConstants");
				}
			    }
			  if(errorCode == CELLML_MODEL_NO_ERROR)
			    {
			      if(numberOfStates > 0)
				{
				  ratesRoutine = (void (*)(double,double*,double*,double*)) dlsym(dsoHandle,"computeRates");
				  if(ratesRoutine == NULL)
				    {
				      errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Error getting method: computeRates");
				    }
				}
			      if(errorCode == CELLML_MODEL_NO_ERROR)
				{
				  if(numberOfAlgebraics > 0)
				    {
				      variablesRoutine = (void (*)(double,double*,double*,double*)) dlsym(dsoHandle,"computeVariables");
				      if(variablesRoutine == NULL)
					{
					  errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Error getting method: computeVariables");
					}
				    }
				  if(errorCode == CELLML_MODEL_NO_ERROR)
				    {
				      if(((numberOfComputedConstants == 0) || computedConstantsRoutine) &&
					 ((numberOfStates == 0) || ratesRoutine) &&
					 ((numberOfAlgebraics == 0) || variablesRoutine))
					{
					  isInstantiated = true;
					}
				      else
					{
					  errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Error getting one or more of the required methods.");
					}
				    }
				}
			    }
			}
		      else
			{
			  errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Error opening shared object ("+dsoFileName+"): "+dlerror());
			}
		    }
		  else
		    {
		      errorCode = CellMLModel::setLastError(systemCode,strerror(systemCode));
		    }
		  if(compCommand) free(compCommand);
		}
	      else
		{
		  errorCode = CellMLModel::setLastError(errno,strerror(errno));
		}
	    }
	  else
	    {
	      errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"Error getting C-code for model.");
	    }
	}
    }
  else
    {
      errorCode = CellMLModel::setLastError(CELLML_MODEL_ERROR,"CellML model is not valid.");
    }
  
  return errorCode;
}

/*
 * Local methods
 */

void CellMLModel::clearLastError()
{
  lastErrorCode = CELLML_MODEL_NO_ERROR;
  lastErrorString = "";
}

int CellMLModel::dumpIssues(
			     char* type,
			     libcellml::LoggerPtr logger
			     )
{
  size_t numberOfIssues = logger->issueCount();
  int numberOfErrors = 0;
  
  if(numberOfIssues != 0)
    {
      std::cout << "The " << type << " has found " << numberOfIssues << " issues!" << std::endl;
      for (size_t e = 0; e < numberOfIssues; ++e)
	{
	  libcellml::IssuePtr issue = logger->issue(e);
	  libcellml::Issue::Level issueLevel = issue->level();
	  if(issueLevel == libcellml::Issue::Level::ERROR) numberOfErrors++;
	  std::string issueSpecificationReference = issue->referenceHeading();
	  
	  std::cout << "  " << type << " issue[" << e << "]:" << std::endl;
	  std::cout << "     Description: " << issue->description()
		    << std::endl;
	  std::cout << "     Type of item stored: " << cellmlElementTypeAsString(issue->item()->type()) << std::endl;
	  std::cout << "     URL: " << issue->url() << std::endl;
	  if (issueSpecificationReference != "")
	    {
	      std::cout << "    See section " << issueSpecificationReference
			<< " in the CellML specification." << std::endl;
	    }
	}	        
    }

  return numberOfErrors;
}

int CellMLModel::setLastError(
			       const int errorCode,
			       std::string errorString
			       )
{
  if(!(errorCode == CELLML_MODEL_NO_ERROR))
    {
      lastErrorCode = errorCode;
      lastErrorString = errorString;
    }

  return errorCode;
}
