#ifndef _CELLMLMODEL_H_
#define _CELLMLMODEL_H_

#include <filesystem>
#include <map>
#include <vector>

#include <libcellml>

#define CELLML_MODEL_NO_ERROR 0
#define CELLML_MODEL_COULD_NOT_CREATE_MODEL_ERROR -9995
#define CELLML_MODEL_INVALID_URI_ERROR -9996
#define CELLML_MODEL_PTR_NOT_NULL_ERROR -9997
#define CELLML_MODEL_PTR_IS_NULL_ERROR -9998
#define CELLML_MODEL_ERROR -9999

/**
 * The primary object used to define a CellML model for use in OpenCMISS.
 *
 * This is the interface object sitting between a CellML description
 * of a mathematical model and the use of that model in OpenCMISS.
 */
class CellMLModel
{
public:
  /**
   * Default constructor
   */
  CellMLModel();
  /**
   * TODO: Need to do a copy constructor
   */
  CellMLModel(
	      const CellMLModel& src
	      )
  {
  }

  /**
   * Construct a model definition from a given source document at the givien URL.
   * @param url The URL of the source document form which to create the model defintion.
   */
  CellMLModel(
	      const char* url
	      );
  /**
   * Destructor.
   */
  ~CellMLModel();
  
  /**
   * Get the initial value (if specified) of the named variable.
   * @param name The name of the model variable to get the initial value of. This string
   * should be in the format of 'component_name/variable_name'.
   * @param value On successful exit, will be the initial value of the named variable; otherwise uninitialised.
   * @return zero if no error occured; otherwise non-zero to indicate and error occured and value is not set.
   */
  int getInitialValue(
		      const char* name,
		      double* value
		      );

  /**
   * Get the initial value (if specified) of the index'th variable of the given type.
   * Will search the local components of this model for a variable flagged with the specified type and assigned the
   * specified index. If that variable (or its corresponding source variable) has an initial value, return that value.
   * \todo Really need to look at using the evaluation type constant to do this, as it is possible to assign initial
   * values and parameters using 'x = value' type equations.
   * @param type The type of the variable to look for.
   * @param index The index of the variable to look for (this is the index of the variable in the array that corresponds
   * to the given type.
   * @param value On successful exit, will be the initial value of the found variable; otherwise uninitialised.
   * @return zero if no error occured; otherwise non-zero to indicate and error occured and value is not set.
   */
  int getInitialValueByIndex(
			     const int type,
			     const int index,
			     double* value
			     );

  /**
   * Get the number of algebraics for a valid model.
   * @param numberOfAlgebraics On exit, the number of algebraics for a valid model.
   * @return the last error code for the model.
   */
  int getNumberOfAlgebraics(
			    int* numberOfAlgebraics
			    );

  /**
   * Get the number of computed constants for a valid model.
   * @param numberOfComputedConstants On exit, the number of computed constants for a valid model.
   * @return the last error code for the model.
   */
  int getNumberOfComputedConstants(
				   int* numberOfComputedConstants
				   );

  /**
   * Get the number of constants for a valid model.
   * @param numberOfConstants On exit, the number of constants for a valid model.
   * @return the last error code for the model.
   */
  int getNumberOfConstants(
			   int* numberOfConstants
			   );

  /**
   * Get the number of states for a valid model.
   * @param numberOfStates On exit, the number of states for a valid model.
   * @return the last error code for the model.
   */
  int getNumberOfStates(
			int* numberOfStates
			);

  /**
   * Get the last error code and message.
   * @param maxErrorStringLength The maximum size of the last error string.
   * @param lastError On exit, the last error message for the mode
   * @return the last error code for the model.
   */
  int getLastError(
		   const int maxErrorStringLength,
		   char* lastError
		   );

  /**
   * Get the current type of the specified variable.
   * @param name The name of the model variable to get the type of. This string should be in the
   * format of 'component_name/variable_name'.
   * @param type On successful exit, will be the type of the named variable; otherwise uninitialised. State=1, known=2, wanted=3, independent=4?
   * @return zero if no error occured; otherwise non-zero to indicate and error occured and variable_type is not set.
   */
  int getVariableType(
		      const char* name,
		      int* type
		      );

  /**
   * Get the current index of the specified variable. C-style index will be returned (starting from 0).
   * @param name The name of the model variable to get the index of. This string
   * should be in the format of 'component_name/variable_name'.
   * @param index On successful exit, will be the C index of the named variable; otherwise uninitialised.
   * @return zero if no error occured; otherwise non-zero to indicate and error occured and variable_index is not set.
   */
  int getVariableIndex(
		       const char* name,
		       int* index
		       );

  /**
   * Flag the specified variable as being 'known' for the purposes of code generation. This implies
   * that the variable will have its value set externally to the CellML model.
   * @param name The name of the model variable to flag. This string should be in the format of 'component_name/variable_name'.
   * @return 0 if no error, non-zero otherwise.
   */
  int setVariableAsKnown(
			 const char* name
			 );

  /**
   * Flag the specified variable as being 'wanted' for the purposes of code generation. This implies
   * that the variable will have its value used externally to the CellML model.
   * @param name The name of the model variable to flag. This string should be in the format of 'component_name/variable_name'.
   * @return 0 if no error, non-zero otherwise.
   */
  int setVariableAsWanted(
			  const char* name
			  );

  /**
   * Instantiate the model definition into simulat-able code.
   * @return 0 if success; non-zero otherwise.
   */
  int instantiate();

  /**
   * Set the compile command to use to compile the generated code into a dynamic shared object.
   * @param command The compile command.
   */
  void setCompileCommand(
			 const std::string& command
			 );
    
  /**
   * Get the current compile command used to compile generated code.
   * @return The current compile command.
   */
  std::string getCompileCommand();

  /**
   * Check model instantiation.
   * @return True if the model is instantiated; false otherwise.
   */
  
  bool instantiated();

  int32_t nBound;
  int32_t nRates;
  int32_t nAlgebraic;
  int32_t nConstants;

  inline void callModelComputedConstantsFunction(
						 double* variables
						 )
  {
    computedConstantsRoutine(variables);
  }

  inline void callModelRatesFunction(
				     double voi,
				     double* states,
				     double* rates,
				     double* variables,
				     double* known
				     )
  {
    ratesRoutine(voi, states, rates, known);
  }
  
  inline void callModelVariablesFunction(
					 double voi,
					 double* states,
					 double* rates,
					 double* variables
					 )
  {
    variablesRoutine(voi, states, rates, variables);
  }

private:
  
  // loaded from the generated and compiled DSO
  
  /* Compute the computed constants of the model
   */
  void (*computedConstantsRoutine)(
				   double* variables
				   );
  
  /* Compute the rates of the model
   */
  void (*ratesRoutine)(
		       double voi,
		       double* states,
		       double* rates,
		       double* variables
		       );

  /* Compute the variables of the model
   */
  void (*variablesRoutine)(
			   double voi,
			   double* states,
			   double* rates,
			   double* variables
			   );

  /**
   * Clear the last error code and string.
   */
  void clearLastError();
  
  /**
   * Dump any libCellML issues.
   * @param type The name of the libCellML issue logger.
   * @param logger The issue logger to dump.
   */
  int dumpIssues(
		 char* type,
		 libcellml::LoggerPtr logger
		 );
  
  /**
   * Set the last error code and string.
   * @param errorCode The error to set code.
   * @param errorString The error string to set.
   * @return The error code.
   */
  int setLastError(
		   const int errorCode,
		   std::string errorString
		   );

  int lastErrorCode; /* The last error code */
  std::string lastErrorString; /* The last error string */
  std::string modelURL; /* The URL of the CellML model */
  std::filesystem::path modelFile; /* The file path for the CellML model file */
  libcellml::ParserPtr parser; /* The pointer to the libCellML parser */
  int numberOfParserErrors;
  libcellml::ModelPtr model; /* The pointer to the libCellML model */
  libcellml::ImporterPtr importer; /* The pointer to the libCellML importer */
  int numberOfImporterErrors;  
  libcellml::ValidatorPtr validator; /* The pointer to the libCellL validator */
  int numberOfValidatorErrors;
  libcellml::AnalyserPtr analyser; /* The pointer to the libCellML analyser */
  int numberOfAnalyserErrors;
  bool isValid; /* True if the model is valid and has no errors */
  libcellml::AnalyserModelPtr analyserModel;
  libcellml::GeneratorPtr generator; /* The pointer to the libCellML generator */
  libcellml::GeneratorProfilePtr generatorProfile; /* The pointer to the libCellML generator profile */
  std::string interfaceCodeString;
  std::string implementationCodeString;
  std::string tmpDirectoryName;
  bool tmpDirectoryExists;
  std::string interfaceCodeFileName;
  bool interfaceCodeFileExists;
  std::string implementationCodeFileName;
  bool implementationCodeFileExists;
  std::string dsoFileName;
  bool dsoFileExists;
  std::string compileCommand;
  void* dsoHandle;
  bool isInstantiated;
  size_t numberOfStateVariables;
  size_t numberOfCellMLVariables;
  size_t totalNumberOfVariables;
  int numberOfStates;
  int numberOfConstants;
  int numberOfComputedConstants;
  int numberOfAlgebraics;
  int numberOfExternals;
  int numberOfParameters;
  int numberOfIntermediates;
  std::map<std::string, int> variableNameToIndexMap;
  std::map<int, std::string> variableIndexToNameMap;
  std::map<int, int> variableIndexToCellMLTypeMap;
  std::map<int, int> variableIndexToCellMLIndexMap;
  std::map<std::pair<int,int>, int> cellMLToVariableIndexMap;
  std::map<int, int> variableIndexToOpenCMISSTypeMap;
  std::map<int, int> variableIndexToOpenCMISSIndexMap;
  std::map<std::pair<int,int>, int> openCMISSToVariableIndexMap;
  std::map<int, double> variableIndexToInitialValueMap;
    
  bool mSaveTempFiles;
  std::map<std::pair<int,int>, double> mInitialValues;
  std::map<std::string, int> mVariableTypes;
  std::map<std::string, int> mVariableIndices;

public:

  void* mModel;
  int mNumberOfWantedVariables;
  int mNumberOfKnownVariables;
  int mNumberOfIndependentVariables;
  int mStateCounter;
  int mIntermediateCounter;
  int mParameterCounter;
};

#endif // _CELLMLMODEL_H_
