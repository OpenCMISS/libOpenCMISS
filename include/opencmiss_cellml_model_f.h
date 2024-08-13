#ifndef _CELLMLMODELF_H_
#define _CELLMLMODELF_H_

typedef void * OpaqueCellMLObject;

#ifdef __cplusplus
extern "C"
{
#else
  typedef struct CellMLModel CellMLModel;
#endif
  
  /**
   * Create a CellMLModel object from the provided URI.
   * @param uri NULL-terminated string containing the URI of the model definition.
   * @param model A newly created CellMLModel object.
   * @return zero if no error occured, otherwise non-zero to indicate an error occured.
   */
  int CellMLModel_CreateF(
			  const char* uri,
			  OpaqueCellMLObject* model
			  );
  
  /**
   * Destroy an allocated CellMLModel object.
   * @param model An existing CellML model definition object
   * @return zero if no error occured, otherwise non-zero to indicate an error occured.
   */
  int CellMLModel_DestroyF(
			   OpaqueCellMLObject* model
			   );
  
  /**
   * Get the initial value (if specified) of the named variable.
   * @param model An existing CellMLModel object.
   * @param variableName NULL-terminated string containing the name of the model variable to get the initial value of.
   * This string should be in the format of 'componentName/variableName'.
   * @param value On successful exit, will be the initial value of the named variable; otherwise uninitialised.
   * @return zero if no error occured; otherwise non-zero to indicate an error occured and value is not set.
   */  
  int CellMLModel_GetInitialValueF(
				   OpaqueCellMLObject model,
				   const char* variableName,
				   double* value
				   );
  
  /**
   * Get the initial value (if specified) of the variableIndex'th variable of the given type.
   * Will search the local components of this model for a variable flagged with the specified type and assigned the
   * specified index. If that variable (or its corresponding source variable) has an initial value, return that value.
   * \todo Really need to look at using the evaluation type constant to do this, as it is possible to assign initial
   * values and parameters using 'x = value' type equations.
   * @param model An existing CellMLModel object.
   * @param variableType The type of the variable to look for.
   * @param variableIndex The index of the variable to look for (this is the index of the variable in the array that
   * corresponds to the given type.
   * @param value On successful exit, will be the initial value of the found variable; otherwise uninitialised.
   * @return zero if no error occured; otherwise non-zero to indicate an error occured and value is not set.
   */  
  int CellMLModel_GetInitialValueByIndexF(
					  OpaqueCellMLObject model,
					  const int variableType,
					  const int variableIndex,
					  double* value
					  );
  
  
  /**
   * Get the number of algebraic variables for the given valid model.
   * @param model The CellMLModel to use.
   * @param numOfAlgebraics On return, the number of algebraic variables in the CellMLModel object.
   * @return zero if success; non-zero otherwise.
   */
  int CellMLModel_GetNumberOfAlgebraicsF(
					 OpaqueCellMLObject model,
					 int* numOfAlgebraics
					 );
  
  /**
   * Get the number of computed constant variables for the given valid model.
   * @param model The CellMLModel to use.
   * @param numOfComputedConstants On return, the number of computed constant variables in the CellMLModel object.
   * @return zero if success; non-zero otherwise.
   */
  int CellMLModel_GetNumberOfComputedConstantsF(
						OpaqueCellMLObject model,
						int* numOfComputedConstants
						);
  
  /**
   * Get the number of constant variables for the given valid model.
   * @param model The CellML model definition to use.
   * @param numberOfConstants On return, the number of constants in the CellMLModel object.
   * @return zero if success; non-zero otherwise.
   */
  int CellMLModel_GetNumberOfConstantsF(
					OpaqueCellMLObject model,
					int* numOfConstants
					);
  
  /**
   * Get the required size of the rates array for the given valid model.
   * @param model The CellMLModel to use.
   * @param numOfRates On return, the number of rates in the CellMLModel object.
   * @return zero if success; non-zero otherwise.
   */
  int CellMLModel_GetNumberOfRatesF(
				    OpaqueCellMLObject model,
				    int* numOfRates
				    );

  /**
   * Get the number of state variables for the given valid model.
   * @param model The CellMLModel to use.
   * @param numOfStates On return, the number of state variables in the CellMLModel object.
   * @return zero if success; non-zero otherwise.
   */
  int CellMLModel_GetNumberOfStatesF(
				     OpaqueCellMLObject model,
				     int* numOfStates
				     );
  
  /**
   * Get the current type of the specified variable.
   * @param model An existing CellMLModel object.
   * @param name NULL-terminated string containing the name of the model variable to get the type of. This string
   * should be in the format of 'component_name/variable_name'.
   * @param variableType On successful exit, will be the type of the named variable; otherwise uninitialised.
   * State=1, known=2, wanted=3, independent=4?
   * @return zero if no error occured; otherwise non-zero to indicate an error occured and variable_type is not set.
   */
  int CellMLModel_GetVariableTypeF(
				   OpaqueCellMLObject model,
				   const char* name,
				   int* variableType
				   );

  /**
   * Get the current index of the specified variable. A C-style index will be returned (i.e., starting from 0).
   * @param model An existing CellMLModel object.
   * @param name NULL-terminated string containing the name of the model variable to get the index of. This string
   * should be in the format of 'component_name/variable_name'.
   * @param variable_index On successful exit, will be the index of the named variable; otherwise uninitialised.
   * @return zero if no error occured; otherwise non-zero to indicate an error occured and variable_index is not set.
   */
  int CellMLModel_GetVariableIndexF(
				    OpaqueCellMLObject model,
				    const char* name,
				    int* variableIndex
				    );

  /**
   * Flag the specified variable as being 'known' for the purposes of code generation. This implies
   * that the variable will have its value set externally to the CellML model.
   * @param model An existing CellMLModel object.
   * @param name NULL-terminated string containing the name of the model variable to flag. This string
   * should be in the format of 'componentName/variableName'.
   * @return zero if no error, non-zero otherwise.
   */
  int CellMLModel_SetVariableAsKnownF(
				      OpaqueCellMLObject model,
				      const char* variableName
				      );

  /**
   * Flag the specified variable as being 'wanted' for the purposes of code generation. This implies
   * that the variable will have its value used externally to the CellML model.
   * @param model An existing CellMLModel object.
   * @param variableName NULL-terminated string containing the name of the model variable to flag. This string
   * should be in the format of 'componentName/variableName'.
   * @return zero if no error, non-zero otherwise.
   */
  int CellMLModel_SetVariableAsWantedF(
				       OpaqueCellMLObject model,
				       const char* name
				       );

  /**
   * Instantiate the CellML model definition into simulat-able code.
   * @return zero if success; non-zero otherwise.
   */
  int CellMLModel_InstantiateF(
			       OpaqueCellMLObject model
			       );

  /**
   * Call the model's compute computed constants method.
   * This method evaluates all required computed constant variables in the model.
   * @param model The CellML model definition to use.
   * @param variables the array to use for computed constant variables -., those model outputs previously
   * flagged as wanted to indicate that their value will be used externally.
   */
  void CellMLModel_CallComputedConstantsRoutineF(
						 OpaqueCellMLObject model,
						 double* variables
						 );
  
  /**
   * Call the model's compute rates method.
   * This method evaluates all required rates in the model.
   * @param model The CellML model definition to use.
   * @param voi The current value of the variable of integration (usually time).
   * @param states The array to use for state variables.
   * @param rates The array to use for rates.
   * @param wanted The array to use for wanted variables - i.e., those model outputs previously
   * flagged as wanted to indicate that their value will be used externally.
   * @param known The array to use for known variables - i.e., those model parameters previously
   * flagged as known to indicate that their value will be set externally.
   */
  void CellMLModel_CallRatesRoutineF(
				     OpaqueCellMLObject model,
				     double voi,
				     double* states,
				     double* rates,
				     double* wanted,
				     double* known
				     );
  
  /**
   * Call the model' compute variables method.
   * This method evaluates all required variables the model.
   * @param model The CellML model definition to use.
   * @param voi The current value of the variable of integration (usually time).
   * @param states The array to use for state variables.
   * @param rates The array to use for rates.
   * @param variables The array to use for variables.
   */
  void CellMLModel_CallVariablesRoutineF(
					 OpaqueCellMLObject model,
					 double voi,
					 double* states,
					 double* rates,
					 double* variables
					 );
  
 /**
   * Return the last CellML model error code and string.
   * @return zero if success; non-zero otherwise.
   */
  int CellMLModel_GetLastErrorF(
				OpaqueCellMLObject model,
				const int maxErrorStringSize,
				char* lastError
			       );
  
  /**
   * Convert a URI to an absolute URI.
   * @param URI The URI to get the absolute URI from.
   * @param asoluteURI On return, the absolute URI.
   * @return zero if success; non-zero otherwise.
   */
  int getAbsoluteURI(
		     const char* URI,
		     char** absoluteURI
		     );
  

#ifdef __cplusplus
} /* extern C */
#endif

#endif // _CELLMLMODELF_H_

