#ifndef _CELLMLMODELF_H_
#define _CELLMLMODELF_H_

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * Create the CellMLModel from the provided URI.
   * @param uri NULL-terminated string containing the URI of the model definition.
   * @return A newly created cellml_model_definition_f object.
   */
  void* cellml_model_create_f(const char* uri);
  
  /**
   * Destroy an allocated CellMLModel.
   * @param model An existing CellML model definition object
   */
  void cellml_model_destroy_f(void** model);
  
  /**
   * Get the initial value (if specified) of the named variable.
   * @param model An existing CellML model definition object.
   * @param name NULL-terminated string containing the name of the model variable to get the initial value of. This string
   * should be in the format of 'component_name/variable_name'.
   * @param value On successful exit, will be the initial value of the named variable; otherwise uninitialised.
   * @return zero if no error occured; otherwise non-zero to indicate and error occured and value is not set.
   */  
  int cellml_model_get_initial_value_f(void* model,
				       const char* name,
				       double* value);
  
  /**
   * Get the initial value (if specified) of the index'th variable of the given type.
   * Will search the local components of this model for a variable flagged with the specified type and assigned the
   * specified index. If that variable (or its corresponding source variable) has an initial value, return that value.
   * \todo Really need to look at using the evaluation type constant to do this, as it is possible to assign initial
   * values and parameters using 'x = value' type equations.
   * @param model An existing CellML model definition object.
   * @param type The type of the variable to look for.
   * @param index The index of the variable to look for (this is the index of the variable in the array that corresponds
   * to the given type.
   * @param value On successful exit, will be the initial value of the found variable; otherwise uninitialised.
   * @return zero if no error occured; otherwise non-zero to indicate and error occured and value is not set.
   */  
  int cellml_model_get_initial_value_by_index_f(void* model,
						const int* const type,
						const int* const index,
						double* value);
  
  /**
   * Get the current type of the specified variable.
   * @param model An existing CellML model definition object.
   * @param name NULL-terminated string containing the name of the model variable to get the type of. This string
   * should be in the format of 'component_name/variable_name'.
   * @param variable_type On successful exit, will be the type of the named variable; otherwise uninitialised. State=1, known=2, wanted=3, independent=4?
   * @return zero if no error occured; otherwise non-zero to indicate and error occured and variable_type is not set.
   */
  int cellml_model_get_variable_type_f(void* model,
				       const char* name,
				       int* variable_type);

  /**
   * Get the current index of the specified variable. A C-style index will be returned (i.e., starting from 0).
   * @param model An existing CellML model definition object.
   * @param name NULL-terminated string containing the name of the model variable to get the index of. This string
   * should be in the format of 'component_name/variable_name'.
   * @param variable_index On successful exit, will be the index of the named variable; otherwise uninitialised.
   * @return zero if no error occured; otherwise non-zero to indicate and error occured and variable_index is not set.
   */
  int cellml_model_get_variable_index_f(void* model,
					const char* name,
					int* variable_index);

  /**
   * Flag the specified variable as being 'known' for the purposes of code generation. This implies
   * that the variable will have its value set externally to the CellML model.
   * @param model An existing CellML model definition object.
   * @param name NULL-terminated string containing the name of the model variable to flag. This string
   * should be in the format of 'component_name/variable_name'.
   * @return 0 if no error, non-zero otherwise.
   */
  int cellml_model_set_variable_as_known_f(void* model,
					   const char* name);

  /**
   * Flag the specified variable as being 'wanted' for the purposes of code generation. This implies
   * that the variable will have its value used externally to the CellML model.
   * @param model An existing CellML model definition object.
   * @param name NULL-terminated string containing the name of the model variable to flag. This string
   * should be in the format of 'component_name/variable_name'.
   * @return 0 if no error, non-zero otherwise.
   */
  int cellml_model_set_variable_as_wanted_f(void* model,
					    const char* name);

  /**
   * Set the save temporary files flag.
   * @param model The CellML model definition to use.
   * @param state Should be set to non-zero to save generated temporary files.
   */
  void cellml_model_set_save_temp_files_f(void* model,
					  int state);

  /**
   * Get the current state of the save temporary files flag.
   * @param model The CellML model definition to use.
   * @return 1 if the current state is true (save), zero if it is false, and -1 in case of error.
   */
  int cellml_model_get_save_temp_files_f(void* model);

  /**
   * Instantiate the CellML model definition into simulat-able code.
   * @return 0 if success; non-zero otherwise.
   */
  int cellml_model_instantiate_f(void* model);

  /**
   * Get the required size of the constants array for the given instantiated model.
   * @param model The CellML model definition to use.
   * @return The required size of the constants array.
   */
  int cellml_model_get_number_constants_f(void* model);
  
  /**
   * Get the required size of the rates array for the given instantiated model.
   * @param model The CellML model definition to use.
   * @return The required size of the rates array.
   */
  int cellml_model_get_number_rates_f(void* model);
  
  /**
   * Get the required size of the algebraic array for the given instantiated model.
   * @param model The CellML model definition to use.
   * @return The required size of the algebraic array.
   */
  int cellml_model_get_number_algebraic_f(void* model);

  /**
   * Call the model's right-hand-side method.
   * This method evaluates all required equations in this model's system of ODEs.
   * @param model The CellML model definition to use.
   * @param voi The current value of the variable of integration (usually time).
   * @param states The array to use for state variables.
   * @param rates The array to use for rates.
   * @param wanted The array to use for wanted variables - i.e., those model outputs previously flagged as wanted to indicate that their value will be used externally.
   * @param known The array to use for known variables - i.e., those model parameters previously flagged as known to indicate that their value will be set externally.
   */
  void cellml_model_call_rhs_routine_f(void* model,
				       double voi,
				       double* states,
				       double* rates,
				       double* wanted,
				       double* known);
  
#ifdef __cplusplus
} /* extern C */
#endif

#endif // _CELLMLMODELF_H_

