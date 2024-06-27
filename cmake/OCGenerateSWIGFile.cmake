file(MAKE_DIRECTORY "${LibOpenCMISS_GEN_BINDINGS_PYTHON_SWIG_DIR}")
execute_process(
  COMMAND "${Python_EXECUTABLE}" generate_bindings "${LibOpenCMISS_ROOT}" SWIG "${LibOpenCMISS_OPENCMISS_GENERATED_I}"
  RESULT_VARIABLE LibOpenCMISS_Python_RESULT_VAR
  OUTPUT_VARIABLE LibOpenCMISS_Python_OUTPUT_VAR
  ERROR_VARIABLE LibOpenCMISS_Python_ERROR_VAR
  WORKING_DIRECTORY ${LibOpenCMISS_BINDINGS_DIR}
)
if(LibOpenCMISS_Python_RESULT_VAR NOT_EQUAL 0)
  message(STATUS "Generate SWIG file failed.")
  message(STAUTS "  Result: '${LibOpenCMISS_Python_RESULT_VAR}'")
  message(STAUTS "  Output: '${LibOpenCMISS_Python_OUTPUT_VAR}'")
  message(STAUTS "  Error: '${LibOpenCMISS_Python_ERROR_VAR}'")
endif()
set_property(SOURCE "${LibOpenCMISS_OPENCMISS_GENERATED_I}" PROPERTY CPLUSPLUS ON)
list(APPEND LibOpenCMISS_CLEANUP_FILES "${LibOpenCMISS_OPENCMISS_GENERATED_I}")
