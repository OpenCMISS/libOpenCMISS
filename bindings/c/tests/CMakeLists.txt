function(oc_add_test TESTNAME)
    add_test(NAME ${TESTNAME} COMMAND ${ARGN})
    add_test_path(${TESTNAME} "${OpenCMISS_BINARY_DIR}" "${OpenCMISS_BINDINGS_BIN_DIR}")
endfunction()

add_executable(ComplexMesh complex_mesh.c)
set_target_properties(ComplexMesh PROPERTIES LINKER_LANGUAGE C)
target_link_libraries(ComplexMesh ${OpenCMISS_SHARED_C_TARGET})
target_include_directories(ComplexMesh PRIVATE "${OpenCMISS_INC_DIR}" "${OpenCMISS_GEN_BINDINGS_C_DIR}" "${CURRENT_CMAKE_BINARY_DIR}")
oc_add_test(C_Bindings_ComplexMesh ComplexMesh)

add_executable(Laplace_C laplace.c)
set_target_properties(Laplace_C PROPERTIES LINKER_LANGUAGE C)
target_link_libraries(Laplace_C ${OpenCMISS_SHARED_C_TARGET})
target_include_directories(Laplace_C PRIVATE "${OpenCMISS_INC_DIR}" "${OpenCMISS_GEN_BINDINGS_C_DIR}" "${CURRENT_CMAKE_BINARY_DIR}")
oc_add_test(C_Bindings_Laplace Laplace_C)
    
