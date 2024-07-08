add_executable(CellMLModelIntegration CellMLModelIntegration.F90)
set_target_properties(CellMLModelIntegration PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(CellMLModelIntegration ${OpenCMISS_SHARED_Fortran_TARGET})
oc_add_test(CellML_ModelIntegration CellMLModelIntegration ${CMAKE_CURRENT_SOURCE_DIR}/n98.xml)

add_executable(Monodomain Monodomain.F90)
set_target_properties(Monodomain PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(Monodomain ${OpenCMISS_SHARED_Fortran_TARGET})
oc_add_test(CellML_Monodomain Monodomain 0.005 0.1001 70 ${CMAKE_CURRENT_SOURCE_DIR}/n98.xml)
