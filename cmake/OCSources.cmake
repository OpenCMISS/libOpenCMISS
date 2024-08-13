set(OpenCMISS_C_SRC
  #binary_file_c.c
  opencmiss_init_c.c
  external_dae_solver_routines.c
  timer_c.c
)
set(OpenCMISS_C_HEADERS
  external_dae_solver_routines.h
)
set(OpenCMISS_CPP_HEADERS
  macros.h
  dllexport.h
  opencmiss_version.h
)
set(OpenCMISS_CXX_SRC )
set(OpenCMISS_CXX_HEADERS )
set(OpenCMISS_Fortran_SRC
  advection_diffusion_equation_routines.F90
  advection_equation_routines.F90
  analytic_analysis_routines.F90
  base_routines.F90
  basis_routines.F90
  basis_access_routines.F90
  #binary_file_f.F90
  biodomain_equation_routines.F90
  bioelectric_finite_elasticity_routines.F90
  bioelectric_routines.F90
  blas.F90
  boundary_condition_routines.F90
  boundary_condition_access_routines.F90
  Burgers_equation_routines.F90
  cellml_access_routines.F90
  characteristic_equation_routines.F90
  classical_field_routines.F90
  computation_routines.F90
  computation_access_routines.F90
  constants.F90
  context_routines.F90
  context_access_routines.F90
  control_loop_routines.F90
  control_loop_access_routines.F90
  coordinate_routines.F90
  coordinate_access_routines.F90
  Darcy_equations_routines.F90
  Darcy_pressure_equations_routines.F90
  data_point_routines.F90
  data_point_access_routines.F90
  data_projection_routines.F90
  data_projection_access_routines.F90
  decomposition_routines.F90
  decomposition_access_routines.F90
  diffusion_advection_diffusion_routines.F90
  diffusion_diffusion_routines.F90
  diffusion_equation_routines.F90
  distributed_matrix_vector_IO.F90
  distributed_matrix_vector.F90
  distributed_matrix_vector_access.F90
  domain_mappings.F90
  elasticity_routines.F90
  electromechanics_routines.F90
  electrophysiology_cell_routines.F90
  equations_routines.F90
  equations_access_routines.F90
  equations_mapping_routines.F90
  equations_mapping_access_routines.F90
  equations_matrices_routines.F90
  equations_matrices_access_routines.F90
  equations_set_routines.F90
  equations_set_access_routines.F90
  field_IO_routines.F90
  field_routines.F90
  field_access_routines.F90
  finite_elasticity_Darcy_routines.F90
  finite_elasticity_fluid_pressure_routines.F90
  finite_elasticity_routines.F90
  finite_elasticity_utility_routines.F90
  fitting_routines.F90
  fluid_mechanics_IO_routines.F90
  fluid_mechanics_routines.F90
  fsi_routines.F90
  generated_mesh_routines.F90
  generated_mesh_access_routines.F90
  Hamilton_Jacobi_equations_routines.F90
  Helmholtz_equations_routines.F90
  #Helmholtz_TEMPLATE_equations_routines.F90
  history_routines.F90
  input_output.F90
  interface_conditions_routines.F90
  interface_condition_access_routines.F90
  interface_equations_routines.F90
  interface_equations_access_routines.F90
  interface_mapping_routines.F90
  interface_mapping_access_routines.F90
  interface_matrices_access_routines.F90
  interface_matrices_routines.F90
  interface_operators_routines.F90
  interface_routines.F90
  interface_access_routines.F90
  iso_varying_string.F90
  kinds.F90
  lapack.F90
  Laplace_equations_routines.F90
  linear_elasticity_routines.F90
  linkedlist_routines.F90
  lists.F90
  maths.F90
  matrix_vector.F90
  matrix_vector_access.F90
  mesh_routines.F90
  mesh_access_routines.F90
  monodomain_equations_routines.F90
  multi_compartment_transport_routines.F90
  multi_physics_routines.F90
  Navier_Stokes_equations_routines.F90
  node_routines.F90
  node_access_routines.F90
  opencmiss.F90
  opencmiss_cellml.F90
  opencmiss_fortran_c.F90
  opencmiss_init.F90
  opencmiss_mpi.F90
  opencmiss_parmetis.F90
  opencmiss_petsc.F90
  opencmiss_petsc_types.F90
  Poiseuille_equations_routines.F90
  Poisson_equations_routines.F90
  problem_routines.F90
  problem_access_routines.F90
  profiling_routines.F90
  reaction_diffusion_equation_routines.F90
  reaction_diffusion_IO_routines.F90
  region_routines.F90
  region_access_routines.F90
  solver_mapping_routines.F90
  solver_mapping_access_routines.F90
  solver_matrices_routines.F90
  solver_matrices_access_routines.F90
  solver_routines.F90
  solver_access_routines.F90
  sorting.F90
  Stokes_equations_routines.F90
  stree_equation_routines.F90
  strings.F90
  test_framework_routines.F90
  timer_f.F90
  trees.F90
  types.F90
  util_array.F90
)
# Add platform dependent files
if(UNIX OR LINUX} MATCHES linux)
  list(APPEND OpenCMISS_Fortran_SRC machine_constants_linux.F90)
elseif(APPLE)
  list(APPEND OpenCMISS_Fortran_SRC machine_constants_linux.F90)
elseif(WIN32)
  list(APPEND OpenCMISS_Fortran_SRC machine_constants_win32.F90)
else()
  message(WARNING "The operating system is not implemented.")
endif()

# Add in CellML files
if(OpenCMISS_WITH_CellML)
  list(APPEND OpenCMISS_C_HEADERS
    opencmiss_cellml_model_f.h
  )
  list(APPEND OpenCMISS_CXX_HEADERS
    opencmiss_cellml_model.hpp
  )
  list(APPEND OpenCMISS_CXX_SRC
    opencmiss_cellml_model.cpp
    opencmiss_cellml_model_f.cpp
  )
  list(APPEND OpenCMISS_Fortran_SRC
    opencmiss_cellml_model.F90
   )
endif()

# Add in FieldML files
if(OpenCMISS_WITH_FieldML)
  list(APPEND OpenCMISS_C_HEADERS
    FieldExport.h
    FieldExportConstants.h
  )
  list(APPEND OpenCMISS_C_SRC
    FieldExport.c
  )
  list(APPEND OpenCMISS_Fortran_SRC
    fieldml_input_routines.F90
    fieldml_output_routines.F90
    fieldml_types.F90
    fieldml_util_routines.F90
  )
endif()

# Fix paths to files
set(FIXPATH_VARS OpenCMISS_C_SRC OpenCMISS_CXX_SRC OpenCMISS_Fortran_SRC)
foreach(varname ${FIXPATH_VARS})
  set(_tmp )
  foreach(filename ${${varname}})
    list(APPEND _tmp ${OpenCMISS_SRC_DIR}/${filename}) 
  endforeach()
  set(${varname} ${_tmp})
endforeach()
set(FIXPATH_VARS OpenCMISS_C_HEADERS OpenCMISS_CXX_HEADERS OpenCMISS_CPP_HEADERS)
foreach(varname ${FIXPATH_VARS})
  set(_tmp )
  foreach(filename ${${varname}})
    list(APPEND _tmp ${OpenCMISS_INC_DIR}/${filename}) 
  endforeach()
  set(${varname} ${_tmp})
endforeach()

# Set combined sources variable
set(OpenCMISS_HEADERS ${OpenCMISS_C_HEADERS} ${OpenCMISS_CXX_HEADERS} ${OpenCMISS_CPP_HEADERS})
set(OpenCMISS_Fortran_SOURCES ${OpenCMISS_Fortran_SRC} ${OpenCMISS_CPP_HEADERS})
set(OpenCMISS_Fortran_C_SOURCES ${OpenCMISS_C_SRC} ${OpenCMISS_C_HEADERS} ${OpenCMISS_CPP_HEADERS})
set(OpenCMISS_C_SOURCES ${OpenCMISS_C_SRC} ${OpenCMISS_C_HEADERS} ${OpenCMISS_CPP_HEADERS})
set(OpenCMISS_CXX_SOURCES ${OpenCMISS_CXX_SRC} ${OpenCMISS_CXX_HEADERS} ${OpenCMISS_CPP_HEADERS})
set(OpenCMISS_SOURCES ${OpenCMISS_C_SRC} ${OpenCMISS_CXX_SRC} ${OpenCMISS_Fortran_SRC} ${OpenCMISS_C_HEADERS} ${OpenCMISS_CXX_HEADERS} ${OpenCMISS_CPP_HEADERS})
