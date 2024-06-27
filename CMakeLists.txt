#[=======================================================================[.rst:
libOpenCMISS CMake Build
========================

The CMake build system to build the OpenCMISS library.


#]=======================================================================]

cmake_minimum_required (VERSION 3.26 FATAL_ERROR)

if(NOT DEFINED LibOpenCMISS_ROOT)
  set(LibOpenCMISS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "LibOpenCMISS build system root" FORCE)
endif()

set(LibOpenCMISS_INC_DIR "${LibOpenCMISS_ROOT}/include")
set(LibOpenCMISS_VERSION_H "${LibOpenCMISS_INC_DIR}/opencmiss_version.h")
if(EXISTS "${LibOpenCMISS_VERSION_H}")
  file(STRINGS "${LibOpenCMISS_VERSION_H}" LibOpenCMISS_VERSION_H_LINES)
  foreach(line ${LibOpenCMISS_VERSION_H_LINES})
    if(line MATCHES "^#define[ \t](LibOpenCMISS_MAJOR_VERSION|LibOpenCMISS_MINOR_VERSION|LibOpenCMISS_PATCH_VERSION)[ \t]*([0-9]+)$")
      set(${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
    endif()
  endforeach()
else()
  message(FATAL_ERROR "OpenCMISS version file not found? Invalid setup?"
endif()
set(LibOpenCMISS_VERSION "${LibOpenCMISS_MAJOR_VERSION}.${LibOpenCMISS_MINOR_VERSION}.${LibOpenCMISS_PATCH_VERSION}")

project(LibOpenCMISS
  VERSION "${LibOpenCMISS_VERSION}"
  LANGUAGES C CXX Fortran
  DESCRIPTION "OpenCMISS computional library"
  HOMEPAGE_URL "https://www.opencmiss.org/")  
)

# Set additional directory variables
set(LibOpenCMISS_BINDINGS_DIR "${LibOpenCMISS_ROOT}/bindings")
set(LibOpenCMISS_CMAKE_DIR "${LibOpenCMISS_ROOT}/cmake")
set(LibOpenCMISS_DOC_DIR "${LibOpenCMISS_ROOT}/doc")
set(LibOpenCMISS_DOXYGEN_DIR "${LibOpenCMISS_DOC_DIR}/doxygen")
set(LibOpenCMISS_GEN_BINDINGS_DIR "${CMAKE_CURRENT_BINARY_DIRECTORY}/bindings")
set(LibOpenCMISS_GEN_INC_DIR "${CMAKE_CURRENT_BINARY_DIRECTORY}/include")
set(LibOpenCMISS_Fortran_MODULES_DIR "${CMAKE_CURRENT_BINARY_DIRECTORY}/fortran_modules")
set(LibOpenCMISS_GEN_SRC_DIR "${CMAKE_CURRENT_BINARY_DIRECTORY}/src")
set(LibOpenCMISS_SRC_DIR "${LibOpenCMISS_ROOT}/src")
set(LibOpenCMISS_TESTS_DIR "${LibOpenCMISS_ROOT}/tests")
set(LibOpenCMISS_UTILS_DIR "${LibOpenCMISS_ROOT}/utils")

# Add in library CMake modules to the module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${LibOpenCMISS_CMAKE_DIR}")

set(_CORRECT_CMAKE_MODULE_PATH FALSE)
# First check if the CMAKE_MODULE_PATH is already set properly.
foreach(_PATH ${CMAKE_MODULE_PATH})
    if (EXISTS "${_PATH}/OpenCMISS/OCMiscFunctions.cmake")
        set(_CORRECT_CMAKE_MODULE_PATH TRUE)
        if (NOT DEFINED IRON_CMAKE_MODULE_PATH)
            set(IRON_CMAKE_MODULE_PATH "${_PATH}" CACHE STRING "Location of OpenCMISS CMake modules.")
        endif ()
        break()
    endif ()
endforeach()

if (NOT _CORRECT_CMAKE_MODULE_PATH)
    # Check bare minimum configuration
    set(IRON_CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}" CACHE STRING "Location of OpenCMISS CMake modules.")
    if (EXISTS "${IRON_CMAKE_MODULE_PATH}/OpenCMISS/OCMiscFunctions.cmake")
        list(APPEND CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}/FindModuleWrappers")
        list(APPEND CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}")
        list(APPEND CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}/OpenCMISS")
        set(_CORRECT_CMAKE_MODULE_PATH TRUE)
    endif ()
endif ()

if (NOT _CORRECT_CMAKE_MODULE_PATH)
    message(STATUS "CMAKE_MODULE_PATH: '${CMAKE_MODULE_PATH}'")
    message(STATUS "IRON_CMAKE_MODULE_PATH: '${IRON_CMAKE_MODULE_PATH}'")
    message(FATAL_ERROR "'OCMiscFunctions.cmake' was not found on the module path.  Check the value of CMAKE_MODULE_PATH and/or IRON_CMAKE_MODULE_PATH.")
endif ()

# Check multi-config environment
get_property(LibOpenCMISS_HAVE_MULTICONFIG_ENV GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(LibOpenCMISS_HAVE_MULTICONFIG_ENV)
  string(TOUPPER "$<$<CONFIG>>" LibOpenCMISS_CONFIG_TYPE)
else()
  string(TOUPPER "${CMAKE_BUILD_TYPE}" LibOpenCMISS_CONFIG_TYPE)
endif()
if(LibOpenCMISS_CONFIG_TYPE MATCHES "^(DEBUG|RELWITHDEBINFO)$")
  set(LibOpenCMISS_IS_DEBUG_BUILD TRUE)
else()
  set(LibOpenCMISS_IS_DEBUG_BUILD FALSE)
endif()

# Options for LibOpenCMISS

# Options for installation
include(GNUInstallDirs)
set(LibOpenCMISS_INSTALL_BIN_DIR "${CMAKE_INSTALL_BIN}" CACHE PATH "Installation directory for OpenCMISS executables.")
set(LibOpenCMISS_INSTALL_LIB_DIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "Installation directory for OpenCMISS libraries.")
set(LibOpenCMISS_INSTALL_INC_DIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "Installation directory for OpenCMISS headers.")
set(LibOpenCMISS_INSTALL_DATA_DIR "${CMAKE_INSTALL_DATAROOTDIR}" CACHE PATH "Installation directory for OpenCMISS data.")
set(LibOpenCMISS_INSTALL_DOC_DIR "${CMAKE_INSTALL_DOCDIR}" CACHE PATH "Installation directory for OpenCMISS documents.")
set(LibOpenCMISS_INSTALL_MAN_DIR "${CMAKE_INSTALL_MANDIR}" CACHE PATH "Installation directory for OpenCMISS manual pages.")
set(LibOpenCMISS_INSTALL_CMAKE_DIR "${LibOpenCMISS_INSTALL_LIB_DIR}/cmake" CACHE PATH "Installation directory for OpenCMISS manual pages.")
set(LibOpenCMISS_INSTALL_PKGCONFIG_DIR "${LibOpenCMISS_INSTALL_DATA_DIR}/pkgconfig" CACHE PATH "Installation directory for OpenCMISS pkgconfig (.pc) files.")
# Options for library type
if(DEFINED BUILD_SHARED_LIBS)
  if(BUILD_SHARED_LIBS)
    set(LibOpenCMISS_BUILD_SHARED_LIBRARY ON CACHE BOOL "Build shared OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_BUILD_SHARED_LIBRARY OFF CACHE BOOL "Build shared OpenCMISS libraries." FORCE)
  endif()
else()
  set(LibOpenCMISS_BUILD_SHARED_LIBRARY ON CACHE BOOL "Build shared OpenCMISS libraries." FORCE)
endif()
set(LibOpenCMISS_BUILD_STATIC_LIBRARY OFF CACHE BOOL "Build static OpenCMISS libraries." FORCE)
# Options for data types
option(LibOpenCMISS_WITH_SINGLE_PRECISION "Build OpenCMISS with single precision real numbers." OFF)
option(LibOpenCMISS_WITH_LARGE_INDICES "Build OpenCMISS with extended integers for larger indicies." OFF)
# Options for language bindings
option(LibOpenCMISS_WITH_C_BINDINGS "Build the C bindings for OpenCMISS." ON)
option(LibOpenCMISS_WITH_Python_BINDINGS "Build the Python bindings for OpenCMISS." ON)
option(LibOpenCMISS_WITH_Python_VIRTUALENV "Build the Python bindings vitual environment for OpenCMISS." OFF)
option(LibOpenCMISS_WITH_DOCKER "Build the docker environment for OpenCMISS." OFF)
# Options for parallelism
option(LibOpenCMISS_WITH_MPI "Build OpenCMISS with MPI support." ON)
option(LibOpenCMISS_WITH_OpenMP "Build OpenCMISS with OpenMP support." ON)
option(LibOpenCMISS_WITH_CUDA "Build OpenCMISS with CUDA support." ON)
# Options for diagnostics, checking, profiling etc. 
if(NOT DEFINED LibOpenCMISS_WITH_DIAGNOSTICS)
  if(LibOpenCMISS_IS_DEBUG_BUILD)
    set(LibOpenCMISS_WITH_DIAGNOSTICS ON)
  else()
    set(LibOpenCMISS_WITH_DIAGNOSTICS OFF)
  endif()
else()
  set(LibOpenCMISS_WITH_DIAGNOSTICS ${LibOpenCMISS_WITH_DIAGNOSTICS})
endif()
if(DEFINED LibOpenCMISS_IS_DEBUG_BUILD)
  opton(LibOpenCMISS_WITH_DIAGNOSTICS "Build OpenCMISS with diagnostics support." ON)
else()
  option(LibOpenCMISS_WITH_DIAGNOSTICS "Build OpenCMISS with diagnostics support." OFF)
endif()
if(NOT DEFINED LibOpenCMISS_WITH_NO_CHECKS)
  set(LibOpenCMISS_WITH_NO_CHECKS OFF)
else()
  set(LibOpenCMISS_WITH_NO_CHECKS ${LibOpenCMISS_WITH_NO_CHECKS})
endif()
option(LibOpenCMISS_WITH_NO_CHECKS "Build OpenCMISS with no pre- or post-checks" OFF)
if(NOT DEFINED LibOpenCMISS_WITH_NO_PRECHECKS)
  set(LibOpenCMISS_WITH_NO_PRECHECKS OFF)
else()
  set(LibOpenCMISS_WITH_NO_PRECHECKS ${LibOpenCMISS_WITH_NO_PRECHECKS})
endif()
option(LibOpenCMISS_WITH_NO_PRECHECKS "Build OpenCMISS with no pre-checks" OFF)
if(NOT DEFINED LibOpenCMISS_WITH_NO_POSTCHECKS)
  set(LibOpenCMISS_WITH_NO_POSTCHECKS OFF)
else()
  set(LibOpenCMISS_WITH_NO_POSTCHECKS ${LibOpenCMISS_WITH_NO_POSTCHECKS})
endif()
option(LibOpenCMISS_WITH_NO_POSTCHECKS "Build OpenCMISS with no post-checks." OFF)
# Options for compiling
if(NOT DEFINED LibOpenCMISS_WITH_CHECK_ALL)
  if(LibOpenCMISS_IS_DEBUG_BUILD)
    set(LibOpenCMISS_WITH_CHECK_ALL ON)
  else()
    set(LibOpenCMISS_WITH_CHECK_ALL OFF)
  endif()
else()
  set(LibOpenCMISS_WITH_CHECK_ALL ${LibOpenCMISS_WITH_CHECK_ALL})
endif()
option(LibOpenCMISS_WITH_CHECK_ALL "Build OpenCMISS with all checks enabled." OFF)
option(LibOpenCMISS_WITH_WARN_ALL "Build OpenCMISS with all warnings enabled." ON)
option(LibOpenCMISS_WITH_PROFILING "Build OpenCMISS with profiling enbabled." OFF)
# Options for testing
option(LibOpenCMISS_WITH_TESTS "Build OpenCMISS with tests enabled." ON)
# Options for external functionality
option(LibOpenCMISS_WITH_CellML "Build OpenCMISS with CellML support." ON)
option(LibOpenCMISS_WITH_FieldML "Build OpenCMISS with FieldML support." ON)
option(LibOpenCMISS_WITH_HYPRE "Build OpenCMISS with HYPRE integration." OFF)
option(LibOpenCMISS_WITH_MUMPS "Build OpenCMISS with MUMPS integration." OFF)
option(LibOpenCMISS_WITH_ParMETIS "Build OpenCMISS with ParMETIS integration." ON)
option(LibOpenCMISS_WITH_PETSc "Build OpenCMISS with PETSc integration." ON)
option(LibOpenCMISS_WITH_SLEPc "Build OpenCMISS with SLEPc integration." ON)
option(LibOpenCMISS_WITH_SUNDIALS "Build OpenCMISS with SUNDIALS integration." OFF)

# Set versions of external libraries/packages
include(OCDefaultVersions)

# Check consistency of options
if(LibOpenCMISS_WITH_Python_BINDINGS)
  # Check we are building a shared library if Python bindings are turned on
  if(NOT LibOpenCMISS_BUILD_SHARED_LIBRARY)
    message(WARNING "Building with Python bindings requires a shared libary. Turning building a shared library on.")
    SET(LibOpenCMISS_BUILD_SHARED_LIBRARY ON)
  endif()
  # Check we have C bindings if Python bindings are turned on
  if(NOT LibOpenCMISS_WITH_C_BINDINGS)
    message(WARNING "Building with Python bindings requires C bindings. Turning C bindings on.")
    SET(LibOpenCMISS_WITH_C_BINDINGS ON)
  endif()
endif()
if(LibOpenCMISS_WITH_Python_VIRTUALENV)
  if(NOT LibOpenCMISS_WITH_Python_BINDINGS)
    # Check we have Python bindings if Python a Python virtual environment is turned on
    message(WARNING "Building with Python virtual environment requires Python bindings. Turning Python bindings on.")
    SET(LibOpenCMISS_WITH_Python_BINDINGS ON)
  endif()
endif()

# Set other variables from options
if(LibOpenCMISS_WITH_SINGLE_PRECISION)
  set(LibOpenCMISS_BUILD_PRECISION "SINGLE")
else()
  set(LibOpenCMISS_BUILD_PRECISION "DOUBLE")
endif()

option(PRINT_CONFIG_SUMMARY "Show a summary of the configuration." TRUE)

string(TOUPPER "${REAL_PRECISION}" _uppercase_REAL_PRECISION)
if(_uppercase_REAL_PRECISION STREQUAL "SINGLE")
  set(WORKING_REAL_PRECISION_SINGLE)
elseif(_uppercase_REAL_PRECISION STREQUAL "DOUBLE")
  set(WORKING_REAL_PRECISION_DOUBLE)
else()
  message(STATUS "The specified real precision of '${REAL_PRECISION}' is invalid. Defaulting to double.")
  set(WORKING_REAL_PRECISION_DOUBLE)
endif()

include(OCMiscFunctions)
include(OCMultiConfigEnvironment)

# Find and check required packages etc.
# C bindings
if(LibOpenCMISS_WITH_C_BINDINGS)
  find_package(Python ${LibOpenCMISS_Python_VERSION} COMPONENTS Interpreter)
  if (NOT Python_Interpreter_FOUND)
    message(FATAL_ERROR "No Python found. A working Python interpreter is required to generate the OpenCMISS C bindings.")
  endif()
endif()
# Python bindings
if(LibOpenCMISS_WITH_Python_BINDINGS)
  find_package(Python ${LibOpenCMISS_Python_VERSION} COMPONENTS Interpreter Development NumPy)
  if (NOT Python_Interpreter_FOUND)
    message(FATAL_ERROR "No Python found. A working Python interpreter is required to generate the OpenCMISS Python bindings.")
  endif()
  if(NOT Python_Development_FOUND)
    message(FATAL_ERROR "No Python development libraries found. Unable to generate OpenCMISS Python bindings.")
  endif()
  if(NOT Python_NumPy_FOUND)
    message(FATAL_ERROR "No Python NumPy module found. Unable to generate OpenCMISS Python bindings.")
  endif()
  find_package(SWIG ${LibOpenCMISS_SWIG_VERSION} COMPONENTS Python QUIET)
  if(NOT SWIG_FOUND)
    message(FATAL_ERROR "SWIG could not be found. Unable to generate OpenCMISS Python bindings.")
  endif()
  set(UseSWIG_MODULE_VERSION 2)
  if(Python_VERSION_MAJOR STREQUAL "3")
    set(LibOpenCMISS_PYTHONLIBS_SWIG_FLAGS -py3)
  endif()
endif()
if(LibOpenCMISS_WITH_Python_VIRTUALENV)
  # Bindings options
  find_program(LibOpenCMISS_VIRTUALENV_EXEC virtualenv)
  mark_as_advanced(LibOpenCMISS_VIRTUALENV_EXEC)
  if(NOT EXISTS "${LibOpenCMISS_VIRTUALENV_EXEC}")
    message(WARNING "The program virtualenv was not found. Disabling Python virtual environment generation.")
    set(LibOpenCMISS_WITH_Python_VIRTUALENV OFF FORCE)
  endif()
endif()
# CUDA
if(LibOpenCMISS_WITH_CUDA)
  if(NOT CUDA_FOUND)
    message(WARNING "CUDA Toolkit was not found. Disabling CUDA support in OpenCMISS.")
    set(LibOpenCMISS_WITH_CUDA OFF)
  endif()
endif()
# MPI
if(LibOpenCMISS_WITH_MPI)
  find_package(MPI ${LibOpenCMISS_MPI_VERSION} COMPONENTS C CXX Fortran REQUIRED)
endif()
# OpenMP
if(LibOpenCMISS_WITH_OpenMP)
endif()
# CellML
if(LibOpenCMISS_WITH_CellML)
  find_package(LibCellML ${LibOpenCMISS_CellML_VERSION} REQUIED)
endif()
# FieldML
if(LibOpenCMISS_WITH_FieldML)
  find_package(FieldML ${LibOpenCMISS_FieldML_VERSION} REQUIRED) 
endif()
# ParMETIS
if(LibOpenCMISS_WITH_ParMETIS)
  find_package(ParMETIS ${LibOpenCMISS_ParMETIS_VERSION} REQUIRED)
endif()
# PETSc
if(LibOpenCMISS_WITH_PETSc)
  find_package(PETSc ${LibOpenCMISS_PETSc_VERSION} REQUIRED)
endif()
# SLEPc
if(LibOpenCMISS_WITH_SLEPc)
  find_package(SLEPc ${LibOpenCMISS_SLEPc_VERSION} REQUIRED)
endif()

# Populates the LibOpenCMISS_XX_SRC variables
include(Sources)

set(BUILD_SHARED_LIBS ${LibOpenCMISS_BUILD_SHARED_LIBRARY})
# See http://cmake.3232098.n2.nabble.com/How-to-install-Fortran-module-files-td7584399.html
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LibOpenCMISS_Fortran_MODULES_DIR})
SET(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)

set(LibOpenCMISS_CLEANUP_FILES )
set(LibOpenCMISS_INCLUDE_INSTALL_FILES )
set(LibOpenCMISS_INSTALL_TARGETS )

# Handle the main Fortran libOpenCMISS library

# Platform-checks
include(CheckSymbolExists)
CHECK_SYMBOL_EXISTS(SA_NODEFER "signal.h" HAVE_SA_NODEFER)
CHECK_SYMBOL_EXISTS(sigaction "signal.h" HAVE_SIGACTION_STRUCT)

set(LibOpenCMISS_CMISS_H "${LibOpenCMISS_GEN_INC_DIR}/cmiss.h")
configure_file(${LibOpenCMISS_CMAKE_DIR}/cmiss.h.in ${LibOpenCMISS_CMISS_H)
list(APPEND LibOpenCMISS_CLEANUP_FILES "${LibOpenCMISS_CMISS_H}")

set(LibOpenCMISS_Fortran_TARGETS )
if(MSVC)
  # Need separate C/Fortran code due to Project natures in VS
  if(LibOpenCMISS_BUILD_SHARED_LIBRARY)
    list(APPEND LibOpenCMISS_Fortran_TARGETS LibOpenCMISS LibOpenCMISSC)
    add_library(LibOpenCMISS SHARED ${LibOpenCMISS_Fortran_SRC})
    add_library(LibOpenCMISSC SHARED ${LibOpenCMISS_C_SRC})
    target_link_libraries(LibOpenCMISS PRIVATE LibOpenCMISSC)
  endif()
  if(LibOpenCMISS_BUILD_STATIC_LIBRARY)
    list(APPEND LibOpenCMISS_Fortran_TARGETS LibOpenCMISS-static LibOpenCMISSC-static)
    add_library(LibOpenCMISS-static STATIC ${LibOpenCMISS_Fortran_SRC})
    add_library(LibOpenCMISSC-static STATIC ${LibOpenCMISS_C_SRC})
    target_link_libraries(LibOpenCMISS-static PRIVATE LibOpenCMISSC-static)
  endif()   
  # IDE stuff
  source_group(headers FILES ${LibOpenCMISS_HEADERS})
  source_group(src FILES ${LibOpenCMISS_C_SRC})
  source_group(src FILES ${LibOpenCMISS_Fortran_SRC})
else()
  if(LibOpenCMISS_BUILD_SHARED_LIBRARY)
    list(APPEND LibOpenCMISS_Fortran_TARGETS LibOpenCMISS)
    add_library(LibOpenCMISS SHARED ${LibOpenCMISS_SRC})
  endif()
  if(LibOpenCMISS_BUILD_STATIC_LIBRARY)
    list(APPEND LibOpenCMISS_Fortran_TARGETS LibOpenCMISS-static)
    add_library(LibOpenCMISS-static STATIC ${LibOpenCMISS_SRC})
  endif()
  source_group(headers FILES ${LibOpenCMISS_HEADERS})
  source_group(src FILES ${LibOpenCMISS_Fortran_SRC} ${LibOpenCMISS_C_SRC})
endif()

# Set the same stuff for either the separate targets (windows) or one target (else)
foreach(_target ${LibOpenCMISS_Fortran_TARGETS})
  target_include_directories(${_target}
    PRIVATE
    $<BUILD_INTERFACE:${LibOpenCMISS_INC_DIR}>
    $<BUILD_INTERFACE:${LibOpenCMISS_GEN_INC_DIR}>
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_INSTALL_PREFIX}/${LibOpenCMISS_INSTALL_INC_DIR}>
  )
  # Dont use "USE MPI" instructions if the precompiled mpi module files are incompatible.
  # This happens when the fortran compiler used to build mpi differs from the currently used one.
  # This holds on windows (only pre-compiled mpi implementations available; MPICH2 and MSMPI)
  # and on e.g. clusters with inappropriately configured toolchain environments.
  if(NOT MPI_Fortran_MODULE_COMPATIBLE)
    target_compile_definitions(${_target} PUBLIC $<BUILD_INTERFACE:NOMPIMOD>)
  endif()
  if(LibOpenCMISS_IS_DEBUG_BUILD)
    target_compile_definitions(${_target} PUBLIC DEBUG)
  endif()
  if(WORKING_REAL_PRECISION_SINGLE)
    target_compile_definitions(${_target} PUBLIC SINGLE_REAL_PRECISION)
  else()
    target_compile_definitions(${_target} PUBLIC DOUBLE_REAL_PRECISION)
  endif()
  if(LibOpenCMISS_WITH_LARGE_INDICES)
    target_compile_definitions(${_target} PUBLIC LARGE_INDICES)
  else()
    target_compile_definitions(${_target} PUBLIC STANDARD_INDICES)
  endif()
  if(LibOpenCMISS_WITH_MPI)
    target_compile_definitions(${_target} PRIVATE WITH_MPI)
    target_link_libraries(${_target} PUBLIC MPI::MPI_Fortran)
  endif()
  if(LibOpenCMISS_WITH_OpenMP)
    target_compile_definitions(${_target} PRIVATE WITH_OpenMP)
  endif()
  if(LibOpenCMISS_WITH_CUDA)
    target_compile_definitions(${_target} PRIVATE WITH_CUDA)
  endif()
  if(LibOpenCMISS_WITH_DIAGNOSTICS)
    target_compile_definitions(${_target} PUBLIC WITH_DIAGNOSTICS)
  endif()
  if(LibOpenCMISS_WITH_NO_CHECKS)
    target_compile_definitions(${_target} PUBLIC WITH_NO_CHECKS)
  endif()
  if(LibOpenCMISS_WITH_NO_PRECHECKS)
    target_compile_definitions(${_target} PUBLIC WITH_NO_PRECHECKS)
  endif()
  if(LibOpenCMISS_WITH_NO_POSTCHECKS)
    target_compile_definitions(${_target} PUBLIC WITH_NO_POSTCHECKS)
  endif()
  if(LibOpenCMISS_WITH_CELLML)
    target_compile_definitions(${_target} PRIVATE WITH_CELLML)
    target_link_libraries(${_target} PUBLIC libCellML::libCellML)
  endif()
  if(LibOpenCMISS_WITH_FIELDML)
    target_compile_definitions(${_target} PRIVATE WITH_FIELDML)
    target_link_libraries(${_target} PUBLIC FieldML::FieldML)
  endif()
  if(LibOpenCMISS_WITH_HYPRE)
    target_compile_definitions(${_target} PRIVATE WITH_HYPRE)
  endif()
  if(LibOpenCMISS_WITH_MUMPS)
    target_compile_definitions(${_target} PRIVATE WITH_MUMPS)
  endif()
  if(LibOpenCMISS_WITH_PETSc)
    target_compile_definitions(${_target} PRIVATE WITH_PETSC)
    target_link_libraries(${_target} PUBLIC PETSc::PETSc)
  endif()
  if(LibOpenCMISS_WITH_SLEPc)
    target_compile_definitions(${_target} PUBLIC WITH_SLEPC)
    target_link_libraries(${_target} PUBLIC SLEPc::SLEPc)
  endif()
  if(LibOpenCMISS_WITH_SUNDIALS)
    target_compile_definitions(${_target} PRIVATE WITH_SUNDIALS)
  endif()
  list(APPEND LibOpenCMISS_INSTALL_TARGETS ${_target})
endforeach()

# Handle bindings generation with target
set(LibOpenCMISS_GEN_BINDINGS_C_DIR "${LibOpenCMISS_GEN_BINDINGS_DIR}/c")
set(LibOpenCMISS_OPENCMISS_H "${LibOpenCMISS_GEN_BINDINGS_C_DIR}/opencmiss.h")
set(LibOpenCMISS_OPENCMISS_C_F90 "${LibOpenCMISS_GEN_SRC_DIR}/opencmiss_c.F90")
set(LibOpemCMISS_GEN_C_BINDINGS_TARGETS LibOpenCMISS_OPENCMISS_H LibOpenCMISS_OPENCMISS_C_F90)

add_custom_command(
  OUTPUT ${LibOpenCMISS_GEN_C_BINDINGS_TARGETS}
  COMMENT "Generating C bindings for OpenCMISS..."
  COMMAND ${CMAKE_EXECUTABLE} -P ${LibOpenCMISS_CMAKE_DIR}/OCGenerateCBindings.cmake
  MAIN_DEPENDENCY ${LibOpenCMISS_SRC_DIR}/opencmiss.F90
  WORKING_DIRECTORY ${LibOpenCMISS_BINDINGS_DIR}
)

# Handle the C bindings and C library
if(LibOpenCMISS_WITH_C_BINDINGS)
  add_custom_target(generate_c_bindings DEPENDS ${LibOpenCMISS_GEN_C_BINDINGS_TARGETS})
  set(LibOpenCMISS_C_TARGETS )
  if(LibOpenCMISS_BUILD_SHARED_LIBRARIES)
    list(APPEND LibOpenCMISS_C_TARGETS LibOpenCMISS_C)
    add_library(LibOpenCMISS_C SHARED ${LibOpenCMISS_OPENCMISS_C_F90})
    target_link_library(LibOpenCMISS_C PUBLIC LibOpenCMISS)
  endif()
  if(LibOpenCMISS_BUILD_STATIC_LIBRARIES)
    list(APPEND LibOpenCMISS_C_TARGETS LibOpenCMISS-static_C)
    add_library(LibOpenCMISS-static_C STATIC ${LibOpenCMISS_OPENCMISS_C_F90})    
    target_link_library(LibOpenCMISS-static_C PUBLIC LibOpenCMISS-static)
  endif()
  include(GenerateExportHeader)
  foreach(_target ${LibOpenCMISS_C_TARGETS})
    target_include_directories(${_target}
      PRIVATE
      $<BUILD_INTERFACE:${LibOpenCMISS_INC_DIR}>
      $<BUILD_INTERFACE:${LibOpenCMISS_GEN_INC_DIR}>
      INTERFACE
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${LibOpenCMISS_INSTALL_INC_DIR}>
    )
    generate_export_header(${_target})
    list(APPEND LibOpenCMISS_INCLUDE_INSTALL_FILES "${CMAKE_CURRENT_BINARY_DIR}/{$_target}_export.h")
    list(APPEND LibOpenCMISS_INSTALL_TARGETS ${_target})
  endforeach()
endif()

set(LibOpenCMISS_BINDINGS_PYTHON_DIR "${LibOpenCMISS_BINDINGS_DIR}/python")
set(LibOpenCMISS_BINDINGS_PYTHON_SWIG_DIR "${LibOpenCMISS_BINDINGS_PYTHON_DIR}/swig")
set(LibOpenCMISS_GEN_BINDINGS_PYTHON_DIR "${LibOpenCMISS_GEN_BINDINGS_DIR}/python")
set(LibOpenCMISS_GEN_BINDINGS_PYTHON_SWIG_DIR "${LibOpenCMISS_GEN_BINDINGS_PYTHON_DIR}/swig")
set(LibOpenCMISS_GEN_BINDINGS_PYTHON_SWIG_WRAPPER_DIR "${LibOpenCMISS_GEN_BINDINGS_PYTHON_SWIG_DIR}/wrapper")
set(LibOpenCMISS_PYTHON_PACKAGE_NAME "opencmiss")
set(LibOpenCMISS_PYTHON_PACKAGE_MODULE_NAME "opencmiss")
set(LibOpenCMISS_PYTHON_MODULE_NAME opencmiss_python)
set(LibOpenCMISS_OPENCMISS_I "${LibOpenCMISS_BINDINGS_PYTHON_SWIG_PYTHON_DIR}/opencmiss.i")
set(LibOpenCMISS_NUMPY_I "${LibOpenCMISS_BINDINGS_PYTHON_SWIG_DIR}/numpy.i")
set(LibOpenCMISS_NUMPY_EXTRA_I "${LibOpenCMISS_BINDINGS_PYTHON_SWIG_DIR}/numpy_extra.i")
set(LibOpenCMISS_OPENCMISS_GENERATED_I "${LibOpenCMISS_GEN_BINDINGS_PYTHON_SWIG_DIR}/opencmiss_generated.i")
set(LibOpenCMISS_OPENCMISS_PY "${LibOpenCMISS_GEN_BINDINGS_PYTHON_DIR}/${LibOpenCMISS_PYTHON_PACKAGE_NAME}/${LibOpenCMISS_PYTHON_PACKAGE_MODULE_NAME}/opencmiss.py")
set(LibOpenCMISS_OPENCMISS_PYTHON_WRAPPER_C "${LibOpenCMISS_GEN_BINDINGS_PYTHON_SWIG_WRAPPER_DIR}/opencmiss_python_wrapper.c")

add_custom_command(
  OUTPUT ${LibOpenCMISS_OPENCMISS_GENERATED_I}
  COMMENT "Generating SWIG file for OpenCMISS..."
  COMMAND ${CMAKE_EXECUTABLE} -P ${LibOpenCMISS_CMAKE_DIR}/OCGenerateSWIGFile.cmake
  MAIN_DEPENDENCY ${LibOpenCMISS_SRC_DIR}/opencmiss.F90
  WORKING_DIRECTORY ${LibOpenCMISS_BINDINGS_DIR}
)

add_custom_command(
  OUTPUT ${LibOpenCMISS_OPENCMISS_PY}
  COMMENT "Generating Python module file for OpenCMISS..."
  COMMAND ${CMAKE_EXECUTABLE} -P ${LibOpenCMISS_CMAKE_DIR}/OCGeneratePythonModuleFile.cmake
  MAIN_DEPENDENCY ${LibOpenCMISS_SRC_DIR}/opencmiss.F90
  WORKING_DIRECTORY ${LibOpenCMISS_BINDINGS_DIR}
)

add_custom_command(
  OUTPUT ${LibOpenCMISS_OPENCMISS_PYTHON_WRAPPER_C}
  COMMENT "Generating SWIG Python wrapper for OpenCMISS..."
  COMMAND ${CMAKE_EXECUTABLE} -P ${LibOpenCMISS_CMAKE_DIR}/OCGenerateSWIGPythonWrapper.cmake
  MAIN_DEPENDENCY ${LibOpenCMISS_OPENCMISS_GENERATED_I}
  DEPENDS ${LibOpenCMISS_OPENCMISS_I} ${LibOpenCMISS_NUMPY_I} ${LibOpenCMISS_NUMPY_EXTRA_I}
  WORKING_DIRECTORY ${LibOpenCMISS_BINDINGS_DIR}
)

# Handle the Python bindings 
if(LibOpenCMISS_WITH_Python_BINDINGS)
  set(LibOpenCMISS_PYTHON_PACKAGE_DESTINATION_PREFIX python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/$<CONFIG>/${LibOpenCMISS_PYTHON_PACKAGE_NAME}/${LibOpenCMISS_PACKAGE_MODULE_NAME})
  file(RELATIVE_PATH _REL_INSTALL_PATH ${CMAKE_INSTALL_PREFIX}/${LibOpenCMISS_INSTALL_LIB_DIR}/${LibOpenCMISS_PYTHON_PACKAGE_DESTINATION_PREFIX} ${CMAKE_INSTALL_PREFIX}/${LibOpenCMISS_INSTALL_LIB_DIR})
  add_custom_target(generate_python_files
    DEPENDS ${LibOpenCMISS_OPENCMISS_PY} ${LibOpenCMISS_OPENCMISS_GENERATED_I} ${LibOpenCMISS_OPENCMISS_PYTHON_WRAPPER_C}
  )
  add_library(${LibOpenCMISS_PYTHON_MODULE} MODULE
    ${LibOpenCMIS_OPENCMISS_PYTHON_WRAPPER_C} ${LibOpenCMISS_OPENCMISS_H} ${LibOpencmiss_OPENCMISS_C_F90}
  )
  add_dependencies(${LibOpenCMISS_PYTHON_MODULE} LibOpenCMISS_C)
  target_include_directories(${LibOpenCMISS_PYTHON_MODULE} PRIVATE LibOpenCMISS_C Python::NumPy Python::Python)
  )
  target_link_libraries(${LibOpenCMISS_PYTHON_MODULE} PUBLIC LibOpenCMISS_C Python::Python)
  set_target_properties(${LibOpenCMISS_PYTHON_MODULE} PROPERTIES
    PREFIX "_"
    DEBUG_POSTFIX ""
    LIBRARY_OUTPUT_DIRECTORY ${LibOpenCMISS_GEN_BINDINGS_PYTHON_DIR}/${LibOpenCMISS_PYTHON_PACKAGE_NAME}/${LibOpenCMISS_PYTHON_PACKGE_MODULE_NAME}
    RUNTIME_OUTPUT_DIRECTORY ${LibOpenCMISS_GEN_BINDINGS_PYTHON_DIR}/${LibOpenCMISS_PYTHON_PACKAGE_NAME}/${LibOpenCMISS_PYTHON_PACKGE_MODULE_NAME}
  )
  if(APPLE)
    set_target_properties(${LibOpenCMISS_PYTHON_MODULE} PROPERTIES
      INSTALL_RPATH "@loader_path/${_REL_INSTALL_PATH};@loader_path/")
  elseif(UNIX)
    set_target_properties(${LibOpenCMISS_PYTHON_MODULE} PROPERTIES
      INSTALL_RPATH "\$ORIGIN/${_REL_INSTALL_PATH};\$ORIGIN/")
  elseif(WIN32 AND NOT CYGWIN)
    set_target_properties(${LibOpenCMISS_PYTHON_MODULE} PROPERTIES
      SUFFIX ".pyd"
    )
  endif()      
  # Generate Python module
  set(LibOpenCMISS_PYTHON_SETUP_PY ${LibOpenCMISS_GEN_BINDINGS_PYTHON_DIR}/setup.py)
  set(LibOpenCMISS_PYTHON_INIT_PY ${LibOpenCMISS_GEN_BINDINGS_PYTHON_DIR}/__init__.py)
  set(LibOpenCMISS_PYTHON_UTILS_PY ${LibOpenCMISS_GEN_BINDINGS_PYTHON_DIR}/_utils.py)
  configure_file(${LibOpenCMISS_BINDINGS_PYTHON_DIR}/setup.py.in ${LibOpenCMISS_PYTHON_SETUP_PY} @ONLY)
  if(WIN32)
    file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/${LibOpenCMISS_INSTALL_BIN_DIR}" NATIVE_LIBRARY_PATH)
  else()
    file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/${LibOpenCMISS_INSTALL_LIB_DIR}" NATIVE_LIBRARY_PATH)
  endif()
  configure_file(${LibOpenCMISS_BINDINGS_PYTHON_DIR}/${LibOpenCMISS_PYTHON_PACKAGE_NAME}/${LibOpenCMISS_PYTHON_PACKAGE_MODULE_NAME}/__init__.py.in ${LibOpenCMISS_PYTHON_INIT_PY} @ONLY)
  configure_file(${LibOpenCMISS_BINDINGS_PYTHON_DIR}/${LibOpenCMISS_PYTHON_PACKAGE_NAME}/${LibOpenCMISS_PYTHON_PACKAGE_MODULE_NAME}/_utils.py.in ${UTILS_PY} @ONLY)
  if(LibOpenCMISS_WITH_PYTHON_VIRTUALENV)
    # Generate virtual environment
  endif()
endif()

if(LibOpenCMISS_WITH_TESTS)
  enable_testing()
  
  # Composes a native PATH-compatible variable to use for DLL/SO finding.
  # Each extra argument is assumed a path to add. Added in the order specified.
  function(get_library_path OUTPUT_VARIABLE)
    if (WIN32)
      set(PSEP "\\;")
      set(LD_VARNAME "PATH")
    elseif(APPLE)
      set(LD_VARNAME "DYLD_LIBRARY_PATH")
      set(PSEP ":")
    elseif(UNIX)
      set(LD_VARNAME "LD_LIBRARY_PATH")
      set(PSEP ":")
    else()
      message(WARNING "get_library_path not implemented for '${CMAKE_HOST_SYSTEM}'")
    endif()
    # Load system environment - on windows its separated by semicolon, so we need to protect those 
    string(REPLACE ";" "\\;" LD_PATH "$ENV{${LD_VARNAME}}")
    foreach(_PATH ${ARGN})
      if (OPENCMISS_HAVE_MULTICONFIG_ENV)
        file(TO_NATIVE_PATH "${_PATH}/$<CONFIG>" _PATH)
      else()
        file(TO_NATIVE_PATH "${_PATH}" _PATH)
      endif()
      set(LD_PATH "${_PATH}${PSEP}${LD_PATH}")
    endforeach()
    set(${OUTPUT_VARIABLE} "${LD_VARNAME}=${LD_PATH}" PARENT_SCOPE)
  endfunction()
  
  function(add_test_path TESTNAME)
    get_library_path(LD_PATH ${ARGN})
    set_tests_properties(${TESTNAME} PROPERTIES
      ENVIRONMENT "${LD_PATH}"
    )
  endfunction()
  
  add_subdirectory(tests)
  if(LibOpenCMISS_WITH_C_BINDINGS)
    add_subdirectory(c/tests)
  endif()
  if(LibOpenCMISS_WITH_Python_BINDINGS)
    add_subdirectory(python/tests)
  endif()
endif()

# Install steps
install(TARGETS ${LibOpenCMISS_INSTALL_TARGETS}
    EXPORT opencmiss-targets
    RUNTIME DESTINATION ${LibOpenCMISS_INSTALL_BIN_DIR}
    COMPONENT Runtime
    ARCHIVE DESTINATION ${LibOpenCMISS_INSTALL_LIB_DIR}
    
    COMPONENT Runtime
    LIBRARY DESTINATION ${LibOpenCMISS_INSTALL_LIB_DIR}
    COMPONENT Development
    INCLUDES DESTINATION ${LibOpenCMISS_INSTALL_INC_DIR}/opencmiss
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION ${Lib_OpenCMISS_INSTALL_INC_DIR}/opencmiss
    COMPONENT Development
    FILES_MATCHING PATTERN "opencmiss*.mod"
)

install(EXPORT opencmiss-targets
  NAMESPACE OpenCMISS
  DESTINATION ${PACKAGE_CONFIG_DIR}
  COMPONENT Development
)
include(CMakePackageConfigHelpers)
configure_package_config_file("${LibOpenCMISS_CMAKE_DIR}/opencmiss-config-template.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/opencmiss-config.cmake"
  INSTALL_DESTINATION "cmake/"
)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/opencmiss-config-version.cmake
  VERSION ${LibOpenCMISS_VERSION}
  COMPATIBILITY AnyNewerVersion
)
install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/opencmiss-config-version.cmake
    DESTINATION ${PACKAGE_CONFIG_DIR}
    COMPONENT Development
)
# Bundle fixup for iron shared library
# We need non-system libraries shipped out with the Runtime libraries.
# This will probably also work on other platforms - need to check
if (WIN32)
    # The wild turkey is out! We need this generator-intermediate step to get a hold of the effective output location
    # of the iron shared library. This in turn in placed inside a directory that is configuration-dependent,
    # and in the cmake_install script the CMAKE_INSTALL_CONFIG_NAME helps us getting the job done.
    file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/wildturkey.cmake"
        CONTENT "set(IRON_LIB_NAME \"$<TARGET_FILE:iron>\")")
    # Hack: The BundleUtilities dont seem to work unless you give them an executable (just the dll wont work),
    # and as we dont always have the tests around we just use the cmake executable, which will always be around.
    install(FILES ${CMAKE_COMMAND}
        DESTINATION bin)
    get_filename_component(CMAKE_NAME "${CMAKE_COMMAND}" NAME)
    install(CODE "
       include(BundleUtilities)
       include(${CMAKE_CURRENT_BINARY_DIR}/\${CMAKE_INSTALL_CONFIG_NAME}/wildturkey.cmake)
       file(COPY \"${CMAKE_COMMAND}\" DESTINATION \"\${CMAKE_INSTALL_PREFIX}/bin\")
       fixup_bundle(\"\${CMAKE_INSTALL_PREFIX}/bin/${CMAKE_NAME}\" \"\${IRON_LIB_NAME}\" \"\")
       file(REMOVE \"\${CMAKE_INSTALL_PREFIX}/bin/${CMAKE_NAME}\")
       "
       COMPONENT Redist)
endif()

if (PRINT_CONFIG_SUMMARY)
	message(STATUS   "")
	message(STATUS   "=====================================================")
	message(STATUS   "Iron has been configured with the following settings:")
	message(STATUS   "")
	message(STATUS   "          Build type: '${CMAKE_BUILD_TYPE}'")
	if (IS_DEBUG_BUILD)
	  message(STATUS "      Is debug build: TRUE")
	else()
	  message(STATUS "      Is debug build: FALSE")
	endif()
        if(WORKING_REAL_PRECISION_SINGLE)
	  message(STATUS "      Real precision: SINGLE")
	else()
	  message(STATUS "      Real precision: DOUBLE")
	endif()
	if (IRON_BUILD_SHARED_LIBRARY)
	  message(STATUS "        Library type: SHARED")
	elseif(IRON_BUILD_STATIC_LIBRARY)
	  message(STATUS "        Library type: STATIC")
	endif()
	message(STATUS   "         Install dir: '${CMAKE_INSTALL_PREFIX}'")
	if (BUILD_TESTS)
	  message(STATUS "         Build tests: ON")
	else()
	  message(STATUS "         Build tests: OFF")
	endif()
	if (WITH_PROFILING)
	  message(STATUS "           Profiling: ON")
	else()
	  message(STATUS "           Profiling: OFF")
	endif()
	if (WITH_C_BINDINGS)
	  message(STATUS "          C Bindings: ON")
	else()
	  message(STATUS "          C Bindings: OFF")
	endif()
	if (WITH_Python_BINDINGS)
	  message(STATUS "     Python Bindings: ON")
	else()
	  message(STATUS "     Python Bindings: OFF")
	endif()
	if (WITH_DIAGNOSTICS)
	  message(STATUS "         Diagnostics: ON")
	else()
	  message(STATUS "         Diagnostics: OFF")
	endif()
	if (WITH_NO_CHECKS)
	  message(STATUS "           No checks: ON")
	else()
	  message(STATUS "           No checks: OFF")
	endif()
	if (WITH_NO_PRECHECKS)
	  message(STATUS "       No pre-checks: ON")
	else()
	  message(STATUS "       No pre-checks: OFF")
	endif()
	if (WITH_NO_POSTCHECKS)
	  message(STATUS "      No post-checks: ON")
	else()
	  message(STATUS "      No post-checks: OFF")
	endif()
	if (WITH_CELLML)
	  message(STATUS "              CellML: ON")
	else()
	  message(STATUS "              CellML: OFF")
	endif()
	if (WITH_FIELDML)
	  message(STATUS "             FieldML: ON")
	else()
	  message(STATUS "             FieldML: OFF")
	endif()
	if (WITH_HYPRE)
	  message(STATUS "               Hypre: ON")
	else()
	  message(STATUS "               Hypre: OFF")
	endif()
	if (WITH_SUNDIALS)
	  message(STATUS "            SUNDIALS: ON")
	else()
	  message(STATUS "            SUNDIALS: OFF")
	endif()
	if (WITH_MUMPS)
	  message(STATUS "               MUMPS: ON")
	else()
	  message(STATUS "               MUMPS: OFF")
	endif()
	if (WITH_SCALAPACK)
	  message(STATUS "           ScaLAPACK: ON")
	else()
	  message(STATUS "           ScaLAPACK: OFF")
	endif()
	if (WITH_PETSC)
	  message(STATUS "               PETSc: ON")
	else()
	  message(STATUS "               PETSc: OFF")
	endif()
	message(STATUS   "=====================================================")
	message(STATUS   "")
endif()

# Set the files for cleaning
set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES "${LibOpenCMISS_CLEANUP_FILES}")

# Summary

# Packaging
