#[=======================================================================[.rst:
OpenCMISS CMake Build
=====================

The CMake build system to build the OpenCMISS library.


#]=======================================================================]

cmake_minimum_required(VERSION 3.26 FATAL_ERROR)

if(NOT DEFINED OpenCMISS_ROOT)
  set(OpenCMISS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "OpenCMISS library root" FORCE)
endif()

set(OpenCMISS_INC_DIR "${OpenCMISS_ROOT}/include")
set(OpenCMISS_VERSION_H "${OpenCMISS_INC_DIR}/opencmiss_version.h")
if(EXISTS "${OpenCMISS_VERSION_H}")
  file(STRINGS "${OpenCMISS_VERSION_H}" OpenCMISS_VERSION_H_LINES)
  foreach(line ${OpenCMISS_VERSION_H_LINES})
    if(line MATCHES "^#define[ \t](OpenCMISS_MAJOR_VERSION|OpenCMISS_MINOR_VERSION|OpenCMISS_PATCH_VERSION)[ \t]*([0-9]+)$")
      set(${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
    endif()
  endforeach()
else()
  message(FATAL_ERROR "OpenCMISS version file not found? Invalid setup?")
endif()
set(OpenCMISS_VERSION "${OpenCMISS_MAJOR_VERSION}.${OpenCMISS_MINOR_VERSION}.${OpenCMISS_PATCH_VERSION}")

project(OpenCMISS
  VERSION "${OpenCMISS_VERSION}"
  LANGUAGES C CXX Fortran
  DESCRIPTION "OpenCMISS computional library"
  HOMEPAGE_URL "https://www.opencmiss.org/"
)

# Set additional directory variables
set(OpenCMISS_BINDINGS_DIR "${OpenCMISS_ROOT}/bindings")
set(OpenCMISS_CMAKE_DIR "${OpenCMISS_ROOT}/cmake")
set(OpenCMISS_DOC_DIR "${OpenCMISS_ROOT}/doc")
set(OpenCMISS_Doxygen_DIR "${OpenCMISS_DOC_DIR}/doxygen")
set(OpenCMISS_SRC_DIR "${OpenCMISS_ROOT}/src")
set(OpenCMISS_TESTS_DIR "${OpenCMISS_ROOT}/tests")
set(OpenCMISS_UTILS_DIR "${OpenCMISS_ROOT}/utils")

# Set additional directory variables for the build tree
set(OpenCMISS_GEN_BINDINGS_DIR "${CMAKE_CURRENT_BINARY_DIR}/bindings")
set(OpenCMISS_GEN_INC_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
set(OpenCMISS_Fortran_MODULES_DIR "${CMAKE_CURRENT_BINARY_DIR}/fortran_modules")
set(OpenCMISS_GEN_SRC_DIR "${CMAKE_CURRENT_BINARY_DIR}/src")

# Add in library cmake modules to the module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${OpenCMISS_CMAKE_DIR}")

# set(_CORRECT_CMAKE_MODULE_PATH FALSE)
# # First check if the CMAKE_MODULE_PATH is already set properly.
# foreach(_PATH ${CMAKE_MODULE_PATH})
#     if (EXISTS "${_PATH}/OpenCMISS/OCMiscFunctions.cmake")
#         set(_CORRECT_CMAKE_MODULE_PATH TRUE)
#         if (NOT DEFINED IRON_CMAKE_MODULE_PATH)
#             set(IRON_CMAKE_MODULE_PATH "${_PATH}" CACHE STRING "Location of OpenCMISS CMake modules.")
#         endif ()
#         break()
#     endif ()
# endforeach()

# if (NOT _CORRECT_CMAKE_MODULE_PATH)
#     # Check bare minimum configuration
#     set(IRON_CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}" CACHE STRING "Location of OpenCMISS CMake modules.")
#     if (EXISTS "${IRON_CMAKE_MODULE_PATH}/OpenCMISS/OCMiscFunctions.cmake")
#         list(APPEND CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}/FindModuleWrappers")
#         list(APPEND CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}")
#         list(APPEND CMAKE_MODULE_PATH "${IRON_CMAKE_MODULE_PATH}/OpenCMISS")
#         set(_CORRECT_CMAKE_MODULE_PATH TRUE)
#     endif ()
# endif ()

# if (NOT _CORRECT_CMAKE_MODULE_PATH)
#     message(STATUS "CMAKE_MODULE_PATH: '${CMAKE_MODULE_PATH}'")
#     message(STATUS "IRON_CMAKE_MODULE_PATH: '${IRON_CMAKE_MODULE_PATH}'")
#     message(FATAL_ERROR "'OCMiscFunctions.cmake' was not found on the module path.  Check the value of CMAKE_MODULE_PATH and/or IRON_CMAKE_MODULE_PATH.")
# endif ()

# Check multi-config environment
get_property(OpenCMISS_HAVE_MULTICONFIG_ENV GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(OpenCMISS_HAVE_MULTICONFIG_ENV)
  string(TOUPPER "$<$<CONFIG>>" OpenCMISS_CONFIG_TYPE)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" OpenCMISS_CONFIG_TYPE)
endif()
if(OpenCMISS_CONFIG_TYPE MATCHES "^(DEBUG|RELWITHDEBINFO)$")
  set(OpenCMISS_IS_DEBUG_BUILD TRUE)
else()
  set(OpenCMISS_IS_DEBUG_BUILD FALSE)
endif()

# Options for OpenCMISS

# Options for installation
include(GNUInstallDirs)
set(OpenCMISS_INSTALL_BIN_SUBDIR "${CMAKE_INSTALL_BINDIR}" CACHE PATH "Sub-directory under CMAKE_INSTALL_PREFIX for OpenCMISS executables.")
set(OpenCMISS_INSTALL_LIB_SUBDIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "Sub-directory under CMAKE_INSTALL_PREFIX for OpenCMISS libraries.")
set(OpenCMISS_INSTALL_INC_SUBDIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "Sub-directory under CMAKE_INSTALL_PREFIX for OpenCMISS headers.")
set(OpenCMISS_INSTALL_INC_SUBSUBDIR "/opencmiss" CACHE PATH "Sub-directory under CMAKE_INSTALL_PREFIX/OpenCMISS_INSTALL_INC_SUBDIR for OpenCMISS headers.")
set(OpenCMISS_INSTALL_DATA_SUBDIR "${CMAKE_INSTALL_DATAROOTDIR}" CACHE PATH "Sub-directory under CMAKE_INSTALL_PREFIX for OpenCMISS data.")
set(OpenCMISS_INSTALL_DOC_SUBDIR "${CMAKE_INSTALL_DOCDIR}" CACHE PATH "Sub-directory under CMAKE_INSTALL_PREFIX for OpenCMISS documents.")
set(OpenCMISS_INSTALL_MAN_SUBDIR "${CMAKE_INSTALL_MANDIR}" CACHE PATH "Sub-directory under CMAKE_INSTALL_PREFIX for OpenCMISS manual pages.")
set(OpenCMISS_INSTALL_CMAKE_CONFIG_SUBDIR "/cmake" CACHE PATH "Sub-directory under CMAKE_INSTALL_PREFIX/OpenCMISS_INSTALL_LIB_SUBDIR for OpenCMISS CMake configuration files.")
set(OpenCMISS_INSTALL_PKGCONFIG_SUBDIR "/pkconfig" CACHE PATH "Sub-directory under CMAKE_INSTALL_PREFIX/OpenCMISS_INSTALL_DATA_SUBDIR for OpenCMISS Pkg-config files.")

# Further install location variables
set(OpenCMISS_INSTALL_INCLUDE_SUBDIR "${OpenCMISS_INSTALL_INC_SUBDIR}${OpenCMISS_INSTALL_INC_SUBSUBDIR}" CACHE PATH "Installation directory for OpenCMISS include files.")
set(OpenCMISS_INSTALL_CMAKE_SUBDIR "${OpenCMISS_INSTALL_LIB_SUBDIR}${OpenCMISS_INSTALL_CMAKE_CONFIG_SUBDIR}" CACHE PATH "Installation directory for OpenCMISS cmake files.")
set(OpenCMISS_INSTALL_CMAKE_PROJECT_CONFIG_SUBDIR "${OpenCMISS_INSTALL_CMAKE_SUBDIR}/${PROJECT_NAME}-${PROJECT_VERSION}" CACHE PATH "Installation directory for OpenCMISS cmake package configure file.")
set(OpenCMISS_INSTALL_PKGCONFIG_SUBDIR "${OpenCMISS_INSTALL_DATA_SUBDIR}${OpenCMISS_INSTALL_PKGCONFIG_SUBDIR}" CACHE PATH "Installation directory for OpenCMISS pkgconfig (.pc) files.")

mark_as_advanced(OpenCMISS_INSTALL_INCLUDE_SUBDIR OpenCMISS_INSTALL_CMAKE_SUBDIR OpenCMISS_INSTALL_CMAKE_PROJECT_CONFIG_SUBDIR OpenCMISS_INSTALL_PKGCONFIG_SUBDIR)

# Options for library type
set(OpenCMISS_BUILD_SHARED_LIBRARIES ON CACHE BOOL "Build shared OpenCMISS libraries." FORCE)
set(OpenCMISS_BUILD_STATIC_LIBRARIES OFF CACHE BOOL "Build static OpenCMISS libraries." FORCE)

# Options for data types
option(OpenCMISS_WITH_SINGLE_PRECISION "Build OpenCMISS with single precision real numbers." OFF)
option(OpenCMISS_WITH_LARGE_INDICES "Build OpenCMISS with extended integers for larger indicies." OFF)

# Options for language bindings
option(OpenCMISS_WITH_C_BINDINGS "Build the C bindings for OpenCMISS." ON)
option(OpenCMISS_WITH_Python_BINDINGS "Build the Python bindings for OpenCMISS." ON)
option(OpenCMISS_WITH_Python_VIRTUALENV "Build the Python bindings vitual environment for OpenCMISS." OFF)

# Options for packaging
option(OpenCMISS_WITH_DOCKER "Build the docker environment for OpenCMISS." OFF)

# Options for parallelism
option(OpenCMISS_WITH_MPI "Build OpenCMISS with MPI support." ON)
option(OpenCMISS_WITH_OpenMP "Build OpenCMISS with OpenMP support." ON)
option(OpenCMISS_WITH_CUDA "Build OpenCMISS with CUDA support." ON)

# Options for diagnostics, checking, profiling etc. 
if(NOT DEFINED OpenCMISS_WITH_DIAGNOSTICS)
  if(OpenCMISS_IS_DEBUG_BUILD)
    set(OpenCMISS_WITH_DIAGNOSTICS ON)
  else()
    set(OpenCMISS_WITH_DIAGNOSTICS OFF)
  endif()
else()
  set(OpenCMISS_WITH_DIAGNOSTICS ${OpenCMISS_WITH_DIAGNOSTICS})
endif()
if(DEFINED OpenCMISS_IS_DEBUG_BUILD)
  option(OpenCMISS_WITH_DIAGNOSTICS "Build OpenCMISS with diagnostics support." ON)
else()
  option(OpenCMISS_WITH_DIAGNOSTICS "Build OpenCMISS with diagnostics support." OFF)
endif()
if(NOT DEFINED OpenCMISS_WITH_NO_CHECKS)
  set(OpenCMISS_WITH_NO_CHECKS OFF)
else()
  set(OpenCMISS_WITH_NO_CHECKS ${OpenCMISS_WITH_NO_CHECKS})
endif()
option(OpenCMISS_WITH_NO_CHECKS "Build OpenCMISS with no pre- or post-checks" OFF)
if(NOT DEFINED OpenCMISS_WITH_NO_PRECHECKS)
  set(OpenCMISS_WITH_NO_PRECHECKS OFF)
else()
  set(OpenCMISS_WITH_NO_PRECHECKS ${OpenCMISS_WITH_NO_PRECHECKS})
endif()
option(OpenCMISS_WITH_NO_PRECHECKS "Build OpenCMISS with no pre-checks" OFF)
if(NOT DEFINED OpenCMISS_WITH_NO_POSTCHECKS)
  set(OpenCMISS_WITH_NO_POSTCHECKS OFF)
else()
  set(OpenCMISS_WITH_NO_POSTCHECKS ${OpenCMISS_WITH_NO_POSTCHECKS})
endif()
option(OpenCMISS_WITH_NO_POSTCHECKS "Build OpenCMISS with no post-checks." OFF)

# Options for compiling
if(NOT DEFINED OpenCMISS_WITH_CHECK_ALL)
  if(OpenCMISS_IS_DEBUG_BUILD)
    set(OpenCMISS_WITH_CHECK_ALL ON)
  else()
    set(OpenCMISS_WITH_CHECK_ALL OFF)
  endif()
else()
  set(OpenCMISS_WITH_CHECK_ALL ${OpenCMISS_WITH_CHECK_ALL})
endif()
option(OpenCMISS_WITH_CHECK_ALL "Build OpenCMISS with all checks enabled." OFF)
option(OpenCMISS_WITH_WARN_ALL "Build OpenCMISS with all warnings enabled." ON)
option(OpenCMISS_WITH_PROFILING "Build OpenCMISS with profiling enbabled." OFF)

# Options for testing
option(OpenCMISS_WITH_TESTS "Build OpenCMISS with tests enabled." ON)

# Options for documentation
option(OpenCMISS_WITH_Doxygen "Build the Doxygen documentation for OpenCMISS." ON)

# Options for external functionality
option(OpenCMISS_WITH_CellML "Build OpenCMISS with CellML support." ON)
option(OpenCMISS_WITH_FieldML "Build OpenCMISS with FieldML support." ON)
option(OpenCMISS_WITH_HYPRE "Build OpenCMISS with HYPRE integration." OFF)
option(OpenCMISS_WITH_MUMPS "Build OpenCMISS with MUMPS integration." OFF)
option(OpenCMISS_WITH_ParMETIS "Build OpenCMISS with ParMETIS integration." ON)
option(OpenCMISS_WITH_PETSc "Build OpenCMISS with PETSc integration." ON)
option(OpenCMISS_WITH_SLEPc "Build OpenCMISS with SLEPc integration." ON)
option(OpenCMISS_WITH_SUNDIALS "Build OpenCMISS with SUNDIALS integration." OFF)

set(OpenCMISS_WITH_CellML OFF)
set(OpenCMISS_WITH_MUMPS OFF)
set(OpenCMISS_WITH_SUNDIALS OFF)

# Set versions of external libraries/packages
include(OCDefaultVersions)

# Check consistency of options
if(OpenCMISS_WITH_Python_BINDINGS)
  # Check we are building a shared library if Python bindings are turned on
  if(NOT OpenCMISS_BUILD_SHARED_LIBRARIES)
    message(WARNING "Building with Python bindings requires a shared libary. Turning building a shared library on.")
    SET(OpenCMISS_BUILD_SHARED_LIBRARIES ON)
  endif()
  # Check we have C bindings if Python bindings are turned on
  if(NOT OpenCMISS_WITH_C_BINDINGS)
    message(WARNING "Building with Python bindings requires C bindings. Turning C bindings on.")
    SET(OpenCMISS_WITH_C_BINDINGS ON)
  endif()
endif()
if(OpenCMISS_WITH_Python_VIRTUALENV)
  if(NOT OpenCMISS_WITH_Python_BINDINGS)
    # Check we have Python bindings if Python a Python virtual environment is turned on
    message(WARNING "Building with Python virtual environment requires Python bindings. Turning Python bindings on.")
    SET(OpenCMISS_WITH_Python_BINDINGS ON)
  endif()
endif()

# Set other variables from options
if(OpenCMISS_WITH_SINGLE_PRECISION)
  set(OpenCMISS_BUILD_PRECISION "SINGLE")
else()
  set(OpenCMISS_BUILD_PRECISION "DOUBLE")
endif()

# Find and check required packages etc.

# C bindings
if(OpenCMISS_WITH_C_BINDINGS)
  find_package(Python ${OpenCMISS_Python_VERSION} COMPONENTS Interpreter)
  if (NOT Python_Interpreter_FOUND)
    message(FATAL_ERROR "No Python found. A working Python interpreter is required to generate the OpenCMISS C bindings.")
  endif()
endif()

# Python bindings
if(OpenCMISS_WITH_Python_BINDINGS)
  find_package(Python ${OpenCMISS_Python_VERSION} COMPONENTS Interpreter Development NumPy)
  if (NOT Python_Interpreter_FOUND)
    message(FATAL_ERROR "No Python found. A working Python interpreter is required to generate the OpenCMISS Python bindings.")
  endif()
  if(NOT Python_Development_FOUND)
    message(FATAL_ERROR "No Python development libraries found. Unable to generate OpenCMISS Python bindings.")
  endif()
  if(NOT Python_NumPy_FOUND)
    message(FATAL_ERROR "No Python NumPy module found. Unable to generate OpenCMISS Python bindings.")
  endif()
  find_package(SWIG ${OpenCMISS_SWIG_VERSION} COMPONENTS python QUIET)
  if(NOT SWIG_FOUND)
    message(FATAL_ERROR "SWIG could not be found. Unable to generate OpenCMISS Python bindings.")
  endif()
  set(UseSWIG_MODULE_VERSION 2)
  #set(OpenCMISS_PYTHONLIBS_SWIG_FLAGS )
endif()
if(OpenCMISS_WITH_Python_VIRTUALENV)
  # Bindings options
  find_program(OpenCMISS_VIRTUALENV_EXEC virtualenv)
  mark_as_advanced(OpenCMISS_VIRTUALENV_EXEC)
  if(NOT EXISTS "${OpenCMISS_VIRTUALENV_EXEC}")
    message(WARNING "The program virtualenv was not found. Disabling Python virtual environment generation.")
    set(OpenCMISS_WITH_Python_VIRTUALENV OFF FORCE)
  endif()
endif()

# CUDA
if(OpenCMISS_WITH_CUDA)
  if(NOT CUDA_FOUND)
    message(WARNING "CUDA Toolkit was not found. Disabling CUDA support in OpenCMISS.")
    set(OpenCMISS_WITH_CUDA OFF)
  endif()
endif()

# MPI
if(OpenCMISS_WITH_MPI)
  find_package(MPI ${OpenCMISS_MPI_VERSION} COMPONENTS C CXX Fortran QUIET)
  if(MPI_FOUND)
    if(MPI_Fortran_FOUND)
      if(MPI_Fortran_HAVE_F08_MODULE)
	set(OpenCMISS_MPI_Fortran_DEFINES "WITH_F08_MPI")
      else()
	if(MPI_Fortran_HAVE_F90_MODULE)
	  set(OpenCMISS_MPI_Fortran_DEFINES "WITH_F90_MPI")
	elseif(MPI_Fortran_HAVE_F77_HEADER)
	  set(OpenCMISS_MPI_Fortran_DEFINES "WITH_F77_MPI")
	else()
	  message(FATAL_ERROR "Modules/Headers for MPI Fortran have not been found. Can not compile OpenCMISS.")
	endif()
      endif()
    else()
      message(FATAL_ERROR "Could not find MPI Fortran. Can not compile OpenCMISS. Please ensure MPI Fortran is installed, can compile programs and can run.")
    endif()
  else()
    message(FATAL_ERROR "Could not find MPI. Check MPI is installed and can compile and run.")
  endif()
endif()

# OpenMP
if(OpenCMISS_WITH_OpenMP)
endif()

# Doxygen
if(OpenCMISS_WITH_Doxygen)
  find_package(Doxygen REQUIRED dot)
  if(NOT DOXYGEN_FOUND)
    message(WARNING "Doxygen was not found. Disabling Doxygen documentation for OpenCMISS.")
    set(OpenCMISS_WITH_Doxygen OFF)
  endif()
endif()

# CellML
if(OpenCMISS_WITH_CellML)
  find_package(libCellML ${OpenCMISS_CellML_VERSION} REQUIRED)
endif()

# FieldML
if(OpenCMISS_WITH_FieldML)
  find_package(FieldML ${OpenCMISS_FieldML_VERSION} REQUIRED) 
endif()

# ParMETIS
if(OpenCMISS_WITH_ParMETIS)
  find_package(ParMETIS ${OpenCMISS_ParMETIS_VERSION} REQUIRED)
endif()

# PETSc
if(OpenCMISS_WITH_PETSc)
  find_package(PETSc ${OpenCMISS_PETSc_VERSION} REQUIRED)
endif()

# SLEPc
if(OpenCMISS_WITH_SLEPc)
  find_package(SLEPc ${OpenCMISS_SLEPc_VERSION} REQUIRED)
endif()

# Populates the OpenCMISS_XX_SRC variables
include(OCSources)

# Standard names, prefixes and suffixes
set(OpenCMISS_TARGET_NAME "${PROJECT_NAME}")
set(OpenCMISS_SHARED_SUFFIX "")
set(OpenCMISS_STATIC_SUFFIX "_static")
set(OpenCMISS_Fortran_SUFFIX "_Fortran")
set(OpenCMISS_C_SUFFIX "_C")
set(OpenCMISS_Python_SUFFIX "_Python")
set(OpenCMISS_Doxygen_SUFFIX "_Doxygen")
set(OpenCMISS_INCLUDE_SUBSUBDIR "opencmiss")
set(OpenCMISS_Python_PACKAGE_NAME "opencmiss")
set(OpenCMISS_Python_PACKAGE_MODULE_NAME "opencmiss")
set(OpenCMISS_Python_MODULE_NAME "${OpenCMISS_TARGET_NAME}${OpenCMISS_Python_SUFFIX}")

# Install directories


# Main target names
set(OpenCMISS_SHARED_Fortran_TARGET "${OpenCMISS_TARGET_NAME}${OpenCMISS_SHARED_SUFFIX}${OpenCMISS_Fortran_SUFFIX}")
set(OpenCMISS_STATIC_Fortran_TARGET "${OpenCMISS_TARGET_NAME}${OpenCMISS_STATIC_SUFFIX}${OpenCMISS_Fortran_SUFFIX}")
if(WIN32)
  set(OpenCMISS_SHARED_Fortran_C_TARGET "${OpenCMISS_TARGET_NAME}${OpenCMISS_SHARED_SUFFIX}${OpenCMISS_Fortran_SUFFIX}${OpenCMISS_C_SUFFIX}")
  set(OpenCMISS_STATIC_Fortran_C_TARGET "${OpenCMISS_TARGET_NAME}${OpenCMISS_STATIC_SUFFIX}${OpenCMISS_Fortran_SUFFIX}${OpenCMISS_C_SUFFIX}")
endif()
if(OpenCMISS_WITH_C_BINDINGS)
  set(OpenCMISS_SHARED_C_TARGET "${OpenCMISS_TARGET_NAME}${OpenCMISS_SHARED_SUFFIX}${OpenCMISS_C_SUFFIX}")
  set(OpenCMISS_STATIC_C_TARGET "${OpenCMISS_TARGET_NAME}${OpenCMISS_STATIC_SUFFIX}${OpenCMISS_C_SUFFIX}")
endif()
if(OpenCMISS_WITH_Python_BINDINGS)
  set(OpenCMISS_Python_MODULE_TARGET "${OpenCMISS_Python_MODULE_NAME}")
endif()

# See http://cmake.3232098.n2.nabble.com/How-to-install-Fortran-module-files-td7584399.html
SET(CMAKE_Fortran_MODULE_DIRECTORY ${OpenCMISS_Fortran_MODULES_DIR})
SET(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)
set(BUILD_SHARED_LIBS ${OpenCMISS_BUILD_SHARED_LIBRARIES})

set(OpenCMISS_CLEANUP_FILES )
set(OpenCMISS_INSTALL_TARGETS )

# Handle the main Fortran OpenCMISS library

# Platform-checks
include(CheckSymbolExists)
CHECK_SYMBOL_EXISTS(SA_NODEFER "signal.h" HAVE_SA_NODEFER)
CHECK_SYMBOL_EXISTS(sigaction "signal.h" HAVE_SIGACTION_STRUCT)

set(OpenCMISS_OPENCMISS_INIT_H "${OpenCMISS_GEN_INC_DIR}/opencmiss_init.h")
configure_file(${OpenCMISS_CMAKE_DIR}/opencmiss_init.h.in ${OpenCMISS_OPENCMISS_INIT_H})
list(APPEND OpenCMISS_CLEANUP_FILES "${OpenCMISS_OPENCMISS_INIT_H}")

set(OpenCMISS_Fortran_TARGETS )
if(MSVC)
  # Need separate C/Fortran code due to Project natures in VS
  if(OpenCMISS_BUILD_SHARED_LIBRARIES)
    add_library(${OpenCMISS_SHARED_Fortran_TARGET} SHARED ${OpenCMISS_Fortran_SOURCES})
    add_library(${OpenCMISS_SHARED_Fortran_C_TARGET} SHARED ${OpenCMISS_Fortran_C_SOURCES})
    target_link_libraries(${OpenCMISS_SHARED_Fortran_TARGET} PRIVATE ${OpenCMISS_SHARED_Fortran_C_TARGET})
    list(APPEND OpenCMISS_Fortran_TARGETS ${OpenCMISS_SHARED_Fortran_TARGET} ${OpenCMISS_SHARED_Fortran_C_TARGET})
  endif()
  if(OpenCMISS_BUILD_STATIC_LIBRARIES)
    add_library(${OpenCMISS_STATIC_Fortran_TARGET} STATIC ${OpenCMISS_Fortran_SOURCES})
    add_library(${OpenCMISS_STATIC_Fortran_C_TARGET} STATIC ${OpenCMISS_Fortran_C_SOURCES})
    target_link_libraries(${OpenCMISS_STATIC_Fortran_TARGET} PRIVATE ${OpenCMISS_STATIC_Fortran_C_TARGET})
    list(APPEND OpenCMISS_Fortran_TARGETS ${OpenCMISS_STATIC_Fortran_TARGET} ${OpenCMISS_STATIC_Fortran_C_TARGET})
  endif()   
  # IDE stuff
  source_group(headers FILES ${OpenCMISS_HEADERS})
  source_group(src FILES ${OpenCMISS_C_SRC})
  source_group(src FILES ${OpenCMISS_Fortran_SRC})
else()
  if(OpenCMISS_BUILD_SHARED_LIBRARIES)
    message(STATUS "Adding library target ${OpenCMISS_SHARED_Fortran_TARGET}...")
    add_library(${OpenCMISS_SHARED_Fortran_TARGET} SHARED ${OpenCMISS_SOURCES})
    target_sources(${OpenCMISS_SHARED_Fortran_TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/opencmiss.F90)
    get_target_property(sourcefiles ${OpenCMISS_SHARED_Fortran_TARGET} SOURCES)
    list(APPEND OpenCMISS_Fortran_TARGETS ${OpenCMISS_SHARED_Fortran_TARGET})
  endif()
  if(OpenCMISS_BUILD_STATIC_LIBRARIES)
    add_library(${OpenCMISS_STATIC_Fortran_TARGET} STATIC ${OpenCMISS_SOURCES})
    list(APPEND OpenCMISS_Fortran_TARGETS ${OpenCMISS_STATIC_Fortran_TARGET})
  endif()
  source_group(headers FILES ${OpenCMISS_HEADERS})
  source_group(src FILES ${OpenCMISS_Fortran_SRC} ${OpenCMISS_C_SRC})
endif()

# Set the same stuff for either the separate targets (windows) or one target (else)
foreach(_target ${OpenCMISS_Fortran_TARGETS})
  target_include_directories(${_target}
    PRIVATE
    ${OpenCMISS_INC_DIR}
    ${OpenCMISS_GEN_INC_DIR}
    INTERFACE
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${OpenCMISS_INSTALL_INCLUDE_SUBDIR}>
  )
  # Dont use "USE MPI" instructions if the precompiled mpi module files are incompatible.
  # This happens when the fortran compiler used to build mpi differs from the currently used one.
  # This holds on windows (only pre-compiled mpi implementations available; MPICH2 and MSMPI)
  # and on e.g. clusters with inappropriately configured toolchain environments.
  if(OpenCMISS_IS_DEBUG_BUILD)
    target_compile_definitions(${_target} PRIVATE DEBUG)
  endif()
  if(OpenCMISS_WITH_SINGLE_PRECISION)
    target_compile_definitions(${_target} PUBLIC SINGLE_REAL_PRECISION)
  else()
    target_compile_definitions(${_target} PUBLIC DOUBLE_REAL_PRECISION)
  endif()
  if(OpenCMISS_WITH_LARGE_INDICES)
    target_compile_definitions(${_target} PUBLIC LARGE_INDICES)
  else()
    target_compile_definitions(${_target} PUBLIC STANDARD_INDICES)
  endif()
  if(OpenCMISS_WITH_MPI)
    target_compile_definitions(${_target} PUBLIC WITH_MPI ${OpenCMISS_MPI_Fortran_DEFINES})    
    target_link_libraries(${_target} PUBLIC MPI::MPI_Fortran)
  endif()
  if(OpenCMISS_WITH_OpenMP)
    target_compile_definitions(${_target} PUBLIC WITH_OpenMP)
  endif()
  if(OpenCMISS_WITH_CUDA)
    target_compile_definitions(${_target} PUBLIC WITH_CUDA)
  endif()
  if(OpenCMISS_WITH_DIAGNOSTICS)
    target_compile_definitions(${_target} PRIVATE WITH_DIAGNOSTICS)
  endif()
  if(OpenCMISS_WITH_NO_CHECKS)
    target_compile_definitions(${_target} PRIVATE WITH_NO_CHECKS)
  endif()
  if(OpenCMISS_WITH_NO_PRECHECKS)
    target_compile_definitions(${_target} PRIVATE WITH_NO_PRECHECKS)
  endif()
  if(OpenCMISS_WITH_NO_POSTCHECKS)
    target_compile_definitions(${_target} PRIVATE WITH_NO_POSTCHECKS)
  endif()
  if(OpenCMISS_WITH_CellML)
    message(STATUS "Adding CellML to target ${_target}...")
    target_compile_definitions(${_target} PUBLIC WITH_CELLML)
    target_link_libraries(${_target} PUBLIC libCellML::libCellML)
  endif()
  if(OpenCMISS_WITH_FieldML)
    message(STATUS "Adding FieldML to target ${_target}...")
    target_compile_definitions(${_target} PUBLIC WITH_FIELDML)
    target_link_libraries(${_target} PUBLIC FieldML::FieldML-API)
  endif()
  if(OpenCMISS_WITH_HYPRE)
    message(STATUS "Adding HYPRE to target ${_target}...")
    target_compile_definitions(${_target} PUBLIC WITH_HYPRE)
  endif()
  if(OpenCMISS_WITH_MUMPS)
    message(STATUS "Adding MUMPS to target ${_target}...")
    target_compile_definitions(${_target} PUBLIC WITH_MUMPS)
  endif()
  if(OpenCMISS_WITH_ParMETIS)
    message(STATUS "Adding ParMETIS to target ${_target}...")
    target_compile_definitions(${_target} PUBLIC WITH_PARMETIS)
    target_link_libraries(${_target} PUBLIC ParMETIS::ParMETIS)
  endif()
  if(OpenCMISS_WITH_PETSc)
    message(STATUS "Adding PETSc to target ${_target}...")
    target_compile_definitions(${_target} PUBLIC WITH_PETSC)
    target_link_libraries(${_target} PUBLIC PETSc::PETSc)
  endif()
  if(OpenCMISS_WITH_SLEPc)
    message(STATUS "Adding SLEPc to target ${_target}...")
    target_compile_definitions(${_target} PUBLIC WITH_SLEPC)
    target_link_libraries(${_target} PUBLIC SLEPc::SLEPc)
  endif()
  if(OpenCMISS_WITH_SUNDIALS)
    message(STATUS "Adding SUNDIALS to target ${_target}...")
    target_compile_definitions(${_target} PUBLIC WITH_SUNDIALS)
  endif()
  list(APPEND OpenCMISS_INSTALL_TARGETS ${_target})
endforeach()

# Handle bindings generation with target
set(OpenCMISS_GEN_BINDINGS_C_DIR "${OpenCMISS_GEN_BINDINGS_DIR}/c")
set(OpenCMISS_OPENCMISS_H "${OpenCMISS_GEN_BINDINGS_C_DIR}${OpenCMISS_INSTALL_INC_SUBSUBDIR}/opencmiss.h")
set(OpenCMISS_OPENCMISS_C_F90 "${OpenCMISS_GEN_SRC_DIR}/opencmiss_c.F90")
set(LibOpemCMISS_GEN_C_BINDINGS_TARGETS OpenCMISS_OPENCMISS_H OpenCMISS_OPENCMISS_C_F90)

add_custom_command(
  OUTPUT "${OpenCMISS_OPENCMISS_H}" "${OpenCMISS_OPENCMISS_C_F90}"
  COMMENT "Generating C bindings for ${PROJECT_NAME}..."
  COMMAND
  ${CMAKE_COMMAND} 
  -DOpenCMISS_ROOT=${OpenCMISS_ROOT}
  -DOpenCMISS_BINDINGS_DIR=${OpenCMISS_BINDINGS_DIR}
  -DOpenCMISS_Python_EXECUTABLE=${Python_EXECUTABLE}
  -DOpenCMISS_OPENCMISS_H=${OpenCMISS_OPENCMISS_H}
  -DOpenCMISS_OPENCMISS_C_F90=${OpenCMISS_OPENCMISS_C_F90}
  -P ${OpenCMISS_CMAKE_DIR}/OCGenerateCBindings.cmake  
  DEPENDS ${OpenCMISS_SRC_DIR}/opencmiss.F90
  WORKING_DIRECTORY ${OpenCMISS_BINDINGS_DIR}
)

# Handle the C bindings and C library
if(OpenCMISS_WITH_C_BINDINGS)
  add_custom_target(generate_c_bindings DEPENDS ${OpenCMISS_GEN_C_BINDINGS_TARGETS})
  set(OpenCMISS_INSTALL_C_INCLUDE_FILES ${OpenCMISS_OPENCMISS_H})
  list(APPEND OpenCMISS_CLEANUP_FILES "${OpenCMISS_OPENCMISS_H}")
  list(APPEND OpenCMISS_CLEANUP_FILES "${OpenCMISS_OPENCMISS_C_F90}")
  set(OpenCMISS_C_TARGETS )
  if(OpenCMISS_BUILD_SHARED_LIBRARIES)
    add_library(${OpenCMISS_SHARED_C_TARGET} SHARED ${OpenCMISS_OPENCMISS_C_F90} ${OpenCMISS_OPENCMISS_H})
    target_link_libraries(${OpenCMISS_SHARED_C_TARGET} PUBLIC ${OpenCMISS_SHARED_Fortran_TARGET})
    list(APPEND OpenCMISS_C_TARGETS ${OpenCMISS_SHARED_C_TARGET})
  endif()
  if(OpenCMISS_BUILD_STATIC_LIBRARIES)
    add_library(${OpenCMISS_STATIC_C_TARGET} STATIC ${OpenCMISS_OPENCMISS_C_F90} ${OpenCMISS_OPENCMISS_H})
    target_link_libraries(${OpenCMISS_STATIC_C_TARGET} PUBLIC ${OpenCMISS_STATIC_Fortran_TARGET})
    list(APPEND OpenCMISS_C_TARGETS ${OpenCMISS_STATIC_C_TARGET})
  endif()
  include(GenerateExportHeader)
  foreach(_target ${OpenCMISS_C_TARGETS})
    target_include_directories(${_target}
      PRIVATE
      $<BUILD_INTERFACE:${OpenCMISS_INC_DIR}>
      $<BUILD_INTERFACE:${OpenCMISS_GEN_BINDINGS_C_DIR}>
      $<BUILD_INTERFACE:${OpenCMISS_GEN_INC_DIR}>
      INTERFACE
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${OpenCMISS_INSTALL_INCLUDE_SUBDIR}>
    )
    generate_export_header(${_target}
      EXPORT_MACRO_NAME OPENCMISS_C_EXPORT
      DEPRECATED_MACRO_NAME OPENCMISS_C_DEPRECATED
      NO_EXPORT_MACRO_NAME OPENCMISS_C_NO_EXPORT
      NO_DEPRECATED_MACRO_NAME OPENCMISS_C_NO_DEPRECATED
      EXPORT_FILE_NAME "${OpenCMISS_GEN_BINDINGS_C_DIR}${OpenCMISS_INSTALL_INC_SUBSUBDIR}/opencmiss_c_export.h"
    )
    list(APPEND OpenCMISS_INSTALL_C_INCLUDE_FILES "${OpenCMISS_GEN_BINDINGS_C_DIR}${OpenCMISS_INSTALL_INC_SUBSUBDIR}/opencmiss_c_export.h")
    list(APPEND OpenCMISS_INSTALL_TARGETS ${_target})
  endforeach()
endif()

set(OpenCMISS_BINDINGS_Python_DIR "${OpenCMISS_BINDINGS_DIR}/python")
set(OpenCMISS_BINDINGS_Python_SWIG_DIR "${OpenCMISS_BINDINGS_Python_DIR}/swig")
set(OpenCMISS_GEN_BINDINGS_Python_DIR "${OpenCMISS_GEN_BINDINGS_DIR}/python")
set(OpenCMISS_GEN_BINDINGS_Python_SWIG_DIR "${OpenCMISS_GEN_BINDINGS_Python_DIR}/swig")
set(OpenCMISS_GEN_BINDINGS_Python_SWIG_WRAPPER_DIR "${OpenCMISS_GEN_BINDINGS_Python_SWIG_DIR}/wrapper")
set(OpenCMISS_OPENCMISS_I "${OpenCMISS_GEN_BINDINGS_Python_SWIG_DIR}/opencmiss.i")
set(OpenCMISS_NUMPY_I "${OpenCMISS_BINDINGS_Python_SWIG_DIR}/numpy.i")
set(OpenCMISS_NUMPY_EXTRA_I "${OpenCMISS_BINDINGS_Python_SWIG_DIR}/numpy_extra.i")
set(OpenCMISS_OPENCMISS_GENERATED_I "${OpenCMISS_GEN_BINDINGS_Python_SWIG_DIR}/opencmiss_generated.i")
set(OpenCMISS_OPENCMISS_PYTHON_PY "${OpenCMISS_GEN_BINDINGS_Python_DIR}/${OpenCMISS_Python_PACKAGE_NAME}/${OpenCMISS_Python_PACKAGE_MODULE_NAME}/${OpenCMISS_Python_MODULE_NAME}.py")
set(OpenCMISS_OPENCMISS_Python_WRAPPER_C "${OpenCMISS_GEN_BINDINGS_Python_SWIG_WRAPPER_DIR}/opencmissPYTHON_wrap.c")

add_custom_command(
  OUTPUT ${OpenCMISS_OPENCMISS_I} ${OpenCMISS_OPENCMISS_GENERATED_I} 
  COMMENT "Generating SWIG files for ${PROJECT_NAME}..."
  COMMAND ${CMAKE_COMMAND}
  -DOpenCMISS_ROOT=${OpenCMISS_ROOT}
  -DOpenCMISS_BINDINGS_DIR=${OpenCMISS_BINDINGS_DIR}
  -DOpenCMISS_BINDINGS_Python_SWIG_DIR=${OpenCMISS_BINDINGS_Python_SWIG_DIR}
  -DOpenCMISS_Python_EXECUTABLE=${Python_EXECUTABLE}
  -DOpenCMISS_OPENCMISS_I=${OpenCMISS_OPENCMISS_I}
  -DOpenCMISS_OPENCMISS_GENERATED_I=${OpenCMISS_OPENCMISS_GENERATED_I}
  -DOpenCMISS_Python_MODULE_NAME=${OpenCMISS_Python_MODULE_NAME}
  -P ${OpenCMISS_CMAKE_DIR}/OCGenerateSWIGFiles.cmake
  DEPENDS ${OpenCMISS_SRC_DIR}/opencmiss.F90
  WORKING_DIRECTORY ${OpenCMISS_BINDINGS_DIR}
)

add_custom_command(
  OUTPUT ${OpenCMISS_OPENCMISS_PYTHON_PY}
  COMMENT "Generating Python module file for ${PROJECT_NAME}..."
  COMMAND ${CMAKE_COMMAND}
  -DOpenCMISS_ROOT=${OpenCMISS_ROOT}
  -DOpenCMISS_BINDINGS_DIR=${OpenCMISS_BINDINGS_DIR}
  -DOpenCMISS_Python_EXECUTABLE=${Python_EXECUTABLE}
  -DOpenCMISS_Python_MODULE_NAME=${OpenCMISS_Python_MODULE_NAME}
  -DOpenCMISS_OPENCMISS_PYTHON_PY=${OpenCMISS_OPENCMISS_PYTHON_PY}
  -P ${OpenCMISS_CMAKE_DIR}/OCGeneratePythonModuleFile.cmake
  DEPENDS ${OpenCMISS_SRC_DIR}/opencmiss.F90
  WORKING_DIRECTORY ${OpenCMISS_BINDINGS_DIR}
)

# Handle the Python bindings 
if(OpenCMISS_WITH_Python_BINDINGS)
  set(OpenCMISS_Python_DESTINATION_PREFIX "python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/$<CONFIG>")
  set(OpenCMISS_Python_PACKAGE_DESTINATION_PREFIX "${OpenCMISS_Python_DESTINATION_PREFIX}/${OpenCMISS_Python_PACKAGE_NAME}")
  set(OpenCMISS_Python_MODULE_DESTINATION_PREFIX "${OpenCMISS_Python_PACKAGE_DESTINATION_PREFIX}/${OpenCMISS_Python_PACKAGE_MODULE_NAME}")
  add_custom_target(generate_python_bindings
   COMMENT "Generating Python bindings for ${PROJECT_NAME}..."
   DEPENDS ${OpenCMISS_OPENCMISS_PYTHON_PY} ${OpenCMISS_OPENCMISS_GENERATED_I}
  )
  list(APPEND OpenCMISS_CLEANUP_FILES "${OpenCMISS_OPENCMISS_PYTHON_PY}" "${OpenCMISS_OPENCMISS_GENERATED_I} ${OpenCMISS_OPENCMISS_Python_WRAPPER_C}")
  file(MAKE_DIRECTORY "${OpenCMISS_GEN_BINDINGS_Python_SWIG_WRAPPER_DIR}")
  include(UseSWIG)
  swig_add_library(${OpenCMISS_Python_MODULE_TARGET}
    TYPE SHARED
    LANGUAGE python
    OUTPUT_DIR ${OpenCMISS_GEN_BINDINGS_Python_SWIG_WRAPPER_DIR}
    SOURCES ${OpenCMISS_OPENCMISS_I}
  )
  set_property(TARGET ${OpenCMISS_Python_MODULE_TARGET} PROPERTY SWIG_INCLUDE_DIRECTORIES ${OpenCMISS_BINDINGS_Python_SWIG_DIR})
  set_property(TARGET ${OpenCMISS_Python_MODULE_TARGET} PROPERTY SWIG_MODULE_NAME "${OpenCMISS_Python_MODULE_NAME}")
  if(DEFINED OpenCMISS_PYTHONLIBS_SWIG_FLAGS)
    set_property(TARGET ${OpenCMISS_Python_MODULE_TARGET} PROPERTY SWIG_COMPILE_OPTIONS ${OpenCMISS_PYTHONLIBS_SWIG_FLAGS})
  endif()  
  add_dependencies(${OpenCMISS_Python_MODULE_TARGET}
    generate_python_bindings
    ${OpenCMISS_SHARED_C_TARGET}
  )
  target_include_directories(${OpenCMISS_Python_MODULE_TARGET}
    PRIVATE
    ${OpenCMISS_GEN_BINDINGS_C_DIR}
    Python::NumPy
    Python::Python)
  target_link_libraries(${OpenCMISS_Python_MODULE_TARGET} PUBLIC ${OpenCMISS_SHARED_C_TARGET} Python::Python)
  set_target_properties(${OpenCMISS_Python_MODULE_TARGET} PROPERTIES
    OUTPUT_NAME ${OpenCMISS_Python_MODULE_NAME}
    LINKER_LANGUAGE C
    PREFIX "_"
    DEBUG_POSTFIX ""
    LIBRARY_OUTPUT_DIRECTORY ${OpenCMISS_GEN_BINDINGS_Python_DIR}/${OpenCMISS_Python_PACKAGE_NAME}/${OpenCMISS_Python_PACKAGE_MODULE_NAME}
    RUNTIME_OUTPUT_DIRECTORY ${OpenCMISS_GEN_BINDINGS_Python_DIR}/${OpenCMISS_Python_PACKAGE_NAME}/${OpenCMISS_Python_PACKAGE_MODULE_NAME}
  )
  list(APPEND OpenCMISS_CLEANUP_FILES "${OpenCMISS_OPENCMISS_Python_WRAPPER_C}")
  # Set RPATH of installed python module to the library directory
  #cmake_path(SET OpenCMISS_Python_MODULE_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}/${OpenCMISS_INSTALL_LIB_SUBDIR}/${OpenCMISS_Python_MODULE_DESTINATION_PREFIX}")
  #cmake_path(RELATIVE_PATH OpenCMISS_Python_MODULE_INSTALL_PATH
  #  BASE_DIRECTORY "${CMAKE_INSTALL_PREFIX}/${OpenCMISS_INSTALL_LIB_SUBDIR}"
  #  OUTPUT_VARIABLE _REL_INSTALL_PATH
  #)
  set(_REL_INSTALL_PATH "../../../../../")
  if(APPLE)
    set_target_properties(${OpenCMISS_Python_MODULE_TARGET} PROPERTIES
      INSTALL_RPATH "@loader_path/${_REL_INSTALL_PATH};@loader_path/"
    )
  elseif(UNIX)
    set_target_properties(${OpenCMISS_Python_MODULE_TARGET} PROPERTIES
      INSTALL_RPATH "\$ORIGIN/${_REL_INSTALL_PATH};\$ORIGIN/${_REL_INSTALL_PATH}/../lib;\$ORIGIN/${_REL_INSTALL_PATH}/../../mpi-none/lib;\$ORIGIN/${_REL_INSTALL_PATH}/../../mpi-none/lib64;\$ORIGIN/"
    )
  elseif(WIN32 AND NOT CYGWIN)
    set_target_properties(${OpenCMISS_Python_MODULE_TARGET} PROPERTIES
      SUFFIX ".pyd"
    )
  endif()      
  # Generate Python module
  set(OpenCMISS_Python_SETUP_PY ${OpenCMISS_GEN_BINDINGS_Python_DIR}/setup.py)
  set(OpenCMISS_Python_PACKAGE_INIT_PY ${OpenCMISS_GEN_BINDINGS_Python_DIR}/${OpenCMISS_Python_PACKAGE_NAME}/__init__.py)
  set(OpenCMISS_Python_MODULE_INIT_PY ${OpenCMISS_GEN_BINDINGS_Python_DIR}/${OpenCMISS_Python_PACKAGE_NAME}/${OpenCMISS_Python_PACKAGE_MODULE_NAME}/__init__.py)
  set(OpenCMISS_Python_MODULE_UTILS_PY ${OpenCMISS_GEN_BINDINGS_Python_DIR}/${OpenCMISS_Python_PACKAGE_NAME}/${OpenCMISS_Python_PACKAGE_MODULE_NAME}/_utils.py)
  configure_file(${OpenCMISS_BINDINGS_Python_DIR}/setup.py.in ${OpenCMISS_Python_SETUP_PY} @ONLY)
  if(WIN32)
    file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/${OpenCMISS_INSTALL_BIN_SUBDIR}" NATIVE_LIBRARY_PATH)
  else()
    file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/${OpenCMISS_INSTALL_LIB_SUBDIR}" NATIVE_LIBRARY_PATH)
  endif()
  configure_file(${OpenCMISS_BINDINGS_Python_DIR}/${OpenCMISS_Python_PACKAGE_NAME}/__init__.py.in ${OpenCMISS_Python_PACKAGE_INIT_PY} @ONLY)
  configure_file(${OpenCMISS_BINDINGS_Python_DIR}/${OpenCMISS_Python_PACKAGE_NAME}/${OpenCMISS_Python_PACKAGE_MODULE_NAME}/__init__.py.in ${OpenCMISS_Python_MODULE_INIT_PY} @ONLY)
  configure_file(${OpenCMISS_BINDINGS_Python_DIR}/${OpenCMISS_Python_PACKAGE_NAME}/${OpenCMISS_Python_PACKAGE_MODULE_NAME}/_utils.py.in ${OpenCMISS_Python_MODULE_UTILS_PY} @ONLY)
  if(OpenCMISS_WITH_Python_VIRTUALENV)
    # Generate virtual environment
  endif()
endif()

if(OpenCMISS_WITH_Doxygen)
  # Generate Doxygen documentation
  set(OpenCMISS_Doxygen_TARGET "${OpenCMISS_TARGET_NAME}${OpenCMISS_Doxygen_SUFFIX}")
  set(OpenCMISS_Doxygen_Fortran_TARGET "${OpenCMISS_TARGET_NAME}${OpenCMISS_Doxygen_SUFFIX}${OpenCMISS_Fortran_SUFFIX}")
  set(OpenCMISS_Doxygen_Programmer_TARGET "${OpenCMISS_TARGET_NAME}${OpenCMISS_Doxygen_SUFFIX}_Programmer")
  doxygen_add_docs(${OpenCMISS_Doxygen_Fortran_TARGET}
    COMMENT "Generating ${PROJECT_NAME} Fortran API Doxygen documentation..."
    WORKING_DIR ${OpenCMISS_Doxygen_DIR}
    CONFIG_FILE ${OpenCMISS_Doxygen_DIR}/Doxyfile_Fortran
  )
  doxygen_add_docs(${OpenCMISS_Doxygen_Programmer_TARGET}
    COMMENT "Generating ${PROJECT_NAME} Programmer Doxygen documentation..."
    WORKING_DIR ${OpenCMISS_Doxygen_DIR}
    CONFIG_FILE ${OpenCMISS_Doxygen_DIR}/Doxyfile_Programmer
  )
  set(OpenCMISS_Doxygen_TARGETS ${OpenCMISS_Doxygen_Fortran_TARGET} ${OpenCMISS_Doxygen_Programmer_TARGET})
  if(OpenCMISS_WITH_C_BINDINGS)
    set(OpenCMISS_Doxygen_C_TARGET "${OpenCMISS_TARGET_NAME}${OpenCMISS_Doxygen_SUFFIX}${OpenCMISS_C_SUFFIX}")
    doxygen_add_docs(${OpenCMISS_Doxygen_C_TARGET}
      COMMENT "Generating ${PROJECT_NAME} C API Doxygen documentation..."
      WORKING_DIR ${OpenCMISS_Doxygen_DIR}
      CONFIG_FILE ${OpenCMISS_Doxygen_DIR}/Doxyfile_C
    )
    list(APPEND OpenCMISS_Doxygen_TARGETS ${OpenCMISS_Doxygen_C_TARGET})
  endif()
  if(OpenCMISS_WITH_Python_BINDINGS)
    set(OpenCMISS_Doxygen_Python_TARGET "${OpenCMISS_TARGET_NAME}${OpenCMISS_Doxygen_SUFFIX}${OpenCMISS_Python_SUFFIX}")
    doxygen_add_docs(${OpenCMISS_Doxygen_Python_TARGET}
      COMMENT "Generating ${PROJECT_NAME} Python API Doxygen documentation..."
      WORKING_DIR ${OpenCMISS_Doxygen_DIR}
      CONFIG_FILE ${OpenCMISS_Doxygen_DIR}/Doxyfile_Python
    )
    list(APPEND OpenCMISS_Doxygen_TARGETS ${OpenCMISS_Doxygen_Python_TARGET})
  endif()
  add_custom_target(${OpenCMISS_Doxygen_TARGET} DEPENDS ${OpenCMISS_Doxygen_TARGETS})
endif()
# Install steps

install(TARGETS ${OpenCMISS_INSTALL_TARGETS}
  EXPORT "${PROJECT_NAME}Targets"
  RUNTIME DESTINATION ${OpenCMISS_INSTALL_BIN_SUBDIR}
  COMPONENT Runtime
  ARCHIVE DESTINATION ${OpenCMISS_INSTALL_LIB_SUBDIR}    
  COMPONENT Runtime
  LIBRARY DESTINATION ${OpenCMISS_INSTALL_LIB_SUBDIR}
  COMPONENT Development
  INCLUDES DESTINATION ${OpenCMISS_INSTALL_INCLUDE_SUBDIR}
)

install(FILES
  "${CMAKE_Fortran_MODULE_DIRECTORY}/opencmiss.mod"
  "${CMAKE_Fortran_MODULE_DIRECTORY}/opencmiss_c.mod"
  "${CMAKE_Fortran_MODULE_DIRECTORY}/opencmissfortranc.mod"
  DESTINATION ${OpenCMISS_INSTALL_INCLUDE_SUBDIR}
  COMPONENT Development
)

if(OpenCMISS_WITH_C_BINDINGS)
  install(FILES ${OpenCMISS_INSTALL_C_INCLUDE_FILES}
    DESTINATION ${OpenCMISS_INSTALL_INCLUDE_SUBDIR}
    COMPONENT CBindings
  )
endif()

if(OpenCMISS_WITH_Python_BINDINGS)
  install(FILES ${OpenCMISS_Python_SETUP_PY}
    DESTINATION ${OpenCMISS_INSTALL_LIB_SUBDIR}/${OpenCMISS_Python_PACKAGE_DESTINATION_PREFIX}
    COMPONENT PythonBindings)
  install(FILES ${OpenCMISS_Python_PACKAGE_INIT_PY}
    DESTINATION ${OpenCMISS_INSTALL_LIB_SUBDIR}/${OpenCMISS_Python_PACKAGE_DESTINATION_PREFIX}${_Python_DESTINATION_PREFIX}/${OpenCMISS_Python_PACKAGE_NAME}
    COMPONENT PythonBindings)
  install(FILES ${OpenCMISS_Python_PACKAGE_INIT_PY} ${OpenCMISS_Python_MODULE_INIT_PY} ${OpenCMISS_OPENCMISS_PYTHON_PY} ${OpenCMISS_Python_MODULE_UTILS_PY}
    DESTINATION ${OpenCMISS_INSTALL_LIB_SUBDIR}/${OpenCMISS_Python_PACKAGE_DESTINATION_PREFIX}${_Python_DESTINATION_PREFIX}/${OpenCMISS_Python_PACKAGE_NAME}/${OpenCMISS_Python_PACKAGE_MODULE_NAME}
    COMPONENT PythonBindings)
  install(TARGETS ${OpenCMISS_Python_MODULE_TARGET}
    DESTINATION ${OpenCMISS_INSTALL_LIB_SUBDIR}/${OpenCMISS_Python_PACKAGE_DESTINATION_PREFIX}${_Python_DESTINATION_PREFIX}/${OpenCMISS_Python_PACKAGE_NAME}/${OpenCMISS_Python_PACKAGE_MODULE_NAME}
    COMPONENT PythonBindings)  
endif()

install(EXPORT "${PROJECT_NAME}Targets"
  NAMESPACE OpenCMISS::
  DESTINATION ${OpenCMISS_INSTALL_CMAKE_PROJECT_CONFIG_SUBDIR}
  COMPONENT Development
)
include(CMakePackageConfigHelpers)
configure_package_config_file("${OpenCMISS_CMAKE_DIR}/${PROJECT_NAME}Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  PATH_VARS
  OpenCMISS_INSTALL_LIB_SUBDIR
  OpenCMISS_INSTALL_INCLUDE_SUBDIR
  OpenCMISS_INSTALL_BIN_SUBDIR
  OpenCMISS_INSTALL_DATA_SUBDIR
  OpenCMISS_INSTALL_DOC_SUBDIR
  OpenCMISS_INSTALL_MAN_SUBDIR
  OpenCMISS_INSTALL_CMAKE_PROJECT_CONFIG_SUBDIR
  OpenCMISS_INSTALL_PKGCONFIG_SUBDIR
  INSTALL_DESTINATION ${OpenCMISS_INSTALL_CMAKE_PROJECT_CONFIG_SUBDIR}
)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${OpenCMISS_VERSION}
  COMPATIBILITY AnyNewerVersion
)
install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${OpenCMISS_INSTALL_CMAKE_PROJECT_CONFIG_SUBDIR}
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


if(OpenCMISS_WITH_TESTS)
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
  if(OpenCMISS_WITH_C_BINDINGS)
    add_subdirectory(bindings/c/tests)
  endif()
  if(OpenCMISS_WITH_Python_BINDINGS)
    add_subdirectory(bindings/python/tests)
  endif()
endif()

message(STATUS "")
message(STATUS "=============================================================")
message(STATUS "${PROJECT_NAME} has been configured with the following settings:")
message(STATUS "")
if(NOT OpenCMISS_HAVE_MULTICONFIG_ENV)
  message(STATUS "                     Build type: '${CMAKE_BUILD_TYPE}'")
endif()
message(STATUS "           Build shared library: ${OpenCMISS_BUILD_SHARED_LIBRARIES}")
message(STATUS "           Build static library: ${OpenCMISS_BUILD_STATIC_LIBRARIES}")
message(STATUS "")
message(STATUS "                 Is debug build: ${OpenCMISS_IS_DEBUG_BUILD}")
if(${OpenCMISS_WITH_SINGLE_PRECISION})
  message(STATUS "                 Real precision: SINGLE")
else()
  message(STATUS "                 Real precision: DOUBLE")
endif()
message(STATUS "             With large indices: ${OpenCMISS_WITH_LARGE_INDICES}")
message(STATUS "")
message(STATUS "                With C bindings: ${OpenCMISS_WITH_C_BINDINGS}")
message(STATUS "           With Python bindings: ${OpenCMISS_WITH_Python_BINDINGS}")
message(STATUS "With Python virtual environment: ${OpenCMISS_WITH_Python_BINDINGS}")
message(STATUS "")
message(STATUS "              With docker image: ${OpenCMISS_WITH_DOCKER}")
message(STATUS "")
message(STATUS "                       With MPI: ${OpenCMISS_WITH_MPI}")
message(STATUS "                    With OpenMP: ${OpenCMISS_WITH_OpenMP}")
message(STATUS "                      With CUDA: ${OpenCMISS_WITH_CUDA}")
message(STATUS "")
message(STATUS "               With diagnostics: ${OpenCMISS_WITH_DIAGNOSTICS}")
message(STATUS "                 With no checks: ${OpenCMISS_WITH_NO_CHECKS}")
message(STATUS "              With no prechecks: ${OpenCMISS_WITH_NO_PRECHECKS}")
message(STATUS "             With no postchecks: ${OpenCMISS_WITH_NO_POSTCHECKS}")
message(STATUS "")
message(STATUS "                With all checks: ${OpenCMISS_WITH_CHECK_ALL}")
message(STATUS "              With all warnings: ${OpenCMISS_WITH_WARN_ALL}")
message(STATUS "                 With profiling: ${OpenCMISS_WITH_PROFILING}")
message(STATUS "")
message(STATUS "                   With testing: ${OpenCMISS_WITH_TESTS}")
message(STATUS "")
message(STATUS "     With Doxygen documentation: ${OpenCMISS_WITH_Doxygen}")
message(STATUS "")
message(STATUS "                    With CellML: ${OpenCMISS_WITH_CellML}")
message(STATUS "                   With FieldML: ${OpenCMISS_WITH_FieldML}")
message(STATUS "                     With HYPRE: ${OpenCMISS_WITH_HYPRE}")
message(STATUS "                     With MUMPS: ${OpenCMISS_WITH_MUMPS}")
message(STATUS "                  With ParMETIS: ${OpenCMISS_WITH_ParMETIS}")
message(STATUS "                     With PETSc: ${OpenCMISS_WITH_PETSc}")
message(STATUS "                     With SLEPc: ${OpenCMISS_WITH_SLEPc}")
message(STATUS "                  With SUNDIALS: ${OpenCMISS_WITH_SUNDIALS}")
message(STATUS   "")
message(STATUS "                 Install prefix: '${CMAKE_INSTALL_PREFIX}'")
message(STATUS   "")
message(STATUS   "=============================================================")
message(STATUS   "")

# Set the files for cleaning
set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES "${OpenCMISS_CLEANUP_FILES}")

# Summary

# Packaging
