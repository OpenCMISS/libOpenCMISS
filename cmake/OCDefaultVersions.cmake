if(LibOpenCMISS_WITH_C_BINDINGS OR LibOpenCMISS_WITH_Python_BINDINGS)
  if(DEFINED LibOpenCMISS_Python_VERSION)
    set(LibOpenCMISS_Python_VERSION "${LibOpenCMISS_Python_VERSION}" CACHE STRING "Python version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_Python_VERSION "3.8" CACHE STRING "Python version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_Python_BINDINGS)
  if(DEFINED LibOpenCMISS_SWIG_VERSION)
    set(LibOpenCMISS_SWIG_VERSION "${LibOpenCMISS_SWIG_VERSION}" CACHE STRING "SWIG version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_SWIG_VERSION "4.0" CACHE STRING "SWIG version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_MPI)
  if(DEFINED LibOpenCMISS_MPI_VERSION)
    set(LibOpenCMISS_MPI_VERSION "${LibOpenCMISS_MPI_VERSION}" CACHE STRING "External MPI library version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_MPI_VERSION "3.0" CACHE STRING "External MPI library version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_OpenMP)
  if(DEFINED LibOpenCMISS_OpenMP_VERSION)
    set(LibOpenCMISS_OpenMP_VERSION "${LibOpenCMISS_OpenMP_VERSION}" CACHE STRING "External OpenMP version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_OpenMP_VERSION "4.0" CACHE STRING "External OpenMP version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_CUDA)
  if(DEFINED LibOpenCMISS_CUDA_VERSION)
    set(LibOpenCMISS_CUDA_VERSION "${LibOpenCMISS_CUDA_VERSION}" CACHE STRING "External CUDA version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_CUDA_VERSION "12.4" CACHE STRING "External CUDA version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_CellML)
  if(DEFINED LibOpenCMISS_CellML_VERSION)
    set(LibOpenCMISS_CellML_VERSION "${LibOpenCMISS_CellML_VERSION}" CACHE STRING "External CellML library version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_CellML_VERSION "0.5.0" CACHE STRING "External CellML library version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_FieldML)
  if(DEFINED LibOpenCMISS_FieldML_VERSION)
    set(LibOpenCMISS_FieldML_VERSION "${LibOpenCMISS_FieldML_VERSION}" CACHE STRING "External FieldML library version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_FieldML_VERSION "0.5.0" CACHE STRING "External FieldML library version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_HYPRE)
  if(DEFINED LibOpenCMISS_HYPRE_VERSION)
    set(LibOpenCMISS_HYPRE_VERSION "${LibOpenCMISS_HYPRE_VERSION}" CACHE STRING "External HYPRE library version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_HYPRE_VERSION "2.23.0" CACHE STRING "External HYPRE library version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_MUMPS)
  if(DEFINED LibOpenCMISS_MUMPS_VERSION)
    set(LibOpenCMISS_MUMPS_VERSION "${LibOpenCMISS_MUMPS_VERSION}" CACHE STRING "External MUMPS library version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_MUMPS_VERSION "5.7.0" CACHE STRING "External MUMPS library version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_ParMETIS)
  if(DEFINED LibOpenCMISS_ParMETIS_VERSION)
    set(LibOpenCMISS_ParMETIS_VERSION "${LibOpenCMISS_ParMETIS_VERSION}" CACHE STRING "External ParMETIS library version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_ParMETIS_VERSION "4.0.3" CACHE STRING "External ParMETIS library version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_PETSc)
  if(DEFINED LibOpenCMISS_PETSc_VERSION)
    set(LibOpenCMISS_PETSc_VERSION "${LibOpenCMISS_PETSc_VERSION}" CACHE STRING "External PETSc library version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_PETSc_VERSION "3.21.0" CACHE STRING "External PETSc library version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_SLEPc)
  if(DEFINED LibOpenCMISS_SLEPc_VERSION)
    set(LibOpenCMISS_SLEPc_VERSION "${LibOpenCMISS_SLEPc_VERSION}" CACHE STRING "External SLEPc library version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_SLEPc_VERSION "3.21.0" CACHE STRING "External SLEPc library version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
if(LibOpenCMISS_WITH_SUNDIALS)
  if(DEFINED LibOpenCMISS_SUNDIALS_VERSION)
    set(LibOpenCMISS_SUNDIALS_VERSION "${LibOpenCMISS_SUNDIALS_VERSION}" CACHE STRING "External SUNDIALS library version number to use when building OpenCMISS libraries." FORCE)
  else()
    set(LibOpenCMISS_SUNDIALS_VERSION "6.0.0" CACHE STRING "External SUNDIALS library version number to use when building OpenCMISS libraries." FORCE)
  endif()
endif()
