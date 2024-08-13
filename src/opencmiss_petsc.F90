!> \file
!> \author Chris Bradley
!> \brief This module is a OpenCMISS buffer module to the PETSc library.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s): Chris Bradley
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \defgroup OpenCMISS_PETSc OpenCMISS::PETSc
!> This module is an OpenCMISS buffer module to the PETSc library.
MODULE OpenCMISSPETSc
  
  USE BaseRoutines
  USE OpenCMISSPetscTypes
  USE Kinds
  USE ISO_VARYING_STRING
#ifdef WITH_MPI
#ifdef WITH_F08_MPI
  USE MPI_F08
#endif  
#endif  
  USE Strings
  USE Types
  
#include "macros.h"
  
#include "petscversion.h"
#include "petsc/finclude/petsc.h"
  USE petsc
#include "petsc/finclude/petsctao.h"
  USE petsctao
  
  IMPLICIT NONE
 
  PRIVATE
  
  !Module parameters

  !> \addtogroup PETSc_PetscMatInsertModes OpenCMISS::PETSc::Constants
  !> \brief OpenCMISS PETSc constants
  !> \see OpenCMISSPETSc
  !>@{
  
  !Insert mode types
  !> \addtogroup OpenCMISSPETSc_PetscMatInsertModes OpenCMISS::PETSc::Constants::PetscMatInsertModes
  !> \brief Types of PETSc matrix insert modes
  !> \see OpenCMISSPETSc
  !>@{
  InsertMode, PARAMETER :: PETSC_NOT_SET_VALUES = NOT_SET_VALUES !<Do not actually use the values
  InsertMode, PARAMETER :: PETSC_INSERT_VALUES = INSERT_VALUES !<Values set in a matrix will overwrite any previous values
  InsertMode, PARAMETER :: PETSC_ADD_VALUES = ADD_VALUES !<Values set in a matrix will add to any previous values
  InsertMode, PARAMETER :: PETSC_MAX_VALUES = MAX_VALUES !<Use the maximum of each current value and the provided value
  InsertMode, PARAMETER :: PETSC_MIN_VALUES = MIN_VALUES !<Use the minimum of each current value and the provided value
  InsertMode, PARAMETER :: PETSC_INSERT_ALL_VALUES = INSERT_ALL_VALUES !<Insert, even if indices that are not marked as constrained by the PetscSection
  InsertMode, PARAMETER :: PETSC_ADD_ALL_VALUES = ADD_ALL_VALUES !<Add, even if indices that are not marked as constrained by the PetscSection
  InsertMode, PARAMETER :: PETSC_INSERT_BC_VALUES = INSERT_BC_VALUES !<Insert, but ignore indices that are marked as constrained by the PetscSection
  InsertMode, PARAMETER :: PETSC_ADD_BC_VALUES = ADD_BC_VALUES !<Add, but ignore indices that are marked as constrained by the PetscSection
  !>@}
  
  !Scatter mode types
  !> \addtogroup OpenCMISSPETSc_PetscScatterModes OpenCMISS::PETSc::Constants::PetscScatterModes
  !> \brief Types of PETSc matrix insert modes
  !> \see OpenCMISSPETSc
  !>@{
  ScatterMode, PARAMETER :: PETSC_SCATTER_FORWARD = SCATTER_FORWARD !<Scatters the values as dictated by the VecScatterCreate call
  ScatterMode, PARAMETER :: PETSC_SCATTER_REVERSE = SCATTER_REVERSE !<Moves the values in the opposite direction than the directions indicated in the VecScatterCreate call
  ScatterMode, PARAMETER :: PETSC_SCATTER_FORWARD_LOCAL = SCATTER_FORWARD_LOCAL !<Scatters the values as dictated by the VecScatterCreate call except NO MPI communication is done
  ScatterMode, PARAMETER :: PETSC_SCATTER_REVERSE_LOCAL = SCATTER_REVERSE_LOCAL !<Moves the values in the opposite direction than the directions indicated in the VecScatterCreate call except NO MPI communication is done.
  !>@}
  
  !> \addtogroup OpenCMISSPETSc_PetscKSP OpenCMISS::PETSc::Constants::KSP
  !> \brief Constants for PETSc KSP (Krylov Subspace) objects
  !> \see OpenCMISSPETSc
  !>@{
  !KSP types
  !> \addtogroup OpenCMISSPETSc_PetscKSPTypes OpenCMISS::PETSc::Constants::KSP::KSPType
  !> \brief Types of PETSc KSP (Krylov Subspace) solvers
  !> \see OpenCMISSPETSc
  !>@{
  KSPType, PARAMETER :: PETSC_KSP_RICHARDSON = KSPRICHARDSON !<Preconditioned Richardson iterative solver
  KSPType, PARAMETER :: PETSC_KSP_CHEBYSHEV = KSPCHEBYSHEV !<Preconditioned Chebyshev iterative solver
  KSPType, PARAMETER :: PETSC_KSP_CG = KSPCG !<Preconditioned conjugate gradient (PCG) solver
  KSPType, PARAMETER :: PETSC_KSP_GROPPCG = KSPGROPPCG !<Gropp's overlapped reduction conjugate gradient
  KSPType, PARAMETER :: PETSC_KSP_PIPECG = KSPPIPECG !<Pipelined conjugate gradient
  KSPType, PARAMETER :: PETSC_KSP_PIPECGRR = KSPPIPECGRR !<Pipelined conjugate gradient with residual replacement
  KSPType, PARAMETER :: PETSC_KSP_PIPELCG = KSPPIPELCG !<Pipelined depth l conjugate gradient
  KSPType, PARAMETER :: PETSC_KSP_PIPECG2 = KSPPIPECG2 !<Pipelined conjugate gradient over iteration pairs
  KSPType, PARAMETER :: PETSC_KSP_CGNE = KSPCGNE !<Conjugate gradient on normal equations
  KSPType, PARAMETER :: PETSC_KSP_NASH = KSPNASH !<Nash conjugate gradient with trust region constraint
  KSPType, PARAMETER :: PETSC_KSP_STCG = KSPSTCG !<Conjugate graident with trust region constraint
  KSPType, PARAMETER :: PETSC_KSP_GLTR = KSPGLTR !<Gould et al conjugate gradient with trust region constraint
  KSPType, PARAMETER :: PETSC_KSP_FCG = KSPFCG !<Flexiable conjugate gradient
  KSPType, PARAMETER :: PETSC_KSP_PIPEFCG = KSPPIPEFCG !<Pipelined flexiable conjugate gradient
  KSPType, PARAMETER :: PETSC_KSP_GMRES = KSPGMRES !<Generalised Minimal Residual solver (GMRES)
  KSPType, PARAMETER :: PETSC_KSP_PIPEFGMRES = KSPPIPEFGMRES !<Pipelined flexiable GMRES
  KSPType, PARAMETER :: PETSC_KSP_FGMRES = KSPFGMRES !<Flexible GMRES
  KSPType, PARAMETER :: PETSC_KSP_LGMRES = KSPLGMRES !<LGMRES
  KSPType, PARAMETER :: PETSC_KSP_DGMRES = KSPDGMRES !<Deflated GMRES
  KSPType, PARAMETER :: PETSC_KSP_PGMRES = KSPPGMRES !<Piplined GMRES
  KSPType, PARAMETER :: PETSC_KSP_TCQMR = KSPTCQMR !<Tony Chan's Quasi Minimal Residual solver
  KSPType, PARAMETER :: PETSC_KSP_BCGS = KSPBCGS !<Stablised BiConjugate Gradient Squared solver
  KSPType, PARAMETER :: PETSC_KSP_IBCGS = KSPIBCGS !<Improved stabilised bi-conjugate gradient
  KSPType, PARAMETER :: PETSC_KSP_QMRCGS = KSPQMRCGS !<QMR stabilised BiCGStab
  KSPType, PARAMETER :: PETSC_KSP_FBCGS = KSPFBCGS !<Flexible stabilised bi-conjugate gradient
  KSPType, PARAMETER :: PETSC_KSP_FBCGSR = KSPFBCGSR !<Flexible stabilised bi-conjugate gradient with fewer reductions
  KSPType, PARAMETER :: PETSC_KSP_BCGSL = KSPBCGSL !<Stabilised bi-conjugate gradients with lenght l recurrence
  KSPType, PARAMETER :: PETSC_KSP_PIPEBCGS = KSPPIPEBCGS !<Pipelined stabilised bi-conjugate gradients
  KSPType, PARAMETER :: PETSC_KSP_CGS = KSPCGS !<Conjugate Gradient Squared solver
  KSPType, PARAMETER :: PETSC_KSP_TFQMR = KSPTFQMR !<Transpose Free Quasi Minimum Residual solver
  KSPType, PARAMETER :: PETSC_KSP_CR = KSPCR !<Conjugate Residuals
  KSPType, PARAMETER :: PETSC_KSP_LSQR = KSPLSQR !LSQR
  KSPType, PARAMETER :: PETSC_KSP_PREONLY = KSPPREONLY !<Stub solver that only applies the preconditioner.
  KSPType, PARAMETER :: PETSC_KSP_QCG = KSPQCG !<Steinhaug conjugate gradient with trust region constraint
  KSPType, PARAMETER :: PETSC_KSP_BICG = KSPBICG !<BiConjugate Gradient solver
  KSPType, PARAMETER :: PETSC_KSP_MINRES = KSPMINRES !<Minimum Residual solver
  KSPType, PARAMETER :: PETSC_KSP_SYMMLQ = KSPSYMMLQ !<SYMMLQ
  KSPType, PARAMETER :: PETSC_KSP_LCD = KSPLCD !<Left Conjugate Direction solver
  KSPType, PARAMETER :: PETSC_KSP_GCR = KSPGCR !<Generalised Conjugate Resdiuals solver
  KSPType, PARAMETER :: PETSC_KSP_PIPEGCR = KSPPIPEGCR !<Pipelined generalised Conjugate Resdiuals solver
  KSPType, PARAMETER :: PETSC_KSP_TSIRM = KSPTSIRM !<Two-stage with least squares residual minimisation
  KSPType, PARAMETER :: PETSC_KSP_CGLS = KSPCGLS !<Conjugate gradient for least squares
  KSPType, PARAMETER :: PETSC_KSP_FETIDP = KSPFETIDP !<FETI-DP (reduction to dual-primal sub-problem)
  KSPType, PARAMETER :: PETSC_KSP_HPDDM = KSPHPDDM !<Generalised conjugate resiudal (with inner normalisation and deflated restarts)
  !>@}
  
  !KSPConvergedReason types
  !> \addtogroup OpenCMISSPETSc_PetscKSPConvergedReasons OpenCMISS::PETSc::Constants::KSP::KSPConvergedReason
  !> \brief Types of KSP converged reasons for PETSc KSP solvers
  !> \see OpenCMISSPETSc
  !>@{
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_RTOL = KSP_CONVERGED_RTOL
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_ATOL = KSP_CONVERGED_ATOL
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_ITS = KSP_CONVERGED_ITS
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_NEG_CURVE = KSP_CONVERGED_NEG_CURVE
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_STEP_LENGTH = KSP_CONVERGED_STEP_LENGTH
  KSPConvergedReason, PARAMETER :: PETSC_KSP_CONVERGED_HAPPY_BREAKDOWN = KSP_CONVERGED_HAPPY_BREAKDOWN
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_NULL = KSP_DIVERGED_NULL
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_ITS = KSP_DIVERGED_ITS
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_DTOL = KSP_DIVERGED_DTOL
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_BREAKDOWN = KSP_DIVERGED_BREAKDOWN
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_BREAKDOWN_BICG = KSP_DIVERGED_BREAKDOWN_BICG
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_NONSYMMETRIC = KSP_DIVERGED_NONSYMMETRIC
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_INDEFINITE_PC = KSP_DIVERGED_INDEFINITE_PC
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_NANORINF = KSP_DIVERGED_NANORINF
  KSPConvergedReason, PARAMETER :: PETSC_KSP_DIVERGED_INDEFINITE_MAT = KSP_DIVERGED_INDEFINITE_MAT
  !>@}

  !KSPNorm types
  !> \addtogroup OpenCMISSPETSc_PetscKSPNormTypes OpenCMISS::PETSc::Constants::KSP::KSPNormType
  !> \brief Types of KSP norm parameters for PETSc KSP solvers
  !> \see OpenCMISSPETSc
  !>@{
  KSPNormType, PARAMETER :: PETSC_KSP_NORM_NONE = KSP_NORM_NONE
  KSPNormType, PARAMETER :: PETSC_KSP_NORM_PRECONDITIONED = KSP_NORM_PRECONDITIONED
  KSPNormType, PARAMETER :: PETSC_KSP_NORM_UNPRECONDITIONED = KSP_NORM_UNPRECONDITIONED
  KSPNormType, PARAMETER :: PETSC_KSP_NORM_NATURAL = KSP_NORM_NATURAL
  !>@}
  !>@}

  !> \addtogroup OpenCMISSPETSc_PetscMatConstants OpenCMISS::PETSc::Constants::Mat
  !> \brief Constants for PETSc matrices
  !> \see OpenCMISSPETSc
  !>@{
  !MatAssembly types
  !> \addtogroup OpenCMISSPETSc_PetscMatAssemblyTypes OpenCMISS::PETSc::Constants::Mat::MatAssemblyType
  !> \brief Types of Matrix assembly parameters for PETSc matrices
  !> \see OpenCMISSPETSc
  !>@{
  MatAssemblyType, PARAMETER :: PETSC_MAT_FLUSH_ASSEMBLY = MAT_FLUSH_ASSEMBLY !<Continue to put values into the matrix
  MatAssemblyType, PARAMETER :: PETSC_MAT_FINAL_ASSEMBLY = MAT_FINAL_ASSEMBLY !<Wish to use the matrix with the values currently in the matrix
  !>@}

  !MatDuplicate types
  !> \addtogroup OpenCMISSPETSc_PetscMatDuplicateOptions OpenCMISS::PETSc::Constants::Mat::MatDuplicateOption
  !> \brief Types of Matrix duplicate options
  !> \see OpenCMISSPETSc
  !>@{
  MatDuplicateOption, PARAMETER :: PETSC_MAT_DO_NOT_COPY_VALUES = MAT_DO_NOT_COPY_VALUES !<Create matrix using the same nonzero pattern but use zero for the numerical values
  MatDuplicateOption, PARAMETER :: PETSC_MAT_COPY_VALUES = MAT_COPY_VALUES !<Create the matrix with the same nonzero pattern and the same numerical values
  MatDuplicateOption, PARAMETER :: PETSC_MAT_SHARE_NONZERO_PATTERN = MAT_SHARE_NONZERO_PATTERN !<Create matrix with the same non-zero pattern but do not copy the values.
  !>@}

  !MatFactor types,
  !> \addtogroup OpenCMISSPETSc_PetscMatFactorTypes OpenCMISS::PETSc::Constants::Mat::MatFactorType
  !> \brief Types of Matrix factor parameters for PETSc matrices
  !> \see OpenCMISSPETSc
  !>@{
  MatFactorType, PARAMETER :: PETSC_MAT_FACTOR_NONE = MAT_FACTOR_NONE !<No matrix factor 
  MatFactorType, PARAMETER :: PETSC_MAT_FACTOR_LU = MAT_FACTOR_LU !<LU factorisation
  MatFactorType, PARAMETER :: PETSC_MAT_FACTOR_CHOLESKY = MAT_FACTOR_CHOLESKY !<Cholesky factorisation
  MatFactorType, PARAMETER :: PETSC_MAT_FACTOR_ILU = MAT_FACTOR_ILU !<ILU factorisation
  MatFactorType, PARAMETER :: PETSC_MAT_FACTOR_ICC = MAT_FACTOR_ICC !<Incomplete Cholesky factorisation
  MatFactorType, PARAMETER :: PETSC_MAT_FACTOR_ILUDT = MAT_FACTOR_ILUDT !<ILU factorisation with drop tolerances
  MatFactorType, PARAMETER :: PETSC_MAT_FACTOR_QR = MAT_FACTOR_QR !<QR factorisation
  !>@}

  !MatInfo types
  !> \addtogroup OpenCMISSPETSc_PetscMatInfos OpenCMISS::PETSc::Constants::Mat::MatInfo
  !> \brief Matrix info parameters for PETSc matrices
  !> \see OpenCMISSPETSc
  !>@{
  MatInfo, PARAMETER :: PETSC_MAT_INFO_SIZE = MAT_INFO_SIZE
  MatInfo, PARAMETER :: PETSC_MAT_INFO_BLOCK_SIZE = MAT_INFO_BLOCK_SIZE
  MatInfo, PARAMETER :: PETSC_MAT_INFO_NZ_ALLOCATED = MAT_INFO_NZ_ALLOCATED
  MatInfo, PARAMETER :: PETSC_MAT_INFO_NZ_USED = MAT_INFO_NZ_USED
  MatInfo, PARAMETER :: PETSC_MAT_INFO_NZ_UNNEEDED = MAT_INFO_NZ_UNNEEDED
  MatInfo, PARAMETER :: PETSC_MAT_INFO_MEMORY = MAT_INFO_MEMORY
  MatInfo, PARAMETER :: PETSC_MAT_INFO_ASSEMBLIES = MAT_INFO_ASSEMBLIES
  MatInfo, PARAMETER :: PETSC_MAT_INFO_MALLOCS = MAT_INFO_MALLOCS
  MatInfo, PARAMETER :: PETSC_MAT_INFO_FILL_RATIO_GIVEN = MAT_INFO_FILL_RATIO_GIVEN
  MatInfo, PARAMETER :: PETSC_MAT_INFO_FILL_RATIO_NEEDED = MAT_INFO_FILL_RATIO_NEEDED
  MatInfo, PARAMETER :: PETSC_MAT_INFO_FACTOR_MALLOCS = MAT_INFO_FACTOR_MALLOCS  
  !>@}

  !MatInfoType types
  !> \addtogroup OpenCMISSPETSc_PetscMatInfoTypes OpenCMISS::PETSc::Constants::Mat::MatInfoType
  !> \brief Types of matrix info for PETSc matrices
  !> \see OpenCMISSPETSc
  !>@{
  MatInfoType, PARAMETER :: PETSC_MAT_LOCAL = MAT_LOCAL !<Values for each MPI process part of the matrix
  MatInfoType, PARAMETER :: PETSC_MAT_GLOBAL_MAX = MAT_GLOBAL_MAX !<Maximum of each value over all MPI processes
  MatInfoType, PARAMETER :: PETSC_MAT_GLOBAL_SUM = MAT_GLOBAL_SUM !<Sum of each value over all MPI processes
  !>@}
  
  !MatOption types
  !> \addtogroup OpenCMISSPETSc_PetscMatOptions OpenCMISS::PETSc::Constants::Mat::MatOption
  !> \brief Types of matrix options for PETSc matrices
  !> \see OpenCMISSPETSc
  !>@{
  MatOption, PARAMETER :: PETSC_MAT_SPD = MAT_SPD !<Matrix is symmetric and positive definite
  MatOption, PARAMETER :: PETSC_MAT_SYMMETRIC = MAT_SYMMETRIC !<Matrix is symmetric in terms of both structure and values.
  MatOption, PARAMETER :: PETSC_MAT_STRUCTURALLY_SYMMETRIC = MAT_STRUCTURALLY_SYMMETRIC !<Matrix has symmetric nonzero structure
  MatOption, PARAMETER :: PETSC_MAT_HERMITIAN = MAT_HERMITIAN !<Hermitian matrix, the transpose is the complex conjugation
  MatOption, PARAMETER :: PETSC_MAT_SYMMETRY_ETERNAL = MAT_SYMMETRY_ETERNAL !<Indicates the symmetry (or Hermitian structure) or its absense will persist through any changes to the matrix
  MatOption, PARAMETER :: PETSC_MAT_STRUCTURAL_SYMMETRY_ETERNAL = MAT_STRUCTURAL_SYMMETRY_ETERNAL !<Indicates the structural symmetry or its absence will persist through any changes to the matrix
  MatOption, PARAMETER :: PETSC_MAT_SPD_ETERNAL = MAT_SPD_ETERNAL !<Indicates the value of MAT_SPD (true or false) will persist through any changes to the matrix
  
  MatOption, PARAMETER :: PETSC_MAT_ROW_ORIENTED = MAT_ROW_ORIENTED !<Matrix will be stored in a row orientated fashion.
  
  MatOption, PARAMETER :: PETSC_MAT_NEW_NONZERO_LOCATIONS = MAT_NEW_NONZERO_LOCATIONS !<Any additions or insertions that would generate a new nonzero location are ignored.
  MatOption, PARAMETER :: PETSC_MAT_FORCE_DIAGONAL_ENTRIES = MAT_FORCE_DIAGONAL_ENTRIES
  MatOption, PARAMETER :: PETSC_MAT_IGNORE_OFF_PROC_ENTRIES = MAT_IGNORE_OFF_PROC_ENTRIES !<Any entries that are for other processors will be dropped
  MatOption, PARAMETER :: PETSC_MAT_NEW_NONZERO_LOCATION_ERR = MAT_NEW_NONZERO_LOCATION_ERR !<Any addition or insertion that will generate a new nonzero location will produce an error
  MatOption, PARAMETER :: PETSC_MAT_USE_HASH_TABLE = MAT_USE_HASH_TABLE !<Use a hash table for matrix assembly in order to improve matrix searches
  MatOption, PARAMETER :: PETSC_MAT_NO_OFF_PROC_ENTRIES = MAT_NO_OFF_PROC_ENTRIES !<Each process will only set values for its own rows.
  MatOption, PARAMETER :: PETSC_SUBSET_OFF_PROC_ENTRIES = MAT_SUBSET_OFF_PROC_ENTRIES !<you know that the first assembly after setting this flag will set a superset of the off-process entries required for all subsequent assemblies.
  
  MatOption, PARAMETER :: PETSC_MAT_NEW_NONZERO_ALLOCATION_ERR = MAT_NEW_NONZERO_ALLOCATION_ERR !<Any addition or insertion tat would generate a new entry that has not been preallocated will produce an error
  MatOption, PARAMETER :: PETSC_MAT_KEEP_NONZERO_PATTERN = MAT_KEEP_NONZERO_PATTERN !<When MatZeroRows is called the zeroed entries are kept in the nonzero strucutre
  MatOption, PARAMETER :: PETSC_MAT_IGNORE_ZERO_ENTRIES = MAT_IGNORE_ZERO_ENTRIES !<Stop zero values from creating a zero location in the matrix 
  MatOption, PARAMETER :: PETSC_MAT_USE_INODES = MAT_USE_INODES !<Matrix will using an inode version of code
  MatOption, PARAMETER :: PETSC_MAT_IGNORE_LOWER_TRIANGULAR = MAT_IGNORE_LOWER_TRIANGULAR !<Ignore any additions or insertions in the lower triangular part of the matrix
  MatOption, PARAMETER :: PETSC_MAT_ERROR_LOWER_TRIANGULAR = MAT_ERROR_LOWER_TRIANGULAR
  MatOption, PARAMETER :: PETSC_MAT_GETROW_UPPERTRIANGULAR = MAT_GETROW_UPPERTRIANGULAR
  MatOption, PARAMETER :: PETSC_MAT_UNUSED_NONZERO_LOCATION_ERR = MAT_UNUSED_NONZERO_LOCATION_ERR
  MatOption, PARAMETER :: PETSC_MAT_OPTION_MAX = MAT_OPTION_MAX
  MatOption, PARAMETER :: PETSC_MAT_OPTION_MIN = MAT_OPTION_MIN
  MatOption, PARAMETER :: PETSC_MAT_NO_OFF_PROC_ZERO_ROWS = MAT_NO_OFF_PROC_ZERO_ROWS
  !>@}
  
  !Matrix solver types
  !> \addtogroup OpenCMISSPETSc_PetscMatSolverTypes OpenCMISS::PETSc::Constants::Mat::MatSolverType
  !> \brief Types of matrix solver for PETSc matrices
  !> \see OpenCMISSPETSc
  !>@{
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_SUPERLU = MATSOLVERSUPERLU !<SuperLU matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_SUPERLU_DIST = MATSOLVERSUPERLU_DIST !<SuperLU_DIST matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_STRUMPACK = MATSOLVERSTRUMPACK !<Strumpack matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_UMFPACK = MATSOLVERUMFPACK !<UMFPack matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_CHOLMOD = MATSOLVERCHOLMOD !<CHOLmod matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_KLU = MATSOLVERKLU !<klu matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_ELEMENTAL = MATSOLVERELEMENTAL !<Elemental matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_SCALAPACK = MATSOLVERSCALAPACK !<ScaLAPACK matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_ESSL = MATSOLVERESSL !<ESSL matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_LUSOL = MATSOLVERLUSOL !LUSol matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_MUMPS = MATSOLVERMUMPS !<MUMPS matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_MKL_PARDISO = MATSOLVERMKL_PARDISO !<MKL Pardiso matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_MKL_CPARDISO = MATSOLVERMKL_CPARDISO !<MKL CPardiso matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_PASTIX = MATSOLVERPASTIX !PaStiX matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_MATLAB = MATSOLVERMATLAB !<Matlab matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_PETSC = MATSOLVERPETSC !PETSc matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_BAS = MATSOLVERBAS !BAS matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_CUSPARSE = MATSOLVERCUSPARSE !<CUDA Sparse matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_CUDA = MATSOLVERCUDA !<CUDA matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_HIPSPARSE = MATSOLVERHIPSPARSE !<HIP Sparse matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_HIP = MATSOLVERHIP !<HIP matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_KOKKOS = MATSOLVERKOKKOS !<KOKKOS matrix solver
  MatSolverType, PARAMETER :: PETSC_MAT_SOLVER_SPQR = MATSOLVERSPQR !<SPQR matrix solver
  !>@}
  
  !MatStructure types
  !> \addtogroup OpenCMISSPETSc_PetscMatStructures OpenCMISS::PETSc::Constants::Mat::MatStructure
  !> \brief The matrix structure parameters for PETSc matrices
  !> \see OpenCMISSPETSc
  !>@{
  MatStructure, PARAMETER :: PETSC_SAME_NONZERO_PATTERN = SAME_NONZERO_PATTERN !<Two matrices have identical nonzero patterns
  MatStructure, PARAMETER :: PETSC_DIFFERENT_NONZERO_PATTERN = DIFFERENT_NONZERO_PATTERN !<Two matrices may have different nonzero patterns
  MatStructure, PARAMETER :: PETSC_SUBSET_NONZERO_PATTERN = SUBSET_NONZERO_PATTERN !<The nonzero pattern of the second matrix is a subset of the first matrix
  MatStructure, PARAMETER :: PETSC_UNKNOWN_NONZERO_PATTERN = UNKNOWN_NONZERO_PATTERN !<There is no known relattionship between the nonzero patterns
  !>@}

  !MatReuse types
  !> \addtogroup OpenCMISSPETSc_PetscMatReuses OpenCMISS::PETSc::Constants::Mat::MatReuse
  !> \brief The matrix reuse parameters for PETSc matrices
  !> \see OpenCMISSPETSc
  !>@{
  MatReuse, PARAMETER :: PETSC_MAT_INITIAL_MATRIX = MAT_INITIAL_MATRIX !<Create a new matrix
  MatReuse, PARAMETER :: PETSC_MAT_REUSE_MATRIX = MAT_REUSE_MATRIX !<Reuse the matrix created with a previous call
  MatReuse, PARAMETER :: PETSC_MAT_IGNORE_MATRIX = MAT_IGNORE_MATRIX !<Replace the first matrix with the new matrix
  MatReuse, PARAMETER :: PETSC_MAT_INPLACE_MATRIX = MAT_INPLACE_MATRIX !<Do not create a new matrix or reuse a given matrix
  !>@}

  !MatColoring types
  !> \addtogroup OpenCMISSPETSc_PetscMatColoringTypes OpenCMISS::PETSc::Constants::Mat::MatColoringType
  !> \brief The matrix colouring types for PETSc matrices
  !> \see OpenCMISSPETSc
  !>@{
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_JP = MATCOLORINGJP
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_POWER = MATCOLORINGPOWER
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_NATURAL = MATCOLORINGNATURAL
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_SL = MATCOLORINGSL
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_LF = MATCOLORINGLF
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_ID = MATCOLORINGID
  MatColoringType, PARAMETER :: PETSC_MATCOLORING_GREEDY = MATCOLORINGGREEDY
  !>@}
  !>@}

  !Norm types
  !> \addtogroup OpenCMISSPETSc_PetscNormTypes OpenCMISS::PETSc::Constants::NormType
  !> \brief Types of norms for PETSc
  !> \see OpenCMISSPETSc
  !>@{
  NormType, PARAMETER :: PETSC_NORM_1 = NORM_1 !<One norm
  NormType, PARAMETER :: PETSC_NORM_2 = NORM_2 !<Two norm
  NormType, PARAMETER :: PETSC_NORM_FROBENIUS = NORM_FROBENIUS !<Frobenius norm
  NormType, PARAMETER :: PETSC_NORM_INFINITY = NORM_INFINITY !<Infinity norm
  NormType, PARAMETER :: PETSC_NORM_1_AND_2 = NORM_1_AND_2 !<Computes both the one and two norm
  !>@}
  
  !PC types
  !> \addtogroup OpenCMISSPETSc_PetscPCypes OpenCMISS::PETSc::Constants::PetscPCTypes
  !> \brief Types of preconditioners for PETSc KSP solvers.
  !> \see OpenCMISSPETSc
  !>@{
  PCType, PARAMETER ::  PETSC_PC_NONE = PCNONE !<No preconditioning
  PCType, PARAMETER ::  PETSC_PC_JACOBI = PCJACOBI !<Jacobi preconditioner
  PCType, PARAMETER ::  PETSC_PC_SOR = PCSOR !<SOR (and SSOR) preconditioner
  PCType, PARAMETER ::  PETSC_PC_LU = PCLU !<LU preconditioner
  PCType, PARAMETER ::  PETSC_PC_QR = PCQR
  PCType, PARAMETER ::  PETSC_PC_SHELL = PCSHELL !<Shell for user defined preconditioner
  PCType, PARAMETER ::  PETSC_PC_AMGX = PCAMGX
  PCType, PARAMETER ::  PETSC_PC_BJACOBI = PCBJACOBI !<Block Jacobi preconditioner
  PCType, PARAMETER ::  PETSC_PC_MG = PCMG
  PCType, PARAMETER ::  PETSC_PC_EISENSTAT = PCEISENSTAT !<SOR with Eisenstat trick
  PCType, PARAMETER ::  PETSC_PC_ILU = PCILU !<Incomplete LU preconditioner
  PCType, PARAMETER ::  PETSC_PC_ICC = PCICC !<Incomplete Cholesky preconditioner
  PCType, PARAMETER ::  PETSC_PC_ASM = PCASM !<Additive Schwarz preconditioner
  PCType, PARAMETER ::  PETSC_PC_GASM = PCGASM !<Generalised additive Schwarz
  PCType, PARAMETER ::  PETSC_PC_KSP = PCKSP !<Linear solver
  PCType, PARAMETER ::  PETSC_PC_COMPOSITE = PCCOMPOSITE !<Combination of preconditioners
  PCType, PARAMETER ::  PETSC_PC_REDUNDANT = PCREDUNDANT
  PCType, PARAMETER ::  PETSC_PC_SPAI = PCSPAI
  PCType, PARAMETER ::  PETSC_PC_NN = PCNN
  PCType, PARAMETER ::  PETSC_PC_CHOLESKY = PCCHOLESKY !<Cholesky preconditioner
  PCType, PARAMETER ::  PETSC_PC_PBJACOBI = PCPBJACOBI
  PCType, PARAMETER ::  PETSC_PC_VPBJACOBI = PCVPBJACOBI
  PCType, PARAMETER ::  PETSC_PC_MAT = PCMAT
  PCType, PARAMETER ::  PETSC_PC_HYPRE = PCHYPRE
  PCType, PARAMETER ::  PETSC_PC_PARMS = PCPARMS
  PCType, PARAMETER ::  PETSC_PC_FIELDSPLIT = PCFIELDSPLIT
  PCType, PARAMETER ::  PETSC_PC_TFS = PCTFS
  PCType, PARAMETER ::  PETSC_PC_ML = PCML
  PCType, PARAMETER ::  PETSC_PC_GALERKIN = PCGALERKIN
  PCType, PARAMETER ::  PETSC_PC_EXOTIC = PCEXOTIC
  PCType, PARAMETER ::  PETSC_PC_CP = PCCP
  PCType, PARAMETER ::  PETSC_PC_BFBT = PCBFBT
  PCType, PARAMETER ::  PETSC_PC_LSC = PCLSC
  PCType, PARAMETER ::  PETSC_PC_PYTHON = PCPYTHON
  PCType, PARAMETER ::  PETSC_PC_PFMG = PCPFMG
  PCType, PARAMETER ::  PETSC_PC_SYSPFMG = PCSYSPFMG
  PCType, PARAMETER ::  PETSC_PC_REDISTRIBUTE = PCREDISTRIBUTE
  PCType, PARAMETER ::  PETSC_PC_SVD = PCSVD
  PCType, PARAMETER ::  PETSC_PC_GAMG = PCGAMG
  PCType, PARAMETER ::  PETSC_PC_CHOWILUVIENNACL = PCCHOWILUVIENNACL
  PCType, PARAMETER ::  PETSC_PC_ROWSCALINGVIENNACL = PCROWSCALINGVIENNACL
  PCType, PARAMETER ::  PETSC_PC_SAVIENNACL = PCSAVIENNACL
  PCType, PARAMETER ::  PETSC_PC_BDDC = PCBDDC
  PCType, PARAMETER ::  PETSC_PC_KACZMARZ = PCKACZMARZ
  PCType, PARAMETER ::  PETSC_PC_TELESCOPE = PCTELESCOPE
  PCType, PARAMETER ::  PETSC_PC_PATCJ = PCPATCH
  PCType, PARAMETER ::  PETSC_PC_LMVM = PCLMVM
  PCType, PARAMETER ::  PETSC_PC_HMG = PCHMG
  PCType, PARAMETER ::  PETSC_PC_DEFLATION = PCDEFLATION
  PCType, PARAMETER ::  PETSC_PC_HPDDM = PCHPDDM
  PCType, PARAMETER ::  PETSC_PC_H2OPUS = PCH2OPUS
  PCType, PARAMETER ::  PETSC_PC_MPI = PCMPI
  !>@}

  !SNES types
  SNESType, PARAMETER :: PETSC_SNES_NEWTONLS = SNESNEWTONLS !<Newton's method
  SNESType, PARAMETER :: PETSC_SNES_NEWTONTR = SNESNEWTONTR !<Newton's method with trust region
  SNESType, PARAMETER :: PETSC_SNES_NRICHARDSON = SNESNRICHARDSON !<Nonlinear Richardson
  !SNESType, PARAMETER :: PETSC_SNES_ANDERSON = SNESANDERSON !<Anderson mixing
  SNESType, PARAMETER :: PETSC_SNES_KSPONLY = SNESKSPONLY !<Single linearisation
  SNESType, PARAMETER :: PETSC_SNES_VINEWTONRSLS = SNESVINEWTONRSLS
  SNESType, PARAMETER :: PETSC_SNES_VINEWTONSSLS = SNESVINEWTONSSLS
  SNESType, PARAMETER :: PETSC_SNES_NGMRES = SNESNGMRES !<Nonlinear GMRES
  SNESType, PARAMETER :: PETSC_SNES_QN = SNESQN !<Quasi-Newton (BFGS)
  SNESType, PARAMETER :: PETSC_SNES_SHELL = SNESSHELL !<Shell (user defined)
  SNESType, PARAMETER :: PETSC_SNES_NCG = SNESNCG !<Nonlinear CG
  SNESType, PARAMETER :: PETSC_SNES_FAS = SNESFAS !<Full Approximation Scheme (nonlinear multigrid)
  !SNESType, PARAMETER :: PETSC_SNES_NGS = SNESNGS !<Nonlinear Gauss-Siedel
  !SNESType, PARAMETER :: PETSC_SNES_NASM = SNESNASM !<Nonlinear additive Schwarz
  !SNESType, PARAMETER :: PETSC_SNES_ASPIN = SNESASPIN !<Nonlinear additive Schwarz preconditioned inexact Newton (ASPIN) method
  !SNESType, PARAMETER :: PETSC_SNES_COMPOSITE = SNESCOMPOSITE !<Composite (comine several nonlinear solvers)

  !SNES converged types
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_FNORM_ABS = SNES_CONVERGED_FNORM_ABS
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_FNORM_RELATIVE = SNES_CONVERGED_FNORM_RELATIVE
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_SNORM_RELATIVE = SNES_CONVERGED_SNORM_RELATIVE
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_ITS = SNES_CONVERGED_ITS
  !SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_INNER_ITER = SNES_CONVERGED_INNER_ITER
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_FUNCTION_DOMAIN = SNES_DIVERGED_FUNCTION_DOMAIN
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_FUNCTION_COUNT = SNES_DIVERGED_FUNCTION_COUNT
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_LINEAR_SOLVE = SNES_DIVERGED_LINEAR_SOLVE
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_FNORM_NAN = SNES_DIVERGED_FNORM_NAN
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_MAX_IT = SNES_DIVERGED_MAX_IT
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_LINE_SEARCH = SNES_DIVERGED_LINE_SEARCH
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_INNER = SNES_DIVERGED_INNER
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_LOCAL_MIN = SNES_DIVERGED_LOCAL_MIN
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_DTOL = SNES_DIVERGED_DTOL
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_JACOBIAN_DOMAIN = SNES_DIVERGED_JACOBIAN_DOMAIN
  SNESConvergedReason, PARAMETER :: PETSC_SNES_DIVERGED_TR_DELTA = SNES_DIVERGED_TR_DELTA
  SNESConvergedReason, PARAMETER :: PETSC_SNES_CONVERGED_ITERATING = SNES_CONVERGED_ITERATING
  
  !SNES line search type
  SNESLineSearchType, PARAMETER :: PETSC_SNES_LINESEARCH_BASIC = SNESLINESEARCHBASIC
  SNESLineSearchType, PARAMETER :: PETSC_SNES_LINESEARCH_BT = SNESLINESEARCHBT
  SNESLineSearchType, PARAMETER :: PETSC_SNES_LINESEARCH_L2 = SNESLINESEARCHL2
  SNESLineSearchType, PARAMETER :: PETSC_SNES_LINESEARCH_CP = SNESLINESEARCHCP
  SNESLineSearchType, PARAMETER :: PETSC_SNES_LINESEARCH_SHELL = SNESLINESEARCHSHELL
  SNESLineSearchType, PARAMETER :: PETSC_SNES_LINESEARCH_NCGLINEAR = SNESLINESEARCHNCGLINEAR
  
  !SNES line search order  
  SNESLineSearchOrder, PARAMETER :: PETSC_SNES_LINESEARCH_ORDER_LINEAR = SNES_LINESEARCH_ORDER_LINEAR
  SNESLineSearchOrder, PARAMETER :: PETSC_SNES_LINESEARCH_ORDER_QUADRATIC = SNES_LINESEARCH_ORDER_QUADRATIC
  SNESLineSearchOrder, PARAMETER :: PETSC_SNES_LINESEARCH_ORDER_CUBIC = SNES_LINESEARCH_ORDER_CUBIC

  !SNES multi-stage smoother type 
  SNESMSType, PARAMETER :: PETSC_SNESMS_M62 = SNESMSM62
  SNESMSType, PARAMETER :: PETSC_SNESMS_EULER = SNESMSEULER
  SNESMSType, PARAMETER :: PETSC_SNESMS_JAMESON83 = SNESMSJAMESON83
  SNESMSType, PARAMETER :: PETSC_SNESMS_VLTP11 = SNESMSVLTP11
  SNESMSType, PARAMETER :: PETSC_SNESMS_VLTP21 = SNESMSVLTP21
  SNESMSType, PARAMETER :: PETSC_SNESMS_VLTP31 = SNESMSVLTP31
  SNESMSType, PARAMETER :: PETSC_SNESMS_VLTP41 = SNESMSVLTP41
  SNESMSType, PARAMETER :: PETSC_SNESMS_VLTP51 = SNESMSVLTP51
  SNESMSType, PARAMETER :: PETSC_SNESMS_VLTP61 = SNESMSVLTP61

  !SNES norm schedules
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_DEFAULT = SNES_NORM_DEFAULT
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_NONE = SNES_NORM_NONE
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_ALWAYS = SNES_NORM_ALWAYS
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_INITIAL_ONLY = SNES_NORM_INITIAL_ONLY
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_FINAL_ONLY = SNES_NORM_FINAL_ONLY
  SNESNormSchedule, PARAMETER :: PETSC_SNES_NORM_INITIAL_FINAL_ONLY = SNES_NORM_INITIAL_FINAL_ONLY

  !SNES QN types
  SNESQNType, PARAMETER :: PETSC_SNES_QN_LBFGS = SNES_QN_LBFGS 
  SNESQNType, PARAMETER :: PETSC_SNES_QN_BROYDEN = SNES_QN_BROYDEN
  SNESQNType, PARAMETER :: PETSC_SNES_QN_BADBROYDEN = SNES_QN_BADBROYDEN
  
  !SNES QN restart types
  SNESQNRestartType, PARAMETER :: PETSC_SNES_QN_RESTART_DEFAULT = SNES_QN_RESTART_DEFAULT
  SNESQNRestartType, PARAMETER :: PETSC_SNES_QN_RESTART_NONE = SNES_QN_RESTART_NONE 
  SNESQNRestartType, PARAMETER :: PETSC_SNES_QN_RESTART_POWELL = SNES_QN_RESTART_POWELL
  SNESQNRestartType, PARAMETER :: PETSC_SNES_QN_RESTART_PERIODIC = SNES_QN_RESTART_PERIODIC 

  !SNES QN scaling types
  SNESQNScaleType, PARAMETER :: PETSC_SNES_QN_SCALE_DEFAULT = SNES_QN_SCALE_DEFAULT
  SNESQNScaleType, PARAMETER :: PETSC_SNES_QN_SCALE_NONE = SNES_QN_SCALE_NONE 
  SNESQNScaleType, PARAMETER :: PETSC_SNES_QN_SCALE_SHANNO = SNES_QN_SCALE_SHANNO 
  SNESQNScaleType, PARAMETER :: PETSC_SNES_QN_SCALE_LINESEARCH = SNES_QN_SCALE_LINESEARCH  
  SNESQNScaleType, PARAMETER :: PETSC_SNES_QN_SCALE_JACOBIAN = SNES_QN_SCALE_JACOBIAN

!Temporary until we have a bug fix in PETSc  
!#define TaoType character*(80)
  
  !TAO types
  TaoType, PARAMETER :: PETSC_TAO_NLS = TAONLS !Newton's method with linesearch
  TaoType, PARAMETER :: PETSC_TAO_NTR = TAONTR !Newton's method with trust region
  TaoType, PARAMETER :: PETSC_TAO_NTL = TAONTL !Newton's method with trust region, lineasearch
  TaoType, PARAMETER :: PETSC_TAO_LMVM = TAOLMVM !Limited Memory Variable Metric (quasi-Newton method)
  TaoType, PARAMETER :: PETSC_TAO_CG = TAOCG !Nonlinear conjugate gradient
  TaoType, PARAMETER :: PETSC_TAO_NM = TAONM !Nelder-Mead solver
  TaoType, PARAMETER :: PETSC_TAO_TRON = TAOTRON !TRON active-set Newton trust region algorithm
  TaoType, PARAMETER :: PETSC_TAO_GPCG = TAOGPCG !Gradient projected conjugate gradient algorithm
  TaoType, PARAMETER :: PETSC_TAO_BLMVM = TAOBLMVM !Bounded limited memory variable metric (quasi-Newton method)
  TaoType, PARAMETER :: PETSC_TAO_OWLQN = TAOOWLQN !Orthant-wise Limited Memory (quasi-Newton method)
  TaoType, PARAMETER :: PETSC_TAO_LCL = TAOLCL !Linearly constrained Lagrangian method for PDE-constrained optimisation
  TaoType, PARAMETER :: PETSC_TAO_POUNDERS = TAOPOUNDERS !POUNDERS algorithm for nonlinear least squares
  TaoType, PARAMETER :: PETSC_TAO_BMRM = TAOBMRM !Bundle method for regularized risk minimisation
  TaoType, PARAMETER :: PETSC_TAO_BQPIP = TAOBQPIP !Bounded quadratic interior point alofrith for quadratic optimisation
  TaoType, PARAMETER :: PETSC_TAO_SSILS = TAOSSILS !Semi-smooth infeasible linesearch algorithm
  TaoType, PARAMETER :: PETSC_TAO_SSFLS = TAOSSFLS !Semi-smooth feasible linesearch algorithm
  TaoType, PARAMETER :: PETSC_TAO_ASILS = TAOASILS !Active-set infeasible linesearch algorithm
  TaoType, PARAMETER :: PETSC_TAO_ASFLS = TAOASFLS !Active-set feasible linesearch algorithm
  TaoType, PARAMETER :: PETSC_TAO_IPM = TAOIPM !Interior point algorithm

  !TAO converged reasons
  TaoConvergedReason, PARAMETER :: PETSC_TAO_CONVERGED_GATOL = TAO_CONVERGED_GATOL !<||g(X)|| < gatol
  TaoConvergedReason, PARAMETER :: PETSC_TAO_CONVERGED_GRTOL = TAO_CONVERGED_GRTOL !<||g(X)|| / f(X)  < grtol 
  TaoConvergedReason, PARAMETER :: PETSC_TAO_CONVERGED_GTTOL = TAO_CONVERGED_GTTOL !<||g(X)|| / ||g(X0)|| < gttol 
  TaoConvergedReason, PARAMETER :: PETSC_TAO_CONVERGED_STEPTOL = TAO_CONVERGED_STEPTOL  !<Step size small
  TaoConvergedReason, PARAMETER :: PETSC_TAO_CONVERGED_MINF = TAO_CONVERGED_MINF !<F < F_min
  TaoConvergedReason, PARAMETER :: PETSC_TAO_CONVERGED_USER = TAO_CONVERGED_USER !<User defined
  TaoConvergedReason, PARAMETER :: PETSC_TAO_DIVERGED_MAXITS = TAO_DIVERGED_MAXITS !<Maximum number of iterations reached
  TaoConvergedReason, PARAMETER :: PETSC_TAO_DIVERGED_NAN = TAO_DIVERGED_NAN !<NaN appeared in the computations
  TaoConvergedReason, PARAMETER :: PETSC_TAO_DIVERGED_MAXFCN = TAO_DIVERGED_MAXFCN !<Maximum number of function evaluations reached
  TaoConvergedReason, PARAMETER :: PETSC_TAO_DIVERGED_LS_FAILURE = TAO_DIVERGED_LS_FAILURE !<Line-search failure
  TaoConvergedReason, PARAMETER :: PETSC_TAO_DIVERGED_TR_REDUCTION = TAO_DIVERGED_TR_REDUCTION !<Trust-region failure
  TaoConvergedReason, PARAMETER :: PETSC_TAO_DIVERGED_USER = TAO_DIVERGED_USER !<User defined
  TaoConvergedReason, PARAMETER :: PETSC_TAO_CONTINUE_ITERATING = TAO_CONTINUE_ITERATING !<The optimisation is still running

  !TS constants
  !> \addtogroup OpenCMISSPETSc_PetscTSConstants OpenCMISS::PETSc::Constants::TS
  !> \brief Constants of PETSc time stepping (TS) objects
  !> \see OpenCMISSPETSc
  !>@{
  
  !TS exact final time options
  !> \addtogroup OpenCMISSPETSc_PetscTSExactFinalTimeOptions OpenCMISS::PETSc::Constants::TS::TSExactFinalTimeOptions
  !> \brief Constants for the handling of the final time step for PETSc TS objects.
  !> \see OpenCMISSPETSc
  !>@{
  TSExactFinalTimeOption, PARAMETER :: PETSC_TS_EXACTFINALTIME_STEPOVER = TS_EXACTFINALTIME_STEPOVER !<Don't do anything if requested final time is exeeded
  TSExactFinalTimeOption, PARAMETER :: PETSC_TS_EXACTFINALTIME_INTERPOLATE = TS_EXACTFINALTIME_INTERPOLATE !<Interpolate back to the final time
  TSExactFinalTimeOption, PARAMETER :: PETSC_TS_EXACTFINALTIME_MATCHSTEP = TS_EXACTFINALTIME_MATCHSTEP !<Adapt the final time step to match the final time requested 
  !>@}
  
  !TS types
  !> \addtogroup OpenCMISSPETSc_PetscTSTypes OpenCMISS::PETSc::Constants::TS::TSType
  !> \brief The types of PETSc time stepping (TS) methods
  !> \see OpenCMISSPETSc
  !>@{
  TSType, PARAMETER :: PETSC_TS_EULER = TSEULER !<Forward Euler time integrator
  TSType, PARAMETER :: PETSC_TS_BEULER = TSBEULER !<Backward Euler time integrator
  TSType, PARAMETER :: PETSC_TS_PSEUDO = TSPSEUDO !<
  TSType, PARAMETER :: PETSC_TS_CN = TSCN !<Crank-Nicolson time integrator
  TSType, PARAMETER :: PETSC_TS_SUNDIALS = TSSUNDIALS
  TSType, PARAMETER :: PETSC_TS_RK = TSRK !<Runge-Kutta time integrator
  TSType, PARAMETER :: PETSC_TS_PYTHON = TSPYTHON
  TSType, PARAMETER :: PETSC_TS_THETA = TSTHETA !<Theta based time integrator
  TSType, PARAMETER :: PETSC_TS_ALPHA = TSALPHA !<Alpha based time integrator
  TSType, PARAMETER :: PETSC_TS_GLLE = TSGLLE
  TSType, PARAMETER :: PETSC_TS_SSP = TSSSP !<Multistage SSP Runge-Kutta time integrator
  TSType, PARAMETER :: PETSC_TS_ARKIMEX = TSARKIMEX !<IMEX Runge-Kutta time integrator
  TSType, PARAMETER :: PETSC_TS_ROSW = TSROSW !<Rosenbrock-W time integrator
  TSType, PARAMETER :: PETSC_TS_EIMEX = TSEIMEX !<Extrapolated IMEX time integrator
  TSType, PARAMETER :: PETSC_TS_RADAU5 = TSRADAU5
  TSType, PARAMETER :: PETSC_TS_MPRK = TSMPRK !<Multirate Partitioned Runge-Kutta time inegrator
  !>@}

  !TS convergence flags
  !> \addtogroup OpenCMISSPETSc_PetscTSConvergedReasons OpenCMISS::PETSc::Constants::TS::TSConvergedReason
  !> \brief The types of PETSc time stepping (TS) methods
  !> \see OpenCMISSPETSc
  !>@{
  TSConvergedReason, PARAMETER :: PETSC_TS_CONVERGED_ITERATING = TS_CONVERGED_ITERATING
  TSConvergedReason, PARAMETER :: PETSC_TS_CONVERGED_TIME = TS_CONVERGED_TIME
  TSConvergedReason, PARAMETER :: PETSC_TS_CONVERGED_ITS = TS_CONVERGED_ITS
  TSConvergedReason, PARAMETER :: PETSC_TS_DIVERGED_NONLINEAR_SOLVE = TS_DIVERGED_NONLINEAR_SOLVE
  TSConvergedReason, PARAMETER :: PETSC_TS_DIVERGED_STEP_REJECTED = TS_DIVERGED_STEP_REJECTED
  !>@}
  
  !TS problem types
  !> \addtogroup OpenCMISSPETSc_PetscTSProblemTypes OpenCMISS::PETSc::Constants::TS::TSProblemType
  !> \brief The types of PETSc time stepping (TS) methods
  !> \see OpenCMISSPETSc
  !>@{
  TSProblemType, PARAMETER :: PETSC_TS_LINEAR = TS_LINEAR
  TSProblemType, PARAMETER :: PETSC_TS_NONLINEAR = TS_NONLINEAR
  !>@}
  
  !TS Sundials types
  !> \addtogroup OpenCMISSPETSc_PetscTSSundialsTypes OpenCMISS::PETSc::Constants::TS::TSSundialsType
  !> \brief The types of PETSc time stepping (TS) methods
  !> \see OpenCMISSPETSc
  !>@{
  TSSundialsType, PARAMETER :: PETSC_SUNDIALS_ADAMS = SUNDIALS_ADAMS
  TSSundialsType, PARAMETER :: PETSC_SUNDIALS_BDF = SUNDIALS_BDF
  !>@}

  !TS Sundials Gram Schmidt Type
  !> \addtogroup OpenCMISSPETSc_PetscTSSundialsGramSchmidtTypes OpenCMISS::PETSc::Constants::TS::TSSundialsGramSchmidtType
  !> \brief The types of PETSc time stepping (TS) methods
  !> \see OpenCMISSPETSc
  !>@{
  TSSundialsGramSchmidtType, PARAMETER :: PETSC_SUNDIALS_MODIFIED_GS = SUNDIALS_MODIFIED_GS
  TSSundialsGramSchmidtType, PARAMETER :: PETSC_SUNDIALS_CLASSICAL_GS = SUNDIALS_CLASSICAL_GS
  !>@}
  !>@}
  
  !Module types

  !Module variables

  LOGICAL, SAVE :: petscHandleError

  !Interfaces

  INTERFACE Petsc_SnesGetJacobian
    MODULE PROCEDURE Petsc_SnesGetJacobianSolver
    MODULE PROCEDURE Petsc_SnesGetJacobianSpecial
  END INTERFACE Petsc_SnesGetJacobian

  INTERFACE Petsc_SnesSetJacobian
    MODULE PROCEDURE Petsc_SnesSetJacobianSolver
  END INTERFACE Petsc_SnesSetJacobian

  !Miscelaneous routines and constants

  PUBLIC PETSC_TRUE,PETSC_FALSE
  
  PUBLIC PETSC_NULL_BOOL,PETSC_NULL_CHARACTER,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,PETSC_NULL_DOUBLE,PETSC_NULL_OBJECT, &
     & PETSC_NULL_SCALAR,PETSC_NULL_REAL

  PUBLIC PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_REAL

  PUBLIC PETSC_DECIDE

  PUBLIC PETSC_COMM_WORLD,PETSC_COMM_SELF

  PUBLIC PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_STDOUT_SELF,PETSC_VIEWER_DRAW_WORLD,PETSC_VIEWER_DRAW_SELF

  PUBLIC Petsc_Initialise,Petsc_Finalise
  
  PUBLIC Petsc_ErrorHandlingSetOn,Petsc_ErrorHandlingSetOff

  PUBLIC Petsc_LogView
  
  PUBLIC Petsc_OptionsSetValue
  
  !Value insert constants
  
  PUBLIC PETSC_NOT_SET_VALUES,PETSC_ADD_VALUES,PETSC_INSERT_VALUES,PETSC_MAX_VALUES,PETSC_MIN_VALUES,PETSC_INSERT_ALL_VALUES, &
    & PETSC_ADD_ALL_VALUES,PETSC_INSERT_BC_VALUES,PETSC_ADD_BC_VALUES
  
  PUBLIC PETSC_SCATTER_FORWARD,PETSC_SCATTER_REVERSE,PETSC_SCATTER_FORWARD_LOCAL,PETSC_SCATTER_REVERSE_LOCAL

  !IS routines

  PUBLIC Petsc_ISInitialise,Petsc_ISFinalise
  
  PUBLIC Petsc_ISDestroy
  
  !IS coloring routines

  PUBLIC Petsc_ISColoringInitialise,Petsc_ISColoringFinalise
  
  PUBLIC Petsc_ISColoringDestroy
  
  !IS local to global mapping routines

  PUBLIC Petsc_ISLocalToGlobalMappingInitialise,Petsc_ISLocalToGlobalMappingFinalise
  
  PUBLIC Petsc_ISLocalToGlobalMappingApply,Petsc_ISLocalToGlobalMappingApplyIS,Petsc_ISLocalToGlobalMappingCreate, &
    & Petsc_ISLocalToGlobalMappingDestroy

  !Norm type

  PUBLIC PETSC_NORM_1,PETSC_NORM_2,PETSC_NORM_FROBENIUS,PETSC_NORM_INFINITY,PETSC_NORM_1_AND_2
    
  !KSP routines and constants
 
  PUBLIC PETSC_KSP_RICHARDSON,PETSC_KSP_CHEBYSHEV,PETSC_KSP_CG,PETSC_KSP_GROPPCG,PETSC_KSP_PIPECG,PETSC_KSP_PIPECGRR, &
    & PETSC_KSP_PIPELCG,PETSC_KSP_PIPECG2,PETSC_KSP_CGNE,PETSC_KSP_NASH,PETSC_KSP_STCG,PETSC_KSP_GLTR,PETSC_KSP_FCG, &
    & PETSC_KSP_PIPEFCG,PETSC_KSP_GMRES,PETSC_KSP_PIPEFGMRES,PETSC_KSP_FGMRES,PETSC_KSP_LGMRES,PETSC_KSP_DGMRES, &
    & PETSC_KSP_PGMRES,PETSC_KSP_TCQMR,PETSC_KSP_BCGS,PETSC_KSP_IBCGS,PETSC_KSP_QMRCGS,PETSC_KSP_FBCGS,PETSC_KSP_FBCGSR, &
    & PETSC_KSP_BCGSL,PETSC_KSP_PIPEBCGS,PETSC_KSP_CGS,PETSC_KSP_LSQR,PETSC_KSP_PREONLY,PETSC_KSP_QCG,PETSC_KSP_BICG, &
    & PETSC_KSP_MINRES,PETSC_KSP_SYMMLQ,PETSC_KSP_LCD,PETSC_KSP_GCR,PETSC_KSP_PIPEGCR,PETSC_KSP_TSIRM,PETSC_KSP_CGLS, &
    & PETSC_KSP_FETIDP,PETSC_KSP_HPDDM

  PUBLIC PETSC_KSP_CONVERGED_RTOL,PETSC_KSP_CONVERGED_ATOL,PETSC_KSP_CONVERGED_ITS,PETSC_KSP_CONVERGED_NEG_CURVE, &
    & PETSC_KSP_CONVERGED_STEP_LENGTH,PETSC_KSP_CONVERGED_HAPPY_BREAKDOWN,PETSC_KSP_DIVERGED_NULL,PETSC_KSP_DIVERGED_ITS, &
    & PETSC_KSP_DIVERGED_DTOL,PETSC_KSP_DIVERGED_BREAKDOWN,PETSC_KSP_DIVERGED_BREAKDOWN_BICG,PETSC_KSP_DIVERGED_NONSYMMETRIC, &
    & PETSC_KSP_DIVERGED_INDEFINITE_PC,PETSC_KSP_DIVERGED_NANORINF,PETSC_KSP_DIVERGED_INDEFINITE_MAT

  PUBLIC PETSC_KSP_NORM_NONE,PETSC_KSP_NORM_PRECONDITIONED,PETSC_KSP_NORM_UNPRECONDITIONED,PETSC_KSP_NORM_NATURAL
  
  PUBLIC Petsc_KSPInitialise,Petsc_KSPFinalise
  
  PUBLIC Petsc_KSPCreate,Petsc_KSPDestroy,Petsc_KSPGetConvergedReason,Petsc_KSPGetIterationNumber,Petsc_KSPGetPC, &
    & Petsc_KSPGetResidualNorm,Petsc_KSPGMRESSetRestart,Petsc_KSPSetFromOptions,Petsc_KSPSetInitialGuessNonZero, &
    & Petsc_KSPSetOperators,Petsc_KSPSetReusePreconditioner,Petsc_KSPSetTolerances,Petsc_KSPSetType,Petsc_KSPSetUp, &
    & Petsc_KSPSolve
  
  !Matrix routines and constants

  PUBLIC PETSC_MAT_FLUSH_ASSEMBLY,PETSC_MAT_FINAL_ASSEMBLY

  PUBLIC PETSC_MAT_DO_NOT_COPY_VALUES,PETSC_MAT_COPY_VALUES,PETSC_MAT_SHARE_NONZERO_PATTERN

  PUBLIC PETSC_MAT_FACTOR_NONE,PETSC_MAT_FACTOR_LU,PETSC_MAT_FACTOR_CHOLESKY,PETSC_MAT_FACTOR_ILU,PETSC_MAT_FACTOR_ICC, &
    & PETSC_MAT_FACTOR_ILUDT,PETSC_MAT_FACTOR_QR
  
  PUBLIC PETSC_MAT_INFO_SIZE,PETSC_MAT_INFO_BLOCK_SIZE,PETSC_MAT_INFO_NZ_ALLOCATED,PETSC_MAT_INFO_NZ_USED, &
    & PETSC_MAT_INFO_NZ_UNNEEDED,PETSC_MAT_INFO_MEMORY,PETSC_MAT_INFO_ASSEMBLIES,PETSC_MAT_INFO_MALLOCS, &
    & PETSC_MAT_INFO_FILL_RATIO_GIVEN,PETSC_MAT_INFO_FILL_RATIO_NEEDED,PETSC_MAT_INFO_FACTOR_MALLOCS

  PUBLIC PETSC_MAT_LOCAL,PETSC_MAT_GLOBAL_MAX,PETSC_MAT_GLOBAL_SUM

  PUBLIC PETSC_MAT_SPD,PETSC_MAT_SYMMETRIC,PETSC_MAT_STRUCTURALLY_SYMMETRIC,PETSC_MAT_HERMITIAN,PETSC_MAT_SYMMETRY_ETERNAL, &
    & PETSC_MAT_STRUCTURAL_SYMMETRY_ETERNAL,PETSC_MAT_SPD_ETERNAL,PETSC_MAT_ROW_ORIENTED,PETSC_MAT_NEW_NONZERO_LOCATIONS, &
    & PETSC_MAT_FORCE_DIAGONAL_ENTRIES,PETSC_MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_MAT_NEW_NONZERO_LOCATION_ERR, &
    & PETSC_MAT_USE_HASH_TABLE,PETSC_MAT_NO_OFF_PROC_ENTRIES,PETSC_SUBSET_OFF_PROC_ENTRIES, &
    & PETSC_MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_MAT_KEEP_NONZERO_PATTERN,PETSC_MAT_IGNORE_ZERO_ENTRIES,PETSC_MAT_USE_INODES, &
    & PETSC_MAT_IGNORE_LOWER_TRIANGULAR,PETSC_MAT_ERROR_LOWER_TRIANGULAR,PETSC_MAT_GETROW_UPPERTRIANGULAR, &
    & PETSC_MAT_UNUSED_NONZERO_LOCATION_ERR,PETSC_MAT_OPTION_MAX,PETSC_MAT_OPTION_MIN,PETSC_MAT_NO_OFF_PROC_ZERO_ROWS
  
  PUBLIC PETSC_MAT_SOLVER_SUPERLU,PETSC_MAT_SOLVER_SUPERLU_DIST,PETSC_MAT_SOLVER_STRUMPACK,PETSC_MAT_SOLVER_UMFPACK, &
    & PETSC_MAT_SOLVER_CHOLMOD,PETSC_MAT_SOLVER_KLU,PETSC_MAT_SOLVER_ELEMENTAL,PETSC_MAT_SOLVER_SCALAPACK, &
    & PETSC_MAT_SOLVER_ESSL,PETSC_MAT_SOLVER_LUSOL,PETSC_MAT_SOLVER_MUMPS,PETSC_MAT_SOLVER_MKL_PARDISO, &
    & PETSC_MAT_SOLVER_MKL_CPARDISO,PETSC_MAT_SOLVER_PASTIX,PETSC_MAT_SOLVER_MATLAB,PETSC_MAT_SOLVER_PETSC, &
    & PETSC_MAT_SOLVER_BAS,PETSC_MAT_SOLVER_CUSPARSE,PETSC_MAT_SOLVER_CUDA,PETSC_MAT_SOLVER_HIPSPARSE,PETSC_MAT_SOLVER_HIP, &
    & PETSC_MAT_SOLVER_KOKKOS,PETSC_MAT_SOLVER_SPQR

  PUBLIC PETSC_SAME_NONZERO_PATTERN,PETSC_DIFFERENT_NONZERO_PATTERN,PETSC_SUBSET_NONZERO_PATTERN,PETSC_UNKNOWN_NONZERO_PATTERN

  PUBLIC PETSC_MAT_INITIAL_MATRIX,PETSC_MAT_REUSE_MATRIX,PETSC_MAT_IGNORE_MATRIX,PETSC_MAT_INPLACE_MATRIX
 
  PUBLIC Petsc_MatInitialise,Petsc_MatFinalise
  
  PUBLIC Petsc_MatAssemblyBegin,Petsc_MatAssemblyEnd,Petsc_MatCreate,Petsc_MatCreateAIJ,Petsc_MatCreateDense, &
    & Petsc_MatCreateSeqAIJ,Petsc_MatCreateSeqDense,Petsc_MatDenseGetArrayF90,Petsc_MatDenseRestoreArrayF90, &
    & Petsc_MatDestroy,Petsc_MatGetInfo,Petsc_MatGetOwnershipRange,Petsc_MatGetRow,Petsc_MatGetValues, &
    & Petsc_MatMumpsSetIcntl,Petsc_MatMumpsSetCntl,Petsc_MatRestoreRow,Petsc_MatSeqAIJGetArrayF90, &
    & Petsc_MatSeqAIJGetMaxRowNonzeros,Petsc_MatSeqAIJRestoreArrayF90,Petsc_MatSetLocalToGlobalMapping, &
    & Petsc_MatSetOption,Petsc_MatSetSizes,Petsc_MatSetValue,Petsc_MatSetValues,Petsc_MatSetValueLocal, &
    & Petsc_MatSetValuesLocal,Petsc_MatView,Petsc_MatZeroEntries
  
  !Matrix coloring routines and constants
  
  PUBLIC PETSC_MATCOLORING_JP,PETSC_MATCOLORING_POWER,PETSC_MATCOLORING_NATURAL,PETSC_MATCOLORING_SL,PETSC_MATCOLORING_LF, &
    & PETSC_MATCOLORING_ID,PETSC_MATCOLORING_GREEDY
  
  PUBLIC Petsc_MatColoringInitialise,Petsc_MatColoringFinalise
  
  PUBLIC Petsc_MatColoringApply,Petsc_MatColoringCreate,Petsc_MatColoringDestroy,Petsc_MatColoringSetFromOptions, &
    & Petsc_MatColoringSetType

  !Matrix FD coloring routines and constants

  PUBLIC Petsc_MatFDColoringInitialise,Petsc_MatFDColoringFinalise
  
  PUBLIC Petsc_MatFDColoringCreate,Petsc_MatFDColoringDestroy,Petsc_MatFDColoringSetFromOptions,Petsc_MatFDColoringSetFunction, &
    & Petsc_MatFDColoringSetParameters,Petsc_MatFDColoringSetup

  !Pre-conditioner routines and constants
  
  PUBLIC PETSC_PC_NONE,PETSC_PC_JACOBI,PETSC_PC_SOR,PETSC_PC_LU,PETSC_PC_QR,PETSC_PC_SHELL,PETSC_PC_AMGX,PETSC_PC_BJACOBI, &
    & PETSC_PC_MG,PETSC_PC_EISENSTAT,PETSC_PC_ILU,PETSC_PC_ICC,PETSC_PC_ASM,PETSC_PC_KSP,PETSC_PC_COMPOSITE, &
    & PETSC_PC_REDUNDANT,PETSC_PC_SPAI,PETSC_PC_NN,PETSC_PC_CHOLESKY,PETSC_PC_PBJACOBI,PETSC_PC_VPBJACOBI,PETSC_PC_MAT, &
    & PETSC_PC_HYPRE,PETSC_PC_PARMS,PETSC_PC_FIELDSPLIT,PETSC_PC_TFS,PETSC_PC_ML,PETSC_PC_GALERKIN,PETSC_PC_EXOTIC, &
    & PETSC_PC_CP,PETSC_PC_BFBT,PETSC_PC_LSC,PETSC_PC_PYTHON,PETSC_PC_PFMG,PETSC_PC_SYSPFMG,PETSC_PC_REDISTRIBUTE, &
    & PETSC_PC_SVD,PETSC_PC_GAMG,PETSC_PC_GASM,PETSC_PC_CHOWILUVIENNACL,PETSC_PC_ROWSCALINGVIENNACL,PETSC_PC_SAVIENNACL, &
    & PETSC_PC_BDDC,PETSC_PC_KACZMARZ,PETSC_PC_TELESCOPE,PETSC_PC_PATCJ,PETSC_PC_LMVM,PETSC_PC_HMG,PETSC_PC_DEFLATION, &
    & PETSC_PC_HPDDM,PETSC_PC_H2OPUS,PETSC_PC_MPI
  
  PUBLIC Petsc_PCInitialise,Petsc_PCFinalise
  
  PUBLIC Petsc_PCFactorGetMatrix,Petsc_PCFactorSetMatSolverType,Petsc_PCFactorSetUpMatSolverType,Petsc_PCSetFromOptions, &
    & Petsc_PCSetReusePreconditioner,Petsc_PCSetType

  !SNES routines and constants
  
  PUBLIC PETSC_SNES_NEWTONLS,PETSC_SNES_NEWTONTR,PETSC_SNES_NRICHARDSON,PETSC_SNES_KSPONLY,PETSC_SNES_VINEWTONRSLS, &
    & PETSC_SNES_VINEWTONSSLS,PETSC_SNES_NGMRES,PETSC_SNES_QN,PETSC_SNES_SHELL,PETSC_SNES_NCG,PETSC_SNES_FAS
  
  PUBLIC PETSC_SNES_CONVERGED_FNORM_ABS,PETSC_SNES_CONVERGED_FNORM_RELATIVE,PETSC_SNES_CONVERGED_SNORM_RELATIVE, &
    & PETSC_SNES_CONVERGED_ITS,PETSC_SNES_DIVERGED_FUNCTION_DOMAIN, &
    & PETSC_SNES_DIVERGED_FUNCTION_COUNT,PETSC_SNES_DIVERGED_LINEAR_SOLVE,PETSC_SNES_DIVERGED_FNORM_NAN, &
    & PETSC_SNES_DIVERGED_MAX_IT,PETSC_SNES_DIVERGED_LINE_SEARCH,PETSC_SNES_DIVERGED_LOCAL_MIN,PETSC_SNES_CONVERGED_ITERATING

  PUBLIC PETSC_SNES_NORM_DEFAULT,PETSC_SNES_NORM_NONE,PETSC_SNES_NORM_ALWAYS,PETSC_SNES_NORM_INITIAL_ONLY, &
    & PETSC_SNES_NORM_FINAL_ONLY,PETSC_SNES_NORM_INITIAL_FINAL_ONLY

  PUBLIC PETSC_SNES_QN_LBFGS,PETSC_SNES_QN_BROYDEN,PETSC_SNES_QN_BADBROYDEN
  
  PUBLIC PETSC_SNES_QN_RESTART_NONE,PETSC_SNES_QN_RESTART_POWELL,PETSC_SNES_QN_RESTART_PERIODIC
  
  PUBLIC PETSC_SNES_QN_SCALE_DEFAULT,PETSC_SNES_QN_SCALE_NONE,PETSC_SNES_QN_SCALE_SHANNO,PETSC_SNES_QN_SCALE_LINESEARCH, &
    & PETSC_SNES_QN_SCALE_JACOBIAN

  PUBLIC Petsc_SnesInitialise,Petsc_SnesFinalise
  
  PUBLIC Petsc_SnesComputeJacobianDefault,Petsc_SnesComputeJacobianDefaultColor
  
  PUBLIC Petsc_SnesCreate,Petsc_SnesDestroy,Petsc_SnesGetApplicationContext,Petsc_SnesGetConvergedReason,Petsc_SnesGetFunction, &
    & Petsc_SnesGetIterationNumber,Petsc_SnesGetJacobian,Petsc_SnesGetKSP,Petsc_SnesGetLineSearch,Petsc_SnesGetSolutionUpdate, &
    & Petsc_SnesMonitorSet,Petsc_SnesQNSetRestartType,Petsc_SnesQNSetScaleType,Petsc_SnesQNSetType, &
    & Petsc_SnesSetApplicationContext,Petsc_SnesSetConvergenceTest,Petsc_SnesSetFromOptions,Petsc_SnesSetFunction, &
    & Petsc_SnesSetJacobian,Petsc_SnesSetKSP,Petsc_SnesSetNormSchedule,Petsc_SnesSetTolerances,Petsc_SnesSetTrustRegionTolerance, &
    & Petsc_SnesSetType,Petsc_SnesSolve
  
  !SNES line search routines and constants

  PUBLIC PETSC_SNES_LINESEARCH_BASIC,PETSC_SNES_LINESEARCH_BT,PETSC_SNES_LINESEARCH_L2,PETSC_SNES_LINESEARCH_CP, &
    & PETSC_SNES_LINESEARCH_SHELL
  
  PUBLIC PETSC_SNES_LINESEARCH_ORDER_LINEAR,PETSC_SNES_LINESEARCH_ORDER_QUADRATIC,PETSC_SNES_LINESEARCH_ORDER_CUBIC

  PUBLIC Petsc_SnesLineSearchInitialise,Petsc_SnesLineSearchFinalise

  PUBLIC Petsc_SnesLineSearchBTSetAlpha,Petsc_SnesLineSearchComputeNorms,Petsc_SnesLineSearchGetNorms,Petsc_SnesLineSearchGetVecs, &
    & Petsc_SnesLineSearchMonitorSet,Petsc_SnesLineSearchSetComputeNorms,Petsc_SnesLineSearchSetNorms, &
    & Petsc_SnesLineSearchSetOrder,Petsc_SnesLineSearchSetTolerances,Petsc_SnesLineSearchSetType
  
  !TAO routines and constants
  
  PUBLIC PETSC_TAO_LMVM,PETSC_TAO_NLS,PETSC_TAO_NTR,PETSC_TAO_NTL,PETSC_TAO_CG,PETSC_TAO_TRON,PETSC_TAO_OWLQN,PETSC_TAO_BMRM, &
    & PETSC_TAO_BLMVM,PETSC_TAO_BQPIP,PETSC_TAO_GPCG,PETSC_TAO_NM,PETSC_TAO_POUNDERS,PETSC_TAO_LCL,PETSC_TAO_SSILS, &
    & PETSC_TAO_SSFLS,PETSC_TAO_ASILS,PETSC_TAO_ASFLS,PETSC_TAO_IPM
  
  PUBLIC PETSC_TAO_CONVERGED_GATOL,PETSC_TAO_CONVERGED_GRTOL,PETSC_TAO_CONVERGED_GTTOL,PETSC_TAO_CONVERGED_STEPTOL, &
    & PETSC_TAO_CONVERGED_MINF,PETSC_TAO_CONVERGED_USER,PETSC_TAO_DIVERGED_MAXITS,PETSC_TAO_DIVERGED_MAXFCN, &
    & PETSC_TAO_DIVERGED_LS_FAILURE,PETSC_TAO_DIVERGED_TR_REDUCTION,PETSC_TAO_DIVERGED_USER,PETSC_TAO_CONTINUE_ITERATING
  
  !TS routines and constants

  PUBLIC PETSC_TS_EXACTFINALTIME_STEPOVER,PETSC_TS_EXACTFINALTIME_INTERPOLATE,PETSC_TS_EXACTFINALTIME_MATCHSTEP
  
  PUBLIC PETSC_TS_EULER,PETSC_TS_BEULER,PETSC_TS_PSEUDO,PETSC_TS_CN,PETSC_TS_SUNDIALS,PETSC_TS_RK,PETSC_TS_PYTHON, &
    & PETSC_TS_THETA,PETSC_TS_ALPHA,PETSC_TS_GLLE,PETSC_TS_SSP,PETSC_TS_ARKIMEX,PETSC_TS_ROSW,PETSC_TS_EIMEX

  PUBLIC PETSC_TS_CONVERGED_ITERATING,PETSC_TS_CONVERGED_TIME,PETSC_TS_CONVERGED_ITS,PETSC_TS_DIVERGED_NONLINEAR_SOLVE, &
    & PETSC_TS_DIVERGED_STEP_REJECTED

  PUBLIC PETSC_TS_LINEAR,PETSC_TS_NONLINEAR

  PUBLIC PETSC_SUNDIALS_ADAMS,PETSC_SUNDIALS_BDF

  PUBLIC PETSC_SUNDIALS_MODIFIED_GS,PETSC_SUNDIALS_CLASSICAL_GS

  PUBLIC Petsc_TSInitialise,Petsc_TSFinalise

  PUBLIC Petsc_TSCreate,Petsc_TSDestroy,Petsc_TSGetSolution,Petsc_TSMonitorSet,Petsc_TSSetDuration,Petsc_TSSetExactFinalTime, &
    & Petsc_TSSetFromOptions,Petsc_TSSetMaxTime,Petsc_TSSetProblemType,Petsc_TSSetRHSFunction, &
    & Petsc_TSSetSolution,Petsc_TSSetTime,Petsc_TSSetTimeStep,Petsc_TSSetType,Petsc_TSSolve,Petsc_TSStep, &
    & Petsc_TSSundialsSetTolerance,Petsc_TSSundialsSetType
  
  !Vector routines and constants

  PUBLIC Petsc_VecInitialise,Petsc_VecFinalise
  
  PUBLIC Petsc_VecAssemblyBegin,Petsc_VecAssemblyEnd,Petsc_VecCopy,Petsc_VecCreate,Petsc_VecCreateGhost, &
    & Petsc_VecCreateGhostWithArray,Petsc_VecCreateMPI,Petsc_VecCreateMPIWithArray,Petsc_VecCreateSeq, &
    & Petsc_VecCreateSeqWithArray,Petsc_VecDestroy,Petsc_VecDuplicate,Petsc_VecDot,Petsc_VecGetArrayF90, &
    & Petsc_VecGetArrayReadF90,Petsc_VecGetLocalSize,Petsc_VecGetOwnershipRange,Petsc_VecGetSize,Petsc_VecGetValues, &
    & Petsc_VecNorm,Petsc_VecRestoreArrayF90,Petsc_VecRestoreArrayReadF90,Petsc_VecScale,Petsc_VecSet,Petsc_VecSetFromOptions, &
    & Petsc_VecSetLocalToGlobalMapping,Petsc_VecSetSizes,Petsc_VecSetValues,Petsc_VecView
  
CONTAINS
  
  !
  !================================================================================================================================
  !
  
  !>Set PETSc error handling on
  SUBROUTINE Petsc_ErrorHandlingSetOff(err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ErrorHandlingSetOff",err,error,*999)

    petscHandleError=.FALSE.
    
    EXITS("Petsc_ErrorHandlingSetOff")
    RETURN
999 ERRORSEXITS("Petsc_ErrorHandlingSetOff",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ErrorHandlingSetOff
    
  !
  !================================================================================================================================
  !

  !>Set PETSc error handling on
  SUBROUTINE Petsc_ErrorHandlingSetOn(err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ErrorHandlingSetOn",err,error,*999)

    petscHandleError=.TRUE.
    
    EXITS("Petsc_ErrorHandlingSetOn")
    RETURN
999 ERRORSEXITS("Petsc_ErrorHandlingSetOn",err,error)
    RETURN 1
  END SUBROUTINE Petsc_ErrorHandlingSetOn
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscFinalize routine
  SUBROUTINE Petsc_Finalise(err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_Finalise",err,error,*999)

    CALL PetscFinalize(err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PetscFinalize.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_Finalise")
    RETURN
999 ERRORSEXITS("Petsc_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_Finalise
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscInitialize routine.
  SUBROUTINE Petsc_Initialise(file,err,error,*)

    !Argument Variables
    CHARACTER(LEN=*), INTENT(IN) :: file !<Filename for PETSc options file
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_Initialise",err,error,*999)

    petscHandleError=.TRUE.
    CALL PetscInitialize(file,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PetscInitialize.",err,error,*999)
    ENDIF
    ! Disable PETSc's signal handler as we have our own OpenCMISS signal handlers in cmiss_c.c
    CALL PetscPopSignalHandler(err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PetscPopSignalHandler.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_Initialise")
    RETURN
999 ERRORSEXITS("Petsc_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_Initialise
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscLogView routine.
  SUBROUTINE Petsc_LogView(viewer,err,error,*)

    !Argument Variables
    PetscViewer, INTENT(IN) :: viewer !<The viewer to print the log to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_LogView",err,error,*999)

    CALL PetscLogView(viewer,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PetscLogView.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_LogView")
    RETURN
999 ERRORSEXITS("Petsc_LogView",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_LogView
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PetscOptionsSetValue routine.
  SUBROUTINE Petsc_OptionsSetValue(name,VALUE,err,error,*)

    !Argument Variables
    CHARACTER(LEN=*), INTENT(IN) :: name !<The name of the option to set
    CHARACTER(LEN=*), INTENT(IN) :: value !<The value of the option to set    
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_OptionSetValue",err,error,*999)

    CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS,name,VALUE,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PetscOptionsSetValue.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_OptionsSetValue")
    RETURN
999 ERRORSEXITS("Petsc_OptionsSetValue",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_OptionsSetValue
  
  !
  !================================================================================================================================
  !

  !Finalise the PETSc IS structure and destroy the IS
  SUBROUTINE Petsc_ISFinalise(is,err,error,*)

    !Argument Variables
    TYPE(PetscISType), INTENT(INOUT) :: is !<The IS to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(is%is)) THEN
      CALL Petsc_ISDestroy(is,err,error,*999)
    ENDIF
    
    EXITS("Petsc_ISFinalise")
    RETURN
999 ERRORSEXITS("Petsc_ISFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISFinalise
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc IS structure
  SUBROUTINE Petsc_ISInitialise(is,err,error,*)

    !Argument Variables
    TYPE(PetscISType), INTENT(INOUT) :: is !<The IS to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISInitialise",err,error,*999)

    is%is%v=-2
    
    EXITS("Petsc_ISInitialise")
    RETURN
999 ERRORSEXITS("Petsc_ISInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISInitialise
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISDestroy routine.
  SUBROUTINE Petsc_ISDestroy(is,err,error,*)

    !Argument Variables
    TYPE(PetscISType), INTENT(INOUT) :: is !<The index set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISDestroy",err,error,*999)

    CALL ISDestroy(is%is,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISDestroy.",err,error,*999)
    ENDIF
        
    EXITS("Petsc_ISDestroy")
    RETURN
999 ERRORSEXITS("Petsc_ISDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISDestroy
    
  !
  !
  !================================================================================================================================
  !

  !Finalise the PETSc ISColoring structure and destroy the ISColoring
  SUBROUTINE Petsc_ISColoringFinalise(iscoloring,err,error,*)

    !Argument Variables
    TYPE(PetscISColoringType), INTENT(INOUT) :: iscoloring !<The ISColoring to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISColoringFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(iscoloring%iscoloring)) THEN
      CALL Petsc_ISColoringDestroy(iscoloring,err,error,*999)
    ENDIF
    
    EXITS("Petsc_ISColoringFinalise")
    RETURN
999 ERRORSEXITS("Petsc_ISColoringFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISColoringFinalise
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc ISColoring structure
  SUBROUTINE Petsc_ISColoringInitialise(iscoloring,err,error,*)

    !Argument Variables
    TYPE(PetscISColoringType), INTENT(INOUT) :: iscoloring !<The ISColoring to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISColoringInitialise",err,error,*999)
    
    iscoloring%iscoloring%v=-2
    
    EXITS("Petsc_ISColoringInitialise")
    RETURN
999 ERRORSEXITS("Petsc_ISColoringInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISColoringInitialise
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISColoringDestroy routine.
  SUBROUTINE Petsc_ISColoringDestroy(iscoloring,err,error,*)

    !Argument Variables
    TYPE(PetscISColoringType), INTENT(INOUT) :: iscoloring !<The index set coloring
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISColoringDestroy",err,error,*999)

    CALL ISColoringDestroy(iscoloring%iscoloring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISColoringDestroy.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_ISColoringDestroy")
    RETURN
999 ERRORSEXITS("Petsc_ISColoringDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISColoringDestroy
  
  !  
  !================================================================================================================================
  !

  !Finalise the PETSc ISLocalToGlobalMapping structure and destroy the ISLocalToGlobalMapping
  SUBROUTINE Petsc_ISLocalToGlobalMappingFinalise(isLocalToGlobalMapping,err,error,*)

    !Argument Variables
    TYPE(PetscISLocalToGloabalMappingType), INTENT(INOUT) :: isLocalToGlobalMapping !<The ISLocalToGlobalMapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISLocalToGlobalMappingFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(isLocalToGlobalMapping%isLocalToGlobalMapping)) THEN
      CALL Petsc_ISLocalToGlobalMappingDestroy(isLocalToGlobalMapping,err,error,*999)
    ENDIF
    
    EXITS("Petsc_ISLocalToGlobalMappingFinalise")
    RETURN
999 ERRORSEXITS("Petsc_ISLocalToGlobalMappingFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISLocalToGlobalMappingFinalise
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc ISLocalToGlobalMapping structure
  SUBROUTINE Petsc_ISLocalToGlobalMappingInitialise(isLocalToGlobalMapping,err,error,*)

    !Argument Variables
    TYPE(PetscISLocalToGloabalMappingType), INTENT(INOUT) :: isLocalToGlobalMapping !<The ISLocalToGlobalMapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISLocalToGlobalMappingInitialise",err,error,*999)

    isLocalToGlobalMapping%isLocalToGlobalMapping%v=-2
    
    EXITS("Petsc_ISLocalToGlobalMappingInitialise")
    RETURN
999 ERRORSEXITS("Petsc_ISLocalToGlobalMappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISLocalToGlobalMappingInitialise
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingApply routine.
  SUBROUTINE Petsc_ISLocalToGlobalMappingApply(isLocalToGlobalMapping,n,idxIn,idxOut,err,error,*)

    !Argument Variables
    TYPE(PetscISLocalToGloabalMappingType), INTENT(INOUT) :: isLocalToGlobalMapping !<The local to global mapping to apply
    INTEGER(INTG), INTENT(IN) :: n !<The number of indicies
    INTEGER(INTG), INTENT(IN) :: idxIn(:)
    INTEGER(INTG), INTENT(OUT) :: idxOut(:)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISLocalToGlobalMappingApply",err,error,*999)

    CALL ISLocalToGlobalMappingApply(isLocalToGlobalMapping%isLocalToGlobalMapping,n,idxIn,idxOut,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISLocalToGlobalMappingApply.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_ISLocalToGlobalMappingApply")
    RETURN
999 ERRORSEXITS("Petsc_ISLocalToGlobalMappingApply",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISLocalToGlobalMappingApply
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingApplyIS routine.
  SUBROUTINE Petsc_ISLocalToGlobalMappingApplyIS(isLocalToGlobalMapping,isIn,isOut,err,error,*)

    !Argument Variables
    TYPE(PetscISLocalToGloabalMappingType), INTENT(IN) :: isLocalToGlobalMapping !<The local to global mapping to apply
    TYPE(PetscISType), INTENT(IN) :: isIn
    TYPE(PetscISType), INTENT(OUT) :: isOut
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISLocalToGlobalMappingApplyIS",err,error,*999)

    CALL ISLocalToGlobalMappingApplyIS(isLocalToGlobalMapping%isLocalToGlobalMapping,isIn%is,isOut%is,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISLocalToGlobalMappingApplyIS.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_ISLocalToGlobalMappingApplyIS")
    RETURN
999 ERRORSEXITS("Petsc_ISLocalToGlobalMappingApplyIS",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISLocalToGlobalMappingApplyIS
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingCreate routine.
  SUBROUTINE Petsc_ISLocalToGlobalMappingCreate(communicator,blockSize,n,indices,mode,isLocalToGlobalMapping,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator
#endif    
    INTEGER(INTG), INTENT(IN) :: blockSize !<The block size
    INTEGER(INTG), INTENT(IN) :: n !<The number of local indices divided by the block size (or the number of block indices)
    INTEGER(INTG), INTENT(IN) :: indices(:) !<The global index for each local element
    PetscCopyMode, INTENT(IN) :: mode !<The copy mode
    TYPE(PetscISLocalToGloabalMappingType), INTENT(INOUT) :: isLocalToGlobalMapping !<On exit, the local to global mapping context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISLocalToGlobalMappingCreate",err,error,*999)

    CALL ISLocalToGlobalMappingCreate(communicator,blockSize,n,indices,mode,isLocalToGlobalMapping%isLocalToGlobalMapping,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISLocalToGlobalMappingCreate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_ISLocalToGlobalMappingCreate")
    RETURN
999 ERRORSEXITS("Petsc_ISLocalToGlobalMappingCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISLocalToGlobalMappingCreate
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc ISLocalToGlobalMappingDestroy routine.
  SUBROUTINE Petsc_ISLocalToGlobalMappingDestroy(isLocalToGlobalMapping,err,error,*)

    !Argument Variables
    TYPE(PetscISLocalToGloabalMappingType), INTENT(INOUT) :: isLocalToGlobalMapping !<The local to global mapping context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_ISLocalToGlobalMappingDestroy",err,error,*999)

    CALL ISLocalToGlobalMappingDestroy(isLocalToGlobalMapping%isLocalToGlobalMapping,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in ISLocalToGlobalMappingDestroy.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_ISLocalToGlobalMappingDestroy")
    RETURN
999 ERRORSEXITS("Petsc_ISLocalToGlobalMappingDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_ISLocalToGlobalMappingDestroy
    
  !
  !================================================================================================================================
  !

  !>Finalise the PETSc KSP structure and destroy the KSP
  SUBROUTINE Petsc_KSPFinalise(ksp,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(ksp%ksp)) THEN
      CALL Petsc_KSPDestroy(ksp,err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPFinalise")
    RETURN
999 ERRORSEXITS("Petsc_KSPFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPFinalise
    
  !
  !================================================================================================================================
  !

  !>Initialise the PETSc KSP structure
  SUBROUTINE Petsc_KSPInitialise(ksp,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPInitialise",err,error,*999)

    ksp%ksp%v=-2
    
    EXITS("Petsc_KSPInitialise")
    RETURN
999 ERRORSEXITS("Petsc_KSPInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPInitialise
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPCreate routine
  SUBROUTINE Petsc_KSPCreate(communicator,ksp,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator for the KSP creation
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator for the KSP creation
#endif    
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<On exit, the KSP object
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPCreate",err,error,*999)

    CALL KSPCreate(communicator,ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPCreate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPCreate")
    RETURN
999 ERRORSEXITS("Petsc_KSPCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPCreate
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPDestroy routine
  SUBROUTINE Petsc_KSPDestroy(ksp,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPDestroy",err,error,*999)

    CALL KSPDestroy(ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPDestroy.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPDestroy")
    RETURN
999 ERRORSEXITS("Petsc_KSPDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPDestroy
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetConvergedReason routine
  SUBROUTINE Petsc_KSPGetConvergedReason(ksp,reason,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP information
    INTEGER(INTG), INTENT(OUT) :: reason !<On exit, the converged reason
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPGetConvergedReason",err,error,*999)

    CALL KSPGetConvergedReason(ksp%ksp,reason,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPGetConvergedReason.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPGetConvergedReason")
    RETURN
999 ERRORSEXITS("Petsc_KSPGetConvergedReason",err,error)
    RETURN 1
  END SUBROUTINE Petsc_KSPGetConvergedReason
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetIterationNumber routine
  SUBROUTINE Petsc_KSPGetIterationNumber(ksp,iterationNumber,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP information
    INTEGER(INTG), INTENT(OUT) :: iterationNumber !<On exit, the number of iterations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPGetIterationNumber",err,error,*999)

    CALL KSPGetIterationNumber(ksp%ksp,iterationNumber,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPGetIterationNumber.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPGetIterationNumber")
    RETURN
999 ERRORSEXITS("Petsc_KSPGetIterationNumber",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPGetIterationNumber
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetPC routine
  SUBROUTINE Petsc_KSPGetPC(ksp,pc,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to get the PC for
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<On exit, the PC associated with the KSP
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPGetPC",err,error,*999)

    CALL KSPGetPC(ksp%ksp,pc%pc,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPGetPC.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPGetPC")
    RETURN
999 ERRORSEXITS("Petsc_KSPGetPC",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPGetPC
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGetResidualNorm routine
  SUBROUTINE Petsc_KSPGetResidualNorm(ksp,residualNorm,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to get the PC for
    REAL(DP), INTENT(OUT) :: residualNorm !<On exit, the residual norm for the KSP
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPGetResidualNorm",err,error,*999)

    CALL KSPGetResidualNorm(ksp%ksp,residualNorm,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPGetResidualNorm.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPGetResidualNorm")
    RETURN
999 ERRORSEXITS("Petsc_KSPGetResidualNorm",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPGetResidualNorm
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPGMRESSetRestart routine
  SUBROUTINE Petsc_KSPGMRESSetRestart(ksp,restart,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to set the GMRES restart for
    INTEGER(INTG), INTENT(OUT) :: restart !<The restart value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPGMRESSetRestart",err,error,*999)

    CALL KSPGMRESSetRestart(ksp%ksp,restart,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPGMRESSetRestart.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPGMRESSetRestart")
    RETURN
999 ERRORSEXITS("Petsc_KSPGMRESSetRestart",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPGMRESSetRestart
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetFromOptions routine
  SUBROUTINE Petsc_KSPSetFromOptions(ksp,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set the options for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPSetFromOptions",err,error,*999)

    CALL KSPSetFromOptions(ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetFromOptions.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPSetFromOptions")
    RETURN
999 ERRORSEXITS("Petsc_KSPSetFromOptions",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPSetFromOptions
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetInitialGuessNonzero routine
  SUBROUTINE Petsc_KSPSetInitialGuessNonZero(ksp,flag,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to set the initial guess non zero for
    LOGICAL, INTENT(IN) :: flag !<If flag is .TRUE. the initial guess is non-zero, if flag is .FALSE. the initial guess is zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPSetInitialGuessNonZero",err,error,*999)

    IF(FLAG) THEN
      CALL KSPSetInitialGuessNonzero(ksp%ksp,PETSC_TRUE,err)
    ELSE
      CALL KSPSetInitialGuessNonzero(ksp%ksp,PETSC_FALSE,err)
    ENDIF
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetInitialGuessNonzero.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPSetInitialGuessNonZero")
    RETURN
999 ERRORSEXITS("Petsc_KSPSetInitialGuessNonZero",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPSetInitialGuessNonZero
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetOperators routine
  SUBROUTINE Petsc_KSPSetOperators(ksp,amat,pmat,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set the operators for
    TYPE(PetscMatType), INTENT(INOUT) :: amat !<The matrix associated with the linear system
    TYPE(PetscMatType), INTENT(INOUT) :: pmat !<The matrix to be used in constructing the preconditioner
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPSetOperators",err,error,*999)

    CALL KSPSetOperators(ksp%ksp,amat%mat,pmat%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetFromOperators.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPSetOperators")
    RETURN
999 ERRORSEXITS("Petsc_KSPSetOperators",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPSetOperators
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetReusePreconditioner routine
  SUBROUTINE Petsc_KSPSetReusePreconditioner(ksp,flag,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to set the options for
    LOGICAL, INTENT(IN) :: flag !<If flag is .TRUE. then the current precondition is reused, if flag is .FALSE. it is not
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPSetReusePreconditioner",err,error,*999)

    IF(flag) THEN
      CALL KSPSetReusePreconditioner(ksp%ksp,PETSC_TRUE,err)
    ELSE
      CALL KSPSetReusePreconditioner(ksp%ksp,PETSC_FALSE,err)
    ENDIF
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetReusePreconditioner.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPSetReusePreconditioner")
    RETURN
999 ERRORSEXITS("Petsc_KSPSetReusePreconditioner",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPSetReusePreconditioner
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetTolerances routine
  SUBROUTINE Petsc_KSPSetTolerances(ksp,rTol,aTol,dTol,maxIterations,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to set the tolerances for
    REAL(DP), INTENT(IN) :: rTol !<The relative tolerance to set
    REAL(DP), INTENT(IN) :: aTol !<The absolution tolerance to set
    REAL(DP), INTENT(IN) :: dTol !<The divergence tolerance to set
    INTEGER(INTG), INTENT(IN) :: maxIterations !<The maximum number of iterations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPSetTolerances",err,error,*999)

    CALL KSPSetTolerances(ksp%ksp,rTol,aTol,dTol,maxIterations,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetTolerances.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPSetTolerances")
    RETURN
999 ERRORSEXITS("Petsc_KSPSetTolerances",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPSetTolerances
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetType routine
  SUBROUTINE Petsc_KSPSetType(ksp,method,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to set the type for
    KSPType, INTENT(IN) :: method !<The KSP method
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPSetType",err,error,*999)

    CALL KSPSetType(ksp%ksp,method,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPSetType")
    RETURN
999 ERRORSEXITS("Petsc_KSPSetType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPSetType
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSetUp routine
  SUBROUTINE Petsc_KSPSetUp(ksp,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to set up
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPSetUp",err,error,*999)

    CALL KSPSetUp(ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSetUp.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPSetUp")
    RETURN
999 ERRORSEXITS("Petsc_KSPSetUp",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPSetUp
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc KSPSolve routine
  SUBROUTINE Petsc_KSPSolve(ksp,b,x,err,error,*)

    !Argument Variables
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The Ksp to set up
    TYPE(PetscVecType), INTENT(INOUT) :: b !<The RHS vector
    TYPE(PetscVecType), INTENT(INOUT)  :: x !<The solution vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_KSPSolve",err,error,*999)

    CALL KSPSolve(ksp%ksp,b%vec,x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in KSPSolve.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_KSPSolve")
    RETURN
999 ERRORSEXITS("Petsc_KSPSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_KSPSolve
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc Mat structure
  SUBROUTINE Petsc_MatFinalise(a,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(a%mat)) THEN
      CALL Petsc_MatDestroy(a,err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFinalise")
    RETURN
999 ERRORSEXITS("Petsc_MatFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFinalise
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc Mat structure
  SUBROUTINE Petsc_MatInitialise(a,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatInitialise",err,error,*999)

    a%mat%v=-2
     
    EXITS("Petsc_MatInitialise")
    RETURN
999 ERRORSEXITS("Petsc_MatInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatInitialise
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatAssemblyBegin routine.
  SUBROUTINE Petsc_MatAssemblyBegin(A,assemblyType,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !The matrix to assemble
    MatAssemblyType, INTENT(IN) :: assemblyType !<The assembly type 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatAssemblyBegin",err,error,*999)

    CALL MatAssemblyBegin(A%mat,assemblyType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatAssemblyBegin.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatAssemblyBegin")
    RETURN
999 ERRORSEXITS("Petsc_MatAssemblyBegin",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatAssemblyBegin
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatAssemblyEnd routine.
  SUBROUTINE Petsc_MatAssemblyEnd(A,assemblyType,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to assemble
    MatAssemblyType, INTENT(IN) :: assemblyType !<The assembly type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatAssemblyEnd",err,error,*999)

    CALL MatAssemblyEnd(A%mat,assemblyType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatAssemblyEnd.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatAssemblyEnd")
    RETURN
999 ERRORSEXITS("Petsc_MatAssemblyEnd",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatAssemblyEnd
       
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatCreate routine.
  SUBROUTINE Petsc_MatCreate(communicator,A,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI Communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI Communicator
#endif    
    TYPE(PetscMatType), INTENT(INOUT) :: A !<On exit, the created matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatCreate",err,error,*999)

    CALL MatCreate(communicator,A%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatCreate")
    RETURN
999 ERRORSEXITS("Petsc_MatCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatCreate
    
  !
  !================================================================================================================================
  !

 !>Buffer routine to the PETSc MatCreateAIJ routine.
  SUBROUTINE Petsc_MatCreateAIJ(communicator,localM,localN,globalM,globalN,diagNumberNonZerosPerRow,diagNumberNonZerosEachRow, &
    & offDiagNumberNonZerosPerRow,offDiagNumberNonZerosEachRow,a,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator
#endif    
    INTEGER(INTG), INTENT(IN) :: localM !<The number of local rows
    INTEGER(INTG), INTENT(IN) :: localN !<The number of local columns
    INTEGER(INTG), INTENT(IN) :: globalM !<The number of global rows
    INTEGER(INTG), INTENT(IN) :: globalN !<The number of global columns
    INTEGER(INTG), INTENT(IN) :: diagNumberNonZerosPerRow !<The maximum number of non-zeros per row in the diagonal part of the matrix
    INTEGER(INTG), INTENT(IN) :: diagNumberNonZerosEachRow(:) !<The number of non-zeros per row in the diagonal part of the matrix
    INTEGER(INTG), INTENT(IN) :: offDiagNumberNonZerosPerRow !<The maximum number of non-zeros per row in the off-diagonal part of the matrix
    INTEGER(INTG), INTENT(IN) :: offDiagNumberNonZerosEachRow(:) !<The number of non-zeros per row in the off-diagonal part of the matrix
    TYPE(PetscMatType), INTENT(INOUT) :: a !<On exit, the created matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatCreateAIJ",err,error,*999)

    CALL MatCreateAIJ(communicator,localM,localN,globalM,globalN,diagNumberNonZerosPerRow,diagNumberNonZerosEachRow, &
      & offDiagNumberNonZerosPerRow,offDiagNumberNonZerosEachRow,a%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreateAIJ.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatCreateAIJ")
    RETURN
999 ERRORSEXITS("Petsc_MatCreateAIJ",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatCreateAIJ
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatCreateDense routine.
  SUBROUTINE Petsc_MatCreateDense(communicator,localM,localN,globalM,globalN,matrixData,a,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator
#endif    
    INTEGER(INTG), INTENT(IN) :: localM !<The number of local rows
    INTEGER(INTG), INTENT(IN) :: localN !<The number of local columns
    INTEGER(INTG), INTENT(IN) :: globalM !<The number of global columns
    INTEGER(INTG), INTENT(IN) :: globalN !<The number of global rows
    REAL(DP), INTENT(IN) :: matrixData !<Optional, the allocated matrix data.
    TYPE(PetscMatType), INTENT(INOUT) :: a !<On exit, the created matrix 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatCreateDense",err,error,*999)

    CALL MatCreateDense(communicator,localM,localN,globalM,globalN,matrixData,a%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreateDense.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatCreateDense")
    RETURN
999 ERRORSEXITS("Petsc_MatCreateDense",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatCreateDense
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatCreateSeqAIJ routine.
  SUBROUTINE Petsc_MatCreateSeqAIJ(communicator,m,n,numberNonZerosPerRow,numberNonZerosEachRow,a,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator
#endif    
    INTEGER(INTG), INTENT(IN) :: m !<The number of rows
    INTEGER(INTG), INTENT(IN) :: n !<The number of columns
    INTEGER(INTG), INTENT(IN) :: numberNonZerosPerRow !<The maximum number of non-zeros per row
    INTEGER(INTG), INTENT(IN) :: numberNonZerosEachRow(:) !<The number of non-zeros in each row
    TYPE(PetscMatType), INTENT(INOUT) :: a !<On exit, the created matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatCreateSeqAIJ",err,error,*999)

    CALL MatCreateSeqAIJ(communicator,m,n,numberNonZerosPerRow,numberNonZerosEachRow,a%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreateSeqAIJ.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatCreateSeqAIJ")
    RETURN
999 ERRORSEXITS("Petsc_MatCreateSeqAIJ",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatCreateSeqAIJ
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatCreateSeqDense routine.
  SUBROUTINE Petsc_MatCreateSeqDense(communicator,m,n,matrixData,a,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI Communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI Communicator
#endif    
    INTEGER(INTG), INTENT(IN) :: m !<The number of rows
    INTEGER(INTG), INTENT(IN) :: n !<The number of columns
    REAL(DP), INTENT(IN) :: matrixData(*) !<Optional, the allocated matrix data
    TYPE(PetscMatType), INTENT(INOUT) :: a !<On exit, the created matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatCreateSeqDense",err,error,*999)

    CALL MatCreateSeqDense(communicator,m,n,matrixData,a%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatCreateSeqDense.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatCreateSeqDense")
    RETURN
999 ERRORSEXITS("Petsc_MatCreateSeqDense",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatCreateSeqDense

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatDenseGetArrayF90 routine.
  SUBROUTINE Petsc_MatDenseGetArrayF90(a,array,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT), TARGET :: a !<The matrix to get the array for
    REAL(DP), POINTER :: array(:,:) !<On exit, a pointer to the matrix array
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatDenseGetArrayF90",err,error,*999)

    IF(ASSOCIATED(array)) THEN
      CALL FlagError("Array is already associated.",err,error,*999)
    ELSE
      CALL MatDenseGetArrayF90(a%mat,array,err)
      IF(err/=0) THEN
        IF(petscHandleError) THEN
          CHKERRQ(err)
        ENDIF
        CALL FlagError("PETSc error in MatDenseGetArrayF90.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Petsc_MatDenseGetArrayF90")
    RETURN
999 ERRORSEXITS("Petsc_MatDenseGetArrayF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatDenseGetArrayF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatDenseRestoreArrayF90 routine.
  SUBROUTINE Petsc_MatDenseRestoreArrayF90(a,array,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to restore the array for
    REAL(DP), POINTER :: array(:,:) !<A pointer to the matrix array to restore
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatDenseRestoreArrayF90",err,error,*999)

    CALL MatDenseRestoreArrayF90(a%mat,array,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatDenseRestoreArrayF90.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatDenseRestoreArrayF90")
    RETURN
999 ERRORSEXITS("Petsc_MatDenseRestoreArrayF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatDenseRestoreArrayF90
    
  !
  !================================================================================================================================
  !    

  !>Buffer routine to the PETSc MatDestroy routine.
  SUBROUTINE Petsc_MatDestroy(a,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatDestroy",err,error,*999)

    CALL MatDestroy(a%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatDestroy.",err,error,*999)
    ENDIF
     
    EXITS("Petsc_MatDestroy")
    RETURN
999 ERRORSEXITS("Petsc_MatDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatDestroy
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetInfo routine.
  SUBROUTINE Petsc_MatGetInfo(a,flag,info,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to get the information of
    MatInfoType, INTENT(IN) :: flag !<The matrix information collective to return
    MatInfo, INTENT(OUT) :: info(MAT_INFO_SIZE) !<On return, the matrix information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatGetInfo",err,error,*999)

    CALL MatGetInfo(a%mat,flag,info,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatGetInfo.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatGetInfo")
    RETURN
999 ERRORSEXITS("Petsc_MatGetInfo",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatGetInfo
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetOwnershipRange routine.
  SUBROUTINE Petsc_MatGetOwnershipRange(a,firstRow,lastRow,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to get the ownership range of
    INTEGER(INTG), INTENT(OUT) :: firstRow !<On exit, the first row for the matrix
    INTEGER(INTG), INTENT(OUT) :: lastRow !<On exit, the last row for the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatGetOwnershipRange",err,error,*999)

    CALL MatGetOwnershipRange(a%mat,firstRow,lastRow,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatGetOwnershipRange.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatGetOwnershipRange")
    RETURN
999 ERRORSEXITS("Petsc_MatGetOwnershipRange",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatGetOwnershipRange
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetRow routine.
  SUBROUTINE Petsc_MatGetRow(A,rowNumber,numberOfColumns,columns,values,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to get the array for
    INTEGER(INTG), INTENT(IN) :: rowNumber !<The row number to get the row values for
    INTEGER(INTG), INTENT(OUT) :: numberOfColumns !<On return, the number of nonzero columns in the row
    INTEGER(INTG), INTENT(OUT) :: columns(:) !<On return, the column numbers for the nonzero columns in the row
    REAL(DP), INTENT(OUT) :: values(:) !<On exit, the nonzero values in the row.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatGetRow",err,error,*999)

    CALL MatGetRow(A%mat,rowNumber,numberOfColumns,columns,values,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatGetRow.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatGetRow")
    RETURN
999 ERRORSEXITS("Petsc_MatGetRow",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatGetRow
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatGetValues routine.
  SUBROUTINE Petsc_MatGetValues(a,m,mIndices,n,nIndices,values,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to get the values of
    INTEGER(INTG), INTENT(IN) :: m !<The number of row indices
    INTEGER(INTG), INTENT(IN) :: mIndices(*) !<The row indices
    INTEGER(INTG), INTENT(IN) :: n !<The number of column indices
    INTEGER(INTG), INTENT(IN) :: nIndices(*) !<The column indices
    REAL(DP), INTENT(OUT) :: values(*) !<The values to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatGetValues",err,error,*999)

    CALL MatGetValues(a%mat,m,mIndices,n,nIndices,values,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatGetValues.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatGetValues")
    RETURN
999 ERRORSEXITS("Petsc_MatGetValues",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatGetValues
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatMumpsSetIcntl routine.
  SUBROUTINE Petsc_MatMumpsSetIcntl(factoredMatrix,icntl,ival,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: factoredMatrix !<The factored matrix from PETSc-MUMPS interface
    INTEGER(INTG), INTENT(IN) :: icntl !<The MUMPS ICNTL integer control parameter
    INTEGER(INTG), INTENT(IN) :: ival !<The MUMPS ICNTL integer value to set: ICNTL(icntl)=ival
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatMumpsSetIcntl",err,error,*999)

#ifdef WITH_MUMPS    
    CALL MatMumpsSetIcntl(factoredMatrix%mat,icntl,ival,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatMumpsSetIcntl.",err,error,*999)
    ENDIF
#endif    
    
    EXITS("Petsc_MatMumpsSetIcntl")
    RETURN
999 ERRORSEXITS("Petsc_MatMumpsSetIcntl",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatMumpsSetIcntl

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatMumpsSetCntl routine.
  SUBROUTINE Petsc_MatMumpsSetCntl(factoredMatrix,icntl,val,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: factoredMatrix !<The factored matrix from PETSc-MUMPS interface
    INTEGER(INTG), INTENT(IN) :: icntl !<The MUMPS CNTL integer control parameter
    REAL(DP), INTENT(IN) :: val !<The MUMPS CNTL real value to set: CNTL(icntl)=val
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatMumpsSetCntl",err,error,*999)

#ifdef WITH_MUMPS    
    CALL MatMumpsSetCntl(factoredMatrix%mat,icntl,val,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatMumpsSetCntl.",err,error,*999)
    ENDIF
#endif    
    
    EXITS("Petsc_MatMumpsSetCntl")
    RETURN
999 ERRORSEXITS("Petsc_MatMumpsSetCntl",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatMumpsSetCntl

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatRestoreRow routine.
  SUBROUTINE Petsc_MatRestoreRow(A,rowNumber,numberOfColumns,columns,values,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: A !<The matrix to restore the row for
    INTEGER(INTG), INTENT(IN) :: rowNumber !<The row number to restore the row values for
    INTEGER(INTG), INTENT(OUT) :: numberOfColumns !<The number of nonzero columns in the row
    INTEGER(INTG), INTENT(OUT) :: columns(:) !<The column numbers for the nonzero columns in the row
    REAL(DP), INTENT(OUT) :: values(:) !<The nonzero values in the row to restore.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatRestoreRow",err,error,*999)

    CALL MatRestoreRow(A%mat,rowNumber,numberOfColumns,columns,values,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatRestoreRow.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatRestoreRow")
    RETURN
999 ERRORSEXITS("Petsc_MatRestoreRow",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatRestoreRow
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSeqAIJGetArrayF90 routine.
  SUBROUTINE Petsc_MatSeqAIJGetArrayF90(a,array,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT), TARGET :: a !<The matrix to get the array for
    REAL(DP), POINTER :: array(:) !<On exit, a pointer to the matrix array
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSeqAIJGetArrayF90",err,error,*999)

    IF(ASSOCIATED(array)) THEN
      CALL FlagError("Array is already associated.",err,error,*999)
    ELSE
      CALL MatSeqAIJGetArrayF90(a%mat,array,err)
      IF(err/=0) THEN
        IF(petscHandleError) THEN
          CHKERRQ(err)
        ENDIF
        CALL FlagError("PETSc error in MatSeqAIJGetArrayF90.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Petsc_MatSeqAIJGetArrayF90")
    RETURN
999 ERRORSEXITS("Petsc_MatSeqAIJGetArrayF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSeqAIJGetArrayF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSeqAIJGetMaxRowNonzeros routine.
  SUBROUTINE Petsc_MatSeqAIJGetMaxRowNonzeros(a,maxNumberNonZeros,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to get the maximum number of non zeros for
    INTEGER(INTG), INTENT(OUT) :: maxNumberNonZeros!<On exit, the maximum number of non zeros in any row of the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSeqAIJGetMaxRowNonzeros",err,error,*999)

    !CALL MatSeqAIJGetMaxRowNonzeros(A%mat,maxNumberNonZeros,err)
    maxNumberNonZeros=0
    CALL FlagError("Not implemented.",err,error,*999)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSeqAIJGetMaxRowNonzeros.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatSeqAIJGetMaxRowNonzeros")
    RETURN
999 ERRORSEXITS("Petsc_MatSeqAIJGetMaxRowNonzeros",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSeqAIJGetMaxRowNonzeros
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatSeqAIJRestoreArrayF90 routine.
  SUBROUTINE Petsc_MatSeqAIJRestoreArrayF90(a,array,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to restore the array for
    REAL(DP), POINTER :: array(:) !<A pointer to the matrix array to restore
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSeqAIJRestoreArrayF90",err,error,*999)

    CALL MatSeqAIJRestoreArrayF90(a%mat,array,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSeqAIJRestoreArrayF90.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatSeqAIJRestoreArrayF90")
    RETURN
999 ERRORSEXITS("Petsc_MatSeqAIJRestoreArrayF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSeqAIJRestoreArrayF90
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetLocalToGlobalMapping routine.
  SUBROUTINE Petsc_MatSetLocalToGlobalMapping(a,rowMapping,columnMapping,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to set the local to global mapping for
    TYPE(PetscISLocalToGloabalMappingType), INTENT(IN) :: rowMapping !<The row local to global mapping 
    TYPE(PetscISLocalToGloabalMappingType), INTENT(IN) :: columnMapping !<The column local to global mapping 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSetLocalToGlobalMapping",err,error,*999)

    CALL MatSetLocalToGlobalMapping(a%mat,rowMapping%isLocalToGlobalMapping,columnMapping%isLocalToGlobalMapping,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetLocalToGlobalMapping.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatSetLocalToGlobalMapping")
    RETURN
999 ERRORSEXITS("Petsc_MatSetLocalToGlobalMapping",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSetLocalToGlobalMapping
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetOption routine.
  SUBROUTINE Petsc_MatSetOption(a,option,flag,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to set the option for
    MatOption, INTENT(IN) :: option !<The option to set \see OpenCMISSPETSc_PetscMatOptionTypes,OpenCMISSPETSc
    LOGICAL, INTENT(IN) :: flag !<The option flag to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSetOption",err,error,*999)

    IF(flag) THEN
      CALL MatSetOption(a%mat,option,PETSC_TRUE,err)
    ELSE
      CALL MatSetOption(a%mat,option,PETSC_FALSE,err)
    ENDIF
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetOption.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatSetOption")
    RETURN
999 ERRORSEXITS("Petsc_MatSetOption",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSetOption
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetSizes routine.
  SUBROUTINE Petsc_MatSetSizes(a,localM,localN,globalM,globalN,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to set the size of
    INTEGER(INTG), INTENT(IN) :: localM !<Number of local rows
    INTEGER(INTG), INTENT(IN) :: localN !<Number of local columns
    INTEGER(INTG), INTENT(IN) :: globalM !<Number of global rows
    INTEGER(INTG), INTENT(IN) :: globalN !<Number of global columns
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSetSizes",err,error,*999)

    CALL MatSetSizes(a%mat,localM,localN,globalM,globalN,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetSizes.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatSetSizes")
    RETURN
999 ERRORSEXITS("Petsc_MatSetSizes",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSetSizes
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatSetType routine.
  SUBROUTINE Petsc_MatSetType(a,matrixType,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to set the type of
    MatType, INTENT(IN) :: matrixType !<The matrix type 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSetType",err,error,*999)

    CALL MatSetType(a%mat,matrixType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatSetType")
    RETURN
999 ERRORSEXITS("Petsc_MatSetType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSetType

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetValue routine.
  SUBROUTINE Petsc_MatSetValue(a,row,col,value,insertMode,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: row !<The row index
    INTEGER(INTG), INTENT(IN) :: col !<The column index
    REAL(DP), INTENT(IN) :: value !<The value to set
    InsertMode, INTENT(IN) :: insertMode !<The insert mode \see OpenCMISSPETSc_PetscMatInsertMode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSetValue",err,error,*999)

    CALL MatSetValue(a%mat,row,col,value,insertMode,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetValue.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatSetValue")
    RETURN
999 ERRORSEXITS("Petsc_MatSetValue",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSetValue
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetValues routine.
  SUBROUTINE Petsc_MatSetValues(a,m,mIndices,n,nIndices,values,insertMode,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: m !<The number of row indices
    INTEGER(INTG), INTENT(IN) :: mIndices(*) !<The row indices
    INTEGER(INTG), INTENT(IN) :: n !<The number of column indices
    INTEGER(INTG), INTENT(IN) :: nIndices(*) !<The column indices
    REAL(DP), INTENT(IN) :: values(*) !<The values to set
    InsertMode, INTENT(IN) :: insertMode !<The insert mode \see OpenCMISSPETSc_PetscMatInsertMode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSetValues",err,error,*999)

    CALL MatSetValues(a%mat,m,mIndices,n,nIndices,values,insertMode,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetValues.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatSetValues")
    RETURN
999 ERRORSEXITS("Petsc_MatSetValues",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSetValues
    
  !
  !================================================================================================================================
  !
  !>Buffer routine to the PETSc MatSetValueLocal routine.
  SUBROUTINE Petsc_MatSetValueLocal(a,row,col,VALUE,insertMode,err,error,*)
    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: row !<The row index
    INTEGER(INTG), INTENT(IN) :: col !<The column index
    REAL(DP), INTENT(IN) :: value !<The value to set
    InsertMode, INTENT(IN) :: insertMode !<The insert mode \see OpenCMISSPETSc_PetscMatInsertMode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSetValueLocal",err,error,*999)

    CALL MatSetValueLocal(a%mat,row,col,value,insertMode,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetValueLocal.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatSetValueLocal")
    RETURN
999 ERRORSEXITS("Petsc_MatSetValueLocal",err,error)
    RETURN 1
  END SUBROUTINE Petsc_MatSetValueLocal

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatSetValuesLocal routine.
  SUBROUTINE Petsc_MatSetValuesLocal(a,m,mIndices,n,nIndices,values,insertMode,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to set the values of
    INTEGER(INTG), INTENT(IN) :: m !<The number of row indices
    INTEGER(INTG), INTENT(IN) :: mIndices(:) !<The row indices
    INTEGER(INTG), INTENT(IN) :: n !<The number of column indices
    INTEGER(INTG), INTENT(IN) :: nIndices(:) !<The column indices
    REAL(DP), INTENT(IN) :: values(:) !<The values to set
    InsertMode, INTENT(IN) :: insertMode !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatSetValuesLocal",err,error,*999)

    CALL MatSetValuesLocal(a%mat,m,mIndices,n,nIndices,values,insertMode,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatSetValuesLocal.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatSetValuesLocal")
    RETURN
999 ERRORSEXITS("Petsc_MatSetValuesLocal",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatSetValuesLocal
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatView routine.
  SUBROUTINE Petsc_MatView(a,viewer,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to view
    PetscViewer, INTENT(IN) :: viewer !<The viewer
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatView",err,error,*999)

    CALL MatView(a%mat,viewer,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatView.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatView")
    RETURN
999 ERRORSEXITS("Petsc_MatView",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatView
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatZeroEntries routine.
  SUBROUTINE Petsc_MatZeroEntries(a,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to zero the entries of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatZeroEntries",err,error,*999)

    CALL MatZeroEntries(a%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatZeroEntries.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatZeroEntries")
    RETURN
999 ERRORSEXITS("Petsc_MatZeroEntries",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatZeroEntries
    
  !
  !================================================================================================================================
  !
  
  !Finalise the PETSc MatColoring structure and destroy the MatColoring
  SUBROUTINE Petsc_MatColoringFinalise(matColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The matColoring to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(matColoring%matColoring)) THEN
      CALL Petsc_MatColoringDestroy(matColoring,err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatColoringFinalise")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringFinalise
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc MatColoring structure
  SUBROUTINE Petsc_MatColoringInitialise(matColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The matColoring to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringInitialise",err,error,*999)

    matColoring%matColoring%v=-2
    
    EXITS("Petsc_MatColoringInitialise")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringInitialise
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatColoringApply routine.
  SUBROUTINE Petsc_MatColoringApply(matColoring,isColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The coloring to apply
    TYPE(PetscISColoringType), INTENT(INOUT) :: isColoring !<The index set coloring to apply to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringApply",err,error,*999)

    CALL MatColoringApply(matColoring%matColoring,isColoring%isColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatColoringApply.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatColoringApply")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringApply",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringApply
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatColoringCreate routine.
  SUBROUTINE Petsc_MatColoringCreate(a,matColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The matrix to create the coloring for
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<On exit, the created coloring
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringCreate",err,error,*999)

    CALL MatColoringCreate(a%mat,matColoring%matColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatColoringCreate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatColoringCreate")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringCreate
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatColoringDestroy routine.
  SUBROUTINE Petsc_MatColoringDestroy(matColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The coloring to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringDestroy",err,error,*999)

    CALL MatColoringDestroy(matColoring%matColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatColoringDestroy.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatColoringDestroy")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringDestroy
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatColoringSetFromOptions routine.
  SUBROUTINE Petsc_MatColoringSetFromOptions(matColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The coloring to set the options for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringSetFromOptions",err,error,*999)

    CALL MatColoringSetFromOptions(matColoring%matColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatColoringSetFromOptions.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatColoringSetFromOptions")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringSetFromOptions",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringSetFromOptions
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc MatColoringSetType routine.
  SUBROUTINE Petsc_MatColoringSetType(matColoring,coloringType,err,error,*)

    !Argument Variables
    TYPE(PetscMatColoringType), INTENT(INOUT) :: matColoring !<The coloring to set the type for
    MatColoringType, INTENT(IN) :: coloringType !<The type of coloring to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatColoringSetType",err,error,*999)

    CALL MatColoringSetType(matColoring%matColoring,coloringType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatColoringSetType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatColoringSetType")
    RETURN
999 ERRORSEXITS("Petsc_MatColoringSetType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatColoringSetType
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc MatFDColoring structure and destroy the MatFDColoring
  SUBROUTINE Petsc_MatFDColoringFinalise(matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The MatFDColoring to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(matFDColoring%matFDColoring)) THEN
      CALL Petsc_MatFDColoringDestroy(matFDColoring,err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringFinalise")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringFinalise
    
  !
  !================================================================================================================================
  !
  
  !Initialise the PETSc MatFDColoring structure
  SUBROUTINE Petsc_MatFDColoringInitialise(matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The MatFDColoring to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringInitialise",err,error,*999)

    matFDColoring%matFDColoring%v=-2
     
    EXITS("Petsc_MatFDColoringInitialise")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringInitialise
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringCreate routine.
  SUBROUTINE Petsc_MatFDColoringCreate(a,isColoring,matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The PETSc matrix to create the FD coloring for
    TYPE(PetscISColoringType), INTENT(IN) :: isColoring !<The index set coloring to create the finite difference coloring for
    TYPE(PetscMatFDColoringType), INTENT(OUT) :: matFDColoring !<On exit, the matrix finite difference coloring
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringCreate",err,error,*999)

    CALL MatFDColoringCreate(a%mat,isColoring%isColoring,matFDColoring%matFDColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringCreate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringCreate")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringCreate
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringDestroy routine.
  SUBROUTINE Petsc_MatFDColoringDestroy(matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The matrix finite difference coloring to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringDestroy",err,error,*999)

    CALL MatFDColoringDestroy(matFDColoring%matFDColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringDestroy.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringDestroy")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringDestroy
    
  !
  !================================================================================================================================2
  !

  !>Buffer routine to the PETSc MatFDColoringSetFromOptions routine.
  SUBROUTINE Petsc_MatFDColoringSetFromOptions(matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The matrix finite difference coloring to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringSetFromOptions",err,error,*999)

    CALL MatFDColoringSetFromOptions(matFDColoring%matFDColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringSetFromOptions.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringSetFromOptions")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringSetFromOptions",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringSetFromOptions
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringSetParameters routine.
  SUBROUTINE Petsc_MatFDColoringSetParameters(matFDColoring,rError,uMin,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The matrix finite difference coloring to set
    REAL(DP) :: rError !<The square root of the relative error
    REAL(DP) :: uMin !<MatFDColoring umin option
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringSetParameters",err,error,*999)

    CALL MatFDColoringSetParameters(matFDColoring%matFDColoring,rError,uMin,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringSetParameters.",err,error,*999)
    ENDIF

    EXITS("Petsc_MatFDColoringSetParameters")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringSetParameters",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringSetParameters

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringSetFunction routine.
  SUBROUTINE Petsc_MatFDColoringSetFunction(matFDColoring,fFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The matrix finite difference coloring to set
    EXTERNAL fFunction !<The external function to call
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringSetFunction",err,error,*999)

    CALL MatFDColoringSetFunction(matFDColoring%matFDColoring,fFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringSetFunction.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringSetFunction")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringSetFunction",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringSetFunction
        
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc MatFDColoringSetup routine.
  SUBROUTINE Petsc_MatFDColoringSetup(a,isColoring,matFDColoring,err,error,*)

    !Argument Variables
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The PETSc matrix to setup the FD coloring for
    TYPE(PetscISColoringType), INTENT(INOUT) :: isColoring !<The index set coloring to setup the finite difference coloring for
    TYPE(PetscMatFDColoringType), INTENT(INOUT) :: matFDColoring !<The matrix finite difference coloring to setup
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_MatFDColoringSetup",err,error,*999)

    CALL MatFDColoringSetup(A%mat,isColoring%isColoring,matFDColoring%matFDColoring,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in MatFDColoringSetup.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_MatFDColoringSetup")
    RETURN
999 ERRORSEXITS("Petsc_MatFDColoringSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_MatFDColoringSetup
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc PC structure 
  SUBROUTINE Petsc_PCFinalise(pc,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The PC to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_PCFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(pc%pc)) THEN      
      !Do nothing - should be destroyed when the KSP is destroyed.
      !CALL Petsc_PCDestroy(pc,err,error,*999)
    ENDIF
    
    EXITS("Petsc_PCFinalise")
    RETURN
999 ERRORSEXITS("Petsc_PCFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_PCFinalise
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc PC structure
  SUBROUTINE Petsc_PCInitialise(pc,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The PC to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_PCInitialise",err,error,*999)

    pc%pc%v=-2
   
    EXITS("Petsc_PCInitialise")
    RETURN
999 ERRORSEXITS("Petsc_PCInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_PCInitialise
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PCDestroy routine.
  SUBROUTINE Petsc_PCDestroy(pc,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The PC to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_PCDestroy",err,error,*999)

    CALL PCDestroy(pc%pc,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCDestroy.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_PCDestroy")
    RETURN
999 ERRORSEXITS("Petsc_PCDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_PCDestroy
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PCFactorGetMatrix routine.
  SUBROUTINE Petsc_PCFactorGetMatrix(pc,factoredMatrix,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The preconditioner to set the solver package for
    TYPE(PetscMatType), INTENT(OUT) :: factoredMatrix !<The factored matrix to get from preconditioner context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_PCFactorGetMatrix",err,error,*999)

    CALL PCFactorGetMatrix(pc%pc,factoredMatrix%mat,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCFactorGetMatrix",err,error,*999)
    ENDIF
    
    EXITS("Petsc_PCFactorGetMatrix")
    RETURN
999 ERRORSEXITS("Petsc_PCFactorGetMatrix",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_PCFactorGetMatrix

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PCFactoSetMatSolverType routine.
  SUBROUTINE Petsc_PCFactorSetMatSolverType(pc,matSolverType,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The preconditioner to set the solver package for
    MatSolverType, INTENT(IN) :: matSolverType !<The matrix solver type to set \see  OpenCMISSPETSc_PetscMatSolverTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_PCFactorSetMatSolverType",err,error,*999)

    CALL PCFactorSetMatSolverType(pc%pc,matSolverType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCFactorSetMatSolverType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_PCFactorSetMatSolverType")
    RETURN
999 ERRORSEXITS("Petsc_PCFactorSetMatSolverType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_PCFactorSetMatSolverType

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PCFactorSetUpMatSolverType routine.
  SUBROUTINE Petsc_PCFactorSetUpMatSolverType(pc,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The preconditioner to set the solver package for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_PCFactorSetUpMatSolverType",err,error,*999)

    CALL PCFactorSetUpMatSolverType(pc%pc,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCFactorSetUpMatSolverType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_PCFactorSetUpMatSolverType")
    RETURN
999 ERRORSEXITS("Petsc_PCFactorSetUpMatSolverType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_PCFactorSetUpMatSolverType

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PCSetReusePreconditioner routine.
  SUBROUTINE Petsc_PCSetReusePreconditioner(pc,flag,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The preconditioner to set the reuse for
    LOGICAL, INTENT(IN) :: flag !<True/false
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_PCSetReusePreconditioner",err,error,*999)

    IF(flag) THEN
      CALL PCSetReusePreconditioner(pc%pc,PETSC_TRUE,err)
    ELSE
      CALL PCSetReusePreconditioner(pc%pc,PETSC_FALSE,err)
    ENDIF
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCSetReusePreconditioner.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_PCSetReusePreconditioner")
    RETURN
999 ERRORSEXITS("Petsc_PCSetReusePreconditioner",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_PCSetReusePreconditioner
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PCSetFromOptions routine.
  SUBROUTINE Petsc_PCSetFromOptions(pc,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The preconditioner to set the options for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_PCSetFromOptions",err,error,*999)

    CALL PCSetFromOptions(pc%pc,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCSetFromOptions.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_PCSetFromOptions")
    RETURN
999 ERRORSEXITS("Petsc_PCSetFromOptions",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_PCSetFromOptions
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc PCSetType routine.
  SUBROUTINE Petsc_PCSetType(pc,method,err,error,*)

    !Argument Variables
    TYPE(PetscPCType), INTENT(INOUT) :: pc !<The preconditioner to set the type of
    PCType, INTENT(IN) :: method !<The preconditioning method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_PCSetType",err,error,*999)

    CALL PCSetType(pc%pc,method,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in PCSetType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_PCSetType")
    RETURN
999 ERRORSEXITS("Petsc_PCSetType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_PCSetType
    
  !
  !================================================================================================================================
  !

  !Finalise the PETSc SNES structure and destroy the SNES
  SUBROUTINE Petsc_SnesFinalise(snes,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(snes%snes)) THEN
      CALL Petsc_SnesDestroy(snes,err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesFinalise")
    RETURN
999 ERRORSEXITS("Petsc_SnesFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesFinalise
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc SNES structure
  SUBROUTINE Petsc_SnesInitialise(snes,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The snes to 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesInitialise",err,error,*999)

    snes%snes%v=-2
     
    EXITS("Petsc_SnesInitialise")
    RETURN
999 ERRORSEXITS("Petsc_SnesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesInitialise
    
!
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESComputeJacobianDefault routine.
  SUBROUTINE Petsc_SnesComputeJacobianDefault(snes,x,j,b,ctx,err,error,*)

    !Argument variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc SNES
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The PETSc x Vec
    TYPE(PetscMatType), INTENT(INOUT) :: j !<The PETSc J Mat
    TYPE(PetscMatType), INTENT(INOUT) :: b !<The PETSc B Mat
    TYPE(SolverType), POINTER :: ctx !<The passed through context
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesComputeJacobianDefault",err,error,*999)
    
    CALL SNESComputeJacobianDefault(snes%snes,x%vec,j%mat,b%mat,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESComputeJacobianDefault.",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesComputeJacobianDefault")    
    RETURN
999 ERRORSEXITS("Petsc_SnesComputeJacobianDefault",err,error)
    RETURN
    
  END SUBROUTINE Petsc_SnesComputeJacobianDefault

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESComputeJacobianDefaultColor routine.
  SUBROUTINE Petsc_SnesComputeJacobianDefaultColor(snes,x,j,b,ctx,err,error,*)

    !Argument variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc SNES
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The PETSc X Vec
    TYPE(PetscMatType), INTENT(INOUT) :: j !<The PETSc J Mat
    TYPE(PetscMatType), INTENT(INOUT) :: b !<The PETSc B Mat
    TYPE(PetscMatFDColoringType) :: ctx !<The passed through context
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesComputeJacobianDefaultColor",err,error,*999)
    
    CALL SNESComputeJacobianDefaultColor(snes%snes,x%vec,j%mat,b%mat,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESComputeJacobianDefaultColor.",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesComputeJacobianDefaultColor")
    RETURN
999 ERRORSEXITS("Petsc_SnesComputeJacobianDefaultColor",err,error)
    RETURN
    
  END SUBROUTINE Petsc_SnesComputeJacobianDefaultColor

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESCreate routine.
  SUBROUTINE Petsc_SnesCreate(communicator,snes,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator for the SNES creation
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator for the SNES creation
#endif    
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<On exit, the SNES information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesCreate",err,error,*999)

    CALL SNESCreate(communicator,snes%snes,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESCreate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesCreate")
    RETURN
999 ERRORSEXITS("Petsc_SnesCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesCreate
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESDestroy routine.
  SUBROUTINE Petsc_SnesDestroy(snes,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesDestroy",err,error,*999)

    CALL SNESDestroy(snes%snes,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESDestroy.",err,error,*999)
    ENDIF
     
    EXITS("Petsc_SnesDestroy")
    RETURN
999 ERRORSEXITS("Petsc_SnesDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesDestroy
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetApplicationContext routine.
  SUBROUTINE Petsc_SnesGetApplicationContext(snes,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the context for
    TYPE(SolverType), POINTER :: ctx !<On exit, the solver data context to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesGetApplicationContext",err,error,*999)

    IF(ASSOCIATED(ctx)) THEN
      CALL FlagError("Context is already associated.",err,error,*999)
    ELSE
      CALL SNESGetApplicationContext(snes%snes,ctx,err)
      IF(err/=0) THEN
        IF(petscHandleError) THEN
          CHKERRQ(err)
        ENDIF
        CALL FlagError("PETSc error in SNESGetApplicationContext.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Petsc_SnesGetApplicationContext")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetApplicationContext",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetApplicationContext

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetConvergedReason routine.
  SUBROUTINE Petsc_SnesGetConvergedReason(snes,reason,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the converged reason for
    INTEGER(INTG), INTENT(OUT) :: reason !<On exit, the converged reason
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesGetConvergedReason",err,error,*999)

    CALL SNESGetConvergedReason(snes%snes,reason,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetConvergedReason.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesGetConvergedReason")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetConvergedReason",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetConvergedReason
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetFunction routine.
  SUBROUTINE Petsc_SnesGetFunction(snes,f,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the function from
    TYPE(PetscVecType), INTENT(OUT) :: f !<The residual vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    ENTERS("Petsc_SnesGetFunction",err,error,*999)

    CALL SNESGetFunction(snes%snes,f%vec,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetFunction.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesGetFunction")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetFunction",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetFunction

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetIterationNumber routine.
  SUBROUTINE Petsc_SnesGetIterationNumber(snes,iterationNumber,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the iteration number for
    INTEGER(INTG), INTENT(OUT) :: iterationNumber !<On exit, the number of iterations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesGetIterationNumber",err,error,*999)

    CALL SNESGetIterationNumber(snes%snes,iterationNumber,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetIterationNumber.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesGetIterationNumber")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetIterationNumber",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetIterationNumber
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetJacobian routine for solver contexts.
  SUBROUTINE Petsc_SnesGetJacobianSolver(snes,a,b,jFunction,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the function for
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The Jacobian matrix
    TYPE(PetscMatType), INTENT(INOUT) :: b !<The Jacobian preconditioning matrix
    EXTERNAL jFunction !<The external function to call
!     TYPE(SolverType), POINTER :: CTX !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesGetJacobianSolver",err,error,*999)

    CALL SNESGetJacobian(snes%snes,a%mat,b%mat,jFunction,PETSC_NULL_INTEGER,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetJacobian.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesGetJacobianSolver")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetJacobianSolver",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetJacobianSolver

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetJacobian routine for solver contexts.
  SUBROUTINE Petsc_SnesGetJacobianSpecial(snes,a,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the function for
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The Jacobian matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesGetJacobianSpecial",err,error,*999)

    CALL SNESGetJacobian(snes%snes,a%mat,PETSC_NULL_MAT,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,err)

    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetJacobian.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesGetJacobianSpecial")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetJacobianSpecial",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetJacobianSpecial

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetKSP routine.
  SUBROUTINE Petsc_SnesGetKSP(snes,ksp,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the iteration number for
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<On exit, the KSP to associated with the SNES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesGetKSP",err,error,*999)

    CALL SNESGetKSP(snes%snes,ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetKSP.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesGetKSP")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetKSP",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetKSP

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetLineSearch routine.
  SUBROUTINE Petsc_SnesGetLineSearch(snes,lineSearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the SNES line search for
    TYPE(PetscSnesLinesearchType), INTENT(OUT) :: lineSearch !<On return, the SNES line search object
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesGetLineSearch",err,error,*999)

    CALL SNESGetLineSearch(snes%snes,lineSearch%snesLineSearch,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetLineSearch.",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesGetLineSearch")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetLineSearch",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetLineSearch

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSC SNESGetSolutionUpdate routine.
  SUBROUTINE Petsc_SnesGetSolutionUpdate(snes,solutionUpdate,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to get the solution update for
    TYPE(PetscVecType), INTENT(INOUT) :: solutionUpdate !<On exit, the solution update
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesGetSolutionUpdate",err,error,*999)

    CALL SNESGetSolutionUpdate(snes%snes,solutionUpdate%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESGetSolutionUpdate.",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesGetSolutionUpdate")
    RETURN
999 ERRORSEXITS("Petsc_SnesGetSolutionUpdate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesGetSolutionUpdate
    
!
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESMonitorSet routine.
  SUBROUTINE Petsc_SnesMonitorSet(snes,mFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set from the command line options
    EXTERNAL :: mFunction !<The external monitor function to set
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the monitor function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesMonitorSet",err,error,*999)

    CALL SNESMonitorSet(snes%snes,mFunction,ctx,PETSC_NULL_FUNCTION,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESMonitorSet.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesMonitorSet")
    RETURN
999 ERRORSEXITS("Petsc_SnesMonitorSet",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesMonitorSet

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESQNSetRestartType routine.
  SUBROUTINE Petsc_SnesQNSetRestartType(snes,rType,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the type for
    SNESQNRestartType, INTENT(IN) :: rType !<The SNES QN restart type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesQNSetRestartType",err,error,*999)

    CALL SNESQNSetRestartType(snes%snes,rType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESQNSetRestartType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesQNSetRestartType")
    RETURN
999 ERRORSEXITS("Petsc_SnesQNSetRestartType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesQNSetRestartType

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESQNSetScaleType routine.
  SUBROUTINE Petsc_SnesQNSetScaleType(snes,sType,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the type for
    SNESQNScaleType, INTENT(IN) :: sType !<The SNES QN scaling type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesQNSetScaleType",err,error,*999)

    CALL SNESQNSetScaleType(snes%snes,sType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESQNSetScaleType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesQNSetScaleType")
    RETURN
999 ERRORSEXITS("Petsc_SnesQNSetScaleType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesQNSetScaleType

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESQNSetType routine.
  SUBROUTINE Petsc_SnesQNSetType(snes,qType,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the type for
    SNESQNType, INTENT(IN) :: qType !<The SNES QN type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesQNSetType",err,error,*999)

    CALL SNESQNSetType(snes%snes,qType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESQNSetType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesQNSetType")
    RETURN
999 ERRORSEXITS("Petsc_SnesQNSetType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesQNSetType

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetApplicationContext routine.
  SUBROUTINE Petsc_SnesSetApplicationContext(snes,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the context for
    TYPE(SolverType), POINTER :: ctx !<The solver data context to set as a context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSetApplicationContext",err,error,*999)

    CALL SNESSetApplicationContext(snes%snes,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetApplicationContext.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSetApplicationContext")
    RETURN
999 ERRORSEXITS("Petsc_SnesSetApplicationContext",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSetApplicationContext

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetFunction routine.
  SUBROUTINE Petsc_SnesSetConvergenceTest(snes,cFunction,ctx,err,error,*)
    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the function for
    EXTERNAL cFunction !<The external function to call (OpenCMISS subroutine to calculate convergence)
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the convergence test function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSetConvergenceTest",err,error,*999)

    CALL SNESSetConvergenceTest(snes%snes,cFunction,ctx,PETSC_NULL_FUNCTION,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetConvergenceTest.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSetConvergenceTest")
    RETURN
999 ERRORSEXITS("Petsc_SnesSetConvergenceTest",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSetConvergenceTest

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetFromOptions routine.
  SUBROUTINE Petsc_SnesSetFromOptions(snes,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set from the command line options
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSetFromOptions",err,error,*999)

    CALL SNESSetFromOptions(snes%snes,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetFromOptions.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSetFromOptions")
    RETURN
999 ERRORSEXITS("Petsc_SnesSetFromOptions",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSetFromOptions
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetFunction routine.
  SUBROUTINE Petsc_SnesSetFunction(snes,f,fFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the function for
    TYPE(PetscVecType), INTENT(INOUT) :: f !<The residual vector
    EXTERNAL fFunction !<The external function to call
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSetFunction",err,error,*999)

    CALL SNESSetFunction(snes%snes,f%vec,fFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetFunction.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSetFunction")
    RETURN
999 ERRORSEXITS("Petsc_SnesSetFunction",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSetFunction

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetJacobian routine for solver contexts.
  SUBROUTINE Petsc_SnesSetJacobianSolver(snes,a,b,jFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the function for
    TYPE(PetscMatType), INTENT(INOUT) :: a !<The Jacobian matrix
    TYPE(PetscMatType), INTENT(INOUT) :: b !<The Jacobian preconditioning matrix
    EXTERNAL jFunction !<The external function to call
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSetJacobianSolver",err,error,*999)

    CALL SNESSetJacobian(snes%snes,a%mat,b%mat,jFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetJacobian.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSetJacobianSolver")
    RETURN
999 ERRORSEXITS("Petsc_SnesSetJacobianSolver",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSetJacobianSolver

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESGetKSP routine.
  SUBROUTINE Petsc_SnesSetKsp(snes,ksp,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the KSP for
    TYPE(PetscKspType), INTENT(INOUT) :: ksp !<The KSP to be associated with the SNES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSetKsp",err,error,*999)

    CALL SNESSetKSP(snes%snes,ksp%ksp,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetKSP.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSetKsp")
    RETURN
999 ERRORSEXITS("Petsc_SnesSetKsp",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSetKsp

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetNormSchedule routine.
  SUBROUTINE Petsc_SnesSetNormSchedule(snes,normSchedule,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the norm type for
    SNESNormSchedule, INTENT(IN) :: normSchedule !<The norm schedule
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSetNormSchedule",err,error,*999)

    CALL SNESSetNormSchedule(snes%snes,normSchedule,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetNormSchedule.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSetNormSchedule")
    RETURN
999 ERRORSEXITS("Petsc_SnesSetNormSchedule",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSetNormSchedule

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetTolerances routine.
  SUBROUTINE Petsc_SnesSetTolerances(snes,absTol,rTol,sTol,maxIterations,maxFunctionEvals,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the tolerances for
    REAL(DP), INTENT(IN) :: absTol !<The absolute convergence tolerance
    REAL(DP), INTENT(IN) :: rTol !<The relative convergence tolerance
    REAL(DP), INTENT(IN) :: sTol !<The convergence tolerance for the change in the solution between steps
    INTEGER(INTG), INTENT(IN) :: maxIterations !<The maximum number of iterations
    INTEGER(INTG), INTENT(IN) :: maxFunctionEvals !<The maximum number of function evaluations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSetTolerances",err,error,*999)

    CALL SNESSetTolerances(snes%snes,absTol,rTol,sTol,maxIterations,maxFunctionEvals,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetTolerances.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSetTolerances")
    RETURN
999 ERRORSEXITS("Petsc_SnesSetTolerances",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSetTolerances
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetTrustRegionTolerance routine.
  SUBROUTINE Petsc_SnesSetTrustRegionTolerance(snes,trTol,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the tolerances for
    REAL(DP), INTENT(IN) :: trTol !<The trust region tolerance
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSetTrustRegionTolerance",err,error,*999)

    CALL SNESSetTrustRegionTolerance(snes%snes,trTol,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetTrustRegionTolerance.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSetTrustRegionTolerance")
    RETURN
999 ERRORSEXITS("Petsc_SnesSetTrustRegionTolerance",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSetTrustRegionTolerance
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSetType routine.
  SUBROUTINE Petsc_SnesSetType(snes,method,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to set the type for
    SNESType, INTENT(IN) :: method !<The SNES type
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSetType",err,error,*999)

    CALL SNESSetType(snes%snes,method,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSetType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSetType")
    RETURN
999 ERRORSEXITS("Petsc_SnesSetType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSetType

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESSolve routine.
  SUBROUTINE Petsc_SnesSolve(snes,b,x,err,error,*)

    !Argument Variables
    TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The SNES to solve
    TYPE(PetscVecType), INTENT(INOUT) :: b !<The constant part of the equation
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The solution vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesSolve",err,error,*999)

    CALL SNESSolve(snes%snes,b%vec,x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESSolve.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesSolve")
    RETURN
999 ERRORSEXITS("Petsc_SnesSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesSolve

  !
  !================================================================================================================================
  !

  !Finalise the PETSc SNES line search structure
  SUBROUTINE Petsc_SnesLineSearchFinalise(lineSearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesLineSearchFinalise",err,error,*999)

    !We don't actually call PETSc's SNESLineSearchDestroy as PETSc accesses and destroys the LineSearch when calling
    !SNESDestroy, so we'll just let PETSc clean it up.

    EXITS("Petsc_SnesLineSearchFinalise")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesLineSearchFinalise

  !
  !================================================================================================================================
  !

  !Initialise the PETSc SNES line search structure
  SUBROUTINE Petsc_SnesLineSearchInitialise(lineSearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesLineSearchInitialise",err,error,*999)
    
    lineSearch%snesLineSearch%v=-2
    
    EXITS("Petsc_SnesLineSearchInitialise")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchInitialise",err,error)
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchInitialise

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESLineSearchBTSetAlpha routine.
  SUBROUTINE Petsc_SnesLineSearchBTSetAlpha(lineSearch,alpha,err,error,*)

    !Argument variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES back-tracking LineSearch to set alpha for
    REAL(DP), INTENT(IN) :: alpha !<The alpha descent parameter to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesLineSearchBTSetAlpha",err,error,*999)

    CALL SNESLineSearchBTSetAlpha(lineSearch%snesLineSearch,alpha,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      END IF
      CALL FlagError("PETSc error in SNESLineSearchBTSetAlpha.",err,error,*999)
    END IF

    EXITS("Petsc_SnesLineSearchBTSetAlpha")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchBTSetAlpha",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesLineSearchBTSetAlpha

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESLineSearchComputeNorms routine.
  SUBROUTINE Petsc_SnesLineSearchComputeNorms(lineSearch,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to compute norms for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesLineSearchComputeNorms",err,error,*999)

    CALL SnesLineSearchComputeNorms(lineSearch%snesLineSearch,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SnesLineSearchComputeNorms.",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesLineSearchComputeNorms")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchComputeNorms",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesLineSearchComputeNorms

  !
  !================================================================================================================================
  !

  !>Buffer routine to the petsc SnesLineSearchGetNorms routine.
  SUBROUTINE Petsc_SnesLineSearchGetNorms(lineSearch,xNorm,fNorm,yNorm,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to get the norms for X, Y, and F from.
    REAL(DP), INTENT(INOUT) :: xNorm !<On exit, the norm of the current solution
    REAL(DP), INTENT(INOUT) :: fNorm !<On exit, the norm of the current function
    REAL(DP), INTENT(INOUT) :: yNorm !<On exit, the norm of the current update
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesLineSearchGetNorms",err,error,*999)

    CALL SnesLineSearchGetNorms(lineSearch%snesLineSearch,xNorm,fNorm,yNorm,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SnesLineSearchGetNorms.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesLineSearchGetNorms")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchGetNorms",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesLineSearchGetNorms

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESLineSearchGetVecs routine.
  SUBROUTINE Petsc_SnesLineSearchGetVecs(lineSearch,x,f,y,w,g,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The PetcsSnesLineSearch to get the vectors from the SNESLineSearch
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The The old solution 
    TYPE(PetscVecType), INTENT(INOUT) :: f !<The old function 
    TYPE(PetscVecType), INTENT(INOUT) :: y !<The search direction 
    TYPE(PetscVecType), INTENT(INOUT) :: w !<The new solution 
    TYPE(PetscVecType), INTENT(INOUT) :: g !<The new function 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesLineSearchGetVecs",err,error,*999)

    CALL SNESLineSearchGetVecs(lineSearch%snesLineSearch,x%vec,f%vec,y%vec,w%vec,g%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchGetVecs.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesLineSearchGetVecs")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchGetVecs",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesLineSearchGetVecs

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESLineSearchMonitorSet routine.
  SUBROUTINE Petsc_SnesLineSearchMonitorSet(lineSearch,mFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set the linesearch monitor routine for
    EXTERNAL :: mFunction !<The external monitor function to set
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the monitor function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesLineSearchMonitorSet",err,error,*999)

    !CALL SnesLineSearchMonitorSet(lineSearch%snesLineSearch,mFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchMonitorSet.",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesLineSearchMonitorSet")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchMonitorSet",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesLineSearchMonitorSet

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESLineSearchSetComputeNorms routine.
  SUBROUTINE Petsc_SnesLineSearchSetComputeNorms(lineSearch,computeNorms,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set whether to compute norms
    LOGICAL, INTENT(IN) :: computeNorms !<Whether to compute norms
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesLineSearchSetComputeNorms",err,error,*999)

    IF(computeNorms) THEN
      CALL SNESLineSearchSetComputeNorms(lineSearch%snesLineSearch,PETSC_TRUE,err)
    ELSE
      CALL SNESLineSearchSetComputeNorms(lineSearch%snesLineSearch,PETSC_FALSE,err)      
    ENDIF
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetComputeNorms.",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesLineSearchSetComputeNorms")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchSetComputeNorms",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesLineSearchSetComputeNorms

  !
  !================================================================================================================================
  !

  !>Buffer routine to the petsc SnesLineSearchSetNorms routine.
  SUBROUTINE Petsc_SnesLineSearchSetNorms(lineSearch,xNorm,fNorm,yNorm,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNESLineSearch to set the computed norms for x, y, and f
    REAL(DP), INTENT(IN) :: xNorm !<The norm of the current solution
    REAL(DP), INTENT(IN) :: fNorm !<The norm of the current function
    REAL(DP), INTENT(IN) :: yNorm !<The norm of the current update
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SnesLineSearchSetNorms",err,error,*999)

    CALL SnesLineSearchSetNorms(lineSearch%snesLineSearch,xNorm,fNorm,yNorm,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("petsc error in SnesLineSearchSetNorms.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SnesLineSearchSetNorms")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchSetNorms",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesLineSearchSetNorms

  !
  !================================================================================================================================
  !
  !>Buffer routine to the PETSc SNESLineSearchSetOrder routine.
  SUBROUTINE Petsc_SnesLineSearchSetOrder(lineSearch,lineSearchOrder,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set the line search order for
    SNESLineSearchOrder, INTENT(IN) :: lineSearchOrder !<The line search order to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Petsc_SnesLineSearchSetOrder",err,error,*999)

    SELECT CASE(lineSearchOrder)
    CASE(PETSC_SNES_LINESEARCH_ORDER_LINEAR)
      CALL SNESLineSearchSetOrder(lineSearch%snesLineSearch,SNES_LINESEARCH_ORDER_LINEAR,err)
    CASE(PETSC_SNES_LINESEARCH_ORDER_QUADRATIC)
      CALL SNESLineSearchSetOrder(lineSearch%snesLineSearch,SNES_LINESEARCH_ORDER_QUADRATIC,err)
    CASE(PETSC_SNES_LINESEARCH_ORDER_CUBIC)
      CALL SNESLineSearchSetOrder(lineSearch%snesLineSearch,SNES_LINESEARCH_ORDER_CUBIC,err)
    CASE DEFAULT
      localError="The specified line search order of "//TRIM(NumberToVString(lineSearchOrder,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetOrder.",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesLineSearchSetOrder")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchSetOrder",err,error)
    RETURN 1
  END SUBROUTINE Petsc_SnesLineSearchSetOrder

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc SNESLineSearchSetTolerances routine.
  SUBROUTINE Petsc_SnesLineSearchSetTolerances(lineSearch,steptol,maxstep,rtol,atol,ltol,maxIt,err,error,*)

    !Argument variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set tolerances for
    REAL(DP), INTENT(IN) :: steptol !<The minimum steplength
    REAL(DP), INTENT(IN) :: maxstep !<The maximum steplength
    REAL(DP), INTENT(IN) :: rtol !<The relative tolerance for iterative line searches
    REAL(DP), INTENT(IN) :: atol !<The absolute tolerance for iterative line searches
    REAL(DP), INTENT(IN) :: ltol !<The change in lambda tolerance for iterative line searches
    INTEGER(INTG), INTENT(IN) :: maxIt !<The maximum number of iterations of the line search
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesLineSearchSetTolerances",err,error,*999)

    CALL SNESLineSearchSetTolerances(lineSearch%snesLineSearch,steptol,maxstep,rtol,atol,ltol,maxIt,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      END IF
      CALL FlagError("PETSc error in SNESLineSearchSetTolerances.",err,error,*999)
    END IF

    EXITS("Petsc_SnesLineSearchSetTolerances")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchSetTolerances",err,error)
    RETURN 1

  END SUBROUTINE Petsc_SnesLineSearchSetTolerances

  !
  !================================================================================================================================
  !
 
  !>Buffer routine to the PETSc SNESLineSearchSetType routine.
  SUBROUTINE Petsc_SnesLineSearchSetType(lineSearch,lineSearchType,err,error,*)

    !Argument Variables
    TYPE(PetscSnesLineSearchType), INTENT(INOUT) :: lineSearch !<The SNES LineSearch to set the line search type for
    SNESLineSearchType, INTENT(IN) :: lineSearchType !<The line search type to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Petsc_SnesLineSearchSetType",err,error,*999)

    CALL SNESLineSearchSetType(lineSearch%snesLineSearch,lineSearchType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in SNESLineSearchSetType.",err,error,*999)
    ENDIF

    EXITS("Petsc_SnesLineSearchSetType")
    RETURN
999 ERRORSEXITS("Petsc_SnesLineSearchSetType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SnesLineSearchSetType

  !
  !================================================================================================================================
  !

  !Finalise the PETSc TAO structure and destroy the TAO
  SUBROUTINE Petsc_TaoFinalise(tao,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The Tao to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(tao%tao)) THEN
      CALL Petsc_TaoDestroy(tao,err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoFinalise")
    RETURN
999 ERRORSEXITS("Petsc_TaoFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoFinalise
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc Tao structure
  SUBROUTINE Petsc_TaoInitialise(tao,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The Tao to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoInitialise",err,error,*999)

    tao%tao%v=-2
     
    EXITS("Petsc_TaoInitialise")
    RETURN
999 ERRORSEXITS("Petsc_TaoInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoInitialise
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TaoCreate routine.
  SUBROUTINE Petsc_TaoCreate(communicator,tao,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(INOUT) :: communicator !<The MPI communicator for the Tao creation
#else    
    MPI_Comm, INTENT(INOUT) :: communicator !<The MPI communicator for the Tao creation
#endif    
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<On exit, the Tao information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoCreate",err,error,*999)

    CALL TaoCreate(communicator,tao%tao,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoCreate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoCreate")
    RETURN
999 ERRORSEXITS("Petsc_TaoCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoCreate
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc TaoDestroy routine.
  SUBROUTINE Petsc_TaoDestroy(tao,err,error,*)

    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The Tao to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoDestroy",err,error,*999)

    CALL TaoDestroy(tao%tao,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoDestroy.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoDestroy")
    RETURN
999 ERRORSEXITS("Petsc_TaoDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoDestroy
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc TaoGetConvergedReason routine
  SUBROUTINE Petsc_TaoGetConvergedReason(tao,reason,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The Tao information
    INTEGER(INTG), INTENT(OUT) :: reason !<On exit, the converged reason
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoGetConvergedReason",err,error,*999)

    CALL TaoGetConvergedReason(tao%tao,reason,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoGetConvergedReason.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoGetConvergedReason")
    RETURN
999 ERRORSEXITS("Petsc_TaoGetConvergedReason",err,error)
    RETURN 1
  END SUBROUTINE Petsc_TaoGetConvergedReason
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc TaoGetTolerances routine
  SUBROUTINE Petsc_TaoGetTolerances(tao,gaTol,grTol,gtTol,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The Tao to get the tolerances for
    REAL(DP), INTENT(OUT) :: gaTol !<On exit, the absolute tolerance of the gradient norm
    REAL(DP), INTENT(OUT) :: grTol !<On exit, the relative tolerance of the gradient norm
    REAL(DP), INTENT(OUT) :: gtTol !<On exit, the factor tolerance of the gradient norm
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoGetTolerances",err,error,*999)

    CALL TaoGetTolerances(tao%tao,gaTol,grTol,gtTol,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoGetTolerances.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoGetTolerances")
    RETURN
999 ERRORSEXITS("Petsc_TaoGetTolerances",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoGetTolerances
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc TaoSetFromOptions routine.
  SUBROUTINE Petsc_TaoSetFromOptions(tao,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The tao to set the options for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoSetFromOptions",err,error,*999)

    CALL TaoSetFromOptions(tao%tao,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoSetFromOptions.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoSetFromOptions")
    RETURN
999 ERRORSEXITS("Petsc_TaoSetFromOptions",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoSetFromOptions
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc TaoSetInequalityBounds routine.
  SUBROUTINE Petsc_TaoSetInequalityBounds(tao,lowerBounds,upperBounds,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The tao to set the inequality bounds for
    TYPE(PetscVecType), INTENT(INOUT) :: lowerBounds !<The vector of lower bounds
    TYPE(PetscVecType), INTENT(INOUT) :: upperBounds !<The vector of upper bounds    
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoSetInequalityBounds",err,error,*999)

    CALL TaoSetInequalityBounds(tao%tao,lowerBounds%vec,upperBounds%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoSetInequalityBounds.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoSetInequalityBounds")
    RETURN
999 ERRORSEXITS("Petsc_TaoSetInequalityBounds",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoSetInequalityBounds
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc TaoSetGradient routine.
  SUBROUTINE Petsc_TaoSetGradient(tao,gFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The tao to set the gradient function for
    EXTERNAL gFunction !<The external function to call
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoSetGradient",err,error,*999)

    CALL TaoSetGradient(tao%tao,gFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoSetGradient.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoSetGradient")
    RETURN
999 ERRORSEXITS("Petsc_TaoSetGradient",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoSetGradient
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc TaoSetHessian routine.
  SUBROUTINE Petsc_TaoSetHessian(tao,H,Hpre,hFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The tao to set the Hessian function for
    TYPE(PetscMatType), INTENT(INOUT) :: H !<The matrix used for the Hessian
    TYPE(PetscMatType), INTENT(INOUT) :: Hpre !<The matrix used for the Hessian preconditioner
    EXTERNAL hFunction !<The external Hessian function to call
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoSetHessian",err,error,*999)

    CALL TaoSetHessian(tao%tao,H%mat,Hpre%mat,hFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoSetHessian.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoSetHessian")
    RETURN
999 ERRORSEXITS("Petsc_TaoSetHessian",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoSetHessian
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc TaoMonitorSet routine.
  SUBROUTINE Petsc_TaoMonitorSet(tao,mFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The tao to set the monitor function for
    EXTERNAL mFunction !<The external monitor function to call
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoMonitorSet",err,error,*999)

    CALL TaoMonitorSet(tao%tao,mFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoMonitorSet.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoMonitorSet")
    RETURN
999 ERRORSEXITS("Petsc_TaoMonitorSet",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoMonitorSet
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc TaoSetObjective routine.
  SUBROUTINE Petsc_TaoSetObjective(tao,oFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The tao to set the objective function for
    EXTERNAL oFunction !<The external objective function to call
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoSetObjective",err,error,*999)

    CALL TaoSetObjective(tao%tao,oFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoSetObjective.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoSetObjective")
    RETURN
999 ERRORSEXITS("Petsc_TaoSetObjective",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoSetObjective
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc TaoSetObjectiveAndGradient routine.
  SUBROUTINE Petsc_TaoSetObjectiveAndGradient(tao,ogFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The tao to set the objective and gradient function for
    EXTERNAL ogFunction !<The external objective and gradient function to call
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoSetObjectiveAndGradient",err,error,*999)

    CALL TaoSetObjectiveAndGradient(tao%tao,ogFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoSetObjectiveAndGradient.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoSetObjectiveAndGradient")
    RETURN
999 ERRORSEXITS("Petsc_TaoSetObjectiveAndGradient",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoSetObjectiveAndGradient
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc TaoSetSolution routine.
  SUBROUTINE Petsc_TaoSetSolution(tao,x0,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The tao to set the solution vector for
    TYPE(PetscVecType), INTENT(INOUT) :: x0 !<The solution vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoSetSolution",err,error,*999)

    CALL TaoSetSolution(tao%tao,x0%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoSetSolution.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoSetSolution")
    RETURN
999 ERRORSEXITS("Petsc_TaoSetSolution",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoSetSolution
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc TaoSetTolerances routine
  SUBROUTINE Petsc_TaoSetTolerances(tao,gaTol,grTol,gtTol,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The Tao to set the tolerances for
    REAL(DP), INTENT(IN) :: gaTol !<The absolute tolerance of the gradient norm to set
    REAL(DP), INTENT(IN) :: grTol !<The relative tolerance of the gradient norm to set
    REAL(DP), INTENT(IN) :: gtTol !<The factor tolerance of the gradient norm to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoSetTolerances",err,error,*999)

    CALL TaoSetTolerances(tao%tao,gaTol,grTol,gtTol,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoSetTolerances.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoSetTolerances")
    RETURN
999 ERRORSEXITS("Petsc_TaoSetTolerances",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoSetTolerances
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc TaoSetVariableBounds routine.
  SUBROUTINE Petsc_TaoSetVariableBounds(tao,lowerBounds,upperBounds,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The tao to set the variable bounds for
    TYPE(PetscVecType), INTENT(INOUT) :: lowerBounds !<The vector of lower bounds
    TYPE(PetscVecType), INTENT(INOUT) :: upperBounds !<The vector of upper bounds    
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoSetVariableBounds",err,error,*999)

    CALL TaoSetVariableBounds(tao%tao,lowerBounds%vec,upperBounds%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoSetVariableBounds.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoSetVariableBounds")
    RETURN
999 ERRORSEXITS("Petsc_TaoSetVariableBounds",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoSetVariableBounds
    
  !
  !================================================================================================================================
  !
  
  !>Buffer routine to the PETSc TaoSolve routine.
  SUBROUTINE Petsc_TaoSolve(tao,err,error,*)

    !Argument Variables
    TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The tao to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TaoSolve",err,error,*999)

    CALL TaoSolve(tao%tao,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TaoSolve.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TaoSolve")
    RETURN
999 ERRORSEXITS("Petsc_TaoSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TaoSolve
    
  !
  !================================================================================================================================
  !
    
  !Finalise the PETSc TS structure and destroy the TS
  SUBROUTINE Petsc_TSFinalise(ts,err,error,*)

    !Argument Variables
    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(ts%ts)) THEN
      CALL Petsc_TSDestroy(ts,err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSFinalise")
    RETURN
999 ERRORSEXITS("Petsc_TSFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSFinalise
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc TS structure
  SUBROUTINE Petsc_TSInitialise(ts,err,error,*)

    !Argument Variables
    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSInitialise",err,error,*999)

    ts%ts%v=-2
     
    EXITS("Petsc_TSInitialise")
    RETURN
999 ERRORSEXITS("Petsc_TSInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSInitialise
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSCreate routine.
  SUBROUTINE Petsc_TSCreate(communicator,ts,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
   TYPE(MPI_Comm), INTENT(INOUT) :: communicator !<The MPI communicator for the TS creation
#else    
    MPI_Comm, INTENT(INOUT) :: communicator !<The MPI communicator for the TS creation
#endif    
    TYPE(PetscTSType), INTENT(INOUT) :: ts !<On exit, the TS information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSCreate",err,error,*999)

    CALL TSCreate(communicator,ts%ts,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSCreate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSCreate")
    RETURN
999 ERRORSEXITS("Petsc_TSCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSCreate
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSDestroy routine.
  SUBROUTINE Petsc_TSDestroy(ts,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSDestroy",err,error,*999)

    CALL TSDestroy(ts%ts,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSDestroy.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSDestroy")
    RETURN
999 ERRORSEXITS("Petsc_TSDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSDestroy
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSGetSolution routine.
  SUBROUTINE Petsc_TSGetSolution(ts,currentSolution,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the time step for
    TYPE(PetscVecType), INTENT(INOUT) :: currentSolution !<The current solution to be set for the TS
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSGetSolution",err,error,*999)

    CALL TSGetSolution(ts%ts,currentSolution%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSGetSolution.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSGetSolution")
    RETURN
999 ERRORSEXITS("Petsc_TSGetSolution",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSGetSolution
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc TSMonitorSet routine.
  SUBROUTINE Petsc_TSMonitorSet(ts,mFunction,ctx,err,error,*)

    !Argument Variables
    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the monitor for
    EXTERNAL :: mFunction !<The external monitor function to set
    TYPE(SolverType), POINTER :: ctx !<The solver data to pass to the monitor function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSMonitorSet",err,error,*999)

    CALL TSMonitorSet(ts%ts,mFunction,ctx,PETSC_NULL_FUNCTION,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSMonitorSet.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSMonitorSet")
    RETURN
999 ERRORSEXITS("Petsc_TSMonitorSet",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSMonitorSet
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetDuration routine.
  SUBROUTINE Petsc_TSSetDuration(ts,maxSteps,maxTime,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set from the options
    INTEGER(INTG), INTENT(IN) :: maxSteps !<The maximum number of steps to use
    REAL(DP), INTENT(IN) :: maxTime !<The maximum time to iteration to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSetDuration",err,error,*999)

    CALL TSSetDuration(ts%ts,maxSteps,maxTime,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetDuration.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSSetDuration")
    RETURN
999 ERRORSEXITS("Petsc_TSSetDuration",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSetDuration
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetExactFinalTime routine.
  SUBROUTINE Petsc_TSSetExactFinalTime(ts,exactFinalTimeOption,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the initial time step for
    TSExactFinalTimeOption, INTENT(IN) :: exactFinalTimeOption !<The option for exact final time to set \see OpenCMISSPETSc_PetscTSExactFinalTimeOptions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSetExactFinalTime",err,error,*999)

    CALL TSSetExactFinalTime(ts%ts,exactFinalTimeOption,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetExactFinalTime.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSSetExactFinalTime")
    RETURN
999 ERRORSEXITS("Petsc_TSSetExactFinalTime",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSetExactFinalTime

  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetFromOptions routine.
  SUBROUTINE Petsc_TSSetFromOptions(ts,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set from the options
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSetFromOptions",err,error,*999)

    CALL TSSetFromOptions(ts%ts,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetFromOptions.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSSetFromOptions")
    RETURN
999 ERRORSEXITS("Petsc_TSSetFromOptions",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSetFromOptions

  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetMaxTime routine.
  SUBROUTINE Petsc_TSSetMaxTime(ts,maxTime,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the time for
    REAL(DP), INTENT(IN) :: maxTime !<The maximum time to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSetMaxTime",err,error,*999)

    CALL TSSetMaxTime(ts%ts,maxTime,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetMaxTime.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSSetMaxTime")
    RETURN
999 ERRORSEXITS("Petsc_TSSetMaxTime",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSetMaxTime

  !
  !================================================================================================================================
  !
      
  !>Buffer routine to the PETSc TSSetProblemType routine.
  SUBROUTINE Petsc_TSSetProblemType(ts,probType,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the problem type for
    INTEGER(INTG), INTENT(IN) :: probType !<The problem type to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSetProblemType",err,error,*999)

    CALL TSSetProblemType(ts%ts,probType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetProblemType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSSetProblemType")
    RETURN
999 ERRORSEXITS("Petsc_TSSetProblemType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSetProblemType
    
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetRHSFunction routine.
  SUBROUTINE Petsc_TSSetRHSFunction(ts,rates,rhsFunction,ctx,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the problem type for
    TYPE(PetscVecType), INTENT(INOUT) :: rates
    EXTERNAL rhsFunction !<The external RHS function to call
    TYPE(CellMLPETScContextType), POINTER :: ctx !<The solver data to pass to the function
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSetRHSFunction",err,error,*999)

    CALL TSSetRHSFunction(ts%ts,rates%vec,rhsFunction,ctx,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetRHSFunction.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSSetRHSFunction")
    RETURN
999 ERRORSEXITS("Petsc_TSSetRHSFunction",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSetRHSFunction

  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetSolution routine.
  SUBROUTINE Petsc_TSSetSolution(ts,initialSolution,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the time step for
    TYPE(PetscVecType), INTENT(IN) :: initialSolution !<The initial solution to be set for the TS
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSetSolution",err,error,*999)

    CALL TSSetSolution(ts%ts,initialSolution%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetSolution.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSSetSolution")
    RETURN
999 ERRORSEXITS("Petsc_TSSetSolution",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSetSolution

  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetTime routine.
  SUBROUTINE Petsc_TSSetTime(ts,time,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the time for
    REAL(DP), INTENT(IN) :: time !<The time to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSetTime",err,error,*999)

    CALL TSSetTime(ts%ts,time,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetTime.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSSetTime")
    RETURN
999 ERRORSEXITS("Petsc_TSSetTime",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSetTime

  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetTimeStep routine.
  SUBROUTINE Petsc_TSSetTimeStep(ts,timeStep,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the time step for
    REAL(DP), INTENT(IN) :: timeStep !<The time step to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSetTimeStep",err,error,*999)

    CALL TSSetTimeStep(ts%ts,timeStep,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetTimeStep.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSSetTimeStep")
    RETURN
999 ERRORSEXITS("Petsc_TSSetTimeStep",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSetTimeStep
  
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSetType routine.
  SUBROUTINE Petsc_TSSetType(ts,method,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to set the type for
    TSType, INTENT(IN) :: method !<The time stepping method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSetType",err,error,*999)

    CALL TSSetType(ts%ts,method,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSetType.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSSetType")
    RETURN
999 ERRORSEXITS("Petsc_TSSetType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSetType
  
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSolve routine.
  SUBROUTINE Petsc_TSSolve(ts,x,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to solve
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The solution vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSolve",err,error,*999)

    CALL TSSolve(ts%ts,x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSolve.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSSolve")
    RETURN
999 ERRORSEXITS("Petsc_TSSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSolve
  
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSStep routine.
  SUBROUTINE Petsc_TSStep(ts,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to step
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSStep",err,error,*999)

    CALL TSStep(ts%ts,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSStep.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_TSStep")
    RETURN
999 ERRORSEXITS("Petsc_TSStep",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSStep

  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSundialsSetType routine.
  SUBROUTINE Petsc_TSSundialsSetType(ts,sundialsType,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to step
    TSSundialsType, INTENT(IN) :: sundialsType !<The type of Sundials solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSundialsSetType",err,error,*999)

#ifdef WITH_SUNDIALS    
    CALL TSSundialsSetType(ts%ts,sundialsType,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSundialsSetType.",err,error,*999)
    ENDIF
#endif    
    
    EXITS("Petsc_TSSundialsSetType")
    RETURN
999 ERRORSEXITS("Petsc_TSSundialsSetType",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSundialsSetType
  !
  !================================================================================================================================
  !
    
  !>Buffer routine to the PETSc TSSundialsSetTolerance routine.
  SUBROUTINE Petsc_TSSundialsSetTolerance(ts,absTol,relTol,err,error,*)

    TYPE(PetscTSType), INTENT(INOUT) :: ts !<The TS to step
    REAL(DP), INTENT(IN) :: absTol !<The absolute tolerance
    REAL(DP), INTENT(IN) :: relTol !<The relative tolerance
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_TSSundialsSetTolerance",err,error,*999)

#ifdef WITH_SUNDIALS    
    CALL TSSundialsSetTolerance(ts%ts,absTol,relTol,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in TSSundialsSetTolerance.",err,error,*999)
    ENDIF
#endif    
    
    EXITS("Petsc_TSSundialsSetTolerance")
    RETURN
999 ERRORSEXITS("Petsc_TSSundialsSetTolerance",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_TSSundialsSetTolerance
  !
  !================================================================================================================================
  !

  !Finalise the PETSc Vec structure and destroy the KSP
  SUBROUTINE Petsc_VecFinalise(x,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The Vec to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecFinalise",err,error,*999)

    IF(.NOT. PetscObjectIsNull(x%vec)) THEN
      CALL Petsc_VecDestroy(x,err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecFinalise")
    RETURN
999 ERRORSEXITS("Petsc_VecFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecFinalise
    
  !
  !================================================================================================================================
  !

  !Initialise the PETSc Vec structure
  SUBROUTINE Petsc_VecInitialise(x,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The Vec to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecInitialise",err,error,*999)

    x%vec%v=-2
    
    EXITS("Petsc_VecInitialise")
    RETURN
999 ERRORSEXITS("Petsc_VecInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecInitialise
  
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecAssemblyBegin routine.
  SUBROUTINE Petsc_VecAssemblyBegin(x,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to begin the assembly of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecAssemblyBegin",err,error,*999)

    CALL VecAssemblyBegin(x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecAssemblyBegin.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecAssemblyBegin")
    RETURN
999 ERRORSEXITS("Petsc_VecAssemblyBegin",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecAssemblyBegin
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecAssemblyEnd routine.
  SUBROUTINE Petsc_VecAssemblyEnd(x,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to end the assembly of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecAssemblyEnd",err,error,*999)

    CALL VecAssemblyEnd(x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecAssemblyEnd.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecAssemblyEnd")
    RETURN
999 ERRORSEXITS("Petsc_VecAssemblyEnd",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecAssemblyEnd
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCopy routine.
  SUBROUTINE Petsc_VecCopy(x,y,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to copy from
    TYPE(PetscVecType), INTENT(INOUT) :: y !<The vector to copy to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecCopy",err,error,*999)

    CALL VecCopy(x%vec,y%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCopy.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecCopy")
    RETURN
999 ERRORSEXITS("Petsc_VecCopy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecCopy
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreate routine.
  SUBROUTINE Petsc_VecCreate(communicator,x,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator
#endif    
    TYPE(PetscVecType), INTENT(INOUT) :: x !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecCreate",err,error,*999)

    CALL VecCreate(communicator,x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecCreate")
    RETURN
999 ERRORSEXITS("Petsc_VecCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecCreate
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateGhost routine.
  SUBROUTINE Petsc_VecCreateGhost(communicator,localN,globalN,numGhosts,ghosts,x,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator
#endif    
    INTEGER(INTG), INTENT(IN) :: localN !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: globalN !<The number of global elements
    INTEGER(INTG), INTENT(IN) :: numGhosts !<The number of ghost elements
    INTEGER(INTG), INTENT(IN) :: ghosts(:) !<The global location of the each ghost element
    TYPE(PetscVecType), INTENT(INOUT) :: x !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecCreateGhost",err,error,*999)

    CALL VecCreateGhost(communicator,localN,globalN,numGhosts,ghosts,x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateGhost.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecCreateGhost")
    RETURN
999 ERRORSEXITS("Petsc_VecCreateGhost",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecCreateGhost
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateGhostWithArray routine.
  SUBROUTINE Petsc_VecCreateGhostWithArray(communicator,localN,globalN,numGhosts,ghosts,array,x,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator
#endif    
    INTEGER(INTG), INTENT(IN) :: localN !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: globalN !<The number of global elements
    INTEGER(INTG), INTENT(IN) :: numGhosts !<The number of ghost elements
    INTEGER(INTG), INTENT(IN) :: ghosts(:) !<The global location of the each ghost element
    REAL(DP), INTENT(OUT) :: array(:) !<The preallocated array of matrix data
    TYPE(PetscVecType), INTENT(INOUT) :: x !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecCreateGhostWithArray",err,error,*999)

    CALL VecCreateGhostWithArray(communicator,localN,globalN,numGhosts,ghosts,array,x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateGhostWithArray.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecCreateGhostWithArray")
    RETURN
999 ERRORSEXITS("Petsc_VecCreateGhostWithArray",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecCreateGhostWithArray
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateMPI routine.
  SUBROUTINE Petsc_VecCreateMPI(communicator,localN,globalN,x,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator
#endif    
    INTEGER(INTG), INTENT(IN) :: localN !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: globalN !<The number of global elements
    TYPE(PetscVecType), INTENT(INOUT) :: x !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecCreateMPI",err,error,*999)

    CALL VecCreateMPI(communicator,localN,globalN,x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateMPI.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecCreateMPI")
    RETURN
999 ERRORSEXITS("Petsc_VecCreateMPI",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecCreateMPI
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateMPIWithArray routine.
  SUBROUTINE Petsc_VecCreateMPIWithArray(communicator,bs,localN,globalN,array,x,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator
#endif    
    INTEGER(INTG), INTENT(IN) :: bs !<The block size
    INTEGER(INTG), INTENT(IN) :: localN !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: globalN !<The number of global elements
    REAL(DP), INTENT(OUT) :: array(:) !<The preallocated array for the vector data
    TYPE(PetscVecType), INTENT(INOUT) :: x !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecCreateMPIWithArray",err,error,*999)

    CALL VecCreateMPIWithArray(communicator,bs,localN,globalN,array,x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateMPIWithArray.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecCreateMPIWithArray")
    RETURN
999 ERRORSEXITS("Petsc_VecCreateMPIWithArray",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecCreateMPIWithArray
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateSeq routine.
  SUBROUTINE Petsc_VecCreateSeq(communicator,n,x,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator
#endif    
    INTEGER(INTG), INTENT(IN) :: n !<The size of the vector
    TYPE(PetscVecType), INTENT(INOUT) :: x !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecCreateSeq",err,error,*999)

    CALL VecCreateSeq(communicator,n,x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateSeq.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecCreateSeq")
    RETURN
999 ERRORSEXITS("Petsc_VecCreateSeq",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecCreateSeq
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecCreateSeqWithArray routine.
  SUBROUTINE Petsc_VecCreateSeqWithArray(communicator,bs,n,array,x,err,error,*)

    !Argument Variables
#ifdef WITH_F08_MPI
    TYPE(MPI_Comm), INTENT(IN) :: communicator !<The MPI communicator
#else    
    MPI_Comm, INTENT(IN) :: communicator !<The MPI communicator
#endif    
    INTEGER(INTG), INTENT(IN) :: bs !<The block size 
    INTEGER(INTG), INTENT(IN) :: n !<The size of the vector
    REAL(DP), INTENT(OUT) :: array(:) !<The preallocated array for the vector data
    TYPE(PetscVecType), INTENT(INOUT) :: x !<On exit, the created vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecCreateSeqWithArray",err,error,*999)

    CALL VecCreateSeqWithArray(communicator,bs,n,array,x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecCreateSeqWithArray.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecCreateSeqWithArray")
    RETURN
999 ERRORSEXITS("Petsc_VecCreateSeqWithArray",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecCreateSeqWithArray
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecDestroy routine.
  SUBROUTINE Petsc_VecDestroy(x,err,error,*)

   !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecDestroy",err,error,*999)

    CALL VecDestroy(x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecDestroy.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecDestroy")
    RETURN
999 ERRORSEXITS("Petsc_VecDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecDestroy
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecDuplicate routine.
  SUBROUTINE Petsc_VecDuplicate(x,y,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to duplicate
    TYPE(PetscVecType), INTENT(OUT) :: y !<On exit, the new duplicated vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecDuplicate",err,error,*999)

    CALL VecDuplicate(x%vec,y%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecDuplicate.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecDuplicate")
    RETURN
999 ERRORSEXITS("Petsc_VecDuplicate",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecDuplicate
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecDot routine.
  SUBROUTINE Petsc_VecDot(x,y,dotProduct,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(IN) :: x !<The vector x
    TYPE(PetscVecType), INTENT(IN) :: y !<The vector y
    REAL(DP), INTENT(OUT) :: dotProduct !<The dot product 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Petsc_VecDot",err,error,*999)
    
    CALL VecDot(x%vec,y%vec,dotProduct,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecDot.",err,error,*999)
    ENDIF

    EXITS("Petsc_VecDot")
    RETURN
999 ERRORSEXITS("Petsc_VecDot",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecDot
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetArrayF90 routine.
  SUBROUTINE Petsc_VecGetArrayF90(x,array,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT), TARGET :: x !<The vector to get the array of
    REAL(DP), POINTER :: array(:) !<On exit, a pointer to the array of the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecGetArrayF90",err,error,*999)

    IF(ASSOCIATED(array)) THEN
      CALL FlagError("Array is already associated",err,error,*999)
    ELSE
      CALL VecGetArrayF90(x%vec,array,err)
      IF(err/=0) THEN
        IF(petscHandleError) THEN
          CHKERRQ(err)
        ENDIF
        CALL FlagError("PETSc error in VecGetArrayF90.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Petsc_VecGetArrayF90")
    RETURN
999 ERRORSEXITS("Petsc_VecGetArrayF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecGetArrayF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetArrayReadF90 routine.
  SUBROUTINE Petsc_VecGetArrayReadF90(x,array,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT), TARGET :: x !<The vector to get the array of
    REAL(DP), POINTER :: array(:) !<On exit, a pointer to the array of the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecGetArrayReadF90",err,error,*999)

    IF(ASSOCIATED(array)) THEN
      CALL FlagError("Array is already associated",err,error,*999)
    ELSE
      CALL VecGetArrayReadF90(x%vec,array,err)
      IF(err/=0) THEN
        IF(petscHandleError) THEN
          CHKERRQ(err)
        ENDIF
        CALL FlagError("PETSc error in VecGetArrayReadF90.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Petsc_VecGetArrayReadF90")
    RETURN
999 ERRORSEXITS("Petsc_VecGetArrayReadF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecGetArrayReadF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetLocalSize routine.
  SUBROUTINE Petsc_VecGetLocalSize(x,n,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to get the local size of
    INTEGER(INTG), INTENT(OUT) :: n !<On exit, the local size of the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecGetLocalSize",err,error,*999)

    CALL VecGetLocalSize(x%vec,n,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGetLocalSize.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecGetLocalSize")
    RETURN
999 ERRORSEXITS("Petsc_VecGetLocalSize",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecGetLocalSize
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetOwnershipRange routine.
  SUBROUTINE Petsc_VecGetOwnershipRange(x,low,high,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to get the ownership range of 
    INTEGER(INTG), INTENT(OUT) :: low !<On exit, the low end of the range
    INTEGER(INTG), INTENT(OUT) :: high !<On exit, the high end of the range
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecGetOwnershipRange",err,error,*999)

    CALL VecGetOwnershipRange(x%vec,low,high,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
        ENDIF
      CALL FlagError("PETSc error in VecGetOwnershipRange.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecGetOwnershipRange")
    RETURN
999 ERRORSEXITS("Petsc_VecGetOwnershipRange",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecGetOwnershipRange

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetSize routine.
  SUBROUTINE Petsc_VecGetSize(x,n,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to get the size of
    INTEGER(INTG), INTENT(OUT) :: n !<On exit, the size of the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecGetSize",err,error,*999)

    CALL VecGetSize(x%vec,n,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGetSize.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecGetSize")
    RETURN
999 ERRORSEXITS("Petsc_VecGetSize",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecGetSize
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGetValues routine.
  SUBROUTINE Petsc_VecGetValues(x,n,indices,values,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to set the values for
    INTEGER(INTG), INTENT(IN) :: n !<The number of indicies to get
    INTEGER(INTG), INTENT(IN) :: indices(:) !<The indices to get
    REAL(DP), INTENT(OUT) :: values(:) !<On return, the values at the specified indicies
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecGetValues",err,error,*999)

    CALL VecGetValues(x%vec,n,indices,values,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGetValues.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecGetValues")
    RETURN
999 ERRORSEXITS("Petsc_VecGetValues",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecGetValues
    
!
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostGetLocalForm routine.
  SUBROUTINE Petsc_VecGhostGetLocalForm(g,l,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: g !<The global form of the vector
    TYPE(PetscVecType), INTENT(INOUT) :: l !<On exit, the local form of the vector with ghosts
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecGhostGetLocalForm",err,error,*999)

    CALL VecGhostGetLocalForm(g%vec,l%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGhostGetLocalForm.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecGhostGetLocalForm")
    RETURN
999 ERRORSEXITS("Petsc_VecGhostGetLocalForm",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecGhostGetLocalForm
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostRestoreLocalForm routine.
  SUBROUTINE Petsc_VecGhostRestoreLocalForm(g,l,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: g !<The global form of the vector
    TYPE(PetscVecType), INTENT(INOUT) :: l !<The local form of the vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecGhostRestoreLocalForm",err,error,*999)

    CALL VecGhostRestoreLocalForm(g%vec,l%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGhostRestoreLocalForm.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecGhostRestoreLocalForm")
    RETURN
999 ERRORSEXITS("Petsc_VecGhostRestoreLocalForm",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecGhostRestoreLocalForm
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostUpdateBegin routine.
  SUBROUTINE Petsc_VecGhostUpdateBegin(x,insertMode,scatterMode,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to begin the ghost update for
    InsertMode, INTENT(IN) :: insertMode !<The insert mode
    ScatterMode, INTENT(IN) :: scatterMode !<The scatter mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecGhostUpdateBegin",err,error,*999)

    CALL VecGhostUpdateBegin(x%vec,insertMode,scatterMode,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGhostUpdateBegin.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecGhostUpdateBegin")
    RETURN
999 ERRORSEXITS("Petsc_VecGhostUpdateBegin",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecGhostUpdateBegin
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecGhostUpdateEnd routine.
  SUBROUTINE Petsc_VecGhostUpdateEnd(x,insertMode,scatterMode,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to end the ghost update for
    InsertMode, INTENT(IN) :: insertMode !<The insert mode
    ScatterMode, INTENT(IN) :: scatterMode !<The scatter mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecGhostUpdateEnd",err,error,*999)

    CALL VecGhostUpdateEnd(x%vec,insertMode,scatterMode,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecGhostUpdateEnd.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecGhostUpdateEnd")
    RETURN
999 ERRORSEXITS("Petsc_VecGhostUpdateEnd",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecGhostUpdateEnd
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecNorm routine.
  SUBROUTINE Petsc_VecNorm(x,normType,norm,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(IN) :: x !<The vector x to find the norm of
    NormType, INTENT(IN) :: normType !<The norm type
    REAL(DP), INTENT(OUT) :: norm !<On exit, the vector norm
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Petsc_VecNorm",err,error,*999)
    
    CALL VecNorm(x%vec,normType,norm,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecNorm.",err,error,*999)
    ENDIF

    EXITS("Petsc_VecNorm")
    RETURN
999 ERRORSEXITS("Petsc_VecNorm",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecNorm

  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecRestoreArrayF90 routine.
  SUBROUTINE Petsc_VecRestoreArrayF90(x,array,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to restore the array of
    REAL(DP), POINTER :: array(:) !<A pointer to the data to restore
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecRestoreArrayF90",err,error,*999)

    CALL VecRestoreArrayF90(x%vec,array,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecRestoreArrayF90.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecRestoreArrayF90")
    RETURN
999 ERRORSEXITS("Petsc_VecRestoreArrayF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecRestoreArrayF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecRestoreArrayReadF90 routine.
  SUBROUTINE Petsc_VecRestoreArrayReadF90(x,array,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT), TARGET :: x !<The vector to restore the array for
    REAL(DP), POINTER :: array(:) !<A pointer to the array of the vector to restore
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecRestoreArrayReadF90",err,error,*999)

    CALL VecRestoreArrayReadF90(x%vec,array,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecRestoreArrayReadF90.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecRestoreArrayReadF90")
    RETURN
999 ERRORSEXITS("Petsc_VecRestoreArrayReadF90",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecRestoreArrayReadF90
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecScale routine.
  SUBROUTINE Petsc_VecScale(x,alpha,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to scale
    REAL(DP), INTENT(IN) :: alpha !<The scaling factor
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecScale",err,error,*999)

    CALL VecScale(x%vec,alpha,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecScale.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecScale")
    RETURN
999 ERRORSEXITS("Petsc_VecScale",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecScale
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSet routine.
  SUBROUTINE Petsc_VecSet(x,VALUE,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to set the value of
    REAL(DP), INTENT(IN) :: value !<The value to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecSet",err,error,*999)

    CALL VecSet(x%vec,value,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSet.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecSet")
    RETURN
999 ERRORSEXITS("Petsc_VecSet",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecSet
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetFromOptions routine.
  SUBROUTINE Petsc_VecSetFromOptions(x,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to set the options for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecSetFromOptions",err,error,*999)

    CALL VecSetFromOptions(x%vec,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSetFromOptions.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecSetFromOptions")
    RETURN
999 ERRORSEXITS("Petsc_VecSetFromOptions",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecSetFromOptions
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetLocalToGlobalMapping routine.
  SUBROUTINE Petsc_VecSetLocalToGlobalMapping(x,isLocalToGlobalMapping,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to set the local to global mapping for
    TYPE(PetscISLocalToGloabalMappingType), INTENT(IN) :: isLocalToGlobalMapping !<The local to global mapping context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecSetLocalToGlobalMapping",err,error,*999)

    CALL VecSetLocalToGlobalMapping(x%vec,isLocalToGlobalMapping%isLocalToGlobalMapping,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSetLocalToGlobalMapping.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecSetLocalToGlobalMapping")
    RETURN
999 ERRORSEXITS("Petsc_VecSetLocalToGlobalMapping",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecSetLocalToGlobalMapping
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetSizes routine.
  SUBROUTINE Petsc_VecSetSizes(x,localN,globalN,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to set the sizes of
    INTEGER(INTG), INTENT(IN) :: localN !<The number of local elements
    INTEGER(INTG), INTENT(IN) :: globalN !<The number of global elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecSetSizes",err,error,*999)

    CALL VecSetSizes(x%vec,localN,globalN,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSetSizes.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecSetSizes")
    RETURN
999 ERRORSEXITS("Petsc_VecSetSizes",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecSetSizes
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetValues routine.
  SUBROUTINE Petsc_VecSetValues(x,n,indices,values,insertMode,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to set the values for
    INTEGER(INTG), INTENT(IN) :: n !<The number of indicies
    INTEGER(INTG), INTENT(IN) :: indices(*) !<The indices
    REAL(DP), INTENT(IN) :: values(*) !<The values to set
    InsertMode, INTENT(IN) :: insertMode !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecSetValues",err,error,*999)

    CALL VecSetValues(x%vec,n,indices,values,insertMode,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSetValues.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecSetValues")
    RETURN
999 ERRORSEXITS("Petsc_VecSetValues",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecSetValues
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecSetValuesLocal routine.
  SUBROUTINE Petsc_SetValuesLocal(x,n,indices,values,insertMode,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to set the values of
    INTEGER(INTG), INTENT(IN) :: n !<The number of indices
    INTEGER(INTG), INTENT(IN) :: indices(*) !<The local indices
    REAL(DP), INTENT(IN) :: values(*) !<The values to set
    InsertMode :: insertMode !<The insert mode
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_SetValuesLocal",err,error,*999)

    CALL VecSetValuesLocal(x%vec,n,indices,values,insertMode,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecSetValuesLocal.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_SetValuesLocal")
    RETURN
999 ERRORSEXITS("Petsc_SetValuesLocal",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_SetValuesLocal
    
  !
  !================================================================================================================================
  !

  !>Buffer routine to the PETSc VecView routine.
  SUBROUTINE Petsc_VecView(x,viewer,err,error,*)

    !Argument Variables
    TYPE(PetscVecType), INTENT(INOUT) :: x !<The vector to view
    PetscViewer, INTENT(IN) :: viewer !<The viewer
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Petsc_VecView",err,error,*999)

    CALL VecView(x%vec,viewer,err)
    IF(err/=0) THEN
      IF(petscHandleError) THEN
        CHKERRQ(err)
      ENDIF
      CALL FlagError("PETSc error in VecView.",err,error,*999)
    ENDIF
    
    EXITS("Petsc_VecView")
    RETURN
999 ERRORSEXITS("Petsc_VecView",err,error)
    RETURN 1
    
  END SUBROUTINE Petsc_VecView

  !
  !================================================================================================================================
  !

END MODULE OpenCMISSPETSc
