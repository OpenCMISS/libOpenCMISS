!> \file
!> \author Chris Bradley
!> \brief This module contains all kind definitions.
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

!> This module contains all kind definitions.
MODULE Kinds

  IMPLICIT NONE

  !Module parameters

  !> \addtogroup Kinds_IntegerKinds Kinds::IntegerKinds
  !> \brief Kind parameters for integer data types
  !> \see Kinds,OpenCMISS_IntegerKinds
  !>@{
  INTEGER, PARAMETER :: INTG=SELECTED_INT_KIND(9) !<Standard integer kind \see Kinds_IntegerKinds,Kinds
  INTEGER, PARAMETER :: SINTG=SELECTED_INT_KIND(4) !<Short integer kind \see Kinds_IntegerKinds,Kinds
  INTEGER, PARAMETER :: LINTG=SELECTED_INT_KIND(18) !<Long integer kind \see Kinds_IntegerKinds,Kinds
  INTEGER, PARAMETER :: PTR=INTG !<Pointer integer kind
  INTEGER, PARAMETER :: IDX=INTG !<Integer index kind
  INTEGER, PARAMETER :: LIDX=LINTG !<Long integer index kind
  !>@}
  
  !> \addtogroup Kinds_RealKinds Kinds::RealKinds
  !> \brief Kind parameters for real data types
  !> \see Kinds,OpenCMISS_RealKinds
  !>@{
  INTEGER, PARAMETER :: SP=SELECTED_REAL_KIND(6,15) !<Single precision real kind \see Kinds_RealKinds,Kinds
  INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(15,307) !<Double precision real kind \see Kinds_RealKinds,Kinds
  INTEGER, PARAMETER :: QP=SELECTED_REAL_KIND(30,1000) !<Quadruple precision real kind \see Kinds_RealKinds,Kinds
#ifdef SINGLE_REAL_PRECISION
  INTEGER, PARAMETER :: RP=SP !<Real working precision kind i.e., single, double, etc. \see Kinds_RealKinds,Kinds
#else
  INTEGER, PARAMETER :: RP=DP !<Real working precision kind i.e., single, double, etc. \see Kinds_RealKinds,Kinds
#endif
  !>@}

  !> \addtogroup Kinds_ComplexKinds Kinds::ComplexKinds
  !> \brief Kind parameters for complex data types
  !> \see Kinds,OpenCMISS_ComplexKinds
  !>@{
  INTEGER, PARAMETER :: SPC=KIND((1.0_SP,1.0_SP)) !<Single precision complex kind \see Kinds_ComplexKinds,Kinds
  INTEGER, PARAMETER :: DPC=KIND((1.0_DP,1.0_DP)) !<Double precision complex kind \see Kinds_ComplexKinds,Kinds
  !INTEGER, PARAMETER :: QPC=KIND((1.0_QP,1.0_QP))
#ifdef SINGLE_REAL_PRECISION
  INTEGER, PARAMETER :: RPC=SPC !<Real working precision complex kind i.e., single, double, etc. \see Kinds_ComplexKinds,Kinds
#else
  INTEGER, PARAMETER :: RPC=DPC !<Real working precision complex kind i.e., single, double, etc. \see Kinds_ComplexKinds,Kinds
#endif
  !>@}

END MODULE Kinds
