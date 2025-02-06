!> \file
!> \author Chris Bradley
!> \brief This module handles all VTK IO related routines.
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

!> This module handles all VTK IO related routines.
MODULE VTKRoutines

  USE BaseRoutines
  USE BasisAccessRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE DecompositionAccessRoutines
  USE ExportAccessRoutines
  USE FieldAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  
  
  IMPLICIT NONE

  
  PRIVATE

  !Module parameters


  !> \addtogroup VTK_Constants OpenCMISS::VTK::Constants
  !> This group contains all OpenCMISS VTK constants
  !>@{
  
  !> \addtogroup VTK_FileFormatTypes OpenCMISS::VTK::Constants::FileFormatTypes
  !> \brief VTK file format types
  !> \see VTKRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: VTK_FILE_FORMAT_LEGAGY = 1 !<Legacy VTK file format type \see VTK_FileFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_FILE_FORMAT_XML = 2 !<XML VTK file format type \see VTK_FileFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_FILE_FORMAT_VTKHDF = 3  !<HDF5 VTK file format type \see VTK_FileFormatTypes,VTKRoutines
  !>@}

  !> \addtogroup VTK_DataFormatTypes OpenCMISS::VTK::Constants::DataFormatTypes
  !> \brief VTK data format types
  !> \see VTKRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: VTK_ASCII_DATA_FORMAT = 1 !<ASCII data format type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BINARY_DATA_FORMAT = 2 !<Binary data format type \see VTK_DataFormatTypes,VTKRoutines
  !>@}

  !> \addtogroup VTK_DatasetTypes OpenCMISS::VTK::Constants::DatasetTypes
  !> \brief VTK data set types
  !> \see VTKRoutines
  !>@{t
  INTEGER(INTG), PARAMETER :: VTK_IMAGE_DATA = 1 !<Image dataset type \see VTK_DatasetTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_RECTLINEAR_DATA = 2 !<Rectlinear dataset type \see VTK_DatasetTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_STRUCTURED_GRID_DATA = 3 !<Structured grid dataset type \see VTK_DatasetTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_POLY_DATA = 4 !<Poly dataset type \see VTK_DatasetTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_UNSTRUCTURED_GRID_DATA = 5 !<Unstructured grid dataset type \see VTK_DatasetTypes,VTKRoutines
  !>@}
  
  !> \addtogroup VTK_DataArrayTypes OpenCMISS::VTK::Constants::DataArrayTypes
  !> \brief VTK data array types
  !> \see VTKRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: VTK_SCALAR_DATA_ARRAY = 1 !<Scalar data array type \see VTK_DataArrayTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_VECTOR_DATA_ARRAY = 2 !<Vector data array type \see VTK_DataArrayTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_NORMALS_DATA_ARRAY = 3 !<Normal data array type \see VTK_DataArrayTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_TENSORS_DATA_ARRAY = 4 !<Tensor data array type \see VTK_DataArrayTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_TCOORDS_DATA_ARRAY = 5 !<TCoord data array type \see VTK_DataArrayTypes,VTKRoutines
  !>@}
  
  !> \addtogroup VTK_DataArrayKind OpenCMISS::VTK::Constants::DataArrayKind
  !> \brief VTK data array kind (integer, float, double, etc.)
  !> \see VTKRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: VTK_INT8_DATA_KIND = 1 !<8-bit integer data kind \see VTK_DataArrayKind,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_UINT8_DATA_KIND = 2 !<Unsigned 8-bit integer data kind \see VTK_DataArrayKind,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_INT16_DATA_KIND = 3 !<16-bit integer data kind \see VTK_DataArrayKind,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_UINT16_DATA_KIND = 4 !<Unsigned 16-bit integer data kind \see VTK_DataArrayKind,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_INT32_DATA_KIND = 5 !<32-bit integer data kind \see VTK_DataArrayKind,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_UINT32_DATA_KIND = 6 !<Unsigned 32-bit integer data kind \see VTK_DataArrayKind,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_INT64_DATA_KIND = 7 !<64-bit integer data kind \see VTK_DataArrayKind,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_UINT64_DATA_KIND = 8 !<Unsigned 64-bit integer data kind \see VTK_DataArrayKind,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_FLOAT32_DATA_KIND = 9 !<32-bit real data kind \see VTK_DataArrayKind,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_FLOAT64_DATA_KIND = 10 !<64-bit real data kind \see VTK_DataArrayKind,VTKRoutines
  !>@}

  
  !>@}
  !> \addtogroup VTK_CellTypes OpenCMISS::VTK::Constants::CellTypes
  !> \brief VTK cell types
  !> \see VTKRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: VTK_VERTEX_CELL_TYPE = 1 !<Vertex cell type \see VTK_Types,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_POLY_VERTEX_CELL_TYPE = 2 !<Poly vertex cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_LINE_CELL_TYPE = 3 !<Line cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_POLY_LINE_CELL_TYPE = 4 !<Poly line cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_TRIANGLE_CELL_TYPE = 5 !<Triangle cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_TRIANGLE_STRIP_CELL_TYPE = 6 !<Triangle strip cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_POLYGON_CELL_TYPE = 7 !<Polygon cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_PIXEL_CELL_TYPE = 8 !<Pixel cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_QUAD_CELL_TYPE = 9 !<Quadrilateral cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_TETRA_CELL_TYPE = 10 !<Tetrahedral cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_VOXEL_CELL_TYPE = 11 !<Voxel cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_HEXAHEDRON_CELL_TYPE = 12 !<Hexahedron cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_WEDGE_CELL_TYPE = 13 !<Wedge cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_PYRAMID_CELL_TYPE = 14 !<Pyramid cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_PENTAGONAL_PRISM_CELL_TYPE = 15 !<Pentagonal pyramid cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_HEXAGONAL_PRISM_CELL_TYPE = 16 !<Hexagonal pyramid cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_QUADRATIC_EDGE_CELL_TYPE = 21 !<Quadratic edge cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_QUADRATIC_TRIANGLE_CELL_TYPE = 22 !<Quadratic triangle cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_QUADRATIC_QUAD_CELL_TYPE = 23 !<Quadratic quadrilateral cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_QUADRATIC_TETRA_CELL_TYPE = 24 !<Quadratic tetrahedron cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_QUADRATIC_HEXAHEDRON_CELL_TYPE = 25 !<Quadratic hexahedron cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_QUADRATIC_WEDGE_CELL_TYPE = 26 !<Quadratic wedge cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_QUADRATIC_PYRAMID_CELL_TYPE = 27 !<Quadratic pyramid cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BIQUADRATIC_QUAD_CELL_TYPE = 28 !<Bi-quadratic quadrilateral cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_TRIQUADRATIC_HEXAGON_CELL_TYPE = 29 !<Tri-quadratic hexagon cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_QUADRATIC_LINEAR_QUAD_CELL_TYPE = 30 !<Quadratic-linear quadrilateral cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_QUADRATIC_LINEAR_WEDGE_CELL_TYPE = 31 !<Quadratic-linear wedge cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BIQUADRATIC_QUADRATIC_WEDGE_CELL_TYPE = 32 !<Bi-quadratic-quadratic wedge cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON_CELL_TYPE = 33 !<Bi-quadratic-quadratic hexahedron cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BIQUADRATIC_TRIANGLE_CELL_TYPE = 34 !<Bi-quadratic triangle cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_CUBIC_LINE_CELL_TYPE = 35 !<Cubic line cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BEZIER_CURVE_CELL_TYPE = 75 !<Bezier curve cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BEZIER_TRIANGLE_CELL_TYPE = 76 !<Bezier triangle cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BEZIER_QUADRILATERAL_CELL_TYPE = 77 !<Bezier quadrilateral cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BEZIER_TETRAHEDRON_CELL_TYPE = 78 !<Bezier tetrahedron cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BEZIER_HEXAHEDRON_CELL_TYPE = 79 !<Bezier hexahedron cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BEZIER_WEDGE_CELL_TYPE = 80 !<Bezier wedge cell type \see VTK_DataFormatTypes,VTKRoutines
  INTEGER(INTG), PARAMETER :: VTK_BEZIER_PYRAMID_CELL_TYPE = 81 !<Bezier pyramid cell type \see VTK_DataFormatTypes,VTKRoutines
  !>@}

  !>@}

  INTEGER(INTG), PARAMETER :: VTKUNIT = 50
  
  !Module types

  

  !Module variables

  !Interfaces

  PUBLIC VTK_FILE_FORMAT_LEGAGY,VTK_FILE_FORMAT_XML,VTK_FILE_FORMAT_VTKHDF

  PUBLIC VTK_ASCII_DATA_FORMAT,VTK_BINARY_DATA_FORMAT

  PUBLIC VTK_Export

  PUBLIC VTK_Import

  PUBLIC VTKExport_CheckVariable

  PUBLIC VTKExport_CreateFinish

CONTAINS

  !
  !================================================================================================================================
  !

  !>Exports a VTK file
  SUBROUTINE VTK_Export(vtkExport,err,error,*)

    !Argument variables
    TYPE(VTKExportType), POINTER, INTENT(IN) :: vtkExport !<The VTK export information object
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: groupNodeNumber,numberOfGroupNodes
    LOGICAL :: isBinary,isParallel
    CHARACTER(LEN=6) :: dataFormat
    TYPE(DomainType), POINTER :: domain
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(ExportType), POINTER :: export
    TYPE(VARYING_STRING) :: filename,localError,vtkFilename
    TYPE(WorkGroupType), POINTER :: workGroup
        
    ENTERS("VTK_Export",err,error,*999)

    IF(.NOT.ASSOCIATED(vtkExport)) CALL FlagError("VTK export is not associated.",err,error,*999)

    NULLIFY(export)
    CALL VTKExport_ExportGet(vtkExport,export,err,error,*999)
   
    NULLIFY(domain)
    CALL VTKExport_DomainGet(vtkExport,domain,err,error,*999)
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(workgroup)
    CALL Decomposition_WorkGroupGet(decomposition,workgroup,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workgroup,numberOfGroupNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workgroup,groupNodeNumber,err,error,*999)
    isParallel = (numberOfGroupNodes > 1)
    isBinary = (vtkExport%vtkDataFormat == VTK_BINARY_DATA_FORMAT)

    SELECT CASE(vtkExport%vtkFormat)
    CASE(VTK_FILE_FORMAT_LEGAGY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(VTK_FILE_FORMAT_XML)
      IF(isParallel) THEN
        filename=export%baseFilename//"_"//TRIM(NumberToVString(groupNodeNumber,"I3.3",err,error))//".vtu"
      ELSE
        filename=export%baseFilename//".vtu"
      ENDIF
      IF(isBinary) THEN
        dataFormat = "binary"
      ELSE
        dataFormat = "ascii"
      ENDIF

      !Open the .vtu file
      OPEN(UNIT=VTKUNIT,FILE=CHAR(filename),ACCESS="stream",ACTION="write",ERR=998)

      !
      ! H E A D E R
      !
      !Write file header tag
      IF(IS_BIG_ENDIAN) THEN
        WRITE(VTKUNIT,'("<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""BigEndian"">")')
      ELSE
        WRITE(VTKUNIT,'("<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">")')
      ENDIF
      !
      ! S T A R T  T A G S
      !
      !Write start unstructured grid tag
      WRITE(VTKUNIT,'("<UnstructuredGrid>")')
      !Write start piece tag
      WRITE(VTKUNIT,'("<Piece> NumberOfPoints=""",I0,""" NumberOfCells=""",I0,""">")') vtkExport%numberOfPoints, &
        & vtkExport%numberOfCells
      !
      ! P O I N T S
      !
      !Write start points tag
      WRITE(VTKUNIT,'("<Points>")')
      !Write end points tag
      WRITE(VTKUNIT,'("</Points>")')
      !
      ! C E L L S
      !
      !Write start cells tag
      WRITE(VTKUNIT,'("<Cells>")')
      !   C e l l s  c o n n e c t i v i t y
      WRITE(VTKUNIT,'("<DataArray type=""Int32"" Name=""connectivity"" format=""",A,""">")')  &
        & dataFormat(1:LEN_TRIM(dataFormat))
      WRITE(VTKUNIT,'("</DataArray>")')
      !   C e l l s  o f f s e t s
      WRITE(VTKUNIT,'("<DataArray type=""Int32"" Name=""offsets"" format=""",A,""">")')  &
        & dataFormat(1:LEN_TRIM(dataFormat))
      WRITE(VTKUNIT,'("</DataArray>")')
      !   C e l l s  t y p e s
      WRITE(VTKUNIT,'("<DataArray type=""UInt8"" Name=""types"" format=""",A,""">")')  &
        & dataFormat(1:LEN_TRIM(dataFormat))
      WRITE(VTKUNIT,'("</DataArray>")')      
      !Write end cells tag
      WRITE(VTKUNIT,'("</Cells>")')
      !
      ! E N D  T A G S
      !
      !Write end piece tag
      WRITE(VTKUNIT,'("</Piece>")')
      !Write end unstructured grid tag      
      WRITE(VTKUNIT,'("</UnstructuredGrid>")')
      !End file tag
      WRITE(VTKUNIT,'("</VTKFile>")')
      !Close file
      CLOSE(VTKUNIT)
      
    CASE(VTK_FILE_FORMAT_VTKHDF)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The VTK export format of "//TRIM(NumberToVString(vtkExport%vtkFormat,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT      
    
    EXITS("VTK_Export")
    RETURN
998 CLOSE(VTKUNIT)
999 ERRORSEXITS("VTK_Export",err,error)
    RETURN 1

  END SUBROUTINE VTK_Export
  
  !
  !================================================================================================================================
  !

  !>Imports a VTK file
  SUBROUTINE VTK_Import(filename,err,error,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: filename !<The VTK file to import
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("VTK_Import",err,error,*999)

 
    EXITS("VTK_Import")
    RETURN
999 ERRORSEXITS("VTK_Import",err,error)
    RETURN 1

  END SUBROUTINE VTK_Import

  !
  !================================================================================================================================
  !

  !>Checks that a field variable can be handled by a VTK export
  SUBROUTINE VTKExport_CheckVariable(vtkExport,fieldVariable,startComponent,endComponent,err,error,*)

    !Argument variables
    TYPE(VTKExportType), POINTER, INTENT(IN) :: vtkExport !<The VTK export information object to check the variable for
    TYPE(FieldVariableType), POINTER, INTENT(IN) :: fieldVariable !<The field variable to check the export for
    INTEGER(INTG), INTENT(IN) :: startComponent !<The start component of the field variable to check
    INTEGER(INTG), INTENT(IN) :: endComponent !<The end component of the field variable to check
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementBasisType,elementIdx,interpolationOrder(4),interpolationType(1),numberOfElements,numberOfXi, &
      & numberOfXiCoordinates,xicIdx
    TYPE(BasisType), POINTER :: basis
    TYPE(DomainType), POINTER :: checkDomain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError
       
    ENTERS("VTKExport_CheckVariable",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vtkExport)) CALL FlagError("VTK export is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
#endif    

    !Check all components have the same domain
    CALL FieldVariable_AssertOneDomain(fieldVariable,startComponent,endComponent,err,error,*999)
    !Check the domain is the same as the current VTK export domain    
    IF(ASSOCIATED(vtkExport%domain)) THEN
      NULLIFY(checkDomain)
      CALL FieldVariable_DomainGet(fieldVariable,startComponent,checkDomain,err,error,*999)
      IF(.NOT.ASSOCIATED(vtkExport%domain,checkDomain))  &
        & CALL FlagError("The specified field variable has a different domain to the VTK export domain.",err,error,*999)
    ELSE
      CALL FieldVariable_DomainGet(fieldVariable,startComponent,vtkExport%domain,err,error,*999)
    ENDIF
    !Check geometric domain???
    !Check that the elements are compatible.
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(vtkExport%domain,domainTopology,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    !Get the number of local elements
    CALL DomainElements_NumberOfElementsGet(domainElements,numberOfElements,err,error,*999)
    !Loop over the local elements
    DO elementIdx=1,numberOfElements
      !Get the basis for the element
      NULLIFY(basis)
      CALL DomainElements_ElementBasisGet(domainElements,elementIdx,basis,err,error,*999)
      !Get the number of xi
      CALL Basis_NumberOfXiGet(basis,numberOfXi,err,error,*999)
      CALL Basis_NumberOfXiCoordinatesGet(basis,numberOfXiCoordinates,err,error,*999)
      !Get the basis type
      CALL Basis_TypeGet(basis,elementBasisType,err,error,*999)
      !Get the basis interpolation types and orders
      CALL Basis_InterpolationTypeGet(basis,interpolationType,err,error,*999)
      CALL Basis_InterpolationOrderGet(basis,interpolationOrder,err,error,*999)
      !For now just do elements that have the same interpolation in all directions
      DO xicIdx=2,numberOfXiCoordinates
        IF(interpolationType(xicIdx)/=interpolationType(1)) THEN
          localError="The interpolation type of "// &
            & TRIM(NumberToVString(interpolationType(xicIdx),"*",err,error))//" in xic direction "// &
            & TRIM(NumberToVString(xicIdx,"*",err,error))//" is different to the interpolation type of "//&
            & TRIM(NumberToVString(interpolationType(1),"*",err,error))//" in xic direction 1 for local element number "// &
            & TRIM(NumberToVString(elementIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(interpolationOrder(xicIdx)/=interpolationOrder(1)) THEN
          localError="The interpolation order of "// &
            & TRIM(NumberToVString(interpolationOrder(xicIdx),"*",err,error))//" in xic direction "// &
            & TRIM(NumberToVString(xicIdx,"*",err,error))//" is different to the interpolation order of "//&
            & TRIM(NumberToVString(interpolationOrder(1),"*",err,error))//" in xic direction 1 for local element number "// &
            & TRIM(NumberToVString(elementIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !xicIdx
      SELECT CASE(elementBasisType)
      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
        !Only Lagrange for now
        IF(interpolationType(1)/=BASIS_LAGRANGE_INTERPOLATION) THEN
          localError="The interpolation type of "//TRIM(NumberToVString(interpolationType(1),"*",err,error))// &
            & " is not implemented for VTK export."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        SELECT CASE(numberOfXi)
        CASE(1)
        CASE(2)
        CASE(3)
        CASE DEFAULT
          localError="The number of xi directions of "//TRIM(NumberToVString(numberOfXi,"*",err,error))// &
            & " for the basis for local element number "//TRIM(NumberToVString(elementIdx,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(BASIS_SIMPLEX_TYPE)
        !Linear and quadratic implemented
        SELECT CASE(interpolationOrder(1))
        CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
          !OK
        CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
          !OK
        CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
          localError="Local element number "//TRIM(NumberToVString(elementIdx,"*",err,error))// &
            & " has cubic Simplex interpolation which is not implemented for VTK export."
          CALL FlagError(localError,err,error,*999)
        CASE DEFAULT
          localError="The interpolation order of "//TRIM(NumberToVString(interpolationOrder(1),"*",err,error))// &
            & " for the basis in local element number "//TRIM(NumberToVString(elementIdx,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The basis type of "//TRIM(NumberToVString(elementBasisType,"*",err,error))// &
          & " is invalid or not implemented for VTK export."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !elementIdx
      
    EXITS("VTKExport_CheckVariable")
    RETURN
999 ERRORSEXITS("VTKExport_CheckVariable",err,error)
    RETURN 1

  END SUBROUTINE VTKExport_CheckVariable
  
  !
  !================================================================================================================================
  !

  !>Finishes the creation of a VTK export
  SUBROUTINE VTKExport_CreateFinish(vtkExport,err,error,*)

    !Argument variables
    TYPE(VTKExportType), POINTER :: vtkExport !<The VTK export to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("VTKExport_CreateFinish",err,error,*999)

 
    EXITS("VTKExport_CreateFinish")
    RETURN
999 ERRORSEXITS("VTKExport_CreateFinish",err,error)
    RETURN 1

  END SUBROUTINE VTKExport_CreateFinish

END MODULE VTKRoutines
