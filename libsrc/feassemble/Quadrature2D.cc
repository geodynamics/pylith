// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "Quadrature2D.hh" // implementation of class methods

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature2D::Quadrature2D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature2D::~Quadrature2D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature2D::Quadrature2D(const Quadrature2D& q) :
  Quadrature(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell.
void
pylith::feassemble::Quadrature2D::_computeGeometry(
		       const ALE::Obj<ALE::Mesh::section_type>& coordinates,
		       const ALE::Mesh::point_type& cell)
{ // _computeGeometry
  assert(2 == _cellDim);
  assert(2 == _spaceDim);
  assert(0 != _basisDeriv);
  assert(0 != _quadPtsRef);
  assert(0 != _quadPts);
  assert(0 != _quadWts);
  assert(0 != _jacobian);
  assert(0 != _jacobianInv);

  _resetGeometry();

  // Get coordinates of cell's vertices
  const ALE::Mesh::topology_type::patch_type patch  = 0;
  const ALE::Mesh::section_type::value_type* vertCoords = 
    coordinates->restrict(patch, cell);
  //assert(2 == coordinates.GetFiberDimensionByDepth(patch,
  //*vertices->begin(), 0));

  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < _numQuadPts; ++iQuadPt) {
    
    // Compute coordinates of quadrature point in cell
    // x = sum[i=0,n-1] (Ni * xi)
    // y = sum[i=0,n-1] (Ni * yi)
    for (int iVertex=0, iB=iQuadPt*_numCorners;
	 iVertex < _numCorners;
	 ++iVertex) {
      const double basis = _basis[iB+iVertex];
      for (int iDim=0, iQ=iQuadPt*_spaceDim, iV=iVertex*_spaceDim;
	   iDim < _spaceDim;
	   ++iDim)
	_quadPts[iQ+iDim] +=  basis * vertCoords[iV+iDim];
    } // for
    
    // Compute Jacobian at quadrature point
    // J = [dx/dp dy/dp]
    //     [dx/dq dy/dq]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    // dx/dq = sum[i=0,n-1] (dNi/dq * xi)
    // dy/dq = sum[i=0,n-1] (dNi/dq * yi)
    for (int iVertex=0;
	 iVertex < _numCorners;
	 ++iVertex)
      for (int iRow=0, 
	     iB=iQuadPt*_numCorners*_spaceDim+iVertex*_cellDim;
	   iRow < _cellDim;
	   ++iRow) {
	const double deriv = _basisDeriv[iB+iRow];
	for (int iCol=0, iJ=iQuadPt*_cellDim*_spaceDim + iRow*_spaceDim;
	     iCol < _spaceDim; ++iCol)
	  _jacobian[iJ + iCol] += deriv * vertCoords[iVertex*_spaceDim+iCol];
      } // for
  
    // Compute determinant of Jacobian at quadrature point
    // |J| = j00*j11-j01*j10
    const int iJ = iQuadPt*_cellDim*_spaceDim;
    const int i00 = iJ + 0*_spaceDim + 0;
    const int i11 = iJ + 1*_spaceDim + 1;
    const int i01 = iJ + 0*_spaceDim + 1;
    const int i10 = iJ + 1*_spaceDim + 0;
    const double det = 
      _jacobian[i00]*_jacobian[i11] - _jacobian[i01]*_jacobian[i10];
    _checkJacobianDet(det);
    _jacobianDet[iQuadPt] = det;

    // Compute inverse of Jacobian at quadrature point
    // Jinv = 1/det*[ j11 -j01]
    //              [-j10  j00]
    _jacobianInv[i00] = _jacobian[i11] / det;
    _jacobianInv[i01] = -_jacobian[i01] / det;
    _jacobianInv[i10] = -_jacobian[i10] / det;
    _jacobianInv[i11] = _jacobian[i00] / det;
  } // for
} // _computeGeometry

// End of file 
