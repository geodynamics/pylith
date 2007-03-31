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

#include "Quadrature1Din3D.hh" // implementation of class methods

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature1Din3D::Quadrature1Din3D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature1Din3D::~Quadrature1Din3D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature1Din3D::Quadrature1Din3D(const Quadrature1Din3D& q) :
  Quadrature(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell.
void
pylith::feassemble::Quadrature1Din3D::computeGeometry(
		       const ALE::Obj<Mesh>& mesh,
		       const ALE::Obj<real_section_type>& coordinates,
		       const Mesh::point_type& cell)
{ // computeGeometry
  assert(1 == _cellDim);
  assert(3 == _spaceDim);
  assert(0 != _basisDeriv);
  assert(0 != _quadPtsRef);
  assert(0 != _quadPts);
  assert(0 != _quadWts);
  assert(0 != _jacobian);
  assert(0 != _jacobianInv);

  _resetGeometry();

  // Get coordinates of cell's vertices
  const real_section_type::value_type* vertCoords = 
    mesh->restrict(coordinates, cell);
  //assert(3 == coordinates.GetFiberDimensionByDepth(patch,
  //*vertices->begin(), 0));

  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < _numQuadPts; ++iQuadPt) {
    
    // Compute coordinates of quadrature point in cell
    // x = sum[i=0,n-1] (Ni * xi)
    // y = sum[i=0,n-1] (Ni * yi)
    // z = sum[i=0,n-1] (Ni * zi)
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
    // J = [dx/dp dy/dp dz/dp]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    // dz/dp = sum[i=0,n-1] (dNi/dp * zi)
    for (int iVertex=0, iB=iQuadPt*_numCorners;
	 iVertex < _numCorners;
	 ++iVertex) {
      const double deriv = _basisDeriv[iB+iVertex];
      for (int iDim=0, iJ=iQuadPt*_spaceDim, iV=iVertex*_spaceDim;
	   iDim < _spaceDim;
	   ++iDim)
	_jacobian[iJ+iDim] += deriv * vertCoords[iV+iDim];
    } // for

    // Compute determinant of Jacobian at quadrature point
    // |J| = sqrt(J transpose(J))
    double det = 0.0;
    for (int iDim=0, iJ=iQuadPt*_spaceDim; iDim < _spaceDim; ++iDim)
      det += _jacobian[iJ+iDim]*_jacobian[iJ+iDim];
    det = sqrt(det);
    _checkJacobianDet(det);
    _jacobianDet[iQuadPt] = det;

    // Compute inverse of Jacobian at quadrature point
    // Jinv = 1.0/[J]
    for (int iDim=0, iJ=iQuadPt*_spaceDim; iDim < _spaceDim; ++iDim)
      _jacobianInv[iJ+iDim] = 1.0/_jacobian[iJ+iDim];
  } // for
} // computeGeometry

// End of file 
