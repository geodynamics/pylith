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

#include "Quadrature1Din2D.hh" // implementation of class methods

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature1Din2D::Quadrature1Din2D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature1Din2D::~Quadrature1Din2D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature1Din2D::Quadrature1Din2D(const Quadrature1Din2D& q) :
  Quadrature(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell.
void
pylith::feassemble::Quadrature1Din2D::computeGeometry(
		       const ALE::Obj<Mesh>& mesh,
		       const ALE::Obj<real_section_type>& coordinates,
		       const Mesh::point_type& cell)
{ // computeGeometry
  assert(1 == _cellDim);
  assert(2 == _spaceDim);

  _resetGeometry();

  // Get coordinates of cell's vertices
  const real_section_type::value_type* vertCoords = 
    mesh->restrict(coordinates, cell);
  //assert(2 == coordinates.GetFiberDimensionByDepth(patch, 
  //*vertices->begin(), 0));

  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < _numQuadPts; ++iQuadPt) {
    
    // Compute coordinates of quadrature point in cell
    // x = sum[i=0,n-1] (Ni * xi)
    // y = sum[i=0,n-1] (Ni * yi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis) {
      const double basis = _basis[iQuadPt*_numBasis+iBasis];
      for (int iDim=0; iDim < _spaceDim; ++iDim)
	_quadPts[iQuadPt*_spaceDim+iDim] +=
	  basis * vertCoords[iBasis*_spaceDim+iDim];
    } // for
    
    // Compute Jacobian at quadrature point
    // J = [dx/dp dy/dp]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis) {
      const double deriv = _basisDeriv[iQuadPt*_numBasis+iBasis];
      for (int iDim=0; iDim < _spaceDim; ++iDim)
	_jacobian[iQuadPt*_spaceDim+iDim] += 
	  deriv * vertCoords[iBasis*_spaceDim+iDim];
    } // for

    // Compute determinant of Jacobian at quadrature point
    // |J| = sqrt(J transpose(J))
    double det = 0.0;
    for (int iDim=0; iDim < _spaceDim; ++iDim)
      det += _jacobian[iQuadPt*_spaceDim+iDim] * 
	_jacobian[iQuadPt*_spaceDim+iDim];
    det = sqrt(det);
    _checkJacobianDet(det);
    _jacobianDet[iQuadPt] = det;

    // Compute inverse of Jacobian at quadrature point
    // Jinv = 1.0/[J]
    for (int iDim=0; iDim < _spaceDim; ++iDim)
      _jacobianInv[iQuadPt*_spaceDim+iDim] = 
	1.0 / _jacobian[iQuadPt*_spaceDim+iDim];
  } // for
} // computeGeometry


// End of file 
