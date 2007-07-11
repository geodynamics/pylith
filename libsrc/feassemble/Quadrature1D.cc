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

#include "Quadrature1D.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature1D::Quadrature1D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature1D::~Quadrature1D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature1D::Quadrature1D(const Quadrature1D& q) :
  Quadrature(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature1D::computeGeometry(
			      const ALE::Obj<Mesh>& mesh,
			      const ALE::Obj<real_section_type>& coordinates,
			      const Mesh::point_type& cell)
{ // computeGeometry
  assert(1 == _cellDim);
  assert(1 == _spaceDim);

  _resetGeometry();
  
  // Get coordinates of cell's vertices
  const real_section_type::value_type* vertCoords = 
    mesh->restrict(coordinates, cell);
  assert(1 == coordinates->getFiberDimension(*mesh->depthStratum(0)->begin()));
  
  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < _numQuadPts; ++iQuadPt) {
    
    // Compute coordinates of quadrature point in cell
    // x = sum[i=0,n-1] (Ni * xi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis)
      _quadPts[iQuadPt] += 
	_basis[iQuadPt*_numBasis+iBasis]*vertCoords[iBasis];

    // Compute Jacobian at quadrature point
    // J = dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis)
      _jacobian[iQuadPt] += 
	_basisDerivRef[iQuadPt*_numBasis+iBasis] * vertCoords[iBasis];

    // Compute determinant of Jacobian at quadrature point
    // |J| = j00
    const double det = _jacobian[iQuadPt];
    _checkJacobianDet(det, cell);
    _jacobianDet[iQuadPt] = _jacobian[iQuadPt];
    
    // Compute inverse of Jacobian at quadrature point
    // Jinv = 1/j00
    _jacobianInv[iQuadPt] = 1.0/_jacobianDet[iQuadPt];

    assert(_numQuadPts*_numBasis*_spaceDim == _basisDeriv.size());
    assert(_numQuadPts*_numBasis*_cellDim == _basisDerivRef.size());
    assert(_numQuadPts*_cellDim*_spaceDim == _jacobianInv.size());

    // Compute derivatives of basis functions with respect to global
    // coordinates
    // dNi/dx = dNi/dp dp/dx + dNi/dq dq/dx + dNi/dr dr/dx
    for (int iBasis=0; iBasis < _numBasis; ++iBasis)
      _basisDeriv[iQuadPt*_numBasis+iBasis] +=
	  _basisDerivRef[iQuadPt*_numBasis+iBasis] *
	  _jacobianInv[iQuadPt];
  } // for
} // computeGeometry


// End of file 
