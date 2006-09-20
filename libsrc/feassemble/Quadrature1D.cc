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

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

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
// Compute geometric quantities for a cell.
void
pylith::feassemble::Quadrature1D::_computeGeometry(
		       const ALE::Obj<ALE::Mesh::section_type>& coordinates,
		       const ALE::Mesh::point_type& cell)
{ // _computeGeometry
  assert(1 == _cellDim);
  assert(1 == _spaceDim);
  assert(0 != _basisDeriv);
  assert(0 != _quadPtsRef);
  assert(0 != _quadPts);
  assert(0 != _quadWts);
  assert(0 != _jacobian);
  assert(0 != _jacobianInv);

  // Get coordinates of cell's vertices
  const ALE::Mesh::topology_type::patch_type patch  = 0;
  const ALE::Mesh::section_type::value_type* vertCoords = 
    coordinates->restrict(patch, cell);
  //assert(1 == coordinates.GetFiberDimension(patch, *vertices->begin()));

  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < _numQuadPts; ++iQuadPt) {

    // Compute coordinates of quadrature point in cell
    // x = sum[i=1,n] (Ni * xi)
    _quadPts[iQuadPt] = 0.0;
    for (int iVertex=0; iVertex < _numCorners; ++iVertex)
      _quadPts[iQuadPt] += 
	_basis[iQuadPt*_numCorners+iVertex]*vertCoords[iVertex];

    // Compute Jacobian at quadrature point
    // J = dx/dp = sum[i=1,n] (dNi/dp * xi)
    _jacobian[iQuadPt] = 0.0;
    for (int iVertex=0; iVertex < _numCorners; ++iVertex)
      _jacobian[iQuadPt] += 
	_basisDeriv[iQuadPt*_numCorners+iVertex] * vertCoords[iVertex];

    // Compute determinant of Jacobian at quadrature point
    // |J| = j00
    const double det = _jacobian[iQuadPt];
    if (det < _jacobianTol) {
      std::ostringstream msg;
      msg << "Determinant of Jacobian (" << det << ") is below minimum\n"
	  << "permissible value (" << _jacobianTol << ")!\n";
      throw std::runtime_error(msg.str());
    } // for
    _jacobianDet[iQuadPt] = _jacobian[iQuadPt];

    // Compute inverse of Jacobian at quadrature point
    // Jinv = 1/j00
    _jacobianInv[iQuadPt] = 1.0/_jacobianDet[iQuadPt];
  } // for
} // _computeGeometry

// End of file 
