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

#include "Quadrature0D.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature0D::Quadrature0D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature0D::~Quadrature0D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature0D::Quadrature0D(const Quadrature0D& q) :
  Quadrature(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature0D::computeGeometry(
			      const ALE::Obj<Mesh>& mesh,
			      const ALE::Obj<real_section_type>& coordinates,
			      const Mesh::point_type& cell)
{ // computeGeometry
  assert(0 == _cellDim);
  assert(1 == _numQuadPts);
  assert(1 == _numBasis);

  _resetGeometry();
  
  // Get coordinates of cell's vertices
  const real_section_type::value_type* vertCoords = 
    mesh->restrict(coordinates, cell);
  assert(1 == coordinates->getFiberDimension(*mesh->depthStratum(0)->begin()));

  for (int i=0; i < _spaceDim; ++i)
    _quadPts[i] = vertCoords[i];

  _jacobian[0] = 1.0;
  _jacobianDet[0] = 1.0;
  _jacobianInv[0] = 1.0;
} // computeGeometry


// End of file 
