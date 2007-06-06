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

#include "Quadrature.hh" // implementation of class methods

#include <string.h> // USES memcpy()
#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature::Quadrature(void) :
  _minJacobian(0),
  _cellDim(0),
  _numBasis(0),
  _numQuadPts(0),
  _spaceDim(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature::~Quadrature(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor
pylith::feassemble::Quadrature::Quadrature(const Quadrature& q) :
  _minJacobian(q._minJacobian),
  _vertices(q._vertices),
  _quadPtsRef(q._quadPtsRef),
  _quadPts(q._quadPts),
  _quadWts(q._quadWts),
  _basis(q._basis),
  _basisDeriv(q._basisDeriv),
  _jacobian(q._jacobian),
  _jacobianInv(q._jacobianInv),
  _jacobianDet(q._jacobianDet),
  _cellDim(q._cellDim),
  _numBasis(q._numBasis),
  _numQuadPts(q._numQuadPts),
  _spaceDim(q._spaceDim)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Set basis functions and their derivatives and coordinates and
//   weights of the quadrature points.
void
pylith::feassemble::Quadrature::initialize(const double* vertices,
					   const double* basis,
					   const double* basisDeriv,
					   const double* quadPtsRef,
					   const double* quadWts,
					   const int cellDim,
					   const int numBasis,
					   const int numQuadPts,
					   const int spaceDim)
{ // initialize
  if (0 == vertices ||
      0 == basis ||
      0 == basisDeriv ||
      0 == quadPtsRef ||
      0 == quadWts ||
      cellDim < 1 || cellDim > 3 ||
      numBasis < 1 ||
      numQuadPts < 1 ||
      spaceDim < 1 || spaceDim > 3) {
    std::ostringstream msg;
    msg << "Incompatible values for quadrature information. Basis functions,\n"
	<< "their derivatives, and coordinates and weights of quadrature\n"
	<< "points must all be specified.\n"
	<< "Values:\n"
	<< "  vertices pointer: " << vertices << "\n"
	<< "  basis pointer: " << basis << "\n"
	<< "  basis derivatites pointer: " << basisDeriv << "\n"
	<< "  quadrature points pointer: " << quadPtsRef << "\n"
	<< "  quadrature weights pointer: " << quadWts << "\n"
	<< "  space dimension: " << spaceDim << "\n"
	<< "  # vertices per cell: " << numBasis << "\n"
	<< "  # quadrature points: " << numQuadPts << "\n"
	<< "  dimension of coordinate space: " << spaceDim << "\n";
    throw std::runtime_error(msg.str());
  } // if

  int size = numBasis * cellDim;
  _vertices.resize(size);
  for (int i=0; i < size; ++i)
    _vertices[i] = vertices[i];

  size = numBasis * numQuadPts;
  _basis.resize(size);
  for (int i=0; i < size; ++i)
    _basis[i] = basis[i];

  size = numBasis * numQuadPts * cellDim;
  _basisDeriv.resize(size);
  for (int i=0; i < size; ++i)
    _basisDeriv[i] = basisDeriv[i];

  size = numQuadPts * cellDim;
  _quadPtsRef.resize(size);
  for (int i=0; i < size; ++i)
    _quadPtsRef[i] = quadPtsRef[i];

  size = numQuadPts;
  _quadWts.resize(size);
  for (int i=0; i < size; ++i)
    _quadWts[i] = quadWts[i];

  _cellDim = cellDim;
  _numBasis = numBasis;
  _numQuadPts = numQuadPts;
  _spaceDim = spaceDim;

  // Allocate for Jacobian and its inverse
  size = numQuadPts * cellDim * spaceDim;
  _jacobian.resize(size);
  _jacobianInv.resize(size);

  // Allocate for Jacobian determinant
  size = numQuadPts;
  _jacobianDet.resize(size);

  // Allocate for quad pts
  size = numQuadPts*spaceDim;
  _quadPts.resize(size);
} // initialize

// ----------------------------------------------------------------------
// Set entries in geometry arrays to zero.
void
pylith::feassemble::Quadrature::_resetGeometry(void)
{ // _resetGeometry
  _jacobian = 0.0;
  _jacobianDet = 0.0;
  _jacobianInv = 0.0;
  _quadPts = 0.0;
} // _resetGeometry

// ----------------------------------------------------------------------
// Check determinant of Jacobian against minimum allowable value
void
pylith::feassemble::Quadrature::_checkJacobianDet(const double det) const
{ // _checkJacobianDet
  if (det < _minJacobian) {
    std::ostringstream msg;
    msg << "Determinant of Jacobian (" << det << ") is below minimum\n"
	<< "permissible value (" << _minJacobian << ")!\n";
    throw std::runtime_error(msg.str());
  } // if
} // _checkJacobianDet

// End of file 
