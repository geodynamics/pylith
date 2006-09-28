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
  _basis(0),
  _basisDeriv(0),
  _quadPtsRef(0),
  _quadPts(0),
  _quadWts(0),
  _jacobian(0),
  _jacobianInv(0),
  _jacobianDet(0),
  _cellDim(0),
  _numCorners(0),
  _numQuadPts(0),
  _spaceDim(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature::~Quadrature(void)
{ // destructor
  delete[] _basis; _basis = 0;
  delete[] _basisDeriv; _basisDeriv = 0;
  delete[] _quadPtsRef; _quadPtsRef = 0;
  delete[] _quadPts; _quadPts = 0;
  delete[] _quadWts; _quadWts = 0;
  delete[] _jacobian; _jacobian = 0;
  delete[] _jacobianInv; _jacobianInv = 0;
  delete[] _jacobianDet; _jacobianDet = 0;
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor
pylith::feassemble::Quadrature::Quadrature(const Quadrature& q) :
  _minJacobian(q._minJacobian),
  _basis(0),
  _basisDeriv(0),
  _quadPtsRef(0),
  _quadPts(0),
  _quadWts(0),
  _jacobian(0),
  _jacobianInv(0),
  _jacobianDet(0),
  _cellDim(q._cellDim),
  _numCorners(q._numCorners),
  _numQuadPts(q._numQuadPts),
  _spaceDim(q._spaceDim)
{ // copy constructor
  if (0 != q._basis) {
    const int size = _numCorners * _numQuadPts;
    _basis = (size > 0) ? new double[size] : 0;
    memcpy(_basis, q._basis, size*sizeof(double));
  } // if

  if (0 != q._basisDeriv) {
    const int size = _numCorners * _numQuadPts * _cellDim;
    _basisDeriv = (size > 0) ? new double[size] : 0;
    memcpy(_basisDeriv, q._basisDeriv, size*sizeof(double));
  } // if

  if (0 != q._quadPtsRef) {
    const int size = _numQuadPts * _cellDim;
    _quadPtsRef = (size > 0) ? new double[size] : 0;
    memcpy(_quadPtsRef, q._quadPtsRef, size*sizeof(double));
  } // if

  if (0 != q._quadPts) {
    const int size = _numQuadPts*_spaceDim;
    _quadPts = (size > 0) ? new double[size] : 0;
    memcpy(_quadPts, q._quadPts, size*sizeof(double));
  } // if

  if (0 != q._quadWts) {
    const int size = _numQuadPts;
    _quadWts = (size > 0) ? new double[size] : 0;
    memcpy(_quadWts, q._quadWts, size*sizeof(double));
  } // if
  
  if (0 != q._jacobian) {
    const int size = _numQuadPts*_cellDim*_spaceDim;
    _jacobian = (size > 0) ? new double[size] : 0;
    memcpy(_jacobian, q._jacobian, size*sizeof(double));
  } // if

  if (0 != q._jacobianInv) {
    const int size = _numQuadPts*_cellDim*_spaceDim;
    _jacobianInv = (size > 0) ? new double[size] : 0;
    memcpy(_jacobianInv, q._jacobianInv, size*sizeof(double));
  } // if

  if (0 != q._jacobianDet) {
    const int size = _numQuadPts;
    _jacobianDet = (size > 0) ? new double[size] : 0;
    memcpy(_jacobianDet, q._jacobianDet, size*sizeof(double));
  } // if
} // copy constructor

// ----------------------------------------------------------------------
// Set basis functions and their derivatives and coordinates and
//   weights of the quadrature points.
void
pylith::feassemble::Quadrature::initialize(const double* basis,
					   const double* basisDeriv,
					   const double* quadPtsRef,
					   const double* quadWts,
					   const int cellDim,
					   const int numCorners,
					   const int numQuadPts,
					   const int spaceDim)
{ // initialize
  if (0 == basis ||
      0 == basisDeriv ||
      0 == quadPtsRef ||
      0 == quadWts ||
      cellDim < 1 || cellDim > 3 ||
      numCorners < 1 ||
      numQuadPts < 1 ||
      spaceDim < 1 || spaceDim > 3) {
    std::ostringstream msg;
    msg << "Incompatible values for quadrature information. Basis functions,\n"
	<< "their derivatives, and coordinates and weights of quadrature\n"
	<< "points must all be specified.\n"
	<< "Values:\n"
	<< "  basis pointer: " << basis << "\n"
	<< "  basis derivatites pointer: " << basisDeriv << "\n"
	<< "  quadrature points pointer: " << quadPtsRef << "\n"
	<< "  quadrature weights pointer: " << quadWts << "\n"
	<< "  space dimension: " << spaceDim << "\n"
	<< "  # vertices per cell: " << numCorners << "\n"
	<< "  # quadrature points: " << numQuadPts << "\n"
	<< "  dimension of coordinate space: " << spaceDim << "\n";
    throw std::runtime_error(msg.str());
  } // if

  int size = numCorners * numQuadPts;
  delete[] _basis; _basis = (size > 0) ? new double[size] : 0;
  for (int i=0; i < size; ++i)
    _basis[i] = basis[i];

  size = numCorners * numQuadPts * cellDim;
  delete[] _basisDeriv; _basisDeriv = (size > 0) ? new double[size] : 0;
  for (int i=0; i < size; ++i)
    _basisDeriv[i] = basisDeriv[i];

  size = numQuadPts * cellDim;
  delete[] _quadPtsRef; _quadPtsRef = (size > 0) ? new double[size] : 0;
  for (int i=0; i < size; ++i)
    _quadPtsRef[i] = quadPtsRef[i];

  size = numQuadPts;
  delete[] _quadWts; _quadWts = (size > 0) ? new double[size] : 0;
  for (int i=0; i < size; ++i)
    _quadWts[i] = quadWts[i];

  _cellDim = cellDim;
  _numCorners = numCorners;
  _numQuadPts = numQuadPts;
  _spaceDim = spaceDim;

  // Allocate for Jacobian and its inverse
  size = numQuadPts*cellDim*spaceDim;
  delete[] _jacobian; _jacobian = (size > 0) ? new double[size] : 0;
  delete[] _jacobianInv; _jacobianInv = (size > 0) ? new double[size] : 0;

  // Allocate for Jacobian determinant
  size = numQuadPts;
  delete[] _jacobianDet; _jacobianDet = (size > 0) ? new double[size] : 0;

  // Allocate for quad pts
  size = numQuadPts*spaceDim;
  delete[] _quadPts; _quadPts = (size > 0) ? new double[size] : 0;
} // initialize

// ----------------------------------------------------------------------
// Set entries in geometry arrays to zero.
void
pylith::feassemble::Quadrature::_resetGeometry(void)
{ // _resetGeometry
  // Zero out Jacobian and its inverse
  int size = _numQuadPts*_cellDim*_spaceDim;
  for (int i=0; i < size; ++i)
    _jacobian[i] = 0.0;
  for (int i=0; i < size; ++i)
    _jacobianInv[i] = 0.0;

  // Zero out Jacobian determinant
  size = _numQuadPts;
  for (int i=0; i < size; ++i)
    _jacobianDet[i] = 0.0;

  // Zero out quad pts
  size = _numQuadPts*_spaceDim;
  for (int i=0; i < size; ++i)
    _quadPts[i] = 0.0;
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
