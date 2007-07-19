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

#include "CellGeometry.hh" // USES CellGeometry

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
  _spaceDim(0),
  _geometry(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature::~Quadrature(void)
{ // destructor
  delete _geometry; _geometry = 0;
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor
pylith::feassemble::Quadrature::Quadrature(const Quadrature& q) :
  _minJacobian(q._minJacobian),
  _quadPtsRef(q._quadPtsRef),
  _quadPts(q._quadPts),
  _quadWts(q._quadWts),
  _basis(q._basis),
  _basisDerivRef(q._basisDerivRef),
  _basisDeriv(q._basisDeriv),
  _jacobian(q._jacobian),
  _jacobianInv(q._jacobianInv),
  _jacobianDet(q._jacobianDet),
  _cellDim(q._cellDim),
  _numBasis(q._numBasis),
  _numQuadPts(q._numQuadPts),
  _spaceDim(q._spaceDim),
  _geometry(0)
{ // copy constructor
  if (0 != q._geometry)
    _geometry = q._geometry->clone();
} // copy constructor

// ----------------------------------------------------------------------
// Set basis functions and their derivatives and coordinates and
//   weights of the quadrature points.
void
pylith::feassemble::Quadrature::initialize(const double* basis,
					   const double* basisDerivRef,
					   const double* quadPtsRef,
					   const double* quadWts,
					   const int cellDim,
					   const int numBasis,
					   const int numQuadPts,
					   const int spaceDim)
{ // initialize
  if (0 == basis ||
      0 == basisDerivRef ||
      0 == quadPtsRef ||
      0 == quadWts ||
      cellDim < 0 || cellDim > 3 ||
      numBasis < 1 ||
      numQuadPts < 1 ||
      spaceDim < 1 || spaceDim > 3) {
    std::ostringstream msg;
    msg << "Incompatible values for quadrature information. Basis functions,\n"
	<< "their derivatives, and coordinates and weights of quadrature\n"
	<< "points must all be specified.\n"
	<< "Values:\n"
	<< "  basis pointer: " << basis << "\n"
	<< "  basis derivatites pointer: " << basisDerivRef << "\n"
	<< "  quadrature points pointer: " << quadPtsRef << "\n"
	<< "  quadrature weights pointer: " << quadWts << "\n"
	<< "  space dimension: " << spaceDim << "\n"
	<< "  # basis functions: " << numBasis << "\n"
	<< "  # quadrature points: " << numQuadPts << "\n"
	<< "  dimension of coordinate space: " << spaceDim << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (cellDim > 0) {
    int size = numBasis * numQuadPts; assert(size > 0);
    _basis.resize(size);
    for (int i=0; i < size; ++i)
      _basis[i] = basis[i];

    size = numQuadPts * numBasis * cellDim; assert(size > 0);
    _basisDerivRef.resize(size);
    for (int i=0; i < size; ++i)
      _basisDerivRef[i] = basisDerivRef[i];

    size = numQuadPts * cellDim; assert(size > 0);
    _quadPtsRef.resize(size);
    for (int i=0; i < size; ++i)
      _quadPtsRef[i] = quadPtsRef[i];

    size = numQuadPts; assert(size > 0);
    _quadWts.resize(size);
    for (int i=0; i < size; ++i)
      _quadWts[i] = quadWts[i];

    _cellDim = cellDim;
    _numBasis = numBasis;
    _numQuadPts = numQuadPts;
    _spaceDim = spaceDim;

    // Allocate for Jacobian and its inverse
    size = numQuadPts * cellDim * spaceDim; assert(size > 0);
    _jacobian.resize(size);
    _jacobianInv.resize(size);

    // Allocate for Jacobian determinant
    size = numQuadPts; assert(size > 0);
    _jacobianDet.resize(size);

    // Allocate for basis derivatives (in global coordinates)
    size = numQuadPts * numBasis * spaceDim; assert(size > 0);
    _basisDeriv.resize(size);

    // Allocate for quad pts
    size = numQuadPts*spaceDim; assert(size > 0);
    _quadPts.resize(size);
  } else {
    if (1 != numBasis ||
	1 != numQuadPts ||
	1 != spaceDim) {
      std::ostringstream msg;
      msg << "0-D quadrature only works in 1-D and is limited to 1 basis "
	  << "function and 1 quadrature point.\n"
	  << "Values:\n"
	  << "  cell dimension: " << cellDim << "\n"
	  << "  spatial dimension: " << spaceDim << "\n"
	  << "  # basis functions: " << numBasis << "\n"
	  << "  # quadrature points: " << numQuadPts << "\n";
      throw std::runtime_error(msg.str());
    } // if

    int size = 1;
    _basis.resize(size);
    for (int i=0; i < size; ++i)
      _basis[i] = basis[i];

    size = 1;
    _basisDerivRef.resize(size);
    for (int i=0; i < size; ++i)
      _basisDerivRef[i] = basisDerivRef[i];

    size = 1;
    _quadPtsRef.resize(size);
    for (int i=0; i < size; ++i)
      _quadPtsRef[i] = quadPtsRef[i];

    size = 1;
    _quadWts.resize(size);
    for (int i=0; i < size; ++i)
      _quadWts[i] = quadWts[i];

    _cellDim = cellDim;
    _numBasis = numBasis;
    _numQuadPts = numQuadPts;
    _spaceDim = spaceDim;

    // Allocate for Jacobian and its inverse
    size = 1;
    _jacobian.resize(size);
    _jacobianInv.resize(size);

    // Allocate for Jacobian determinant
    size = 1;
    _jacobianDet.resize(size);

    // Allocate for basis derivatives (in global coordinates)
    size = numQuadPts * numBasis * spaceDim; assert(size > 0);
    _basisDeriv.resize(size);

    // Allocate for quad pts
    size = spaceDim; assert(size > 0);
    _quadPts.resize(size);
  } // else
} // initialize

// ----------------------------------------------------------------------
// Set geometry associated with reference cell.
void
pylith::feassemble::Quadrature::refGeometry(CellGeometry* const geometry)
{ // refGeometry
  delete _geometry; _geometry = (0 != geometry) ? geometry->clone() : 0;
} // refGeometry

// ----------------------------------------------------------------------
// Get geometry associated with reference cell.
const pylith::feassemble::CellGeometry&
pylith::feassemble::Quadrature::refGeometry(void) const
{ // refGeometry
  assert(0 != _geometry);
  return *_geometry;
} // refGeometry

// ----------------------------------------------------------------------
// Set entries in geometry arrays to zero.
void
pylith::feassemble::Quadrature::_resetGeometry(void)
{ // _resetGeometry
  _jacobian = 0.0;
  _jacobianDet = 0.0;
  _jacobianInv = 0.0;
  _basisDeriv = 0.0;
  _quadPts = 0.0;
} // _resetGeometry

// ----------------------------------------------------------------------
// Check determinant of Jacobian against minimum allowable value
void
pylith::feassemble::Quadrature::_checkJacobianDet(const double det,
					   const Mesh::point_type& cell) const
{ // _checkJacobianDet
  if (det < _minJacobian) {
    std::ostringstream msg;
    msg << "Determinant of Jacobian (" << det << ") for cell " << cell
	<< " is smaller than minimum permissible value (" << _minJacobian
	<< ")!\n";
    throw std::runtime_error(msg.str());
  } // if
} // _checkJacobianDet


// End of file 
