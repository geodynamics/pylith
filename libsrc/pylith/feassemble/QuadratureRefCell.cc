// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "QuadratureRefCell.hh" // implementation of class methods

#include "CellGeometry.hh" // USES CellGeometry

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::QuadratureRefCell::QuadratureRefCell(void) :
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
pylith::feassemble::QuadratureRefCell::~QuadratureRefCell(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::QuadratureRefCell::deallocate(void)
{ // deallocate
  delete _geometry; _geometry = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor
pylith::feassemble::QuadratureRefCell::QuadratureRefCell(const QuadratureRefCell& q) :
  _minJacobian(q._minJacobian),
  _quadPtsRef(q._quadPtsRef),
  _quadWts(q._quadWts),
  _basis(q._basis),
  _basisDerivRef(q._basisDerivRef),
  _cellDim(q._cellDim),
  _numBasis(q._numBasis),
  _numQuadPts(q._numQuadPts),
  _spaceDim(q._spaceDim),
  _geometry(0)
{ // copy constructor
  if (q._geometry)
    _geometry = q._geometry->clone();
} // copy constructor

// ----------------------------------------------------------------------
// Set basis functions and their derivatives and coordinates and
//   weights of the quadrature points.
void
pylith::feassemble::QuadratureRefCell::initialize(const PylithScalar* basis,
						  const int numQuadPts1,
						  const int numBasis1,
						  const PylithScalar* basisDerivRef,
						  const int numQuadPts2,
						  const int numBasis2,
						  const int cellDim2,
						  const PylithScalar* quadPtsRef,
						  const int numQuadPts3,
						  const int cellDim3,
						  const PylithScalar* quadWts,
						  const int numQuadPts4,
						  const int spaceDim)
{ // initialize
  PYLITH_METHOD_BEGIN;

  const int numQuadPts = numQuadPts1;
  const int numBasis = numBasis1;
  const int cellDim = (numBasis != 1) ? cellDim2 : 0;

  assert(numQuadPts == numQuadPts2);
  assert(numQuadPts == numQuadPts3);
  assert(numQuadPts == numQuadPts4);
  assert(numBasis == numBasis2);
  assert(cellDim2 == cellDim3);

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

  } // else

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Set geometry associated with reference cell.
void
pylith::feassemble::QuadratureRefCell::refGeometry(CellGeometry* const geometry)
{ // refGeometry
  PYLITH_METHOD_BEGIN;

  delete _geometry; _geometry = (geometry) ? geometry->clone() : 0;

  PYLITH_METHOD_END;
} // refGeometry

// ----------------------------------------------------------------------
// Get geometry associated with reference cell.
const pylith::feassemble::CellGeometry&
pylith::feassemble::QuadratureRefCell::refGeometry(void) const
{ // refGeometry
  assert(_geometry);
  return *_geometry;
} // refGeometry


// End of file 
