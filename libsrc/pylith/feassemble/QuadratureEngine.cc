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

#include "QuadratureEngine.hh" // implementation of class methods

#include "CellGeometry.hh" // USES CellGeometry
#include "QuadratureRefCell.hh" // QuadratureRefCell

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor.
pylith::feassemble::QuadratureEngine::QuadratureEngine(const QuadratureRefCell& q) :
  _quadRefCell(q)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::QuadratureEngine::~QuadratureEngine(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::QuadratureEngine::deallocate(void)
{ // deallocate
} // deallocate
  
// ----------------------------------------------------------------------
// Allocate cell buffers.
void
pylith::feassemble::QuadratureEngine::initialize(void)
{ // initialize
  PYLITH_METHOD_BEGIN;

  const int numQuadPts = _quadRefCell.numQuadPts();
  const int numBasis = _quadRefCell.numBasis();
  const int cellDim = _quadRefCell.cellDim();
  const int spaceDim = _quadRefCell.spaceDim();

  if (cellDim > 0) {
    _quadPts.resize(numQuadPts*spaceDim);
    _jacobian.resize(numQuadPts*cellDim*spaceDim);
    _jacobianInv.resize(numQuadPts*cellDim*spaceDim);
    _jacobianDet.resize(numQuadPts);
    _basisDeriv.resize(numQuadPts*numBasis*spaceDim);
  } else {
    _quadPts.resize(numQuadPts*spaceDim);
    _jacobian.resize(1);
    _jacobianInv.resize(1);
    _jacobianDet.resize(1);
    _basisDeriv.resize(numQuadPts*numBasis*spaceDim);
  } // if/else

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Fill cell buffers with zeros.
void
pylith::feassemble::QuadratureEngine::zero(void)
{ // zero
  _quadPts = 0.0;
  _jacobian = 0.0;
  _jacobianDet = 0.0;
  _jacobianInv = 0.0;
  _basisDeriv = 0.0;
} // zero

// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::QuadratureEngine::QuadratureEngine(const QuadratureEngine& q) :
  _quadPts(q._quadPts),
  _jacobian(q._jacobian),
  _jacobianDet(q._jacobianDet),
  _jacobianInv(q._jacobianInv),
  _basisDeriv(q._basisDeriv),
  _quadRefCell(q._quadRefCell)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Check determinant of Jacobian against minimum allowable value
void
pylith::feassemble::QuadratureEngine::_checkJacobianDet(const PylithScalar det,
							const int cell) const
{ // _checkJacobianDet
  const PylithScalar minJacobian = _quadRefCell.minJacobian();
  if (det < minJacobian) {
    std::ostringstream msg;
    msg << "Determinant of Jacobian (" << det << ") for cell " << cell
	<< " is smaller than minimum permissible value (" << minJacobian
	<< ")!\n"
	<< "The two most likely causes of this are highly distorted cells "
	<< "and nondimensionalization with a length scale that is much larger "
	<< "than the dimensions of the cells.\n";
    throw std::runtime_error(msg.str());
  } // if
} // _checkJacobianDet


// End of file 
