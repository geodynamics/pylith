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

#include "QuadratureEngine.hh" // implementation of class methods

#include "CellGeometry.hh" // USES CellGeometry
#include "QuadratureBase.hh" // QuadratureBase

// ----------------------------------------------------------------------
// Constructor.
pylith::feassemble::QuadratureEngine::QuadratureEngine(const QuadratureBase& q) :
  _quadRefCell(q)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::QuadratureEngine::~QuadratureEngine(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Allocate cell buffers.
void
pylith::feassemble::QuadratureEngine::initialize(void)
{ // initialize
  const int numQuadPts = _quadRefCell.numQuadPts();
  const int numBasis = _quadRefCell.numBasis();
  const int cellDim = _quadRefCell.cellDim();
  const int spaceDim = _quadRefCell.spaceDim();

  _quadPts.resize(numQuadPts*spaceDim);
  _jacobian.resize(numQuadPts*cellDim*spaceDim);
  _jacobianInv.resize(numQuadPts*cellDim*spaceDim);
  _jacobianDet.resize(numQuadPts);
  _basisDeriv.resize(numQuadPts*numBasis*spaceDim);
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


// End of file 
