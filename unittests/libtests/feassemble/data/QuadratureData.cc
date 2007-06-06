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

#include "QuadratureData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::QuadratureData::QuadratureData(void) :
  numVertices(0),
  spaceDim(0),
  numCells(0),
  cellDim(0),
  numBasis(0),
  numQuadPts(0),
  vertices(0),
  cells(0),
  verticesRef(0),
  quadPtsRef(0),
  quadWts(0),
  quadPts(0),
  basis(0),
  basisDeriv(0),
  jacobian(0),
  jacobianDet(0),
  jacobianInv(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::QuadratureData::~QuadratureData(void)
{ // destructor
} // destructor

// End of file
