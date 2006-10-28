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

#include "IntegratorData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorData::IntegratorData(void) :
  numVertices(0),
  spaceDim(0),
  numCells(0),
  cellDim(0),
  numCorners(0),
  numQuadPts(0),
  fiberDim(0),
  vertices(0),
  cells(0),
  quadPts(0),
  quadWts(0),
  basis(0),
  basisDeriv(0),
  fieldIn(0),
  valsAction(0),
  valsMatrix(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorData::~IntegratorData(void)
{ // destructor
} // destructor

// End of file
