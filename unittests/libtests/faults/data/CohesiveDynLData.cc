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

#include "CohesiveDynLData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::CohesiveDynLData::CohesiveDynLData(void) :
  meshFilename(0),
  spaceDim(0),
  cellDim(0),
  numBasis(0),
  numQuadPts(0),
  quadPts(0),
  quadWts(0),
  basis(0),
  basisDeriv(0),
  verticesRef(0),
  id(0),
  label(0),
  initialTractFilename(0),
  fieldT(0),
  fieldIncr(0),
  jacobian(0),
  orientation(0),
  area(0),
  fieldIncrSlip(0),
  slipSlip(0),
  fieldIncrOpen(0),
  slipOpen(0),
  constraintVertices(0),
  constraintCells(0),
  numConstraintVert(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::faults::CohesiveDynLData::~CohesiveDynLData(void)
{ // destructor
} // destructor


// End of file
