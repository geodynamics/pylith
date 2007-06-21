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

#include "CohesiveKinData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::CohesiveKinData::CohesiveKinData(void) :
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
  finalSlipFilename(0),
  slipTimeFilename(0),
  peakRateFilename(0),
  matPropsFilename(0),
  fieldT(0),
  orientation(0),
  constraintVertices(0),
  constraintCells(0),
  valsResidual(0),
  valsJacobian(0),
  pseudoStiffness(0),
  numConstraintVert(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::faults::CohesiveKinData::~CohesiveKinData(void)
{ // destructor
} // destructor


// End of file
