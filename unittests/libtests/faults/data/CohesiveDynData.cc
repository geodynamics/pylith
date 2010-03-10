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

#include "CohesiveDynData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::CohesiveDynData::CohesiveDynData(void) :
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
  fieldIncrStick(0),
  fieldIncrSlip(0),
  fieldIncrOpen(0),
  jacobian(0),
  orientation(0),
  forcesInitial(0),
  area(0),
  fieldIncrSlipE(0),
  slipSlipE(0),
  fieldIncrOpenE(0),
  slipOpenE(0),
  constraintVertices(0),
  numConstraintVert(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::faults::CohesiveDynData::~CohesiveDynData(void)
{ // destructor
} // destructor


// End of file
