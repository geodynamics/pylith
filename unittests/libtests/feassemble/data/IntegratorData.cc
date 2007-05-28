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
  meshFilename(0),
  spaceDim(0),
  cellDim(0),
  numBasis(0),
  numQuadPts(0),
  quadPts(0),
  quadWts(0),
  basis(0),
  basisDeriv(0),
  matType(0),
  matDBFilename(0),
  matId(0),
  matlabel(0),
  fieldTpdt(0),
  fieldT(0),
  fieldTmdt(0),
  valsResidual(0),
  valsJacobian(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorData::~IntegratorData(void)
{ // destructor
} // destructor


// End of file
