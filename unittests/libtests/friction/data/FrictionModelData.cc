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

#include "FrictionModelData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::friction::FrictionModelData::FrictionModelData(void) :
  numLocs(0),
  numProperties(0),
  numStateVars(0),
  numDBProperties(0),
  numDBStateVars(0),
  numPropsVertex(0),
  numVarsVertex(0),
  numPropertyValues(0),
  numStateVarValues(0),
  dbPropertyValues(0),
  dbStateVarValues(0),
  dt(0.0),
  dbProperties(0),
  dbStateVars(0),
  properties(0),
  stateVars(0),
  propertiesNondim(0),
  stateVarsNondim(0),
  friction(0),
  slip(0),
  slipRate(0),
  normalTraction(0),
  stateVarsUpdated(0),
  lengthScale(0),
  timeScale(0),
  pressureScale(0),
  densityScale(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::friction::FrictionModelData::~FrictionModelData(void)
{ // destructor
} // destructor

// End of file
