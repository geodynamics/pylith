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

#include "MaterialData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::materials::MaterialData::MaterialData(void) :
  dimension(0),
  numLocs(0),
  numProperties(0),
  numStateVars(0),
  numDBProperties(0),
  numDBStateVars(0),
  numPropsQuadPt(0),
  numVarsQuadPt(0),
  numPropertyValues(0),
  numStateVarValues(0),
  dbPropertyValues(0),
  dbStateVarValues(0),
  dbProperties(0),
  dbStateVars(0),
  properties(0),
  stateVars(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::materials::MaterialData::~MaterialData(void)
{ // destructor
} // destructor

// End of file
