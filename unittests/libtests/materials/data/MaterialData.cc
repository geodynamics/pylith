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
  numDBValues(0),
  numParameters(0),
  numParamValues(0),
  dbValues(0),
  parameterNames(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::materials::MaterialData::~MaterialData(void)
{ // destructor
} // destructor

// End of file
