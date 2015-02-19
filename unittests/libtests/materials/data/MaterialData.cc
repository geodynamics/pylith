// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
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
  stateVars(0),
  propertiesNondim(0),
  stateVarsNondim(0),
  lengthScale(1.0e+3),
  timeScale(2.0),
  pressureScale(2.25e+10),
  densityScale(0)
{ // constructor
  const PylithScalar velScale = lengthScale / timeScale;
  densityScale = pressureScale / (velScale*velScale);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::materials::MaterialData::~MaterialData(void)
{ // destructor
} // destructor

// End of file
