// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
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
  frictionDeriv(0),
  slip(0),
  slipRate(0),
  normalTraction(0),
  stateVarsUpdated(0),
  setLengthScale(0),
  setTimeScale(0),
  setPressureScale(0),
  setDensityScale(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::friction::FrictionModelData::~FrictionModelData(void)
{ // destructor
} // destructor

// End of file
