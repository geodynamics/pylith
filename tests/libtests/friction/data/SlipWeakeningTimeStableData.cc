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

#include "SlipWeakeningTimeStableData.hh"

const int pylith::friction::SlipWeakeningTimeStableData::_numLocs = 6;

const int pylith::friction::SlipWeakeningTimeStableData::_numProperties = 6;

const int pylith::friction::SlipWeakeningTimeStableData::_numStateVars = 2;

const int pylith::friction::SlipWeakeningTimeStableData::_numDBProperties = 6;

const int pylith::friction::SlipWeakeningTimeStableData::_numDBStateVars = 2;

const int pylith::friction::SlipWeakeningTimeStableData::_numPropsVertex = 6;

const int pylith::friction::SlipWeakeningTimeStableData::_numVarsVertex = 2;

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_lengthScale =   1.00000000e+03;

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_timeScale =   1.00000000e+01;

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_pressureScale =   2.25000000e+10;

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_densityScale =   1.00000000e+03;

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_dt = 0.01;

const int pylith::friction::SlipWeakeningTimeStableData::_numPropertyValues[6] = {
  1,
  1,
  1,
  1,
  1,
  1,
};

const int pylith::friction::SlipWeakeningTimeStableData::_numStateVarValues[2] = {
  1,
  1,
};

const char* pylith::friction::SlipWeakeningTimeStableData::_dbPropertyValues[6] = {
  "static-coefficient",
  "dynamic-coefficient",
  "slip-weakening-parameter",
  "cohesion",
  "time-weakening-time",
  "time-weakening-parameter",
};

const char* pylith::friction::SlipWeakeningTimeStableData::_dbStateVarValues[2] = {
  "cumulative-slip",
  "previous-slip",
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_dbProperties[6*6] = {
  0.6,
  0.5,
  0.8,
  1000000,
  2.0,
  0.5,

  0.6,
  0.5,
  0.8,
  1000000,
  1.25,
  0.5,

  0.6,
  0.5,
  0.8,
  1000000,
  0.6,
  0.5,

  0.6,
  0.5,
  0.4,
  1000000,
  2.0,
  0.5,

  0.6,
  0.5,
  0.4,
  1000000,
  1.25,
  0.5,

  0.6,
  0.5,
  0.4,
  1000000,
  0.6,
  0.5,
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_dbStateVars[6*2] = {
  0.15,
  0.1,

  0.15,
  0.1,

  0.15,
  0.1,

  0.8,
  0.4,

  0.8,
  0.4,

  0.8,
  0.4,
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_properties[6*6] = {
  0.6,
  0.5,
  0.8,
  1000000,
  2.0,
  0.5,

  0.6,
  0.5,
  0.8,
  1000000,
  1.25,
  0.5,

  0.6,
  0.5,
  0.8,
  1000000,
  0.6,
  0.5,

  0.6,
  0.5,
  0.4,
  1000000,
  2.0,
  0.5,

  0.6,
  0.5,
  0.4,
  1000000,
  1.25,
  0.5,

  0.6,
  0.5,
  0.4,
  1000000,
  0.6,
  0.5,
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_stateVars[6*2] = {
  0.15,
  0.1,

  0.15,
  0.1,

  0.15,
  0.1,

  0.8,
  0.4,

  0.8,
  0.4,

  0.8,
  0.4,
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_propertiesNondim[6*6] = {
  0.6,
  0.5,
  0.0008,
  0.000044444444,
  0.2,
  0.05,

  0.6,
  0.5,
  0.0008,
  0.000044444444,
  0.125,
  0.05,

  0.6,
  0.5,
  0.0008,
  0.000044444444,
  0.06,
  0.05,

  0.6,
  0.5,
  0.0004,
  0.000044444444,
  0.2,
  0.05,

  0.6,
  0.5,
  0.0004,
  0.000044444444,
  0.125,
  0.05,

  0.6,
  0.5,
  0.0004,
  0.000044444444,
  0.06,
  0.05,
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_stateVarsNondim[6*2] = {
  0.00015,
  0.0001,

  0.00015,
  0.0001,

  0.00015,
  0.0001,

  0.0008,
  0.0004,

  0.0008,
  0.0004,

  0.0008,
  0.0004,
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_friction[6] = {
  11.265e+5,
  11.210e+5,
  11.100e+5,
  11.15e+5,
  11.15e+5,
  11.15e+5,
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_frictionDeriv[6] = {
  -2.2e+5*0.1/0.8,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_slip[6] = {
  0.15,
  0.15,
  0.15,
  0.6,
  0.6,
  0.6,
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_slipRate[6] = {
  0.74,
  0.74,
  0.74,
  0.64,
  0.64,
  0.64,
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_normalTraction[6] = {
  -2.2e+5,
  -2.2e+5,
  -2.2e+5,
  -2.3e+5,
  -2.3e+5,
  -2.3e+5,
};

const PylithScalar pylith::friction::SlipWeakeningTimeStableData::_stateVarsUpdated[6*2] = {
  0.20,
  0.15,

  0.20,
  0.15,

  0.20,
  0.15,

  1.0,
  0.6,

  1.0,
  0.6,

  1.0,
  0.6,
};

pylith::friction::SlipWeakeningTimeStableData::SlipWeakeningTimeStableData(void)
{ // constructor
  numLocs = _numLocs;
  numProperties = _numProperties;
  numStateVars = _numStateVars;
  numDBProperties = _numDBProperties;
  numDBStateVars = _numDBStateVars;
  numPropsVertex = _numPropsVertex;
  numVarsVertex = _numVarsVertex;
  lengthScale = _lengthScale;
  timeScale = _timeScale;
  pressureScale = _pressureScale;
  densityScale = _densityScale;
  numPropertyValues = const_cast<int*>(_numPropertyValues);
  numStateVarValues = const_cast<int*>(_numStateVarValues);
  dbPropertyValues = const_cast<char**>(_dbPropertyValues);
  dbStateVarValues = const_cast<char**>(_dbStateVarValues);
  dbProperties = const_cast<PylithScalar*>(_dbProperties);
  dbStateVars = const_cast<PylithScalar*>(_dbStateVars);
  dt = _dt;
  properties = const_cast<PylithScalar*>(_properties);
  stateVars = const_cast<PylithScalar*>(_stateVars);
  propertiesNondim = const_cast<PylithScalar*>(_propertiesNondim);
  stateVarsNondim = const_cast<PylithScalar*>(_stateVarsNondim);
  friction = const_cast<PylithScalar*>(_friction);
  frictionDeriv = const_cast<PylithScalar*>(_frictionDeriv);
  slip = const_cast<PylithScalar*>(_slip);
  slipRate = const_cast<PylithScalar*>(_slipRate);
  normalTraction = const_cast<PylithScalar*>(_normalTraction);
  stateVarsUpdated = const_cast<PylithScalar*>(_stateVarsUpdated);
} // constructor

pylith::friction::SlipWeakeningTimeStableData::~SlipWeakeningTimeStableData(void)
{}


// End of file
