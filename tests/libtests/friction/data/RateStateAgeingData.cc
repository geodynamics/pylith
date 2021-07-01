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

#include "RateStateAgeingData.hh"

const int pylith::friction::RateStateAgeingData::_numLocs = 2;

const int pylith::friction::RateStateAgeingData::_numProperties = 6;

const int pylith::friction::RateStateAgeingData::_numStateVars = 1;

const int pylith::friction::RateStateAgeingData::_numDBProperties = 6;

const int pylith::friction::RateStateAgeingData::_numDBStateVars = 1;

const int pylith::friction::RateStateAgeingData::_numPropsVertex = 6;

const int pylith::friction::RateStateAgeingData::_numVarsVertex = 1;

const PylithScalar pylith::friction::RateStateAgeingData::_lengthScale =   1.00000000e+03;

const PylithScalar pylith::friction::RateStateAgeingData::_timeScale =   1.00000000e+00;

const PylithScalar pylith::friction::RateStateAgeingData::_pressureScale =   2.25000000e+10;

const PylithScalar pylith::friction::RateStateAgeingData::_densityScale =   1.00000000e+03;

const PylithScalar pylith::friction::RateStateAgeingData::_dt = 0.01;

const int pylith::friction::RateStateAgeingData::_numPropertyValues[] = {
  1,
  1,
  1,
  1,
  1,
  1,
};

const int pylith::friction::RateStateAgeingData::_numStateVarValues[] = {
  1,
};

const char* pylith::friction::RateStateAgeingData::_dbPropertyValues[] = {
  "reference-friction-coefficient",
  "reference-slip-rate",
  "characteristic-slip-distance",
  "constitutive-parameter-a",
  "constitutive-parameter-b",
  "cohesion",
};

const char* pylith::friction::RateStateAgeingData::_dbStateVarValues[] = {
  "state-variable",
};

const PylithScalar pylith::friction::RateStateAgeingData::_dbProperties[] = {
  0.6,
  0.000001,
  0.0370,
  0.0125,
  0.0172,
  1000000,
  0.5,
  0.000002,
  0.0470,
  0.0225,
  0.0272,
  1000000,
};

const PylithScalar pylith::friction::RateStateAgeingData::_dbStateVars[] = {
  92.7,
  93.7,
};

const PylithScalar pylith::friction::RateStateAgeingData::_properties[] = {
  0.6,
  0.000001,
  0.0370,
  0.0125,
  0.0172,
  1000000,
  0.5,
  0.000002,
  0.0470,
  0.0225,
  0.0272,
  1000000,
};

const PylithScalar pylith::friction::RateStateAgeingData::_stateVars[] = {
  92.7,
  93.7,
};

const PylithScalar pylith::friction::RateStateAgeingData::_propertiesNondim[] = {
  0.6,
  0.000000001,
  0.0000370,
  0.0125,
  0.0172,
  0.000044444444,
  0.5,
  0.000000002,
  0.0000470,
  0.0225,
  0.0272,
  0.000044444444,
};

const PylithScalar pylith::friction::RateStateAgeingData::_stateVarsNondim[] = {
  92.7,
  93.7,
};

const PylithScalar pylith::friction::RateStateAgeingData::_friction[] = {
  1000001.285949009547604,
  1000001.164378652801948,
};

const PylithScalar pylith::friction::RateStateAgeingData::_frictionDeriv[] = {
  2.2*0.0125/(0.0011*0.01),
  2.3*0.0225/(0.0021*0.01),
};

const PylithScalar pylith::friction::RateStateAgeingData::_slip[] = {
  0.12,
  0.22,
};

const PylithScalar pylith::friction::RateStateAgeingData::_slipRate[] = {
  0.0011,
  0.0021,
};

const PylithScalar pylith::friction::RateStateAgeingData::_normalTraction[] = {
  -2.2,
  -2.3,
};

const PylithScalar pylith::friction::RateStateAgeingData::_stateVarsUpdated[] = {
  92.682443150471812,
  93.668141160483529,
};

pylith::friction::RateStateAgeingData::RateStateAgeingData(void)
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

pylith::friction::RateStateAgeingData::~RateStateAgeingData(void)
{}


// End of file
