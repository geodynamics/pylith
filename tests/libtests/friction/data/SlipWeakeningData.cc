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

#include "SlipWeakeningData.hh"

const int pylith::friction::SlipWeakeningData::_numLocs = 2;

const int pylith::friction::SlipWeakeningData::_numProperties = 4;

const int pylith::friction::SlipWeakeningData::_numStateVars = 2;

const int pylith::friction::SlipWeakeningData::_numDBProperties = 4;

const int pylith::friction::SlipWeakeningData::_numDBStateVars = 2;

const int pylith::friction::SlipWeakeningData::_numPropsVertex = 4;

const int pylith::friction::SlipWeakeningData::_numVarsVertex = 2;

const PylithScalar pylith::friction::SlipWeakeningData::_lengthScale =   1.00000000e+03;

const PylithScalar pylith::friction::SlipWeakeningData::_timeScale =   1.00000000e+00;

const PylithScalar pylith::friction::SlipWeakeningData::_pressureScale =   2.25000000e+10;

const PylithScalar pylith::friction::SlipWeakeningData::_densityScale =   1.00000000e+03;

const PylithScalar pylith::friction::SlipWeakeningData::_dt = 0.01;

const int pylith::friction::SlipWeakeningData::_numPropertyValues[] = {
  1,
  1,
  1,
  1,
};

const int pylith::friction::SlipWeakeningData::_numStateVarValues[] = {
  1,
  1,
};

const char* pylith::friction::SlipWeakeningData::_dbPropertyValues[] = {
  "static-coefficient",
  "dynamic-coefficient",
  "slip-weakening-parameter",
  "cohesion",
};

const char* pylith::friction::SlipWeakeningData::_dbStateVarValues[] = {
  "cumulative-slip",
  "previous-slip",
};

const PylithScalar pylith::friction::SlipWeakeningData::_dbProperties[] = {
  0.6,
  0.5,
  0.8,
  1000000,
  0.6,
  0.5,
  0.4,
  1000000,
};

const PylithScalar pylith::friction::SlipWeakeningData::_dbStateVars[] = {
  0.4,
  0.2,
  0.5,
  0.1,
};

const PylithScalar pylith::friction::SlipWeakeningData::_properties[] = {
  0.6,
  0.5,
  0.8,
  1000000,
  0.6,
  0.5,
  0.4,
  1000000,
};

const PylithScalar pylith::friction::SlipWeakeningData::_stateVars[] = {
  0.4,
  0.2,
  0.5,
  0.1,
};

const PylithScalar pylith::friction::SlipWeakeningData::_propertiesNondim[] = {
  0.6,
  0.5,
  0.0008,
  0.000044444444,
  0.6,
  0.5,
  0.0004,
  0.000044444444,
};

const PylithScalar pylith::friction::SlipWeakeningData::_stateVarsNondim[] = {
  0.0004,
  0.0002,
  0.0005,
  0.0001,
};

const PylithScalar pylith::friction::SlipWeakeningData::_friction[] = {
  1000001.21,
  1000001.15,
};

const PylithScalar pylith::friction::SlipWeakeningData::_frictionDeriv[] = {
  -2.2*0.1/0.8,
  0.0,
};

const PylithScalar pylith::friction::SlipWeakeningData::_slip[] = {
  0.12,
  0.25,
};

const PylithScalar pylith::friction::SlipWeakeningData::_slipRate[] = {
  0.74,
  0.64,
};

const PylithScalar pylith::friction::SlipWeakeningData::_normalTraction[] = {
  -2.2,
  -2.3,
};

const PylithScalar pylith::friction::SlipWeakeningData::_stateVarsUpdated[] = {
  0.48,
  0.12,
  0.65,
  0.25,
};

pylith::friction::SlipWeakeningData::SlipWeakeningData(void)
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

pylith::friction::SlipWeakeningData::~SlipWeakeningData(void)
{}


// End of file
