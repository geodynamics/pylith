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

#include "SlipWeakeningTimeData.hh"

const int pylith::friction::SlipWeakeningTimeData::_numLocs = 2;

const int pylith::friction::SlipWeakeningTimeData::_numProperties = 5;

const int pylith::friction::SlipWeakeningTimeData::_numStateVars = 2;

const int pylith::friction::SlipWeakeningTimeData::_numDBProperties = 5;

const int pylith::friction::SlipWeakeningTimeData::_numDBStateVars = 2;

const int pylith::friction::SlipWeakeningTimeData::_numPropsVertex = 5;

const int pylith::friction::SlipWeakeningTimeData::_numVarsVertex = 2;

const PylithScalar pylith::friction::SlipWeakeningTimeData::_lengthScale =   1.00000000e+03;

const PylithScalar pylith::friction::SlipWeakeningTimeData::_timeScale =   1.00000000e+01;

const PylithScalar pylith::friction::SlipWeakeningTimeData::_pressureScale =   2.25000000e+10;

const PylithScalar pylith::friction::SlipWeakeningTimeData::_densityScale =   1.00000000e+03;

const PylithScalar pylith::friction::SlipWeakeningTimeData::_dt = 0.01;

const int pylith::friction::SlipWeakeningTimeData::_numPropertyValues[] = {
  1,
  1,
  1,
  1,
  1,
};

const int pylith::friction::SlipWeakeningTimeData::_numStateVarValues[] = {
  1,
  1,
};

const char* pylith::friction::SlipWeakeningTimeData::_dbPropertyValues[] = {
  "static-coefficient",
  "dynamic-coefficient",
  "slip-weakening-parameter",
  "cohesion",
  "weakening-time",
};

const char* pylith::friction::SlipWeakeningTimeData::_dbStateVarValues[] = {
  "cumulative-slip",
  "previous-slip",
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_dbProperties[] = {
  0.6,
  0.5,
  0.8,
  1000000,
  2.0,
  0.6,
  0.5,
  0.4,
  1000000,
  1.0,
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_dbStateVars[] = {
  0.4,
  0.2,
  0.5,
  0.1,
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_properties[] = {
  0.6,
  0.5,
  0.8,
  1000000,
  2.0,
  0.6,
  0.5,
  0.4,
  1000000,
  1.0,
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_stateVars[] = {
  0.4,
  0.2,
  0.5,
  0.1,
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_propertiesNondim[] = {
  0.6,
  0.5,
  0.0008,
  0.000044444444,
  0.2,
  0.6,
  0.5,
  0.0004,
  0.000044444444,
  0.1,
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_stateVarsNondim[] = {
  0.0004,
  0.0002,
  0.0005,
  0.0001,
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_friction[] = {
  1000001.21,
  1000001.15,
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_frictionDeriv[] = {
  -2.2*0.1/0.8,
  0.0,
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_slip[] = {
  0.12,
  0.25,
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_slipRate[] = {
  0.74,
  0.64,
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_normalTraction[] = {
  -2.2,
  -2.3,
};

const PylithScalar pylith::friction::SlipWeakeningTimeData::_stateVarsUpdated[] = {
  0.48,
  0.12,
  0.65,
  0.25,
};

pylith::friction::SlipWeakeningTimeData::SlipWeakeningTimeData(void)
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

pylith::friction::SlipWeakeningTimeData::~SlipWeakeningTimeData(void)
{}


// End of file
