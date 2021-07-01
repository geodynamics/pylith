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

#include "TimeWeakeningData.hh"

const int pylith::friction::TimeWeakeningData::_numLocs = 2;

const int pylith::friction::TimeWeakeningData::_numProperties = 4;

const int pylith::friction::TimeWeakeningData::_numStateVars = 1;

const int pylith::friction::TimeWeakeningData::_numDBProperties = 4;

const int pylith::friction::TimeWeakeningData::_numDBStateVars = 1;

const int pylith::friction::TimeWeakeningData::_numPropsVertex = 4;

const int pylith::friction::TimeWeakeningData::_numVarsVertex = 1;

const PylithScalar pylith::friction::TimeWeakeningData::_lengthScale =   1.00000000e+03;

const PylithScalar pylith::friction::TimeWeakeningData::_timeScale =   1.00000000e+00;

const PylithScalar pylith::friction::TimeWeakeningData::_pressureScale =   2.25000000e+10;

const PylithScalar pylith::friction::TimeWeakeningData::_densityScale =   1.00000000e+03;

const PylithScalar pylith::friction::TimeWeakeningData::_dt = 0.01;

const int pylith::friction::TimeWeakeningData::_numPropertyValues[] = {
  1,
  1,
  1,
  1,
};

const int pylith::friction::TimeWeakeningData::_numStateVarValues[] = {
  1,
};

const char* pylith::friction::TimeWeakeningData::_dbPropertyValues[] = {
  "static-coefficient",
  "dynamic-coefficient",
  "time-weakening-parameter",
  "cohesion",
};

const char* pylith::friction::TimeWeakeningData::_dbStateVarValues[] = {
  "elapsed-time",
};

const PylithScalar pylith::friction::TimeWeakeningData::_dbProperties[] = {
  0.6,
  0.5,
  0.8,
  1000000,
  0.6,
  0.5,
  0.4,
  1000000,
};

const PylithScalar pylith::friction::TimeWeakeningData::_dbStateVars[] = {
  0.4,
  0.5,
};

const PylithScalar pylith::friction::TimeWeakeningData::_properties[] = {
  0.6,
  0.5,
  0.8,
  1000000,
  0.6,
  0.5,
  0.4,
  1000000,
};

const PylithScalar pylith::friction::TimeWeakeningData::_stateVars[] = {
  0.4,
  0.5,
};

const PylithScalar pylith::friction::TimeWeakeningData::_propertiesNondim[] = {
  0.6,
  0.5,
  0.8,
  0.000044444444,
  0.6,
  0.5,
  0.4,
  0.000044444444,
};

const PylithScalar pylith::friction::TimeWeakeningData::_stateVarsNondim[] = {
  0.4,
  0.5,
};

const PylithScalar pylith::friction::TimeWeakeningData::_friction[] = {
  1000001.21,
  1000001.15,
};

const PylithScalar pylith::friction::TimeWeakeningData::_frictionDeriv[] = {
  0.0,
  0.0,
};

const PylithScalar pylith::friction::TimeWeakeningData::_slip[] = {
  0.12,
  0.25,
};

const PylithScalar pylith::friction::TimeWeakeningData::_slipRate[] = {
  0.74,
  0.64,
};

const PylithScalar pylith::friction::TimeWeakeningData::_normalTraction[] = {
  -2.2,
  -2.3,
};

const PylithScalar pylith::friction::TimeWeakeningData::_stateVarsUpdated[] = {
  0.41,
  0.51,
};

pylith::friction::TimeWeakeningData::TimeWeakeningData(void)
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

pylith::friction::TimeWeakeningData::~TimeWeakeningData(void)
{}


// End of file
