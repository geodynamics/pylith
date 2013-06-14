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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "SlipWeakeningStressData.hh"

const int pylith::friction::SlipWeakeningStressData::_numLocs = 2;

const int pylith::friction::SlipWeakeningStressData::_numProperties = 4;

const int pylith::friction::SlipWeakeningStressData::_numStateVars = 2;

const int pylith::friction::SlipWeakeningStressData::_numDBProperties = 4;

const int pylith::friction::SlipWeakeningStressData::_numDBStateVars = 2;

const int pylith::friction::SlipWeakeningStressData::_numPropsVertex = 4;

const int pylith::friction::SlipWeakeningStressData::_numVarsVertex = 2;

const PylithScalar pylith::friction::SlipWeakeningStressData::_lengthScale =   1.00000000e+03;

const PylithScalar pylith::friction::SlipWeakeningStressData::_timeScale =   1.00000000e+01;

const PylithScalar pylith::friction::SlipWeakeningStressData::_pressureScale =   2.5000000e+10;

const PylithScalar pylith::friction::SlipWeakeningStressData::_densityScale =   1.00000000e+03;

const PylithScalar pylith::friction::SlipWeakeningStressData::_dt = 0.01;

const int pylith::friction::SlipWeakeningStressData::_numPropertyValues[] = {
  1,
  1,
  1,
  1,
};

const int pylith::friction::SlipWeakeningStressData::_numStateVarValues[] = {
  1,
  1,
};

const char* pylith::friction::SlipWeakeningStressData::_dbPropertyValues[] = {
  "static-stress",
  "dynamic-stress",
  "slip-weakening-parameter",
  "weakening-time",
};

const char* pylith::friction::SlipWeakeningStressData::_dbStateVarValues[] = {
  "cumulative-slip",
  "previous-slip",
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_dbProperties[] = {
  2.0e+6,
  1.0e+6,
  0.8,
  2.0,
  2.5e+6,
  1.25e+6,
  0.4,
  1.0,
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_dbStateVars[] = {
  0.4,
  0.2,
  0.5,
  0.1,
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_properties[] = {
  2.0e+6,
  1.0e+6,
  0.8,
  2.0,
  2.5e+6,
  1.25e+6,
  0.4,
  1.0,
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_stateVars[] = {
  0.4,
  0.2,
  0.5,
  0.1,
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_propertiesNondim[] = {
  2.0e+6/2.5e+10,
  1.0e+6/2.5e+10,
  0.0008,
  0.2,
  2.5e+6/2.5e+10,
  1.25e+6/2.5e+10,
  0.0004,
  0.1,
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_stateVarsNondim[] = {
  0.0004,
  0.0002,
  0.0005,
  0.0001,
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_friction[] = {
  1.5e+6,
  1.25e+6,
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_frictionDeriv[] = {
  -1.0e+6 / 0.8,
  0.0,
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_slip[] = {
  0.12,
  0.25,
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_slipRate[] = {
  0.74,
  0.64,
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_normalTraction[] = {
  -2.2e+5,
  -2.3e+5,
};

const PylithScalar pylith::friction::SlipWeakeningStressData::_stateVarsUpdated[] = {
  0.48,
  0.12,
  0.65,
  0.25,
};

pylith::friction::SlipWeakeningStressData::SlipWeakeningStressData(void)
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

pylith::friction::SlipWeakeningStressData::~SlipWeakeningStressData(void)
{}


// End of file
