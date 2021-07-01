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

#include "StaticFrictionData.hh"

const int pylith::friction::StaticFrictionData::_numLocs = 2;

const int pylith::friction::StaticFrictionData::_numProperties = 2;

const int pylith::friction::StaticFrictionData::_numStateVars = 0;

const int pylith::friction::StaticFrictionData::_numDBProperties = 2;

const int pylith::friction::StaticFrictionData::_numDBStateVars = 0;

const int pylith::friction::StaticFrictionData::_numPropsVertex = 2;

const int pylith::friction::StaticFrictionData::_numVarsVertex = 0;

const PylithScalar pylith::friction::StaticFrictionData::_lengthScale =   1.00000000e+03;

const PylithScalar pylith::friction::StaticFrictionData::_timeScale =   1.00000000e+00;

const PylithScalar pylith::friction::StaticFrictionData::_pressureScale =   2.25000000e+10;

const PylithScalar pylith::friction::StaticFrictionData::_densityScale =   1.00000000e+03;

const PylithScalar pylith::friction::StaticFrictionData::_dt = 0.01;

const int pylith::friction::StaticFrictionData::_numPropertyValues[] = {
1,
1,
};

const int* pylith::friction::StaticFrictionData::_numStateVarValues = 0;

const char* pylith::friction::StaticFrictionData::_dbPropertyValues[] = {
"friction-coefficient",
"cohesion",
};

const char** pylith::friction::StaticFrictionData::_dbStateVarValues = 0;

const PylithScalar pylith::friction::StaticFrictionData::_dbProperties[] = {
  0.6,
  1000000,
  0.6,
  1000000,
};

const PylithScalar* pylith::friction::StaticFrictionData::_dbStateVars = 0;

const PylithScalar pylith::friction::StaticFrictionData::_properties[] = {
  0.6,
  1000000,
  0.6,
  1000000,
};

const PylithScalar* pylith::friction::StaticFrictionData::_stateVars = 0;

const PylithScalar pylith::friction::StaticFrictionData::_propertiesNondim[] = {
   0.6,
   0.000044444444,
   0.6,
   0.000044444444,
};

const PylithScalar* pylith::friction::StaticFrictionData::_stateVarsNondim = 0;

const PylithScalar pylith::friction::StaticFrictionData::_friction[] = {
  1000001.32,
  1.0e+6,
};

const PylithScalar pylith::friction::StaticFrictionData::_frictionDeriv[] = {
  0.0,
  0.0,
};

const PylithScalar pylith::friction::StaticFrictionData::_slip[] = {
  0.12,
  0.25,
};

const PylithScalar pylith::friction::StaticFrictionData::_slipRate[] = {
  0.74,
  0.64,
};

const PylithScalar pylith::friction::StaticFrictionData::_normalTraction[] = {
  -2.2,
  0.8,
};

const PylithScalar* pylith::friction::StaticFrictionData::_stateVarsUpdated = 0;

pylith::friction::StaticFrictionData::StaticFrictionData(void)
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

pylith::friction::StaticFrictionData::~StaticFrictionData(void)
{}


// End of file
