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

#include "TimeWeakeningData.hh"

const int pylith::friction::TimeWeakeningData::_numLocs = 2;

const int pylith::friction::TimeWeakeningData::_numProperties = 4;

const int pylith::friction::TimeWeakeningData::_numStateVars = 1;

const int pylith::friction::TimeWeakeningData::_numDBProperties = 4;

const int pylith::friction::TimeWeakeningData::_numDBStateVars = 1;

const int pylith::friction::TimeWeakeningData::_numPropsVertex = 4;

const int pylith::friction::TimeWeakeningData::_numVarsVertex = 1;

const double pylith::friction::TimeWeakeningData::_lengthScale =   1.00000000e+03;

const double pylith::friction::TimeWeakeningData::_timeScale =   1.00000000e+00;

const double pylith::friction::TimeWeakeningData::_pressureScale =   2.25000000e+10;

const double pylith::friction::TimeWeakeningData::_densityScale =   1.00000000e+03;

const double pylith::friction::TimeWeakeningData::_dt = 0.01;

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
  "time-weakeneing-parameter",
  "cohesion",
};

const char* pylith::friction::TimeWeakeningData::_dbStateVarValues[] = {
  "Elapsed-slip",
};

const double pylith::friction::TimeWeakeningData::_dbProperties[] = {
  0.6,
  0.5,
  0.8,
  1000000,
  0.6,
  0.5,
  0.4,
  1000000,
};

const double pylith::friction::TimeWeakeningData::_dbStateVars[] = {
  0.4,
  0.5,
};

const double pylith::friction::TimeWeakeningData::_properties[] = {
  0.6,
  0.5,
  0.8,
  1000000,
  0.6,
  0.5,
  0.4,
  1000000,
};

const double pylith::friction::TimeWeakeningData::_stateVars[] = {
  0.4,
  0.5,
};

const double pylith::friction::TimeWeakeningData::_propertiesNondim[] = {
  0.6,
  0.5,
  0.8,
  0.000044444444,
  0.6,
  0.5,
  0.4,
  0.000044444444,
};

const double pylith::friction::TimeWeakeningData::_stateVarsNondim[] = {
  0.4,
  0.5,
};

const double pylith::friction::TimeWeakeningData::_friction[] = {
  1000001.21,
  1000001.15,
};

const double pylith::friction::TimeWeakeningData::_slip[] = {
  0.12,
  0.25,
};

const double pylith::friction::TimeWeakeningData::_slipRate[] = {
  0.74,
  0.64,
};

const double pylith::friction::TimeWeakeningData::_normalTraction[] = {
  -2.2,
  -2.3,
};

const double pylith::friction::TimeWeakeningData::_stateVarsUpdated[] = {
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
  dbProperties = const_cast<double*>(_dbProperties);
  dbStateVars = const_cast<double*>(_dbStateVars);
  dt = _dt;
  properties = const_cast<double*>(_properties);
  stateVars = const_cast<double*>(_stateVars);
  propertiesNondim = const_cast<double*>(_propertiesNondim);
  stateVarsNondim = const_cast<double*>(_stateVarsNondim);
  friction = const_cast<double*>(_friction);
  slip = const_cast<double*>(_slip);
  slipRate = const_cast<double*>(_slipRate);
  normalTraction = const_cast<double*>(_normalTraction);
  stateVarsUpdated = const_cast<double*>(_stateVarsUpdated);
} // constructor

pylith::friction::TimeWeakeningData::~TimeWeakeningData(void)
{}


// End of file
