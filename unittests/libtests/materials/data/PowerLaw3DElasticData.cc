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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application powerlaw3delastic.

#include "PowerLaw3DElasticData.hh"

const int pylith::materials::PowerLaw3DElasticData::_dimension = 3;

const int pylith::materials::PowerLaw3DElasticData::_numLocs = 2;

const int pylith::materials::PowerLaw3DElasticData::_numProperties = 6;

const int pylith::materials::PowerLaw3DElasticData::_numStateVars = 2;

const int pylith::materials::PowerLaw3DElasticData::_numDBProperties = 6;

const int pylith::materials::PowerLaw3DElasticData::_numDBStateVars = 12;

const int pylith::materials::PowerLaw3DElasticData::_numPropsQuadPt = 6;

const int pylith::materials::PowerLaw3DElasticData::_numVarsQuadPt = 12;

const double pylith::materials::PowerLaw3DElasticData::_lengthScale =   1.00000000e+03;

const double pylith::materials::PowerLaw3DElasticData::_timeScale =   1.00000000e+00;

const double pylith::materials::PowerLaw3DElasticData::_pressureScale =   2.25000000e+10;

const double pylith::materials::PowerLaw3DElasticData::_densityScale =   1.00000000e+03;

const double pylith::materials::PowerLaw3DElasticData::_dtStableImplicit =   4.44444444e+06;

const int pylith::materials::PowerLaw3DElasticData::_numPropertyValues[] = {
1,
1,
1,
1,
1,
1,
};

const int pylith::materials::PowerLaw3DElasticData::_numStateVarValues[] = {
6,
6,
};

const char* pylith::materials::PowerLaw3DElasticData::_dbPropertyValues[] = {
"density",
"vs",
"vp",
"reference-strain-rate",
"reference-stress",
"power-law-exponent",
};

const char* pylith::materials::PowerLaw3DElasticData::_dbStateVarValues[] = {
"viscous-strain-xx",
"viscous-strain-yy",
"viscous-strain-zz",
"viscous-strain-xy",
"viscous-strain-yz",
"viscous-strain-xz",
"stress-xx",
"stress-yy",
"stress-zz",
"stress-xy",
"stress-yz",
"stress-xz",
};

const double pylith::materials::PowerLaw3DElasticData::_dbProperties[] = {
  2.50000000e+03,
  3.00000000e+03,
  5.19615242e+03,
  1.00000000e-06,
  2.00000000e+12,
  1.00000000e+00,
  2.00000000e+03,
  1.20000000e+03,
  2.07846097e+03,
  1.00000000e-06,
  1.25992105e+08,
  3.00000000e+00,
};

const double pylith::materials::PowerLaw3DElasticData::_dbStateVars[] = {
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

const double pylith::materials::PowerLaw3DElasticData::_properties[] = {
  2.50000000e+03,
  2.25000000e+10,
  2.25000000e+10,
  1.00000000e-06,
  2.00000000e+12,
  1.00000000e+00,
  2.00000000e+03,
  2.88000000e+09,
  2.88000000e+09,
  1.00000000e-06,
  1.25992105e+08,
  3.00000000e+00,
};

const double pylith::materials::PowerLaw3DElasticData::_stateVars[] = {
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

const double pylith::materials::PowerLaw3DElasticData::_propertiesNondim[] = {
  2.50000000e+00,
  1.00000000e+00,
  1.00000000e+00,
  1.00000000e-06,
  8.88888889e+01,
  1.00000000e+00,
  2.00000000e+00,
  1.28000000e-01,
  1.28000000e-01,
  1.00000000e-06,
  5.59964911e-03,
  3.00000000e+00,
};

const double pylith::materials::PowerLaw3DElasticData::_stateVarsNondim[] = {
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

const double pylith::materials::PowerLaw3DElasticData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const double pylith::materials::PowerLaw3DElasticData::_strain[] = {
  1.10000000e-04,
  1.20000000e-04,
  1.30000000e-04,
  1.40000000e-04,
  1.50000000e-04,
  1.60000000e-04,
  4.10000000e-04,
  4.20000000e-04,
  4.30000000e-04,
  4.40000000e-04,
  4.50000000e-04,
  4.60000000e-04,
};

const double pylith::materials::PowerLaw3DElasticData::_stress[] = {
 -2.24790000e+07,
 -2.24780000e+07,
 -2.24770000e+07,
 -8.97600000e+06,
 -8.97500000e+06,
 -8.97400000e+06,
 -2.82900000e+06,
 -2.82800000e+06,
 -2.82700000e+06,
 -1.09800000e+06,
 -1.09700000e+06,
 -1.09600000e+06,
};

const double pylith::materials::PowerLaw3DElasticData::_elasticConsts[] = {
  6.75000000e+10,
  2.25000000e+10,
  2.25000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.25000000e+10,
  6.75000000e+10,
  2.25000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.25000000e+10,
  2.25000000e+10,
  6.75000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  8.64000000e+09,
  2.88000000e+09,
  2.88000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.88000000e+09,
  8.64000000e+09,
  2.88000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.88000000e+09,
  2.88000000e+09,
  8.64000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
};

const double pylith::materials::PowerLaw3DElasticData::_initialStress[] = {
  2.10000000e+04,
  2.20000000e+04,
  2.30000000e+04,
  2.40000000e+04,
  2.50000000e+04,
  2.60000000e+04,
  5.10000000e+04,
  5.20000000e+04,
  5.30000000e+04,
  5.40000000e+04,
  5.50000000e+04,
  5.60000000e+04,
};

const double pylith::materials::PowerLaw3DElasticData::_initialStrain[] = {
  3.10000000e-04,
  3.20000000e-04,
  3.30000000e-04,
  3.40000000e-04,
  3.50000000e-04,
  3.60000000e-04,
  6.10000000e-04,
  6.20000000e-04,
  6.30000000e-04,
  6.40000000e-04,
  6.50000000e-04,
  6.60000000e-04,
};

const double pylith::materials::PowerLaw3DElasticData::_stateVarsUpdated[] = {
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
 -2.24790000e+07,
 -2.24780000e+07,
 -2.24770000e+07,
 -8.97600000e+06,
 -8.97500000e+06,
 -8.97400000e+06,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
 -2.82900000e+06,
 -2.82800000e+06,
 -2.82700000e+06,
 -1.09800000e+06,
 -1.09700000e+06,
 -1.09600000e+06,
};

pylith::materials::PowerLaw3DElasticData::PowerLaw3DElasticData(void)
{ // constructor
  dimension = _dimension;
  numLocs = _numLocs;
  numProperties = _numProperties;
  numStateVars = _numStateVars;
  numDBProperties = _numDBProperties;
  numDBStateVars = _numDBStateVars;
  numPropsQuadPt = _numPropsQuadPt;
  numVarsQuadPt = _numVarsQuadPt;
  lengthScale = _lengthScale;
  timeScale = _timeScale;
  pressureScale = _pressureScale;
  densityScale = _densityScale;
  dtStableImplicit = _dtStableImplicit;
  numPropertyValues = const_cast<int*>(_numPropertyValues);
  numStateVarValues = const_cast<int*>(_numStateVarValues);
  dbPropertyValues = const_cast<char**>(_dbPropertyValues);
  dbStateVarValues = const_cast<char**>(_dbStateVarValues);
  dbProperties = const_cast<double*>(_dbProperties);
  dbStateVars = const_cast<double*>(_dbStateVars);
  properties = const_cast<double*>(_properties);
  stateVars = const_cast<double*>(_stateVars);
  propertiesNondim = const_cast<double*>(_propertiesNondim);
  stateVarsNondim = const_cast<double*>(_stateVarsNondim);
  density = const_cast<double*>(_density);
  strain = const_cast<double*>(_strain);
  stress = const_cast<double*>(_stress);
  elasticConsts = const_cast<double*>(_elasticConsts);
  initialStress = const_cast<double*>(_initialStress);
  initialStrain = const_cast<double*>(_initialStrain);
  stateVarsUpdated = const_cast<double*>(_stateVarsUpdated);
} // constructor

pylith::materials::PowerLaw3DElasticData::~PowerLaw3DElasticData(void)
{}


// End of file
