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

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticisotropic3d.

#include "ElasticIsotropic3DData.hh"

const int pylith::materials::ElasticIsotropic3DData::_dimension = 3;

const int pylith::materials::ElasticIsotropic3DData::_numLocs = 2;

const int pylith::materials::ElasticIsotropic3DData::_numProperties = 3;

const int pylith::materials::ElasticIsotropic3DData::_numStateVars = 0;

const int pylith::materials::ElasticIsotropic3DData::_numDBProperties = 3;

const int pylith::materials::ElasticIsotropic3DData::_numDBStateVars = 0;

const int pylith::materials::ElasticIsotropic3DData::_numPropsQuadPt = 3;

const int pylith::materials::ElasticIsotropic3DData::_numVarsQuadPt = 0;

const double pylith::materials::ElasticIsotropic3DData::_lengthScale =   1.00000000e+03;

const double pylith::materials::ElasticIsotropic3DData::_timeScale =   1.00000000e+00;

const double pylith::materials::ElasticIsotropic3DData::_pressureScale =   2.25000000e+10;

const double pylith::materials::ElasticIsotropic3DData::_densityScale =   1.00000000e+03;

const double pylith::materials::ElasticIsotropic3DData::_dtStableImplicit =   1.00000000e+30;

const int pylith::materials::ElasticIsotropic3DData::_numPropertyValues[] = {
1,
1,
1,
};

const int* pylith::materials::ElasticIsotropic3DData::_numStateVarValues = 0;

const char* pylith::materials::ElasticIsotropic3DData::_dbPropertyValues[] = {
"density",
"vs",
"vp",
};

const char** pylith::materials::ElasticIsotropic3DData::_dbStateVarValues = 0;

const double pylith::materials::ElasticIsotropic3DData::_dbProperties[] = {
  2.50000000e+03,
  3.00000000e+03,
  5.19615242e+03,
  2.00000000e+03,
  1.20000000e+03,
  2.07846097e+03,
};

const double* pylith::materials::ElasticIsotropic3DData::_dbStateVars = 0;

const double pylith::materials::ElasticIsotropic3DData::_properties[] = {
  2.50000000e+03,
  2.25000000e+10,
  2.25000000e+10,
  2.00000000e+03,
  2.88000000e+09,
  2.88000000e+09,
};

const double* pylith::materials::ElasticIsotropic3DData::_stateVars = 0;

const double pylith::materials::ElasticIsotropic3DData::_propertiesNondim[] = {
  2.50000000e+00,
  1.00000000e+00,
  1.00000000e+00,
  2.00000000e+00,
  1.28000000e-01,
  1.28000000e-01,
};

const double* pylith::materials::ElasticIsotropic3DData::_stateVarsNondim = 0;

const double pylith::materials::ElasticIsotropic3DData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const double pylith::materials::ElasticIsotropic3DData::_strain[] = {
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

const double pylith::materials::ElasticIsotropic3DData::_stress[] = {
 -2.24790000e+07,
 -2.24780000e+07,
 -2.24770000e+07,
 -8.97600000e+06,
  6.61750000e+06,
 -8.97400000e+06,
 -2.82900000e+06,
 -2.82800000e+06,
 -2.82700000e+06,
 -1.09800000e+06,
  2.60956000e+06,
 -1.09600000e+06,
};

const double pylith::materials::ElasticIsotropic3DData::_elasticConsts[] = {
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

const double pylith::materials::ElasticIsotropic3DData::_initialStress[] = {
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

const double pylith::materials::ElasticIsotropic3DData::_initialStrain[] = {
  3.10000000e-04,
  3.20000000e-04,
  3.30000000e-04,
  3.40000000e-04,
  3.50000000e-06,
  3.60000000e-04,
  6.10000000e-04,
  6.20000000e-04,
  6.30000000e-04,
  6.40000000e-04,
  6.50000000e-06,
  6.60000000e-04,
};

const double* pylith::materials::ElasticIsotropic3DData::_stateVarsUpdated = 0;

pylith::materials::ElasticIsotropic3DData::ElasticIsotropic3DData(void)
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

pylith::materials::ElasticIsotropic3DData::~ElasticIsotropic3DData(void)
{}


// End of file
