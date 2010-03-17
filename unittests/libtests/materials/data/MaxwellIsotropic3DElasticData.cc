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
// This file was generated from python application maxwellisotropic3delastic.

#include "MaxwellIsotropic3DElasticData.hh"

const int pylith::materials::MaxwellIsotropic3DElasticData::_dimension = 3;

const int pylith::materials::MaxwellIsotropic3DElasticData::_numLocs = 2;

const int pylith::materials::MaxwellIsotropic3DElasticData::_numProperties = 4;

const int pylith::materials::MaxwellIsotropic3DElasticData::_numStateVars = 2;

const int pylith::materials::MaxwellIsotropic3DElasticData::_numDBProperties = 4;

const int pylith::materials::MaxwellIsotropic3DElasticData::_numDBStateVars = 12;

const int pylith::materials::MaxwellIsotropic3DElasticData::_numPropsQuadPt = 4;

const int pylith::materials::MaxwellIsotropic3DElasticData::_numVarsQuadPt = 12;

const double pylith::materials::MaxwellIsotropic3DElasticData::_lengthScale =   1.00000000e+03;

const double pylith::materials::MaxwellIsotropic3DElasticData::_timeScale =   1.00000000e+00;

const double pylith::materials::MaxwellIsotropic3DElasticData::_pressureScale =   2.25000000e+10;

const double pylith::materials::MaxwellIsotropic3DElasticData::_densityScale =   1.00000000e+03;

const double pylith::materials::MaxwellIsotropic3DElasticData::_dtStableImplicit =   8.88888889e+06;

const int pylith::materials::MaxwellIsotropic3DElasticData::_numPropertyValues[] = {
1,
1,
1,
1,
};

const int pylith::materials::MaxwellIsotropic3DElasticData::_numStateVarValues[] = {
6,
6,
};

const char* pylith::materials::MaxwellIsotropic3DElasticData::_dbPropertyValues[] = {
"density",
"vs",
"vp",
"viscosity",
};

const char* pylith::materials::MaxwellIsotropic3DElasticData::_dbStateVarValues[] = {
"total-strain-xx",
"total-strain-yy",
"total-strain-zz",
"total-strain-xy",
"total-strain-yz",
"total-strain-xz",
"viscous-strain-xx",
"viscous-strain-yy",
"viscous-strain-zz",
"viscous-strain-xy",
"viscous-strain-yz",
"viscous-strain-xz",
};

const double pylith::materials::MaxwellIsotropic3DElasticData::_dbProperties[] = {
  2.50000000e+03,
  3.00000000e+03,
  5.19615242e+03,
  1.00000000e+18,
  2.00000000e+03,
  1.20000000e+03,
  2.07846097e+03,
  1.00000000e+18,
};

const double pylith::materials::MaxwellIsotropic3DElasticData::_dbStateVars[] = {
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

const double pylith::materials::MaxwellIsotropic3DElasticData::_properties[] = {
  2.50000000e+03,
  2.25000000e+10,
  2.25000000e+10,
  4.44444444e+07,
  2.00000000e+03,
  2.88000000e+09,
  2.88000000e+09,
  3.47222222e+08,
};

const double pylith::materials::MaxwellIsotropic3DElasticData::_stateVars[] = {
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

const double pylith::materials::MaxwellIsotropic3DElasticData::_propertiesNondim[] = {
  2.50000000e+00,
  1.00000000e+00,
  1.00000000e+00,
  4.44444444e+07,
  2.00000000e+00,
  1.28000000e-01,
  1.28000000e-01,
  3.47222222e+08,
};

const double pylith::materials::MaxwellIsotropic3DElasticData::_stateVarsNondim[] = {
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

const double pylith::materials::MaxwellIsotropic3DElasticData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const double pylith::materials::MaxwellIsotropic3DElasticData::_strain[] = {
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

const double pylith::materials::MaxwellIsotropic3DElasticData::_stress[] = {
  9.51600000e+06,
  9.92200000e+06,
  1.03280000e+07,
  4.79400000e+06,
  5.20000000e+06,
  5.60600000e+06,
  5.15436000e+06,
  5.20720000e+06,
  5.26004000e+06,
  2.21976000e+06,
  2.27260000e+06,
  2.32544000e+06,
};

const double pylith::materials::MaxwellIsotropic3DElasticData::_elasticConsts[] = {
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

const double pylith::materials::MaxwellIsotropic3DElasticData::_initialStress[] = {
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

const double pylith::materials::MaxwellIsotropic3DElasticData::_initialStrain[] = {
  3.10000000e-05,
  3.20000000e-05,
  3.30000000e-05,
  3.40000000e-05,
  3.50000000e-05,
  3.60000000e-05,
  6.10000000e-05,
  6.20000000e-05,
  6.30000000e-05,
  6.40000000e-05,
  6.50000000e-05,
  6.60000000e-05,
};

const double pylith::materials::MaxwellIsotropic3DElasticData::_stateVarsUpdated[] = {
  1.10000000e-04,
  1.20000000e-04,
  1.30000000e-04,
  1.40000000e-04,
  1.50000000e-04,
  1.60000000e-04,
 -1.00000000e-05,
  1.35525272e-20,
  1.00000000e-05,
  1.40000000e-04,
  1.50000000e-04,
  1.60000000e-04,
  4.10000000e-04,
  4.20000000e-04,
  4.30000000e-04,
  4.40000000e-04,
  4.50000000e-04,
  4.60000000e-04,
 -1.00000000e-05,
  0.00000000e+00,
  1.00000000e-05,
  4.40000000e-04,
  4.50000000e-04,
  4.60000000e-04,
};

pylith::materials::MaxwellIsotropic3DElasticData::MaxwellIsotropic3DElasticData(void)
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

pylith::materials::MaxwellIsotropic3DElasticData::~MaxwellIsotropic3DElasticData(void)
{}


// End of file
