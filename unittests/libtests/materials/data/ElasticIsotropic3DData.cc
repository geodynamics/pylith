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

const int pylith::materials::ElasticIsotropic3DData::_numDBValues = 3;

const int pylith::materials::ElasticIsotropic3DData::_numParameters = 3;

const int pylith::materials::ElasticIsotropic3DData::_numLocs = 2;

const int pylith::materials::ElasticIsotropic3DData::_numParamValues[] = {
1,
1,
1,
};

const char* pylith::materials::ElasticIsotropic3DData::_dbValues[] = {
"density",
"vs",
"vp",
};

const char* pylith::materials::ElasticIsotropic3DData::_parameterNames[] = {
"density",
"mu",
"lambda",
};

const double pylith::materials::ElasticIsotropic3DData::_dbData[] = {
  2.50000000e+03,
  3.00000000e+03,
  5.19615242e+03,
  2.00000000e+03,
  1.20000000e+03,
  2.07846097e+03,
};

const double pylith::materials::ElasticIsotropic3DData::_parameterData[] = {
  2.50000000e+03,
  2.25000000e+10,
  2.25000000e+10,
  2.00000000e+03,
  2.88000000e+09,
  2.88000000e+09,
};

const double pylith::materials::ElasticIsotropic3DData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const double pylith::materials::ElasticIsotropic3DData::_strain[] = {
  1.10000000e-04,
  2.20000000e-04,
  3.30000000e-04,
  4.40000000e-04,
  5.50000000e-04,
  6.60000000e-04,
  1.20000000e-04,
  2.30000000e-04,
  3.40000000e-04,
  4.50000000e-04,
  5.60000000e-04,
  6.70000000e-04,
};

const double pylith::materials::ElasticIsotropic3DData::_stress[] = {
  1.98000000e+07,
  2.47500000e+07,
  2.97000000e+07,
  1.98000000e+07,
  2.47500000e+07,
  2.97000000e+07,
  2.67840000e+06,
  3.31200000e+06,
  3.94560000e+06,
  2.59200000e+06,
  3.22560000e+06,
  3.85920000e+06,
};

const double pylith::materials::ElasticIsotropic3DData::_elasticConsts[] = {
  6.75000000e+10,
  2.25000000e+10,
  2.25000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  6.75000000e+10,
  2.25000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  6.75000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  0.00000000e+00,
  4.50000000e+10,
  8.64000000e+09,
  2.88000000e+09,
  2.88000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  8.64000000e+09,
  2.88000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  8.64000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
  0.00000000e+00,
  5.76000000e+09,
};

pylith::materials::ElasticIsotropic3DData::ElasticIsotropic3DData(void)
{ // constructor
  dimension = _dimension;
  numDBValues = _numDBValues;
  numParameters = _numParameters;
  numLocs = _numLocs;
  numParamValues = const_cast<int*>(_numParamValues);
  dbValues = const_cast<char**>(_dbValues);
  parameterNames = const_cast<char**>(_parameterNames);
  dbData = const_cast<double*>(_dbData);
  parameterData = const_cast<double*>(_parameterData);
  density = const_cast<double*>(_density);
  strain = const_cast<double*>(_strain);
  stress = const_cast<double*>(_stress);
  elasticConsts = const_cast<double*>(_elasticConsts);
} // constructor

pylith::materials::ElasticIsotropic3DData::~ElasticIsotropic3DData(void)
{}


// End of file
