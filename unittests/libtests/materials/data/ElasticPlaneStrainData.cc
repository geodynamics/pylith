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
// This file was generated from python application elasticplanestrain.

#include "ElasticPlaneStrainData.hh"

const int pylith::materials::ElasticPlaneStrainData::_dimension = 2;

const int pylith::materials::ElasticPlaneStrainData::_numDBValues = 3;

const int pylith::materials::ElasticPlaneStrainData::_numParameters = 3;

const int pylith::materials::ElasticPlaneStrainData::_numLocs = 2;

const char* pylith::materials::ElasticPlaneStrainData::_dbValues[] = {
"density",
"vs",
"vp",
};

const char* pylith::materials::ElasticPlaneStrainData::_parameterNames[] = {
"density",
"mu",
"lambda",
};

const double pylith::materials::ElasticPlaneStrainData::_dbData[] = {
  2.50000000e+03,
  3.00000000e+03,
  5.19615242e+03,
  2.00000000e+03,
  1.20000000e+03,
  2.07846097e+03,
};

const double pylith::materials::ElasticPlaneStrainData::_parameterData[] = {
  2.50000000e+03,
  2.25000000e+10,
  2.25000000e+10,
  2.00000000e+03,
  2.88000000e+09,
  2.88000000e+09,
};

const double pylith::materials::ElasticPlaneStrainData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const double pylith::materials::ElasticPlaneStrainData::_strain[] = {
  1.10000000e-04,
  2.20000000e-04,
  3.30000000e-04,
  1.20000000e-04,
  2.30000000e-04,
  3.40000000e-04,
};

const double pylith::materials::ElasticPlaneStrainData::_stress[] = {
  1.23750000e+07,
  1.73250000e+07,
  1.48500000e+07,
  1.69920000e+06,
  2.33280000e+06,
  1.95840000e+06,
};

const double pylith::materials::ElasticPlaneStrainData::_elasticConsts[] = {
  6.75000000e+10,
  2.25000000e+10,
  0.00000000e+00,
  6.75000000e+10,
  0.00000000e+00,
  4.50000000e+10,
  8.64000000e+09,
  2.88000000e+09,
  0.00000000e+00,
  8.64000000e+09,
  0.00000000e+00,
  5.76000000e+09,
};

pylith::materials::ElasticPlaneStrainData::ElasticPlaneStrainData(void)
{ // constructor
  dimension = _dimension;
  numDBValues = _numDBValues;
  numParameters = _numParameters;
  numLocs = _numLocs;
  dbValues = const_cast<char**>(_dbValues);
  parameterNames = const_cast<char**>(_parameterNames);
  dbData = const_cast<double*>(_dbData);
  parameterData = const_cast<double*>(_parameterData);
  density = const_cast<double*>(_density);
  strain = const_cast<double*>(_strain);
  stress = const_cast<double*>(_stress);
  elasticConsts = const_cast<double*>(_elasticConsts);
} // constructor

pylith::materials::ElasticPlaneStrainData::~ElasticPlaneStrainData(void)
{}


// End of file
