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

#include "ElasticIsotropic3DData.hh"

const int pylith::materials::ElasticIsotropic3DData::_numDBValues = 3;

const int pylith::materials::ElasticIsotropic3DData::_numParameters = 1;

const char* pylith::materials::ElasticIsotropic3DData::_dbValues[] =
  { "density", "vp", "vs" };

const char* pylith::materials::ElasticIsotropic3DData::_parameterNames[] =
  { "density", "mu", "lambda" };


pylith::materials::ElasticIsotropic3DData::ElasticIsotropic3DData(void)
{ // constructor
  numDBValues = _numDBValues;
  numParameters = _numParameters;
  dbValues = const_cast<char**>(_dbValues);
  parameterNames = const_cast<char**>(_parameterNames);
} // constructor

pylith::materials::ElasticIsotropic3DData::~ElasticIsotropic3DData(void)
{}


// End of file
