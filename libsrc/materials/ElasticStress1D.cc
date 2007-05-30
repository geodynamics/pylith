// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "ElasticStress1D.hh" // implementation of object methods

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticStress1D {

      /// Number of entries in stress tensor.
      const int tensorSize = 1;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 1;

      /// Values expected in spatial database
      const int numDBValues = 3;
      const char* namesDBValues[] =
	{"density", "vs", "vp" };

      /// Indices (order) of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;

      /// Parameters
      const int numParameters = 3;
      const char* namesParameters[] = {"density", "mu", "lambda" };
      const int numParamValues[] = { 1, 1, 1 };
      
      /// Indices (order) of parameters
      const int pidDensity = 0;
      const int pidMu = 1;
      const int pidLambda = 2;

    } // _ElasticStress1D
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticStress1D::ElasticStress1D(void)
{ // constructor
  _dimension = 1;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticStress1D::~ElasticStress1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get names of values expected to be in database of parameters for
const char**
pylith::materials::ElasticStress1D::_dbValues(void) const
{ // _dbValues
  return _ElasticStress1D::namesDBValues;
} // _dbValues

// ----------------------------------------------------------------------
// Get number of values expected to be in database of parameters for
int
pylith::materials::ElasticStress1D::_numDBValues(void) const
{ // _numDBValues
  return _ElasticStress1D::numDBValues;
} // _numDBValues

// ----------------------------------------------------------------------
// Get names of parameters for physical properties.
const char**
pylith::materials::ElasticStress1D::_parameterNames(void) const
{ // _parameterNames
  return _ElasticStress1D::namesParameters;
} // _parameterNames

// ----------------------------------------------------------------------
// Get number of values for each parameter for physical properties.
void
pylith::materials::ElasticStress1D::_numParamValues(int_array* numValues) const
{ // _numParamValues
  assert(0 != numValues);

  const int numParams = _ElasticStress1D::numParameters;
  numValues->resize(numParams);
  for (int i=0; i < numParams; ++i)
    (*numValues)[i] = _ElasticStress1D::numParamValues[i];
} // _numParamValues

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticStress1D::_dbToParameters(std::vector<double_array>* const paramVals,
					   const double_array& dbValues) const
{ // computeParameters
  assert(0 != paramVals);
  const int numParams = paramVals->size();
  assert(_ElasticStress1D::numParameters == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticStress1D::numDBValues == numDBValues);
  for (int i=0; i < numParams; ++i)
    assert(_ElasticStress1D::numParamValues[i] == (*paramVals)[i].size());

  const double density = dbValues[_ElasticStress1D::didDensity];
  const double vs = dbValues[_ElasticStress1D::didVs];
  const double vp = dbValues[_ElasticStress1D::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  (*paramVals)[_ElasticStress1D::pidDensity][0] = density;
  (*paramVals)[_ElasticStress1D::pidMu][0] = mu;
  (*paramVals)[_ElasticStress1D::pidLambda][0] = lambda;
} // computeParameters

// ----------------------------------------------------------------------
// Get number of entries in stress tensor.
int
pylith::materials::ElasticStress1D::_tensorSize(void) const
{ // _tensorSize
  return _ElasticStress1D::tensorSize;
} // _tensorSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
pylith::materials::ElasticStress1D::_numElasticConsts(void) const
{ // _numElasticConsts
  return _ElasticStress1D::numElasticConsts;
} // _numElasticConsts

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticStress1D::_calcDensity(double_array* const density,
				  const std::vector<double_array>& parameters)
{ // _calcDensity
  assert(0 != density);
  assert(1 == density->size());
  assert(_ElasticStress1D::numParameters == parameters.size());
  assert(1 == parameters[_ElasticStress1D::pidDensity].size());

  (*density)[0] = parameters[_ElasticStress1D::pidDensity][0];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticStress1D::_calcStress(double_array* const stress,
				 const std::vector<double_array>& parameters,
				 const double_array& totalStrain)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticStress1D::tensorSize == stress->size());
  assert(_ElasticStress1D::numParameters == parameters.size());
  assert(_ElasticStress1D::tensorSize == totalStrain.size());

  const double density = parameters[_ElasticStress1D::pidDensity][0];
  const double mu = parameters[_ElasticStress1D::pidMu][0];
  const double lambda = parameters[_ElasticStress1D::pidLambda][0];

  const double e11 = totalStrain[0];
  (*stress)[0] = mu * (3.0*lambda+2.0*mu) / (lambda + mu) * e11;
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::ElasticStress1D::_calcElasticConsts(
				  double_array* const elasticConsts,
				  const std::vector<double_array>& parameters,
				  const double_array& totalStrain)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticStress1D::numElasticConsts == elasticConsts->size());
  assert(_ElasticStress1D::numParameters == parameters.size());
  assert(_ElasticStress1D::tensorSize == totalStrain.size());
 
  const double density = parameters[_ElasticStress1D::pidDensity][0];
  const double mu = parameters[_ElasticStress1D::pidMu][0];
  const double lambda = parameters[_ElasticStress1D::pidLambda][0];

  (*elasticConsts)[0] = mu * (3.0*lambda+2.0*mu) / (lambda + mu);
} // _calcElasticConsts


// End of file 
