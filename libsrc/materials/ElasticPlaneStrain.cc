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

#include "ElasticPlaneStrain.hh" // implementation of object methods

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticPlaneStrain {

      /// Number of entries in stress tensor.
      const int tensorSize = 3;

      /// Number of elastic constants (for general 3-D elastic material)
      const int numElasticConsts = 6;

      /// Values expected in spatial database
      const int numDBValues = 3;
      const int numParamValues[] = { 1, 1, 1 };
      const char* namesDBValues[] = { "density", "vs", "vp" };
      
      /// Indices (order) of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;
      
      /// Parameters
      const int numParameters = 3;
      const char* namesParameters[] = { "density", "mu", "lambda" };
      
      /// Indices (order) of parameters
      const int pidDensity = 0;
      const int pidMu = 1;
      const int pidLambda = 2;

    } // _ElasticPlaneStrain
  } // materials
} // pylith


// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticPlaneStrain::ElasticPlaneStrain(void)
{ // constructor
  _dimension = 2;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticPlaneStrain::~ElasticPlaneStrain(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get names of values expected to be in database of parameters for
const char**
pylith::materials::ElasticPlaneStrain::_dbValues(void) const
{ // _dbValues
  return _ElasticPlaneStrain::namesDBValues;
} // _dbValues

// ----------------------------------------------------------------------
// Get number of values expected to be in database of parameters for
int
pylith::materials::ElasticPlaneStrain::_numDBValues(void) const
{ // _numDBValues
  return _ElasticPlaneStrain::numDBValues;
} // _numDBValues

// ----------------------------------------------------------------------
// Get names of parameters for physical properties.
const char**
pylith::materials::ElasticPlaneStrain::_parameterNames(void) const
{ // _parameterNames
  return _ElasticPlaneStrain::namesParameters;
} // _parameterNames

// ----------------------------------------------------------------------
// Get number of values for each parameter for physical properties.
void
pylith::materials::ElasticPlaneStrain::_numParamValues(int_array* numValues) const
{ // _numParamValues
  assert(0 != numValues);

  const int numParams = _ElasticPlaneStrain::numParameters;
  numValues->resize(numParams);
  for (int i=0; i < numParams; ++i)
    (*numValues)[i] = _ElasticPlaneStrain::numParamValues[i];
} // _numParamValues

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticPlaneStrain::_dbToParameters(
				   std::vector<double_array>* const paramVals,
				   const double_array& dbValues) const
{ // computeParameters
  assert(0 != paramVals);
  const int numParams = paramVals->size();
  assert(_ElasticPlaneStrain::numParameters == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticPlaneStrain::numDBValues == numDBValues);
  for (int i=0; i < numParams; ++i)
    assert(_ElasticPlaneStrain::numParamValues[i] == (*paramVals)[i].size());

  const double density = dbValues[_ElasticPlaneStrain::didDensity];
  const double vs = dbValues[_ElasticPlaneStrain::didVs];
  const double vp = dbValues[_ElasticPlaneStrain::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  (*paramVals)[_ElasticPlaneStrain::pidDensity][0] = density;
  (*paramVals)[_ElasticPlaneStrain::pidMu][0] = mu;
  (*paramVals)[_ElasticPlaneStrain::pidLambda][0] = lambda;
} // computeParameters

// ----------------------------------------------------------------------
// Get number of entries in stress tensor.
int
pylith::materials::ElasticPlaneStrain::_tensorSize(void) const
{ // _tensorSize
  return _ElasticPlaneStrain::tensorSize;
} // _tensorSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
pylith::materials::ElasticPlaneStrain::_numElasticConsts(void) const
{ // _numElasticConsts
  return _ElasticPlaneStrain::numElasticConsts;
} // _numElasticConsts

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStrain::_calcDensity(
				  double_array* const density,
				  const std::vector<double_array>& parameters)
{ // calcDensity
  assert(0 != density);
  assert(1 == density->size());
  assert(_ElasticPlaneStrain::numParameters == parameters.size());
  assert(1 == parameters[_ElasticPlaneStrain::pidDensity].size());

  (*density)[0] = parameters[_ElasticPlaneStrain::pidDensity][0];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticPlaneStrain::_calcStress(
				  double_array* const stress,
				  const std::vector<double_array>& parameters,
				  const double_array& totalStrain)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticPlaneStrain::tensorSize == stress->size());
  assert(_ElasticPlaneStrain::numParameters == parameters.size());
  assert(_ElasticPlaneStrain::tensorSize == totalStrain.size());
  assert(1 == parameters[_ElasticPlaneStrain::pidDensity].size());
  assert(1 == parameters[_ElasticPlaneStrain::pidMu].size());
  assert(1 == parameters[_ElasticPlaneStrain::pidLambda].size());

  const double density = parameters[_ElasticPlaneStrain::pidDensity][0];
  const double mu = parameters[_ElasticPlaneStrain::pidMu][0];
  const double lambda = parameters[_ElasticPlaneStrain::pidLambda][0];
  
  const double lambda2mu = lambda + 2.0 * mu;
  
  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e12 = totalStrain[2];

  const double s12 = lambda * (e11 + e22);

  (*stress)[0] = s12 + 2.0*mu*e11;
  (*stress)[1] = s12 + 2.0*mu*e22;
  (*stress)[2] = 2.0 * mu * e12;
} // _calcStress

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStrain::_calcElasticConsts(
				   double_array* const elasticConsts,
				   const std::vector<double_array>& parameters,
				   const double_array& totalStrain)
{ // calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticPlaneStrain::numElasticConsts == elasticConsts->size());
  assert(_ElasticPlaneStrain::numParameters == parameters.size());
  assert(_ElasticPlaneStrain::tensorSize == totalStrain.size());
  assert(1 == parameters[_ElasticPlaneStrain::pidDensity].size());
  assert(1 == parameters[_ElasticPlaneStrain::pidMu].size());
  assert(1 == parameters[_ElasticPlaneStrain::pidLambda].size());

  const double density = parameters[_ElasticPlaneStrain::pidDensity][0];
  const double mu = parameters[_ElasticPlaneStrain::pidMu][0];
  const double lambda = parameters[_ElasticPlaneStrain::pidLambda][0];

  const double lambda2mu = lambda + 2.0*mu;
  
  (*elasticConsts)[0] = lambda2mu; // C1111
  (*elasticConsts)[1] = lambda; // C1122
  (*elasticConsts)[2] = 0; // C1112
  (*elasticConsts)[3] = lambda2mu; // C2222
  (*elasticConsts)[4] = 0; // C2212
  (*elasticConsts)[5] = 2.0*mu; // C1212
} // calcElasticConsts


// End of file 
