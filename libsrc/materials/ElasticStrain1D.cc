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

#include "ElasticStrain1D.hh" // implementation of object methods

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticStrain1D {

      /// Number of entries in stress tensor.
      const int tensorSize = 1;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 1;

      /// Values expected in spatial database
      const int numDBValues = 2;
      const char* namesDBValues[] = { "density", "vp" };      
      
      /// Indices (order) of database values
      const int didDensity = 0;
      const int didVp = 1;
      
      /// Parameters
      const int numParameters = 2;
      const int numParamValues[] = { 1, 1 };
      const char* namesParameters[] = { "density", "lambda2mu" };
      
      /// Indices (order) of parameters
      const int pidDensity = 0;
      const int pidLambda2mu = 1;
    } // _ElasticStrain1D
  } // materials
} // pylith


// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticStrain1D::ElasticStrain1D(void)
{ // constructor
  _dimension = 1;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticStrain1D::~ElasticStrain1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::ElasticStrain1D::ElasticStrain1D(
						const ElasticStrain1D& m) :
  ElasticMaterial(m)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Get names of values expected to be in database of parameters for
const char**
pylith::materials::ElasticStrain1D::_dbValues(void) const
{ // _dbValues
  return _ElasticStrain1D::namesDBValues;
} // _dbValues

// ----------------------------------------------------------------------
// Get number of values expected to be in database of parameters for
int
pylith::materials::ElasticStrain1D::_numDBValues(void) const
{ // _numDBValues
  return _ElasticStrain1D::numDBValues;
} // _numDBValues

// ----------------------------------------------------------------------
// Get names of parameters for physical properties.
const char**
pylith::materials::ElasticStrain1D::_parameterNames(void) const
{ // _parameterNames
  return _ElasticStrain1D::namesParameters;
} // _parameterNames

// ----------------------------------------------------------------------
// Get number of values for each parameter for physical properties.
void
pylith::materials::ElasticStrain1D::_numParamValues(int_array* numValues) const
{ // _numParamValues
  assert(0 != numValues);
  
  const int numParams = _ElasticStrain1D::numParameters;
  numValues->resize(numParams);
  for (int i=0; i < numParams; ++i)
    (*numValues)[i] = _ElasticStrain1D::numParamValues[i];
} // _numParamValues

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticStrain1D::_dbToParameters(std::vector<double_array>* const paramVals,
					   const double_array& dbValues) const
{ // _dbToParameters
  assert(0 != paramVals);
  const int numParams = paramVals->size();
  assert(_ElasticStrain1D::numParameters == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticStrain1D::numDBValues == numDBValues);
  for (int i=0; i < numParams; ++i)
    assert(_ElasticStrain1D::numParamValues[i] == (*paramVals)[i].size());

  const double density = dbValues[_ElasticStrain1D::didDensity];
  const double vp = dbValues[_ElasticStrain1D::didVp];
  const double lambda2mu = density * vp*vp;

  (*paramVals)[_ElasticStrain1D::pidDensity][0] = density;
  (*paramVals)[_ElasticStrain1D::pidLambda2mu][0] = lambda2mu;
} // _dbToParameters

// ----------------------------------------------------------------------
// Get number of entries in stress tensor.
int
pylith::materials::ElasticStrain1D::_tensorSize(void) const
{ // _tensorSize
  return _ElasticStrain1D::tensorSize;
} // _tensorSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
pylith::materials::ElasticStrain1D::_numElasticConsts(void) const
{ // _numElasticConsts
  return _ElasticStrain1D::numElasticConsts;
} // _numElasticConsts

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticStrain1D::_calcDensity(double_array* const density,
				 const std::vector<double_array>& parameters)
{ // _calcDensity
  assert(0 != density);
  assert(1 == density->size());
  assert(_ElasticStrain1D::numParameters == parameters.size());
  assert(1 == parameters[_ElasticStrain1D::pidDensity].size());

  (*density)[0] = parameters[_ElasticStrain1D::pidDensity][0];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticStrain1D::_calcStress(
				   double_array* const stress,
				   const std::vector<double_array>& parameters,
				   const double_array& totalStrain)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticStrain1D::tensorSize == stress->size());
  assert(_ElasticStrain1D::numParameters == parameters.size());
  assert(_ElasticStrain1D::tensorSize == totalStrain.size());
  assert(1 == parameters[_ElasticStrain1D::pidDensity].size());
  assert(1 == parameters[_ElasticStrain1D::pidLambda2mu].size());

  const double density = parameters[_ElasticStrain1D::pidDensity][0];
  const double lambda2mu = parameters[_ElasticStrain1D::pidLambda2mu][0];

  const double e11 = totalStrain[0];
  (*stress)[0] = lambda2mu * e11;
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::ElasticStrain1D::_calcElasticConsts(
				   double_array* const elasticConsts,
			           const std::vector<double_array>& parameters,
				   const double_array& totalStrain)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticStrain1D::numElasticConsts == elasticConsts->size());
  assert(_ElasticStrain1D::numParameters == parameters.size());
  assert(_ElasticStrain1D::tensorSize == totalStrain.size());
  assert(1 == parameters[_ElasticStrain1D::pidDensity].size());
  assert(1 == parameters[_ElasticStrain1D::pidLambda2mu].size());
 
  const double density = parameters[_ElasticStrain1D::pidDensity][0];
  const double lambda2mu = parameters[_ElasticStrain1D::pidLambda2mu][0];

  (*elasticConsts)[0] = lambda2mu;
} // _calcElasticConsts


// End of file 
