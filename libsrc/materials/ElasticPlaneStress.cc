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

#include "ElasticPlaneStress.hh" // implementation of object methods

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticPlaneStress {

      /// Number of entries in stress tensor.
      const int tensorSize = 3;

      /// Number of elastic constants (for general 3-D elastic material)
      const int numElasticConsts = 6;

      /// Values expected in spatial database
      const int numDBValues = 3;
      const char* namesDBValues[] = { "density", "vs", "vp" };
      
      /// Indices (order) of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;
	    
      /// Parameters
      const int numParameters = 3;
      const int numParamValues[] = { 1, 1, 1 };
      const char* namesParameters[] = { "density", "mu", "lambda" };
	      
      /// Indices (order) of parameters
      const int pidDensity = 0;
      const int pidMu = 1;
      const int pidLambda = 2;

    } // _ElasticPlaneStress
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticPlaneStress::ElasticPlaneStress(void)
{ // constructor
  _dimension = 2;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticPlaneStress::~ElasticPlaneStress(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get names of values expected to be in database of parameters for
const char**
pylith::materials::ElasticPlaneStress::_dbValues(void) const
{ // _dbValues
  return _ElasticPlaneStress::namesDBValues;
} // _dbValues

// ----------------------------------------------------------------------
// Get number of values expected to be in database of parameters for
int
pylith::materials::ElasticPlaneStress::_numDBValues(void) const
{ // _numDBValues
  return _ElasticPlaneStress::numDBValues;
} // _numDBValues

// ----------------------------------------------------------------------
// Get names of parameters for physical properties.
const char**
pylith::materials::ElasticPlaneStress::_parameterNames(void) const
{ // _parameterNames
  return _ElasticPlaneStress::namesParameters;
} // _parameterNames

// ----------------------------------------------------------------------
// Get number of values for each parameter for physical properties.
void
pylith::materials::ElasticPlaneStress::_numParamValues(int_array* numValues) const
{ // _numParamValues
  assert(0 != numValues);

  const int numParams = _ElasticPlaneStress::numParameters;
  numValues->resize(numParams);
  for (int i=0; i < numParams; ++i)
    (*numValues)[i] = _ElasticPlaneStress::numParamValues[i];
} // _numParamValues

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticPlaneStress::_dbToParameters(
				      std::vector<double_array>* paramVals,
				      const double_array& dbValues) const
{ // _dbToParameters
  assert(0 != paramVals);
  const int numParams = paramVals->size();
  assert(_ElasticPlaneStress::numParameters == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticPlaneStress::numDBValues == numDBValues);
  for (int i=0; i < numParams; ++i)
    assert(_ElasticPlaneStress::numParamValues[i] == (*paramVals)[i].size());

  const double density = dbValues[_ElasticPlaneStress::didDensity];
  const double vs = dbValues[_ElasticPlaneStress::didVs];
  const double vp = dbValues[_ElasticPlaneStress::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;
  if (lambda < 0.0) {
    std::ostringstream msg;
    msg << "Attempted to set Lame's constant lambda to negative value.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if

  (*paramVals)[_ElasticPlaneStress::pidDensity][0] = density;
  (*paramVals)[_ElasticPlaneStress::pidMu][0] = mu;
  (*paramVals)[_ElasticPlaneStress::pidLambda][0] = lambda;
} // _dbToParameters

// ----------------------------------------------------------------------
// Get number of entries in stress tensor.
int
pylith::materials::ElasticPlaneStress::_tensorSize(void) const
{ // _tensorSize
  return _ElasticPlaneStress::tensorSize;
} // _tensorSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
pylith::materials::ElasticPlaneStress::_numElasticConsts(void) const
{ // _numElasticConsts
  return _ElasticPlaneStress::numElasticConsts;
} // _numElasticConsts

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStress::_calcDensity(double_array* const density,
						    const std::vector<double_array>& parameters)
{ // calcDensity
  assert(0 != density);
  assert(1 == density->size());
  assert(_ElasticPlaneStress::numParameters == parameters.size());
  assert(1 == parameters[_ElasticPlaneStress::pidDensity].size());

  (*density)[0] = parameters[_ElasticPlaneStress::pidDensity][0];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticPlaneStress::_calcStress(
				 double_array* const stress,
				 const std::vector<double_array>& parameters,
				 const double_array& totalStrain)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticPlaneStress::tensorSize == stress->size());
  assert(_ElasticPlaneStress::numParameters == parameters.size());
  assert(_ElasticPlaneStress::tensorSize == totalStrain.size());
  assert(1 == parameters[_ElasticPlaneStress::pidDensity].size());
  assert(1 == parameters[_ElasticPlaneStress::pidMu].size());
  assert(1 == parameters[_ElasticPlaneStress::pidLambda].size());

  const double density = parameters[_ElasticPlaneStress::pidDensity][0];
  const double mu = parameters[_ElasticPlaneStress::pidMu][0];
  const double lambda = parameters[_ElasticPlaneStress::pidLambda][0];

  const double lambda2mu = lambda + 2.0 * mu;
  const double lambdamu = lambda + mu;
  
  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e12 = totalStrain[2];

  (*stress)[0] = (4.0*mu*lambdamu * e11 + 2.0*mu*lambda * e22)/lambda2mu;
  (*stress)[1] = (2.0*mu*lambda * e11 + 4.0*mu*lambdamu * e22)/lambda2mu;
  (*stress)[2] = 2.0 * mu * e12;
} // _calcStress

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStress::_calcElasticConsts(
				  double_array* const elasticConsts,
				  const std::vector<double_array>& parameters,
				  const double_array& totalStrain)
{ // calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticPlaneStress::numElasticConsts == elasticConsts->size());
  assert(_ElasticPlaneStress::numParameters == parameters.size());
  assert(_ElasticPlaneStress::tensorSize == totalStrain.size());
  assert(1 == parameters[_ElasticPlaneStress::pidDensity].size());
  assert(1 == parameters[_ElasticPlaneStress::pidMu].size());
  assert(1 == parameters[_ElasticPlaneStress::pidLambda].size());

  const double density = parameters[_ElasticPlaneStress::pidDensity][0];
  const double mu = parameters[_ElasticPlaneStress::pidMu][0];
  const double lambda = parameters[_ElasticPlaneStress::pidLambda][0];
  
  const double lambda2mu = lambda + 2.0 * mu;
  const double c11 = 4.0 * mu * (lambda + mu) / lambda2mu;

  (*elasticConsts)[0] = c11; // C1111
  (*elasticConsts)[1] = 2.0 * mu * lambda / lambda2mu; // C1122
  (*elasticConsts)[2] = 0; // C1112
  (*elasticConsts)[3] = c11; // C2222
  (*elasticConsts)[4] = 0; // C2212
  (*elasticConsts)[5] = 2.0 * mu; // C1212
} // calcElasticConsts


// End of file 
