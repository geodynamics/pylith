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

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    class _ElasticPlaneStress;
  } // materials
} // pylith

/** _ElasticPlaneStress is a helper class for ElasticPlaneStress. We define
 * it in this implementation file to insulate other objects from these
 * details.
 */
class pylith::materials::_ElasticPlaneStress {
public:
  /// Number of entries in stress tensor.
  static const int tensorSize;

  /// Number of elastic constants (for general 3-D elastic material)
  static const int numElasticConsts;

  /// Values expected in spatial database
  static const int numDBValues;
  static const char* namesDBValues[];

  /// Indices (order) of database values
  static const int didDensity;
  static const int didVs;
  static const int didVp;

  /// Parameters
  static const int numParameters;
  static const char* namesParameters[];

  /// Indices (order) of parameters
  static const int pidDensity;
  static const int pidMu;
  static const int pidLambda;
}; // _ElasticPlaneStress

const int pylith::materials::_ElasticPlaneStress::tensorSize = 3;
const int pylith::materials::_ElasticPlaneStress::numElasticConsts = 6;
const int pylith::materials::_ElasticPlaneStress::numDBValues = 3;
const char* pylith::materials::_ElasticPlaneStress::namesDBValues[] =
  {"density", "vs", "vp" };
const int pylith::materials::_ElasticPlaneStress::numParameters = 3;
const char* pylith::materials::_ElasticPlaneStress::namesParameters[] =
  {"density", "mu", "lambda" };
const int pylith::materials::_ElasticPlaneStress::didDensity = 0;
const int pylith::materials::_ElasticPlaneStress::didVs = 1;
const int pylith::materials::_ElasticPlaneStress::didVp = 2;
const int pylith::materials::_ElasticPlaneStress::pidDensity = 0;
const int pylith::materials::_ElasticPlaneStress::pidMu = 1;
const int pylith::materials::_ElasticPlaneStress::pidLambda = 2;

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
// Copy constructor.
pylith::materials::ElasticPlaneStress::ElasticPlaneStress(
						const ElasticPlaneStress& m) :
  ElasticMaterial(m)
{ // copy constructor
} // copy constructor

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
// Get number of parameters for physical properties.
int
pylith::materials::ElasticPlaneStress::_numParameters(void) const
{ // _numParameters
  return _ElasticPlaneStress::numParameters;
} // _numParameters

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticPlaneStress::_dbToParameters(double_array* paramVals,
					  const double_array& dbValues) const
{ // _dbToParameters
  assert(0 != paramVals);
  const int numParams = paramVals->size();
  assert(_ElasticPlaneStress::numParameters == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticPlaneStress::numDBValues == numDBValues);

  const double density = dbValues[_ElasticPlaneStress::didDensity];
  const double vs = dbValues[_ElasticPlaneStress::didVs];
  const double vp = dbValues[_ElasticPlaneStress::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  (*paramVals)[_ElasticPlaneStress::pidDensity] = density;
  (*paramVals)[_ElasticPlaneStress::pidMu] = mu;
  (*paramVals)[_ElasticPlaneStress::pidLambda] = lambda;
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
						    const double_array& parameters)
{ // calcDensity
  assert(0 != density);
  assert(1 == density->size());
  assert(_ElasticPlaneStress::numParameters == parameters.size());

  (*density)[0] = parameters[_ElasticPlaneStress::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticPlaneStress::_calcStress(double_array* const stress,
					       const double_array& parameters,
					       const double_array& totalStrain)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticPlaneStress::tensorSize == stress->size());
  assert(_ElasticPlaneStress::numParameters == parameters.size());
  assert(_ElasticPlaneStress::tensorSize == totalStrain.size());

  const double density = parameters[_ElasticPlaneStress::pidDensity];
  const double mu = parameters[_ElasticPlaneStress::pidMu];
  const double lambda = parameters[_ElasticPlaneStress::pidLambda];

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
pylith::materials::ElasticPlaneStress::_calcElasticConsts(double_array* const elasticConsts,
					     const double_array& parameters,
					     const double_array& totalStrain)
{ // calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticPlaneStress::numElasticConsts == elasticConsts->size());
  assert(_ElasticPlaneStress::numParameters == parameters.size());
  assert(_ElasticPlaneStress::tensorSize == totalStrain.size());

  const double density = parameters[_ElasticPlaneStress::pidDensity];
  const double mu = parameters[_ElasticPlaneStress::pidMu];
  const double lambda = parameters[_ElasticPlaneStress::pidLambda];
  
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
