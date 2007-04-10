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
    class _ElasticPlaneStrain;
  } // materials
} // pylith

/** _ElasticPlaneStrain is a helper class for ElasticPlaneStrain. We define
 * it in this implementation file to insulate other objects from these
 * details.
 */
class pylith::materials::_ElasticPlaneStrain {
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
}; // _ElasticPlaneStrain

const int pylith::materials::_ElasticPlaneStrain::tensorSize = 3;
const int pylith::materials::_ElasticPlaneStrain::numElasticConsts = 6;
const int pylith::materials::_ElasticPlaneStrain::numDBValues = 3;
const char* pylith::materials::_ElasticPlaneStrain::namesDBValues[] =
  {"density", "vs", "vp" };
const int pylith::materials::_ElasticPlaneStrain::numParameters = 3;
const char* pylith::materials::_ElasticPlaneStrain::namesParameters[] =
  {"density", "mu", "lambda" };
const int pylith::materials::_ElasticPlaneStrain::didDensity = 0;
const int pylith::materials::_ElasticPlaneStrain::didVs = 1;
const int pylith::materials::_ElasticPlaneStrain::didVp = 2;
const int pylith::materials::_ElasticPlaneStrain::pidDensity = 0;
const int pylith::materials::_ElasticPlaneStrain::pidMu = 1;
const int pylith::materials::_ElasticPlaneStrain::pidLambda = 2;

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
// Copy constructor.
pylith::materials::ElasticPlaneStrain::ElasticPlaneStrain(
						const ElasticPlaneStrain& m) :
  ElasticMaterial(m)
{ // copy constructor
} // copy constructor

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
// Get number of parameters for physical properties.
int
pylith::materials::ElasticPlaneStrain::_numParameters(void) const
{ // _numParameters
  return _ElasticPlaneStrain::numParameters;
} // _numParameters

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticPlaneStrain::_dbToParameters(double_array* const paramVals,
					   const double_array& dbValues) const
{ // computeParameters
  assert(0 != paramVals);
  const int numParams = paramVals->size();
  assert(_ElasticPlaneStrain::numParameters == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticPlaneStrain::numDBValues == numDBValues);

  const double density = dbValues[_ElasticPlaneStrain::didDensity];
  const double vs = dbValues[_ElasticPlaneStrain::didVs];
  const double vp = dbValues[_ElasticPlaneStrain::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  (*paramVals)[_ElasticPlaneStrain::pidDensity] = density;
  (*paramVals)[_ElasticPlaneStrain::pidMu] = mu;
  (*paramVals)[_ElasticPlaneStrain::pidLambda] = lambda;
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
pylith::materials::ElasticPlaneStrain::_calcDensity(double_array* const density,
						    const double_array& parameters)
{ // calcDensity
  assert(0 != density);
  assert(1 == density->size());
  assert(_ElasticPlaneStrain::numParameters == parameters.size());

  (*density)[0] = parameters[_ElasticPlaneStrain::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticPlaneStrain::_calcStress(double_array* const stress,
					       const double_array& parameters,
					       const double_array& totalStrain)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticPlaneStrain::tensorSize == stress->size());
  assert(_ElasticPlaneStrain::numParameters == parameters.size());
  assert(_ElasticPlaneStrain::tensorSize == totalStrain.size());

  const double density = parameters[_ElasticPlaneStrain::pidDensity];
  const double mu = parameters[_ElasticPlaneStrain::pidMu];
  const double lambda = parameters[_ElasticPlaneStrain::pidLambda];
  
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
pylith::materials::ElasticPlaneStrain::_calcElasticConsts(double_array* const elasticConsts,
					     const double_array& parameters,
					     const double_array& totalStrain)
{ // calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticPlaneStrain::numElasticConsts == elasticConsts->size());
  assert(_ElasticPlaneStrain::numParameters == parameters.size());
  assert(_ElasticPlaneStrain::tensorSize == totalStrain.size());

  const double density = parameters[_ElasticPlaneStrain::pidDensity];
  const double mu = parameters[_ElasticPlaneStrain::pidMu];
  const double lambda = parameters[_ElasticPlaneStrain::pidLambda];

  const double lambda2mu = lambda + 2.0*mu;
  
  (*elasticConsts)[0] = lambda2mu; // C1111
  (*elasticConsts)[1] = lambda; // C1122
  (*elasticConsts)[2] = 0; // C1112
  (*elasticConsts)[3] = lambda2mu; // C2222
  (*elasticConsts)[4] = 0; // C2212
  (*elasticConsts)[5] = 2.0*mu; // C1212
} // calcElasticConsts


// End of file 
