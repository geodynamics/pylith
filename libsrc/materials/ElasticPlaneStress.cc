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
  static const int stressSize;

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

const int pylith::materials::_ElasticPlaneStress::stressSize = 3;
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
// Get number of entries in stress tensor.
int
pylith::materials::ElasticPlaneStress::stressSize(void) const
{ // stressSize
  return _ElasticPlaneStress::stressSize;
} // stressSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
pylith::materials::ElasticPlaneStress::numElasticConsts(void) const
{ // numElasticConsts
  return _ElasticPlaneStress::numElasticConsts;
} // numElasticConsts

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
pylith::materials::ElasticPlaneStress::_dbToParameters(double* paramVals,
						       const int numParams,
						       const double* dbValues,
						       const int numValues) const
{ // computeParameters
  assert(0 != paramVals);
  assert(_ElasticPlaneStress::numParameters == numParams);
  assert(0 != dbValues);
  assert(_ElasticPlaneStress::numDBValues == numValues);

  const double density = dbValues[_ElasticPlaneStress::didDensity];
  const double vs = dbValues[_ElasticPlaneStress::didVs];
  const double vp = dbValues[_ElasticPlaneStress::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  paramVals[_ElasticPlaneStress::pidDensity] = density;
  paramVals[_ElasticPlaneStress::pidMu] = mu;
  paramVals[_ElasticPlaneStress::pidLambda] = lambda;
} // computeParameters

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStress::_calcDensity(double* const density,
						    const int size,
						    const double* parameters,
						    const int numParameters)
{ // calcDensity
  assert(0 != density);
  assert(1 == size);
  assert(0 != parameters);

  *density = parameters[_ElasticPlaneStress::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticPlaneStress::_calcStress(double* const stress,
						   const int size,
						   const double* parameters,
						   const int numParameters,
						   const double* totalStrain,
						   const int spaceDim)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticPlaneStress::stressSize == size);
  assert(0 != parameters);
  assert(_ElasticPlaneStress::numParameters == numParameters);
  assert(spaceDim == _dimension);

  const double density = parameters[_ElasticPlaneStress::pidDensity];
  const double mu = parameters[_ElasticPlaneStress::pidMu];
  const double lambda = parameters[_ElasticPlaneStress::pidLambda];

  const double lambda2mu = lambda + 2.0 * mu;
  const double lambdamu = lambda + mu;
  
  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e12 = totalStrain[2];

  stress[0] = (4.0*mu*lambdamu * e11 + 2.0*mu*lambda * e22)/lambda2mu;
  stress[1] = (2.0*mu*lambda * e11 + 4.0*mu*lambdamu * e22)/lambda2mu;
  stress[2] = 2.0 * mu * e12;
} // _calcStress

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStress::_calcElasticConsts(double* const elasticConsts,
							  const int size,
							  const double* parameters,
							  const int numParameters,
							  const double* totalStrain,
							  const int spaceDim)
{ // calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticPlaneStress::numElasticConsts == size);
  assert(0 != parameters);
  assert(_ElasticPlaneStress::numParameters == numParameters);

  const double density = parameters[_ElasticPlaneStress::pidDensity];
  const double mu = parameters[_ElasticPlaneStress::pidMu];
  const double lambda = parameters[_ElasticPlaneStress::pidLambda];
  
  const double lambda2mu = lambda + 2.0 * mu;
  const double c11 = 4.0 * mu * (lambda + mu) / lambda2mu;

  elasticConsts[0] = c11; // C1111
  elasticConsts[1] = 2.0 * mu * lambda / lambda2mu; // C1122
  elasticConsts[2] = 0; // C1112
  elasticConsts[3] = c11; // C2222
  elasticConsts[4] = 0; // C2212
  elasticConsts[5] = 2.0 * mu; // C1212
} // calcElasticConsts


// End of file 
