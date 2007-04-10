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
    class _ElasticStress1D;
  } // materials
} // pylith

/** _ElasticStress1D is a helper class for ElasticStress1D. We define
 * it in this implementation file to insulate other objects from these
 * details.
 */
class pylith::materials::_ElasticStress1D {
public:
  /// Number of entries in stress tensor.
  static const int tensorSize;

  /// Number of entries in derivative of elasticity matrix.
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
}; // _ElasticStress1D

const int pylith::materials::_ElasticStress1D::tensorSize = 1;
const int pylith::materials::_ElasticStress1D::numElasticConsts = 1;
const int pylith::materials::_ElasticStress1D::numDBValues = 3;
const char* pylith::materials::_ElasticStress1D::namesDBValues[] =
  {"density", "vs", "vp" };
const int pylith::materials::_ElasticStress1D::numParameters = 3;
const char* pylith::materials::_ElasticStress1D::namesParameters[] =
  {"density", "mu", "lambda" };
const int pylith::materials::_ElasticStress1D::didDensity = 0;
const int pylith::materials::_ElasticStress1D::didVs = 1;
const int pylith::materials::_ElasticStress1D::didVp = 2;
const int pylith::materials::_ElasticStress1D::pidDensity = 0;
const int pylith::materials::_ElasticStress1D::pidMu = 1;
const int pylith::materials::_ElasticStress1D::pidLambda = 2;

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
// Copy constructor.
pylith::materials::ElasticStress1D::ElasticStress1D(
						const ElasticStress1D& m) :
  ElasticMaterial(m)
{ // copy constructor
} // copy constructor

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
// Get number of parameters for physical properties.
int
pylith::materials::ElasticStress1D::_numParameters(void) const
{ // _numParameters
  return _ElasticStress1D::numParameters;
} // _numParameters

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticStress1D::_dbToParameters(double_array* const paramVals,
					   const double_array& dbValues) const
{ // computeParameters
  assert(0 != paramVals);
  const int numParams = paramVals->size();
  assert(_ElasticStress1D::numParameters == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticStress1D::numDBValues == numDBValues);

  const double density = dbValues[_ElasticStress1D::didDensity];
  const double vs = dbValues[_ElasticStress1D::didVs];
  const double vp = dbValues[_ElasticStress1D::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  (*paramVals)[_ElasticStress1D::pidDensity] = density;
  (*paramVals)[_ElasticStress1D::pidMu] = mu;
  (*paramVals)[_ElasticStress1D::pidLambda] = lambda;
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
					      const double_array& parameters)
{ // _calcDensity
  assert(0 != density);
  assert(1 == density->size());
  assert(_ElasticStress1D::numParameters == parameters.size());

  (*density)[0] = parameters[_ElasticStress1D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticStress1D::_calcStress(double_array* const stress,
					       const double_array& parameters,
					       const double_array& totalStrain)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticStress1D::tensorSize == stress->size());
  assert(_ElasticStress1D::numParameters == parameters.size());
  assert(_ElasticStress1D::tensorSize == totalStrain.size());

  const double density = parameters[_ElasticStress1D::pidDensity];
  const double mu = parameters[_ElasticStress1D::pidMu];
  const double lambda = parameters[_ElasticStress1D::pidLambda];

  const double e11 = totalStrain[0];
  (*stress)[0] = mu * (3.0*lambda+2.0*mu) / (lambda + mu) * e11;
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::ElasticStress1D::_calcElasticConsts(
				       double_array* const elasticConsts,
				       const double_array& parameters,
				       const double_array& totalStrain)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticStress1D::numElasticConsts == elasticConsts->size());
  assert(_ElasticStress1D::numParameters == parameters.size());
  assert(_ElasticStress1D::tensorSize == totalStrain.size());
 
  const double density = parameters[_ElasticStress1D::pidDensity];
  const double mu = parameters[_ElasticStress1D::pidMu];
  const double lambda = parameters[_ElasticStress1D::pidLambda];

  (*elasticConsts)[0] = mu * (3.0*lambda+2.0*mu) / (lambda + mu);
} // _calcElasticConsts


// End of file 
