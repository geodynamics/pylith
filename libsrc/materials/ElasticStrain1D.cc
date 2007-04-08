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

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    class _ElasticStrain1D;
  } // materials
} // pylith

/** _ElasticStrain1D is a helper class for ElasticStrain1D. We define
 * it in this implementation file to insulate other objects from these
 * details.
 */
class pylith::materials::_ElasticStrain1D {
public:
  /// Number of entries in stress tensor.
  static const int stressSize;

  /// Number of entries in derivative of elasticity matrix.
  static const int numElasticConsts;

  /// Values expected in spatial database
  static const int numDBValues;
  static const char* namesDBValues[];

  /// Indices (order) of database values
  static const int didDensity;
  static const int didVp;

  /// Parameters
  static const int numParameters;
  static const char* namesParameters[];

  /// Indices (order) of parameters
  static const int pidDensity;
  static const int pidLambda2mu;
}; // _ElasticStrain1D

const int pylith::materials::_ElasticStrain1D::stressSize = 1;
const int pylith::materials::_ElasticStrain1D::numElasticConsts = 1;
const int pylith::materials::_ElasticStrain1D::numDBValues = 2;
const char* pylith::materials::_ElasticStrain1D::namesDBValues[] =
  {"density", "vp" };
const int pylith::materials::_ElasticStrain1D::numParameters = 2;
const char* pylith::materials::_ElasticStrain1D::namesParameters[] =
  {"density", "lambda2mu" };
const int pylith::materials::_ElasticStrain1D::didDensity = 0;
const int pylith::materials::_ElasticStrain1D::didVp = 1;
const int pylith::materials::_ElasticStrain1D::pidDensity = 0;
const int pylith::materials::_ElasticStrain1D::pidLambda2mu = 1;

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
// Get number of entries in stress tensor.
int
pylith::materials::ElasticStrain1D::stressSize(void) const
{ // stressSize
  return _ElasticStrain1D::stressSize;
} // stressSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
pylith::materials::ElasticStrain1D::numElasticConsts(void) const
{ // numElasticConsts
  return _ElasticStrain1D::numElasticConsts;
} // numElasticConsts

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
// Get number of parameters for physical properties.
int
pylith::materials::ElasticStrain1D::_numParameters(void) const
{ // _numParameters
  return _ElasticStrain1D::numParameters;
} // _numParameters

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticStrain1D::_dbToParameters(double* paramVals,
						    const int numParams,
						    const double* dbValues,
						    const int numValues) const
{ // computeParameters
  assert(0 != paramVals);
  assert(_ElasticStrain1D::numParameters == numParams);
  assert(0 != dbValues);
  assert(_ElasticStrain1D::numDBValues == numValues);

  const double density = dbValues[_ElasticStrain1D::didDensity];
  const double vp = dbValues[_ElasticStrain1D::didVp];
  const double lambda2mu = density * vp*vp;

  paramVals[_ElasticStrain1D::pidDensity] = density;
  paramVals[_ElasticStrain1D::pidLambda2mu] = lambda2mu;
} // computeParameters

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticStrain1D::_calcDensity(double* const density,
						 const int size,
						 const double* parameters,
						 const int numParameters)
{ // _calcDensity
  assert(0 != density);
  assert(1 == size);
  assert(0 != parameters);

  *density = parameters[_ElasticStrain1D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticStrain1D::_calcStress(double* const stress,
						const int size,
						const double* parameters,
						const int numParameters,
						const double* totalStrain,
						const int spaceDim)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticStrain1D::stressSize == size);
  assert(0 != parameters);
  assert(_ElasticStrain1D::numParameters == numParameters);
  assert(spaceDim == _dimension);

  const double density = parameters[_ElasticStrain1D::pidDensity];
  const double lambda2mu = parameters[_ElasticStrain1D::pidLambda2mu];

  const double e11 = totalStrain[0];
  stress[0] = lambda2mu * e11;
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::ElasticStrain1D::_calcElasticConsts(double* const elasticConsts,
						       const int size,
						       const double* parameters,
						       const int numParameters,
						       const double* totalStrain,
						       const int spaceDim)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticStrain1D::numElasticConsts == size);
  assert(0 != parameters);
  assert(_ElasticStrain1D::numParameters == numParameters);
  assert(spaceDim == _dimension);
 
  const double density = parameters[_ElasticStrain1D::pidDensity];
  const double lambda2mu = parameters[_ElasticStrain1D::pidLambda2mu];

  elasticConsts[0] = lambda2mu;
} // _calcElasticConsts


// End of file 
