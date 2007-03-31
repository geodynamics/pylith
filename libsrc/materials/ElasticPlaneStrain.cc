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

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// _ElasticPlaneStrain is a helper class for ElasticPlaneStrain. We define
// it in this implementation file to insulate other objects from these
// details.
namespace pylith {
  namespace materials {
    class _ElasticPlaneStrain;
  } // materials
} // pylith

class pylith::materials::_ElasticPlaneStrain {
public:
  // Number of entries in stress tensor.
  static const int stressSize;

  // Number of elastic constants (for general 3-D elastic material)
  static const int numElasticConsts;

  // Values expected in spatial database
  static const int numDBValues;
  static const char* namesDBValues[];

  // Indices (order) of database values
  static const int didDensity;
  static const int didVs;
  static const int didVp;

  // Parameters
  static const int numParameters;
  static const char* namesParameters[];

  // Indices (order) of parameters
  static const int pidDensity;
  static const int pidMu;
  static const int pidLambda;
}; // _ElasticPlaneStrain

const int pylith::materials::_ElasticPlaneStrain::stressSize = 3;
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
// Get number of entries in stress tensor.
int
pylith::materials::ElasticPlaneStrain::stressSize(void) const
{ // stressSize
  return _ElasticPlaneStrain::stressSize;
} // stressSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
pylith::materials::ElasticPlaneStrain::numElasticConsts(void) const
{ // numElasticConsts
  return _ElasticPlaneStrain::numElasticConsts;
} // numElasticConsts

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
pylith::materials::ElasticPlaneStrain::_dbToParameters(double* paramVals,
						       const int numParams,
						       const double* dbValues,
						       const int numValues) const
{ // computeParameters
  assert(0 != paramVals);
  assert(_ElasticPlaneStrain::numParameters == numParams);
  assert(0 != dbValues);
  assert(_ElasticPlaneStrain::numDBValues == numValues);

  const double density = dbValues[_ElasticPlaneStrain::didDensity];
  const double vs = dbValues[_ElasticPlaneStrain::didVs];
  const double vp = dbValues[_ElasticPlaneStrain::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  paramVals[_ElasticPlaneStrain::pidDensity] = density;
  paramVals[_ElasticPlaneStrain::pidMu] = mu;
  paramVals[_ElasticPlaneStrain::pidLambda] = lambda;
} // computeParameters

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStrain::_calcDensity(const double* parameters,
						    const int numParameters,
						    const int numLocs)
{ // calcDensity
  assert(0 != _density);
  assert(0 != parameters);

  for (int iLoc=0, index=0; iLoc < numLocs; ++iLoc, index+=numParameters)
    _density[iLoc] = 
      parameters[index+_ElasticPlaneStrain::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticPlaneStrain::_calcStress(const double* parameters,
						   const int numParameters,
						   const double* totalStrain,
						   const int numLocs,
						   const int spaceDim)
{ // _calcStress
  assert(0 != _stress);
  assert(0 != parameters);
  assert(_ElasticPlaneStrain::numParameters == numParameters);
  assert(spaceDim == _dimension);

  for (int iLoc=0, indexP=0, indexS=0;
       iLoc < numLocs; 
       ++iLoc, 
	 indexP+=_ElasticPlaneStrain::numParameters,
	 indexS+=_ElasticPlaneStrain::stressSize) {
    const double density = parameters[indexP+_ElasticPlaneStrain::pidDensity];
    const double mu = parameters[indexP+_ElasticPlaneStrain::pidMu];
    const double lambda = parameters[indexP+_ElasticPlaneStrain::pidLambda];

    const double lambda2mu = lambda + 2.0 * mu;
  
    const double e11 = totalStrain[indexS  ];
    const double e22 = totalStrain[indexS+1];
    const double e12 = totalStrain[indexS+3];

    const double s12 = lambda * (e11 + e22);

    _stress[indexS  ] = s12 + 2.0*mu*e11;
    _stress[indexS+1] = s12 + 2.0*mu*e22;
    _stress[indexS+2] = 2.0 * mu * e12;
  } // for
} // _calcStress

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStrain::_calcElasticConsts(const double* parameters,
							  const int numParameters,
							  const double* totalStrain,
							  const int numLocs,
							  const int spaceDim)
{ // calcElasticConsts
  assert(0 != _elasticConsts);
  assert(0 != parameters);
  assert(_ElasticPlaneStrain::numParameters == numParameters);

  for (int iLoc=0, indexP=0, indexC=0;
       iLoc < numLocs; 
       ++iLoc, 
	 indexP+=_ElasticPlaneStrain::numParameters,
	 indexC+=_ElasticPlaneStrain::numElasticConsts) {
    const double density = parameters[indexP+_ElasticPlaneStrain::pidDensity];
    const double mu = parameters[indexP+_ElasticPlaneStrain::pidMu];
    const double lambda = parameters[indexP+_ElasticPlaneStrain::pidLambda];

    const double lambda2mu = lambda + 2.0*mu;
  
    _elasticConsts[indexC+ 0] = lambda2mu; // C1111
    _elasticConsts[indexC+ 1] = lambda; // C1122
    _elasticConsts[indexC+ 3] = 0; // C1112
    _elasticConsts[indexC+ 6] = lambda2mu; // C2222
    _elasticConsts[indexC+ 8] = 0; // C2212
    _elasticConsts[indexC+15] = 2.0*mu; // C1212
  } // for
} // calcElasticConsts


// End of file 
