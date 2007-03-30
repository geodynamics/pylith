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
// _ElasticPlaneStress is a helper class for ElasticPlaneStress. We define
// it in this implementation file to insulate other objects from these
// details.
namespace pylith {
  namespace materials {
    class _ElasticPlaneStress;
  } // materials
} // pylith

class pylith::materials::_ElasticPlaneStress {
public:
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
}; // _ElasticPlaneStress

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
// Get number of elastic constants for material.
const int
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
pylith::materials::ElasticPlaneStress::_calcDensity(const double* parameters,
						    const int numParameters,
						    const int numLocs)
{ // calcDensity
  assert(0 != _density);
  assert(0 != parameters);

  for (int iLoc=0, index=0; iLoc < numLocs; ++iLoc, index+=numParameters)
    _density[iLoc] = 
      parameters[index+_ElasticPlaneStress::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStress::_calcElasticConsts(const double* parameters,
							  const int numParameters,
							  const int numLocs)
{ // calcElasticConsts
  assert(0 != _elasticConsts);
  assert(0 != parameters);
  assert(_ElasticPlaneStress::numParameters == numParameters);

  for (int iLoc=0, indexP=0, indexC=0;
       iLoc < numLocs; 
       ++iLoc, 
	 indexP+=_ElasticPlaneStress::numParameters,
	 indexC+=_ElasticPlaneStress::numElasticConsts) {
    const double density = parameters[indexP+_ElasticPlaneStress::pidDensity];
    const double mu = parameters[indexP+_ElasticPlaneStress::pidMu];
    const double lambda = parameters[indexP+_ElasticPlaneStress::pidLambda];
  
    _elasticConsts[indexC+ 0] = lambda + 2.0*mu; // C1111
    _elasticConsts[indexC+ 1] = lambda; // C1122
    _elasticConsts[indexC+ 2] = 0; // C1112
    _elasticConsts[indexC+ 3] = lambda + 2.0*mu; // C2222
    _elasticConsts[indexC+ 4] = 0; // C2212
    _elasticConsts[indexC+15] = mu; // C1212
  } // for
} // calcElasticConsts


// End of file 
