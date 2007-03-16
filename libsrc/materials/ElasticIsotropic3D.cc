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

#include "ElasticIsotropic3D.hh" // implementation of object methods

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// _ElasticIsotropic is a helper class for ElasticIsotropic. We define
// it in this implementation file to insulate other objects from these
// details.
namespace pylith {
  namespace materials {
    class _ElasticIsotropic3D;
  } // materials
} // pylith

class pylith::materials::_ElasticIsotropic3D {
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
}; // _ElasticIsotropic3D

const int pylith::materials::_ElasticIsotropic3D::numElasticConsts = 21;
const int pylith::materials::_ElasticIsotropic3D::numDBValues = 3;
const char* pylith::materials::_ElasticIsotropic3D::namesDBValues[] =
  {"density", "vp", "vs" };
const int pylith::materials::_ElasticIsotropic3D::numParameters = 3;
const char* pylith::materials::_ElasticIsotropic3D::namesParameters[] =
  {"density", "mu", "lambda" };
const int pylith::materials::_ElasticIsotropic3D::didDensity = 0;
const int pylith::materials::_ElasticIsotropic3D::didVs = 1;
const int pylith::materials::_ElasticIsotropic3D::didVp = 2;
const int pylith::materials::_ElasticIsotropic3D::pidDensity = 0;
const int pylith::materials::_ElasticIsotropic3D::pidMu = 1;
const int pylith::materials::_ElasticIsotropic3D::pidLambda = 2;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticIsotropic3D::ElasticIsotropic3D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticIsotropic3D::~ElasticIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::ElasticIsotropic3D::ElasticIsotropic3D(
						const ElasticIsotropic3D& m) :
  ElasticMaterial(m)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
const int
pylith::materials::ElasticIsotropic3D::numElasticConsts(void) const
{ // numElasticConsts
  return _ElasticIsotropic3D::numElasticConsts;
} // numElasticConsts

// ----------------------------------------------------------------------
// Get names of values expected to be in database of parameters for
const char**
pylith::materials::ElasticIsotropic3D::_dbValues(void) const
{ // _dbValues
  return _ElasticIsotropic3D::namesDBValues;
} // _dbValues

// ----------------------------------------------------------------------
// Get number of values expected to be in database of parameters for
int
pylith::materials::ElasticIsotropic3D::_numDBValues(void) const
{ // _numDBValues
  return _ElasticIsotropic3D::numDBValues;
} // _numDBValues

// ----------------------------------------------------------------------
// Get names of parameters for physical properties.
const char**
pylith::materials::ElasticIsotropic3D::_parameterNames(void) const
{ // _parameterNames
  return _ElasticIsotropic3D::namesParameters;
} // _parameterNames

// ----------------------------------------------------------------------
// Get number of parameters for physical properties.
int
pylith::materials::ElasticIsotropic3D::_numParameters(void) const
{ // _numParameters
  return _ElasticIsotropic3D::numParameters;
} // _numParameters

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticIsotropic3D::_dbToParameters(double* paramVals,
						       const int numParams,
						       double* dbValues,
						       const int numValues) const
{ // computeParameters
  assert(0 != paramVals);
  assert(_ElasticIsotropic3D::numParameters == numParams);
  assert(0 != dbValues);
  assert(_ElasticIsotropic3D::numDBValues == numValues);

  const double density = dbValues[_ElasticIsotropic3D::didDensity];
  const double vs = dbValues[_ElasticIsotropic3D::didVs];
  const double vp = dbValues[_ElasticIsotropic3D::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  paramVals[_ElasticIsotropic3D::pidDensity] = density;
  paramVals[_ElasticIsotropic3D::pidMu] = mu;
  paramVals[_ElasticIsotropic3D::pidLambda] = lambda;
} // computeParameters

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticIsotropic3D::_calcDensity(const double* parameters,
						    const int numParameters,
						    const int numLocs)
{ // calcDensity
  assert(0 != _density);
  assert(0 != parameters);

  for (int iLoc=0, index=0; iLoc < numLocs; ++iLoc, index+=numParameters)
    _density[iLoc] = 
      parameters[index+_ElasticIsotropic3D::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticIsotropic3D::_calcElasticConsts(const double* parameters,
							  const int numParameters,
							  const int numLocs)
{ // calcElasticConsts
  assert(0 != _elasticConsts);
  assert(0 != parameters);
  assert(_ElasticIsotropic3D::numParameters == numParameters);

  for (int iLoc=0, indexP=0, indexC=0;
       iLoc < numLocs; 
       ++iLoc, 
	 indexP+=_ElasticIsotropic3D::numParameters,
	 indexC+=_ElasticIsotropic3D::numElasticConsts) {
    const double density = parameters[indexP+_ElasticIsotropic3D::pidDensity];
    const double mu = parameters[indexP+_ElasticIsotropic3D::pidMu];
    const double lambda = parameters[indexP+_ElasticIsotropic3D::pidLambda];
  
    _elasticConsts[indexC+ 0] = lambda + 2.0*mu; // C1111
    _elasticConsts[indexC+ 1] = lambda; // C1122
    _elasticConsts[indexC+ 2] = lambda; // C1133
    _elasticConsts[indexC+ 3] = 0; // C1112
    _elasticConsts[indexC+ 4] = 0; // C1123
    _elasticConsts[indexC+ 5] = 0; // C1113
    _elasticConsts[indexC+ 6] = lambda + 2.0*mu; // C2222
    _elasticConsts[indexC+ 7] = lambda; // C2233
    _elasticConsts[indexC+ 8] = 0; // C2212
    _elasticConsts[indexC+ 9] = 0; // C2223
    _elasticConsts[indexC+10] = 0; // C2212
    _elasticConsts[indexC+11] = lambda + 2.0*mu; // C3333
    _elasticConsts[indexC+12] = 0; // C3312
    _elasticConsts[indexC+13] = 0; // C3323
    _elasticConsts[indexC+14] = 0; // C3312
    _elasticConsts[indexC+15] = mu; // C1212
    _elasticConsts[indexC+16] = 0; // C1223
    _elasticConsts[indexC+17] = 0; // C1213
    _elasticConsts[indexC+18] = mu; // C2323
    _elasticConsts[indexC+19] = 0; // C2313
    _elasticConsts[indexC+20] = mu; // C1313
  } // for
} // calcElasticConsts


// End of file 
