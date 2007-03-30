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

#include "ElasticIsotropic1D.hh" // implementation of object methods

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// _ElasticIsotropic1D is a helper class for ElasticIsotropic1D. We define
// it in this implementation file to insulate other objects from these
// details.
namespace pylith {
  namespace materials {
    class _ElasticIsotropic1D;
  } // materials
} // pylith

class pylith::materials::_ElasticIsotropic1D {
public:
  // Number of elastic constants (for general 3-D elastic material)
  static const int numElasticConsts;

  // Values expected in spatial database
  static const int numDBValues;
  static const char* namesDBValues[];

  // Indices (order) of database values
  static const int didDensity;
  static const int didVp;

  // Parameters
  static const int numParameters;
  static const char* namesParameters[];

  // Indices (order) of parameters
  static const int pidDensity;
  static const int pidLambda2Mu;
}; // _ElasticIsotropic1D

const int pylith::materials::_ElasticIsotropic1D::numElasticConsts = 1;
const int pylith::materials::_ElasticIsotropic1D::numDBValues = 2;
const char* pylith::materials::_ElasticIsotropic1D::namesDBValues[] =
  {"density", "vp" };
const int pylith::materials::_ElasticIsotropic1D::numParameters = 3;
const char* pylith::materials::_ElasticIsotropic1D::namesParameters[] =
  {"density", "lambda+2mu" };
const int pylith::materials::_ElasticIsotropic1D::didDensity = 0;
const int pylith::materials::_ElasticIsotropic1D::didVp = 1;
const int pylith::materials::_ElasticIsotropic1D::pidDensity = 0;
const int pylith::materials::_ElasticIsotropic1D::pidLambda2Mu = 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticIsotropic1D::ElasticIsotropic1D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticIsotropic1D::~ElasticIsotropic1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::ElasticIsotropic1D::ElasticIsotropic1D(
						const ElasticIsotropic1D& m) :
  ElasticMaterial(m)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
const int
pylith::materials::ElasticIsotropic1D::numElasticConsts(void) const
{ // numElasticConsts
  return _ElasticIsotropic1D::numElasticConsts;
} // numElasticConsts

// ----------------------------------------------------------------------
// Get names of values expected to be in database of parameters for
const char**
pylith::materials::ElasticIsotropic1D::_dbValues(void) const
{ // _dbValues
  return _ElasticIsotropic1D::namesDBValues;
} // _dbValues

// ----------------------------------------------------------------------
// Get number of values expected to be in database of parameters for
int
pylith::materials::ElasticIsotropic1D::_numDBValues(void) const
{ // _numDBValues
  return _ElasticIsotropic1D::numDBValues;
} // _numDBValues

// ----------------------------------------------------------------------
// Get names of parameters for physical properties.
const char**
pylith::materials::ElasticIsotropic1D::_parameterNames(void) const
{ // _parameterNames
  return _ElasticIsotropic1D::namesParameters;
} // _parameterNames

// ----------------------------------------------------------------------
// Get number of parameters for physical properties.
int
pylith::materials::ElasticIsotropic1D::_numParameters(void) const
{ // _numParameters
  return _ElasticIsotropic1D::numParameters;
} // _numParameters

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticIsotropic1D::_dbToParameters(double* paramVals,
						       const int numParams,
						       const double* dbValues,
						       const int numValues) const
{ // computeParameters
  assert(0 != paramVals);
  assert(_ElasticIsotropic1D::numParameters == numParams);
  assert(0 != dbValues);
  assert(_ElasticIsotropic1D::numDBValues == numValues);

  const double density = dbValues[_ElasticIsotropic1D::didDensity];
  const double vp = dbValues[_ElasticIsotropic1D::didVp];
 
  const double lambda2mu = density * vp*vp;

  paramVals[_ElasticIsotropic1D::pidDensity] = density;
  paramVals[_ElasticIsotropic1D::pidLambda2Mu] = lambda2mu;
} // computeParameters

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticIsotropic1D::_calcDensity(const double* parameters,
						    const int numParameters,
						    const int numLocs)
{ // calcDensity
  assert(0 != _density);
  assert(0 != parameters);

  for (int iLoc=0, index=0; iLoc < numLocs; ++iLoc, index+=numParameters)
    _density[iLoc] = 
      parameters[index+_ElasticIsotropic1D::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticIsotropic1D::_calcElasticConsts(const double* parameters,
							  const int numParameters,
							  const int numLocs)
{ // calcElasticConsts
  assert(0 != _elasticConsts);
  assert(0 != parameters);
  assert(_ElasticIsotropic1D::numParameters == numParameters);

  for (int iLoc=0, indexP=0, indexC=0;
       iLoc < numLocs; 
       ++iLoc, 
	 indexP+=_ElasticIsotropic1D::numParameters,
	 indexC+=_ElasticIsotropic1D::numElasticConsts) {
    const double density = parameters[indexP+_ElasticIsotropic1D::pidDensity];
    const double lambda2mu = 
      parameters[indexP+_ElasticIsotropic1D::pidLambda2Mu];
  
    _elasticConsts[indexC+ 0] = lambda2mu; // C1111
  } // for
} // calcElasticConsts


// End of file 
