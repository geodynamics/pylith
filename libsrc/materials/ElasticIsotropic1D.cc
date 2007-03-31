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
  // Number of entries in stress tensor.
  static const int stressSize;

  // Number of entries in derivative of elasticity matrix.
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
}; // _ElasticIsotropic1D

const int pylith::materials::_ElasticIsotropic1D::stressSize = 1;
const int pylith::materials::_ElasticIsotropic1D::numElasticConsts = 1;
const int pylith::materials::_ElasticIsotropic1D::numDBValues = 3;
const char* pylith::materials::_ElasticIsotropic1D::namesDBValues[] =
  {"density", "vs", "vp" };
const int pylith::materials::_ElasticIsotropic1D::numParameters = 3;
const char* pylith::materials::_ElasticIsotropic1D::namesParameters[] =
  {"density", "mu", "lambda" };
const int pylith::materials::_ElasticIsotropic1D::didDensity = 0;
const int pylith::materials::_ElasticIsotropic1D::didVs = 1;
const int pylith::materials::_ElasticIsotropic1D::didVp = 2;
const int pylith::materials::_ElasticIsotropic1D::pidDensity = 0;
const int pylith::materials::_ElasticIsotropic1D::pidMu = 1;
const int pylith::materials::_ElasticIsotropic1D::pidLambda = 2;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticIsotropic1D::ElasticIsotropic1D(void)
{ // constructor
  _dimension = 1;
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
// Get number of entries in stress tensor.
int
pylith::materials::ElasticIsotropic1D::stressSize(void) const
{ // stressSize
  return _ElasticIsotropic1D::stressSize;
} // stressSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
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
  const double vs = dbValues[_ElasticIsotropic1D::didVs];
  const double vp = dbValues[_ElasticIsotropic1D::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  paramVals[_ElasticIsotropic1D::pidDensity] = density;
  paramVals[_ElasticIsotropic1D::pidMu] = mu;
  paramVals[_ElasticIsotropic1D::pidLambda] = lambda;
} // computeParameters

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticIsotropic1D::_calcDensity(const double* parameters,
						    const int numParameters,
						    const int numLocs)
{ // _calcDensity
  assert(0 != _density);
  assert(0 != parameters);

  for (int iLoc=0, index=0; iLoc < numLocs; ++iLoc, index+=numParameters)
    _density[iLoc] = 
      parameters[index+_ElasticIsotropic1D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticIsotropic1D::_calcStress(const double* parameters,
						   const int numParameters,
						   const double* totalStrain,
						   const int numLocs,
						   const int spaceDim)
{ // _calcStress
  assert(0 != _stress);
  assert(0 != parameters);
  assert(_ElasticIsotropic1D::numParameters == numParameters);
  assert(spaceDim == _dimension);

  for (int iLoc=0, indexP=0, indexS=0;
       iLoc < numLocs; 
       ++iLoc, 
	 indexP+=_ElasticIsotropic1D::numParameters,
	 indexS+=_ElasticIsotropic1D::stressSize) {
    const double density = parameters[indexP+_ElasticIsotropic1D::pidDensity];
    const double mu = parameters[indexP+_ElasticIsotropic1D::pidMu];
    const double lambda = parameters[indexP+_ElasticIsotropic1D::pidLambda];

    const double e11 = totalStrain[indexS  ];
    _stress[indexS  ] = mu*(3.0*lambda+2.0*mu)/(lambda + mu) * e11;
  } // for
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::ElasticIsotropic1D::_calcElasticConsts(const double* parameters,
							  const int numParameters,
							  const double* totalStrain,
							  const int numLocs,
							  const int spaceDim)
{ // _calcElasticConsts
  assert(0 != _elasticConsts);
  assert(0 != parameters);
  assert(_ElasticIsotropic1D::numParameters == numParameters);
  assert(spaceDim == _dimension);
 
  for (int iLoc=0, indexP=0, indexC=0;
       iLoc < numLocs; 
       ++iLoc, 
       indexP+=_ElasticIsotropic1D::numParameters,
       indexC+=_ElasticIsotropic1D::numElasticConsts) {
    const double density = parameters[indexP+_ElasticIsotropic1D::pidDensity];
    const double mu = parameters[indexP+_ElasticIsotropic1D::pidMu];
    const double lambda = parameters[indexP+_ElasticIsotropic1D::pidLambda];

    _elasticConsts[indexC  ] = mu*(3.0*lambda+2.0*mu)/(lambda + mu);
  } // for
} // _calcElasticConsts


// End of file 
