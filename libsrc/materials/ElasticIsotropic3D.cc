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
}; // _ElasticIsotropic3D

const int pylith::materials::_ElasticIsotropic3D::stressSize = 6;
const int pylith::materials::_ElasticIsotropic3D::numElasticConsts = 21;
const int pylith::materials::_ElasticIsotropic3D::numDBValues = 3;
const char* pylith::materials::_ElasticIsotropic3D::namesDBValues[] =
  {"density", "vs", "vp" };
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
  _dimension = 3;
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
// Get number of entries in stress tensor.
int
pylith::materials::ElasticIsotropic3D::stressSize(void) const
{ // stressSize
  return _ElasticIsotropic3D::stressSize;
} // stressSize

// ----------------------------------------------------------------------
// Get number of entries in elasticity matrix for material.
int
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
						       const double* dbValues,
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
{ // _calcDensity
  assert(0 != _density);
  assert(0 != parameters);

  for (int iLoc=0, index=0; iLoc < numLocs; ++iLoc, index+=numParameters)
    _density[iLoc] = 
      parameters[index+_ElasticIsotropic3D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticIsotropic3D::_calcStress(const double* parameters,
						   const int numParameters,
						   const double* totalStrain,
						   const int numLocs,
						   const int spaceDim)
{ // _calcStress
  assert(0 != _stress);
  assert(0 != parameters);
  assert(_ElasticIsotropic3D::numParameters == numParameters);
  assert(spaceDim == _dimension);

  for (int iLoc=0, indexP=0, indexS=0;
       iLoc < numLocs; 
       ++iLoc, 
	 indexP+=_ElasticIsotropic3D::numParameters,
	 indexS+=_ElasticIsotropic3D::stressSize) {
    const double density = parameters[indexP+_ElasticIsotropic3D::pidDensity];
    const double mu = parameters[indexP+_ElasticIsotropic3D::pidMu];
    const double lambda = parameters[indexP+_ElasticIsotropic3D::pidLambda];

    const double lambda2mu = lambda + 2.0 * mu;
  
    const double e11 = totalStrain[indexS  ];
    const double e22 = totalStrain[indexS+1];
    const double e33 = totalStrain[indexS+2];
    const double e12 = totalStrain[indexS+3];
    const double e23 = totalStrain[indexS+4];
    const double e13 = totalStrain[indexS+5];

    const double s123 = lambda * (e11 + e22 + e33);

    _stress[indexS  ] = s123 + 2.0*mu*e11;
    _stress[indexS+1] = s123 + 2.0*mu*e22;
    _stress[indexS+2] = s123 + 2.0*mu*e33;
    _stress[indexS+3] = 2.0 * mu * e12;
    _stress[indexS+4] = 2.0 * mu * e23;
    _stress[indexS+5] = 2.0 * mu * e13;
  } // for
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::ElasticIsotropic3D::_calcElasticConsts(const double* parameters,
							  const int numParameters,
							  const double* totalStrain,
							  const int numLocs,
							  const int spaceDim)
{ // _calcElasticConsts
  assert(0 != _elasticConsts);
  assert(0 != parameters);
  assert(_ElasticIsotropic3D::numParameters == numParameters);
  assert(spaceDim == _dimension);
 
  for (int iLoc=0, indexP=0, indexC=0;
       iLoc < numLocs; 
       ++iLoc, 
       indexP+=_ElasticIsotropic3D::numParameters,
       indexC+=_ElasticIsotropic3D::numElasticConsts) {
    const double density = parameters[indexP+_ElasticIsotropic3D::pidDensity];
    const double mu = parameters[indexP+_ElasticIsotropic3D::pidMu];
    const double lambda = parameters[indexP+_ElasticIsotropic3D::pidLambda];

    const double lambda2mu = lambda + 2.0 * mu;
   
    _elasticConsts[indexC+ 0] = lambda2mu; // C1111
    _elasticConsts[indexC+ 1] = lambda; // C1122
    _elasticConsts[indexC+ 2] = lambda; // C1133
    _elasticConsts[indexC+ 3] = 0; // C1112
    _elasticConsts[indexC+ 4] = 0; // C1123
    _elasticConsts[indexC+ 5] = 0; // C1113
    _elasticConsts[indexC+ 6] = lambda2mu; // C2222
    _elasticConsts[indexC+ 7] = lambda; // C2233
    _elasticConsts[indexC+ 8] = 0; // C2212
    _elasticConsts[indexC+ 9] = 0; // C2223
    _elasticConsts[indexC+10] = 0; // C2213
    _elasticConsts[indexC+11] = lambda2mu; // C3333
    _elasticConsts[indexC+12] = 0; // C3312
    _elasticConsts[indexC+13] = 0; // C3323
    _elasticConsts[indexC+14] = 0; // C3313
    _elasticConsts[indexC+15] = 2.0 * mu; // C1212
    _elasticConsts[indexC+16] = 0; // C1223
    _elasticConsts[indexC+17] = 0; // C1213
    _elasticConsts[indexC+18] = 2.0 * mu; // C2323
    _elasticConsts[indexC+19] = 0; // C2313
    _elasticConsts[indexC+20] = 2.0 * mu; // C1313
  } // for
} // _calcElasticConsts


// End of file 
