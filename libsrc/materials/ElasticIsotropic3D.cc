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
namespace pylith {
  namespace materials {
    class _ElasticIsotropic3D;
  } // materials
} // pylith

class pylith::materials::_ElasticIsotropic3D {
public:
  static const int numDBValues;
  static const char* namesDBValues[];
  static const int numParameters;
  static const char* namesParameters[];
  static const int didDensity;
  static const int didVs;
  static const int didVp;
  static const int pidDensity;
  static const int pidMu;
  static const int pidLambda;
}; // _ElasticIsotropic3D

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
  ElasticMaterial3D(m)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticIsotropic3D::dbToParams(double* paramVals,
						  const int numParams,
						  double* dbValues,
						  const int numValues) const
{ // computeParameters
  assert(0 != paramVals);
  assert(_NUMPARAMETERS == numParams);
  assert(0 != dbValues);
  assert(_NUMDBVALUES == numValues);

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
pylith::materials::ElasticIsotropic3D::calcDensity(double* density,
						   const double* parameters,
						   const int numParameters)
{ // calcDensity
  assert(0 != density);
  *density = parameters[_ElasticIsotropic3D::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticIsotropic3D::calcElasticConsts(double* elasticConsts,
							 const int numConsts,
							 const double* parameters,
							 const int numParameters)
{ // calcElasticConsts
  assert(0 != elasticConsts);
  assert(NUMELASTCONSTS == numConsts);
  assert(0 != parameters);
  assert(_ElasticIsotropic3D::numParameters == numParameters);

  const double density = parameters[_ElasticIsotropic3D::pidDensity];
  const double mu = parameters[_ElasticIsotropic3D::pidMu];
  const double lambda = parameters[_ElasticIsotropic3D::pidLambda];
  
  elasticConsts[ 0] = lambda + 2.0*mu; // C1111
  elasticConsts[ 1] = lambda; // C1122
  elasticConsts[ 2] = lambda; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = 0; // C1123
  elasticConsts[ 5] = 0; // C1113
  elasticConsts[ 6] = lambda + 2.0*mu; // C2222
  elasticConsts[ 7] = lambda; // C2233
  elasticConsts[ 8] = 0; // C2212
  elasticConsts[ 9] = 0; // C2223
  elasticConsts[10] = 0; // C2212
  elasticConsts[11] = lambda + 2.0*mu; // C3333
  elasticConsts[12] = 0; // C3312
  elasticConsts[13] = 0; // C3323
  elasticConsts[14] = 0; // C3312
  elasticConsts[15] = mu; // C1212
  elasticConsts[16] = 0; // C1223
  elasticConsts[17] = 0; // C1213
  elasticConsts[18] = mu; // C2323
  elasticConsts[19] = 0; // C2313
  elasticConsts[20] = mu; // C1313
} // calcElasticConsts


// End of file 
