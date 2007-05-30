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

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticIsotropic3D {

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 6;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 21;

      /// Values expected in spatial database
      const int numDBValues = 3;
      const char* namesDBValues[] = { "density", "vs", "vp" };

      /// Indices (order) of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;

      /// Parameters
      const int numParameters = 3;
      const int numParamValues[] = { 1, 1, 1 };
      const char* namesParameters[] = { "density", "mu", "lambda" };

      /// Indices (order) of parameters
      const int pidDensity = 0;
      const int pidMu = 1;
      const int pidLambda = 2;

    } // _ElasticIsotropic3D
  } // materials
} // pylith


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
// Get number of values for each parameter for physical properties.
void
pylith::materials::ElasticIsotropic3D::_numParamValues(int_array* numValues) const
{ // _numParamValues
  assert(0 != numValues);

  const int numParams = _ElasticIsotropic3D::numParameters;
  numValues->resize(numParams);
  for (int i=0; i < numParams; ++i)
    (*numValues)[i] = _ElasticIsotropic3D::numParamValues[i];
} // _numParamValues

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticIsotropic3D::_dbToParameters(std::vector<double_array>* paramVals,
					  const double_array& dbValues) const
{ // _dbToParameters
  assert(0 != paramVals);
  const int numParams = paramVals->size();
  assert(_ElasticIsotropic3D::numParameters == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticIsotropic3D::numDBValues == numDBValues);
  for (int i=0; i < numParams; ++i)
    assert(_ElasticIsotropic3D::numParamValues[i] == (*paramVals)[i].size());

  const double density = dbValues[_ElasticIsotropic3D::didDensity];
  const double vs = dbValues[_ElasticIsotropic3D::didVs];
  const double vp = dbValues[_ElasticIsotropic3D::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;
  if (lambda < 0.0) {
    std::ostringstream msg;
    msg << "Attempted to set Lame's constant lambda to negative value.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if

  (*paramVals)[_ElasticIsotropic3D::pidDensity][0] = density;
  (*paramVals)[_ElasticIsotropic3D::pidMu][0] = mu;
  (*paramVals)[_ElasticIsotropic3D::pidLambda][0] = lambda;
} // _dbToParameters

// ----------------------------------------------------------------------
// Get number of entries in stress tensor.
int
pylith::materials::ElasticIsotropic3D::_tensorSize(void) const
{ // _tensorSize
  return _ElasticIsotropic3D::tensorSize;
} // _tensorSize

// ----------------------------------------------------------------------
// Get number of entries in elasticity matrix for material.
int
pylith::materials::ElasticIsotropic3D::_numElasticConsts(void) const
{ // numElasticConsts
  return _ElasticIsotropic3D::numElasticConsts;
} // numElasticConsts

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticIsotropic3D::_calcDensity(
				  double_array* const density,
				  const std::vector<double_array>& parameters)
{ // _calcDensity
  assert(0 != density);
  assert(1 == density->size());
  assert(_ElasticIsotropic3D::numParameters == parameters.size());
  assert(1 == parameters[_ElasticIsotropic3D::pidDensity].size());

  (*density)[0] = parameters[_ElasticIsotropic3D::pidDensity][0];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticIsotropic3D::_calcStress(
				  double_array* const stress,
				  const std::vector<double_array>& parameters,
				  const double_array& totalStrain)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticIsotropic3D::tensorSize == stress->size());
  assert(_ElasticIsotropic3D::numParameters == parameters.size());
  assert(_ElasticIsotropic3D::tensorSize == totalStrain.size());
  assert(1 == parameters[_ElasticIsotropic3D::pidDensity].size());
  assert(1 == parameters[_ElasticIsotropic3D::pidMu].size());
  assert(1 == parameters[_ElasticIsotropic3D::pidLambda].size());

  const double density = parameters[_ElasticIsotropic3D::pidDensity][0];
  const double mu = parameters[_ElasticIsotropic3D::pidMu][0];
  const double lambda = parameters[_ElasticIsotropic3D::pidLambda][0];

  const double lambda2mu = lambda + 2.0 * mu;
  
  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double s123 = lambda * (e11 + e22 + e33);

  (*stress)[0] = s123 + 2.0*mu*e11;
  (*stress)[1] = s123 + 2.0*mu*e22;
  (*stress)[2] = s123 + 2.0*mu*e33;
  (*stress)[3] = 2.0 * mu * e12;
  (*stress)[4] = 2.0 * mu * e23;
  (*stress)[5] = 2.0 * mu * e13;
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::ElasticIsotropic3D::_calcElasticConsts(
				  double_array* const elasticConsts,
				  const std::vector<double_array>& parameters,
				  const double_array& totalStrain)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticIsotropic3D::numElasticConsts == elasticConsts->size());
  assert(_ElasticIsotropic3D::numParameters == parameters.size());
  assert(_ElasticIsotropic3D::tensorSize == totalStrain.size());
  assert(1 == parameters[_ElasticIsotropic3D::pidDensity].size());
  assert(1 == parameters[_ElasticIsotropic3D::pidMu].size());
  assert(1 == parameters[_ElasticIsotropic3D::pidLambda].size());
 
  const double density = parameters[_ElasticIsotropic3D::pidDensity][0];
  const double mu = parameters[_ElasticIsotropic3D::pidMu][0];
  const double lambda = parameters[_ElasticIsotropic3D::pidLambda][0];

  const double lambda2mu = lambda + 2.0 * mu;
   
  (*elasticConsts)[ 0] = lambda2mu; // C1111
  (*elasticConsts)[ 1] = lambda; // C1122
  (*elasticConsts)[ 2] = lambda; // C1133
  (*elasticConsts)[ 3] = 0; // C1112
  (*elasticConsts)[ 4] = 0; // C1123
  (*elasticConsts)[ 5] = 0; // C1113
  (*elasticConsts)[ 6] = lambda2mu; // C2222
  (*elasticConsts)[ 7] = lambda; // C2233
  (*elasticConsts)[ 8] = 0; // C2212
  (*elasticConsts)[ 9] = 0; // C2223
  (*elasticConsts)[10] = 0; // C2213
  (*elasticConsts)[11] = lambda2mu; // C3333
  (*elasticConsts)[12] = 0; // C3312
  (*elasticConsts)[13] = 0; // C3323
  (*elasticConsts)[14] = 0; // C3313
  (*elasticConsts)[15] = 2.0 * mu; // C1212
  (*elasticConsts)[16] = 0; // C1223
  (*elasticConsts)[17] = 0; // C1213
  (*elasticConsts)[18] = 2.0 * mu; // C2323
  (*elasticConsts)[19] = 0; // C2313
  (*elasticConsts)[20] = 2.0 * mu; // C1313
} // _calcElasticConsts


// End of file 
