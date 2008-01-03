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

#include "petsc.h" // USES PetscLogFlopsNoCheck

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticPlaneStrain {

      /// Number of entries in stress tensor.
      const int tensorSize = 3;

      /// Number of elastic constants (for general 3-D elastic material)
      const int numElasticConsts = 6;

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
      
      /// Indices (order) of parameters
      const int pidDensity = 0;
      const int pidMu = pidDensity + 1;
      const int pidLambda = pidMu + 1;

    } // _ElasticPlaneStrain
  } // materials
} // pylith


// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticPlaneStrain::ElasticPlaneStrain(void) :
  ElasticMaterial(_ElasticPlaneStrain::numParamValues,
		  _ElasticPlaneStrain::numParameters)
{ // constructor
  _dimension = 2;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticPlaneStrain::~ElasticPlaneStrain(void)
{ // destructor
} // destructor

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
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticPlaneStrain::_dbToParameters(
				          double* const paramVals,
				          const int numParams,
                                          const double_array& dbValues) const
{ // computeParameters
  assert(0 != paramVals);
  assert(_numParamsQuadPt == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticPlaneStrain::numDBValues == numDBValues);

  const double density = dbValues[_ElasticPlaneStrain::didDensity];
  const double vs = dbValues[_ElasticPlaneStrain::didVs];
  const double vp = dbValues[_ElasticPlaneStrain::didVp];
 
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
  
  paramVals[_ElasticPlaneStrain::pidDensity] = density;
  paramVals[_ElasticPlaneStrain::pidMu] = mu;
  paramVals[_ElasticPlaneStrain::pidLambda] = lambda;

  PetscLogFlopsNoCheck(6);
} // computeParameters

// ----------------------------------------------------------------------
// Get number of entries in stress tensor.
int
pylith::materials::ElasticPlaneStrain::_tensorSize(void) const
{ // _tensorSize
  return _ElasticPlaneStrain::tensorSize;
} // _tensorSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
pylith::materials::ElasticPlaneStrain::_numElasticConsts(void) const
{ // _numElasticConsts
  return _ElasticPlaneStrain::numElasticConsts;
} // _numElasticConsts

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStrain::_calcDensity(
				  double* const density,
				  const double* parameters,
				  const int numParams)
{ // calcDensity
  assert(0 != density);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);

  density[0] = parameters[_ElasticPlaneStrain::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticPlaneStrain::_calcStress(
				  double* const stress,
				  const int stressSize,
				  const double* parameters,
				  const int numParams,
				  const double* totalStrain,
				  const int strainSize)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticPlaneStrain::tensorSize == stressSize);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStrain::tensorSize == strainSize);

  const double density = parameters[_ElasticPlaneStrain::pidDensity];
  const double mu = parameters[_ElasticPlaneStrain::pidMu];
  const double lambda = parameters[_ElasticPlaneStrain::pidLambda];
  
  const double mu2 = 2.0 * mu;
  
  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e12 = totalStrain[2];

  const double s12 = lambda * (e11 + e22);

  stress[0] = s12 + mu2*e11;
  stress[1] = s12 + mu2*e22;
  stress[2] = mu2 * e12;

  PetscLogFlopsNoCheck(8);
} // _calcStress

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStrain::_calcElasticConsts(
					       double* const elasticConsts,
					       const int numElasticConsts,
					       const double* parameters,
					       const int numParams,
					       const double* totalStrain,
					       const int strainSize)
{ // calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticPlaneStrain::numElasticConsts == numElasticConsts);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStrain::tensorSize == strainSize);

  const double density = parameters[_ElasticPlaneStrain::pidDensity];
  const double mu = parameters[_ElasticPlaneStrain::pidMu];
  const double lambda = parameters[_ElasticPlaneStrain::pidLambda];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
  
  elasticConsts[0] = lambda2mu; // C1111
  elasticConsts[1] = lambda; // C1122
  elasticConsts[2] = 0; // C1112
  elasticConsts[3] = lambda2mu; // C2222
  elasticConsts[4] = 0; // C2212
  elasticConsts[5] = mu2; // C1212

  PetscLogFlopsNoCheck(2);
} // calcElasticConsts


// End of file 
