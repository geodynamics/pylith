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

#include "petsc.h" // USES PetscLogFlopsNoCheck

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticPlaneStress {

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

    } // _ElasticPlaneStress
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticPlaneStress::ElasticPlaneStress(void) :
  ElasticMaterial(_ElasticPlaneStress::numParamValues,
		  _ElasticPlaneStress::numParameters)
{ // constructor
  _dimension = 2;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticPlaneStress::~ElasticPlaneStress(void)
{ // destructor
} // destructor

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
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticPlaneStress::_dbToParameters(
				      double* paramVals,
				      const int numParams,
				      const double_array& dbValues) const
{ // _dbToParameters
  assert(0 != paramVals);
  assert(_numParamsQuadPt == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticPlaneStress::numDBValues == numDBValues);

  const double density = dbValues[_ElasticPlaneStress::didDensity];
  const double vs = dbValues[_ElasticPlaneStress::didVs];
  const double vp = dbValues[_ElasticPlaneStress::didVp];
 
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

  paramVals[_ElasticPlaneStress::pidDensity] = density;
  paramVals[_ElasticPlaneStress::pidMu] = mu;
  paramVals[_ElasticPlaneStress::pidLambda] = lambda;

  PetscLogFlopsNoCheck(6);
} // _dbToParameters

// ----------------------------------------------------------------------
// Get number of entries in stress tensor.
int
pylith::materials::ElasticPlaneStress::_tensorSize(void) const
{ // _tensorSize
  return _ElasticPlaneStress::tensorSize;
} // _tensorSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
pylith::materials::ElasticPlaneStress::_numElasticConsts(void) const
{ // _numElasticConsts
  return _ElasticPlaneStress::numElasticConsts;
} // _numElasticConsts

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStress::_calcDensity(double* const density,
						    const double* parameters,
						    const int numParams)
{ // calcDensity
  assert(0 != density);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);

  density[0] = parameters[_ElasticPlaneStress::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticPlaneStress::_calcStress(double* const stress,
						   const int stressSize,
						   const double* parameters,
						   const int numParams,
						   const double* totalStrain,
						   const int strainSize)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticPlaneStress::tensorSize == stressSize);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStress::tensorSize == strainSize);

  const double density = parameters[_ElasticPlaneStress::pidDensity];
  const double mu = parameters[_ElasticPlaneStress::pidMu];
  const double lambda = parameters[_ElasticPlaneStress::pidLambda];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
  const double lambdamu = lambda + mu;
  
  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e12 = totalStrain[2];

  stress[0] = (2.0*mu2*lambdamu * e11 + mu2*lambda * e22) / lambda2mu;
  stress[1] = (mu2*lambda * e11 + 2.0*mu2*lambdamu * e22) / lambda2mu;
  stress[2] = mu2 * e12;

  PetscLogFlopsNoCheck(18);
} // _calcStress

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticPlaneStress::_calcElasticConsts(
						  double* const elasticConsts,
						  const int numElasticConsts,
						  const double* parameters,
						  const int numParams,
						  const double* totalStrain,
						  const int strainSize)
{ // calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticPlaneStress::numElasticConsts == numElasticConsts);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStress::tensorSize == strainSize);

  const double density = parameters[_ElasticPlaneStress::pidDensity];
  const double mu = parameters[_ElasticPlaneStress::pidMu];
  const double lambda = parameters[_ElasticPlaneStress::pidLambda];
  
  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
  const double c11 = 2.0 * mu2 * (lambda + mu) / lambda2mu;

  elasticConsts[0] = c11; // C1111
  elasticConsts[1] = mu2 * lambda / lambda2mu; // C1122
  elasticConsts[2] = 0; // C1112
  elasticConsts[3] = c11; // C2222
  elasticConsts[4] = 0; // C2212
  elasticConsts[5] = mu2; // C1212

  PetscLogFlopsNoCheck(8);
} // calcElasticConsts


// End of file 
