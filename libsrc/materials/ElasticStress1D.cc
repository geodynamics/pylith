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

#include "ElasticStress1D.hh" // implementation of object methods

#include "petsc.h" // USES PetscLogFlopsNoCheck

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticStress1D {

      /// Number of entries in stress tensor.
      const int tensorSize = 1;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 1;

      /// Values expected in spatial database
      const int numDBValues = 3;
      const char* namesDBValues[] =
	{"density", "vs", "vp" };

      /// Indices (order) of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;

      /// Parameters
      const int numParameters = 3;
      const int numParamValues[] = { 1, 1, 1 };

      /// Indices of parameters
      const int pidDensity = 0;
      const int pidMu = pidDensity + 1;
      const int pidLambda = pidMu + 1;

    } // _ElasticStress1D
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticStress1D::ElasticStress1D(void) :
  ElasticMaterial(_ElasticStress1D::numParamValues,
		  _ElasticStress1D::numParameters)
{ // constructor
  _dimension = 1;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticStress1D::~ElasticStress1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get names of values expected to be in database of parameters for
const char**
pylith::materials::ElasticStress1D::_dbValues(void) const
{ // _dbValues
  return _ElasticStress1D::namesDBValues;
} // _dbValues

// ----------------------------------------------------------------------
// Get number of values expected to be in database of parameters for
int
pylith::materials::ElasticStress1D::_numDBValues(void) const
{ // _numDBValues
  return _ElasticStress1D::numDBValues;
} // _numDBValues

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticStress1D::_dbToParameters(
					      double* const paramVals,
					      const int numParams,
					      const double_array& dbValues) const
{ // computeParameters
  assert(0 != paramVals);
  assert(_numParamsQuadPt == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticStress1D::numDBValues == numDBValues);

  const double density = dbValues[_ElasticStress1D::didDensity];
  const double vs = dbValues[_ElasticStress1D::didVs];
  const double vp = dbValues[_ElasticStress1D::didVp];
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  paramVals[_ElasticStress1D::pidDensity] = density;
  paramVals[_ElasticStress1D::pidMu] = mu;
  paramVals[_ElasticStress1D::pidLambda] = lambda;

  PetscLogFlopsNoCheck(6);
} // computeParameters

// ----------------------------------------------------------------------
// Get number of entries in stress tensor.
int
pylith::materials::ElasticStress1D::_tensorSize(void) const
{ // _tensorSize
  return _ElasticStress1D::tensorSize;
} // _tensorSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
pylith::materials::ElasticStress1D::_numElasticConsts(void) const
{ // _numElasticConsts
  return _ElasticStress1D::numElasticConsts;
} // _numElasticConsts

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticStress1D::_calcDensity(double* const density,
						 const double* parameters,
						 const int numParams)
{ // _calcDensity
  assert(0 != density);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);

  density[0] = parameters[_ElasticStress1D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticStress1D::_calcStress(double* const stress,
						const int stressSize,
						const double* parameters,
						const int numParams,
						const double* totalStrain,
						const int strainSize)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticStress1D::tensorSize == stressSize);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_ElasticStress1D::tensorSize == strainSize);

  const double density = parameters[_ElasticStress1D::pidDensity];
  const double mu = parameters[_ElasticStress1D::pidMu];
  const double lambda = parameters[_ElasticStress1D::pidLambda];

  const double e11 = totalStrain[0];
  stress[0] = mu * (3.0*lambda+2.0*mu) / (lambda + mu) * e11;

  PetscLogFlopsNoCheck(7);
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::ElasticStress1D::_calcElasticConsts(
				  double* const elasticConsts,
				  const int numElasticConsts,
				  const double* parameters,
				  const int numParams,
				  const double* totalStrain,
				  const int strainSize)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticStress1D::numElasticConsts == numElasticConsts);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_ElasticStress1D::tensorSize == strainSize);
 
  const double density = parameters[_ElasticStress1D::pidDensity];
  const double mu = parameters[_ElasticStress1D::pidMu];
  const double lambda = parameters[_ElasticStress1D::pidLambda];

  elasticConsts[0] = mu * (3.0*lambda+2.0*mu) / (lambda + mu);

  PetscLogFlopsNoCheck(6);
} // _calcElasticConsts


// End of file 
