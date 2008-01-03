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

#include "ElasticStrain1D.hh" // implementation of object methods

#include "petsc.h" // USES PetscLogFlopsNoCheck

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticStrain1D {

      /// Number of entries in stress tensor.
      const int tensorSize = 1;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 1;

      /// Values expected in spatial database
      const int numDBValues = 2;
      const char* namesDBValues[] = { "density", "vp" };      
      
      /// Indices of database values
      const int didDensity = 0;
      const int didVp = 1;
      
      /// Parameters
      const int numParameters = 2;
      const int numParamValues[] = { 1, 1 };
      
      /// Indices of parameters
      const int pidDensity = 0;
      const int pidLambda2mu = pidDensity + 1;
    } // _ElasticStrain1D
  } // materials
} // pylith


// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticStrain1D::ElasticStrain1D(void) :
  ElasticMaterial(_ElasticStrain1D::numParamValues,
		  _ElasticStrain1D::numParameters)
{ // constructor
  _dimension = 1;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticStrain1D::~ElasticStrain1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get names of values expected to be in database of parameters for
const char**
pylith::materials::ElasticStrain1D::_dbValues(void) const
{ // _dbValues
  return _ElasticStrain1D::namesDBValues;
} // _dbValues

// ----------------------------------------------------------------------
// Get number of values expected to be in database of parameters for
int
pylith::materials::ElasticStrain1D::_numDBValues(void) const
{ // _numDBValues
  return _ElasticStrain1D::numDBValues;
} // _numDBValues

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticStrain1D::_dbToParameters(
					    double* const paramVals,
					    const int numParams,
					    const double_array& dbValues) const
{ // _dbToParameters
  assert(0 != paramVals);
  assert(_numParamsQuadPt == numParams);
  const int numDBValues = dbValues.size();
  assert(_ElasticStrain1D::numDBValues == numDBValues);

  const double density = dbValues[_ElasticStrain1D::didDensity];
  const double vp = dbValues[_ElasticStrain1D::didVp];
  const double lambda2mu = density * vp*vp;

  paramVals[_ElasticStrain1D::pidDensity] = density;
  paramVals[_ElasticStrain1D::pidLambda2mu] = lambda2mu;

  PetscLogFlopsNoCheck(2);
} // _dbToParameters

// ----------------------------------------------------------------------
// Get number of entries in stress tensor.
int
pylith::materials::ElasticStrain1D::_tensorSize(void) const
{ // _tensorSize
  return _ElasticStrain1D::tensorSize;
} // _tensorSize

// ----------------------------------------------------------------------
// Get number of elastic constants for material.
int
pylith::materials::ElasticStrain1D::_numElasticConsts(void) const
{ // _numElasticConsts
  return _ElasticStrain1D::numElasticConsts;
} // _numElasticConsts

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::ElasticStrain1D::_calcDensity(double* const density,
						 const double* parameters,
						 const int numParams)
{ // _calcDensity
  assert(0 != density);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);

  density[0] = parameters[_ElasticStrain1D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::ElasticStrain1D::_calcStress(
				   double* const stress,
				   const int stressSize,
				   const double* parameters,
				   const int numParams,
				   const double* totalStrain,
				   const int strainSize)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticStrain1D::tensorSize == stressSize);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_ElasticStrain1D::tensorSize == strainSize);

  const double density = parameters[_ElasticStrain1D::pidDensity];
  const double lambda2mu = parameters[_ElasticStrain1D::pidLambda2mu];

  const double e11 = totalStrain[0];
  stress[0] = lambda2mu * e11;

  PetscLogFlopsNoCheck(1);
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::ElasticStrain1D::_calcElasticConsts(
				   double* const elasticConsts,
				   const int numElasticConsts,
				   const double* parameters,
				   const int numParams,
				   const double* totalStrain,
				   const int strainSize)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticStrain1D::numElasticConsts == numElasticConsts);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_ElasticStrain1D::tensorSize == strainSize);
 
  const double density = parameters[_ElasticStrain1D::pidDensity];
  const double lambda2mu = parameters[_ElasticStrain1D::pidLambda2mu];

  elasticConsts[0] = lambda2mu;
} // _calcElasticConsts


// End of file 
