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

#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
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

      /// Number of physical properties.
      const int numProperties = 3;

      /// Physical properties.
      const Material::PropMetaData properties[] = {
	{ "density", 1, OTHER_FIELD },
	{ "mu", 1, OTHER_FIELD },
	{ "lambda", 1, OTHER_FIELD },
      };
      /// Indices of physical properties
      const int pidDensity = 0;
      const int pidMu = pidDensity + 1;
      const int pidLambda = pidMu + 1;

      /// Values expected in spatial database
      const int numDBValues = 3;
      const char* namesDBValues[] = { "density", "vs", "vp" };      
      
      /// Indices of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;

      /// Initial state values expected in spatial database
      const int numInitialStateDBValues = tensorSize;
      const char* namesInitialStateDBValues[] = { "stress_xx", "stress_yy",
						  "stress_xy" };

    } // _ElasticPlaneStress
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticPlaneStress::ElasticPlaneStress(void) :
  ElasticMaterial(_ElasticPlaneStress::tensorSize,
		  _ElasticPlaneStress::numElasticConsts,
		  _ElasticPlaneStress::namesDBValues,
		  _ElasticPlaneStress::namesInitialStateDBValues,
		  _ElasticPlaneStress::numDBValues,
		  _ElasticPlaneStress::properties,
		  _ElasticPlaneStress::numProperties)
{ // constructor
  _dimension = 2;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticPlaneStress::~ElasticPlaneStress(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticPlaneStress::_dbToProperties(
					  double* propValues,
					  const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_ElasticPlaneStress::numDBValues == numDBValues);

  const double density = dbValues[_ElasticPlaneStress::didDensity];
  const double vs = dbValues[_ElasticPlaneStress::didVs];
  const double vp = dbValues[_ElasticPlaneStress::didVp];
 
  if (density <= 0.0 || vs <= 0.0 || vp <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for physical "
	<< "properties.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if

  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  if (lambda <= 0.0) {
    std::ostringstream msg;
    msg << "Attempted to set Lame's constant lambda to nonpositive value.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[_ElasticPlaneStress::pidDensity] = density;
  propValues[_ElasticPlaneStress::pidMu] = mu;
  propValues[_ElasticPlaneStress::pidLambda] = lambda;

  PetscLogFlops(6);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::ElasticPlaneStress::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticPlaneStress::numProperties);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  values[_ElasticPlaneStress::pidDensity] = 
    _normalizer->nondimensionalize(values[_ElasticPlaneStress::pidDensity],
				   densityScale);
  values[_ElasticPlaneStress::pidMu] = 
    _normalizer->nondimensionalize(values[_ElasticPlaneStress::pidMu],
				   pressureScale);
  values[_ElasticPlaneStress::pidLambda] = 
    _normalizer->nondimensionalize(values[_ElasticPlaneStress::pidLambda],
				   pressureScale);

  PetscLogFlops(3);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::ElasticPlaneStress::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticPlaneStress::numProperties);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  values[_ElasticPlaneStress::pidDensity] = 
    _normalizer->dimensionalize(values[_ElasticPlaneStress::pidDensity],
				   densityScale);
  values[_ElasticPlaneStress::pidMu] = 
    _normalizer->dimensionalize(values[_ElasticPlaneStress::pidMu],
				   pressureScale);
  values[_ElasticPlaneStress::pidLambda] = 
    _normalizer->dimensionalize(values[_ElasticPlaneStress::pidLambda],
				   pressureScale);

  PetscLogFlops(3);
} // _dimProperties

// ----------------------------------------------------------------------
// Nondimensionalize initial state.
void
pylith::materials::ElasticPlaneStress::_nondimInitState(double* const values,
							const int nvalues) const
{ // _nondimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticPlaneStress::numInitialStateDBValues);

  const double pressureScale = _normalizer->pressureScale();
  _normalizer->nondimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _nondimInitState

// ----------------------------------------------------------------------
// Dimensionalize initial state.
void
pylith::materials::ElasticPlaneStress::_dimInitState(double* const values,
						     const int nvalues) const
{ // _dimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticPlaneStress::numInitialStateDBValues);
  
  const double pressureScale = _normalizer->pressureScale();
  _normalizer->dimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _dimInitState

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::ElasticPlaneStress::_calcDensity(double* const density,
						    const double* properties,
						    const int numProperties)
{ // calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  density[0] = properties[_ElasticPlaneStress::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties.
void
pylith::materials::ElasticPlaneStress::_calcStress(double* const stress,
						   const int stressSize,
						   const double* properties,
						   const int numProperties,
						   const double* totalStrain,
						   const int strainSize,
						   const double* initialState,
						   const int initialStateSize,
						   const bool computeStateVars)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticPlaneStress::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStress::tensorSize == strainSize);
  assert(_ElasticPlaneStress::tensorSize == initialStateSize);

  const double density = properties[_ElasticPlaneStress::pidDensity];
  const double mu = properties[_ElasticPlaneStress::pidMu];
  const double lambda = properties[_ElasticPlaneStress::pidLambda];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
  const double lambdamu = lambda + mu;
  
  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e12 = totalStrain[2];

  stress[0] = (2.0*mu2*lambdamu * e11 + mu2*lambda * e22) / lambda2mu + initialState[0];
  stress[1] = (mu2*lambda * e11 + 2.0*mu2*lambdamu * e22) / lambda2mu + initialState[1];
  stress[2] = mu2 * e12 + initialState[2];

  PetscLogFlops(18);
} // _calcStress

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::ElasticPlaneStress::_calcElasticConsts(
						  double* const elasticConsts,
						  const int numElasticConsts,
						  const double* properties,
						  const int numProperties,
						  const double* totalStrain,
						  const int strainSize,
						  const double* initialState,
						  const int initialStateSize)
{ // calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticPlaneStress::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStress::tensorSize == strainSize);
  assert(_ElasticPlaneStress::tensorSize == initialStateSize);

  const double density = properties[_ElasticPlaneStress::pidDensity];
  const double mu = properties[_ElasticPlaneStress::pidMu];
  const double lambda = properties[_ElasticPlaneStress::pidLambda];
  
  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
  const double c11 = 2.0 * mu2 * (lambda + mu) / lambda2mu;

  elasticConsts[0] = c11; // C1111
  elasticConsts[1] = mu2 * lambda / lambda2mu; // C1122
  elasticConsts[2] = 0; // C1112
  elasticConsts[3] = c11; // C2222
  elasticConsts[4] = 0; // C2212
  elasticConsts[5] = mu2; // C1212

  PetscLogFlops(8);
} // calcElasticConsts

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::ElasticPlaneStress::_stableTimeStepImplicit(const double* properties,
				 const int numProperties) const
{ // _stableTimeStepImplicit
  return pylith::PYLITH_MAXDOUBLE;
} // _stableTimeStepImplicit


// End of file 
