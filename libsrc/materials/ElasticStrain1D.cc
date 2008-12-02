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

#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticStrain1D {

      /// Number of entries in stress tensor.
      const int tensorSize = 1;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 1;

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

      /// Initial state values expected in spatial database.
      const int numInitialStateDBValues = tensorSize;
      const char* namesInitialStateDBValues[] = { "stress_xx" };
      
    } // _ElasticStrain1D
  } // materials
} // pylith


// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticStrain1D::ElasticStrain1D(void) :
  ElasticMaterial(_ElasticStrain1D::tensorSize,
		  _ElasticStrain1D::numElasticConsts,
		  _ElasticStrain1D::namesDBValues,
		  _ElasticStrain1D::namesInitialStateDBValues,
		  _ElasticStrain1D::numDBValues,
		  _ElasticStrain1D::properties,
		  _ElasticStrain1D::numProperties)
{ // constructor
  _dimension = 1;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticStrain1D::~ElasticStrain1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticStrain1D::_dbToProperties(
					    double* const propValues,
					    const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_ElasticStrain1D::numDBValues == numDBValues);

  const double density = dbValues[_ElasticStrain1D::didDensity];
  const double vp = dbValues[_ElasticStrain1D::didVp];
  const double vs = dbValues[_ElasticStrain1D::didVs];

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

  propValues[_ElasticStrain1D::pidDensity] = density;
  propValues[_ElasticStrain1D::pidMu] = mu;
  propValues[_ElasticStrain1D::pidLambda] = lambda;

  PetscLogFlops(6);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::ElasticStrain1D::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticStrain1D::numProperties);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  values[_ElasticStrain1D::pidDensity] = 
    _normalizer->nondimensionalize(values[_ElasticStrain1D::pidDensity],
				   densityScale);
  values[_ElasticStrain1D::pidMu] = 
    _normalizer->nondimensionalize(values[_ElasticStrain1D::pidMu],
				   pressureScale);
  values[_ElasticStrain1D::pidLambda] = 
    _normalizer->nondimensionalize(values[_ElasticStrain1D::pidLambda],
				   pressureScale);

  PetscLogFlops(3);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::ElasticStrain1D::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticStrain1D::numProperties);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  values[_ElasticStrain1D::pidDensity] = 
    _normalizer->dimensionalize(values[_ElasticStrain1D::pidDensity],
				   densityScale);
  values[_ElasticStrain1D::pidMu] = 
    _normalizer->dimensionalize(values[_ElasticStrain1D::pidMu],
				   pressureScale);
  values[_ElasticStrain1D::pidLambda] = 
    _normalizer->dimensionalize(values[_ElasticStrain1D::pidLambda],
				   pressureScale);

  PetscLogFlops(3);
} // _dimProperties

// ----------------------------------------------------------------------
// Nondimensionalize initial state.
void
pylith::materials::ElasticStrain1D::_nondimInitState(double* const values,
							const int nvalues) const
{ // _nondimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticStrain1D::numInitialStateDBValues);

  const double pressureScale = _normalizer->pressureScale();
  _normalizer->nondimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _nondimInitState

// ----------------------------------------------------------------------
// Dimensionalize initial state.
void
pylith::materials::ElasticStrain1D::_dimInitState(double* const values,
						     const int nvalues) const
{ // _dimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticStrain1D::numInitialStateDBValues);
  
  const double pressureScale = _normalizer->pressureScale();
  _normalizer->dimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _dimInitState

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::ElasticStrain1D::_calcDensity(double* const density,
						 const double* properties,
						 const int numProperties)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  density[0] = properties[_ElasticStrain1D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties.
void
pylith::materials::ElasticStrain1D::_calcStress(
				   double* const stress,
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
  assert(_ElasticStrain1D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_ElasticStrain1D::tensorSize == strainSize);
  assert(_ElasticStrain1D::tensorSize == initialStateSize);

  const double density = properties[_ElasticStrain1D::pidDensity];
  const double lambda = properties[_ElasticStrain1D::pidLambda];
  const double mu = properties[_ElasticStrain1D::pidMu];

  const double e11 = totalStrain[0];
  stress[0] = (lambda + 2.0*mu) * e11 + initialState[0];

  PetscLogFlops(4);
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::ElasticStrain1D::_calcElasticConsts(
				   double* const elasticConsts,
				   const int numElasticConsts,
				   const double* properties,
				   const int numProperties,
				   const double* totalStrain,
				   const int strainSize,
				   const double* initialState,
				   const int initialStateSize)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticStrain1D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_ElasticStrain1D::tensorSize == strainSize);
  assert(_ElasticStrain1D::tensorSize == initialStateSize);
 
  const double density = properties[_ElasticStrain1D::pidDensity];
  const double lambda = properties[_ElasticStrain1D::pidLambda];
  const double mu = properties[_ElasticStrain1D::pidMu];

  elasticConsts[0] = lambda + 2.0*mu;

  PetscLogFlops(2);
} // _calcElasticConsts

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::ElasticStrain1D::_stableTimeStepImplicit(const double* properties,
				 const int numProperties) const
{ // _stableTimeStepImplicit
  return pylith::PYLITH_MAXDOUBLE;
} // _stableTimeStepImplicit


// End of file 
