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
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticIsotropic3D {

      /// Number of entries in stress tensor.
      const int tensorSize = 6;

      /// Number of elastic constants (for general 3-D elastic material)
      const int numElasticConsts = 21;

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
						  "stress_zz", "stress_xy",
						  "stress_yz", "stress_xz" };

    } // _ElasticIsotropic3D
  } // materials
} // pylith


// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticIsotropic3D::ElasticIsotropic3D(void) :
  ElasticMaterial(_ElasticIsotropic3D::tensorSize,
		  _ElasticIsotropic3D::numElasticConsts,
		  _ElasticIsotropic3D::namesDBValues,
		  _ElasticIsotropic3D::namesInitialStateDBValues,
		  _ElasticIsotropic3D::numDBValues,
		  _ElasticIsotropic3D::properties,
		  _ElasticIsotropic3D::numProperties)
{ // constructor
  _dimension = 3;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticIsotropic3D::~ElasticIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::ElasticIsotropic3D::_dbToProperties(
		   double* const propValues,
		   const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_ElasticIsotropic3D::numDBValues == numDBValues);

  const double density = dbValues[_ElasticIsotropic3D::didDensity];
  const double vs = dbValues[_ElasticIsotropic3D::didVs];
  const double vp = dbValues[_ElasticIsotropic3D::didVp];
 
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

  propValues[_ElasticIsotropic3D::pidDensity] = density;
  propValues[_ElasticIsotropic3D::pidMu] = mu;
  propValues[_ElasticIsotropic3D::pidLambda] = lambda;

  PetscLogFlops(6);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::ElasticIsotropic3D::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticIsotropic3D::numProperties);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  values[_ElasticIsotropic3D::pidDensity] = 
    _normalizer->nondimensionalize(values[_ElasticIsotropic3D::pidDensity],
				   densityScale);
  values[_ElasticIsotropic3D::pidMu] = 
    _normalizer->nondimensionalize(values[_ElasticIsotropic3D::pidMu],
				   pressureScale);
  values[_ElasticIsotropic3D::pidLambda] = 
    _normalizer->nondimensionalize(values[_ElasticIsotropic3D::pidLambda],
				   pressureScale);

  PetscLogFlops(3);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::ElasticIsotropic3D::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticIsotropic3D::numProperties);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  values[_ElasticIsotropic3D::pidDensity] = 
    _normalizer->dimensionalize(values[_ElasticIsotropic3D::pidDensity],
				   densityScale);
  values[_ElasticIsotropic3D::pidMu] = 
    _normalizer->dimensionalize(values[_ElasticIsotropic3D::pidMu],
				   pressureScale);
  values[_ElasticIsotropic3D::pidLambda] = 
    _normalizer->dimensionalize(values[_ElasticIsotropic3D::pidLambda],
				   pressureScale);

  PetscLogFlops(3);
} // _dimProperties

// ----------------------------------------------------------------------
// Nondimensionalize initial state.
void
pylith::materials::ElasticIsotropic3D::_nondimInitState(double* const values,
							const int nvalues) const
{ // _nondimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticIsotropic3D::numInitialStateDBValues);

  const double pressureScale = _normalizer->pressureScale();
  _normalizer->nondimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _nondimInitState

// ----------------------------------------------------------------------
// Dimensionalize initial state.
void
pylith::materials::ElasticIsotropic3D::_dimInitState(double* const values,
						     const int nvalues) const
{ // _dimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticIsotropic3D::numInitialStateDBValues);
  
  const double pressureScale = _normalizer->pressureScale();
  _normalizer->dimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _dimInitState

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::ElasticIsotropic3D::_calcDensity(
				  double* const density,
				  const double* properties,
				  const int numProperties)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  density[0] = properties[_ElasticIsotropic3D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties.
void
pylith::materials::ElasticIsotropic3D::_calcStress(
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
  assert(_ElasticIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_ElasticIsotropic3D::tensorSize == strainSize);
  assert(_ElasticIsotropic3D::tensorSize == initialStateSize);

  const double density = properties[_ElasticIsotropic3D::pidDensity];
  const double mu = properties[_ElasticIsotropic3D::pidMu];
  const double lambda = properties[_ElasticIsotropic3D::pidLambda];

  const double mu2 = 2.0*mu;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double s123 = lambda * (e11 + e22 + e33);

  stress[0] = s123 + mu2*e11 + initialState[0];
  stress[1] = s123 + mu2*e22 + initialState[1];
  stress[2] = s123 + mu2*e33 + initialState[2];
  stress[3] = mu2 * e12 + initialState[3];
  stress[4] = mu2 * e23 + initialState[4];
  stress[5] = mu2 * e13 + initialState[5];

  PetscLogFlops(19);
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::ElasticIsotropic3D::_calcElasticConsts(
				  double* const elasticConsts,
				  const int numElasticConsts,
				  const double* properties,
				  const int numProperties,
				  const double*  totalStrain,
				  const int strainSize,
				  const double* initialState,
				  const int initialStateSize)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_ElasticIsotropic3D::tensorSize == strainSize);
  assert(_ElasticIsotropic3D::tensorSize == initialStateSize);
 
  const double density = properties[_ElasticIsotropic3D::pidDensity];
  const double mu = properties[_ElasticIsotropic3D::pidMu];
  const double lambda = properties[_ElasticIsotropic3D::pidLambda];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
   
  elasticConsts[ 0] = lambda2mu; // C1111
  elasticConsts[ 1] = lambda; // C1122
  elasticConsts[ 2] = lambda; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = 0; // C1123
  elasticConsts[ 5] = 0; // C1113
  elasticConsts[ 6] = lambda2mu; // C2222
  elasticConsts[ 7] = lambda; // C2233
  elasticConsts[ 8] = 0; // C2212
  elasticConsts[ 9] = 0; // C2223
  elasticConsts[10] = 0; // C2213
  elasticConsts[11] = lambda2mu; // C3333
  elasticConsts[12] = 0; // C3312
  elasticConsts[13] = 0; // C3323
  elasticConsts[14] = 0; // C3313
  elasticConsts[15] = mu2; // C1212
  elasticConsts[16] = 0; // C1223
  elasticConsts[17] = 0; // C1213
  elasticConsts[18] = mu2; // C2323
  elasticConsts[19] = 0; // C2313
  elasticConsts[20] = mu2; // C1313

  PetscLogFlops(2);
} // _calcElasticConsts

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::ElasticIsotropic3D::_stableTimeStepImplicit(const double* properties,
				 const int numProperties) const
{ // _stableTimeStepImplicit
  return pylith::PYLITH_MAXDOUBLE;
} // _stableTimeStepImplicit


// End of file 
