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

#include "MaxwellPlaneStrain.hh" // implementation of object methods

#include "ViscoelasticMaxwell.hh" // USES computeVisStrain

#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <cstring> // USES memcpy()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _MaxwellPlaneStrain{

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 4;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 6;

      /// Number of physical properties.
      const int numProperties = 6;

      /// Physical properties.
      const Material::PropMetaData properties[] = {
	{ "density", 1, OTHER_FIELD },
	{ "mu", 1, OTHER_FIELD },
	{ "lambda", 1, OTHER_FIELD },
	{ "maxwell_time", 1, OTHER_FIELD },
	{ "total_strain", 4, OTHER_FIELD },
	{ "viscous_strain", 4, OTHER_FIELD },
      };
      /// Indices (order) of properties.
      const int pidDensity = 0;
      const int pidMu = pidDensity + 1;
      const int pidLambda = pidMu + 1;
      const int pidMaxwellTime = pidLambda + 1;
      const int pidStrainT = pidMaxwellTime + 1;
      const int pidVisStrain = pidStrainT + tensorSize;

      /// Values expected in spatial database
      const int numDBValues = 4;
      const char* namesDBValues[] = {"density", "vs", "vp" , "viscosity"};

      /// Indices (order) of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;
      const int didViscosity = 3;

      /// Initial state values expected in spatial database
      const int numInitialStateDBValues = tensorSize;
      const char* namesInitialStateDBValues[] = { "stress_xx", "stress_yy",
                                                  "stress_zz", "stress_xy" };

    } // _MaxwellPlaneStrain
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::MaxwellPlaneStrain::MaxwellPlaneStrain(void) :
  ElasticMaterial(_MaxwellPlaneStrain::tensorSize,
		  _MaxwellPlaneStrain::numElasticConsts,
		  _MaxwellPlaneStrain::namesDBValues,
		  _MaxwellPlaneStrain::namesInitialStateDBValues,
		  _MaxwellPlaneStrain::numDBValues,
		  _MaxwellPlaneStrain::properties,
		  _MaxwellPlaneStrain::numProperties),
  _calcElasticConstsFn(&pylith::materials::MaxwellPlaneStrain::_calcElasticConstsElastic),
  _calcStressFn(&pylith::materials::MaxwellPlaneStrain::_calcStressElastic),
  _updatePropertiesFn(&pylith::materials::MaxwellPlaneStrain::_updatePropertiesElastic)
{ // constructor
  _dimension = 2;
  _visStrain.resize(_MaxwellPlaneStrain::tensorSize);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::MaxwellPlaneStrain::~MaxwellPlaneStrain(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::MaxwellPlaneStrain::_dbToProperties(
					    double* const propValues,
					    const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_MaxwellPlaneStrain::numDBValues == numDBValues);

  const double density = dbValues[_MaxwellPlaneStrain::didDensity];
  const double vs = dbValues[_MaxwellPlaneStrain::didVs];
  const double vp = dbValues[_MaxwellPlaneStrain::didVp];
  const double viscosity = dbValues[_MaxwellPlaneStrain::didViscosity];
 
  if (density <= 0.0 || vs <= 0.0 || vp <= 0.0 || viscosity <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for physical "
	<< "properties.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n"
	<< "viscosity: " << viscosity << "\n";
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
  assert(mu > 0);

  const double maxwelltime = viscosity / mu;
  assert(maxwelltime > 0.0);

  propValues[_MaxwellPlaneStrain::pidDensity] = density;
  propValues[_MaxwellPlaneStrain::pidMu] = mu;
  propValues[_MaxwellPlaneStrain::pidLambda] = lambda;
  propValues[_MaxwellPlaneStrain::pidMaxwellTime] = maxwelltime;

  PetscLogFlops(7);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::MaxwellPlaneStrain::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _totalPropsQuadPt);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();
  values[_MaxwellPlaneStrain::pidDensity] = 
    _normalizer->nondimensionalize(values[_MaxwellPlaneStrain::pidDensity],
				   densityScale);
  values[_MaxwellPlaneStrain::pidMu] = 
    _normalizer->nondimensionalize(values[_MaxwellPlaneStrain::pidMu],
				   pressureScale);
  values[_MaxwellPlaneStrain::pidLambda] = 
    _normalizer->nondimensionalize(values[_MaxwellPlaneStrain::pidLambda],
				   pressureScale);
  values[_MaxwellPlaneStrain::pidMaxwellTime] = 
    _normalizer->nondimensionalize(values[_MaxwellPlaneStrain::pidMaxwellTime],
				   timeScale);

  PetscLogFlops(4);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::MaxwellPlaneStrain::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _totalPropsQuadPt);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();
  values[_MaxwellPlaneStrain::pidDensity] = 
    _normalizer->dimensionalize(values[_MaxwellPlaneStrain::pidDensity],
				densityScale);
  values[_MaxwellPlaneStrain::pidMu] = 
    _normalizer->dimensionalize(values[_MaxwellPlaneStrain::pidMu],
				pressureScale);
  values[_MaxwellPlaneStrain::pidLambda] = 
    _normalizer->dimensionalize(values[_MaxwellPlaneStrain::pidLambda],
				pressureScale);
  values[_MaxwellPlaneStrain::pidMaxwellTime] = 
    _normalizer->dimensionalize(values[_MaxwellPlaneStrain::pidMaxwellTime],
				timeScale);

  PetscLogFlops(4);
} // _dimProperties

// ----------------------------------------------------------------------
// Nondimensionalize initial state.
void
pylith::materials::MaxwellPlaneStrain::_nondimInitState(double* const values,
							const int nvalues) const
{ // _nondimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _MaxwellPlaneStrain::numInitialStateDBValues);

  const double pressureScale = _normalizer->pressureScale();
  _normalizer->nondimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _nondimInitState

// ----------------------------------------------------------------------
// Dimensionalize initial state.
void
pylith::materials::MaxwellPlaneStrain::_dimInitState(double* const values,
						     const int nvalues) const
{ // _dimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _MaxwellPlaneStrain::numInitialStateDBValues);
  
  const double pressureScale = _normalizer->pressureScale();
  _normalizer->dimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _dimInitState

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::MaxwellPlaneStrain::_calcDensity(double* const density,
						    const double* properties,
						    const int numProperties)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  density[0] = properties[_MaxwellPlaneStrain::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute viscous strain for current time step.
// material.
void
pylith::materials::MaxwellPlaneStrain::_computeStateVars(
				         const double* properties,
					 const int numProperties,
					 const double* totalStrain,
					 const int strainSize,
					 const double* initialState,
					 const int initialStateSize)
{ // _computeStateVars
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(_MaxwellPlaneStrain::tensorSize == initialStateSize);

  const int tensorSize = _MaxwellPlaneStrain::tensorSize;
  const double maxwelltime = properties[_MaxwellPlaneStrain::pidMaxwellTime];

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  
  const double meanStrainTpdt = (e11 + e22 + e33)/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0 };

  const double meanStrainT =
    (properties[_MaxwellPlaneStrain::pidStrainT+0] +
     properties[_MaxwellPlaneStrain::pidStrainT+1] +
     properties[_MaxwellPlaneStrain::pidStrainT+2])/3.0;
  
  // Time integration.
  double dq = ViscoelasticMaxwell::computeVisStrain(_dt, maxwelltime);
  const double expFac = exp(-_dt/maxwelltime);

  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;

  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp] * meanStrainTpdt;
    devStrainT = properties[_MaxwellPlaneStrain::pidStrainT+iComp] -
      diag[iComp] * meanStrainT;
    _visStrain[iComp] = expFac *
      properties[_MaxwellPlaneStrain::pidVisStrain + iComp] +
      dq * (devStrainTpdt - devStrainT);
  } // for

  PetscLogFlops(8 + 7 * tensorSize);
} // _computeStateVars

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as an elastic
// material.
void
pylith::materials::MaxwellPlaneStrain::_calcStressElastic(
				         double* const stress,
					 const int stressSize,
					 const double* properties,
					 const int numProperties,
					 const double* totalStrain,
					 const int strainSize,
					 const double* initialState,
					 const int initialStateSize,
					 const bool computeStateVars)
{ // _calcStressElastic
  assert(0 != stress);
  assert(_MaxwellPlaneStrain::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(_MaxwellPlaneStrain::tensorSize == initialStateSize);

  const double mu = properties[_MaxwellPlaneStrain::pidMu];
  const double lambda = properties[_MaxwellPlaneStrain::pidLambda];
  const double mu2 = 2.0 * mu;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  
  const double traceStrainTpdt = e11 + e22 + e33;
  const double s123 = lambda * traceStrainTpdt;

  stress[0] = s123 + mu2*e11 + initialState[0];
  stress[1] = s123 + mu2*e22 + initialState[1];
  stress[2] = s123 +           initialState[2];
  stress[3] = mu2 * e12 + initialState[3];

  PetscLogFlops(13);
} // _calcStressElastic

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as a viscoelastic
// material.
void
pylith::materials::MaxwellPlaneStrain::_calcStressViscoelastic(
				         double* const stress,
					 const int stressSize,
					 const double* properties,
					 const int numProperties,
					 const double* totalStrain,
					 const int strainSize,
					 const double* initialState,
					 const int initialStateSize,
					 const bool computeStateVars)
{ // _calcStressViscoelastic
  assert(0 != stress);
  assert(_MaxwellPlaneStrain::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(_MaxwellPlaneStrain::tensorSize == initialStateSize);

  const int tensorSize = _MaxwellPlaneStrain::tensorSize;

  const double mu = properties[_MaxwellPlaneStrain::pidMu];
  const double lambda = properties[_MaxwellPlaneStrain::pidLambda];
  const double maxwelltime = properties[_MaxwellPlaneStrain::pidMaxwellTime];

  const double mu2 = 2.0 * mu;
  const double bulkModulus = lambda + mu2/3.0;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  
  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;
  const double meanStressTpdt = bulkModulus * traceStrainTpdt;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0 };

  // Get viscous strains
  if (computeStateVars) {
    pylith::materials::MaxwellPlaneStrain::_computeStateVars(properties,
							     numProperties,
							     totalStrain,
							     strainSize,
							     initialState,
							     initialStateSize);
  } else {
    memcpy(&_visStrain[0], &properties[_MaxwellPlaneStrain::pidVisStrain],
	   tensorSize * sizeof(double));
  } // else

  // Compute new stresses
  double devStressTpdt = 0.0;

  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStressTpdt = mu2 * _visStrain[iComp];

    stress[iComp] = diag[iComp] * meanStressTpdt + devStressTpdt +
	    initialState[iComp];
  } // for

  PetscLogFlops(7 + 4 * tensorSize);
} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::MaxwellPlaneStrain::_calcElasticConstsElastic(
				         double* const elasticConsts,
					 const int numElasticConsts,
					 const double* properties,
					 const int numProperties,
					 const double* totalStrain,
					 const int strainSize,
					 const double* initialState,
					 const int initialStateSize)
{ // _calcElasticConstsElastic
  assert(0 != elasticConsts);
  assert(_MaxwellPlaneStrain::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(_MaxwellPlaneStrain::tensorSize == initialStateSize);
 
  const double mu = properties[_MaxwellPlaneStrain::pidMu];
  const double lambda = properties[_MaxwellPlaneStrain::pidLambda];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;

  elasticConsts[ 0] = lambda2mu; // C1111
  elasticConsts[ 1] = lambda; // C1122
  elasticConsts[ 2] = lambda; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = lambda2mu; // C2222
  elasticConsts[ 5] = lambda; // C2233
  elasticConsts[ 6] = 0; // C2212
  elasticConsts[ 7] = lambda2mu; // C3333
  elasticConsts[ 8] = 0; // C3312
  elasticConsts[ 9] = mu2; // C1212

  PetscLogFlops(2);
} // _calcElasticConstsElastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties
// as an elastic material.
void
pylith::materials::MaxwellPlaneStrain::_calcElasticConstsViscoelastic(
				         double* const elasticConsts,
					 const int numElasticConsts,
					 const double* properties,
					 const int numProperties,
					 const double* totalStrain,
					 const int strainSize,
					 const double* initialState,
					 const int initialStateSize)
{ // _calcElasticConstsViscoelastic
  assert(0 != elasticConsts);
  assert(_MaxwellPlaneStrain::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(_MaxwellPlaneStrain::tensorSize == initialStateSize);
 
  const double mu = properties[_MaxwellPlaneStrain::pidMu];
  const double lambda = properties[_MaxwellPlaneStrain::pidLambda];
  const double maxwelltime = properties[_MaxwellPlaneStrain::pidMaxwellTime];

  const double mu2 = 2.0 * mu;
  const double bulkModulus = lambda + mu2/3.0;

  double dq = ViscoelasticMaxwell::computeVisStrain(_dt, maxwelltime);

  const double visFac = mu*dq/3.0;
  elasticConsts[ 0] = bulkModulus + 4.0*visFac; // C1111
  elasticConsts[ 1] = bulkModulus - 2.0*visFac; // C1122
  elasticConsts[ 2] = elasticConsts[1]; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = elasticConsts[0]; // C2222
  elasticConsts[ 5] = elasticConsts[1]; // C2233
  elasticConsts[ 6] = 0; // C2212
  elasticConsts[ 7] = elasticConsts[0]; // C3333
  elasticConsts[ 8] = 0; // C3312
  elasticConsts[ 9] = 6.0 * visFac; // C1212

  PetscLogFlops(10);
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::MaxwellPlaneStrain::_stableTimeStepImplicit(const double* properties,
				 const int numProperties) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  const double maxwellTime = 
    properties[_MaxwellPlaneStrain::pidMaxwellTime];
  const double dtStable = 0.1*maxwellTime;

  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::MaxwellPlaneStrain::_updatePropertiesElastic(
				         double* const properties,
					 const int numProperties,
					 const double* totalStrain,
					 const int strainSize,
					 const double* initialState,
					 const int initialStateSize)
{ // _updatePropertiesElastic
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);

  const double maxwelltime = properties[_MaxwellPlaneStrain::pidMaxwellTime];

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0 };

  for (int iComp=0; iComp < _MaxwellPlaneStrain::tensorSize; ++iComp) {
    properties[_MaxwellPlaneStrain::pidStrainT+iComp] = totalStrain[iComp];
    properties[_MaxwellPlaneStrain::pidVisStrain+iComp] =
      totalStrain[iComp] - diag[iComp] * meanStrainTpdt;
  } // for
  PetscLogFlops(3 + 2 * _MaxwellPlaneStrain::tensorSize);

  _needNewJacobian = true;
} // _updatePropertiesElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::MaxwellPlaneStrain::_updatePropertiesViscoelastic(
						 double* const properties,
						 const int numProperties,
						 const double* totalStrain,
						 const int strainSize,
						 const double* initialState,
						 const int initialStateSize)
{ // _updatePropertiesViscoelastic
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(_MaxwellPlaneStrain::tensorSize == initialStateSize);

  const int tensorSize = _MaxwellPlaneStrain::tensorSize;

  pylith::materials::MaxwellPlaneStrain::_computeStateVars(properties,
							   numProperties,
							   totalStrain,
							   strainSize,
							   initialState,
							   initialStateSize);

  memcpy(&properties[_MaxwellPlaneStrain::pidVisStrain],
	 &_visStrain[0], 
	 tensorSize * sizeof(double));
  memcpy(&properties[_MaxwellPlaneStrain::pidStrainT],
	 &totalStrain[0], 
	 tensorSize * sizeof(double));

  _needNewJacobian = false;

} // _updatePropertiesViscoelastic


// End of file 
