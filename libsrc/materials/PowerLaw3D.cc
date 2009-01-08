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

#include "PowerLaw3D.hh" // implementation of object methods

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
    namespace _PowerLaw3D{

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 6;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 21;

      /// Number of physical properties.
      const int numProperties = 7;

      /// Physical properties.
      const Material::PropMetaData properties[] = {
	{ "density", 1, OTHER_FIELD },
	{ "mu", 1, OTHER_FIELD },
	{ "lambda", 1, OTHER_FIELD },
	{ "viscosity_coeff", 1, OTHER_FIELD },
	{ "power_law_exponent", 1, OTHER_FIELD },
	{ "total_strain", 6, OTHER_FIELD },
	{ "viscous_strain", 6, OTHER_FIELD },
      };
      /// Indices (order) of properties.
      const int pidDensity = 0;
      const int pidMu = pidDensity + 1;
      const int pidLambda = pidMu + 1;
      const int pidViscosityCoeff = pidLambda + 1;
      const int pidPowerLawExp = pidViscosityCoeff + 1;
      const int pidStrainT = pidPowerLawExp + 1;
      const int pidVisStrain = pidStrainT + tensorSize;

      /// Values expected in spatial database
      const int numDBValues = 5;
      const char* namesDBValues[] = {"density", "vs", "vp" , "viscosity_coeff",
				     "power_law_exponent"};

      /// Indices (order) of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;
      const int didViscosityCoeff = 3;
      const int didPowerLawExp = 4;

      /// Initial state values expected in spatial database
      const int numInitialStateDBValues = tensorSize;
      const char* namesInitialStateDBValues[] = { "stress_xx", "stress_yy",
                                                  "stress_zz", "stress_xy",
                                                  "stress_yz", "stress_xz" };

    } // _PowerLaw3D
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::PowerLaw3D::PowerLaw3D(void) :
  ElasticMaterial(_PowerLaw3D::tensorSize,
		  _PowerLaw3D::numElasticConsts,
		  _PowerLaw3D::namesDBValues,
		  _PowerLaw3D::namesInitialStateDBValues,
		  _PowerLaw3D::numDBValues,
		  _PowerLaw3D::properties,
		  _PowerLaw3D::numProperties),
  _calcElasticConstsFn(&pylith::materials::PowerLaw3D::_calcElasticConstsElastic),
  _calcStressFn(&pylith::materials::PowerLaw3D::_calcStressElastic),
  _updatePropertiesFn(&pylith::materials::PowerLaw3D::_updatePropertiesElastic)
{ // constructor
  _dimension = 3;
  _visStrain.resize(_PowerLaw3D::tensorSize);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::PowerLaw3D::~PowerLaw3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::PowerLaw3D::_dbToProperties(
					    double* const propValues,
					    const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_PowerLaw3D::numDBValues == numDBValues);

  const double density = dbValues[_PowerLaw3D::didDensity];
  const double vs = dbValues[_PowerLaw3D::didVs];
  const double vp = dbValues[_PowerLaw3D::didVp];
  const double viscosityCoeff = dbValues[_PowerLaw3D::didViscosityCoeff];
  const double powerLawExp = dbValues[_PowerLaw3D::didPowerLawExp];
 
  if (density <= 0.0 || vs <= 0.0 || vp <= 0.0 || viscosityCoeff <= 0.0
      || powerLawExp < 1.0) {
    std::ostringstream msg;
    msg << "Spatial database returned illegal value for physical "
	<< "properties.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n"
	<< "viscosityCoeff: " << viscosityCoeff << "\n"
	<< "powerLawExp: " << powerLawExp << "\n";
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

  propValues[_PowerLaw3D::pidDensity] = density;
  propValues[_PowerLaw3D::pidMu] = mu;
  propValues[_PowerLaw3D::pidLambda] = lambda;
  propValues[_PowerLaw3D::pidViscosityCoeff] = viscosityCoeff;
  propValues[_PowerLaw3D::pidPowerLawExp] = powerLawExp;

  PetscLogFlops(6);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::PowerLaw3D::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _totalPropsQuadPt);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();
  // **** NOTE:  Make sure scaling is correct for viscosity coefficient.
  const double powerLawExp = _PowerLaw3D::pidPowerLawExp;
  const double viscosityCoeffScale = timeScale *
    pressureScale^(1.0/powerLawExp);
  const double powerLawExpScale = 1.0;
  values[_PowerLaw3D::pidDensity] = 
    _normalizer->nondimensionalize(values[_PowerLaw3D::pidDensity],
				   densityScale);
  values[_PowerLaw3D::pidMu] = 
    _normalizer->nondimensionalize(values[_PowerLaw3D::pidMu],
				   pressureScale);
  values[_PowerLaw3D::pidLambda] = 
    _normalizer->nondimensionalize(values[_PowerLaw3D::pidLambda],
				   pressureScale);
  values[_PowerLaw3D::pidViscosityCoeff] = 
    _normalizer->nondimensionalize(values[_PowerLaw3D::pidViscosityCoeff],
				   viscosityCoeffScale);
  values[_PowerLaw3D::pidPowerLawExp] = 
    _normalizer->nondimensionalize(values[_PowerLaw3D::pidPowerLawExp],
				   powerLawExpScale);

  PetscLogFlops(8);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::PowerLaw3D::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _totalPropsQuadPt);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();
  // **** NOTE:  Make sure scaling is correct for viscosity coefficient.
  const double powerLawExp = _PowerLaw3D::pidPowerLawExp;
  const double viscosityCoeffScale = timeScale *
    pressureScale^(1.0/powerLawExp);
  const double powerLawExpScale = 1.0;
  values[_PowerLaw3D::pidDensity] = 
    _normalizer->dimensionalize(values[_PowerLaw3D::pidDensity],
				densityScale);
  values[_PowerLaw3D::pidMu] = 
    _normalizer->dimensionalize(values[_PowerLaw3D::pidMu],
				pressureScale);
  values[_PowerLaw3D::pidLambda] = 
    _normalizer->dimensionalize(values[_PowerLaw3D::pidLambda],
				pressureScale);
  values[_PowerLaw3D::pidViscosityCoeff] = 
    _normalizer->dimensionalize(values[_PowerLaw3D::pidViscosityCoeff],
				viscosityCoeffScale);
  values[_PowerLaw3D::pidPowerLawExp] = 
    _normalizer->dimensionalize(values[_PowerLaw3D::pidPowerLawExp],
				powerLawExpScale);

  PetscLogFlops(8);
} // _dimProperties

// ----------------------------------------------------------------------
// Nondimensionalize initial state.
void
pylith::materials::PowerLaw3D::_nondimInitState(double* const values,
							const int nvalues) const
{ // _nondimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _PowerLaw3D::numInitialStateDBValues);

  const double pressureScale = _normalizer->pressureScale();
  _normalizer->nondimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _nondimInitState

// ----------------------------------------------------------------------
// Dimensionalize initial state.
void
pylith::materials::PowerLaw3D::_dimInitState(double* const values,
						     const int nvalues) const
{ // _dimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _PowerLaw3D::numInitialStateDBValues);
  
  const double pressureScale = _normalizer->pressureScale();
  _normalizer->dimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _dimInitState

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::PowerLaw3D::_calcDensity(double* const density,
						    const double* properties,
						    const int numProperties)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  density[0] = properties[_PowerLaw3D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute viscous strain for current time step.
// material.
void
pylith::materials::PowerLaw3D::_computeStateVars(
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
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(_PowerLaw3D::tensorSize == initialStateSize);

  const int tensorSize = _PowerLaw3D::tensorSize;
  const double maxwelltime = properties[_PowerLaw3D::pidMaxwellTime];

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double meanStrainTpdt = (e11 + e22 + e33)/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  const double meanStrainT =
    (properties[_PowerLaw3D::pidStrainT+0] +
     properties[_PowerLaw3D::pidStrainT+1] +
     properties[_PowerLaw3D::pidStrainT+2])/3.0;
  
  // Time integration.
  double dq = ViscoelasticMaxwell::computeVisStrain(_dt, maxwelltime);
  const double expFac = exp(-_dt/maxwelltime);

  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;

  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp] * meanStrainTpdt;
    devStrainT = properties[_PowerLaw3D::pidStrainT+iComp] -
      diag[iComp] * meanStrainT;
    _visStrain[iComp] = expFac *
      properties[_PowerLaw3D::pidVisStrain + iComp] +
      dq * (devStrainTpdt - devStrainT);
  } // for

  PetscLogFlops(8 + 7 * tensorSize);
} // _computeStateVars

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as an elastic
// material.
void
pylith::materials::PowerLaw3D::_calcStressElastic(
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
  assert(_PowerLaw3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(_PowerLaw3D::tensorSize == initialStateSize);

  const double mu = properties[_PowerLaw3D::pidMu];
  const double lambda = properties[_PowerLaw3D::pidLambda];
  const double mu2 = 2.0 * mu;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double traceStrainTpdt = e11 + e22 + e33;
  const double s123 = lambda * traceStrainTpdt;

  stress[0] = s123 + mu2*e11 + initialState[0];
  stress[1] = s123 + mu2*e22 + initialState[1];
  stress[2] = s123 + mu2*e33 + initialState[2];
  stress[3] = mu2 * e12 + initialState[3];
  stress[4] = mu2 * e23 + initialState[4];
  stress[5] = mu2 * e13 + initialState[5];

  PetscLogFlops(19);
} // _calcStressElastic

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as a viscoelastic
// material.
void
pylith::materials::PowerLaw3D::_calcStressViscoelastic(
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
  assert(_PowerLaw3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(_PowerLaw3D::tensorSize == initialStateSize);

  const int tensorSize = _PowerLaw3D::tensorSize;

  const double mu = properties[_PowerLaw3D::pidMu];
  const double lambda = properties[_PowerLaw3D::pidLambda];
  const double viscosityCoeff = properties[_PowerLaw3D::pidViscosityCoeff];
  const double powerLawExp = properties[_PowerLaw3D::pidPowerLawExp];

  const double mu2 = 2.0 * mu;
  const double bulkModulus = lambda + mu2/3.0;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  
  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;
  const double meanStressTpdt = bulkModulus * traceStrainTpdt;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Get viscous strains
  if (computeStateVars) {
    pylith::materials::PowerLaw3D::_computeStateVars(properties,
							     numProperties,
							     totalStrain,
							     strainSize,
							     initialState,
							     initialStateSize);
  } else {
    memcpy(&_visStrain[0], &properties[_PowerLaw3D::pidVisStrain],
	   tensorSize * sizeof(double));
  } // else

  // Compute new stresses
  double devStressTpdt = 0.0;

  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStressTpdt = mu2 * _visStrain[iComp];

    // Later I will want to put in initial stresses.
    stress[iComp] = diag[iComp] * meanStressTpdt + devStressTpdt +
	    initialState[iComp];
  } // for

  PetscLogFlops(7 + 4 * tensorSize);
} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::PowerLaw3D::_calcElasticConstsElastic(
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
  assert(_PowerLaw3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(_PowerLaw3D::tensorSize == initialStateSize);
 
  const double mu = properties[_PowerLaw3D::pidMu];
  const double lambda = properties[_PowerLaw3D::pidLambda];

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

  PetscLogFlops(4);
} // _calcElasticConstsElastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties
// as an elastic material.
void
pylith::materials::PowerLaw3D::_calcElasticConstsViscoelastic(
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
  assert(_PowerLaw3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(_PowerLaw3D::tensorSize == initialStateSize);
 
  const double mu = properties[_PowerLaw3D::pidMu];
  const double lambda = properties[_PowerLaw3D::pidLambda];
  const double maxwelltime = properties[_PowerLaw3D::pidMaxwellTime];

  const double mu2 = 2.0 * mu;
  const double bulkModulus = lambda + mu2/3.0;

  double dq = ViscoelasticMaxwell::computeVisStrain(_dt, maxwelltime);

  const double visFac = mu*dq/3.0;
  elasticConsts[ 0] = bulkModulus + 4.0*visFac; // C1111
  elasticConsts[ 1] = bulkModulus - 2.0*visFac; // C1122
  elasticConsts[ 2] = elasticConsts[1]; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = 0; // C1123
  elasticConsts[ 5] = 0; // C1113
  elasticConsts[ 6] = elasticConsts[0]; // C2222
  elasticConsts[ 7] = elasticConsts[1]; // C2233
  elasticConsts[ 8] = 0; // C2212
  elasticConsts[ 9] = 0; // C2223
  elasticConsts[10] = 0; // C2213
  elasticConsts[11] = elasticConsts[0]; // C3333
  elasticConsts[12] = 0; // C3312
  elasticConsts[13] = 0; // C3323
  elasticConsts[14] = 0; // C3313
  elasticConsts[15] = 6.0 * visFac; // C1212
  elasticConsts[16] = 0; // C1223
  elasticConsts[17] = 0; // C1213
  elasticConsts[18] = elasticConsts[15]; // C2323
  elasticConsts[19] = 0; // C2313
  elasticConsts[20] = elasticConsts[15]; // C1313

  PetscLogFlops(10);
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::PowerLaw3D::_stableTimeStepImplicit(const double* properties,
				 const int numProperties) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  const double maxwellTime = 
    properties[_PowerLaw3D::pidMaxwellTime];
  const double dtStable = 0.1*maxwellTime;

  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::PowerLaw3D::_updatePropertiesElastic(
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
  assert(_PowerLaw3D::tensorSize == strainSize);

  const double maxwelltime = properties[_PowerLaw3D::pidMaxwellTime];

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  for (int iComp=0; iComp < _PowerLaw3D::tensorSize; ++iComp) {
    properties[_PowerLaw3D::pidStrainT+iComp] = totalStrain[iComp];
    properties[_PowerLaw3D::pidVisStrain+iComp] =
      totalStrain[iComp] - diag[iComp] * meanStrainTpdt;
  } // for
  PetscLogFlops(3 + 2 * _PowerLaw3D::tensorSize);

  _needNewJacobian = true;
} // _updatePropertiesElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::PowerLaw3D::_updatePropertiesViscoelastic(
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
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(_PowerLaw3D::tensorSize == initialStateSize);

  const int tensorSize = _PowerLaw3D::tensorSize;

  pylith::materials::PowerLaw3D::_computeStateVars(properties,
							   numProperties,
							   totalStrain,
							   strainSize,
							   initialState,
							   initialStateSize);

  memcpy(&properties[_PowerLaw3D::pidVisStrain],
	 &_visStrain[0], 
	 tensorSize * sizeof(double));
  memcpy(&properties[_PowerLaw3D::pidStrainT],
	 &totalStrain[0], 
	 tensorSize * sizeof(double));

  _needNewJacobian = false;

} // _updatePropertiesViscoelastic


// End of file 
