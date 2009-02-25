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

#include "MaxwellIsotropic3D.hh" // implementation of object methods

#include "ViscoelasticMaxwell.hh" // USES computeVisStrain
#include "Metadata.hh" // USES Metadata

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
    namespace _MaxwellIsotropic3D{

      // Dimension of material.
      const int dimension = 3;

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 6;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 21;

      /// Number of physical properties.
      const int numProperties = 4;
      
      /// Physical properties.
      const Metadata::ParamDescription properties[] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
	{ "maxwell_time", 1, pylith::topology::FieldBase::SCALAR },
      };
	
      // Values expected in properties spatial database
      const int numDBProperties = 4;
      const char* dbProperties[] = {"density", "vs", "vp" , "viscosity"};

      /// Number of state variables.
      const int numStateVars = 6;
      
      /// State variables.
      const Metadata::ParamDescription stateVars[] = {
	{ "total_strain", 6, pylith::topology::FieldBase::TENSOR },
	{ "viscous_strain", 6, pylith::topology::FieldBase::TENSOR },
      };

      // Values expected in state variables spatial database
      const int numDBStateVars = 12;
      const char* dbStateVars[] = {"total-strain-xx",
				   "total-strain-yy",
				   "total-strain-zz",
				   "total-strain-xy",
				   "total-strain-yz",
				   "total-strain-xz",
				   "viscous-strain-xx",
				   "viscous-strain-yy",
				   "viscous-strain-zz",
				   "viscous-strain-xy",
				   "viscous-strain-yz",
				   "viscous-strain-xz",
      };

    } // _MaxwellIsotropic3D
  } // materials
} // pylith

// Indices of physical properties
const int pylith::materials::MaxwellIsotropic3D::p_density = 0;

const int pylith::materials::MaxwellIsotropic3D::p_mu = 
  pylith::materials::MaxwellIsotropic3D::p_density + 1;

const int pylith::materials::MaxwellIsotropic3D::p_lambda = 
  pylith::materials::MaxwellIsotropic3D::p_mu + 1;

const int pylith::materials::MaxwellIsotropic3D::p_maxwellTime = 
  pylith::materials::MaxwellIsotropic3D::p_maxwellTime + 1;

// Indices of database values (order must match dbProperties)
const int pylith::materials::MaxwellIsotropic3D::db_density = 0;

const int pylith::materials::MaxwellIsotropic3D::db_vs = 
  pylith::materials::MaxwellIsotropic3D::db_density + 1;

const int pylith::materials::MaxwellIsotropic3D::db_vp = 
  pylith::materials::MaxwellIsotropic3D::db_vs + 1;

const int pylith::materials::MaxwellIsotropic3D::db_viscosity = 
  pylith::materials::MaxwellIsotropic3D::db_vp + 1;

// Indices of state variables
const int pylith::materials::MaxwellIsotropic3D::s_totalStrain = 0;

const int pylith::materials::MaxwellIsotropic3D::s_viscousStrain = 
  pylith::materials::MaxwellIsotropic3D::s_totalStrain + 6;

// Indices of database values (order must match dbStateVars)
const int pylith::materials::MaxwellIsotropic3D::db_totalStrain = 0;

const int pylith::materials::MaxwellIsotropic3D::db_viscousStrain = 
  pylith::materials::MaxwellIsotropic3D::db_totalStrain + 6;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::MaxwellIsotropic3D::MaxwellIsotropic3D(void) :
  ElasticMaterial(_MaxwellIsotropic3D::dimension,
		  _MaxwellIsotropic3D::tensorSize,
		  _MaxwellIsotropic3D::numElasticConsts,
		  Metadata(_MaxwellIsotropic3D::properties,
			   _MaxwellIsotropic3D::numProperties,
			   _MaxwellIsotropic3D::dbProperties,
			   _MaxwellIsotropic3D::numDBProperties,
			   _MaxwellIsotropic3D::stateVars,
			   _MaxwellIsotropic3D::numStateVars,
			   _MaxwellIsotropic3D::dbStateVars,
			   _MaxwellIsotropic3D::numDBStateVars)),
  _calcElasticConstsFn(0),
  _calcStressFn(0),
  _updateStateVarsFn(0)
{ // constructor
  useElasticBehavior(true);
  _viscousStrain.resize(_tensorSize);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::MaxwellIsotropic3D::~MaxwellIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set whether elastic or inelastic constitutive relations are used.
void
pylith::materials::MaxwellIsotropic3D::useElasticBehavior(const bool flag)
{ // useElasticBehavior
  if (flag) {
    _calcStressFn = 
      &pylith::materials::MaxwellIsotropic3D::_calcStressElastic;
    _calcElasticConstsFn = 
      &pylith::materials::MaxwellIsotropic3D::_calcElasticConstsElastic;
    _updateStateVarsFn = 
      &pylith::materials::MaxwellIsotropic3D::_updateStateVarsElastic;

  } else {
    _calcStressFn = 
      &pylith::materials::MaxwellIsotropic3D::_calcStressViscoelastic;
    _calcElasticConstsFn = 
      &pylith::materials::MaxwellIsotropic3D::_calcElasticConstsViscoelastic;
    _updateStateVarsFn = 
      &pylith::materials::MaxwellIsotropic3D::_updateStateVarsViscoelastic;
  } // if/else
} // useElasticBehavior

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::MaxwellIsotropic3D::_dbToProperties(
					    double* const propValues,
					    const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_MaxwellIsotropic3D::numDBProperties == numDBValues);

  const double density = dbValues[db_density];
  const double vs = dbValues[db_vs];
  const double vp = dbValues[db_vp];
  const double viscosity = dbValues[db_viscosity];
 
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

  const double maxwellTime = viscosity / mu;
  assert(maxwellTime > 0.0);

  propValues[p_density] = density;
  propValues[p_mu] = mu;
  propValues[p_lambda] = lambda;
  propValues[p_maxwellTime] = maxwellTime;

  PetscLogFlops(7);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::MaxwellIsotropic3D::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numPropsQuadPt);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();
  values[p_density] = 
    _normalizer->nondimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->nondimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->nondimensionalize(values[p_lambda], pressureScale);
  values[p_maxwellTime] = 
    _normalizer->nondimensionalize(values[p_maxwellTime], timeScale);

  PetscLogFlops(4);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::MaxwellIsotropic3D::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numPropsQuadPt);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();
  values[p_density] = 
    _normalizer->dimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->dimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->dimensionalize(values[p_lambda], pressureScale);
  values[p_maxwellTime] = 
    _normalizer->dimensionalize(values[p_maxwellTime], timeScale);

  PetscLogFlops(4);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute initial state variables from values in spatial database.
void
pylith::materials::MaxwellIsotropic3D::_dbToStateVars(
					double* const stateValues,
					const double_array& dbValues) const
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_MaxwellIsotropic3D::numDBStateVars == numDBValues);

  const int totalSize = 2 * _tensorSize;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  memcpy(stateValues, &dbValues[0], totalSize*sizeof(double));

  PetscLogFlops(0);
} // _dbToStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::MaxwellIsotropic3D::_calcDensity(double* const density,
						    const double* properties,
						    const int numProperties,
						    const double* stateVars,
						    const int numStateVars)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);

  density[0] = properties[p_density];
} // _calcDensity

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::MaxwellIsotropic3D::_stableTimeStepImplicit(
					   const double* properties,
					   const int numProperties,
					   const double* stateVars,
					   const int numStateVars) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);

  const double maxwellTime = properties[p_maxwellTime];
  const double dtStable = 0.1 * maxwellTime;

  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Compute viscous strain for current time step.
// material.
void
pylith::materials::MaxwellIsotropic3D::_computeStateVars(
					       const double* stateVars,
					       const int numStateVars,
					       const double* properties,
					       const int numProperties,
					       const double* totalStrain,
					       const int strainSize,
					       const double* initialStress,
					       const int initialStressSize,
					       const double* initialStrain,
					       const int initialStrainSize)
{ // _computeStateVars
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;
  const double maxwellTime = properties[p_maxwellTime];

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double meanStrainTpdt = (e11 + e22 + e33)/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  const double meanStrainT =
    ( stateVars[s_totalStrain+0] +
      stateVars[s_totalStrain+1] +
      stateVars[s_totalStrain+2] ) / 3.0;
  
  // Time integration.
  double dq = ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTime);
  const double expFac = exp(-_dt/maxwellTime);

  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;

  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp] * meanStrainTpdt;
    devStrainT = stateVars[s_totalStrain+iComp] - diag[iComp] * meanStrainT;
    _viscousStrain[iComp] = expFac * stateVars[s_viscousStrain+iComp] +
      dq * (devStrainTpdt - devStrainT);
  } // for

  PetscLogFlops(8 + 7 * tensorSize);
} // _computeStateVars

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as an elastic
// material.
void
pylith::materials::MaxwellIsotropic3D::_calcStressElastic(
					     double* const stress,
					     const int stressSize,
					     const double* properties,
					     const int numProperties,
					     const double* stateVars,
					     const int numStateVars,
					     const double* totalStrain,
					     const int strainSize,
					     const double* initialStress,
					     const int initialStressSize,
					     const double* initialStrain,
					     const int initialStrainSize,
					     const bool computeStateVars)
{ // _calcStressElastic
  assert(0 != stress);
  assert(_MaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];
  const double mu2 = 2.0 * mu;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double traceStrainTpdt = e11 + e22 + e33;
  const double s123 = lambda * traceStrainTpdt;

  stress[0] = s123 + mu2*e11 + initialStress[0];
  stress[1] = s123 + mu2*e22 + initialStress[1];
  stress[2] = s123 + mu2*e33 + initialStress[2];
  stress[3] = mu2 * e12 + initialStress[3];
  stress[4] = mu2 * e23 + initialStress[4];
  stress[5] = mu2 * e13 + initialStress[5];

  PetscLogFlops(19);
} // _calcStressElastic

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as a viscoelastic
// material.
void
pylith::materials::MaxwellIsotropic3D::_calcStressViscoelastic(
					     double* const stress,
					     const int stressSize,
					     const double* properties,
					     const int numProperties,
					     const double* stateVars,
					     const int numStateVars,
					     const double* totalStrain,
					     const int strainSize,
					     const double* initialStress,
					     const int initialStressSize,
					     const double* initialStrain,
					     const int initialStrainSize,
					     const bool computeStateVars)
{ // _calcStressViscoelastic
  assert(0 != stress);
  assert(_MaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _MaxwellIsotropic3D::tensorSize;

  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];
  const double maxwellTime = properties[p_maxwellTime];

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
  if (computeStateVars)
    _computeStateVars(stateVars, numStateVars,
		      properties, numProperties,
		      totalStrain, strainSize,
		      initialStress, initialStressSize,
		      initialStrain, initialStrainSize);
  else
    memcpy(&_viscousStrain[0], &stateVars[s_viscousStrain],
	   tensorSize*sizeof(double));

  // Compute new stresses
  double devStressTpdt = 0.0;

  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStressTpdt = mu2 * _viscousStrain[iComp];

    stress[iComp] = diag[iComp] * meanStressTpdt + devStressTpdt +
	    initialStress[iComp];
  } // for

  PetscLogFlops(7 + 4 * tensorSize);
} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::MaxwellIsotropic3D::_calcElasticConstsElastic(
					double* const elasticConsts,
					const int numElasticConsts,
					const double* properties,
					const int numProperties,
					const double* stateVars,
					const int numStateVars,
					const double* totalStrain,
					const int strainSize,
					const double* initialStress,
					const int initialStressSize,
					const double* initialStrain,
					const int initialStrainSize)
{ // _calcElasticConstsElastic
  assert(0 != elasticConsts);
  assert(_MaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);
 
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];

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
pylith::materials::MaxwellIsotropic3D::_calcElasticConstsViscoelastic(
					double* const elasticConsts,
					const int numElasticConsts,
					const double* properties,
					const int numProperties,
					const double* stateVars,
					const int numStateVars,
					const double* totalStrain,
					const int strainSize,
					const double* initialStress,
					const int initialStressSize,
					const double* initialStrain,
					const int initialStrainSize)
{ // _calcElasticConstsViscoelastic
  assert(0 != elasticConsts);
  assert(_MaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];
  const double maxwellTime = properties[p_maxwellTime];

  const double mu2 = 2.0 * mu;
  const double bulkModulus = lambda + mu2 / 3.0;

  double dq = ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTime);

  const double visFac = mu * dq / 3.0;
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
// Update state variables.
void
pylith::materials::MaxwellIsotropic3D::_updateStateVarsElastic(
					    double* const stateVars,
					    const int numStateVars,
					    const double* properties,
					    const int numProperties,
					    const double* totalStrain,
					    const int strainSize,
					    const double* initialStress,
					    const int initialStressSize,
					    const double* initialStrain,
					    const int initialStrainSize)
{ // _updateStateVarsElastic
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;
  const double maxwellTime = properties[p_maxwellTime];

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt / 3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // :TODO: Need to account for initial values for state variables
  // and the initial strain??
  for (int iComp=0; iComp < tensorSize; ++iComp) {
    stateVars[s_totalStrain+iComp] = totalStrain[iComp];
    stateVars[s_viscousStrain+iComp] =
      totalStrain[iComp] - diag[iComp] * meanStrainTpdt;
  } // for
  PetscLogFlops(3 + 2 * _tensorSize);

  _needNewJacobian = true;
} // _updateStateVarsElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::MaxwellIsotropic3D::_updateStateVarsViscoelastic(
					    double* const stateVars,
					    const int numStateVars,
					    const double* properties,
					    const int numProperties,
					    const double* totalStrain,
					    const int strainSize,
					    const double* initialStress,
					    const int initialStressSize,
					    const double* initialStrain,
					    const int initialStrainSize)
{ // _updateStateVarsViscoelastic
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;

  _computeStateVars(stateVars, numStateVars,
		    properties, numProperties,
		    totalStrain, strainSize,
		    initialStress, initialStressSize,
		    initialStrain, initialStrainSize);

  memcpy(&stateVars[s_totalStrain], totalStrain, tensorSize*sizeof(double));

  memcpy(&stateVars[s_viscousStrain], &_viscousStrain[0], 
	 tensorSize*sizeof(double));

  _needNewJacobian = false;
} // _updateStateVarsViscoelastic


// End of file 
