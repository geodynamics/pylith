// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "MaxwellPlaneStrain.hh" // implementation of object methods

#include "ViscoelasticMaxwell.hh" // USES computeVisStrain
#include "Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

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

      // Dimension of material
      const int dimension = 2;

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 3;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 9;

      /// Number of physical properties.
      const int numProperties = 4;

      /// Physical properties.
      const Metadata::ParamDescription properties[] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
	{ "maxwell_time", 1, pylith::topology::FieldBase::SCALAR },
      };

      /// Values expected in properties spatial database
      const int numDBProperties = 4;
      const char* dbProperties[] = {"density", "vs", "vp" , "viscosity"};

      /// Number of state variables.
      const int numStateVars = 3;

      /// State variables.
      const Metadata::ParamDescription stateVars[] = {
	{ "stress_zz_initial", 1, pylith::topology::FieldBase::SCALAR },
	{ "total_strain", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "viscous_strain", 4, pylith::topology::FieldBase::OTHER },
      };

      /// Values expected in state variables spatial database
      const int numDBStateVars = 1 + tensorSize + 4;
      const char* dbStateVars[] = { "stress-zz-initial",
				    "total-strain-xx",
				    "total-strain-yy",
				    "total-strain-xy",
				    "viscous-strain-xx",
				    "viscous-strain-yy",
				    "viscous-strain-zz",
				    "viscous-strain-xy" };

    } // _MaxwellPlaneStrain
  } // materials
} // pylith

// Indices of physical properties
const int pylith::materials::MaxwellPlaneStrain::p_density = 0;

const int pylith::materials::MaxwellPlaneStrain::p_mu = 
  pylith::materials::MaxwellPlaneStrain::p_density + 1;

const int pylith::materials::MaxwellPlaneStrain::p_lambda = 
  pylith::materials::MaxwellPlaneStrain::p_mu + 1;

const int pylith::materials::MaxwellPlaneStrain::p_maxwellTime = 
  pylith::materials::MaxwellPlaneStrain::p_lambda + 1;

// Indices of database values (order must match dbProperties)
const int pylith::materials::MaxwellPlaneStrain::db_density = 0;

const int pylith::materials::MaxwellPlaneStrain::db_vs = 
  pylith::materials::MaxwellPlaneStrain::db_density + 1;

const int pylith::materials::MaxwellPlaneStrain::db_vp = 
  pylith::materials::MaxwellPlaneStrain::db_vs + 1;

const int pylith::materials::MaxwellPlaneStrain::db_viscosity = 
  pylith::materials::MaxwellPlaneStrain::db_vp + 1;

// Indices of state variables
const int pylith::materials::MaxwellPlaneStrain::s_stressZZInitial = 0;

const int pylith::materials::MaxwellPlaneStrain::s_totalStrain =
  pylith::materials::MaxwellPlaneStrain::s_stressZZInitial + 1;

const int pylith::materials::MaxwellPlaneStrain::s_viscousStrain = 
  pylith::materials::MaxwellPlaneStrain::s_totalStrain + 3;

// Indices of database values (order must match dbStateVars)
const int pylith::materials::MaxwellPlaneStrain::db_stressZZInitial = 0;

const int pylith::materials::MaxwellPlaneStrain::db_totalStrain =
  pylith::materials::MaxwellPlaneStrain::db_stressZZInitial + 1;

const int pylith::materials::MaxwellPlaneStrain::db_viscousStrain = 
  pylith::materials::MaxwellPlaneStrain::db_totalStrain + 3;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::MaxwellPlaneStrain::MaxwellPlaneStrain(void) :
  ElasticMaterial(_MaxwellPlaneStrain::dimension,
		  _MaxwellPlaneStrain::tensorSize,
		  _MaxwellPlaneStrain::numElasticConsts,
		  Metadata(_MaxwellPlaneStrain::properties,
			   _MaxwellPlaneStrain::numProperties,
			   _MaxwellPlaneStrain::dbProperties,
			   _MaxwellPlaneStrain::numDBProperties,
			   _MaxwellPlaneStrain::stateVars,
			   _MaxwellPlaneStrain::numStateVars,
			   _MaxwellPlaneStrain::dbStateVars,
			   _MaxwellPlaneStrain::numDBStateVars)),
  _calcElasticConstsFn(0),
  _calcStressFn(0),
  _updateStateVarsFn(0)
{ // constructor
  useElasticBehavior(true);
  _viscousStrain.resize(4);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::MaxwellPlaneStrain::~MaxwellPlaneStrain(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set whether elastic or inelastic constitutive relations are used.
void
pylith::materials::MaxwellPlaneStrain::useElasticBehavior(const bool flag)
{ // useElasticBehavior
  if (flag) {
    _calcStressFn = 
      &pylith::materials::MaxwellPlaneStrain::_calcStressElastic;
    _calcElasticConstsFn = 
      &pylith::materials::MaxwellPlaneStrain::_calcElasticConstsElastic;
    _updateStateVarsFn = 
      &pylith::materials::MaxwellPlaneStrain::_updateStateVarsElastic;

  } else {
    _calcStressFn = 
      &pylith::materials::MaxwellPlaneStrain::_calcStressViscoelastic;
    _calcElasticConstsFn = 
      &pylith::materials::MaxwellPlaneStrain::_calcElasticConstsViscoelastic;
    _updateStateVarsFn = 
      &pylith::materials::MaxwellPlaneStrain::_updateStateVarsViscoelastic;
  } // if/else
} // useElasticBehavior

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::MaxwellPlaneStrain::_dbToProperties(
					    double* const propValues,
					    const double_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_MaxwellPlaneStrain::numDBProperties == numDBValues);

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

  const double mu = density * vs * vs;
  const double lambda = density * vp * vp - 2.0 * mu;

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
pylith::materials::MaxwellPlaneStrain::_nondimProperties(double* const values,
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
pylith::materials::MaxwellPlaneStrain::_dimProperties(double* const values,
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
pylith::materials::MaxwellPlaneStrain::_dbToStateVars(
					double* const stateValues,
					const double_array& dbValues)
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_MaxwellPlaneStrain::numDBStateVars == numDBValues);

  const int totalSize = 1 + _tensorSize + 4;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  memcpy(stateValues, &dbValues[0], totalSize*sizeof(double));

  PetscLogFlops(0);
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::materials::MaxwellPlaneStrain::_nondimStateVars(double* const values,
							const int nvalues) const
{ // _nondimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numVarsQuadPt);

  const double pressureScale = _normalizer->pressureScale();
  _normalizer->nondimensionalize(&values[s_stressZZInitial], 1, pressureScale);

  PetscLogFlops(1);
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::materials::MaxwellPlaneStrain::_dimStateVars(double* const values,
						     const int nvalues) const
{ // _dimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numVarsQuadPt);

  const double pressureScale = _normalizer->pressureScale();
  _normalizer->dimensionalize(&values[s_stressZZInitial], 1, pressureScale);

  PetscLogFlops(1);
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::MaxwellPlaneStrain::_calcDensity(double* const density,
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
// Compute stress tensor at location from properties as an elastic
// material.
void
pylith::materials::MaxwellPlaneStrain::_calcStressElastic(
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
  assert(_MaxwellPlaneStrain::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellPlaneStrain::tensorSize == initialStrainSize);

  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];
  const double mu2 = 2.0 * mu;

  const double e11 = totalStrain[0] - initialStrain[0];
  const double e22 = totalStrain[1] - initialStrain[1];
  const double e12 = totalStrain[2] - initialStrain[2];
  
  const double s12 = lambda * (e11 + e22);

  stress[0] = s12 + mu2 * e11 + initialStress[0];
  stress[1] = s12 + mu2 * e22 + initialStress[1];
  stress[2] = mu2 * e12 + initialStress[2];

  PetscLogFlops(14);
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
  assert(_MaxwellPlaneStrain::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellPlaneStrain::tensorSize == initialStrainSize);

  const int tensorSize = _MaxwellPlaneStrain::tensorSize;

  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];
  const double maxwellTime = properties[p_maxwellTime];
  const double stressZZInitial = stateVars[s_stressZZInitial];

  const double mu2 = 2.0 * mu;
  const double bulkModulus = lambda + mu2 / 3.0;

  // Initial stress and strain values
  const double meanStrainInitial = (initialStrain[0] + initialStrain[1]) / 3.0;
  const double meanStressInitial = (initialStress[0] + initialStress[1] +
				    stressZZInitial) / 3.0;
  const double devStrainInitial[] = {initialStrain[0] - meanStrainInitial,
				     initialStrain[1] - meanStrainInitial,
				     initialStrain[2]};
  const double devStressInitial[] = {initialStress[0] - meanStressInitial,
				     initialStress[1] - meanStressInitial,
				     initialStress[2]};

  // Mean stress and strain for t + dt
  const double meanStrainTpdt = (totalStrain[0] + totalStrain[1]) / 3.0;
  const double meanStressTpdt = 3.0 * bulkModulus *
    (meanStrainTpdt - meanStrainInitial) + meanStressInitial;

  // Get viscous strains
  if (computeStateVars) {
    _computeStateVars(stateVars, numStateVars,
		      properties, numProperties,
		      totalStrain, strainSize,
		      initialStress, initialStressSize,
		      initialStrain, initialStrainSize);
  } else {
    memcpy(&_viscousStrain[0], &stateVars[s_viscousStrain],
	   4 * sizeof(double));
  } // else

  // Compute new stresses
  stress[0] = meanStressTpdt + mu2 * (_viscousStrain[0] - devStrainInitial[0]) +
    devStressInitial[0];
  stress[1] = meanStressTpdt + mu2 * (_viscousStrain[1] - devStrainInitial[1]) +
    devStressInitial[1];
  stress[2] = mu2 * (_viscousStrain[3] - devStrainInitial[2]) +
    devStressInitial[2];

  PetscLogFlops(30);
} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::MaxwellPlaneStrain::_calcElasticConstsElastic(
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
  assert(_MaxwellPlaneStrain::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellPlaneStrain::tensorSize == initialStrainSize);
 
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;

  elasticConsts[ 0] = lambda2mu; // C1111
  elasticConsts[ 1] = lambda; // C1122
  elasticConsts[ 2] = 0; // C1112
  elasticConsts[ 3] = lambda; // C2211
  elasticConsts[ 4] = lambda2mu; // C2222
  elasticConsts[ 5] = 0; // C2212
  elasticConsts[ 6] = 0; // C1211
  elasticConsts[ 7] = 0; // C1222
  elasticConsts[ 8] = mu2; // C1212

  PetscLogFlops(2);
} // _calcElasticConstsElastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties
// as a viscoelastic material.
void
pylith::materials::MaxwellPlaneStrain::_calcElasticConstsViscoelastic(
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
  assert(_MaxwellPlaneStrain::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellPlaneStrain::tensorSize == initialStrainSize);
 
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];
  const double maxwellTime = properties[p_maxwellTime];

  const double mu2 = 2.0 * mu;
  const double bulkModulus = lambda + mu2 / 3.0;

  double dq = ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTime);

  const double visFac = mu * dq / 3.0;
  elasticConsts[ 0] = bulkModulus + 4.0 * visFac; // C1111
  elasticConsts[ 1] = bulkModulus - 2.0 * visFac; // C1122
  elasticConsts[ 2] = 0; // C1112
  elasticConsts[ 3] = elasticConsts[1]; // C2211
  elasticConsts[ 4] = elasticConsts[0]; // C2222
  elasticConsts[ 5] = 0; // C2212
  elasticConsts[ 6] = 0; // C1211
  elasticConsts[ 7] = 0; // C1222
  elasticConsts[ 8] = 6.0 * visFac; // C1212

  PetscLogFlops(10);
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Update state variables as an elastic material.
void
pylith::materials::MaxwellPlaneStrain::_updateStateVarsElastic(
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
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellPlaneStrain::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;
  const double maxwellTime = properties[p_maxwellTime];

  const double strainTpdt[] = {totalStrain[0] - initialStrain[0],
			       totalStrain[1] - initialStrain[1],
			       0.0,
			       totalStrain[2] - initialStrain[2]};

  const double meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0 };

  stateVars[s_totalStrain] = totalStrain[0];
  stateVars[s_totalStrain + 1] = totalStrain[1];
  stateVars[s_totalStrain + 2] = totalStrain[2];

  for (int iComp=0; iComp < 4; ++iComp) {
    stateVars[s_viscousStrain + iComp] =
      strainTpdt[iComp] - diag[iComp] * meanStrainTpdt;
  } // for
  PetscLogFlops(13);

  _needNewJacobian = true;
} // _updateStateVarsElastic

// ----------------------------------------------------------------------
// Update state variables as a viscoelastic material.
void
pylith::materials::MaxwellPlaneStrain::_updateStateVarsViscoelastic(
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
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellPlaneStrain::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;

  _computeStateVars(stateVars, numStateVars,
		    properties, numProperties,
		    totalStrain, strainSize,
		    initialStress, initialStressSize,
		    initialStrain, initialStrainSize);

  memcpy(&stateVars[s_totalStrain], totalStrain, tensorSize * sizeof(double));

  memcpy(&stateVars[s_viscousStrain], &_viscousStrain[0], 4 * sizeof(double));

  _needNewJacobian = false;

} // _updateStateVarsViscoelastic

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::MaxwellPlaneStrain::_stableTimeStepImplicit(
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
  const double dtStable = 0.2*maxwellTime;

  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Compute viscous strain for current time step.
void
pylith::materials::MaxwellPlaneStrain::_computeStateVars(
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
  assert(_MaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellPlaneStrain::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;
  const double maxwellTime = properties[p_maxwellTime];

  const double strainTpdt[] = {totalStrain[0],
			       totalStrain[1],
			       0.0,
			       totalStrain[2]};
  const double strainT[] = {stateVars[s_totalStrain+0],
			    stateVars[s_totalStrain+1],
			    0.0,
			    stateVars[s_totalStrain+2]};
  
  const double meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;
  const double meanStrainT = (strainT[0] + strainT[1])/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0 };

  // Time integration.
  double dq = ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTime);
  const double expFac = exp(-_dt/maxwellTime);

  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;

  for (int iComp=0; iComp < 4; ++iComp) {
    devStrainTpdt = strainTpdt[iComp] - diag[iComp] * meanStrainTpdt;
    devStrainT = strainT[iComp] - diag[iComp] * meanStrainT;
    _viscousStrain[iComp] = expFac * stateVars[s_viscousStrain+iComp] +
      dq * (devStrainTpdt - devStrainT);
  } // for

  PetscLogFlops(39);
} // _computeStateVars


// End of file 
