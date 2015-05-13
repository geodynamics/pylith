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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "PowerLawPlaneStrain.hh" // implementation of object methods

#include "Metadata.hh" // USES Metadata
#include "EffectiveStress.hh" // USES EffectiveStress

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cmath> // USES fabs()
#include <cassert> // USES assert()
#include <cstring> // USES memcpy()
#include <sstream> // USES std::ostringstream
#include <iostream> // USES std::cout
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _PowerLawPlaneStrain{

      /// Dimension of material.
      const int dimension = 2;

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 3;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 9;

      /// Number of physical properties.
      const int numProperties = 6;

      /// Physical properties.
      const Metadata::ParamDescription properties[6] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
	{ "reference_strain_rate", 1, pylith::topology::FieldBase::SCALAR },
	{ "reference_stress", 1, pylith::topology::FieldBase::SCALAR },
	{ "power_law_exponent", 1, pylith::topology::FieldBase::SCALAR }
      };

      // Values expected in properties spatial database
      const int numDBProperties = 6;
      const char* dbProperties[6] = {"density", "vs", "vp" ,
				     "reference-strain-rate",
				     "reference-stress",
				     "power-law-exponent"};

      /// Number of state variables.
      const int numStateVars = 3;

      /// State variables.
      const Metadata::ParamDescription stateVars[3] = {
	{ "stress_zz_initial", 1, pylith::topology::FieldBase::SCALAR },
	{ "viscous_strain", 4, pylith::topology::FieldBase::OTHER },
	{ "stress4", 4, pylith::topology::FieldBase::OTHER }
      };

      // Values expected in state variables spatial database.
      const int numDBStateVars = 9;
      const char* dbStateVars[9] = { "stress-zz-initial",
				     "viscous-strain-xx",
				     "viscous-strain-yy",
				     "viscous-strain-zz",
				     "viscous-strain-xy",
				     "stress4-xx",
				     "stress4-yy",
				     "stress4-zz",
				     "stress4-xy"
      };

    } // _PowerLawPlaneStrain
  } // materials
} // pylith

// Indices of physical properties.
const int pylith::materials::PowerLawPlaneStrain::p_density = 0;

const int pylith::materials::PowerLawPlaneStrain::p_mu = 
  pylith::materials::PowerLawPlaneStrain::p_density + 1;

const int pylith::materials::PowerLawPlaneStrain::p_lambda = 
  pylith::materials::PowerLawPlaneStrain::p_mu + 1;

const int pylith::materials::PowerLawPlaneStrain::p_referenceStrainRate = 
  pylith::materials::PowerLawPlaneStrain::p_lambda + 1;

const int pylith::materials::PowerLawPlaneStrain::p_referenceStress = 
  pylith::materials::PowerLawPlaneStrain::p_referenceStrainRate + 1;

const int pylith::materials::PowerLawPlaneStrain::p_powerLawExponent = 
  pylith::materials::PowerLawPlaneStrain::p_referenceStress + 1;

// Indices of property database values (order must match dbProperties).
const int pylith::materials::PowerLawPlaneStrain::db_density = 0;

const int pylith::materials::PowerLawPlaneStrain::db_vs = 
  pylith::materials::PowerLawPlaneStrain::db_density + 1;

const int pylith::materials::PowerLawPlaneStrain::db_vp = 
  pylith::materials::PowerLawPlaneStrain::db_vs + 1;

const int pylith::materials::PowerLawPlaneStrain::db_referenceStrainRate = 
  pylith::materials::PowerLawPlaneStrain::db_vp + 1;

const int pylith::materials::PowerLawPlaneStrain::db_referenceStress = 
  pylith::materials::PowerLawPlaneStrain::db_referenceStrainRate + 1;

const int pylith::materials::PowerLawPlaneStrain::db_powerLawExponent = 
  pylith::materials::PowerLawPlaneStrain::db_referenceStress + 1;

// Indices of state variables.
const int pylith::materials::PowerLawPlaneStrain::s_stressZZInitial = 0;

const int pylith::materials::PowerLawPlaneStrain::s_viscousStrain =
  s_stressZZInitial + 1;

const int pylith::materials::PowerLawPlaneStrain::s_stress4 = 
  pylith::materials::PowerLawPlaneStrain::s_viscousStrain + 4;

// Indices of state variable database values (order must match dbStateVars).
const int pylith::materials::PowerLawPlaneStrain::db_stressZZInitial = 0;

const int pylith::materials::PowerLawPlaneStrain::db_viscousStrain =
  db_stressZZInitial + 1;

const int pylith::materials::PowerLawPlaneStrain::db_stress4 = 
  pylith::materials::PowerLawPlaneStrain::db_viscousStrain + 4;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::PowerLawPlaneStrain::PowerLawPlaneStrain(void) :
  ElasticMaterial(_PowerLawPlaneStrain::dimension,
		  _PowerLawPlaneStrain::tensorSize,
		  _PowerLawPlaneStrain::numElasticConsts,
		  Metadata(_PowerLawPlaneStrain::properties,
			   _PowerLawPlaneStrain::numProperties,
			   _PowerLawPlaneStrain::dbProperties,
			   _PowerLawPlaneStrain::numDBProperties,
			   _PowerLawPlaneStrain::stateVars,
			   _PowerLawPlaneStrain::numStateVars,
			   _PowerLawPlaneStrain::dbStateVars,
			   _PowerLawPlaneStrain::numDBStateVars)),
  _calcElasticConstsFn(0),
  _calcStressFn(0),
  _updateStateVarsFn(0)
{ // constructor
  useElasticBehavior(false);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::PowerLawPlaneStrain::~PowerLawPlaneStrain(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set whether elastic or inelastic constitutive relations are used.
void
pylith::materials::PowerLawPlaneStrain::useElasticBehavior(const bool flag)
{ // useElasticBehavior
  if (flag) {
    _calcStressFn = 
      &pylith::materials::PowerLawPlaneStrain::_calcStressElastic;
    _calcElasticConstsFn = 
      &pylith::materials::PowerLawPlaneStrain::_calcElasticConstsElastic;
    _updateStateVarsFn = 
      &pylith::materials::PowerLawPlaneStrain::_updateStateVarsElastic;

  } else {
    _calcStressFn = 
      &pylith::materials::PowerLawPlaneStrain::_calcStressViscoelastic;
    _calcElasticConstsFn = 
      &pylith::materials::PowerLawPlaneStrain::_calcElasticConstsViscoelastic;
    _updateStateVarsFn = 
      &pylith::materials::PowerLawPlaneStrain::_updateStateVarsViscoelastic;
  } // if/else
} // useElasticBehavior

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::PowerLawPlaneStrain::_dbToProperties(
				PylithScalar* const propValues,
				const scalar_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_PowerLawPlaneStrain::numDBProperties == numDBValues);

  const PylithScalar density = dbValues[db_density];
  const PylithScalar vs = dbValues[db_vs];
  const PylithScalar vp = dbValues[db_vp];
  const PylithScalar referenceStrainRate = dbValues[db_referenceStrainRate];
  const PylithScalar referenceStress = dbValues[db_referenceStress];
  const PylithScalar powerLawExponent = dbValues[db_powerLawExponent];
 
  if (density <= 0.0 || vs <= 0.0 || vp <= 0.0 || referenceStrainRate <= 0.0
      || referenceStress <= 0.0 || powerLawExponent < 1.0) {
    std::ostringstream msg;
    msg << "Spatial database returned illegal value for physical "
	<< "properties.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n"
	<< "referenceStrainRate: " << referenceStrainRate << "\n"
	<< "referenceStress: " << referenceStress << "\n"
	<< "powerLawExponent: " << powerLawExponent << "\n";
    throw std::runtime_error(msg.str());
  } // if

  const PylithScalar mu = density * vs * vs;
  const PylithScalar lambda = density * vp * vp - 2.0 * mu;

  if (lambda <= 0.0) {
    std::ostringstream msg;
    msg << "Attempted to set Lame's constant lambda to nonpositive value.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if
  assert(mu > 0);

  propValues[p_density] = density;
  propValues[p_mu] = mu;
  propValues[p_lambda] = lambda;
  propValues[p_referenceStrainRate] = referenceStrainRate;
  propValues[p_referenceStress] = referenceStress;
  propValues[p_powerLawExponent] = powerLawExponent;

  PetscLogFlops(6);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::PowerLawPlaneStrain::_nondimProperties(PylithScalar* const values,
					         const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numPropsQuadPt);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar strainRateScale = 1.0/timeScale;

  values[p_density] = 
    _normalizer->nondimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->nondimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->nondimensionalize(values[p_lambda], pressureScale);
  values[p_referenceStrainRate] = 
    _normalizer->nondimensionalize(values[p_referenceStrainRate],
				   strainRateScale);
  values[p_referenceStress] = 
    _normalizer->nondimensionalize(values[p_referenceStress], pressureScale);

  PetscLogFlops(6);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::PowerLawPlaneStrain::_dimProperties(PylithScalar* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numPropsQuadPt);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar strainRateScale = 1.0/timeScale;

  values[p_density] = 
    _normalizer->dimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->dimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->dimensionalize(values[p_lambda], pressureScale);
  values[p_referenceStrainRate] = 
    _normalizer->dimensionalize(values[p_referenceStrainRate], strainRateScale);
  values[p_referenceStress] = 
    _normalizer->dimensionalize(values[p_referenceStress], pressureScale);

  PetscLogFlops(6);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute initial state variables from values in spatial database.
void
pylith::materials::PowerLawPlaneStrain::_dbToStateVars(
				PylithScalar* const stateValues,
				const scalar_array& dbValues)
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_PowerLawPlaneStrain::numDBStateVars == numDBValues);

  const int totalSize = 1  + 2 * 4;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  memcpy(stateValues, &dbValues[0], totalSize*sizeof(PylithScalar));

  PetscLogFlops(0);
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::materials::PowerLawPlaneStrain::_nondimStateVars(PylithScalar* const values,
						const int nvalues) const
{ // _nondimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numVarsQuadPt);

  const PylithScalar pressureScale = _normalizer->pressureScale();
  _normalizer->nondimensionalize(&values[s_stress4], 4, pressureScale);
  _normalizer->nondimensionalize(&values[s_stressZZInitial], 1, pressureScale);

  PetscLogFlops(5);
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::materials::PowerLawPlaneStrain::_dimStateVars(PylithScalar* const values,
					     const int nvalues) const
{ // _dimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numVarsQuadPt);

  const PylithScalar pressureScale = _normalizer->pressureScale();
  _normalizer->dimensionalize(&values[s_stress4], 4, pressureScale);
  _normalizer->dimensionalize(&values[s_stressZZInitial], 1, pressureScale);

  PetscLogFlops(_tensorSize);
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::PowerLawPlaneStrain::_calcDensity(PylithScalar* const density,
					    const PylithScalar* properties,
					    const int numProperties,
					    const PylithScalar* stateVars,
					    const int numStateVars)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);

  density[0] = properties[p_density];
} // _calcDensity

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::PowerLawPlaneStrain::_stableTimeStepImplicit(
				  const PylithScalar* properties,
				  const int numProperties,
				  const PylithScalar* stateVars,
				  const int numStateVars) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);

  const int tensorSizePS = 4;

  const PylithScalar mu = properties[p_mu];
  const PylithScalar referenceStrainRate = properties[p_referenceStrainRate];
  const PylithScalar referenceStress = properties[p_referenceStress];
  const PylithScalar powerLawExp = properties[p_powerLawExponent];

  const PylithScalar stress4[tensorSizePS] = {stateVars[s_stress4],
					      stateVars[s_stress4 + 1],
					      stateVars[s_stress4 + 2],
					      stateVars[s_stress4 + 3]};
  const PylithScalar meanStress = (stress4[0] + stress4[1] + stress4[2])/3.0;
  const PylithScalar devStress[tensorSizePS] = {stress4[0] - meanStress,
						stress4[1] - meanStress,
						stress4[2] - meanStress,
						stress4[3]};
  const PylithScalar devStressProd = scalarProduct2DPS(devStress, devStress);
  const PylithScalar effStress = sqrt(0.5 * devStressProd);
  PylithScalar dtTest = 0.0;
  if (effStress <= 0.0) {
    dtTest = pylith::PYLITH_MAXSCALAR;
  } else {
    dtTest = pow((referenceStress/effStress), (powerLawExp - 1.0)) *
      (referenceStress/mu)/(referenceStrainRate * 30.0);
  } //else
  const PylithScalar dtStable = dtTest;

#if 0 // DEBUGGING
  PylithScalar maxwellTime = 10.0 * dtStable;
  std::cout << "Maxwell time:  " << maxwellTime << std::endl;
#endif
  PetscLogFlops(23);
  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Get stable time step for explicit time integration.
PylithScalar
pylith::materials::PowerLawPlaneStrain::_stableTimeStepExplicit(const PylithScalar* properties,
								const int numProperties,
								const PylithScalar* stateVars,
								const int numStateVars,
								const double minCellWidth) const
{ // _stableTimeStepExplicit
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
 
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar density = properties[p_density];

  assert(density > 0.0);
  const PylithScalar vp = sqrt((lambda + 2*mu) / density);

  const PylithScalar dtStable = minCellWidth / vp;
  return dtStable;
} // _stableTimeStepExplicit


// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as an elastic
// material.
void
pylith::materials::PowerLawPlaneStrain::_calcStressElastic(
				         PylithScalar* const stress,
					 const int stressSize,
					 const PylithScalar* properties,
					 const int numProperties,
					 const PylithScalar* stateVars,
					 const int numStateVars,
					 const PylithScalar* totalStrain,
					 const int strainSize,
					 const PylithScalar* initialStress,
					 const int initialStressSize,
					 const PylithScalar* initialStrain,
					 const int initialStrainSize,
					 const bool computeStateVars)
{ // _calcStressElastic
  assert(0 != stress);
  assert(_PowerLawPlaneStrain::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLawPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLawPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLawPlaneStrain::tensorSize == initialStrainSize);

  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar mu2 = 2.0 * mu;

  const PylithScalar e11 = totalStrain[0] - initialStrain[0];
  const PylithScalar e22 = totalStrain[1] - initialStrain[1];
  const PylithScalar e12 = totalStrain[2] - initialStrain[2];
  
  const PylithScalar s12 = lambda * (e11 + e22);

  stress[0] = s12 + mu2 * e11 + initialStress[0];
  stress[1] = s12 + mu2 * e22 + initialStress[1];
  stress[2] = mu2 * e12 + initialStress[2];

  PetscLogFlops(14);
} // _calcStressElastic

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as a viscoelastic
// material.
void
pylith::materials::PowerLawPlaneStrain::_calcStressViscoelastic(
					PylithScalar* const stress,
					const int stressSize,
					const PylithScalar* properties,
					const int numProperties,
					const PylithScalar* stateVars,
					const int numStateVars,
					const PylithScalar* totalStrain,
					const int strainSize,
					const PylithScalar* initialStress,
					const int initialStressSize,
					const PylithScalar* initialStrain,
					const int initialStrainSize,
					const bool computeStateVars)
{ // _calcStressViscoelastic
  assert(0 != stress);
  assert(_PowerLawPlaneStrain::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLawPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLawPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLawPlaneStrain::tensorSize == initialStrainSize);

  const int tensorSize = 3;
  const int tensorSizePS = 4;
  assert(_tensorSize == tensorSize);
    
  // We need to do root-finding method if state variables are from previous
  // time step.
  if (computeStateVars) {

    const PylithScalar mu = properties[p_mu];
    const PylithScalar lambda = properties[p_lambda];
    const PylithScalar referenceStrainRate = properties[p_referenceStrainRate];
    const PylithScalar referenceStress = properties[p_referenceStress];
    const PylithScalar powerLawExp = properties[p_powerLawExponent];

    const PylithScalar stressZZInitial = stateVars[s_stressZZInitial];
    const PylithScalar visStrainT[tensorSizePS] = {
      stateVars[s_viscousStrain],
      stateVars[s_viscousStrain + 1],
      stateVars[s_viscousStrain + 2],
      stateVars[s_viscousStrain + 3]};
    const PylithScalar stressT[tensorSizePS] = {
      stateVars[s_stress4],
      stateVars[s_stress4 + 1],
      stateVars[s_stress4 + 2],
      stateVars[s_stress4 + 3]};

    const PylithScalar mu2 = 2.0 * mu;
    const PylithScalar bulkModulus = lambda + mu2/3.0;
    const PylithScalar ae = 1.0/mu2;
    const PylithScalar diagg[tensorSizePS] = { 1.0, 1.0, 1.0, 0.0 };

    // Need to figure out how time integration parameter alpha is going to be
    // specified.  It should probably be specified in the problem definition and
    // then used only by the material types that use it.  For now we are setting
    // it to 0.5, which should probably be the default value.
    const PylithScalar alpha = 0.5;
    const PylithScalar timeFac = _dt * (1.0 - alpha);

    // Initial stress values
    const PylithScalar meanStressInitial = (initialStress[0] +
					    initialStress[1] +
					    stressZZInitial)/3.0;
    const PylithScalar devStressInitial[tensorSizePS] = {
      initialStress[0] - meanStressInitial,
      initialStress[1] - meanStressInitial,
      stressZZInitial - meanStressInitial,
      initialStress[2]
    };
    const PylithScalar stressInvar2Initial =
      0.5 * scalarProduct2DPS(devStressInitial, devStressInitial);

    // Initial strain values
    const PylithScalar meanStrainInitial = (initialStrain[0] +
					    initialStrain[1])/3.0;

    // Values for current time step
    const PylithScalar meanStrainTpdt = (totalStrain[0] + totalStrain[1])/3.0 -
      meanStrainInitial;
    const PylithScalar meanStressTpdt = 3.0 * bulkModulus * meanStrainTpdt;

    // Note that I use the initial strain rather than the deviatoric initial
    // strain since otherwise the initial mean strain would get used twice.
    const PylithScalar strainPPTpdt[tensorSizePS] =
      { totalStrain[0] - meanStrainTpdt - visStrainT[0] - initialStrain[0],
	totalStrain[1] - meanStrainTpdt - visStrainT[1] - initialStrain[1],
	- meanStrainTpdt - visStrainT[2],
	totalStrain[2] - visStrainT[3] - initialStrain[2]
      };
    const PylithScalar strainPPInvar2Tpdt =
      0.5 * scalarProduct2DPS(strainPPTpdt, strainPPTpdt);

    // Values for previous time step
    const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const PylithScalar devStressT[tensorSizePS] = {
      stressT[0] - meanStressT,
      stressT[1] - meanStressT,
      stressT[2] - meanStressT,
      stressT[3]
    };
    const PylithScalar stressInvar2T =
      0.5 * scalarProduct2DPS(devStressT, devStressT);
    const PylithScalar effStressT = sqrt(stressInvar2T);

    // Finish defining parameters needed for root-finding algorithm.
    const PylithScalar b = strainPPInvar2Tpdt + ae *
      scalarProduct2DPS(strainPPTpdt, devStressInitial) +
      ae * ae * stressInvar2Initial;
    const PylithScalar c =
      (scalarProduct2DPS(strainPPTpdt, devStressT) +
       ae * scalarProduct2DPS(devStressT, devStressInitial)) * timeFac;
    const PylithScalar d = timeFac * effStressT;

    PetscLogFlops(94);

    // If b, c, and d are all zero, then the effective stress is zero and we
    // don't need a root-finding algorithm. Otherwise, use the algorithm to
    // find the effective stress.
    PylithScalar effStressTpdt = 0.0;
    if (b != 0.0 || c != 0.0 || d != 0.0) {
      const PylithScalar stressScale = mu;

      // Put parameters into a struct and call root-finding algorithm.
      _effStressParams.ae = ae;
      _effStressParams.b = b;
      _effStressParams.c = c;
      _effStressParams.d = d;
      _effStressParams.alpha = alpha;
      _effStressParams.dt = _dt;
      _effStressParams.effStressT = effStressT;
      _effStressParams.powerLawExp = powerLawExp;
      _effStressParams.referenceStrainRate = referenceStrainRate;
      _effStressParams.referenceStress = referenceStress;
      
      const PylithScalar effStressInitialGuess = effStressT;

      effStressTpdt =
	EffectiveStress::calculate<PowerLawPlaneStrain>(effStressInitialGuess,
					       stressScale, this);
    } // if

    // Compute stresses from effective stress.
    const PylithScalar effStressTau = (1.0 - alpha) * effStressT +
      alpha * effStressTpdt;
    const PylithScalar gammaTau = referenceStrainRate *
      pow((effStressTau/referenceStress),
	  (powerLawExp - 1.0))/referenceStress;
    const PylithScalar factor1 = 1.0/(ae + alpha * _dt * gammaTau);
    const PylithScalar factor2 = timeFac * gammaTau;
    PylithScalar devStressTpdt = 0.0;
    PylithScalar totalStress[tensorSizePS];

    for (int iComp=0; iComp < tensorSizePS; ++iComp) {
      devStressTpdt = factor1 *
	(strainPPTpdt[iComp] - factor2 * devStressT[iComp] +
	 ae * devStressInitial[iComp]);
      totalStress[iComp] = devStressTpdt + diagg[iComp] *
	(meanStressTpdt + meanStressInitial);
    } // for
    stress[0] = totalStress[0];
    stress[1] = totalStress[1];
    stress[2] = totalStress[3];
    PetscLogFlops(14 + 8 * tensorSizePS);

    // If state variables have already been updated, current stress is already
    // contained in stress.
  } else {
    stress[0] = stateVars[s_stress4];
    stress[1] = stateVars[s_stress4 + 1];
    stress[2] = stateVars[s_stress4 + 3];
  } // else

} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Effective stress function that computes effective stress function only
// (no derivative).
PylithScalar
pylith::materials::PowerLawPlaneStrain::effStressFunc(const PylithScalar effStressTpdt)
{ // effStressFunc
  const PylithScalar ae = _effStressParams.ae;
  const PylithScalar b = _effStressParams.b;
  const PylithScalar c = _effStressParams.c;
  const PylithScalar d = _effStressParams.d;
  const PylithScalar alpha = _effStressParams.alpha;
  const PylithScalar dt = _effStressParams.dt;
  const PylithScalar effStressT = _effStressParams.effStressT;
  const PylithScalar powerLawExp = _effStressParams.powerLawExp;
  const PylithScalar referenceStrainRate = _effStressParams.referenceStrainRate;
  const PylithScalar referenceStress = _effStressParams.referenceStress;
  const PylithScalar factor1 = 1.0-alpha;
  const PylithScalar effStressTau = factor1 * effStressT +
    alpha * effStressTpdt;
  const PylithScalar gammaTau = referenceStrainRate * 
    pow((effStressTau/referenceStress), (powerLawExp - 1.0))/referenceStress;
  const PylithScalar a = ae + alpha * dt * gammaTau;
  const PylithScalar y = a * a * effStressTpdt * effStressTpdt - b +
    c * gammaTau - d * d * gammaTau * gammaTau;

  PetscLogFlops(22);

  return y;
} // effStressFunc

// ----------------------------------------------------------------------
// Effective stress function that computes effective stress function
// derivative only (no function value).
PylithScalar
pylith::materials::PowerLawPlaneStrain::effStressDerivFunc(const PylithScalar effStressTpdt)
{ // effStressDFunc
  const PylithScalar ae = _effStressParams.ae;
  const PylithScalar c = _effStressParams.c;
  const PylithScalar d = _effStressParams.d;
  const PylithScalar alpha = _effStressParams.alpha;
  const PylithScalar dt = _effStressParams.dt;
  const PylithScalar effStressT = _effStressParams.effStressT;
  const PylithScalar powerLawExp = _effStressParams.powerLawExp;
  const PylithScalar referenceStrainRate = _effStressParams.referenceStrainRate;
  const PylithScalar referenceStress = _effStressParams.referenceStress;
  const PylithScalar factor1 = 1.0-alpha;
  const PylithScalar effStressTau = factor1 * effStressT +
    alpha * effStressTpdt;
  const PylithScalar gammaTau = referenceStrainRate *
    pow((effStressTau/referenceStress), (powerLawExp - 1.0))/referenceStress;
  const PylithScalar a = ae + alpha * dt * gammaTau;
  const PylithScalar dGammaTau =
    referenceStrainRate * alpha * (powerLawExp - 1.0) *
    pow((effStressTau/referenceStress), (powerLawExp - 2.0))/
    (referenceStress * referenceStress);
  const PylithScalar dy = 2.0 * a * a * effStressTpdt + dGammaTau *
    (2.0 * a * alpha * dt * effStressTpdt * effStressTpdt +
     c - 2.0 * d * d * gammaTau);
  PetscLogFlops(36);

  return dy;
} // effStressDFunc

// ----------------------------------------------------------------------
// Effective stress function that computes effective stress function
// and derivative.
void
pylith::materials::PowerLawPlaneStrain::effStressFuncDerivFunc( PylithScalar* func,
								PylithScalar* dfunc,
								const PylithScalar effStressTpdt)
{ // effStressFuncDFunc
  PylithScalar y = *func;
  PylithScalar dy = *dfunc;

  const PylithScalar ae = _effStressParams.ae;
  const PylithScalar b = _effStressParams.b;
  const PylithScalar c = _effStressParams.c;
  const PylithScalar d = _effStressParams.d;
  const PylithScalar alpha = _effStressParams.alpha;
  const PylithScalar dt = _effStressParams.dt;
  const PylithScalar effStressT = _effStressParams.effStressT;
  const PylithScalar powerLawExp = _effStressParams.powerLawExp;
  const PylithScalar referenceStrainRate = _effStressParams.referenceStrainRate;
  const PylithScalar referenceStress = _effStressParams.referenceStress;
  const PylithScalar factor1 = 1.0-alpha;
  const PylithScalar effStressTau = factor1 * effStressT + alpha *
    effStressTpdt;
  const PylithScalar gammaTau = referenceStrainRate *
    pow((effStressTau/referenceStress), (powerLawExp - 1.0))/referenceStress;
  const PylithScalar dGammaTau =
    referenceStrainRate * alpha * (powerLawExp - 1.0) *
    pow((effStressTau/referenceStress), (powerLawExp - 2.0))/
    (referenceStress * referenceStress);
  const PylithScalar a = ae + alpha * dt * gammaTau;
  y = a * a * effStressTpdt * effStressTpdt -
    b +
    c * gammaTau -
    d * d * gammaTau * gammaTau;
  dy = 2.0 * a * a * effStressTpdt +
    dGammaTau *
    (2.0 * a * alpha * dt * effStressTpdt * effStressTpdt +
     c - 2.0 * d * d * gammaTau);
  
  *func = y;
  *dfunc = dy;

  PetscLogFlops(46);
} // effStressFuncDFunc

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::PowerLawPlaneStrain::_calcElasticConstsElastic(
				         PylithScalar* const elasticConsts,
					 const int numElasticConsts,
					 const PylithScalar* properties,
					 const int numProperties,
					 const PylithScalar* stateVars,
					 const int numStateVars,
					 const PylithScalar* totalStrain,
					 const int strainSize,
					 const PylithScalar* initialStress,
					 const int initialStressSize,
					 const PylithScalar* initialStrain,
					 const int initialStrainSize)
{ // _calcElasticConstsElastic
  assert(0 != elasticConsts);
  assert(_PowerLawPlaneStrain::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLawPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLawPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLawPlaneStrain::tensorSize == initialStrainSize);
 
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar lambda2mu = lambda + mu2;

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
pylith::materials::PowerLawPlaneStrain::_calcElasticConstsViscoelastic(
				         PylithScalar* const elasticConsts,
					 const int numElasticConsts,
					 const PylithScalar* properties,
					 const int numProperties,
					 const PylithScalar* stateVars,
					 const int numStateVars,
					 const PylithScalar* totalStrain,
					 const int strainSize,
					 const PylithScalar* initialStress,
					 const int initialStressSize,
					 const PylithScalar* initialStrain,
					 const int initialStrainSize)
{ // _calcElasticConstsViscoelastic
  assert(0 != elasticConsts);
  assert(_PowerLawPlaneStrain::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLawPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLawPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLawPlaneStrain::tensorSize == initialStrainSize);

  const int tensorSize = 3;
  const int tensorSizePS = 4;

  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar referenceStrainRate = properties[p_referenceStrainRate];
  const PylithScalar referenceStress = properties[p_referenceStress];
  const PylithScalar powerLawExp = properties[p_powerLawExponent];
    
  // State variables.
  const PylithScalar stressZZInitial = stateVars[s_stressZZInitial];
  const PylithScalar visStrainT[tensorSizePS] = {
    stateVars[s_viscousStrain],
    stateVars[s_viscousStrain + 1],
    stateVars[s_viscousStrain + 2],
    stateVars[s_viscousStrain + 3]
  };
  const PylithScalar stressT[tensorSizePS] = {stateVars[s_stress4],
					      stateVars[s_stress4 + 1],
					      stateVars[s_stress4 + 2],
					      stateVars[s_stress4 + 3]};

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2/3.0;
  const PylithScalar ae = 1.0/mu2;
    
  // Need to figure out how time integration parameter alpha is going to be
  // specified.  It should probably be specified in the problem definition and
  // then used only by the material types that use it.  For now we are setting
  // it to 0.5, which should probably be the default value.
  const PylithScalar alpha = 0.5;
  const PylithScalar explicitFac = 1.0 - alpha;
  const PylithScalar timeFac = _dt * explicitFac;
    
  /// Initial state.
  // Initial stress values.
  const PylithScalar meanStressInitial = (initialStress[0] +
					  initialStress[1] +
					  stressZZInitial)/3.0;
  const PylithScalar devStressInitial[tensorSizePS] = {
    initialStress[0] - meanStressInitial,
    initialStress[1] - meanStressInitial,
    stressZZInitial - meanStressInitial,
    initialStress[2]
  };
  const PylithScalar stressInvar2Initial =
    0.5 * scalarProduct2DPS(devStressInitial, devStressInitial);

  // Initial strain values.
  const PylithScalar meanStrainInitial = (initialStrain[0] +
					  initialStrain[1])/3.0;
  
  /// Values for current time step
  const PylithScalar meanStrainTpdt = (totalStrain[0] + totalStrain[1])/3.0 -
    meanStrainInitial;
  
  // Note that I use the initial strain rather than the deviatoric initial
  // strain since otherwise the initial mean strain would get used twice.
  
  const PylithScalar strainPPTpdt[tensorSizePS] =
    { totalStrain[0] - meanStrainTpdt - visStrainT[0] - initialStrain[0],
      totalStrain[1] - meanStrainTpdt - visStrainT[1] - initialStrain[1],
      - meanStrainTpdt - visStrainT[2],
      totalStrain[2] - visStrainT[3] - initialStrain[2]
    };
  const PylithScalar strainPPInvar2Tpdt =
    0.5 * scalarProduct2DPS(strainPPTpdt, strainPPTpdt);
  
  // Values for previous time step
  const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
  const PylithScalar devStressT[tensorSizePS] = { stressT[0] - meanStressT,
						  stressT[1] - meanStressT,
						  stressT[2] - meanStressT,
						  stressT[3] };
  const PylithScalar stressInvar2T =
    0.5 * scalarProduct2DPS(devStressT, devStressT);
  const PylithScalar effStressT = sqrt(stressInvar2T);
    
  // Finish defining parameters needed for root-finding algorithm.
  const PylithScalar b = strainPPInvar2Tpdt +
    ae * scalarProduct2DPS(strainPPTpdt, devStressInitial) +
    ae * ae * stressInvar2Initial;
  const PylithScalar c =
    (scalarProduct2DPS(strainPPTpdt, devStressT) +
     ae * scalarProduct2DPS(devStressT, devStressInitial)) * timeFac;
  const PylithScalar d = timeFac * effStressT;

  PetscLogFlops(96);

  // If b = c = d = 0, the effective stress is zero and the elastic constants
  // will be the same as for the elastic case. Otherwise, compute the tangent
  // matrix using the effective stress function algorithm.
  if (b == 0.0 && c == 0.0 && d == 0.0) {
    _calcElasticConstsElastic(elasticConsts,
			      numElasticConsts,
			      properties,
			      numProperties,
			      stateVars,
			      numStateVars,
			      totalStrain,
			      strainSize,
			      initialStress,
			      initialStressSize,
			      initialStrain,
			      initialStrainSize);
  } else {
    const PylithScalar stressScale = mu;
  
    // Put parameters into a struct and call root-finding algorithm.
    _effStressParams.ae = ae;
    _effStressParams.b = b;
    _effStressParams.c = c;
    _effStressParams.d = d;
    _effStressParams.alpha = alpha;
    _effStressParams.dt = _dt;
    _effStressParams.effStressT = effStressT;
    _effStressParams.powerLawExp = powerLawExp;
    _effStressParams.referenceStrainRate = referenceStrainRate;
    _effStressParams.referenceStress = referenceStress;
    
    const PylithScalar effStressInitialGuess = effStressT;
    
    const PylithScalar effStressTpdt =
      EffectiveStress::calculate<PowerLawPlaneStrain>(effStressInitialGuess,
					     stressScale, this);
  
    // Compute quantities at intermediate time tau used to compute values at
    // end of time step.
    const PylithScalar effStressTau = (1.0 - alpha) * effStressT +
      alpha * effStressTpdt;
    const PylithScalar gammaTau = referenceStrainRate *
      pow((effStressTau/referenceStress),
	  (powerLawExp - 1.0))/referenceStress;
    const PylithScalar a = ae + alpha * _dt * gammaTau;
    const PylithScalar factor1 = 1.0/a;
    const PylithScalar factor2 = timeFac * gammaTau;
    const PylithScalar devStressTpdt[tensorSize] = {
      factor1 *
      (strainPPTpdt[0] - factor2 * devStressT[0] + ae * devStressInitial[0]),
      factor1 *
      (strainPPTpdt[1] - factor2 * devStressT[1] + ae * devStressInitial[1]),
      factor1 *
      (strainPPTpdt[3] - factor2 * devStressT[3] + ae * devStressInitial[3])
    };
    const PylithScalar devStressTau[tensorSize] = {
      alpha * devStressT[0] + explicitFac * devStressTpdt[0],
      alpha * devStressT[1] + explicitFac * devStressTpdt[1],
      alpha * devStressT[2] + explicitFac * devStressTpdt[2]
    };
    const PylithScalar factor3 = 0.5 * referenceStrainRate * _dt * alpha *
      (powerLawExp - 1.0) *
      pow((effStressTau/referenceStress), (powerLawExp - 2.0))/
      (referenceStress * referenceStress * effStressTpdt);

    // Compute deviatoric derivatives
    const PylithScalar dStress11dStrain11 = 1.0/
      (a + devStressTau[0] * devStressTpdt[0] * factor3);
    const PylithScalar dStress22dStrain22 = 1.0/
      (a + devStressTau[1] * devStressTpdt[1] * factor3);
    const PylithScalar dStress12dStrain12 = 1.0/
      (a + 2.0 * devStressTau[2] * devStressTpdt[2] * factor3);
    
    /// Compute tangent matrix.
    // Form elastic constants.
    elasticConsts[ 0] = bulkModulus + 2.0 * dStress11dStrain11/3.0; // C1111
    elasticConsts[ 1] = bulkModulus -       dStress11dStrain11/3.0; // C1122
    elasticConsts[ 2] = 0.0; // C1112
    elasticConsts[ 3] = elasticConsts[ 1]; // C2211
    elasticConsts[ 4] = bulkModulus + 2.0 * dStress22dStrain22/3.0; // C2222
    elasticConsts[ 5] = 0.0; // C2212
    elasticConsts[ 6] = 0.0; // C1211
    elasticConsts[ 7] = 0.0; // C1222
    elasticConsts[ 8] = dStress12dStrain12; // C1212
    PetscLogFlops(71);
  } // else
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::PowerLawPlaneStrain::_updateStateVarsElastic(
				    PylithScalar* const stateVars,
				    const int numStateVars,
				    const PylithScalar* properties,
				    const int numProperties,
				    const PylithScalar* totalStrain,
				    const int strainSize,
				    const PylithScalar* initialStress,
				    const int initialStressSize,
				    const PylithScalar* initialStrain,
				    const int initialStrainSize)
{ // _updateStateVarsElastic
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_PowerLawPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLawPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLawPlaneStrain::tensorSize == initialStrainSize);

  const int tensorSize = 3;
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar stressZZInitial = stateVars[s_stressZZInitial];

  const bool computeStateVars = true;
  PylithScalar stress[tensorSize] = {0.0, 0.0, 0.0};
  const int stressSize = strainSize;
  _calcStressElastic(stress, stressSize,
		     properties, numProperties,
		     stateVars, numStateVars,
		     totalStrain, strainSize,
		     initialStress, initialStressSize,
		     initialStrain, initialStrainSize,
		     computeStateVars);

  stateVars[s_viscousStrain] = 0.0;
  stateVars[s_viscousStrain + 1] = 0.0;
  stateVars[s_viscousStrain + 2] = 0.0;
  stateVars[s_viscousStrain + 3] = 0.0;
  stateVars[s_stress4] = stress[0];
  stateVars[s_stress4 + 1] = stress[1];
  stateVars[s_stress4 + 2] = lambda * (totalStrain[0] + totalStrain[1]) +
    stressZZInitial;
  stateVars[s_stress4 + 3] = stress[2];

  _needNewJacobian = true;
} // _updateStateVarsElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::PowerLawPlaneStrain::_updateStateVarsViscoelastic(
				    PylithScalar* const stateVars,
				    const int numStateVars,
				    const PylithScalar* properties,
				    const int numProperties,
				    const PylithScalar* totalStrain,
				    const int strainSize,
				    const PylithScalar* initialStress,
				    const int initialStressSize,
				    const PylithScalar* initialStrain,
				    const int initialStrainSize)
{ // _updateStateVarsViscoelastic
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_PowerLawPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLawPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLawPlaneStrain::tensorSize == initialStrainSize);

  const int tensorSizePS = 4;

  // For now, we are duplicating the functionality of _calcStressViscoelastic,
  // since otherwise we would have to redo a lot of calculations.
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar referenceStrainRate = properties[p_referenceStrainRate];
  const PylithScalar referenceStress = properties[p_referenceStress];
  const PylithScalar powerLawExp = properties[p_powerLawExponent];

  const PylithScalar stressZZInitial = stateVars[s_stressZZInitial];
  const PylithScalar visStrainT[tensorSizePS] = {
    stateVars[s_viscousStrain],
    stateVars[s_viscousStrain + 1],
    stateVars[s_viscousStrain + 2],
    stateVars[s_viscousStrain + 3]
  };

  const PylithScalar stressT[tensorSizePS] = {stateVars[s_stress4],
					      stateVars[s_stress4 + 1],
					      stateVars[s_stress4 + 2],
					      stateVars[s_stress4 + 3] };
  
  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2/3.0;
  const PylithScalar ae = 1.0/mu2;
  const PylithScalar diagg[tensorSizePS] = { 1.0, 1.0, 1.0, 0.0 };

  // Need to figure out how time integration parameter alpha is going to be
  // specified.  It should probably be specified in the problem definition and
  // then used only by the material types that use it.  For now we are setting
  // it to 0.5, which should probably be the default value.
  const PylithScalar alpha = 0.5;
  const PylithScalar timeFac = _dt * (1.0 - alpha);

  // Initial stress values
  const PylithScalar meanStressInitial = (initialStress[0] + initialStress[1] +
					  stressZZInitial)/3.0;
  const PylithScalar devStressInitial[tensorSizePS] = {
    initialStress[0] - meanStressInitial,
    initialStress[1] - meanStressInitial,
    stressZZInitial - meanStressInitial,
    initialStress[2]
  };
  const PylithScalar stressInvar2Initial =
    0.5 * scalarProduct2DPS(devStressInitial, devStressInitial);

  // Initial strain values
  const PylithScalar meanStrainInitial = (initialStrain[0] +
					  initialStrain[1])/3.0;
  
  // Values for current time step
  const PylithScalar meanStrainTpdt = (totalStrain[0] + totalStrain[1])/3.0
    - meanStrainInitial;
  const PylithScalar meanStressTpdt = 3.0 * bulkModulus * meanStrainTpdt;

  // Note that I use the initial strain rather than the deviatoric initial
  // strain since otherwise the initial mean strain would get used twice.
  const PylithScalar strainPPTpdt[] = {
    totalStrain[0] - meanStrainTpdt - visStrainT[0] - initialStrain[0],
    totalStrain[1] - meanStrainTpdt - visStrainT[1] - initialStrain[1],
    - meanStrainTpdt - visStrainT[2],
    totalStrain[2] - visStrainT[3] - initialStrain[2]
  };
  const PylithScalar strainPPInvar2Tpdt =
    0.5 * scalarProduct2DPS(strainPPTpdt, strainPPTpdt);

  // Values for previous time step
  const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
  const PylithScalar devStressT[tensorSizePS] = {
    stressT[0] - meanStressT,
    stressT[1] - meanStressT,
    stressT[2] - meanStressT,
    stressT[3]
  };
  const PylithScalar stressInvar2T =
    0.5 * scalarProduct2DPS(devStressT, devStressT);
  const PylithScalar effStressT = sqrt(stressInvar2T);

  // Finish defining parameters needed for root-finding algorithm.
  const PylithScalar b = strainPPInvar2Tpdt +
    ae * scalarProduct2DPS(strainPPTpdt, devStressInitial) +
    ae * ae * stressInvar2Initial;
  const PylithScalar c =
    (scalarProduct2DPS(strainPPTpdt, devStressT) +
     ae * scalarProduct2DPS(devStressT, devStressInitial)) * timeFac;
  const PylithScalar d = timeFac * effStressT;
  PetscLogFlops(96);

  // If b, c, and d are all zero, then the effective stress is zero and we
  // don't need a root-finding algorithm. Otherwise, use the algorithm to
  // find the effective stress.
  PylithScalar effStressTpdt = 0.0;
  if (b != 0.0 || c != 0.0 || d != 0.0) {
    const PylithScalar stressScale = mu;

    // Put parameters into a struct and call root-finding algorithm.
    _effStressParams.ae = ae;
    _effStressParams.b = b;
    _effStressParams.c = c;
    _effStressParams.d = d;
    _effStressParams.alpha = alpha;
    _effStressParams.dt = _dt;
    _effStressParams.effStressT = effStressT;
    _effStressParams.powerLawExp = powerLawExp;
    _effStressParams.referenceStrainRate = referenceStrainRate;
    _effStressParams.referenceStress = referenceStress;

    const PylithScalar effStressInitialGuess = effStressT;

    effStressTpdt =
      EffectiveStress::calculate<PowerLawPlaneStrain>(effStressInitialGuess,
					     stressScale, this);

  } // if

  // Compute stress and viscous strain and update appropriate state variables.
  const PylithScalar effStressTau = (1.0 - alpha) * effStressT +
    alpha * effStressTpdt;
  const PylithScalar gammaTau = referenceStrainRate *
    pow((effStressTau/referenceStress),
	(powerLawExp - 1.0))/referenceStress;
  const PylithScalar factor1 = 1.0/(ae + alpha * _dt * gammaTau);
  const PylithScalar factor2 = timeFac * gammaTau;
  PylithScalar devStressTpdt = 0.0;
  PylithScalar devStressTau = 0.0;

  for (int iComp=0; iComp < tensorSizePS; ++iComp) {
    devStressTpdt = factor1 *
      (strainPPTpdt[iComp] - factor2 * devStressT[iComp] +
       ae * devStressInitial[iComp]);
    stateVars[s_stress4 + iComp] = devStressTpdt + diagg[iComp] *
      (meanStressTpdt + meanStressInitial);
    devStressTau = (1.0 - alpha) * devStressT[iComp] + alpha * devStressTpdt;
    stateVars[s_viscousStrain+iComp] += _dt * gammaTau * devStressTau;
  } // for

  _needNewJacobian = true;
  PetscLogFlops(14 + tensorSizePS * 15);

} // _updateStateVarsViscoelastic

// End of file 
