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

#include "DruckerPragerPlaneStrain.hh" // implementation of object methods

#include "Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES MAXSCALAR

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
    namespace _DruckerPragerPlaneStrain{

      /// Dimension of material.
      const int dimension = 2;

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 3;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 9;

      /// Number of physical properties.
      const int numProperties = 6;

      /// Physical properties.
      const Metadata::ParamDescription properties[numProperties] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
	{ "alpha_yield", 1, pylith::topology::FieldBase::SCALAR },
	{ "beta", 1, pylith::topology::FieldBase::SCALAR },
	{ "alpha_flow", 1, pylith::topology::FieldBase::SCALAR }
      };

      // Values expected in properties spatial database
      const int numDBProperties = 6;
      const char* dbProperties[6] = {
	"density",
	"vs",
	"vp" ,
	"friction-angle",
	"cohesion",
	"dilatation-angle"
      };

      /// Number of state variables.
      const int numStateVars = 2;

      /// State variables.
      const Metadata::ParamDescription stateVars[numStateVars] = {
	{ "stress_zz_initial", 1, pylith::topology::FieldBase::SCALAR },
	{ "plastic_strain", 4, pylith::topology::FieldBase::OTHER }
      };

      // Values expected in state variables spatial database.
      const int numDBStateVars = 5;
      const char* dbStateVars[5] = {
	"stress-zz-initial",
	"plastic-strain-xx",
	"plastic-strain-yy",
	"plastic-strain-zz",
	"plastic-strain-xy"
      };

    } // _DruckerPragerPlaneStrain
  } // materials
} // pylith

// Indices of physical properties.
const int pylith::materials::DruckerPragerPlaneStrain::p_density = 0;

const int pylith::materials::DruckerPragerPlaneStrain::p_mu = 
  pylith::materials::DruckerPragerPlaneStrain::p_density + 1;

const int pylith::materials::DruckerPragerPlaneStrain::p_lambda = 
  pylith::materials::DruckerPragerPlaneStrain::p_mu + 1;

const int pylith::materials::DruckerPragerPlaneStrain::p_alphaYield = 
  pylith::materials::DruckerPragerPlaneStrain::p_lambda + 1;

const int pylith::materials::DruckerPragerPlaneStrain::p_beta = 
  pylith::materials::DruckerPragerPlaneStrain::p_alphaYield + 1;

const int pylith::materials::DruckerPragerPlaneStrain::p_alphaFlow = 
  pylith::materials::DruckerPragerPlaneStrain::p_beta + 1;

// Indices of property database values (order must match dbProperties).
const int pylith::materials::DruckerPragerPlaneStrain::db_density = 0;

const int pylith::materials::DruckerPragerPlaneStrain::db_vs = 
  pylith::materials::DruckerPragerPlaneStrain::db_density + 1;

const int pylith::materials::DruckerPragerPlaneStrain::db_vp = 
  pylith::materials::DruckerPragerPlaneStrain::db_vs + 1;

const int pylith::materials::DruckerPragerPlaneStrain::db_frictionAngle = 
  pylith::materials::DruckerPragerPlaneStrain::db_vp + 1;

const int pylith::materials::DruckerPragerPlaneStrain::db_cohesion = 
  pylith::materials::DruckerPragerPlaneStrain::db_frictionAngle + 1;

const int pylith::materials::DruckerPragerPlaneStrain::db_dilatationAngle = 
  pylith::materials::DruckerPragerPlaneStrain::db_cohesion + 1;

// Indices of state variables.
const int pylith::materials::DruckerPragerPlaneStrain::s_stressZZInitial = 0;

const int pylith::materials::DruckerPragerPlaneStrain::s_plasticStrain = 
  pylith::materials::DruckerPragerPlaneStrain::s_stressZZInitial + 1;

// Indices of state variable database values (order must match dbStateVars).
const int pylith::materials::DruckerPragerPlaneStrain::db_stressZZInitial = 0;

const int pylith::materials::DruckerPragerPlaneStrain::db_plasticStrain = 
  pylith::materials::DruckerPragerPlaneStrain::db_stressZZInitial + 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::DruckerPragerPlaneStrain::DruckerPragerPlaneStrain(void) :
  ElasticMaterial(_DruckerPragerPlaneStrain::dimension,
		  _DruckerPragerPlaneStrain::tensorSize,
		  _DruckerPragerPlaneStrain::numElasticConsts,
		  Metadata(_DruckerPragerPlaneStrain::properties,
			   _DruckerPragerPlaneStrain::numProperties,
			   _DruckerPragerPlaneStrain::dbProperties,
			   _DruckerPragerPlaneStrain::numDBProperties,
			   _DruckerPragerPlaneStrain::stateVars,
			   _DruckerPragerPlaneStrain::numStateVars,
			   _DruckerPragerPlaneStrain::dbStateVars,
			   _DruckerPragerPlaneStrain::numDBStateVars)),
  _calcElasticConstsFn(0),
  _calcStressFn(0),
  _updateStateVarsFn(0),
  _fitMohrCoulomb(MOHR_COULOMB_INSCRIBED),
  _allowTensileYield(false)
{ // constructor
  useElasticBehavior(false);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::DruckerPragerPlaneStrain::~DruckerPragerPlaneStrain(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set boolean for whether tensile yield is allowed.
void
pylith::materials::DruckerPragerPlaneStrain::allowTensileYield(const bool flag)
{ // allowTensileYield
  _allowTensileYield = flag;
} // allowTensileYield

// ----------------------------------------------------------------------
// Set fit to Mohr-Coulomb surface.
void
pylith::materials::DruckerPragerPlaneStrain::fitMohrCoulomb(
						FitMohrCoulombEnum value)
{ // fitMohrCoulomb
  _fitMohrCoulomb = value;
} // fitMohrCoulomb

// ----------------------------------------------------------------------
// Set whether elastic or inelastic constitutive relations are used.
void
pylith::materials::DruckerPragerPlaneStrain::useElasticBehavior(const bool flag)
{ // useElasticBehavior
  if (flag) {
    _calcStressFn = 
      &pylith::materials::DruckerPragerPlaneStrain::_calcStressElastic;
    _calcElasticConstsFn = 
      &pylith::materials::DruckerPragerPlaneStrain::_calcElasticConstsElastic;
    _updateStateVarsFn = 
      &pylith::materials::DruckerPragerPlaneStrain::_updateStateVarsElastic;

  } else {
    _calcStressFn = 
      &pylith::materials::DruckerPragerPlaneStrain::_calcStressElastoplastic;
    _calcElasticConstsFn = 
      &pylith::materials::DruckerPragerPlaneStrain::_calcElasticConstsElastoplastic;
    _updateStateVarsFn = 
      &pylith::materials::DruckerPragerPlaneStrain::_updateStateVarsElastoplastic;
  } // if/else
} // useElasticBehavior

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::DruckerPragerPlaneStrain::_dbToProperties(
				PylithScalar* const propValues,
				const scalar_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_DruckerPragerPlaneStrain::numDBProperties == numDBValues);

  const PylithScalar density = dbValues[db_density];
  const PylithScalar vs = dbValues[db_vs];
  const PylithScalar vp = dbValues[db_vp];
  const PylithScalar frictionAngle = dbValues[db_frictionAngle];
  const PylithScalar cohesion = dbValues[db_cohesion];
  const PylithScalar dilatationAngle = dbValues[db_dilatationAngle];

  const PylithScalar pi = M_PI;
 
  if (density <= 0.0 || vs <= 0.0 || vp <= 0.0 || cohesion < 0.0 ||
      frictionAngle < 0.0 || frictionAngle > pi/2 ||
      dilatationAngle < 0.0 || dilatationAngle > pi/2 ||
      frictionAngle < dilatationAngle) {
    std::ostringstream msg;
    msg << "Spatial database returned illegal value for physical "
	<< "properties.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n"
	<< "frictionAngle: " << frictionAngle << "\n"
	<< "cohesion: " << cohesion << "\n"
	<< "dilatationAngle: " << dilatationAngle << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (fabs(frictionAngle - dilatationAngle) > 1.0e-6)
    _isJacobianSymmetric = false;

  const PylithScalar mu = density * vs*vs;
  const PylithScalar lambda = density * vp*vp - 2.0*mu;

  PylithScalar alphaYield = 0.0;
  PylithScalar beta = 0.0;
  PylithScalar alphaFlow = 0.0;
  switch (_fitMohrCoulomb) { // switch
  case MOHR_COULOMB_INSCRIBED: {
    const PylithScalar denomFriction = sqrt(3.0) * (3.0 - sin(frictionAngle));
    const PylithScalar denomDilatation = sqrt(3.0) *
      (3.0 - sin(dilatationAngle));
    alphaYield = 2.0 * sin(frictionAngle)/denomFriction;
    beta = 6.0 * cohesion * cos(frictionAngle)/denomFriction;
    alphaFlow = 2.0 * sin(dilatationAngle)/denomDilatation;
    break;
  } // MOHR_COULOMB_INSCRIBED
  case MOHR_COULOMB_MIDDLE: {
    alphaYield = sin(frictionAngle)/3.0;
    beta = cohesion * cos(frictionAngle);
    alphaFlow = sin(dilatationAngle)/3.0;
    break;
  } // MOHR_COULOMB_MIDDLE
  case MOHR_COULOMB_CIRCUMSCRIBED: {
    const PylithScalar denomFriction = sqrt(3.0) * (3.0 + sin(frictionAngle));
    const PylithScalar denomDilatation = sqrt(3.0) *
      (3.0 + sin(dilatationAngle));
    alphaYield = 2.0 * sin(frictionAngle)/denomFriction;
    beta = 6.0 * cohesion * cos(frictionAngle)/denomFriction;
    alphaFlow = 2.0 * sin(dilatationAngle)/denomDilatation;
    break;
  } // MOHR_COULOMB_CIRCUMSCRIBED
  default :
    assert(0);
    throw std::logic_error("Unknown Mohr-Coulomb fit.");
    break;
  } // switch

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
  propValues[p_alphaYield] = alphaYield;
  propValues[p_beta] = beta;
  propValues[p_alphaFlow] = alphaFlow;

  PetscLogFlops(28);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::DruckerPragerPlaneStrain::_nondimProperties(
						PylithScalar* const values,
						const int nvalues) const
{ // _nondimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _numPropsQuadPt);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_density] = 
    _normalizer->nondimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->nondimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->nondimensionalize(values[p_lambda], pressureScale);
  values[p_beta] = 
    _normalizer->nondimensionalize(values[p_beta], pressureScale);

  PetscLogFlops(4);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::DruckerPragerPlaneStrain::_dimProperties(
						PylithScalar* const values,
						const int nvalues) const
{ // _dimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _numPropsQuadPt);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_density] = 
    _normalizer->dimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->dimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->dimensionalize(values[p_lambda], pressureScale);
  values[p_beta] = 
    _normalizer->dimensionalize(values[p_beta], pressureScale);

  PetscLogFlops(4);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute initial state variables from values in spatial database.
void
pylith::materials::DruckerPragerPlaneStrain::_dbToStateVars(
				PylithScalar* const stateValues,
				const scalar_array& dbValues)
{ // _dbToStateVars
  assert(stateValues);
  const int numDBValues = dbValues.size();
  assert(_DruckerPragerPlaneStrain::numDBStateVars == numDBValues);

  const int totalSize = 5;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  memcpy(stateValues, &dbValues[0], totalSize*sizeof(PylithScalar));

} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::materials::DruckerPragerPlaneStrain::_nondimStateVars(
						PylithScalar* const values,
						const int nvalues) const
{ // _nondimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _numVarsQuadPt);

  const PylithScalar pressureScale = _normalizer->pressureScale();
  _normalizer->nondimensionalize(&values[s_stressZZInitial], 1, pressureScale);

  PetscLogFlops(1);
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::materials::DruckerPragerPlaneStrain::_dimStateVars(
						PylithScalar* const values,
						const int nvalues) const
{ // _dimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _numVarsQuadPt);

  const PylithScalar pressureScale = _normalizer->pressureScale();
  _normalizer->dimensionalize(&values[s_stressZZInitial], 1, pressureScale);

  PetscLogFlops(1);
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::DruckerPragerPlaneStrain::_calcDensity(
						PylithScalar* const density,
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
pylith::materials::DruckerPragerPlaneStrain::_stableTimeStepImplicit(
				  const PylithScalar* properties,
				  const int numProperties,
				  const PylithScalar* stateVars,
				  const int numStateVars) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  // It's unclear what to do for an elasto-plastic material, which has no
  // inherent time scale. For now, just set dtStable to a large value.
  const PylithScalar dtStable = pylith::PYLITH_MAXSCALAR;

  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Get stable time step for explicit time integration.
PylithScalar
pylith::materials::DruckerPragerPlaneStrain::_stableTimeStepExplicit(const PylithScalar* properties,
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
pylith::materials::DruckerPragerPlaneStrain::_calcStressElastic(
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
  assert(stress);
  assert(_DruckerPragerPlaneStrain::tensorSize == stressSize);
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
  assert(stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(totalStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == strainSize);
  assert(initialStress);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStressSize);
  assert(initialStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStrainSize);

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
// Compute stress tensor at location from properties as an elastoplastic
// material.
void
pylith::materials::DruckerPragerPlaneStrain::_calcStressElastoplastic(
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
{ // _calcStressElastoplastic
  assert(stress);
  assert(_DruckerPragerPlaneStrain::tensorSize == stressSize);
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
  assert(stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(totalStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == strainSize);
  assert(initialStress);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStressSize);
  assert(initialStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStrainSize);

  const int tensorSize = 3;
  const int tensorSizePS = 4;
  assert(_tensorSize == tensorSize);
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
    
  // We need to compute the plastic strain increment if state variables are
  // from previous time step.
  if (computeStateVars) {

    const PylithScalar alphaYield = properties[p_alphaYield];
    const PylithScalar beta = properties[p_beta];
    const PylithScalar alphaFlow = properties[p_alphaFlow];
    const PylithScalar mu2 = 2.0 * mu;
    const PylithScalar bulkModulus = lambda + mu2/3.0;
    const PylithScalar ae = 1.0/mu2;
    const PylithScalar am = 1.0/(3.0 * bulkModulus);

    const PylithScalar stressZZInitial = stateVars[s_stressZZInitial];

    const PylithScalar plasticStrainT[tensorSizePS] = {
      stateVars[s_plasticStrain    ],
      stateVars[s_plasticStrain + 1],
      stateVars[s_plasticStrain + 2],
      stateVars[s_plasticStrain + 3]
    };
    const PylithScalar meanPlasticStrainT = (plasticStrainT[0] +
					     plasticStrainT[1] +
					     plasticStrainT[2])/3.0;
    PylithScalar devPlasticStrainT[tensorSizePS];
    calcDeviatoric2DPS(devPlasticStrainT, plasticStrainT, meanPlasticStrainT);

    const PylithScalar diag[tensorSizePS] = { 1.0, 1.0, 1.0, 0.0 };

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

    // Initial strain values
    const PylithScalar meanStrainInitial = (initialStrain[0] +
					    initialStrain[1])/3.0;
    const PylithScalar devStrainInitial[tensorSizePS] = {
      initialStrain[0] - meanStrainInitial,
      initialStrain[1] - meanStrainInitial,
                       - meanStrainInitial,
      initialStrain[2],
    };

    // Values for current time step
    const PylithScalar meanStrainTpdt = (totalStrain[0] + totalStrain[1])/3.0;
    const PylithScalar meanStrainPPTpdt =
      meanStrainTpdt - meanPlasticStrainT - meanStrainInitial;

    // devStrainPPTpdt
    const PylithScalar strainPPTpdt[tensorSizePS] = {
      totalStrain[0] - meanStrainTpdt - devPlasticStrainT[0] -
      devStrainInitial[0],
      totalStrain[1] - meanStrainTpdt - devPlasticStrainT[1] -
      devStrainInitial[1],
      - meanStrainTpdt - devPlasticStrainT[2] - devStrainInitial[2],
      totalStrain[2] - devPlasticStrainT[3] - devStrainInitial[3]
    };

    // Compute trial elastic stresses and yield function to see if yield should
    // occur.
    const PylithScalar trialDevStress[tensorSizePS] = {
      strainPPTpdt[0]/ae + devStressInitial[0],
      strainPPTpdt[1]/ae + devStressInitial[1],
      strainPPTpdt[2]/ae + devStressInitial[2],
      strainPPTpdt[3]/ae + devStressInitial[3]
    };
    const PylithScalar trialMeanStress =
      meanStrainPPTpdt/am + meanStressInitial;
    const PylithScalar stressInvar2 =
      sqrt(0.5 * scalarProduct2DPS(trialDevStress, trialDevStress));
    const PylithScalar yieldFunction =
      3.0 * alphaYield * trialMeanStress + stressInvar2 - beta;
#if 0 // DEBUGGING
    std::cout << "Function _calcStressElastoPlastic: elastic" << std::endl;
    std::cout << "  alphaYield:       " << alphaYield << std::endl;
    std::cout << "  beta:             " << beta << std::endl;
    std::cout << "  trialMeanStress:  " << trialMeanStress << std::endl;
    std::cout << "  stressInvar2:     " << stressInvar2 << std::endl;
    std::cout << "  yieldFunction:    " << yieldFunction << std::endl;
#endif
    PetscLogFlops(62);

    // If yield function is greater than zero, compute elastoplastic stress.
    if (yieldFunction >= 0.0) {
      const PylithScalar devStressInitialProd =
	scalarProduct2DPS(devStressInitial, devStressInitial);
      const PylithScalar strainPPTpdtProd =
	scalarProduct2DPS(strainPPTpdt, strainPPTpdt);
      const PylithScalar d =
	sqrt(ae * ae * devStressInitialProd + 2.0 * ae *
	     scalarProduct2DPS(devStressInitial, strainPPTpdt) +
	     strainPPTpdtProd);

      PylithScalar plasticMult = 0.0;
      if(_allowTensileYield) {
	plasticMult = 
	  std::min(sqrt(2.0)*d,
		   2.0 * ae * am *
		   (3.0 * alphaYield * (meanStrainPPTpdt/am + meanStressInitial)
		    + d/(sqrt(2.0) * ae) 
		    - beta)/ (6.0 * alphaYield * alphaFlow * ae + am));
      } else {
	plasticMult = 
	  2.0 * ae * am *
	  (3.0 * alphaYield * (meanStrainPPTpdt/am + meanStressInitial) +
	   d/(sqrt(2.0) * ae) - beta)/
	  (6.0 * alphaYield * alphaFlow * ae + am);
      } // if/else

      const PylithScalar meanStressTpdt =
	(meanStrainPPTpdt - plasticMult * alphaFlow) / am + meanStressInitial;
      PylithScalar deltaDevPlasticStrain = 0.0;
      PylithScalar devStressTpdt = 0.0;
      PylithScalar totalStress[tensorSizePS];

#if 0 // DEBUGGING
      std::cout << "YIELDING\n";
      std::cout << "  totalStrain\n";
      for (int i=0; i < tensorSize; ++i)
	std::cout << "    " << totalStrain[i] << "\n";
      std::cout << "  plasticStrainT\n";
      for (int i=0; i < tensorSizePS; ++i)
	std::cout << "    " << plasticStrainT[i] << "\n";
      std::cout << "  devPlasticStrainT\n";
      for (int i=0; i < tensorSizePS; ++i)
	std::cout << "    " << devPlasticStrainT[i] << "\n";
      std::cout << "  strainPPTpdt\n";
      for (int i=0; i < tensorSizePS; ++i)
	std::cout << "    " << strainPPTpdt[i] << "\n";
      std::cout << "  meanPlasticStrainT: " << meanPlasticStrainT << "\n"
		<< "  meanStressInitial: " << meanStressInitial << "\n"
		<< "  meanStrainInitial: " << meanStrainInitial << "\n"
		<< "  meanStrainTpdt: " << meanStrainTpdt << "\n"
		<< "  meanStrainPPTpdt: " << meanStrainPPTpdt << "\n"
		<< "  plasticMult: " << plasticMult
		<< std::endl;
#endif

      if (d > 0.0 || !_allowTensileYield) {
	for (int iComp=0; iComp < tensorSizePS; ++iComp) {
	  deltaDevPlasticStrain =
	    plasticMult * (strainPPTpdt[iComp] + ae * devStressInitial[iComp])/
	    (sqrt(2.0) * d);
	  devStressTpdt =
	    (strainPPTpdt[iComp] - deltaDevPlasticStrain)/ae +
	    devStressInitial[iComp];
	  totalStress[iComp] = devStressTpdt + diag[iComp] * meanStressTpdt;
	} // for
      } else {
	for (int iComp=0; iComp < tensorSizePS; ++iComp) {
	  devStressTpdt = (strainPPTpdt[iComp])/ae + devStressInitial[iComp];
	  totalStress[iComp] = devStressTpdt + diag[iComp] * meanStressTpdt;
	} // for
      } // if/else
      stress[0] = totalStress[0];
      stress[1] = totalStress[1];
      stress[2] = totalStress[3];

#if 0 // DEBUGGING
      std::cout << "  totalStress\n";
      for (int i=0; i < tensorSizePS; ++i)
	std::cout << "    " << totalStress[i] << "\n";
#endif

    PetscLogFlops(51 + 11 * tensorSizePS);

    } else {
      // No plastic strain.
      const PylithScalar meanStressTpdt = meanStrainPPTpdt/am +
	meanStressInitial;
      stress[0] = strainPPTpdt[0]/ae + devStressInitial[0] + meanStressTpdt; 
      stress[1] = strainPPTpdt[1]/ae + devStressInitial[1] + meanStressTpdt; 
      stress[2] = strainPPTpdt[3]/ae + devStressInitial[3]; 

      PetscLogFlops(10);
    } // if

    // If state variables have already been updated, the plastic strain for the
    // time step has already been computed.
  } else {
    const PylithScalar mu2 = 2.0 * mu;
    const PylithScalar plasticStrainTpdt[tensorSizePS] = {
      stateVars[s_plasticStrain],
      stateVars[s_plasticStrain + 1],
      stateVars[s_plasticStrain + 2],
      stateVars[s_plasticStrain + 3]
    };

    const PylithScalar e11 = totalStrain[0] - plasticStrainTpdt[0] -
      initialStrain[0];
    const PylithScalar e22 = totalStrain[1] - plasticStrainTpdt[1] -
      initialStrain[1];
    const PylithScalar e33 = -plasticStrainTpdt[2];
    const PylithScalar e12 = totalStrain[2] - plasticStrainTpdt[3] -
      initialStrain[2];

    const PylithScalar traceStrainTpdt = e11 + e22 + e33;
    const PylithScalar s12 = lambda * traceStrainTpdt;

    stress[0] = s12 + mu2 * e11 + initialStress[0];
    stress[1] = s12 + mu2 * e22 + initialStress[1];
    stress[2] = mu2 * e12 + initialStress[2];

    PetscLogFlops(17);

  } // else
#if 0 // DEBUGGING
  //**  This debugging code won't work for plane strain since we're missing
  // stressZZ. Not sure if this can be fixed.
  const PylithScalar alphaYield = properties[p_alphaYield];
  const PylithScalar beta = properties[p_beta];
  const PylithScalar alphaFlow = properties[p_alphaFlow];
  const PylithScalar meanStressTest = (stress[0] + stress[1])/3.0;
  const PylithScalar devStressTest[] = { stress[0] - meanStressTest,
				   stress[1] - meanStressTest,
				   stress[2] - meanStressTest,
				   stress[3],
				   stress[4],
				   stress[5]};
  const PylithScalar stressInvar2Test =
    sqrt(0.5 *
	 pylith::materials::ElasticMaterial::scalarProduct3D(devStressTest,
							     devStressTest));
  
  const PylithScalar yieldFunctionTest = 3.0 * alphaYield * meanStressTest +
      stressInvar2Test - beta;
  std::cout << "Function _calcStressElastoPlastic: end" << std::endl;
  std::cout << "  alphaYield:        " << alphaYield << std::endl;
  std::cout << "  beta:              " << beta << std::endl;
  std::cout << "  meanStressTest:    " << meanStressTest << std::endl;
  std::cout << "  stressInvar2Test:  " << stressInvar2Test << std::endl;
  std::cout << "  yieldFunctionTest: " << yieldFunctionTest << std::endl;
#endif

} // _calcStressElastoplastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::DruckerPragerPlaneStrain::_calcElasticConstsElastic(
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
  assert(elasticConsts);
  assert(_DruckerPragerPlaneStrain::numElasticConsts == numElasticConsts);
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
  assert(stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(totalStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == strainSize);
  assert(initialStress);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStressSize);
  assert(initialStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStrainSize);
 
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
// as an elastoplastic material.
void
pylith::materials::DruckerPragerPlaneStrain::_calcElasticConstsElastoplastic(
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
{ // _calcElasticConstsElastoplastic
  assert(elasticConsts);
  assert(_DruckerPragerPlaneStrain::numElasticConsts == numElasticConsts);
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
  assert(stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(totalStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == strainSize);
  assert(initialStress);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStressSize);
  assert(initialStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStrainSize);

  // Duplicate functionality of _calcStressElastoplastic
  // Get properties
  const int tensorSize = 3;
  const int tensorSizePS = 4;
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar alphaYield = properties[p_alphaYield];
  const PylithScalar beta = properties[p_beta];
  const PylithScalar alphaFlow = properties[p_alphaFlow];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2/3.0;
  const PylithScalar ae = 1.0/mu2;
  const PylithScalar am = 1.0/(3.0 * bulkModulus);

  const PylithScalar stressZZInitial = stateVars[s_stressZZInitial];
  
  // Get state variables from previous time step
  const PylithScalar plasticStrainT[tensorSizePS] = {
    stateVars[s_plasticStrain    ],
    stateVars[s_plasticStrain + 1],
    stateVars[s_plasticStrain + 2],
    stateVars[s_plasticStrain + 3]};
  const PylithScalar meanPlasticStrainT = (plasticStrainT[0] +
					   plasticStrainT[1] +
					   plasticStrainT[2])/3.0;
  PylithScalar devPlasticStrainT[tensorSizePS];
  calcDeviatoric2DPS(devPlasticStrainT, plasticStrainT, meanPlasticStrainT);

  const PylithScalar diag[tensorSize] = { 1.0, 1.0, 0.0 };

  // Initial stress values
  const PylithScalar meanStressInitial = (initialStress[0] +
					  initialStress[1] +
					  stressZZInitial)/3.0;
  const PylithScalar devStressInitial[tensorSizePS] = {
    initialStress[0] - meanStressInitial,
    initialStress[1] - meanStressInitial,
    stressZZInitial - meanStressInitial,
    initialStress[2]};

  // Initial strain values
  const PylithScalar meanStrainInitial = (initialStrain[0] +
					  initialStrain[1])/3.0;
  const PylithScalar devStrainInitial[tensorSizePS] = {
    initialStrain[0] - meanStrainInitial,
    initialStrain[1] - meanStrainInitial,
    -meanStrainInitial,
    initialStrain[2]};

  // Values for current time step
  const PylithScalar meanStrainTpdt = (totalStrain[0] + totalStrain[1])/3.0;
  const PylithScalar meanStrainPPTpdt = meanStrainTpdt - meanPlasticStrainT -
    meanStrainInitial;
  
  const PylithScalar strainPPTpdt[tensorSizePS] = {
    totalStrain[0] - meanStrainTpdt - devPlasticStrainT[0] -
    devStrainInitial[0],
    totalStrain[1] - meanStrainTpdt - devPlasticStrainT[1] -
    devStrainInitial[1],
                   - meanStrainTpdt - devPlasticStrainT[2] -
    devStrainInitial[2],
    totalStrain[2] - devPlasticStrainT[3] - devStrainInitial[3],
  };
  
  // Compute trial elastic stresses and yield function to see if yield should
  // occur.
  const PylithScalar trialDevStress[tensorSizePS] = {
    strainPPTpdt[0]/ae + devStressInitial[0],
    strainPPTpdt[1]/ae + devStressInitial[1],
    strainPPTpdt[2]/ae + devStressInitial[2],
    strainPPTpdt[3]/ae + devStressInitial[3],
  };
  const PylithScalar trialMeanStress = meanStrainPPTpdt/am + meanStressInitial;
  const PylithScalar stressInvar2 =
    sqrt(0.5 * scalarProduct2DPS(trialDevStress, trialDevStress));
  const PylithScalar yieldFunction = 3.0 * alphaYield * trialMeanStress + stressInvar2 - beta;
#if 0 // DEBUGGING
  std::cout << "Function _calcElasticConstsElastoPlastic:" << std::endl;
  std::cout << "  alphaYield:       " << alphaYield << std::endl;
  std::cout << "  beta:             " << beta << std::endl;
  std::cout << "  trialMeanStress:  " << trialMeanStress << std::endl;
  std::cout << "  stressInvar2:     " << stressInvar2 << std::endl;
  std::cout << "  yieldFunction:    " << yieldFunction << std::endl;
#endif
  PetscLogFlops(62);
  
  // If yield function is greater than zero, compute elastoplastic stress and
  // corresponding tangent matrix.
  if (yieldFunction >= 0.0) {
    const PylithScalar devStressInitialProd = 
      scalarProduct2DPS(devStressInitial, devStressInitial);
    const PylithScalar strainPPTpdtProd =
      scalarProduct2DPS(strainPPTpdt, strainPPTpdt);
    const PylithScalar d = 
      sqrt(ae * ae * devStressInitialProd + 2.0 * ae *
	   scalarProduct2DPS(devStressInitial, strainPPTpdt) +
	   strainPPTpdtProd);
    const PylithScalar plasticFac = 2.0 * ae * am/
      (6.0 * alphaYield * alphaFlow * ae + am);
    const PylithScalar meanStrainFac = 3.0 * alphaYield;
    const PylithScalar dFac = 1.0/(sqrt(2.0) * ae);

    PylithScalar plasticMult = 0.0;
    bool tensileYield = false;
    PylithScalar dFac2 = 0.0;
    if (_allowTensileYield) {
      const PylithScalar testMult = plasticFac *
	(meanStrainFac * (meanStrainPPTpdt/am + meanStressInitial) +
	 dFac * d - beta);
      tensileYield = (sqrt(2.0) * d < testMult) ? true: false;
      plasticMult = tensileYield ? sqrt(2.0) * d : testMult;
      dFac2 = (d > 0.0) ? 1.0/(sqrt(2.0) * d) : 0.0;
    } else {
      plasticMult = plasticFac *
	(meanStrainFac * (meanStrainPPTpdt/am + meanStressInitial) +
	 dFac * d - beta);
      dFac2 = 1.0/(sqrt(2.0) * d);
    } // if/else

    // Define some constants, vectors, and matrices.
    const PylithScalar third = 1.0/3.0;
    const PylithScalar dEdEpsilon[3][3] = {
      { 2.0 * third,      -third,      0.0},
      {      -third, 2.0 * third,      0.0},
      {         0.0,         0.0,      1.0}};
    const PylithScalar vec1[3] = {
      strainPPTpdt[0] + ae * devStressInitial[0],
      strainPPTpdt[1] + ae * devStressInitial[1],
      strainPPTpdt[3] + ae * devStressInitial[3],
    };
    
    PylithScalar dDeltaEdEpsilon = 0.0;

    // Compute elasticity matrix.
    if (d > 0.0) {
      const PylithScalar dDdEpsilon[3] = {vec1[0]/d,
					  vec1[1]/d,
					  2.0 * vec1[2]/d};

      PylithScalar dLambdadEpsilon[tensorSize];
      if (tensileYield) {
	dLambdadEpsilon[0] = sqrt(2.0) * dDdEpsilon[0];
	dLambdadEpsilon[1] = sqrt(2.0) * dDdEpsilon[1];
	dLambdadEpsilon[2] = sqrt(2.0) * dDdEpsilon[2];
      } else {
	dLambdadEpsilon[0] = plasticFac *
	  (alphaYield/am + dFac * dDdEpsilon[0]);
	dLambdadEpsilon[1] = plasticFac *
	  (alphaYield/am + dFac * dDdEpsilon[1]);
	dLambdadEpsilon[2] = plasticFac * dFac * dDdEpsilon[2];
      } // else
      for (int iComp=0; iComp < tensorSize; ++iComp) {
	for (int jComp=0; jComp < tensorSize; ++jComp) {
	  int iCount = jComp + tensorSize * iComp;
	  dDeltaEdEpsilon = dFac2 * (vec1[iComp] *
				     (dLambdadEpsilon[jComp] -
				      plasticMult * dDdEpsilon[jComp]/d) +
				     plasticMult * dEdEpsilon[iComp][jComp]);
	  elasticConsts[iCount] = (dEdEpsilon[iComp][jComp] -
				   dDeltaEdEpsilon)/ae +
	    diag[iComp] * (third * diag[jComp] -
			   alphaFlow * dLambdadEpsilon[jComp])/am;
	} // for
      } // for
    } else {
      const PylithScalar dLambdadEpsilon[3] = {0.0,	
					       0.0,	
					       0.0};
      for (int iComp=0; iComp < tensorSize; ++iComp) {
	for (int jComp=0; jComp < tensorSize; ++jComp) {
	  int iCount = jComp + tensorSize * iComp;
	  elasticConsts[iCount] = (dEdEpsilon[iComp][jComp])/ae +
	    diag[iComp] * (third * diag[jComp] -
			   alphaFlow * dLambdadEpsilon[jComp])/am;
	} // for
      } // for
    } // if/else

    PetscLogFlops(76 + tensorSize * tensorSize * 15);

  } else {
    // No plastic strain.
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

    PetscLogFlops(1);
  } // else

} // _calcElasticConstsElastoplastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::DruckerPragerPlaneStrain::_updateStateVarsElastic(
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
  assert(stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
  assert(totalStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == strainSize);
  assert(initialStress);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStressSize);
  assert(initialStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStrainSize);

  for (int iComp=0; iComp < 4; ++iComp) {
    stateVars[s_plasticStrain+iComp] = 0.0;
  } // for

  _needNewJacobian = true;
} // _updateStateVarsElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::DruckerPragerPlaneStrain::_updateStateVarsElastoplastic(
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
{ // _updateStateVarsElastoplastic
  assert(stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
  assert(totalStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == strainSize);
  assert(initialStress);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStressSize);
  assert(initialStrain);
  assert(_DruckerPragerPlaneStrain::tensorSize == initialStrainSize);

  // For now, we are duplicating the functionality of _calcStressElastoplastic,
  // since otherwise we would have to redo a lot of calculations.

  const int tensorSizePS = 4;
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar alphaYield = properties[p_alphaYield];
  const PylithScalar beta = properties[p_beta];
  const PylithScalar alphaFlow = properties[p_alphaFlow];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2/3.0;
  const PylithScalar ae = 1.0/mu2;
  const PylithScalar am = 1.0/(3.0 * bulkModulus);

  const PylithScalar stressZZInitial = stateVars[s_stressZZInitial];

  const PylithScalar plasticStrainT[tensorSizePS] = {
    stateVars[s_plasticStrain    ],
    stateVars[s_plasticStrain + 1],
    stateVars[s_plasticStrain + 2],
    stateVars[s_plasticStrain + 3]};
  const PylithScalar meanPlasticStrainT = (plasticStrainT[0] +
					   plasticStrainT[1] +
					   plasticStrainT[2])/3.0;
  PylithScalar devPlasticStrainT[tensorSizePS];
  calcDeviatoric2DPS(devPlasticStrainT, plasticStrainT, meanPlasticStrainT);

  const PylithScalar diag[tensorSizePS] = { 1.0, 1.0, 1.0, 0.0 };

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

  // Initial strain values
  const PylithScalar meanStrainInitial = (initialStrain[0] +
					  initialStrain[1])/3.0;
  const PylithScalar devStrainInitial[tensorSizePS] = {
    initialStrain[0] - meanStrainInitial,
    initialStrain[1] - meanStrainInitial,
    0.0 - meanStrainInitial,
    initialStrain[2]};

  // Values for current time step
  const PylithScalar meanStrainTpdt = (totalStrain[0] + totalStrain[1])/3.0;
  const PylithScalar meanStrainPPTpdt = meanStrainTpdt - meanPlasticStrainT -
    meanStrainInitial;

  const PylithScalar strainPPTpdt[tensorSizePS] = {
    totalStrain[0] - meanStrainTpdt - devPlasticStrainT[0] -
    devStrainInitial[0],
    totalStrain[1] - meanStrainTpdt - devPlasticStrainT[1] -
    devStrainInitial[1],
    - meanStrainTpdt - devPlasticStrainT[2] - devStrainInitial[2],
    totalStrain[2] - devPlasticStrainT[3] - devStrainInitial[3]
  };

  // Compute trial elastic stresses and yield function to see if yield should
  // occur.
  const PylithScalar trialDevStress[tensorSizePS] = {
    strainPPTpdt[0]/ae + devStressInitial[0],
    strainPPTpdt[1]/ae + devStressInitial[1],
    strainPPTpdt[2]/ae + devStressInitial[2],
    strainPPTpdt[3]/ae + devStressInitial[3]};
  const PylithScalar trialMeanStress = meanStrainPPTpdt/am + meanStressInitial;
  const PylithScalar stressInvar2 =
    sqrt(0.5 * scalarProduct2DPS(trialDevStress, trialDevStress));
  const PylithScalar yieldFunction =
    3.0 * alphaYield * trialMeanStress + stressInvar2 - beta;

#if 0 // DEBUGGING
  std::cout << "Function _updateStateVarsElastoPlastic:" << std::endl;
  std::cout << "  alphaYield:       " << alphaYield << std::endl;
  std::cout << "  beta:             " << beta << std::endl;
  std::cout << "  trialMeanStress:  " << trialMeanStress << std::endl;
  std::cout << "  stressInvar2:     " << stressInvar2 << std::endl;
  std::cout << "  yieldFunction:    " << yieldFunction << std::endl;
#endif
  PetscLogFlops(62);

  // If yield function is greater than zero, compute plastic strains.
  // Otherwise, plastic strains remain the same.
  if (yieldFunction >= 0.0) {
    const PylithScalar devStressInitialProd = 
      scalarProduct2DPS(devStressInitial, devStressInitial);
    const PylithScalar strainPPTpdtProd =
      scalarProduct2DPS(strainPPTpdt, strainPPTpdt);
    const PylithScalar d =
      sqrt(ae * ae * devStressInitialProd + 2.0 * ae *
	   scalarProduct2DPS(devStressInitial, strainPPTpdt) +
	   strainPPTpdtProd);
    const PylithScalar plasticFac = 2.0 * ae * am/
      (6.0 * alphaYield * alphaFlow * ae + am);
    const PylithScalar meanStrainFac = 3.0 * alphaYield;
    const PylithScalar dFac = 1.0/(sqrt(2.0) * ae);

    PylithScalar plasticMult =  0.0;
    const PylithScalar plasticMultNormal = plasticFac *
      (meanStrainFac * (meanStrainPPTpdt/am + meanStressInitial) +
       dFac * d - beta);
  const PylithScalar plasticMultTensile = sqrt(2.0) * d;
    if (_allowTensileYield) {
      plasticMult = std::min(plasticMultNormal, plasticMultTensile);
    } else if (plasticMultNormal <= plasticMultTensile) {
      plasticMult = plasticMultNormal;
    } else {
      std::ostringstream msg;
      const PylithScalar stressInvar2Comp = 0.5 *
	(sqrt(2.0) * d - plasticMultNormal)/ae;
      msg << "Infeasible stress state. Cannot project back to yield surface.\n"
	  << "  alphaYield:         " << alphaYield << "\n"
	  << "  alphaFlow:          " << alphaFlow << "\n"
	  << "  beta:               " << beta << "\n"
	  << "  d:                  " << d << "\n"
	  << "  plasticMultNormal:  " << plasticMultNormal << "\n"
	  << "  plasticMultTensile: " << plasticMultTensile << "\n"
	  << "  yieldFunction:      " << yieldFunction << "\n"
	  << "  stressInvar2Comp:   " << stressInvar2Comp << "\n";
      throw std::runtime_error(msg.str());
    } // if/else

    const PylithScalar deltaMeanPlasticStrain = plasticMult * alphaFlow;
    PylithScalar deltaDevPlasticStrain = 0.0;
    if (d > 0.0 || !_allowTensileYield) {
      for (int iComp=0; iComp < tensorSizePS; ++iComp) {
	deltaDevPlasticStrain = plasticMult *(strainPPTpdt[iComp] +
					      ae * devStressInitial[iComp])/
	  (sqrt(2.0) * d);
	stateVars[s_plasticStrain+iComp] += deltaDevPlasticStrain +
	  diag[iComp] * deltaMeanPlasticStrain;
      } // for
    } else {
      for (int iComp=0; iComp < tensorSizePS; ++iComp) {
	stateVars[s_plasticStrain+iComp] +=
	  diag[iComp] * deltaMeanPlasticStrain;
      } // for
    } // if/else
    
    PetscLogFlops(48 + 9 * tensorSizePS);

  } // if

  _needNewJacobian = true;

} // _updateStateVarsElastoplastic

// End of file 
