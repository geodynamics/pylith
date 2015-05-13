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

#include "DruckerPrager3D.hh" // implementation of object methods

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
    namespace _DruckerPrager3D{

      /// Dimension of material.
      const int dimension = 3;

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 6;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 36;

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
      const int numStateVars = 1;

      /// State variables.
      const Metadata::ParamDescription stateVars[numStateVars] = {
	{ "plastic_strain", tensorSize, pylith::topology::FieldBase::TENSOR }
      };

      // Values expected in state variables spatial database.
      const int numDBStateVars = 6;
      const char* dbStateVars[6] = {
	"plastic-strain-xx",
	"plastic-strain-yy",
	"plastic-strain-zz",
	"plastic-strain-xy",
	"plastic-strain-yz",
	"plastic-strain-xz",
      };

    } // _DruckerPrager3D
  } // materials
} // pylith

// Indices of physical properties.
const int pylith::materials::DruckerPrager3D::p_density = 0;

const int pylith::materials::DruckerPrager3D::p_mu = 
  pylith::materials::DruckerPrager3D::p_density + 1;

const int pylith::materials::DruckerPrager3D::p_lambda = 
  pylith::materials::DruckerPrager3D::p_mu + 1;

const int pylith::materials::DruckerPrager3D::p_alphaYield = 
  pylith::materials::DruckerPrager3D::p_lambda + 1;

const int pylith::materials::DruckerPrager3D::p_beta = 
  pylith::materials::DruckerPrager3D::p_alphaYield + 1;

const int pylith::materials::DruckerPrager3D::p_alphaFlow = 
  pylith::materials::DruckerPrager3D::p_beta + 1;

// Indices of property database values (order must match dbProperties).
const int pylith::materials::DruckerPrager3D::db_density = 0;

const int pylith::materials::DruckerPrager3D::db_vs = 
  pylith::materials::DruckerPrager3D::db_density + 1;

const int pylith::materials::DruckerPrager3D::db_vp = 
  pylith::materials::DruckerPrager3D::db_vs + 1;

const int pylith::materials::DruckerPrager3D::db_frictionAngle = 
  pylith::materials::DruckerPrager3D::db_vp + 1;

const int pylith::materials::DruckerPrager3D::db_cohesion = 
  pylith::materials::DruckerPrager3D::db_frictionAngle + 1;

const int pylith::materials::DruckerPrager3D::db_dilatationAngle = 
  pylith::materials::DruckerPrager3D::db_cohesion + 1;

// Indices of state variables.
const int pylith::materials::DruckerPrager3D::s_plasticStrain = 0;

// Indices of state variable database values (order must match dbStateVars).
const int pylith::materials::DruckerPrager3D::db_plasticStrain = 0;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::DruckerPrager3D::DruckerPrager3D(void) :
  ElasticMaterial(_DruckerPrager3D::dimension,
		  _DruckerPrager3D::tensorSize,
		  _DruckerPrager3D::numElasticConsts,
		  Metadata(_DruckerPrager3D::properties,
			   _DruckerPrager3D::numProperties,
			   _DruckerPrager3D::dbProperties,
			   _DruckerPrager3D::numDBProperties,
			   _DruckerPrager3D::stateVars,
			   _DruckerPrager3D::numStateVars,
			   _DruckerPrager3D::dbStateVars,
			   _DruckerPrager3D::numDBStateVars)),
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
pylith::materials::DruckerPrager3D::~DruckerPrager3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set boolean for whether tensile yield is allowed.
void
pylith::materials::DruckerPrager3D::allowTensileYield(const bool flag)
{ // allowTensileYield
  _allowTensileYield = flag;
} // allowTensileYield

// ----------------------------------------------------------------------
// Set fit to Mohr-Coulomb surface.
void
pylith::materials::DruckerPrager3D::fitMohrCoulomb(FitMohrCoulombEnum value)
{ // fitMohrCoulomb
  _fitMohrCoulomb = value;
} // fitMohrCoulomb

// ----------------------------------------------------------------------
// Set whether elastic or inelastic constitutive relations are used.
void
pylith::materials::DruckerPrager3D::useElasticBehavior(const bool flag)
{ // useElasticBehavior
  if (flag) {
    _calcStressFn = 
      &pylith::materials::DruckerPrager3D::_calcStressElastic;
    _calcElasticConstsFn = 
      &pylith::materials::DruckerPrager3D::_calcElasticConstsElastic;
    _updateStateVarsFn = 
      &pylith::materials::DruckerPrager3D::_updateStateVarsElastic;

  } else {
    _calcStressFn = 
      &pylith::materials::DruckerPrager3D::_calcStressElastoplastic;
    _calcElasticConstsFn = 
      &pylith::materials::DruckerPrager3D::_calcElasticConstsElastoplastic;
    _updateStateVarsFn = 
      &pylith::materials::DruckerPrager3D::_updateStateVarsElastoplastic;
  } // if/else
} // useElasticBehavior

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::DruckerPrager3D::_dbToProperties(
				PylithScalar* const propValues,
				const scalar_array& dbValues)
{ // _dbToProperties
  assert(propValues);
  const int numDBValues = dbValues.size();
  assert(_DruckerPrager3D::numDBProperties == numDBValues);

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
pylith::materials::DruckerPrager3D::_nondimProperties(PylithScalar* const values,
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
    _normalizer->nondimensionalize(values[p_beta],
				   pressureScale);

  PetscLogFlops(4);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::DruckerPrager3D::_dimProperties(PylithScalar* const values,
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
pylith::materials::DruckerPrager3D::_dbToStateVars(
				PylithScalar* const stateValues,
				const scalar_array& dbValues)
{ // _dbToStateVars
  assert(stateValues);
  const int numDBValues = dbValues.size();
  assert(_DruckerPrager3D::numDBStateVars == numDBValues);

  const int totalSize = _tensorSize;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  memcpy(&stateValues[s_plasticStrain], &dbValues[db_plasticStrain],
	 _tensorSize*sizeof(PylithScalar));
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::materials::DruckerPrager3D::_nondimStateVars(PylithScalar* const values,
						const int nvalues) const
{ // _nondimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _numVarsQuadPt);
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::materials::DruckerPrager3D::_dimStateVars(PylithScalar* const values,
					     const int nvalues) const
{ // _dimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _numVarsQuadPt);
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::DruckerPrager3D::_calcDensity(PylithScalar* const density,
					    const PylithScalar* properties,
					    const int numProperties,
					    const PylithScalar* stateVars,
					    const int numStateVars)
{ // _calcDensity
  assert(density);
  assert(properties);
  assert(_numPropsQuadPt == numProperties);

  density[0] = properties[p_density];
} // _calcDensity

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::DruckerPrager3D::_stableTimeStepImplicit(
				  const PylithScalar* properties,
				  const int numProperties,
				  const PylithScalar* stateVars,
				  const int numStateVars) const
{ // _stableTimeStepImplicit
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
  assert(stateVars);
  assert(_numVarsQuadPt == numStateVars);

  // It's unclear what to do for an elasto-plastic material, which has no
  // inherent time scale. For now, just set dtStable to a large value.
  const PylithScalar dtStable = pylith::PYLITH_MAXSCALAR;

  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Get stable time step for explicit time integration.
PylithScalar
pylith::materials::DruckerPrager3D::_stableTimeStepExplicit(const PylithScalar* properties,
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
pylith::materials::DruckerPrager3D::_calcStressElastic(
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
  assert(_DruckerPrager3D::tensorSize == stressSize);
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
  assert(stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(totalStrain);
  assert(_DruckerPrager3D::tensorSize == strainSize);
  assert(initialStress);
  assert(_DruckerPrager3D::tensorSize == initialStressSize);
  assert(initialStrain);
  assert(_DruckerPrager3D::tensorSize == initialStrainSize);

  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar mu2 = 2.0 * mu;

  const PylithScalar e11 = totalStrain[0] - initialStrain[0];
  const PylithScalar e22 = totalStrain[1] - initialStrain[1];
  const PylithScalar e33 = totalStrain[2] - initialStrain[2];
  const PylithScalar e12 = totalStrain[3] - initialStrain[3];
  const PylithScalar e23 = totalStrain[4] - initialStrain[4];
  const PylithScalar e13 = totalStrain[5] - initialStrain[5];
  
  const PylithScalar traceStrainTpdt = e11 + e22 + e33;
  const PylithScalar s123 = lambda * traceStrainTpdt;

  stress[0] = s123 + mu2*e11 + initialStress[0];
  stress[1] = s123 + mu2*e22 + initialStress[1];
  stress[2] = s123 + mu2*e33 + initialStress[2];
  stress[3] = mu2 * e12 + initialStress[3];
  stress[4] = mu2 * e23 + initialStress[4];
  stress[5] = mu2 * e13 + initialStress[5];

  PetscLogFlops(25);
} // _calcStressElastic

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as an elastoplastic
// material.
void
pylith::materials::DruckerPrager3D::_calcStressElastoplastic(
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
  assert(_DruckerPrager3D::tensorSize == stressSize);
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
  assert(stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(totalStrain);
  assert(_DruckerPrager3D::tensorSize == strainSize);
  assert(initialStress);
  assert(_DruckerPrager3D::tensorSize == initialStressSize);
  assert(initialStrain);
  assert(_DruckerPrager3D::tensorSize == initialStrainSize);

  const int tensorSize = 6;
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

    const PylithScalar plasticStrainT[tensorSize] = {
      stateVars[s_plasticStrain  ],
      stateVars[s_plasticStrain+1],
      stateVars[s_plasticStrain+2],
      stateVars[s_plasticStrain+3],
      stateVars[s_plasticStrain+4],
      stateVars[s_plasticStrain+5],
    };
    const PylithScalar meanPlasticStrainT = (plasticStrainT[0] +
					     plasticStrainT[1] +
					     plasticStrainT[2])/3.0;
    PylithScalar devPlasticStrainT[tensorSize];
    calcDeviatoric3D(devPlasticStrainT, plasticStrainT, meanPlasticStrainT);
    
    const PylithScalar diag[tensorSize] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
    
    // Initial stress values
    const PylithScalar meanStressInitial = (initialStress[0] +
					    initialStress[1] +
					    initialStress[2])/3.0;
    PylithScalar devStressInitial[tensorSize];
    calcDeviatoric3D(devStressInitial, initialStress, meanStressInitial);
    
    // Initial strain values
    const PylithScalar meanStrainInitial = (initialStrain[0] +
					    initialStrain[1] +
					    initialStrain[2])/3.0;
    PylithScalar devStrainInitial[tensorSize];
    calcDeviatoric3D(devStrainInitial, initialStrain, meanStrainInitial);
    
    // Values for current time step
    const PylithScalar e11 = totalStrain[0];
    const PylithScalar e22 = totalStrain[1];
    const PylithScalar e33 = totalStrain[2];
    const PylithScalar meanStrainTpdt = (e11 + e22 + e33)/3.0;
    const PylithScalar meanStrainPPTpdt =
      meanStrainTpdt - meanPlasticStrainT - meanStrainInitial;

    const PylithScalar strainPPTpdt[tensorSize] = {
      totalStrain[0] - meanStrainTpdt - devPlasticStrainT[0] -
      devStrainInitial[0],
      totalStrain[1] - meanStrainTpdt - devPlasticStrainT[1] -
      devStrainInitial[1],
      totalStrain[2] - meanStrainTpdt - devPlasticStrainT[2] -
      devStrainInitial[2],
      totalStrain[3] - devPlasticStrainT[3] - devStrainInitial[3],
      totalStrain[4] - devPlasticStrainT[4] - devStrainInitial[4],
      totalStrain[5] - devPlasticStrainT[5] - devStrainInitial[5],
    };

    // Compute trial elastic stresses and yield function to see if yield should
    // occur.
    const PylithScalar trialDevStress[tensorSize] = { 
      strainPPTpdt[0]/ae + devStressInitial[0],
      strainPPTpdt[1]/ae + devStressInitial[1],
      strainPPTpdt[2]/ae + devStressInitial[2],
      strainPPTpdt[3]/ae + devStressInitial[3],
      strainPPTpdt[4]/ae + devStressInitial[4],
      strainPPTpdt[5]/ae + devStressInitial[5],
    };
    const PylithScalar trialMeanStress =
      meanStrainPPTpdt/am + meanStressInitial;
    const PylithScalar stressInvar2 =
      sqrt(0.5 * scalarProduct3D(trialDevStress, trialDevStress));
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
    PetscLogFlops(76);

    // If yield function is greater than zero, compute elastoplastic stress.
    if (yieldFunction >= 0.0) {
      const PylithScalar devStressInitialProd = 
	scalarProduct3D(devStressInitial, devStressInitial);
      const PylithScalar strainPPTpdtProd = 
	scalarProduct3D(strainPPTpdt, strainPPTpdt);
      const PylithScalar d =
	sqrt(ae * ae * devStressInitialProd + 2.0 * ae *
	   scalarProduct3D(devStressInitial, strainPPTpdt) + strainPPTpdtProd);
      
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
	  (3.0 * alphaYield * (meanStrainPPTpdt/am + meanStressInitial)
	   + d/(sqrt(2.0) * ae) - beta)/
	  (6.0 * alphaYield * alphaFlow * ae + am);
      } // if/else

      const PylithScalar meanStressTpdt =
	(meanStrainPPTpdt - plasticMult * alphaFlow)/am + meanStressInitial;
      PylithScalar deltaDevPlasticStrain = 0.0;
      PylithScalar devStressTpdt = 0.0;
      if (d > 0.0 || !_allowTensileYield) {
	for (int iComp=0; iComp < tensorSize; ++iComp) {
	  deltaDevPlasticStrain =
	    plasticMult * (strainPPTpdt[iComp] + ae * devStressInitial[iComp])/
	    (sqrt(2.0) * d);
	  devStressTpdt =
	    (strainPPTpdt[iComp] - deltaDevPlasticStrain)/ae +
	    devStressInitial[iComp];
	  stress[iComp] = devStressTpdt + diag[iComp] * meanStressTpdt;
	} // for
      } else {
	for (int iComp=0; iComp < tensorSize; ++iComp) {
	  devStressTpdt = (strainPPTpdt[iComp])/ae + devStressInitial[iComp];
	  stress[iComp] = devStressTpdt + diag[iComp] * meanStressTpdt;
	} // for
      } // if/else

    PetscLogFlops(62 + 11 * tensorSize);

    } else {
      // No plastic strain.
      const PylithScalar meanStressTpdt = meanStrainPPTpdt/am + meanStressInitial;
      stress[0] = strainPPTpdt[0]/ae + devStressInitial[0] + meanStressTpdt; 
      stress[1] = strainPPTpdt[1]/ae + devStressInitial[1] + meanStressTpdt; 
      stress[2] = strainPPTpdt[2]/ae + devStressInitial[2] + meanStressTpdt; 
      stress[3] = strainPPTpdt[3]/ae + devStressInitial[3]; 
      stress[4] = strainPPTpdt[4]/ae + devStressInitial[4]; 
      stress[5] = strainPPTpdt[5]/ae + devStressInitial[5]; 
    } // if

  } else {
    // If state variables have already been updated, the plastic strain for the
    // time step has already been computed.
    const PylithScalar mu2 = 2.0 * mu;
    const PylithScalar plasticStrainTpdt[tensorSize] = {
      stateVars[s_plasticStrain  ],
      stateVars[s_plasticStrain+1],
      stateVars[s_plasticStrain+2],
      stateVars[s_plasticStrain+3],
      stateVars[s_plasticStrain+4],
      stateVars[s_plasticStrain+5],
    };

    const PylithScalar e11 =
      totalStrain[0] - plasticStrainTpdt[0] - initialStrain[0];
    const PylithScalar e22 =
      totalStrain[1] - plasticStrainTpdt[1] - initialStrain[1];
    const PylithScalar e33 =
      totalStrain[2] - plasticStrainTpdt[2] - initialStrain[2];
    const PylithScalar e12 =
      totalStrain[3] - plasticStrainTpdt[3] - initialStrain[3];
    const PylithScalar e23 =
      totalStrain[4] - plasticStrainTpdt[4] - initialStrain[4];
    const PylithScalar e13 =
      totalStrain[5] - plasticStrainTpdt[5] - initialStrain[5];

    const PylithScalar traceStrainTpdt = e11 + e22 + e33;
    const PylithScalar s123 = lambda * traceStrainTpdt;

    stress[0] = s123 + mu2 * e11 + initialStress[0];
    stress[1] = s123 + mu2 * e22 + initialStress[1];
    stress[2] = s123 + mu2 * e33 + initialStress[2];
    stress[3] = mu2 * e12 + initialStress[3];
    stress[4] = mu2 * e23 + initialStress[4];
    stress[5] = mu2 * e13 + initialStress[5];

    PetscLogFlops(31);

  } // if/else
#if 0 // DEBUGGING
  const PylithScalar alphaYield = properties[p_alphaYield];
  const PylithScalar beta = properties[p_beta];
  const PylithScalar alphaFlow = properties[p_alphaFlow];
  const PylithScalar meanStressTest = (stress[0] + stress[1] + stress[2])/3.0;
  const PylithScalar devStressTest[tensorSize] = {
    stress[0] - meanStressTest,
    stress[1] - meanStressTest,
    stress[2] - meanStressTest,
    stress[3],
    stress[4],
    stress[5]
  };
  const PylithScalar stressInvar2Test =
    sqrt(0.5 * scalarProduct3D(devStressTest, devStressTest));
  
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
pylith::materials::DruckerPrager3D::_calcElasticConstsElastic(
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
  assert(_DruckerPrager3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_DruckerPrager3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_DruckerPrager3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_DruckerPrager3D::tensorSize == initialStrainSize);
 
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar lambda2mu = lambda + mu2;

  elasticConsts[ 0] = lambda2mu; // C1111
  elasticConsts[ 1] = lambda; // C1122
  elasticConsts[ 2] = lambda; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = 0; // C1123
  elasticConsts[ 5] = 0; // C1113
  elasticConsts[ 6] = lambda; // C2211
  elasticConsts[ 7] = lambda2mu; // C2222
  elasticConsts[ 8] = lambda; // C2233
  elasticConsts[ 9] = 0; // C2212
  elasticConsts[10] = 0; // C2223
  elasticConsts[11] = 0; // C2213
  elasticConsts[12] = lambda; // C3311
  elasticConsts[13] = lambda; // C3322
  elasticConsts[14] = lambda2mu; // C3333
  elasticConsts[15] = 0; // C3312
  elasticConsts[16] = 0; // C3323
  elasticConsts[17] = 0; // C3313
  elasticConsts[18] = 0; // C1211
  elasticConsts[19] = 0; // C1222
  elasticConsts[20] = 0; // C1233
  elasticConsts[21] = mu2; // C1212
  elasticConsts[22] = 0; // C1223
  elasticConsts[23] = 0; // C1213
  elasticConsts[24] = 0; // C2311
  elasticConsts[25] = 0; // C2322
  elasticConsts[26] = 0; // C2333
  elasticConsts[27] = 0; // C2312
  elasticConsts[28] = mu2; // C2323
  elasticConsts[29] = 0; // C2313
  elasticConsts[30] = 0; // C1311
  elasticConsts[31] = 0; // C1322
  elasticConsts[32] = 0; // C1333
  elasticConsts[33] = 0; // C1312
  elasticConsts[34] = 0; // C1323
  elasticConsts[35] = mu2; // C1313

  PetscLogFlops(2);
} // _calcElasticConstsElastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties
// as an elastoplastic material.
void
pylith::materials::DruckerPrager3D::_calcElasticConstsElastoplastic(
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
  assert(0 != elasticConsts);
  assert(_DruckerPrager3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_DruckerPrager3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_DruckerPrager3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_DruckerPrager3D::tensorSize == initialStrainSize);

  // Duplicate functionality of _calcStressElastoplastic
  // Get properties
  const int tensorSize = 6;
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar alphaYield = properties[p_alphaYield];
  const PylithScalar beta = properties[p_beta];
  const PylithScalar alphaFlow = properties[p_alphaFlow];
  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2/3.0;
  const PylithScalar ae = 1.0/mu2;
  const PylithScalar am = 1.0/(3.0 * bulkModulus);
  
  // Get state variables from previous time step
  const PylithScalar plasticStrainT[tensorSize] = {
    stateVars[s_plasticStrain  ],
    stateVars[s_plasticStrain+1],
    stateVars[s_plasticStrain+2],
    stateVars[s_plasticStrain+3],
    stateVars[s_plasticStrain+4],
    stateVars[s_plasticStrain+5],
  };
  const PylithScalar meanPlasticStrainT = (plasticStrainT[0] +
					   plasticStrainT[1] +
					   plasticStrainT[2])/3.0;
  PylithScalar devPlasticStrainT[tensorSize];
  calcDeviatoric3D(devPlasticStrainT, plasticStrainT, meanPlasticStrainT);
  
  const PylithScalar diag[tensorSize] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
  
  // Initial stress values
  const PylithScalar meanStressInitial = (initialStress[0] +
					  initialStress[1] +
					  initialStress[2])/3.0;
  PylithScalar devStressInitial[tensorSize];
  calcDeviatoric3D(devStressInitial, initialStress, meanStressInitial);
  
  // Initial strain values
  const PylithScalar meanStrainInitial = (initialStrain[0] +
					  initialStrain[1] +
					  initialStrain[2])/3.0;
  PylithScalar devStrainInitial[tensorSize];
  calcDeviatoric3D(devStrainInitial, initialStrain, meanStrainInitial);

  // Values for current time step
  const PylithScalar meanStrainTpdt = (totalStrain[0] +
				       totalStrain[1] +
				       totalStrain[2])/3.0;
  const PylithScalar meanStrainPPTpdt = meanStrainTpdt - meanPlasticStrainT -
    meanStrainInitial;
  
  const PylithScalar strainPPTpdt[tensorSize] = {
    totalStrain[0] - meanStrainTpdt - devPlasticStrainT[0] -
    devStrainInitial[0],
    totalStrain[1] - meanStrainTpdt - devPlasticStrainT[1] -
    devStrainInitial[1],
    totalStrain[2] - meanStrainTpdt - devPlasticStrainT[2] -
    devStrainInitial[2],
    totalStrain[3] - devPlasticStrainT[3] - devStrainInitial[3],
    totalStrain[4] - devPlasticStrainT[4] - devStrainInitial[4],
    totalStrain[5] - devPlasticStrainT[5] - devStrainInitial[5],
  };
  
  // Compute trial elastic stresses and yield function to see if yield should
  // occur.
  const PylithScalar trialDevStress[tensorSize] = {
    strainPPTpdt[0]/ae + devStressInitial[0],
    strainPPTpdt[1]/ae + devStressInitial[1],
    strainPPTpdt[2]/ae + devStressInitial[2],
    strainPPTpdt[3]/ae + devStressInitial[3],
    strainPPTpdt[4]/ae + devStressInitial[4],
    strainPPTpdt[5]/ae + devStressInitial[5]
  };
  const PylithScalar trialMeanStress = meanStrainPPTpdt/am + meanStressInitial;
  const PylithScalar stressInvar2 =
    sqrt(0.5 * scalarProduct3D(trialDevStress, trialDevStress));
  const PylithScalar yieldFunction = 3.0 * alphaYield * trialMeanStress +
    stressInvar2 - beta;
#if 0 // DEBUGGING
  std::cout << "Function _calcElasticConstsElastoPlastic:" << std::endl;
  std::cout << "  alphaYield:       " << alphaYield << std::endl;
  std::cout << "  beta:             " << beta << std::endl;
  std::cout << "  trialMeanStress:  " << trialMeanStress << std::endl;
  std::cout << "  stressInvar2:     " << stressInvar2 << std::endl;
  std::cout << "  yieldFunction:    " << yieldFunction << std::endl;
#endif
  PetscLogFlops(76);
  
  // If yield function is greater than zero, compute elastoplastic stress and
  // corresponding tangent matrix.
  if (yieldFunction >= 0.0) {
    const PylithScalar devStressInitialProd = 
      scalarProduct3D(devStressInitial, devStressInitial);
    const PylithScalar strainPPTpdtProd =
      scalarProduct3D(strainPPTpdt, strainPPTpdt);
    const PylithScalar d = 
      sqrt(ae * ae * devStressInitialProd + 2.0 * ae *
	   scalarProduct3D(devStressInitial, strainPPTpdt) + strainPPTpdtProd);
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
      tensileYield = (sqrt(2.0)*d < testMult) ? true : false;
      plasticMult = (tensileYield) ? sqrt(2.0)*d : testMult;
      dFac2 = (d > 0.0) ? 1.0/(sqrt(2.0) * d) : 0.0;
    } else {
      plasticMult = plasticFac *
	(meanStrainFac * (meanStrainPPTpdt/am + meanStressInitial) +
	 dFac * d - beta);
      dFac2 = 1.0/(sqrt(2.0) * d);
    } // if/else

    // Define some constants, vectors, and matrices.
    const PylithScalar third = 1.0/3.0;
    const PylithScalar dEdEpsilon[6][6] = {
      { 2.0 * third,      -third,      -third, 0.0, 0.0, 0.0},
      {      -third, 2.0 * third,      -third, 0.0, 0.0, 0.0},
      {      -third,      -third, 2.0 * third, 0.0, 0.0, 0.0},
      {         0.0,         0.0,         0.0, 1.0, 0.0, 0.0},
      {         0.0,         0.0,         0.0, 0.0, 1.0, 0.0},
      {         0.0,         0.0,         0.0, 0.0, 0.0, 1.0}};
    const PylithScalar vec1[tensorSize] = {
      strainPPTpdt[0] + ae * devStressInitial[0],
      strainPPTpdt[1] + ae * devStressInitial[1],
      strainPPTpdt[2] + ae * devStressInitial[2],
      strainPPTpdt[3] + ae * devStressInitial[3],
      strainPPTpdt[4] + ae * devStressInitial[4],
      strainPPTpdt[5] + ae * devStressInitial[5]
    };
    PylithScalar dDeltaEdEpsilon = 0.0;

    // Compute elasticity matrix.
    if (d > 0.0) {
      const PylithScalar dDdEpsilon[tensorSize] = {vec1[0]/d,
						   vec1[1]/d,
						   vec1[2]/d,
						   2.0 * vec1[3]/d,
						   2.0 * vec1[4]/d,
						   2.0 * vec1[5]/d};
      PylithScalar dLambdadEpsilon[tensorSize];
      if (tensileYield) {
	dLambdadEpsilon[0] = sqrt(2.0) * dDdEpsilon[0];
	dLambdadEpsilon[1] = sqrt(2.0) * dDdEpsilon[1];
	dLambdadEpsilon[2] = sqrt(2.0) * dDdEpsilon[2];
	dLambdadEpsilon[3] = sqrt(2.0) * dDdEpsilon[3];
	dLambdadEpsilon[4] = sqrt(2.0) * dDdEpsilon[4];
	dLambdadEpsilon[5] = sqrt(2.0) * dDdEpsilon[5];
      } else {
	dLambdadEpsilon[0] = plasticFac *
	  (alphaYield/am + dFac * dDdEpsilon[0]);
	dLambdadEpsilon[1] = plasticFac *
	  (alphaYield/am + dFac * dDdEpsilon[1]);
	dLambdadEpsilon[2] = plasticFac *
	  (alphaYield/am + dFac * dDdEpsilon[2]);
	dLambdadEpsilon[3] = plasticFac * dFac * dDdEpsilon[3];
	dLambdadEpsilon[4] = plasticFac * dFac * dDdEpsilon[4];
	dLambdadEpsilon[5] = plasticFac * dFac * dDdEpsilon[5];
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
      const PylithScalar dLambdadEpsilon[tensorSize] = {0.0,
							0.0,	
							0.0,	
							0.0,
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

    PetscLogFlops(109 + tensorSize * tensorSize * 15);

  } else {
    // No plastic strain.
    const PylithScalar lambda2mu = lambda + mu2;
    elasticConsts[ 0] = lambda2mu; // C1111
    elasticConsts[ 1] = lambda; // C1122
    elasticConsts[ 2] = lambda; // C1133
    elasticConsts[ 3] = 0; // C1112
    elasticConsts[ 4] = 0; // C1123
    elasticConsts[ 5] = 0; // C1113
    elasticConsts[ 6] = lambda; // C2211
    elasticConsts[ 7] = lambda2mu; // C2222
    elasticConsts[ 8] = lambda; // C2233
    elasticConsts[ 9] = 0; // C2212
    elasticConsts[10] = 0; // C2223
    elasticConsts[11] = 0; // C2213
    elasticConsts[12] = lambda; // C3311
    elasticConsts[13] = lambda; // C3322
    elasticConsts[14] = lambda2mu; // C3333
    elasticConsts[15] = 0; // C3312
    elasticConsts[16] = 0; // C3323
    elasticConsts[17] = 0; // C3313
    elasticConsts[18] = 0; // C1211
    elasticConsts[19] = 0; // C1222
    elasticConsts[20] = 0; // C1233
    elasticConsts[21] = mu2; // C1212
    elasticConsts[22] = 0; // C1223
    elasticConsts[23] = 0; // C1213
    elasticConsts[24] = 0; // C2311
    elasticConsts[25] = 0; // C2322
    elasticConsts[26] = 0; // C2333
    elasticConsts[27] = 0; // C2312
    elasticConsts[28] = mu2; // C2323
    elasticConsts[29] = 0; // C2313
    elasticConsts[30] = 0; // C1311
    elasticConsts[31] = 0; // C1322
    elasticConsts[32] = 0; // C1333
    elasticConsts[33] = 0; // C1312
    elasticConsts[34] = 0; // C1323
    elasticConsts[35] = mu2; // C1313

    PetscLogFlops(1);
  } // else

} // _calcElasticConstsElastoplastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::DruckerPrager3D::_updateStateVarsElastic(
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
  assert(_DruckerPrager3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_DruckerPrager3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_DruckerPrager3D::tensorSize == initialStrainSize);

  for (int iComp=0; iComp < _tensorSize; ++iComp) {
    stateVars[s_plasticStrain+iComp] = 0.0;
  } // for

  _needNewJacobian = true;
} // _updateStateVarsElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::DruckerPrager3D::_updateStateVarsElastoplastic(
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
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_DruckerPrager3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_DruckerPrager3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_DruckerPrager3D::tensorSize == initialStrainSize);

  // For now, we are duplicating the functionality of _calcStressElastoplastic,
  // since otherwise we would have to redo a lot of calculations.

  const int tensorSize = 6;
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar alphaYield = properties[p_alphaYield];
  const PylithScalar beta = properties[p_beta];
  const PylithScalar alphaFlow = properties[p_alphaFlow];
  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2/3.0;
  const PylithScalar ae = 1.0/mu2;
  const PylithScalar am = 1.0/(3.0 * bulkModulus);

  const PylithScalar plasticStrainT[tensorSize] = {
    stateVars[s_plasticStrain],
    stateVars[s_plasticStrain + 1],
    stateVars[s_plasticStrain + 2],
    stateVars[s_plasticStrain + 3],
    stateVars[s_plasticStrain + 4],
    stateVars[s_plasticStrain + 5]
  };
  const PylithScalar meanPlasticStrainT = (plasticStrainT[0] +
					   plasticStrainT[1] +
					   plasticStrainT[2])/3.0;
  PylithScalar devPlasticStrainT[tensorSize];
  calcDeviatoric3D(devPlasticStrainT, plasticStrainT, meanPlasticStrainT);

  const PylithScalar diag[tensorSize] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Initial stress values
  const PylithScalar meanStressInitial = (initialStress[0] +
					  initialStress[1] +
					  initialStress[2])/3.0;
  PylithScalar devStressInitial[tensorSize];
  calcDeviatoric3D(devStressInitial, initialStress, meanStressInitial);

  // Initial strain values
  const PylithScalar meanStrainInitial = (initialStrain[0] +
					  initialStrain[1] +
					  initialStrain[2])/3.0;
  PylithScalar devStrainInitial[tensorSize];
  calcDeviatoric3D(devStrainInitial, initialStrain, meanStrainInitial);

  // Values for current time step
  const PylithScalar e11 = totalStrain[0];
  const PylithScalar e22 = totalStrain[1];
  const PylithScalar e33 = totalStrain[2];
  const PylithScalar meanStrainTpdt = (e11 + e22 + e33)/3.0;
  const PylithScalar meanStrainPPTpdt = meanStrainTpdt - meanPlasticStrainT -
    meanStrainInitial;

  const PylithScalar strainPPTpdt[tensorSize] = {
    totalStrain[0] - meanStrainTpdt - devPlasticStrainT[0] -
    devStrainInitial[0],
    totalStrain[1] - meanStrainTpdt - devPlasticStrainT[1] -
    devStrainInitial[1],
    totalStrain[2] - meanStrainTpdt - devPlasticStrainT[2] -
    devStrainInitial[2],
    totalStrain[3] - devPlasticStrainT[3] - devStrainInitial[3],
    totalStrain[4] - devPlasticStrainT[4] - devStrainInitial[4],
    totalStrain[5] - devPlasticStrainT[5] - devStrainInitial[5]
  };

  // Compute trial elastic stresses and yield function to see if yield should
  // occur.
  const PylithScalar trialDevStress[tensorSize] = {
    strainPPTpdt[0]/ae + devStressInitial[0],
    strainPPTpdt[1]/ae + devStressInitial[1],
    strainPPTpdt[2]/ae + devStressInitial[2],
    strainPPTpdt[3]/ae + devStressInitial[3],
    strainPPTpdt[4]/ae + devStressInitial[4],
    strainPPTpdt[5]/ae + devStressInitial[5]
  };
  const PylithScalar trialMeanStress = meanStrainPPTpdt/am + meanStressInitial;
  const PylithScalar stressInvar2 =
    sqrt(0.5 * scalarProduct3D(trialDevStress, trialDevStress));
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
  PetscLogFlops(76);

  // If yield function is greater than zero, compute plastic strains.
  // Otherwise, plastic strains remain the same.
  if (yieldFunction >= 0.0) {
    const PylithScalar devStressInitialProd = 
      scalarProduct3D(devStressInitial, devStressInitial);
    const PylithScalar strainPPTpdtProd =
      scalarProduct3D(strainPPTpdt, strainPPTpdt);
    const PylithScalar d =
      sqrt(ae * ae * devStressInitialProd + 2.0 * ae *
	   scalarProduct3D(devStressInitial, strainPPTpdt) + strainPPTpdtProd);
    const PylithScalar plasticFac = 2.0 * ae * am/
      (6.0 * alphaYield * alphaFlow * ae + am);
    const PylithScalar meanStrainFac = 3.0 * alphaYield;
    const PylithScalar dFac = 1.0/(sqrt(2.0) * ae);

    PylithScalar plasticMult = 0.0;
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
      for (int iComp=0; iComp < tensorSize; ++iComp) {
	deltaDevPlasticStrain = plasticMult *(strainPPTpdt[iComp] +
					      ae * devStressInitial[iComp])/
	  (sqrt(2.0) * d);
	stateVars[s_plasticStrain+iComp] += deltaDevPlasticStrain +
	  diag[iComp] * deltaMeanPlasticStrain;
      } // for
    } else {
      for (int iComp=0; iComp < tensorSize; ++iComp) {
	stateVars[s_plasticStrain+iComp] +=
	  diag[iComp] * deltaMeanPlasticStrain;
      } // for
    } // if/else

    PetscLogFlops(60 + 9 * tensorSize);

  } // if

  _needNewJacobian = true;

} // _updateStateVarsElastoplastic

// End of file 
