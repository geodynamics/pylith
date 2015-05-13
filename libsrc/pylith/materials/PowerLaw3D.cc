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

#include "PowerLaw3D.hh" // implementation of object methods

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
    namespace _PowerLaw3D{

      /// Dimension of material.
      const int dimension = 3;

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 6;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 36;

      /// Number of physical properties.
      const int numProperties = 6;

      /// Physical properties.
      const Metadata::ParamDescription properties[] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
	{ "reference_strain_rate", 1, pylith::topology::FieldBase::SCALAR },
	{ "reference_stress", 1, pylith::topology::FieldBase::SCALAR },
	{ "power_law_exponent", 1, pylith::topology::FieldBase::SCALAR }
      };

      // Values expected in properties spatial database
      const int numDBProperties = 6;
      const char* dbProperties[] = {"density", "vs", "vp" ,
				    "reference-strain-rate",
				    "reference-stress",
				    "power-law-exponent"};

      /// Number of state variables.
      const int numStateVars = 2;

      /// State variables.
      const Metadata::ParamDescription stateVars[] = {
	{ "viscous_strain", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "stress", tensorSize, pylith::topology::FieldBase::TENSOR }
      };

      // Values expected in state variables spatial database.
      const int numDBStateVars = 12;
      const char* dbStateVars[] = { "viscous-strain-xx",
				    "viscous-strain-yy",
				    "viscous-strain-zz",
				    "viscous-strain-xy",
				    "viscous-strain-yz",
				    "viscous-strain-xz",
				    "stress-xx",
				    "stress-yy",
				    "stress-zz",
				    "stress-xy",
				    "stress-yz",
				    "stress-xz"
      };

    } // _PowerLaw3D
  } // materials
} // pylith

// Indices of physical properties.
const int pylith::materials::PowerLaw3D::p_density = 0;

const int pylith::materials::PowerLaw3D::p_mu = 
  pylith::materials::PowerLaw3D::p_density + 1;

const int pylith::materials::PowerLaw3D::p_lambda = 
  pylith::materials::PowerLaw3D::p_mu + 1;

const int pylith::materials::PowerLaw3D::p_referenceStrainRate = 
  pylith::materials::PowerLaw3D::p_lambda + 1;

const int pylith::materials::PowerLaw3D::p_referenceStress = 
  pylith::materials::PowerLaw3D::p_referenceStrainRate + 1;

const int pylith::materials::PowerLaw3D::p_powerLawExponent = 
  pylith::materials::PowerLaw3D::p_referenceStress + 1;

// Indices of property database values (order must match dbProperties).
const int pylith::materials::PowerLaw3D::db_density = 0;

const int pylith::materials::PowerLaw3D::db_vs = 
  pylith::materials::PowerLaw3D::db_density + 1;

const int pylith::materials::PowerLaw3D::db_vp = 
  pylith::materials::PowerLaw3D::db_vs + 1;

const int pylith::materials::PowerLaw3D::db_referenceStrainRate = 
  pylith::materials::PowerLaw3D::db_vp + 1;

const int pylith::materials::PowerLaw3D::db_referenceStress = 
  pylith::materials::PowerLaw3D::db_referenceStrainRate + 1;

const int pylith::materials::PowerLaw3D::db_powerLawExponent = 
  pylith::materials::PowerLaw3D::db_referenceStress + 1;

// Indices of state variables.
const int pylith::materials::PowerLaw3D::s_viscousStrain = 0;

const int pylith::materials::PowerLaw3D::s_stress = 
  pylith::materials::PowerLaw3D::s_viscousStrain + 
  _PowerLaw3D::tensorSize;

// Indices of state variable database values (order must match dbStateVars).
const int pylith::materials::PowerLaw3D::db_viscousStrain = 0;

const int pylith::materials::PowerLaw3D::db_stress = 
  pylith::materials::PowerLaw3D::db_viscousStrain +
  _PowerLaw3D::tensorSize;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::PowerLaw3D::PowerLaw3D(void) :
  ElasticMaterial(_PowerLaw3D::dimension,
		  _PowerLaw3D::tensorSize,
		  _PowerLaw3D::numElasticConsts,
		  Metadata(_PowerLaw3D::properties,
			   _PowerLaw3D::numProperties,
			   _PowerLaw3D::dbProperties,
			   _PowerLaw3D::numDBProperties,
			   _PowerLaw3D::stateVars,
			   _PowerLaw3D::numStateVars,
			   _PowerLaw3D::dbStateVars,
			   _PowerLaw3D::numDBStateVars)),
  _calcElasticConstsFn(0),
  _calcStressFn(0),
  _updateStateVarsFn(0)
{ // constructor
  useElasticBehavior(false);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::PowerLaw3D::~PowerLaw3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set whether elastic or inelastic constitutive relations are used.
void
pylith::materials::PowerLaw3D::useElasticBehavior(const bool flag)
{ // useElasticBehavior
  if (flag) {
    _calcStressFn = 
      &pylith::materials::PowerLaw3D::_calcStressElastic;
    _calcElasticConstsFn = 
      &pylith::materials::PowerLaw3D::_calcElasticConstsElastic;
    _updateStateVarsFn = 
      &pylith::materials::PowerLaw3D::_updateStateVarsElastic;

  } else {
    _calcStressFn = 
      &pylith::materials::PowerLaw3D::_calcStressViscoelastic;
    _calcElasticConstsFn = 
      &pylith::materials::PowerLaw3D::_calcElasticConstsViscoelastic;
    _updateStateVarsFn = 
      &pylith::materials::PowerLaw3D::_updateStateVarsViscoelastic;
  } // if/else
} // useElasticBehavior

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::PowerLaw3D::_dbToProperties(
				PylithScalar* const propValues,
				const scalar_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_PowerLaw3D::numDBProperties == numDBValues);

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

  const PylithScalar mu = density * vs*vs;
  const PylithScalar lambda = density * vp*vp - 2.0*mu;

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
pylith::materials::PowerLaw3D::_nondimProperties(PylithScalar* const values,
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
pylith::materials::PowerLaw3D::_dimProperties(PylithScalar* const values,
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
pylith::materials::PowerLaw3D::_dbToStateVars(
				PylithScalar* const stateValues,
				const scalar_array& dbValues)
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_PowerLaw3D::numDBStateVars == numDBValues);

  const int totalSize = 2 * _tensorSize;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  memcpy(&stateValues[s_viscousStrain], &dbValues[db_viscousStrain],
	 _tensorSize*sizeof(PylithScalar));
  memcpy(&stateValues[s_stress], &dbValues[db_stress],
	 _tensorSize*sizeof(PylithScalar));

  PetscLogFlops(0);
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::materials::PowerLaw3D::_nondimStateVars(PylithScalar* const values,
						const int nvalues) const
{ // _nondimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numVarsQuadPt);

  const PylithScalar pressureScale = _normalizer->pressureScale();
  _normalizer->nondimensionalize(&values[s_stress], _tensorSize, pressureScale);

  PetscLogFlops(_tensorSize);
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::materials::PowerLaw3D::_dimStateVars(PylithScalar* const values,
					     const int nvalues) const
{ // _dimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numVarsQuadPt);

  const PylithScalar pressureScale = _normalizer->pressureScale();
  _normalizer->dimensionalize(&values[s_stress], _tensorSize, pressureScale);

  PetscLogFlops(_tensorSize);
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::PowerLaw3D::_calcDensity(PylithScalar* const density,
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
pylith::materials::PowerLaw3D::_stableTimeStepImplicit(
				  const PylithScalar* properties,
				  const int numProperties,
				  const PylithScalar* stateVars,
				  const int numStateVars) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  const PylithScalar mu = properties[p_mu];
  const PylithScalar referenceStrainRate = properties[p_referenceStrainRate];
  const PylithScalar referenceStress = properties[p_referenceStress];
  const PylithScalar powerLawExp = properties[p_powerLawExponent];

  const PylithScalar stress[] = {stateVars[s_stress],
			   stateVars[s_stress + 1],
			   stateVars[s_stress + 2],
			   stateVars[s_stress + 3],
			   stateVars[s_stress + 4],
			   stateVars[s_stress + 5]};
  const PylithScalar meanStress = (stress[0] + stress[1] + stress[2])/3.0;
  const PylithScalar devStress[] = {stress[0] - meanStress,
			      stress[1] - meanStress,
			      stress[2] - meanStress,
			      stress[3],
			      stress[4],
			      stress[5] };
  const PylithScalar devStressProd = scalarProduct3D(devStress, devStress);
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
  PetscLogFlops(27);
  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Get stable time step for explicit time integration.
PylithScalar
pylith::materials::PowerLaw3D::_stableTimeStepExplicit(const PylithScalar* properties,
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
pylith::materials::PowerLaw3D::_calcStressElastic(
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
  assert(_PowerLaw3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLaw3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLaw3D::tensorSize == initialStrainSize);

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
// Compute stress tensor at location from properties as a viscoelastic
// material.
void
pylith::materials::PowerLaw3D::_calcStressViscoelastic(
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
  assert(_PowerLaw3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLaw3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLaw3D::tensorSize == initialStrainSize);

  const int tensorSize = 6;
  assert(_tensorSize == tensorSize);
    
  // We need to do root-finding method if state variables are from previous
  // time step.
  if (computeStateVars) {

    const PylithScalar mu = properties[p_mu];
    const PylithScalar lambda = properties[p_lambda];
    const PylithScalar referenceStrainRate = properties[p_referenceStrainRate];
    const PylithScalar referenceStress = properties[p_referenceStress];
    const PylithScalar powerLawExp = properties[p_powerLawExponent];
    const PylithScalar visStrainT[tensorSize] = {
      stateVars[s_viscousStrain],
      stateVars[s_viscousStrain + 1],
      stateVars[s_viscousStrain + 2],
      stateVars[s_viscousStrain + 3],
      stateVars[s_viscousStrain + 4],
      stateVars[s_viscousStrain + 5]
    };
    const PylithScalar stressT[tensorSize] = {stateVars[s_stress],
					      stateVars[s_stress + 1],
					      stateVars[s_stress + 2],
					      stateVars[s_stress + 3],
					      stateVars[s_stress + 4],
					      stateVars[s_stress + 5]};

    const PylithScalar mu2 = 2.0 * mu;
    const PylithScalar bulkModulus = lambda + mu2/3.0;
    const PylithScalar ae = 1.0/mu2;
    const PylithScalar diag[tensorSize] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

    // Need to figure out how time integration parameter alpha is going to be
    // specified.  It should probably be specified in the problem definition and
    // then used only by the material types that use it.  For now we are setting
    // it to 0.5, which should probably be the default value.
    const PylithScalar alpha = 0.5;
    const PylithScalar timeFac = _dt * (1.0 - alpha);

    // Initial stress values
    const PylithScalar meanStressInitial = (initialStress[0] +
					    initialStress[1] +
					    initialStress[2])/3.0;
    const PylithScalar devStressInitial[tensorSize] = {
      initialStress[0] - meanStressInitial,
      initialStress[1] - meanStressInitial,
      initialStress[2] - meanStressInitial,
      initialStress[3],
      initialStress[4],
      initialStress[5]
    };
    const PylithScalar stressInvar2Initial =
      0.5 * scalarProduct3D(devStressInitial, devStressInitial);

    // Initial strain values
    const PylithScalar meanStrainInitial = (initialStrain[0] +
					    initialStrain[1] +
					    initialStrain[2])/3.0;

    // Values for current time step
    const PylithScalar e11 = totalStrain[0];
    const PylithScalar e22 = totalStrain[1];
    const PylithScalar e33 = totalStrain[2];
    const PylithScalar meanStrainTpdt = (e11 + e22 + e33)/3.0 -
      meanStrainInitial;
    const PylithScalar meanStressTpdt = 3.0 * bulkModulus * meanStrainTpdt;

    // Note that I use the initial strain rather than the deviatoric initial
    // strain since otherwise the initial mean strain would get used twice.
    const PylithScalar strainPPTpdt[tensorSize] = {
      totalStrain[0] - meanStrainTpdt - visStrainT[0] - initialStrain[0],
      totalStrain[1] - meanStrainTpdt - visStrainT[1] - initialStrain[1],
      totalStrain[2] - meanStrainTpdt - visStrainT[2] - initialStrain[2],
      totalStrain[3] - visStrainT[3] - initialStrain[3],
      totalStrain[4] - visStrainT[4] - initialStrain[4],
      totalStrain[5] - visStrainT[5] - initialStrain[5]
    };
    const PylithScalar strainPPInvar2Tpdt =
      0.5 * scalarProduct3D(strainPPTpdt, strainPPTpdt);

    // Values for previous time step
    const PylithScalar meanStressT = (stressT[0] +
				      stressT[1] +
				      stressT[2])/3.0;
    const PylithScalar devStressT[tensorSize] = { stressT[0] - meanStressT,
						  stressT[1] - meanStressT,
						  stressT[2] - meanStressT,
						  stressT[3],
						  stressT[4],
						  stressT[5] };
    const PylithScalar stressInvar2T =
      0.5 * scalarProduct3D(devStressT, devStressT);
    const PylithScalar effStressT = sqrt(stressInvar2T);

    // Finish defining parameters needed for root-finding algorithm.
    const PylithScalar b = strainPPInvar2Tpdt + ae *
      scalarProduct3D(strainPPTpdt, devStressInitial) +
      ae * ae * stressInvar2Initial;
    const PylithScalar c =
      (scalarProduct3D(strainPPTpdt, devStressT) +
       ae * scalarProduct3D(devStressT, devStressInitial)) * timeFac;
    const PylithScalar d = timeFac * effStressT;

    PetscLogFlops(92);

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
	EffectiveStress::calculate<PowerLaw3D>(effStressInitialGuess,
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

    for (int iComp=0; iComp < tensorSize; ++iComp) {
      devStressTpdt = factor1 *
	(strainPPTpdt[iComp] - factor2 * devStressT[iComp] +
	 ae * devStressInitial[iComp]);
      stress[iComp] = devStressTpdt + diag[iComp] *
	(meanStressTpdt + meanStressInitial);
    } // for
    PetscLogFlops(14 + 8 * tensorSize);

    // If state variables have already been updated, current stress is already
    // contained in stress.
  } else {
    memcpy(&stress[0], &stateVars[s_stress], tensorSize * sizeof(PylithScalar));
  } // else

} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Effective stress function that computes effective stress function only
// (no derivative).
PylithScalar
pylith::materials::PowerLaw3D::effStressFunc(const PylithScalar effStressTpdt)
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

  PetscLogFlops(21);

  return y;
} // effStressFunc

// ----------------------------------------------------------------------
// Effective stress function that computes effective stress function
// derivative only (no function value).
PylithScalar
pylith::materials::PowerLaw3D::effStressDerivFunc(const PylithScalar effStressTpdt)
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
  const PylithScalar dGammaTau = referenceStrainRate * alpha *
    (powerLawExp - 1.0) *
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
pylith::materials::PowerLaw3D::effStressFuncDerivFunc(PylithScalar* func,
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
  const PylithScalar effStressTau = factor1 * effStressT +
    alpha * effStressTpdt;
  const PylithScalar gammaTau = referenceStrainRate *
    pow((effStressTau/referenceStress), (powerLawExp - 1.0))/referenceStress;
  const PylithScalar dGammaTau = referenceStrainRate * alpha *
    (powerLawExp - 1.0) *
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
pylith::materials::PowerLaw3D::_calcElasticConstsElastic(
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
  assert(_PowerLaw3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLaw3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLaw3D::tensorSize == initialStrainSize);
 
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
// as a viscoelastic material.
void
pylith::materials::PowerLaw3D::_calcElasticConstsViscoelastic(
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
  assert(_PowerLaw3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLaw3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLaw3D::tensorSize == initialStrainSize);

  const int tensorSize = 6;
  assert(_tensorSize == tensorSize);

  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar referenceStrainRate = properties[p_referenceStrainRate];
  const PylithScalar referenceStress = properties[p_referenceStress];
  const PylithScalar powerLawExp = properties[p_powerLawExponent];
    
  // State variables.
  const PylithScalar visStrainT[tensorSize] = {stateVars[s_viscousStrain],
					       stateVars[s_viscousStrain + 1],
					       stateVars[s_viscousStrain + 2],
					       stateVars[s_viscousStrain + 3],
					       stateVars[s_viscousStrain + 4],
					       stateVars[s_viscousStrain + 5]};
  const PylithScalar stressT[tensorSize] = {stateVars[s_stress],
					    stateVars[s_stress + 1],
					    stateVars[s_stress + 2],
					    stateVars[s_stress + 3],
					    stateVars[s_stress + 4],
					    stateVars[s_stress + 5]};

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
					  initialStress[2])/3.0;
  const PylithScalar devStressInitial[tensorSize] = {
    initialStress[0] - meanStressInitial,
    initialStress[1] - meanStressInitial,
    initialStress[2] - meanStressInitial,
    initialStress[3],
    initialStress[4],
    initialStress[5]
  };
  const PylithScalar stressInvar2Initial =
    0.5 * scalarProduct3D(devStressInitial, devStressInitial);

  // Initial strain values.
  const PylithScalar meanStrainInitial = (initialStrain[0] +
					  initialStrain[1] +
					  initialStrain[2])/3.0;
  
  /// Values for current time step
  const PylithScalar e11 = totalStrain[0];
  const PylithScalar e22 = totalStrain[1];
  const PylithScalar e33 = totalStrain[2];
  const PylithScalar meanStrainTpdt = (e11 + e22 + e33)/3.0 - meanStrainInitial;
  
  // Note that I use the initial strain rather than the deviatoric initial
  // strain since otherwise the initial mean strain would get used twice.
  
  const PylithScalar strainPPTpdt[tensorSize] = {
    totalStrain[0] - meanStrainTpdt - visStrainT[0] - initialStrain[0],
    totalStrain[1] - meanStrainTpdt - visStrainT[1] - initialStrain[1],
    totalStrain[2] - meanStrainTpdt - visStrainT[2] - initialStrain[2],
    totalStrain[3] - visStrainT[3] - initialStrain[3],
    totalStrain[4] - visStrainT[4] - initialStrain[4],
    totalStrain[5] - visStrainT[5] - initialStrain[5]
  };
  const PylithScalar strainPPInvar2Tpdt =
    0.5 * scalarProduct3D(strainPPTpdt, strainPPTpdt);
  
  // Values for previous time step
  const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
  const PylithScalar devStressT[tensorSize] = { stressT[0] - meanStressT,
						stressT[1] - meanStressT,
						stressT[2] - meanStressT,
						stressT[3],
						stressT[4],
						stressT[5] };
  const PylithScalar stressInvar2T =
    0.5 * scalarProduct3D(devStressT, devStressT);
  const PylithScalar effStressT = sqrt(stressInvar2T);
    
  // Finish defining parameters needed for root-finding algorithm.
  const PylithScalar b = strainPPInvar2Tpdt +
    ae * scalarProduct3D(strainPPTpdt, devStressInitial) +
    ae * ae * stressInvar2Initial;
  const PylithScalar c =
    (scalarProduct3D(strainPPTpdt, devStressT) +
     ae * scalarProduct3D(devStressT, devStressInitial)) * timeFac;
  const PylithScalar d = timeFac * effStressT;

  PetscLogFlops(92);

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
      EffectiveStress::calculate<PowerLaw3D>(effStressInitialGuess,
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
      (strainPPTpdt[2] - factor2 * devStressT[2] + ae * devStressInitial[2]),
      factor1 *
      (strainPPTpdt[3] - factor2 * devStressT[3] + ae * devStressInitial[3]),
      factor1 *
      (strainPPTpdt[4] - factor2 * devStressT[4] + ae * devStressInitial[4]),
      factor1 *
      (strainPPTpdt[5] - factor2 * devStressT[5] + ae * devStressInitial[5])
    };
    const PylithScalar devStressTau[tensorSize] = {
      alpha * devStressT[0] + explicitFac * devStressTpdt[0],
      alpha * devStressT[1] + explicitFac * devStressTpdt[1],
      alpha * devStressT[2] + explicitFac * devStressTpdt[2],
      alpha * devStressT[3] + explicitFac * devStressTpdt[3],
      alpha * devStressT[4] + explicitFac * devStressTpdt[4],
      alpha * devStressT[5] + explicitFac * devStressTpdt[5]
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
    const PylithScalar dStress33dStrain33 = 1.0/
      (a + devStressTau[2] * devStressTpdt[2] * factor3);
    const PylithScalar dStress12dStrain12 = 1.0/
      (a + 2.0 * devStressTau[3] * devStressTpdt[3] * factor3);
    const PylithScalar dStress23dStrain23 = 1.0/
      (a + 2.0 * devStressTau[4] * devStressTpdt[4] * factor3);
    const PylithScalar dStress13dStrain13 = 1.0/
      (a + 2.0 * devStressTau[5] * devStressTpdt[5] * factor3);
    
    /// Compute tangent matrix.
    // Form elastic constants.
    elasticConsts[ 0] = bulkModulus + 2.0 * dStress11dStrain11/3.0; // C1111
    elasticConsts[ 1] = bulkModulus -       dStress11dStrain11/3.0; // C1122
    elasticConsts[ 2] = elasticConsts[ 1]; // C1133
    elasticConsts[ 3] = 0.0; // C1112
    elasticConsts[ 4] = 0.0; // C1123
    elasticConsts[ 5] = 0.0; // C1113
    elasticConsts[ 6] = elasticConsts[ 1]; // C2211
    elasticConsts[ 7] = bulkModulus + 2.0 * dStress22dStrain22/3.0; // C2222
    elasticConsts[ 8] = bulkModulus -       dStress22dStrain22/3.0; // C2233
    elasticConsts[ 9] = 0.0; // C2212
    elasticConsts[10] = 0.0; // C2223
    elasticConsts[11] = 0.0; // C2213
    elasticConsts[12] = elasticConsts[ 1]; // C3311
    elasticConsts[13] = elasticConsts[ 8]; // C3322
    elasticConsts[14] = bulkModulus + 2.0 * dStress33dStrain33/3.0; // C3333
    elasticConsts[15] = 0; // C3312
    elasticConsts[16] = 0; // C3323
    elasticConsts[17] = 0; // C3313
    elasticConsts[18] = 0; // C1211
    elasticConsts[19] = 0; // C1222
    elasticConsts[20] = 0; // C1233
    elasticConsts[21] = dStress12dStrain12; // C1212
    elasticConsts[22] = 0; // C1223
    elasticConsts[23] = 0; // C1213
    elasticConsts[24] = 0; // C2311
    elasticConsts[25] = 0; // C2322
    elasticConsts[26] = 0; // C2333
    elasticConsts[27] = 0; // C2312
    elasticConsts[28] = dStress23dStrain23; // C2323
    elasticConsts[29] = 0; // C2313
    elasticConsts[30] = 0; // C1311
    elasticConsts[31] = 0; // C1322
    elasticConsts[32] = 0; // C1333
    elasticConsts[33] = 0; // C1312
    elasticConsts[34] = 0; // C1323
    elasticConsts[35] = dStress13dStrain13; // C1313
    PetscLogFlops(114);
  } // else
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::PowerLaw3D::_updateStateVarsElastic(
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
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLaw3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLaw3D::tensorSize == initialStrainSize);

  const int tensorSize = 6;
  assert(_tensorSize == tensorSize);

  const bool computeStateVars = true;
  PylithScalar stress[tensorSize] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const int stressSize = strainSize;
  _calcStressElastic(stress, stressSize,
		     properties, numProperties,
		     stateVars, numStateVars,
		     totalStrain, strainSize,
		     initialStress, initialStressSize,
		     initialStrain, initialStrainSize,
		     computeStateVars);

  for (int iComp=0; iComp < tensorSize; ++iComp) {
    stateVars[s_viscousStrain+iComp] = 0.0;
    stateVars[s_stress+iComp] = stress[iComp];
  } // for

  _needNewJacobian = true;
} // _updateStateVarsElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::PowerLaw3D::_updateStateVarsViscoelastic(
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
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLaw3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLaw3D::tensorSize == initialStrainSize);

  const int tensorSize = 6;
  assert(_tensorSize == tensorSize);

  // For now, we are duplicating the functionality of _calcStressViscoelastic,
  // since otherwise we would have to redo a lot of calculations.
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar referenceStrainRate = properties[p_referenceStrainRate];
  const PylithScalar referenceStress = properties[p_referenceStress];
  const PylithScalar powerLawExp = properties[p_powerLawExponent];

  const PylithScalar visStrainT[tensorSize] = {stateVars[s_viscousStrain],
					       stateVars[s_viscousStrain + 1],
					       stateVars[s_viscousStrain + 2],
					       stateVars[s_viscousStrain + 3],
					       stateVars[s_viscousStrain + 4],
					       stateVars[s_viscousStrain + 5]};

  const PylithScalar stressT[tensorSize] = {stateVars[s_stress],
					    stateVars[s_stress + 1],
					    stateVars[s_stress + 2],
					    stateVars[s_stress + 3],
					    stateVars[s_stress + 4],
					    stateVars[s_stress + 5]};
  
  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2/3.0;
  const PylithScalar ae = 1.0/mu2;
  const PylithScalar diag[tensorSize] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Need to figure out how time integration parameter alpha is going to be
  // specified.  It should probably be specified in the problem definition and
  // then used only by the material types that use it.  For now we are setting
  // it to 0.5, which should probably be the default value.
  const PylithScalar alpha = 0.5;
  const PylithScalar timeFac = _dt * (1.0 - alpha);

  // Initial stress values
  const PylithScalar meanStressInitial = (initialStress[0] + initialStress[1] +
					  initialStress[2])/3.0;
  const PylithScalar devStressInitial[tensorSize] = {
    initialStress[0] - meanStressInitial,
    initialStress[1] - meanStressInitial,
    initialStress[2] - meanStressInitial,
    initialStress[3],
    initialStress[4],
    initialStress[5]
  };
  const PylithScalar stressInvar2Initial =
    0.5 * scalarProduct3D(devStressInitial, devStressInitial);

  // Initial strain values
  const PylithScalar meanStrainInitial = (initialStrain[0] + initialStrain[1] +
					  initialStrain[2])/3.0;
  
  // Values for current time step
  const PylithScalar e11 = totalStrain[0];
  const PylithScalar e22 = totalStrain[1];
  const PylithScalar e33 = totalStrain[2];
  const PylithScalar meanStrainTpdt = (e11 + e22 + e33)/3.0 - meanStrainInitial;
  const PylithScalar meanStressTpdt = 3.0 * bulkModulus * meanStrainTpdt;

  // Note that I use the initial strain rather than the deviatoric initial
  // strain since otherwise the initial mean strain would get used twice.
  const PylithScalar strainPPTpdt[tensorSize] = {
    totalStrain[0] - meanStrainTpdt - visStrainT[0] - initialStrain[0],
    totalStrain[1] - meanStrainTpdt - visStrainT[1] - initialStrain[1],
    totalStrain[2] - meanStrainTpdt - visStrainT[2] - initialStrain[2],
    totalStrain[3] - visStrainT[3] - initialStrain[3],
    totalStrain[4] - visStrainT[4] - initialStrain[4],
    totalStrain[5] - visStrainT[5] - initialStrain[5]
  };
  const PylithScalar strainPPInvar2Tpdt =
    0.5 * scalarProduct3D(strainPPTpdt, strainPPTpdt);

  // Values for previous time step
  const PylithScalar meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
  const PylithScalar devStressT[tensorSize] = { stressT[0] - meanStressT,
						stressT[1] - meanStressT,
						stressT[2] - meanStressT,
						stressT[3],
						stressT[4],
						stressT[5] };
  const PylithScalar stressInvar2T =
    0.5 * scalarProduct3D(devStressT, devStressT);
  const PylithScalar effStressT = sqrt(stressInvar2T);

  // Finish defining parameters needed for root-finding algorithm.
  const PylithScalar b = strainPPInvar2Tpdt +
    ae * scalarProduct3D(strainPPTpdt, devStressInitial) +
    ae * ae * stressInvar2Initial;
  const PylithScalar c =
    (scalarProduct3D(strainPPTpdt, devStressT) +
     ae * scalarProduct3D(devStressT, devStressInitial)) * timeFac;
  const PylithScalar d = timeFac * effStressT;
  PetscLogFlops(92);

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
      EffectiveStress::calculate<PowerLaw3D>(effStressInitialGuess,
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

  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStressTpdt = factor1 *
      (strainPPTpdt[iComp] - factor2 * devStressT[iComp] +
       ae * devStressInitial[iComp]);
    stateVars[s_stress+iComp] = devStressTpdt + diag[iComp] *
      (meanStressTpdt + meanStressInitial);
    devStressTau = (1.0 - alpha) * devStressT[iComp] + alpha * devStressTpdt;
    stateVars[s_viscousStrain+iComp] += _dt * gammaTau * devStressTau;
  } // for

  _needNewJacobian = true;
  PetscLogFlops(14 + _tensorSize * 15);

} // _updateStateVarsViscoelastic

// End of file 
