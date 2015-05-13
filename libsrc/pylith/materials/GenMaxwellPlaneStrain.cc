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

#include "GenMaxwellPlaneStrain.hh" // implementation of object methods

#include "ViscoelasticMaxwell.hh" // USES computeVisStrain
#include "Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES MAXSCALAR

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _GenMaxwellPlaneStrain{

      /// Number of Maxwell models in parallel.
      const int numMaxwellModels = 3;

      // Dimension of material.
      const int dimension = 2;

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 3;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 9;

      /// Number of physical properties.
      const int numProperties = 5;
      
      /// Physical properties.
      const Metadata::ParamDescription properties[numProperties] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
	{ "shear_ratio", numMaxwellModels, pylith::topology::FieldBase::OTHER },
	{ "maxwell_time", numMaxwellModels, pylith::topology::FieldBase::OTHER },
      };
      // Values expected in properties spatial database.  :KLUDGE: Not
      // generalized over number of models.
      const int numDBProperties = 3 + 2 * numMaxwellModels;
      const char* dbProperties[numDBProperties] = {
	"density", "vs", "vp",
	"shear-ratio-1",
	"shear-ratio-2",
	"shear-ratio-3",
	"viscosity-1",
	"viscosity-2",
	"viscosity-3",
      };
      
      /// Number of state variables.
      const int numStateVars = 2 + numMaxwellModels;
      
      /// State variables. :KLUDGE: Not generalized over number of models.
      const Metadata::ParamDescription stateVars[numStateVars] = {
	{ "stress_zz_initial", 1, pylith::topology::FieldBase::SCALAR },
	{ "total_strain", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "viscous_strain_1", 4, pylith::topology::FieldBase::OTHER },
	{ "viscous_strain_2", 4, pylith::topology::FieldBase::OTHER },
	{ "viscous_strain_3", 4, pylith::topology::FieldBase::OTHER },
      };

      // Values expected in state variables spatial database
      const int numDBStateVars = 1 + tensorSize + numMaxwellModels * 4;
      const char* dbStateVars[numDBStateVars] = {"stress-zz-initial",
						 "total-strain-xx",
						 "total-strain-yy",
						 "total-strain-xy",
						 "viscous-strain-1-xx",
						 "viscous-strain-1-yy",
						 "viscous-strain-1-zz",
						 "viscous-strain-1-xy",
						 "viscous-strain-2-xx",
						 "viscous-strain-2-yy",
						 "viscous-strain-2-zz",
						 "viscous-strain-2-xy",
						 "viscous-strain-3-xx",
						 "viscous-strain-3-yy",
						 "viscous-strain-3-zz",
						 "viscous-strain-3-xy",
      };

    } // _GenMaxwellPlaneStrain
  } // materials
} // pylith

// Indices of physical properties
const int pylith::materials::GenMaxwellPlaneStrain::p_density = 0;

const int pylith::materials::GenMaxwellPlaneStrain::p_muEff =
  pylith::materials::GenMaxwellPlaneStrain::p_density + 1;

const int pylith::materials::GenMaxwellPlaneStrain::p_lambdaEff =
  pylith::materials::GenMaxwellPlaneStrain::p_muEff + 1;

const int pylith::materials::GenMaxwellPlaneStrain::p_shearRatio =
  pylith::materials::GenMaxwellPlaneStrain::p_lambdaEff + 1;

const int pylith::materials::GenMaxwellPlaneStrain::p_maxwellTime =
  pylith::materials::GenMaxwellPlaneStrain::p_shearRatio +
  pylith::materials::_GenMaxwellPlaneStrain::numMaxwellModels;

// Indices of database values (order must match dbProperties)
const int pylith::materials::GenMaxwellPlaneStrain::db_density = 0;

const int pylith::materials::GenMaxwellPlaneStrain::db_vs =
  pylith::materials::GenMaxwellPlaneStrain::db_density + 1;

const int pylith::materials::GenMaxwellPlaneStrain::db_vp =
  pylith::materials::GenMaxwellPlaneStrain::db_vs + 1;

const int pylith::materials::GenMaxwellPlaneStrain::db_shearRatio =
  pylith::materials::GenMaxwellPlaneStrain::db_vp + 1;

const int pylith::materials::GenMaxwellPlaneStrain::db_viscosity =
  pylith::materials::GenMaxwellPlaneStrain::db_shearRatio + 
  pylith::materials::_GenMaxwellPlaneStrain::numMaxwellModels;

// Indices of state variables
const int pylith::materials::GenMaxwellPlaneStrain::s_stressZZInitial = 0;

const int pylith::materials::GenMaxwellPlaneStrain::s_totalStrain =
  pylith::materials::GenMaxwellPlaneStrain::s_stressZZInitial + 1;

const int pylith::materials::GenMaxwellPlaneStrain::s_viscousStrain1 = 
  pylith::materials::GenMaxwellPlaneStrain::s_totalStrain +
  pylith::materials::_GenMaxwellPlaneStrain::tensorSize;

const int pylith::materials::GenMaxwellPlaneStrain::s_viscousStrain2 = 
  pylith::materials::GenMaxwellPlaneStrain::s_viscousStrain1 + 4;

const int pylith::materials::GenMaxwellPlaneStrain::s_viscousStrain3 = 
  pylith::materials::GenMaxwellPlaneStrain::s_viscousStrain2 + 4;

// Indices of database values (order must match dbStateVars)
const int pylith::materials::GenMaxwellPlaneStrain::db_stressZZInitial = 0;

const int pylith::materials::GenMaxwellPlaneStrain::db_totalStrain =
  pylith::materials::GenMaxwellPlaneStrain::db_stressZZInitial + 1;

const int pylith::materials::GenMaxwellPlaneStrain::db_viscousStrain1 =
  pylith::materials::GenMaxwellPlaneStrain::db_totalStrain +
  pylith::materials::_GenMaxwellPlaneStrain::tensorSize;

const int pylith::materials::GenMaxwellPlaneStrain::db_viscousStrain2 =
  pylith::materials::GenMaxwellPlaneStrain::db_viscousStrain1 + 4;

const int pylith::materials::GenMaxwellPlaneStrain::db_viscousStrain3 =
  pylith::materials::GenMaxwellPlaneStrain::db_viscousStrain2 + 4;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::GenMaxwellPlaneStrain::GenMaxwellPlaneStrain(void) :
  ElasticMaterial(_GenMaxwellPlaneStrain::dimension,
		  _GenMaxwellPlaneStrain::tensorSize,
		  _GenMaxwellPlaneStrain::numElasticConsts,
		  Metadata(_GenMaxwellPlaneStrain::properties,
			   _GenMaxwellPlaneStrain::numProperties,
			   _GenMaxwellPlaneStrain::dbProperties,
			   _GenMaxwellPlaneStrain::numDBProperties,
			   _GenMaxwellPlaneStrain::stateVars,
			   _GenMaxwellPlaneStrain::numStateVars,
			   _GenMaxwellPlaneStrain::dbStateVars,
			   _GenMaxwellPlaneStrain::numDBStateVars)),
  _calcElasticConstsFn(0),
  _calcStressFn(0),
  _updateStateVarsFn(0)  
{ // constructor
  useElasticBehavior(false);
  _viscousStrain.resize(_GenMaxwellPlaneStrain::numMaxwellModels * 4);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::GenMaxwellPlaneStrain::~GenMaxwellPlaneStrain(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set whether elastic or inelastic constitutive relations are used.
void
pylith::materials::GenMaxwellPlaneStrain::useElasticBehavior(const bool flag)
{ // useElasticBehavior
  if (flag) {
    _calcStressFn = 
      &pylith::materials::GenMaxwellPlaneStrain::_calcStressElastic;
    _calcElasticConstsFn = 
      &pylith::materials::GenMaxwellPlaneStrain::_calcElasticConstsElastic;
    _updateStateVarsFn = 
      &pylith::materials::GenMaxwellPlaneStrain::_updateStateVarsElastic;

  } else {
    _calcStressFn = 
      &pylith::materials::GenMaxwellPlaneStrain::_calcStressViscoelastic;
    _calcElasticConstsFn = 
      &pylith::materials::GenMaxwellPlaneStrain::_calcElasticConstsViscoelastic;
    _updateStateVarsFn = 
      &pylith::materials::GenMaxwellPlaneStrain::_updateStateVarsViscoelastic;
  } // if/else
} // useElasticBehavior

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::GenMaxwellPlaneStrain::_dbToProperties(
					    PylithScalar* const propValues,
					    const scalar_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_GenMaxwellPlaneStrain::numDBProperties == numDBValues);

  const PylithScalar density = dbValues[db_density];
  const PylithScalar vs = dbValues[db_vs];
  const PylithScalar vp = dbValues[db_vp];
  const int numMaxwellModels = _GenMaxwellPlaneStrain::numMaxwellModels;
 
  if (density <= 0.0 || vs <= 0.0 || vp <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for physical "
	<< "properties.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
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
  propValues[p_muEff] = mu;
  propValues[p_lambdaEff] = lambda;

  PylithScalar visFrac = 0.0;
  for (int imodel = 0; imodel < numMaxwellModels; ++imodel) 
    visFrac += dbValues[db_shearRatio + imodel];
  if (visFrac > 1.0) {
    std::ostringstream msg;
    msg << "Shear modulus ratios sum to a value greater than 1.0 for\n"
	<< "Generalized Maxwell model.\n"
	<< "Ratio 1: " << propValues[db_shearRatio  ] << "\n"
	<< "Ratio 2: " << propValues[db_shearRatio+1] << "\n"
	<< "Ratio 3: " << propValues[db_shearRatio+2] << "\n"
	<< "Total:   " << visFrac << "\n";
    throw std::runtime_error(msg.str());
  } // if

  // Loop over number of Maxwell models.
  for (int imodel = 0; imodel < numMaxwellModels; ++imodel) {
    PylithScalar muRatio = dbValues[db_shearRatio + imodel];
    PylithScalar viscosity = dbValues[db_viscosity + imodel];
    PylithScalar muFac = muRatio * mu;
    PylithScalar maxwellTime = pylith::PYLITH_MAXSCALAR;
    if (muFac > 0.0)
      maxwellTime = viscosity / muFac;
    if (muRatio < 0.0 || viscosity < 0.0 || muFac < 0.0 || maxwellTime < 0.0) {
      std::ostringstream msg;
      msg << "Found negative value(s) for physical properties.\n"
	  << "muRatio: " << muRatio << "\n"
	  << "viscosity: " << viscosity << "\n"
	  << "muFac: " << muFac << "\n"
	  << "maxwellTime: " << maxwellTime << "\n";
      throw std::runtime_error(msg.str());
    } // if
    propValues[p_shearRatio + imodel] = muRatio;
    propValues[p_maxwellTime + imodel] = maxwellTime;
  } // for

  PetscLogFlops(6 + 3 * numMaxwellModels);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::GenMaxwellPlaneStrain::_nondimProperties(PylithScalar* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numPropsQuadPt);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  values[p_density] = 
    _normalizer->nondimensionalize(values[p_density], densityScale);
  values[p_muEff] = 
    _normalizer->nondimensionalize(values[p_muEff], pressureScale);
  values[p_lambdaEff] = 
    _normalizer->nondimensionalize(values[p_lambdaEff], pressureScale);
  const int numMaxwellModels = _GenMaxwellPlaneStrain::numMaxwellModels;
  _normalizer->nondimensionalize(&values[p_maxwellTime],
				 numMaxwellModels, timeScale);
  
  PetscLogFlops(3 + 1 * numMaxwellModels);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::GenMaxwellPlaneStrain::_dimProperties(PylithScalar* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numPropsQuadPt);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  values[p_density] = 
    _normalizer->dimensionalize(values[p_density], densityScale);
  values[p_muEff] = 
    _normalizer->dimensionalize(values[p_muEff], pressureScale);
  values[p_lambdaEff] = 
    _normalizer->dimensionalize(values[p_lambdaEff], pressureScale);
  const int numMaxwellModels = _GenMaxwellPlaneStrain::numMaxwellModels;
  _normalizer->dimensionalize(&values[p_maxwellTime],
			      numMaxwellModels, timeScale);

  PetscLogFlops(3 + 1 * numMaxwellModels);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute initial state variables from values in spatial database.
void
pylith::materials::GenMaxwellPlaneStrain::_dbToStateVars(
					PylithScalar* const stateValues,
					const scalar_array& dbValues)
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_GenMaxwellPlaneStrain::numDBStateVars == numDBValues);

  const int numMaxwellModels = _GenMaxwellPlaneStrain::numMaxwellModels;
  const int totalSize = 1 + _tensorSize + 4 * numMaxwellModels;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  for (int i=0; i < totalSize; ++i)
    stateValues[i] = dbValues[i];

  PetscLogFlops(0);
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::materials::GenMaxwellPlaneStrain::_nondimStateVars(
		                            PylithScalar* const values,
		                            const int nvalues) const
{ // _nondimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numVarsQuadPt);

  const PylithScalar pressureScale = _normalizer->pressureScale();
  _normalizer->nondimensionalize(&values[s_stressZZInitial], 1, pressureScale);

  PetscLogFlops(1);
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::materials::GenMaxwellPlaneStrain::_dimStateVars(
						PylithScalar* const values,
						const int nvalues) const
{ // _dimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numVarsQuadPt);

  const PylithScalar pressureScale = _normalizer->pressureScale();
  _normalizer->dimensionalize(&values[s_stressZZInitial], 1, pressureScale);

  PetscLogFlops(1);
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::GenMaxwellPlaneStrain::_calcDensity(PylithScalar* const density,
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
// Compute stress tensor at location from properties as an elastic
// material.
void
pylith::materials::GenMaxwellPlaneStrain::_calcStressElastic(
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
  assert(_GenMaxwellPlaneStrain::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStrainSize);

  const PylithScalar mu = properties[p_muEff];
  const PylithScalar lambda = properties[p_lambdaEff];
  const PylithScalar mu2 = 2.0 * mu;

  // :TODO: Need to consider initial state variables????
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
pylith::materials::GenMaxwellPlaneStrain::_calcStressViscoelastic(
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
  assert(_GenMaxwellPlaneStrain::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStrainSize);

  const int numMaxwellModels = _GenMaxwellPlaneStrain::numMaxwellModels;
  const int tensorSize = _GenMaxwellPlaneStrain::tensorSize;

  const PylithScalar mu = properties[p_muEff];
  const PylithScalar lambda = properties[p_lambdaEff];
  const PylithScalar muRatio[numMaxwellModels] = {
    properties[p_shearRatio  ],
    properties[p_shearRatio+1],
    properties[p_shearRatio+2]
  };
  const PylithScalar stressZZInitial = stateVars[s_stressZZInitial];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2/3.0;

  // Initial stress and strain values
  const PylithScalar meanStrainInitial = (initialStrain[0] + initialStrain[1]) / 3.0;
  const PylithScalar meanStressInitial = (initialStress[0] + initialStress[1] +
				    stressZZInitial) / 3.0;
  const PylithScalar devStrainInitial[] = {initialStrain[0] - meanStrainInitial,
				     initialStrain[1] - meanStrainInitial,
				     initialStrain[2]};
  const PylithScalar devStressInitial[] = {initialStress[0] - meanStressInitial,
				     initialStress[1] - meanStressInitial,
				     initialStress[2]};

  // Mean stress and strain for t + dt
  const PylithScalar meanStrainTpdt = (totalStrain[0] + totalStrain[1]) / 3.0;
  const PylithScalar meanStressTpdt = 3.0 * bulkModulus *
    (meanStrainTpdt - meanStrainInitial) + meanStressInitial;

  const PylithScalar diag[] = { 1.0, 1.0, 0.0 };
  
  PylithScalar visFrac = 0.0;
  for (int imodel=0; imodel < numMaxwellModels; ++imodel) 
    visFrac += muRatio[imodel];
  assert(visFrac <= 1.0);
  const PylithScalar elasFrac = 1.0 - visFrac;

  PetscLogFlops(18 + numMaxwellModels);

  // Get viscous strains
  if (computeStateVars) {
    _computeStateVars(stateVars, numStateVars,
		      properties, numProperties,
		      totalStrain, strainSize,
		      initialStress, initialStressSize,
		      initialStrain, initialStrainSize);
  } else {
    int index = 0;
    for (int iComp=0; iComp < 4; ++iComp)
      _viscousStrain[index++] = stateVars[s_viscousStrain1+iComp];
    for (int iComp=0; iComp < 4; ++iComp)
      _viscousStrain[index++] = stateVars[s_viscousStrain2+iComp];
    for (int iComp=0; iComp < 4; ++iComp)
      _viscousStrain[index++] = stateVars[s_viscousStrain3+iComp];
  } // else

  // Compute new stresses
  PylithScalar devStrainTpdt = 0.0;
  PylithScalar devStressTpdt = 0.0;
  const int visIndex[] = {0, 1, 3};
  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp] * meanStrainTpdt -
      devStrainInitial[iComp];
    devStressTpdt = elasFrac * devStrainTpdt;
    for (int model=0; model < numMaxwellModels; ++model) {
      devStressTpdt += muRatio[model] *
	_viscousStrain[4 * model + visIndex[iComp]];
    } // for

    devStressTpdt = mu2 * devStressTpdt + devStressInitial[iComp];
    stress[iComp] = diag[iComp] * meanStressTpdt + devStressTpdt;
  } // for

  PetscLogFlops((9 + 2 * numMaxwellModels) * tensorSize);
} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::GenMaxwellPlaneStrain::_calcElasticConstsElastic(
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
  assert(_GenMaxwellPlaneStrain::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStrainSize);
 
  const PylithScalar mu = properties[p_muEff];
  const PylithScalar lambda = properties[p_lambdaEff];

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
pylith::materials::GenMaxwellPlaneStrain::_calcElasticConstsViscoelastic(
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
  assert(_GenMaxwellPlaneStrain::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStrainSize);

  const int numMaxwellModels = _GenMaxwellPlaneStrain::numMaxwellModels;

  const PylithScalar mu = properties[p_muEff];
  const PylithScalar lambda = properties[p_lambdaEff];
  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2 / 3.0;

  // Compute viscous contribution.
  PylithScalar visFac = 0.0;
  PylithScalar visFrac = 0.0;
  PylithScalar shearRatio = 0.0;
  for (int imodel = 0; imodel < numMaxwellModels; ++imodel) {
    shearRatio = properties[p_shearRatio + imodel];
    PylithScalar maxwellTime = pylith::PYLITH_MAXSCALAR;
    visFrac += shearRatio;
    if (shearRatio != 0.0) {
      maxwellTime = properties[p_maxwellTime + imodel];
      visFac +=
	shearRatio * ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTime);
    } // if
  } // for
  PylithScalar elasFrac = 1.0 - visFrac;
  PylithScalar shearFac = elasFrac + visFac;

  elasticConsts[ 0] = bulkModulus + 4.0 * mu / 3.0 * shearFac; // C1111
  elasticConsts[ 1] = bulkModulus - 2.0 * mu / 3.0 * shearFac; // C1122
  elasticConsts[ 2] = 0; // C1112
  elasticConsts[ 3] = elasticConsts[1]; // C2211
  elasticConsts[ 4] = elasticConsts[0]; // C2222
  elasticConsts[ 5] = 0; // C2212
  elasticConsts[ 6] = 0; // C1211
  elasticConsts[ 7] = 0; // C1222
  elasticConsts[ 8] = 2.0 * mu * shearFac; // C1212

  PetscLogFlops(15 + 3 * numMaxwellModels);
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::GenMaxwellPlaneStrain::_updateStateVarsElastic(
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
  assert(_GenMaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;

  const PylithScalar strainTpdt[] = {totalStrain[0] - initialStrain[0],
			       totalStrain[1] - initialStrain[1],
			       0.0,
			       totalStrain[2] - initialStrain[2]};
  const PylithScalar meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;

  const PylithScalar diag[] = { 1.0, 1.0, 1.0, 0.0};

  // Update total strain
  for (int iComp=0; iComp < tensorSize; ++iComp)
    stateVars[s_totalStrain+iComp] = totalStrain[iComp];

  // Initialize all viscous strains to deviatoric elastic strains.
  PylithScalar devStrain = 0.0;
  for (int iComp=0; iComp < 4; ++iComp) {
    devStrain = strainTpdt[iComp] - diag[iComp] * meanStrainTpdt;
    // Maxwell model 1
    stateVars[s_viscousStrain1+iComp] = devStrain;
    // Maxwell model 2
    stateVars[s_viscousStrain2+iComp] = devStrain;
    // Maxwell model 3
    stateVars[s_viscousStrain3+iComp] = devStrain;
  } // for
  PetscLogFlops(5 + 2 * 4);

  _needNewJacobian = true;
} // _updateStateVarsElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::GenMaxwellPlaneStrain::_updateStateVarsViscoelastic(
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
  assert(_GenMaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStrainSize);

  _computeStateVars(stateVars, numStateVars,
		    properties, numProperties,
		    totalStrain, strainSize,
		    initialStress, initialStressSize,
		    initialStrain, initialStrainSize);

  const int tensorSize = _tensorSize;

  // Total strain
  for (int iComp=0; iComp < tensorSize; ++iComp)
    stateVars[s_totalStrain+iComp] = totalStrain[iComp];

  // Viscous strain
  int index = 0;
  for (int iComp=0; iComp < 4; ++iComp)
    stateVars[s_viscousStrain1+iComp] = _viscousStrain[index++];
  for (int iComp=0; iComp < 4; ++iComp)
    stateVars[s_viscousStrain2+iComp] = _viscousStrain[index++];
  for (int iComp=0; iComp < 4; ++iComp)
    stateVars[s_viscousStrain3+iComp] = _viscousStrain[index++];

  _needNewJacobian = false;
} // _updateStateVarsViscoelastic

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::GenMaxwellPlaneStrain::_stableTimeStepImplicit(
					   const PylithScalar* properties,
					   const int numProperties,
					   const PylithScalar* stateVars,
					   const int numStateVars) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);

  PylithScalar dtStable = pylith::PYLITH_MAXSCALAR;

  const int numMaxwellModels = _GenMaxwellPlaneStrain::numMaxwellModels;
  for (int i=0; i < numMaxwellModels; ++i) {
    const PylithScalar maxwellTime = properties[p_maxwellTime+i];
    const PylithScalar dt = 0.2*maxwellTime;
    if (dt < dtStable)
      dtStable = dt;
  } // for

  return dtStable;
} // _stableTimeStepImplicit


// ----------------------------------------------------------------------
// Get stable time step for explicit time integration.
PylithScalar
pylith::materials::GenMaxwellPlaneStrain::_stableTimeStepExplicit(const PylithScalar* properties,
								  const int numProperties,
								  const PylithScalar* stateVars,
								  const int numStateVars,
								  const double minCellWidth) const
{ // _stableTimeStepExplicit
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
 
  const PylithScalar mu = properties[p_muEff];
  const PylithScalar lambda = properties[p_lambdaEff];
  const PylithScalar density = properties[p_density];

  assert(density > 0.0);
  const PylithScalar vp = sqrt((lambda + 2*mu) / density);

  const PylithScalar dtStable = minCellWidth / vp;
  return dtStable;
} // _stableTimeStepExplicit


// ----------------------------------------------------------------------
// Compute viscous strain for current time step.
void
pylith::materials::GenMaxwellPlaneStrain::_computeStateVars(
					       const PylithScalar* stateVars,
					       const int numStateVars,
					       const PylithScalar* properties,
					       const int numProperties,
					       const PylithScalar* totalStrain,
					       const int strainSize,
					       const PylithScalar* initialStress,
					       const int initialStressSize,
					       const PylithScalar* initialStrain,
					       const int initialStrainSize)
{ // _computeStateVars
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellPlaneStrain::tensorSize == initialStrainSize);

  const int numMaxwellModels = _GenMaxwellPlaneStrain::numMaxwellModels;

  const PylithScalar muRatio[numMaxwellModels] = {
    properties[p_shearRatio  ],
    properties[p_shearRatio+1],
    properties[p_shearRatio+2]
  };
  const PylithScalar maxwellTime[numMaxwellModels] = {
    properties[p_maxwellTime  ],
    properties[p_maxwellTime+1],
    properties[p_maxwellTime+2]
  };

  const PylithScalar strainTpdt[] = {totalStrain[0],
			       totalStrain[1],
			       0.0,
			       totalStrain[2]};
  const PylithScalar strainT[] = {stateVars[s_totalStrain+0],
			    stateVars[s_totalStrain+1],
			    0.0,
			    stateVars[s_totalStrain+2]};
  
  const PylithScalar meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;
  const PylithScalar meanStrainT = (strainT[0] + strainT[1]) / 3.0;

  const PylithScalar diag[] = { 1.0, 1.0, 1.0, 0.0 };

  PetscLogFlops(4);

  // Compute Prony series terms
  scalar_array dq(numMaxwellModels);
  dq = 0.0;
  for (int i=0; i < numMaxwellModels; ++i)
    if (muRatio[i] != 0.0)
      dq[i] = ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTime[i]);

  // Compute new viscous strains
  PylithScalar devStrainTpdt = 0.0;
  PylithScalar devStrainT = 0.0;
  PylithScalar deltaStrain = 0.0;
  
  for (int iComp=0; iComp < 4; ++iComp) {
    devStrainTpdt = strainTpdt[iComp] - diag[iComp] * meanStrainTpdt;
    devStrainT = strainT[iComp] - diag[iComp] * meanStrainT;
    deltaStrain = devStrainTpdt - devStrainT;

    // Maxwell model 1
    int imodel = 0;
    if (0.0 != muRatio[imodel]) {
      _viscousStrain[imodel * 4 + iComp] = exp(-_dt / maxwellTime[imodel]) *
	stateVars[s_viscousStrain1 + iComp] + dq[imodel] * deltaStrain;
      PetscLogFlops(6);
    } // if

    // Maxwell model 2
    imodel = 1;
    if (0.0 != muRatio[imodel]) {
      _viscousStrain[imodel * 4 + iComp] = exp(-_dt / maxwellTime[imodel]) *
	stateVars[s_viscousStrain2 + iComp] + dq[imodel] * deltaStrain;
      PetscLogFlops(6);
    } // if

    // Maxwell model 3
    imodel = 2;
    if (0.0 != muRatio[imodel]) {
      _viscousStrain[imodel * 4 + iComp] = exp(-_dt / maxwellTime[imodel]) *
	stateVars[s_viscousStrain3 + iComp] + dq[imodel] * deltaStrain;
      PetscLogFlops(6);
    } // if

  } // for

  PetscLogFlops(5 * 4);
} // _computeStateVars


// End of file 
