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

#include "GenMaxwellQpQsIsotropic3D.hh" // implementation of object methods

#include "ViscoelasticMaxwell.hh" // USES computeVisStrain
#include "Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES MAXSCALAR

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

#include <iostream> // TEMPORARY
// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _GenMaxwellQpQsIsotropic3D{

      /// Number of Maxwell models in parallel.
      const int numMaxwellModels = 3;

      // Dimension of material.
      const int dimension = 3;

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 6;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 36;

      /// Number of physical properties.
      const int numProperties = 7;
      
      /// Physical properties.
      const Metadata::ParamDescription properties[numProperties] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "k", 1, pylith::topology::FieldBase::SCALAR },
	{ "shear_ratio", 
	  numMaxwellModels, pylith::topology::FieldBase::OTHER },
	{ "maxwell_time_shear", 
	  numMaxwellModels, pylith::topology::FieldBase::OTHER },
	{ "bulk_ratio", numMaxwellModels, pylith::topology::FieldBase::OTHER },
	{ "maxwell_time_bulk", 
	  numMaxwellModels, pylith::topology::FieldBase::OTHER },
      };
      // Values expected in properties spatial database.
      // :KLUDGE: Not generalized over number of models.
      const int numDBProperties = 3 + 4*numMaxwellModels;
      const char* dbProperties[numDBProperties] = {
	"density", "vs", "vp",
	"shear-ratio-1",
	"shear-ratio-2",
	"shear-ratio-3",
	"shear-viscosity-1",
	"shear-viscosity-2",
	"shear-viscosity-3",
	"bulk-ratio-1",
	"bulk-ratio-2",
	"bulk-ratio-3",
	"bulk-viscosity-1",
	"bulk-viscosity-2",
	"bulk-viscosity-3",
      };
      
      /// Number of state variables.
      const int numStateVars = 3;
      
      /// State variables.
      const Metadata::ParamDescription stateVars[numStateVars] = {
	{ "total_strain", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "viscous_deviatoric_strain", 
	  numMaxwellModels*tensorSize, pylith::topology::FieldBase::OTHER },
	{ "viscous_mean_strain", 
	  numMaxwellModels, pylith::topology::FieldBase::OTHER },
      };

      // Values expected in state variables spatial database
      const int numDBStateVars = 
	tensorSize + numMaxwellModels*tensorSize + numMaxwellModels;
      const char* dbStateVars[numDBStateVars] = {"total-strain-xx",
						 "total-strain-yy",
						 "total-strain-zz",
						 "total-strain-xy",
						 "total-strain-yz",
						 "total-strain-xz",
						 "viscous-deviatoric-strain-1-xx",
						 "viscous-deviatoric-strain-1-yy",
						 "viscous-deviatoric-strain-1-zz",
						 "viscous-deviatoric-strain-1-xy",
						 "viscous-deviatoric-strain-1-yz",
						 "viscous-deviatoric-strain-1-xz",
						 "viscous-deviatoric-strain-2-xx",
						 "viscous-deviatoric-strain-2-yy",
						 "viscous-deviatoric-strain-2-zz",
						 "viscous-deviatoric-strain-2-xy",
						 "viscous-deviatoric-strain-2-yz",
						 "viscous-deviatoric-strain-2-xz",
						 "viscous-deviatoric-strain-3-xx",
						 "viscous-deviatoric-strain-3-yy",
						 "viscous-deviatoric-strain-3-zz",
						 "viscous-deviatoric-strain-3-xy",
						 "viscous-deviatoric-strain-3-yz",
						 "viscous-deviatoric-strain-3-xz",
						 "viscous-mean-strain-1",
						 "viscous-mean-strain-2",
						 "viscous-mean-strain-3",
      };

    } // _GenMaxwellQpQsIsotropic3D
  } // materials
} // pylith

// Indices of physical properties
const int pylith::materials::GenMaxwellQpQsIsotropic3D::p_density = 0;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::p_muEff =
  pylith::materials::GenMaxwellQpQsIsotropic3D::p_density + 1;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::p_kEff =
  pylith::materials::GenMaxwellQpQsIsotropic3D::p_muEff + 1;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::p_shearRatio =
  pylith::materials::GenMaxwellQpQsIsotropic3D::p_kEff + 1;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::p_maxwellTimeShear =
  pylith::materials::GenMaxwellQpQsIsotropic3D::p_shearRatio +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::numMaxwellModels;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::p_bulkRatio =
  pylith::materials::GenMaxwellQpQsIsotropic3D::p_maxwellTimeShear +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::numMaxwellModels;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::p_maxwellTimeBulk =
  pylith::materials::GenMaxwellQpQsIsotropic3D::p_bulkRatio +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::numMaxwellModels;

// Indices of database values (order must match dbProperties)
const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_density = 0;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_vs =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_density + 1;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_vp =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_vs + 1;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_shearRatio =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_vp + 1;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_shearViscosity =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_shearRatio + 
  pylith::materials::_GenMaxwellQpQsIsotropic3D::numMaxwellModels;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_bulkRatio =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_shearViscosity + 
  pylith::materials::_GenMaxwellQpQsIsotropic3D::numMaxwellModels;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_bulkViscosity =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_bulkRatio + 
  pylith::materials::_GenMaxwellQpQsIsotropic3D::numMaxwellModels;

// Indices of state variables
const int pylith::materials::GenMaxwellQpQsIsotropic3D::s_totalStrain = 0;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousDevStrain = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::s_totalStrain +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousMeanStrain = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousDevStrain +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::numMaxwellModels *
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

// Indices of database values (order must match dbStateVars)
const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_totalStrain = 0;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousDevStrain1 =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_totalStrain +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousDevStrain2 =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousDevStrain1 +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousDevStrain3 =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousDevStrain2 +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousMeanStrain1 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousDevStrain3 +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousMeanStrain2 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousMeanStrain1 + 1;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousMeanStrain3 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousMeanStrain2 + 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::GenMaxwellQpQsIsotropic3D::GenMaxwellQpQsIsotropic3D(void) :
  ElasticMaterial(_GenMaxwellQpQsIsotropic3D::dimension,
		  _GenMaxwellQpQsIsotropic3D::tensorSize,
		  _GenMaxwellQpQsIsotropic3D::numElasticConsts,
		  Metadata(_GenMaxwellQpQsIsotropic3D::properties,
			   _GenMaxwellQpQsIsotropic3D::numProperties,
			   _GenMaxwellQpQsIsotropic3D::dbProperties,
			   _GenMaxwellQpQsIsotropic3D::numDBProperties,
			   _GenMaxwellQpQsIsotropic3D::stateVars,
			   _GenMaxwellQpQsIsotropic3D::numStateVars,
			   _GenMaxwellQpQsIsotropic3D::dbStateVars,
			   _GenMaxwellQpQsIsotropic3D::numDBStateVars)),
  _calcElasticConstsFn(0),
  _calcStressFn(0),
  _updateStateVarsFn(0)  
{ // constructor
  useElasticBehavior(false);
  _viscousDevStrain.resize(_GenMaxwellQpQsIsotropic3D::numMaxwellModels*_tensorSize);
  _viscousMeanStrain.resize(_GenMaxwellQpQsIsotropic3D::numMaxwellModels);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::GenMaxwellQpQsIsotropic3D::~GenMaxwellQpQsIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set whether elastic or inelastic constitutive relations are used.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::useElasticBehavior(const bool flag)
{ // useElasticBehavior
  if (flag) {
    _calcStressFn = 
      &pylith::materials::GenMaxwellQpQsIsotropic3D::_calcStressElastic;
    _calcElasticConstsFn = 
      &pylith::materials::GenMaxwellQpQsIsotropic3D::_calcElasticConstsElastic;
    _updateStateVarsFn = 
      &pylith::materials::GenMaxwellQpQsIsotropic3D::_updateStateVarsElastic;

  } else {
    _calcStressFn = 
      &pylith::materials::GenMaxwellQpQsIsotropic3D::_calcStressViscoelastic;
    _calcElasticConstsFn = 
      &pylith::materials::GenMaxwellQpQsIsotropic3D::_calcElasticConstsViscoelastic;
    _updateStateVarsFn = 
      &pylith::materials::GenMaxwellQpQsIsotropic3D::_updateStateVarsViscoelastic;
  } // if/else
} // useElasticBehavior

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_dbToProperties(
					    PylithScalar* const propValues,
					    const scalar_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_GenMaxwellQpQsIsotropic3D::numDBProperties == numDBValues);

  const PylithScalar density = dbValues[db_density];
  const PylithScalar vs = dbValues[db_vs];
  const PylithScalar vp = dbValues[db_vp];
  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;
 
  if (density <= 0.0 || vs <= 0.0 || vp <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for physical "
	<< "properties.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if

  const PylithScalar mu = density * vs*vs;
  const PylithScalar k = density * vp*vp - 4.0/3.0*mu;

  if (k <= 0.0) {
    std::ostringstream msg;
    msg << "Attempted to set bulk modulus to nonpositive value.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if
  assert(mu > 0);

  propValues[p_density] = density;
  propValues[p_muEff] = mu;
  propValues[p_kEff] = k;

  PylithScalar visFrac = 0.0;
  for (int imodel = 0; imodel < numMaxwellModels; ++imodel) 
    visFrac += dbValues[db_shearRatio + imodel];
  if (visFrac > 1.0) {
    std::ostringstream msg;
    msg << "Shear modulus ratios sum to a value greater than 1.0 for\n"
	<< "Generalized Maxwell shear model.\n"
	<< "Ratio 1: " << dbValues[db_shearRatio  ] << "\n"
	<< "Ratio 2: " << dbValues[db_shearRatio+1] << "\n"
	<< "Ratio 3: " << dbValues[db_shearRatio+2] << "\n"
	<< "Total:   " << visFrac << "\n";
    throw std::runtime_error(msg.str());
  } // if

  PylithScalar meanFrac = 0.0;
  for (int imodel = 0; imodel < numMaxwellModels; ++imodel) 
    meanFrac += dbValues[db_bulkRatio + imodel];
  if (meanFrac > 1.0) {
    std::ostringstream msg;
    msg << "Bulk modulus ratios sum to a value greater than 1.0 for\n"
	<< "Generalized Maxwell bulk model.\n"
	<< "Ratio 1: " << dbValues[db_bulkRatio  ] << "\n"
	<< "Ratio 2: " << dbValues[db_bulkRatio+1] << "\n"
	<< "Ratio 3: " << dbValues[db_bulkRatio+2] << "\n"
	<< "Total:   " << meanFrac << "\n";
    throw std::runtime_error(msg.str());
  } // if

  // Loop over number of Maxwell models.
  for (int imodel=0; imodel < numMaxwellModels; ++imodel) {
    PylithScalar shearRatio = dbValues[db_shearRatio + imodel];
    PylithScalar bulkRatio = dbValues[db_bulkRatio + imodel];
    PylithScalar shearViscosity = dbValues[db_shearViscosity + imodel];
    PylithScalar bulkViscosity = dbValues[db_bulkViscosity + imodel];
    PylithScalar maxwellTimeShear = pylith::PYLITH_MAXSCALAR;
    PylithScalar maxwellTimeBulk = pylith::PYLITH_MAXSCALAR;
    maxwellTimeShear = shearViscosity / mu;
    maxwellTimeBulk = bulkViscosity / k;
    if (shearRatio < 0.0 || shearViscosity < 0.0 || maxwellTimeShear < 0.0 || 
	bulkRatio < 0.0 || bulkViscosity < 0.0 || maxwellTimeBulk < 0.0) {
      std::ostringstream msg;
      msg << "Found negative value(s) for physical properties.\n"
	  << "shearRatio: " << shearRatio << "\n"
	  << "dev viscosity: " << shearViscosity << "\n"
	  << "maxwellTimeShear: " << maxwellTimeShear << "\n"
	  << "bulkRatio: " << bulkRatio << "\n"
	  << "bulk viscosity: " << bulkViscosity << "\n"
	  << "maxwellTimeBulk: " << maxwellTimeBulk << "\n";
      throw std::runtime_error(msg.str());
    } // if
    propValues[p_shearRatio + imodel] = shearRatio;
    propValues[p_bulkRatio + imodel] = bulkRatio;
    propValues[p_maxwellTimeShear + imodel] = maxwellTimeShear;
    propValues[p_maxwellTimeBulk + imodel] = maxwellTimeBulk;
  } // for

  PetscLogFlops(6+2*numMaxwellModels);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_nondimProperties(PylithScalar* const values,
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
  values[p_kEff] = 
    _normalizer->nondimensionalize(values[p_kEff], pressureScale);

  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;
  _normalizer->nondimensionalize(&values[p_maxwellTimeShear],
				 numMaxwellModels, timeScale);
  _normalizer->nondimensionalize(&values[p_maxwellTimeBulk],
				 numMaxwellModels, timeScale);
  
  PetscLogFlops(3+2*numMaxwellModels);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_dimProperties(PylithScalar* const values,
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
  values[p_kEff] = 
    _normalizer->dimensionalize(values[p_kEff], pressureScale);

  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;
  _normalizer->dimensionalize(&values[p_maxwellTimeShear],
			      numMaxwellModels, timeScale);
  _normalizer->dimensionalize(&values[p_maxwellTimeBulk],
			      numMaxwellModels, timeScale);

  PetscLogFlops(3+2*numMaxwellModels);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute initial state variables from values in spatial database.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_dbToStateVars(
					PylithScalar* const stateValues,
					const scalar_array& dbValues)
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_GenMaxwellQpQsIsotropic3D::numDBStateVars == numDBValues);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == _tensorSize);

  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;
  const int tensorSize = _GenMaxwellQpQsIsotropic3D::tensorSize;

  const int totalSize = 
    tensorSize + numMaxwellModels*tensorSize + numMaxwellModels;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  for (int i=0; i < totalSize; ++i)
    stateValues[i] = dbValues[i];
} // _dbToStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_calcDensity(PylithScalar* const density,
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
pylith::materials::GenMaxwellQpQsIsotropic3D::_calcStressElastic(
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
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStrainSize);

  const PylithScalar mu = properties[p_muEff];
  const PylithScalar k = properties[p_kEff];
  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar lambda = k - 2.0/3.0 * mu;

  // :TODO: Need to consider initial state variables????
  const PylithScalar e11 = totalStrain[0] - initialStrain[0];
  const PylithScalar e22 = totalStrain[1] - initialStrain[1];
  const PylithScalar e33 = totalStrain[2] - initialStrain[2];
  const PylithScalar e12 = totalStrain[3] - initialStrain[3];
  const PylithScalar e23 = totalStrain[4] - initialStrain[4];
  const PylithScalar e13 = totalStrain[5] - initialStrain[5];
  
  const PylithScalar s123 = lambda * (e11 + e22 + e33);

  stress[0] = s123 + mu2*e11 + initialStress[0];
  stress[1] = s123 + mu2*e22 + initialStress[1];
  stress[2] = s123 + mu2*e33 + initialStress[2];
  stress[3] = mu2 * e12 + initialStress[3];
  stress[4] = mu2 * e23 + initialStress[4];
  stress[5] = mu2 * e13 + initialStress[5];

  PetscLogFlops(28);
} // _calcStressElastic


// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as a viscoelastic
// material.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_calcStressViscoelastic(
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
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStrainSize);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == _tensorSize);

  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;
  const int tensorSize = _GenMaxwellQpQsIsotropic3D::tensorSize;

  const PylithScalar mu = properties[p_muEff];
  const PylithScalar bulkModulus = properties[p_kEff];

  const PylithScalar mu2 = 2.0 * mu;

  const PylithScalar e11 = totalStrain[0] - initialStrain[0];
  const PylithScalar e22 = totalStrain[1] - initialStrain[1];
  const PylithScalar e33 = totalStrain[2] - initialStrain[2];
  
  const PylithScalar e123 = e11 + e22 + e33;
  const PylithScalar meanStrainTpdt = e123 / 3.0;

  PylithScalar elasFracShear = 1.0;  // deviatoric (shear) component
  PylithScalar elasFracBulk = 1.0; // mean (bulk)  component
  for (int imodel=0; imodel < numMaxwellModels; ++imodel) {
    elasFracShear -= properties[p_shearRatio+imodel];
    elasFracBulk -= properties[p_bulkRatio+imodel];
  } // for
  const PylithScalar tolerance = 1.0e-6;
  assert(elasFracShear >= -tolerance);
  assert(elasFracBulk >= -tolerance);
  
  PetscLogFlops(7 + 2*numMaxwellModels);

  // Get viscous strains (deviatoric + mean).
  // Current viscous strains are in _viscousDevStrain and _viscousMeanStrain.
  if (computeStateVars) {
    _computeStateVars(stateVars, numStateVars,
		      properties, numProperties,
		      totalStrain, strainSize,
		      initialStress, initialStressSize,
		      initialStrain, initialStrainSize);
  } else {
    for (int iModel=0; iModel < numMaxwellModels; ++iModel) {
      for (int i=0; i < tensorSize; ++i)
	_viscousDevStrain[iModel*tensorSize+i] = 
	  stateVars[s_viscousDevStrain+iModel*tensorSize+i];

      _viscousMeanStrain[iModel] = stateVars[s_viscousMeanStrain+iModel];
    } // for
  } // if/else


  // Compute mean stresses.
  PylithScalar meanStrain = elasFracBulk * meanStrainTpdt;
  for (int iModel=0; iModel < numMaxwellModels; ++iModel)
    meanStrain += _viscousMeanStrain[iModel];
  const PylithScalar meanStressTpdt = 3.0 * bulkModulus * meanStrain;
  
  // Compute stresses (mean + deviatoric)
  assert(6 == tensorSize);
  const PylithScalar diag[6] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
  for (int i=0; i < tensorSize; ++i) {
    PylithScalar devStrain = elasFracShear * 
      (totalStrain[i] - initialStrain[i] - diag[i]*meanStrainTpdt);
    for (int iModel=0; iModel < numMaxwellModels; ++iModel)
      devStrain += _viscousDevStrain[iModel*tensorSize+i];
    stress[i] = diag[i]*meanStressTpdt + mu2 * devStrain + initialStress[i];
  } // for

  PetscLogFlops(3 + numMaxwellModels*1 + tensorSize*(6 + numMaxwellModels));
} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_calcElasticConstsElastic(
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
  assert(_GenMaxwellQpQsIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStrainSize);
 
  const PylithScalar mu = properties[p_muEff];
  const PylithScalar k = properties[p_kEff];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar lambda2mu = k + 4.0/3.0*mu;
  const PylithScalar lambda = k - 2.0/3.0*mu;

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

  PetscLogFlops(7);
} // _calcElasticConstsElastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties
// as a viscoelastic material.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_calcElasticConstsViscoelastic(
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
  assert(_GenMaxwellQpQsIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStrainSize);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == _tensorSize);

  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;

  const PylithScalar mu = properties[p_muEff];
  const PylithScalar bulkModulus = properties[p_kEff];

  // Compute viscous contribution. (deviatoric + mean)
  PylithScalar elasFracShear = 1.0;  // deviatoric (shear) component
  PylithScalar visFactorDev = 0.0;
  PylithScalar elasFracBulk = 1.0; // mean (bulk) component
  PylithScalar visFactorBulk = 0.0;
  for (int iModel=0; iModel < numMaxwellModels; ++iModel) {
    const PylithScalar shearRatio = properties[p_shearRatio+iModel];
    const PylithScalar bulkRatio = properties[p_bulkRatio+iModel];
    elasFracShear -= shearRatio;
    elasFracBulk -= bulkRatio;

    const PylithScalar maxwellTimeShear = properties[p_maxwellTimeShear+iModel];
    visFactorDev +=
      shearRatio*ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTimeShear);
    const PylithScalar maxwellTimeBulk = properties[p_maxwellTimeBulk+iModel];
    visFactorBulk +=
      bulkRatio*ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTimeBulk);
  } // for
  const PylithScalar tolerance = 1.0e-6;
  assert(elasFracShear >= -tolerance);
  assert(elasFracBulk >= -tolerance);

  const PylithScalar muEff = mu * (elasFracShear + visFactorDev);
  const PylithScalar kEff = bulkModulus * (elasFracBulk + visFactorBulk);
  const PylithScalar mu2Eff = 2.0*muEff;
  const PylithScalar lambda2muEff = kEff + 4.0/3.0*muEff;
  const PylithScalar lambdaEff = kEff - 2.0/3.0*muEff;

  elasticConsts[ 0] = lambda2muEff; // C1111
  elasticConsts[ 1] = lambdaEff; // C1122
  elasticConsts[ 2] = lambdaEff; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = 0; // C1123
  elasticConsts[ 5] = 0; // C1113
  elasticConsts[ 6] = lambdaEff; // C2211
  elasticConsts[ 7] = lambda2muEff; // C2222
  elasticConsts[ 8] = lambdaEff; // C2233
  elasticConsts[ 9] = 0; // C2212
  elasticConsts[10] = 0; // C2223
  elasticConsts[11] = 0; // C2213
  elasticConsts[12] = lambdaEff; // C3311
  elasticConsts[13] = lambdaEff; // C3322
  elasticConsts[14] = lambda2muEff; // C3333
  elasticConsts[15] = 0; // C3312
  elasticConsts[16] = 0; // C3323
  elasticConsts[17] = 0; // C3313
  elasticConsts[18] = 0; // C1211
  elasticConsts[19] = 0; // C1222
  elasticConsts[20] = 0; // C1233
  elasticConsts[21] = mu2Eff; // C1212 // ??
  elasticConsts[22] = 0; // C1223
  elasticConsts[23] = 0; // C1213
  elasticConsts[24] = 0; // C2311
  elasticConsts[25] = 0; // C2322
  elasticConsts[26] = 0; // C2333
  elasticConsts[27] = 0; // C2312
  elasticConsts[28] = mu2Eff; // C2323
  elasticConsts[29] = 0; // C2313
  elasticConsts[30] = 0; // C1311
  elasticConsts[31] = 0; // C1322
  elasticConsts[32] = 0; // C1333
  elasticConsts[33] = 0; // C1312
  elasticConsts[34] = 0; // C1323
  elasticConsts[35] = mu2Eff; // C1313

#if 0 // DEBUGGING
  std::cout << "_calcElasticConstsViscoelastic" << std::endl
	    << "  mu: " << muEff
	    << ", muEff: " << muEff
	    << ", k: " << bulkModulus
	    << ", kEff: " << kEff
	    << ", lambda2muEff: " << lambda2muEff
	    << ", lambdaEff: " << lambdaEff
	    << std::endl;
#endif

  PetscLogFlops(1 + numMaxwellModels*6 + 11);
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_updateStateVarsElastic(
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
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _GenMaxwellQpQsIsotropic3D::tensorSize;
  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;

  const PylithScalar e11 = totalStrain[0];
  const PylithScalar e22 = totalStrain[1];
  const PylithScalar e33 = totalStrain[2];
  const PylithScalar meanStrainTpdt = (e11 + e22 + e33) / 3.0;

  // Update total strain
  for (int i=0; i < tensorSize; ++i)
    stateVars[s_totalStrain+i] = totalStrain[i];

  // Initialize viscous strains to deviatoric elastic strains.

  // Deviatoric viscous strains.
  assert(tensorSize == 6);
  const PylithScalar diag[6] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
  for (int iModel=0; iModel < numMaxwellModels; ++iModel)
    for (int i=0; i < tensorSize; ++i) {
      const PylithScalar devStrain = totalStrain[i] - diag[i]*meanStrainTpdt;
      stateVars[s_viscousDevStrain +iModel*tensorSize+i] = 
	properties[p_shearRatio+iModel] * devStrain;
    } // for

  // Mean viscous strains.
  for (int iModel=0; iModel < numMaxwellModels; ++iModel)  
    stateVars[s_viscousMeanStrain+iModel] = 
      properties[p_bulkRatio+iModel] * meanStrainTpdt;

  PetscLogFlops(3 + 2 * tensorSize);

  _needNewJacobian = true;
} // _updateStateVarsElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_updateStateVarsViscoelastic(
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
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStrainSize);

  _computeStateVars(stateVars, numStateVars,
		    properties, numProperties,
		    totalStrain, strainSize,
		    initialStress, initialStressSize,
		    initialStrain, initialStrainSize);

  const int tensorSize = _GenMaxwellQpQsIsotropic3D::tensorSize;
  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;

  // Total strain
  for (int i=0; i < tensorSize; ++i)
    stateVars[s_totalStrain+i] = totalStrain[i];

  // Viscous deviatoric strains.
  for (int iModel=0; iModel < numMaxwellModels; ++iModel)
    for (int i=0; i < tensorSize; ++i)
      stateVars[s_viscousDevStrain+iModel*tensorSize+i] = 
	_viscousDevStrain[iModel*tensorSize+i];

  for (int iModel=0; iModel < numMaxwellModels; ++iModel)
    stateVars[s_viscousMeanStrain+iModel] = 
      _viscousMeanStrain[iModel];

  _needNewJacobian = false;
} // _updateStateVarsViscoelastic

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::GenMaxwellQpQsIsotropic3D::_stableTimeStepImplicit(
					   const PylithScalar* properties,
					   const int numProperties,
					   const PylithScalar* stateVars,
					   const int numStateVars) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);

  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;

  PylithScalar dtStable = pylith::PYLITH_MAXSCALAR;

  for (int i=0; i < numMaxwellModels; ++i) {
    const PylithScalar maxwellTime = properties[p_maxwellTimeShear+i];
    const PylithScalar dt = 0.2*maxwellTime;
    if (dt < dtStable)
      dtStable = dt;
  } // for

  for (int i=0; i < numMaxwellModels; ++i) {
    const PylithScalar maxwellTime = properties[p_maxwellTimeBulk+i];
    const PylithScalar dt = 0.2*maxwellTime;
    if (dt < dtStable)
      dtStable = dt;
  } // for

  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Get stable time step for explicit time integration.
PylithScalar
pylith::materials::GenMaxwellQpQsIsotropic3D::_stableTimeStepExplicit(const PylithScalar* properties,
								      const int numProperties,
								      const PylithScalar* stateVars,
								      const int numStateVars,
								      const double minCellWidth) const
{ // _stableTimeStepExplicit
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
 
  const PylithScalar mu = properties[p_muEff];
  const PylithScalar kappa = properties[p_kEff];
  const PylithScalar density = properties[p_density];

  assert(density > 0.0);
  const PylithScalar vp = sqrt((kappa + 4.0/3.0*mu) / density);

  const PylithScalar dtStable = minCellWidth / vp;
  return dtStable;
} // _stableTimeStepExplicit


// ----------------------------------------------------------------------
// Compute viscous strain for current time step.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_computeStateVars(
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
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _GenMaxwellQpQsIsotropic3D::tensorSize;
  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;

  // :TODO: Need to account for initial values for state variables
  // and the initial strain??

  const PylithScalar e11 = totalStrain[0];
  const PylithScalar e22 = totalStrain[1];
  const PylithScalar e33 = totalStrain[2];
  const PylithScalar meanStrainTpdt = (e11 + e22 + e33) / 3.0;

  const PylithScalar meanStrainT = 
    ( stateVars[s_totalStrain+0] +
      stateVars[s_totalStrain+1] +
      stateVars[s_totalStrain+2] ) / 3.0;
  
  PetscLogFlops(6);

  // Deviatoric viscous strains.

  assert(6 == tensorSize);
  const PylithScalar diag[6] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
  for (int iModel=0; iModel < numMaxwellModels; ++iModel) {

    const PylithScalar dq = 
      ViscoelasticMaxwell::viscousStrainParam(_dt, properties[p_maxwellTimeShear+iModel]);

    for (int i=0; i < tensorSize; ++i) {
      const PylithScalar devStrainTpdt = totalStrain[i] - diag[i]*meanStrainTpdt;
      const PylithScalar devStrainT = stateVars[s_totalStrain+i] - diag[i]*meanStrainT;
      const PylithScalar deltaStrain = devStrainTpdt - devStrainT;
      
      _viscousDevStrain[iModel*tensorSize+i] = 
	exp(-_dt/properties[p_maxwellTimeShear+iModel]) * 
	stateVars[s_viscousDevStrain+iModel*tensorSize+i] + 
	properties[p_shearRatio+iModel] * dq * deltaStrain;
    } // for
  } // for


  // Mean viscous strains.

  // Compute Prony series terms
  for (int iModel=0; iModel < numMaxwellModels; ++iModel) {

    const PylithScalar dq = 
      ViscoelasticMaxwell::viscousStrainParam(_dt, properties[p_maxwellTimeBulk+iModel]);

    const PylithScalar deltaStrain = meanStrainTpdt - meanStrainT;

    _viscousMeanStrain[iModel] =  
      exp(-_dt/properties[p_maxwellTimeBulk+iModel]) * 
      stateVars[s_viscousMeanStrain+iModel] + 
      properties[p_bulkRatio+iModel] * dq * deltaStrain;
  } // for

} // _computeStateVars


// End of file 
