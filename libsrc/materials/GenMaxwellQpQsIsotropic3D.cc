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

#include "GenMaxwellQpQsIsotropic3D.hh" // implementation of object methods

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
	{ "shear_ratio", numMaxwellModels, pylith::topology::FieldBase::OTHER },
	{ "bulk_ratio", numMaxwellModels, pylith::topology::FieldBase::OTHER },
	{ "maxwell_time_shear", numMaxwellModels, pylith::topology::FieldBase::OTHER },
	{ "maxwell_time_bulk", numMaxwellModels, pylith::topology::FieldBase::OTHER },
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
      const int numStateVars = 1+2*numMaxwellModels;
      
      /// State variables. :KLUDGE: Not generalized over number of models.
      const Metadata::ParamDescription stateVars[numStateVars] = {
	{ "total_strain", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "viscous_strain_1", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "viscous_strain_2", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "viscous_strain_3", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "viscous_strain_1_bulk", 1, pylith::topology::FieldBase::SCALAR },
	{ "viscous_strain_2_bulk", 1, pylith::topology::FieldBase::SCALAR },
	{ "viscous_strain_3_bulk", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in state variables spatial database
      const int numDBStateVars = tensorSize * (numMaxwellModels+1) + numMaxwellModels;
      const char* dbStateVars[numDBStateVars] = {"total-strain-xx",
						 "total-strain-yy",
						 "total-strain-zz",
						 "total-strain-xy",
						 "total-strain-yz",
						 "total-strain-xz",
						 "viscous-strain-1-xx",
						 "viscous-strain-1-yy",
						 "viscous-strain-1-zz",
						 "viscous-strain-1-xy",
						 "viscous-strain-1-yz",
						 "viscous-strain-1-xz",
						 "viscous-strain-2-xx",
						 "viscous-strain-2-yy",
						 "viscous-strain-2-zz",
						 "viscous-strain-2-xy",
						 "viscous-strain-2-yz",
						 "viscous-strain-2-xz",
						 "viscous-strain-3-xx",
						 "viscous-strain-3-yy",
						 "viscous-strain-3-zz",
						 "viscous-strain-3-xy",
						 "viscous-strain-3-yz",
						 "viscous-strain-3-xz",
						 "viscous-strain-1-bulk",
						 "viscous-strain-3-bulk",
						 "viscous-strain-3-bulk",
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

const int pylith::materials::GenMaxwellQpQsIsotropic3D::p_bulkRatio =
  pylith::materials::GenMaxwellQpQsIsotropic3D::p_shearRatio +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::numMaxwellModels;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::p_maxwellTimeShear =
  pylith::materials::GenMaxwellQpQsIsotropic3D::p_bulkRatio +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::numMaxwellModels;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::p_maxwellTimeBulk =
  pylith::materials::GenMaxwellQpQsIsotropic3D::p_maxwellTimeShear +
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

const int pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousStrain1 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::s_totalStrain +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousStrain2 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousStrain1 +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousStrain3 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousStrain2 +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousStrainBulk1 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousStrain3 +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousStrainBulk2 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousStrainBulk1 + 1;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousStrainBulk3 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::s_viscousStrainBulk2 + 1;

// Indices of database values (order must match dbStateVars)
const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_totalStrain = 0;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousStrain1 =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_totalStrain +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousStrain2 =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousStrain1 +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousStrain3 =
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousStrain2 +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousStrainBulk1 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousStrain3 +
  pylith::materials::_GenMaxwellQpQsIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousStrainBulk2 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousStrainBulk1 + 1;

const int pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousStrainBulk3 = 
  pylith::materials::GenMaxwellQpQsIsotropic3D::db_viscousStrainBulk2 + 1;

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
  useElasticBehavior(true);
  _viscousStrain.resize(_GenMaxwellQpQsIsotropic3D::numMaxwellModels*_tensorSize);
  _viscousStrainBulk.resize(_GenMaxwellQpQsIsotropic3D::numMaxwellModels);
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
					    double* const propValues,
					    const double_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_GenMaxwellQpQsIsotropic3D::numDBProperties == numDBValues);

  const double density = dbValues[db_density];
  const double vs = dbValues[db_vs];
  const double vp = dbValues[db_vp];
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

  const double mu = density * vs*vs;
  const double k = density * vp*vp - 4.0/3.0*mu;

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

  double visFrac = 0.0;
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

  double bulkFrac = 0.0;
  for (int imodel = 0; imodel < numMaxwellModels; ++imodel) 
    bulkFrac += dbValues[db_bulkRatio + imodel];
  if (bulkFrac > 1.0) {
    std::ostringstream msg;
    msg << "Bulk modulus ratios sum to a value greater than 1.0 for\n"
	<< "Generalized Maxwell bulk model.\n"
	<< "Ratio 1: " << dbValues[db_bulkRatio  ] << "\n"
	<< "Ratio 2: " << dbValues[db_bulkRatio+1] << "\n"
	<< "Ratio 3: " << dbValues[db_bulkRatio+2] << "\n"
	<< "Total:   " << bulkFrac << "\n";
    throw std::runtime_error(msg.str());
  } // if

  // Loop over number of Maxwell models.
  for (int imodel=0; imodel < numMaxwellModels; ++imodel) {
    double shearRatio = dbValues[db_shearRatio + imodel];
    double bulkRatio = dbValues[db_bulkRatio + imodel];
    double shearViscosity = dbValues[db_shearViscosity + imodel];
    double bulkViscosity = dbValues[db_bulkViscosity + imodel];
    double muFac = shearRatio*mu;
    double kFac = bulkRatio*k;
    double maxwellTimeShear = pylith::PYLITH_MAXDOUBLE;
    double maxwellTimeBulk = pylith::PYLITH_MAXDOUBLE;
    if (muFac > 0.0)
      maxwellTimeShear = shearViscosity / muFac;
    if (kFac > 0.0)
      maxwellTimeBulk = bulkViscosity / kFac;
    if (muFac < 0.0 || shearRatio < 0.0 || shearViscosity < 0.0 || maxwellTimeShear < 0.0 || 
	kFac < 0.0 || bulkRatio < 0.0 || bulkViscosity < 0.0 || maxwellTimeBulk < 0.0) {
      std::ostringstream msg;
      msg << "Found negative value(s) for physical properties.\n"
	  << "shearRatio: " << shearRatio << "\n"
	  << "muFac: " << muFac << "\n"
	  << "shear viscosity: " << shearViscosity << "\n"
	  << "maxwellTimeShear: " << maxwellTimeShear << "\n"
	  << "bulkRatio: " << bulkRatio << "\n"
	  << "kFac: " << kFac << "\n"
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
pylith::materials::GenMaxwellQpQsIsotropic3D::_nondimProperties(double* const values,
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
pylith::materials::GenMaxwellQpQsIsotropic3D::_dimProperties(double* const values,
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
					double* const stateValues,
					const double_array& dbValues)
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_GenMaxwellQpQsIsotropic3D::numDBStateVars == numDBValues);

  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;
  const int totalSize = _tensorSize * (numMaxwellModels+1) + numMaxwellModels;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  for (int i=0; i < totalSize; ++i)
    stateValues[i] = dbValues[i];

  PetscLogFlops(0);
} // _dbToStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_calcDensity(double* const density,
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
pylith::materials::GenMaxwellQpQsIsotropic3D::_calcStressElastic(
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

  const double mu = properties[p_muEff];
  const double k = properties[p_kEff];
  const double mu2 = 2.0 * mu;
  const double lambda = k - 2.0/3.0 * mu;

  // :TODO: Need to consider initial state variables????
  const double e11 = totalStrain[0] - initialStrain[0];
  const double e22 = totalStrain[1] - initialStrain[1];
  const double e33 = totalStrain[2] - initialStrain[2];
  const double e12 = totalStrain[3] - initialStrain[3];
  const double e23 = totalStrain[4] - initialStrain[4];
  const double e13 = totalStrain[5] - initialStrain[5];
  
  const double s123 = lambda * (e11 + e22 + e33);

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
pylith::materials::GenMaxwellQpQsIsotropic3D::_calcStressViscoelastic(
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

  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;
  const int tensorSize = _tensorSize;

  const double mu = properties[p_muEff];
  const double bulkModulus = properties[p_kEff];
  const double shearRatio[numMaxwellModels] = {
    properties[p_shearRatio  ],
    properties[p_shearRatio+1],
    properties[p_shearRatio+2]
  };
  const double bulkRatio[numMaxwellModels] = {
    properties[p_bulkRatio  ],
    properties[p_bulkRatio+1],
    properties[p_bulkRatio+2]
  };

  const double mu2 = 2.0 * mu;

  // :TODO: Need to determine how to incorporate initial strain and
  // state variables
  const double e11 = totalStrain[0] - initialStrain[0];
  const double e22 = totalStrain[1] - initialStrain[1];
  const double e33 = totalStrain[2] - initialStrain[2];
  
  const double e123 = e11 + e22 + e33;
  const double meanStrainTpdt = e123 / 3.0;

  double elasFracShear = 1.0;  // deviatoric component
  double elasFracBulk = 1.0; // volumetric component
  for (int imodel=0; imodel < numMaxwellModels; ++imodel) {
    elasFracShear -= shearRatio[imodel];
    elasFracBulk -= bulkRatio[imodel];
  } // for
  const double tolerance = 1.0e-6;
  assert(elasFracShear >= -tolerance);
  assert(elasFracBulk >= -tolerance);
  
  PetscLogFlops(7 + 2*numMaxwellModels);

  // Get viscous strains ( deviatoric and volumetric )
  // Current viscous strains are in _viscousStrain and _viscousStrainBulk
  if (computeStateVars) {
    _computeStateVars(stateVars, numStateVars,
		      properties, numProperties,
		      totalStrain, strainSize,
		      initialStress, initialStressSize,
		      initialStrain, initialStrainSize);
  } else {
    for (int i=0; i < tensorSize; ++i)
      _viscousStrain[0*tensorSize+i] = stateVars[s_viscousStrain1+i];
    for (int i=0; i < tensorSize; ++i)
      _viscousStrain[1*tensorSize+i] = stateVars[s_viscousStrain2+i];
    for (int i=0; i < tensorSize; ++i)
      _viscousStrain[2*tensorSize+i] = stateVars[s_viscousStrain3+i];
    _viscousStrainBulk[0] = stateVars[s_viscousStrainBulk1];
    _viscousStrainBulk[1] = stateVars[s_viscousStrainBulk2];
    _viscousStrainBulk[2] = stateVars[s_viscousStrainBulk3];
  } // if/else


  // Compute new stresses ( volumetric + deviatoric )
  const double volStressTpdt = 3.0 * bulkModulus * 
    (elasFracBulk * meanStrainTpdt + 
     bulkRatio[0] * _viscousStrainBulk[0] +
     bulkRatio[1] * _viscousStrainBulk[1] +
     bulkRatio[2] * _viscousStrainBulk[2]);
  
  // Axial components
  for (int i=0; i < 3; ++i) {
    const double devStressTpdt = mu2 * 
      (elasFracShear * (totalStrain[i] - meanStrainTpdt) +
       shearRatio[0] * _viscousStrain[0*tensorSize+i] +
       shearRatio[1] * _viscousStrain[1*tensorSize+i] +
       shearRatio[2] * _viscousStrain[2*tensorSize+i]);
    stress[i] = volStressTpdt + devStressTpdt + initialStress[i];
  } // for

  // shear components
  for (int i=3; i < 6; ++i) {
    const double devStressTpdt = mu2 * 
      (elasFracShear * totalStrain[i] +
       shearRatio[0] * _viscousStrain[0*tensorSize+i] +
       shearRatio[1] * _viscousStrain[1*tensorSize+i] +
       shearRatio[2] * _viscousStrain[2*tensorSize+i]);
    stress[i] = devStressTpdt + initialStress[i];
  } // for

  PetscLogFlops(8 + 3*12 + 3*9);
} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_calcElasticConstsElastic(
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
 
  const double mu = properties[p_muEff];
  const double k = properties[p_kEff];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = k + 4.0/3.0*mu;
  const double lambda = k - 2.0/3.0*mu;

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

  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;
  const int tensorSize = _tensorSize;

  const double mu = properties[p_muEff];
  const double bulkModulus = properties[p_kEff];
  const double mu2 = 2.0 * mu;

  // Compute viscous contribution. (deviatoric + volumetric)
  double elasFracShear = 1.0;  // deviatoric component
  double visFactorShear = 0.0;
  double elasFracBulk = 1.0; // volumetric component
  double visFactorBulk = 0.0;
  for (int imodel = 0; imodel < numMaxwellModels; ++imodel) {
    const double shearRatio = properties[p_shearRatio+imodel];
    const double bulkRatio = properties[p_bulkRatio+imodel];
    elasFracShear -= shearRatio;
    elasFracBulk -= bulkRatio;

    double maxwellTimeShear = pylith::PYLITH_MAXDOUBLE;
    double maxwellTimeBulk = pylith::PYLITH_MAXDOUBLE;
    if (shearRatio != 0.0) {
      maxwellTimeShear = properties[p_maxwellTimeShear+imodel];
      visFactorShear +=
	shearRatio*ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTimeShear);
    } // if
    if (bulkRatio != 0.0) {
      maxwellTimeBulk = properties[p_maxwellTimeBulk+imodel];
      visFactorBulk +=
	bulkRatio*ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTimeBulk);
    } // if
  } // for
  const double tolerance = 1.0e-6;
  assert(elasFracShear >= -tolerance);
  assert(elasFracBulk >= -tolerance);

  double shearFactor = elasFracShear + visFactorShear;
  double bulkFactor = elasFracBulk + visFactorBulk;

  elasticConsts[ 0] = bulkModulus*bulkFactor + 4.0*mu/3.0*shearFactor; // C1111
  elasticConsts[ 1] = bulkModulus*bulkFactor - 2.0*mu/3.0*shearFactor; // C1122
  elasticConsts[ 2] = elasticConsts[1]; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = 0; // C1123
  elasticConsts[ 5] = 0; // C1113
  elasticConsts[ 6] = elasticConsts[1]; // C2211
  elasticConsts[ 7] = elasticConsts[0]; // C2222
  elasticConsts[ 8] = elasticConsts[1]; // C2233
  elasticConsts[ 9] = 0; // C2212
  elasticConsts[10] = 0; // C2223
  elasticConsts[11] = 0; // C2213
  elasticConsts[12] = elasticConsts[1]; // C3311
  elasticConsts[13] = elasticConsts[1]; // C3322
  elasticConsts[14] = elasticConsts[0]; // C3333
  elasticConsts[15] = 0; // C3312
  elasticConsts[16] = 0; // C3323
  elasticConsts[17] = 0; // C3313
  elasticConsts[18] = 0; // C1211
  elasticConsts[19] = 0; // C1222
  elasticConsts[20] = 0; // C1233
  elasticConsts[21] = 2.0*mu*shearFactor; // C1212 // ??
  elasticConsts[22] = 0; // C1223
  elasticConsts[23] = 0; // C1213
  elasticConsts[24] = 0; // C2311
  elasticConsts[25] = 0; // C2322
  elasticConsts[26] = 0; // C2333
  elasticConsts[27] = 0; // C2312
  elasticConsts[28] = elasticConsts[21]; // C2323
  elasticConsts[29] = 0; // C2313
  elasticConsts[30] = 0; // C1311
  elasticConsts[31] = 0; // C1322
  elasticConsts[32] = 0; // C1333
  elasticConsts[33] = 0; // C1312
  elasticConsts[34] = 0; // C1323
  elasticConsts[35] = elasticConsts[21]; // C1313

#if 0 // DEBUGGING
  std::cout << "_calcElasticConstsViscoelastic" << std::endl;
  std::cout << elasticConsts[0] << "  " << elasticConsts[1] << "  " << elasticConsts[2] << std::endl;
  std::cout << elasticConsts[6] << "  " << elasticConsts[7] << std::endl;
  std::cout << elasticConsts[11] << std::endl;
  std::cout << elasticConsts[15] << std::endl;
  std::cout << elasticConsts[18] << std::endl;
  std::cout << elasticConsts[20] << std::endl;
#endif

  PetscLogFlops(15 + 6 * numMaxwellModels);
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_updateStateVarsElastic(
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
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double meanStrainTpdt = (e11 + e22 + e33) / 3.0;

  // Update total strain
  for (int i=0; i < tensorSize; ++i)
    stateVars[s_totalStrain+i] = totalStrain[i];

  // Initialize viscous strains to deviatoric elastic strains.

  // Axial strains
  for (int i=0; i < 3; ++i) {
    const double devStrain = totalStrain[i] - meanStrainTpdt;
    stateVars[s_viscousStrain1+i] = devStrain;
    stateVars[s_viscousStrain2+i] = devStrain;
    stateVars[s_viscousStrain3+i] = devStrain;
  } // for

  // Shear strains
  for (int i=3; i < 6; ++i) {
    const double devStrain = totalStrain[i];
    stateVars[s_viscousStrain1+i] = devStrain;
    stateVars[s_viscousStrain2+i] = devStrain;
    stateVars[s_viscousStrain3+i] = devStrain;
  } // for

  // Bulk strains
  stateVars[s_viscousStrainBulk1] = meanStrainTpdt;
  stateVars[s_viscousStrainBulk2] = meanStrainTpdt;
  stateVars[s_viscousStrainBulk3] = meanStrainTpdt;

  PetscLogFlops(3 + 2 * tensorSize);

  _needNewJacobian = true;
} // _updateStateVarsElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_updateStateVarsViscoelastic(
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

  const int tensorSize = _tensorSize;

  // Total strain
  for (int i=0; i < tensorSize; ++i)
    stateVars[s_totalStrain+i] = totalStrain[i];

  for (int i=0; i < tensorSize; ++i)
    stateVars[s_viscousStrain1+i] = _viscousStrain[0*tensorSize+i];
  for (int i=0; i < tensorSize; ++i)
    stateVars[s_viscousStrain2+i] = _viscousStrain[1*tensorSize+i];
  for (int i=0; i < tensorSize; ++i)
    stateVars[s_viscousStrain3+i] = _viscousStrain[2*tensorSize+i];

  stateVars[s_viscousStrainBulk1] = _viscousStrainBulk[0];
  stateVars[s_viscousStrainBulk2] = _viscousStrainBulk[1];
  stateVars[s_viscousStrainBulk3] = _viscousStrainBulk[2];

  _needNewJacobian = false;
} // _updateStateVarsViscoelastic

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::GenMaxwellQpQsIsotropic3D::_stableTimeStepImplicit(
					   const double* properties,
					   const int numProperties,
					   const double* stateVars,
					   const int numStateVars) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);

  double dtStable = pylith::PYLITH_MAXDOUBLE;

  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;
  for (int i=0; i < numMaxwellModels; ++i) {
    const double maxwellTime = properties[p_maxwellTimeShear+i];
    const double dt = 0.2*maxwellTime;
    if (dt < dtStable)
      dtStable = dt;
  } // for

  for (int i=0; i < numMaxwellModels; ++i) {
    const double maxwellTime = properties[p_maxwellTimeBulk+i];
    const double dt = 0.2*maxwellTime;
    if (dt < dtStable)
      dtStable = dt;
  } // for

  return dtStable;
} // _stableTimeStepImplicit

#include <iostream>
// ----------------------------------------------------------------------
// Compute viscous strain for current time step.
void
pylith::materials::GenMaxwellQpQsIsotropic3D::_computeStateVars(
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
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellQpQsIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;
  const int numMaxwellModels = _GenMaxwellQpQsIsotropic3D::numMaxwellModels;

  const double shearRatio[numMaxwellModels] = {
    properties[p_shearRatio  ],
    properties[p_shearRatio+1],
    properties[p_shearRatio+2]
  };
  const double bulkRatio[numMaxwellModels] = {
    properties[p_bulkRatio  ],
    properties[p_bulkRatio+1],
    properties[p_bulkRatio+2]
  };
  const double maxwellTimeShear[numMaxwellModels] = {
    properties[p_maxwellTimeShear  ],
    properties[p_maxwellTimeShear+1],
    properties[p_maxwellTimeShear+2]
  };
  const double maxwellTimeBulk[numMaxwellModels] = {
    properties[p_maxwellTimeBulk  ],
    properties[p_maxwellTimeBulk+1],
    properties[p_maxwellTimeBulk+2]
  };

  // :TODO: Need to account for initial values for state variables
  // and the initial strain??

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  
  const double meanStrainTpdt = (e11 + e22 + e33) / 3.0;

  const double meanStrainT = 
    ( stateVars[s_totalStrain+0] +
      stateVars[s_totalStrain+1] +
      stateVars[s_totalStrain+2] ) / 3.0;
  
  PetscLogFlops(6);

  // Compute Prony series terms
  double dq[numMaxwellModels];
  for (int i=0; i < numMaxwellModels; ++i) {
    dq[i] = 0.0;
    if (shearRatio[i] != 0.0)
      dq[i] = ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTimeShear[i]);
  } // if

  // Compute new viscous strains (deviatoric)
  const double diag[6] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
  for (int i=0; i < tensorSize; ++i) {
    const double devStrainTpdt = totalStrain[i] - diag[i]*meanStrainTpdt;
    const double devStrainT = stateVars[s_totalStrain+i] - diag[i]*meanStrainT;
    const double deltaStrain = devStrainTpdt - devStrainT;

    // Maxwell model 1
    int imodel = 0;
    if (0.0 != shearRatio[imodel]) {
      _viscousStrain[imodel*tensorSize+i] = exp(-_dt/maxwellTimeShear[imodel])*
	stateVars[s_viscousStrain1+i] + dq[imodel] * deltaStrain;
      PetscLogFlops(6);
    } // if

    // Maxwell model 2
    imodel = 1;
    if (0.0 != shearRatio[imodel]) {
      _viscousStrain[imodel*tensorSize+i] = exp(-_dt/maxwellTimeShear[imodel])*
	stateVars[s_viscousStrain2+i] + dq[imodel] * deltaStrain;
      PetscLogFlops(6);
    } // if

    // Maxwell model 3
    imodel = 2;
    if (0.0 != shearRatio[imodel]) {
      _viscousStrain[imodel*tensorSize+i] = exp(-_dt/maxwellTimeShear[imodel])*
	stateVars[s_viscousStrain3+i] + dq[imodel] * deltaStrain;
      PetscLogFlops(6);
    } // if
  } // for

  PetscLogFlops(5*tensorSize);

  // Compute new viscous strains (volumetric)
  // Compute Prony series terms
  for (int i=0; i < numMaxwellModels; ++i) {
    dq[i] = 0.0;
    if (bulkRatio[i] != 0.0)
      dq[i] = ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTimeBulk[i]);
  } // for

  const double deltaStrain = meanStrainTpdt - meanStrainT;
  PetscLogFlops(1);

  // Maxwell model 1
  int imodel = 0;
  if (0.0 != bulkRatio[imodel]) {
    _viscousStrainBulk[imodel] = exp(-_dt/maxwellTimeBulk[imodel])*
      stateVars[s_viscousStrainBulk1] + dq[imodel] * deltaStrain;
    PetscLogFlops(6);
  } // if

    // Maxwell model 2
  imodel = 1;
  if (0.0 != bulkRatio[imodel]) {
    _viscousStrainBulk[imodel] = exp(-_dt/maxwellTimeBulk[imodel])*
      stateVars[s_viscousStrainBulk2] + dq[imodel] * deltaStrain;
    PetscLogFlops(6);
  } // if

    // Maxwell model 3
  imodel = 2;
  if (0.0 != bulkRatio[imodel]) {
    _viscousStrainBulk[imodel] = exp(-_dt/maxwellTimeBulk[imodel])*
      stateVars[s_viscousStrainBulk3] + dq[imodel] * deltaStrain;
    PetscLogFlops(6);
  } // if


} // _computeStateVars


// End of file 
