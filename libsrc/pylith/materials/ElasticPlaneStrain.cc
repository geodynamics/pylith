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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "ElasticPlaneStrain.hh" // implementation of object methods

#include "Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticPlaneStrain {

      // Dimension of material.
      const int dimension = 2;

      // Number of entries in stress tensor.
      const int tensorSize = 3;

      // Number of elastic constants (for general 3-D elastic material)
      const int numElasticConsts = 9;

      // Number of physical properties.
      const int numProperties = 3;

      // Physical properties.
      const Metadata::ParamDescription properties[] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 3;
      const char* dbProperties[] = { "density", "vs", "vp" };      
      
    } // _ElasticPlaneStrain
  } // materials
} // pylith

// Indices of physical properties
const int pylith::materials::ElasticPlaneStrain::p_density = 0;

const int pylith::materials::ElasticPlaneStrain::p_mu = 
  pylith::materials::ElasticPlaneStrain::p_density + 1;

const int pylith::materials::ElasticPlaneStrain::p_lambda = 
  pylith::materials::ElasticPlaneStrain::p_mu + 1;

// Indices of database values (order must match dbProperties)
const int pylith::materials::ElasticPlaneStrain::db_density = 0;

const int pylith::materials::ElasticPlaneStrain::db_vs = 
  pylith::materials::ElasticPlaneStrain::db_density + 1;

const int pylith::materials::ElasticPlaneStrain::db_vp = 
  pylith::materials::ElasticPlaneStrain::db_vs + 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticPlaneStrain::ElasticPlaneStrain(void) :
  ElasticMaterial(_ElasticPlaneStrain::dimension,
		  _ElasticPlaneStrain::tensorSize,
		  _ElasticPlaneStrain::numElasticConsts,
		  Metadata(_ElasticPlaneStrain::properties,
			   _ElasticPlaneStrain::numProperties,
			   _ElasticPlaneStrain::dbProperties,
			   _ElasticPlaneStrain::numDBProperties,
			   0, 0,
			   0, 0))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticPlaneStrain::~ElasticPlaneStrain(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticPlaneStrain::_dbToProperties(
				          double* const propValues,
                                          const double_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_ElasticPlaneStrain::numDBProperties == numDBValues);

  const double density = dbValues[db_density];
  const double vs = dbValues[db_vs];
  const double vp = dbValues[db_vp];
 
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
  const double lambda = density * vp*vp - 2.0*mu;

  if (lambda <= 0.0) {
    std::ostringstream msg;
    msg << "Attempted to set Lame's constant lambda to nonpositive value.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[p_density] = density;
  propValues[p_mu] = mu;
  propValues[p_lambda] = lambda;

  PetscLogFlops(6);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::ElasticPlaneStrain::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticPlaneStrain::numProperties);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();

  values[p_density] = 
    _normalizer->nondimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->nondimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->nondimensionalize(values[p_lambda], pressureScale);

  PetscLogFlops(3);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::ElasticPlaneStrain::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticPlaneStrain::numProperties);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();

  values[p_density] = 
    _normalizer->dimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->dimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->dimensionalize(values[p_lambda], pressureScale);

  PetscLogFlops(3);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::ElasticPlaneStrain::_calcDensity(double* const density,
						    const double* properties,
						    const int numProperties,
						    const double* stateVars,
						    const int numStateVars)
{ // calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);

  density[0] = properties[p_density];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties.
void
pylith::materials::ElasticPlaneStrain::_calcStress(double* const stress,
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
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticPlaneStrain::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_ElasticPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_ElasticPlaneStrain::tensorSize == initialStrainSize);

  const double density = properties[p_density];
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];

  const double mu2 = 2.0*mu;

  const double e11 = totalStrain[0] - initialStrain[0];
  const double e22 = totalStrain[1] - initialStrain[1];
  const double e12 = totalStrain[2] - initialStrain[2];

  const double s12 = lambda * (e11 + e22);

  stress[0] = s12 + mu2*e11 + initialStress[0];
  stress[1] = s12 + mu2*e22 + initialStress[1];
  stress[2] = mu2 * e12 + initialStress[2];

  PetscLogFlops(14);
} // _calcStress

// ----------------------------------------------------------------------
// Compute elastic constants at location from properties.
void
pylith::materials::ElasticPlaneStrain::_calcElasticConsts(
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
{ // calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticPlaneStrain::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStrain::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_ElasticPlaneStrain::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_ElasticPlaneStrain::tensorSize == initialStrainSize);
 
  const double density = properties[p_density];
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
   
  elasticConsts[0] = lambda2mu; // C1111
  elasticConsts[1] = lambda; // C1122
  elasticConsts[2] = 0; // C1112
  elasticConsts[3] = lambda; // C2211
  elasticConsts[4] = lambda2mu; // C2222
  elasticConsts[5] = 0; // C2212
  elasticConsts[6] = 0; // C1211
  elasticConsts[7] = 0; // C1222
  elasticConsts[8] = mu2; // C1212

  PetscLogFlops(2);
} // calcElasticConsts

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::ElasticPlaneStrain::stableTimeStepImplicit(
					const topology::Mesh& mesh) {
  return pylith::PYLITH_MAXDOUBLE;
}

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::ElasticPlaneStrain::_stableTimeStepImplicit(
				     const double* properties,
				     const int numProperties,
				     const double* stateVars,
				     const int numStateVars) const
{ // _stableTimeStepImplicit
  return pylith::PYLITH_MAXDOUBLE;
} // _stableTimeStepImplicit


// End of file 
