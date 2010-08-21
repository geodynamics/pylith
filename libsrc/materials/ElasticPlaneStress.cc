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

#include "ElasticPlaneStress.hh" // implementation of object methods

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
    namespace _ElasticPlaneStress {

      // Dimension of material.
      const int dimension = 2;

      // Number of entries in stress tensor.
      const int tensorSize = 3;

      // Number of elastic constants (for general 2-D elastic material)
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
      
    } // _ElasticPlaneStress
  } // materials
} // pylith

// Indices of physical properties
const int pylith::materials::ElasticPlaneStress::p_density = 0;

const int pylith::materials::ElasticPlaneStress::p_mu = 
  pylith::materials::ElasticPlaneStress::p_density + 1;

const int pylith::materials::ElasticPlaneStress::p_lambda = 
  pylith::materials::ElasticPlaneStress::p_mu + 1;

// Indices of database values (order must match dbProperties)
const int pylith::materials::ElasticPlaneStress::db_density = 0;

const int pylith::materials::ElasticPlaneStress::db_vs = 
  pylith::materials::ElasticPlaneStress::db_density + 1;

const int pylith::materials::ElasticPlaneStress::db_vp = 
  pylith::materials::ElasticPlaneStress::db_vs + 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticPlaneStress::ElasticPlaneStress(void) :
  ElasticMaterial(_ElasticPlaneStress::dimension,
		  _ElasticPlaneStress::tensorSize,
		  _ElasticPlaneStress::numElasticConsts,
		  Metadata(_ElasticPlaneStress::properties,
			   _ElasticPlaneStress::numProperties,
			   _ElasticPlaneStress::dbProperties,
			   _ElasticPlaneStress::numDBProperties,
			   0, 0,
			   0, 0))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticPlaneStress::~ElasticPlaneStress(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticPlaneStress::_dbToProperties(
					  double* propValues,
					  const double_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_ElasticPlaneStress::numDBProperties == numDBValues);

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
pylith::materials::ElasticPlaneStress::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticPlaneStress::numProperties);

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
pylith::materials::ElasticPlaneStress::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticPlaneStress::numProperties);

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
pylith::materials::ElasticPlaneStress::_calcDensity(double* const density,
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
pylith::materials::ElasticPlaneStress::_calcStress(double* const stress,
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
  assert(_ElasticPlaneStress::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStress::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_ElasticPlaneStress::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_ElasticPlaneStress::tensorSize == initialStrainSize);

  const double density = properties[p_density];
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
  const double lambdamu = lambda + mu;
  
  const double e11 = totalStrain[0] - initialStrain[0];
  const double e22 = totalStrain[1] - initialStrain[1];
  const double e12 = totalStrain[2] - initialStrain[2];

  stress[0] = 
    (2.0*mu2*lambdamu * e11 + mu2*lambda * e22) / lambda2mu + initialStress[0];
  stress[1] =
    (mu2*lambda * e11 + 2.0*mu2*lambdamu * e22) / lambda2mu + initialStress[1];
  stress[2] = mu2 * e12 + initialStress[2];

  PetscLogFlops(21);
} // _calcStress

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::ElasticPlaneStress::_calcElasticConsts(
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
  assert(_ElasticPlaneStress::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStress::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_ElasticPlaneStress::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_ElasticPlaneStress::tensorSize == initialStrainSize);
 
  const double density = properties[p_density];
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
  const double c11 = 2.0 * mu2 * (lambda + mu) / lambda2mu;

  elasticConsts[0] = c11; // C1111
  elasticConsts[1] = mu2 * lambda / lambda2mu; // C1122
  elasticConsts[2] = 0; // C1112
  elasticConsts[3] = elasticConsts[1]; // C2211
  elasticConsts[4] = c11; // C2222
  elasticConsts[5] = 0; // C2212
  elasticConsts[6] = 0; // C1211
  elasticConsts[7] = 0; // C1222
  elasticConsts[8] = mu2; // C1212

  PetscLogFlops(8);
} // calcElasticConsts

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::ElasticPlaneStress::stableTimeStepImplicit(
					const topology::Mesh& mesh) {
  return pylith::PYLITH_MAXDOUBLE;
}

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::ElasticPlaneStress::_stableTimeStepImplicit(
				     const double* properties,
				     const int numProperties,
				     const double* stateVars,
				     const int numStateVars) const
{ // _stableTimeStepImplicit
  return pylith::PYLITH_MAXDOUBLE;
} // _stableTimeStepImplicit


// End of file 
