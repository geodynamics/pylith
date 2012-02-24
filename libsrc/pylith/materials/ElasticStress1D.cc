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

#include "ElasticStress1D.hh" // implementation of object methods

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
    namespace _ElasticStress1D {

      // Dimension of material.
      const int dimension = 1;

      /// Number of entries in stress tensor.
      const int tensorSize = 1;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 1;

      /// Number of physical properties.
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
      
    } // _ElasticStress1D
  } // materials
} // pylith

// Indices of physical properties
const int pylith::materials::ElasticStress1D::p_density = 0;

const int pylith::materials::ElasticStress1D::p_mu = 
  pylith::materials::ElasticStress1D::p_density + 1;

const int pylith::materials::ElasticStress1D::p_lambda = 
  pylith::materials::ElasticStress1D::p_mu + 1;

// Indices of database values (order must match dbProperties)
const int pylith::materials::ElasticStress1D::db_density = 0;

const int pylith::materials::ElasticStress1D::db_vs = 
  pylith::materials::ElasticStress1D::db_density + 1;

const int pylith::materials::ElasticStress1D::db_vp = 
  pylith::materials::ElasticStress1D::db_vs + 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticStress1D::ElasticStress1D(void) :
  ElasticMaterial(_ElasticStress1D::dimension,
		  _ElasticStress1D::tensorSize,
		  _ElasticStress1D::numElasticConsts,
		  Metadata(_ElasticStress1D::properties,
			   _ElasticStress1D::numProperties,
			   _ElasticStress1D::dbProperties,
			   _ElasticStress1D::numDBProperties,
			   0, 0,
			   0, 0))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticStress1D::~ElasticStress1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticStress1D::_dbToProperties(
					    PylithScalar* const propValues,
					    const scalar_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_ElasticStress1D::numDBProperties == numDBValues);

  const PylithScalar density = dbValues[db_density];
  const PylithScalar vs = dbValues[db_vs];
  const PylithScalar vp = dbValues[db_vp];
 
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
  const PylithScalar lambda = density * vp*vp - 2.0*mu;

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
pylith::materials::ElasticStress1D::_nondimProperties(PylithScalar* const values,
						      const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticStress1D::numProperties);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

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
pylith::materials::ElasticStress1D::_dimProperties(PylithScalar* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticStress1D::numProperties);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

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
pylith::materials::ElasticStress1D::_calcDensity(PylithScalar* const density,
						 const PylithScalar* properties,
						 const int numProperties,
						 const PylithScalar* stateVars,
						 const int numStateVars)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);

  density[0] = properties[p_density];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties.
void
pylith::materials::ElasticStress1D::_calcStress(PylithScalar* const stress,
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
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticStress1D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);
  assert(0 != totalStrain);
  assert(_ElasticStress1D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_ElasticStress1D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_ElasticStress1D::tensorSize == initialStrainSize);

  const PylithScalar density = properties[p_density];
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];

  const PylithScalar e11 = totalStrain[0] - initialStrain[0];
  stress[0] = mu * (3.0*lambda+2.0*mu) / (lambda + mu) * e11 + initialStress[0];

  PetscLogFlops(9);
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::ElasticStress1D::_calcElasticConsts(
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
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticStress1D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);
  assert(0 != totalStrain);
  assert(_ElasticStress1D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_ElasticStress1D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_ElasticStress1D::tensorSize == initialStrainSize);
 
  const PylithScalar density = properties[p_density];
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];

  elasticConsts[0] = mu * (3.0*lambda+2.0*mu) / (lambda + mu);

  PetscLogFlops(6);
} // _calcElasticConsts

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::ElasticStress1D::stableTimeStepImplicit(
					const topology::Mesh& mesh) {
  return pylith::PYLITH_MAXSCALAR;
}

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::ElasticStress1D::_stableTimeStepImplicit(
					      const PylithScalar* properties,
					      const int numProperties,
					      const PylithScalar* stateVars,
					      const int numStateVars) const
{ // _stableTimeStepImplicit
  return pylith::PYLITH_MAXSCALAR;
} // _stableTimeStepImplicit


// End of file 
