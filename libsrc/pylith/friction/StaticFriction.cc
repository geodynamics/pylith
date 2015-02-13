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

#include "StaticFriction.hh" // implementation of object methods

#include "pylith/materials/Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES MAXSCALAR

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace friction {
    namespace _StaticFriction {

      // Number of physical properties.
      const int numProperties = 2;

      // Physical properties.
      const pylith::materials::Metadata::ParamDescription properties[] = {
	{ "friction_coefficient", 1, pylith::topology::FieldBase::SCALAR },
	{ "cohesion", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 2;
      const char* dbProperties[] = { "friction-coefficient", 
				     "cohesion"
};      
      
    } // _StaticFriction
  } // friction
} // pylith

// Indices of physical properties
const int pylith::friction::StaticFriction::p_coef = 0;
const int pylith::friction::StaticFriction::p_cohesion =
  pylith::friction::StaticFriction::p_coef + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::StaticFriction::db_coef = 0;
const int pylith::friction::StaticFriction::db_cohesion =
  pylith::friction::StaticFriction::db_coef + 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::StaticFriction::StaticFriction(void) :
  FrictionModel(materials::Metadata(_StaticFriction::properties,
				    _StaticFriction::numProperties,
				    _StaticFriction::dbProperties,
				    _StaticFriction::numDBProperties,
				    0, 0,
				    0, 0))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::StaticFriction::~StaticFriction(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::friction::StaticFriction::_dbToProperties(PylithScalar* const propValues,
						  const scalar_array& dbValues) const
{ // _dbToProperties
  assert(propValues);
  const int numDBValues = dbValues.size();
  assert(_StaticFriction::numDBProperties == numDBValues);

  const PylithScalar coef = dbValues[db_coef];
  const PylithScalar cohesion = dbValues[db_cohesion];
 
  if (coef < 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned negative value for coefficient "
	<< "of friction.\n"
	<< "coefficient of friction: " << coef << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[p_coef] = coef;
  propValues[p_cohesion] = cohesion;
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::StaticFriction::_nondimProperties(PylithScalar* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _StaticFriction::numProperties);

  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_cohesion] /= pressureScale;
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::StaticFriction::_dimProperties(PylithScalar* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _StaticFriction::numProperties);

  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_cohesion] *= pressureScale;
} // _dimProperties

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
PylithScalar
pylith::friction::StaticFriction::_calcFriction(const PylithScalar t,
						const PylithScalar slip,
						const PylithScalar slipRate,
						const PylithScalar normalTraction,
						const PylithScalar* properties,
						const int numProperties,
						const PylithScalar* stateVars,
						const int numStateVars)
{ // _calcFriction
  assert(properties);
  assert(_StaticFriction::numProperties == numProperties);
  assert(0 == numStateVars);

  const PylithScalar friction = (normalTraction <= 0.0) ?
    properties[p_cohesion] - properties[p_coef] * normalTraction :
    properties[p_cohesion];

  PetscLogFlops(2);

  return friction;
} // _calcFriction


// ----------------------------------------------------------------------
// Compute derivative of friction with slip from properties and
// state variables.
PylithScalar
pylith::friction::StaticFriction::_calcFrictionDeriv(const PylithScalar t,
						     const PylithScalar slip,
						     const PylithScalar slipRate,
						     const PylithScalar normalTraction,
						     const PylithScalar* properties,
						     const int numProperties,
						     const PylithScalar* stateVars,
						     const int numStateVars)
{ // _calcFrictionDeriv
  return 0.0;
} // _calcFrictionDeriv


// End of file 
