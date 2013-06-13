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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

// See ViscousFriction.hh for a description of each C++ function and
// its arguments.

#include <portinfo> // machine specific info generated by configure

#include "ViscousFriction.hh" // implementation of object methods

#include "pylith/materials/Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES scaary_array
#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Create a local namespace to use for local constants and other
// information. This insulates all other classes from this information
// while preventing clashes with other local constants and data (as
// long as no other object use the same _ViscousFriction namespace in
// the  namespace.
namespace contrib {
  namespace friction {
    namespace _ViscousFriction {

      // These are the fault constitutive parameters stored during the
      // simulation and need not coincide with the physical properties
      // provided by the user.

      // Number of physical properties.
      const int numProperties = 3;

      // Friction model parameters.
      const pylith::materials::Metadata::ParamDescription properties[numProperties] = {
	{ "static_coefficient", 1, pylith::topology::FieldBase::SCALAR },
	{ "reference_slip_rate", 1, pylith::topology::FieldBase::SCALAR },
	{ "cohesion", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Number of state variables.
      const int numStateVars = 1;

      // State variables.
      const pylith::materials::Metadata::ParamDescription stateVars[numStateVars] = {
	{ "slip_rate", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 3;
      const char* dbProperties[numDBProperties] = { 
	"static-coefficient",
	"reference-slip-rate",
	"cohesion",
      };

      // These are the state variables stored during the
      // simulation. Usually, we store only the time-dependent values
      // needed to compute the behavior at a given point in time. In
      // this example, however, for illustration purposes we store the
      // slip rate.

      const int numDBStateVars = 1;
      const char* dbStateVars[numDBStateVars] = { 
	"slip-rate",
      };      
      
    } // _ViscousFriction
  } // friction
} // contrib

// Indices of fault constitutive parameters.
const int contrib::friction::ViscousFriction::p_coefS = 0;
const int contrib::friction::ViscousFriction::p_v0 = 
  contrib::friction::ViscousFriction::p_coefS + 1;
const int contrib::friction::ViscousFriction::p_cohesion =
  contrib::friction::ViscousFriction::p_v0 + 1;

// Indices of database values (order must match dbProperties)
const int contrib::friction::ViscousFriction::db_coefS = 0;
const int contrib::friction::ViscousFriction::db_v0 =
  contrib::friction::ViscousFriction::db_coefS + 1;
const int contrib::friction::ViscousFriction::db_cohesion =
  contrib::friction::ViscousFriction::db_v0 + 1;

// Indices of state variables.
const int contrib::friction::ViscousFriction::s_slipRate = 0;

// Indices of database values (order must match dbProperties)
const int contrib::friction::ViscousFriction::db_slipRate = 0;

// ----------------------------------------------------------------------
// Default constructor.
contrib::friction::ViscousFriction::ViscousFriction(void) :
  pylith::friction::FrictionModel(pylith::materials::Metadata(_ViscousFriction::properties,
				    _ViscousFriction::numProperties,
				    _ViscousFriction::dbProperties,
				    _ViscousFriction::numDBProperties,
				    _ViscousFriction::stateVars,
				    _ViscousFriction::numStateVars,
				    _ViscousFriction::dbStateVars,
				    _ViscousFriction::numDBStateVars))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
contrib::friction::ViscousFriction::~ViscousFriction(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
contrib::friction::ViscousFriction::_dbToProperties(PylithScalar* const propValues,
						    const pylith::scalar_array& dbValues) const
{ // _dbToProperties
  // Check consistency of arguments
  assert(propValues);
  const int numDBValues = dbValues.size();
  assert(_ViscousFriction::numDBProperties == numDBValues);

  // Extract values from array using our defined indices.
  const PylithScalar coefS = dbValues[db_coefS];
  const PylithScalar v0 = dbValues[db_v0];
  const PylithScalar cohesion = dbValues[db_cohesion];

  // Check for reasonable values. If user supplied unreasonable values
  // throw an exception.
  if (coefS <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for static "
	<< "coefficient of friction.\n"
	<< "Static coefficient of friction: " << v0 << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (cohesion < 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned negative value for cohesion.\n"
	<< "Cohesion: " << cohesion << "\n";
    throw std::runtime_error(msg.str());
  } // if

  // Compute parameters that we store from the user-supplied parameters.
  propValues[p_coefS] = coefS;
  propValues[p_v0] = db_v0;
  propValues[p_cohesion] = cohesion;

} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
contrib::friction::ViscousFriction::_nondimProperties(PylithScalar* const values,
						      const int nvalues) const
{ // _nondimProperties
  // Check consistency of arguments.
  assert(_normalizer);
  assert(values);
  assert(nvalues == _ViscousFriction::numProperties);

  // Get scales needed to nondimensional parameters from the
  // Nondimensional object.
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar velocityScale = lengthScale / timeScale;

  // Use the Nondimensional::nondimensionalize() function to
  // nondimensionalize the quantities using the appropriate scale.
  values[p_v0] = _normalizer->nondimensionalize(values[p_v0], velocityScale);
  values[p_cohesion] = 
    _normalizer->nondimensionalize(values[p_cohesion], pressureScale);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
contrib::friction::ViscousFriction::_dimProperties(PylithScalar* const values,
						   const int nvalues) const
{ // _dimProperties
  // Check consistency of arguments.
  assert(_normalizer);
  assert(values);
  assert(nvalues == _ViscousFriction::numProperties);

  // Get scales needed to dimensional parameters from the
  // Nondimensional object.
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar velocityScale = lengthScale / timeScale;

  // Use the Nondimensional::dimensionalize() function to
  // dimensionalize the quantities using the appropriate scale.
  values[p_v0] = _normalizer->dimensionalize(values[p_v0], velocityScale);
  values[p_cohesion] = 
    _normalizer->dimensionalize(values[p_cohesion], pressureScale);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
contrib::friction::ViscousFriction::_dbToStateVars(PylithScalar* const stateValues,
						   const pylith::scalar_array& dbValues) const
{ // _dbToStateVars
  // Check consistency of arguments.
  assert(stateValues);
  const int numDBValues = dbValues.size();
  assert(_ViscousFriction::numDBStateVars == numDBValues);

  // Compute friction parameters that we store from the user-supplied
  // friction parameters.
  const PylithScalar slipRate = dbValues[db_slipRate];
 
  // Store computed friction parameters in the properties array.
  stateValues[s_slipRate] = slipRate;
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
contrib::friction::ViscousFriction::_nondimStateVars(PylithScalar* const values,
						     const int nvalues) const
{ // _nondimStateVars
  // Check consistency of arguments.
  assert(_normalizer);
  assert(values);
  assert(nvalues == _ViscousFriction::numStateVars);

  // Get scales needed to nondimensional parameters from the
  // Nondimensional object.
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar velocityScale = lengthScale / timeScale;

  // Use the Nondimensional::dimensionalize() function to
  // dimensionalize the quantities using the appropriate scale.
  values[s_slipRate] = 
    _normalizer->nondimensionalize(values[s_slipRate], velocityScale);
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
contrib::friction::ViscousFriction::_dimStateVars(PylithScalar* const values,
						  const int nvalues) const
{ // _dimStateVars
  // Check consistency of arguments.
  assert(_normalizer);
  assert(values);
  assert(nvalues == _ViscousFriction::numStateVars);

  // Get scales needed to dimensional parameters from the
  // Nondimensional object.
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar velocityScale = lengthScale / timeScale;

  // Use the Nondimensional::dimensionalize() function to
  // dimensionalize the quantities using the appropriate scale.
  values[s_slipRate] = 
    _normalizer->dimensionalize(values[s_slipRate], velocityScale);
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
PylithScalar
contrib::friction::ViscousFriction::_calcFriction(const PylithScalar t,
						  const PylithScalar slip,
						  const PylithScalar slipRate,
						  const PylithScalar normalTraction,
						  const PylithScalar* properties,
						  const int numProperties,
						  const PylithScalar* stateVars,
						  const int numStateVars)
{ // _calcFriction
  // Check consistency of arguments.
  assert(properties);
  assert(_ViscousFriction::numProperties == numProperties);
  assert(numStateVars);
  assert(_ViscousFriction::numStateVars == numStateVars);

  // Compute friction traction.
  PylithScalar friction = 0.0;
  PylithScalar mu_f = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression
    mu_f = properties[p_coefS] * (1.0 + fabs(slipRate) / properties[p_v0]);
    friction = - mu_f * normalTraction + properties[p_cohesion];
  } // if

  return friction;
} // _calcFriction

// ----------------------------------------------------------------------
// Compute derivative of friction with slip from properties and state variables.
PylithScalar
contrib::friction::ViscousFriction::_calcFrictionDeriv(const PylithScalar t,
						       const PylithScalar slip,
						       const PylithScalar slipRate,
						       const PylithScalar normalTraction,
						       const PylithScalar* properties,
						       const int numProperties,
						       const PylithScalar* stateVars,
						       const int numStateVars)
{ // _calcFrictionDeriv
  // Check consistency of arguments.
  assert(properties);
  assert(_ViscousFriction::numProperties == numProperties);
  assert(numStateVars);
  assert(_ViscousFriction::numStateVars == numStateVars);

  // Compute friction traction.
  PylithScalar frictionDeriv = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression

    // We want the derivative of friction with respect to
    // slip. Because the friction model depends on slip rate, we
    // approximate the derivative with respect to slip by taking the
    // derivative with respect to slip rate and multiplying by the time step (dt).
    frictionDeriv = -normalTraction * properties[p_coefS] / (properties[p_v0] * _dt);
  } // if

  return frictionDeriv;
} // _calcFrictionDeriv

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
contrib::friction::ViscousFriction::_updateStateVars(const PylithScalar t,
						     const PylithScalar slip,
						     const PylithScalar slipRate,
						     const PylithScalar normalTraction,
						     PylithScalar* const stateVars,
						     const int numStateVars,
						     const PylithScalar* properties,
						     const int numProperties)
{ // _updateStateVars
  // Check consistency of arguments.
  assert(properties);
  assert(_ViscousFriction::numProperties == numProperties);
  assert(numStateVars);
  assert(_ViscousFriction::numStateVars == numStateVars);

  // Store state variables.
  stateVars[s_slipRate] = stateVars[s_slipRate]; 
} // _updateStateVars


// End of file 
