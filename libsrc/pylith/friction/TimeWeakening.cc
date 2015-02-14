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

#include "TimeWeakening.hh" // implementation of object methods

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
    namespace _TimeWeakening {

      // Number of physical properties.
      const int numProperties = 4;

      // Physical properties.
      const pylith::materials::Metadata::ParamDescription properties[] = {
	{ "static_coefficient", 1, pylith::topology::FieldBase::SCALAR },
	{ "dynamic_coefficient", 1, pylith::topology::FieldBase::SCALAR },
        { "time_weakening_parameter", 1, pylith::topology::FieldBase::SCALAR },
	{ "cohesion", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Number of State Variables.
      const int numStateVars = 1;

      // State Variables.
      const pylith::materials::Metadata::ParamDescription stateVars[] = {
	{ "elapsed_time", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 4;
      const char* dbProperties[4] = { "static-coefficient",
				      "dynamic-coefficient",
				      "time-weakening-parameter",
				      "cohesion",
      };

      const int numDBStateVars = 1;
      const char* dbStateVars[1] = { "elapsed-time",
      };      
      
    } // _TimeWeakening
  } // friction
} // pylith

// Indices of physical properties
const int pylith::friction::TimeWeakening::p_coefS = 0;
const int pylith::friction::TimeWeakening::p_coefD = 
  pylith::friction::TimeWeakening::p_coefS + 1;
const int pylith::friction::TimeWeakening::p_Tc = 
  pylith::friction::TimeWeakening::p_coefD + 1;
const int pylith::friction::TimeWeakening::p_cohesion =
  pylith::friction::TimeWeakening::p_Tc + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::TimeWeakening::db_coefS = 0;
const int pylith::friction::TimeWeakening::db_coefD = 
  pylith::friction::TimeWeakening::db_coefS + 1;
const int pylith::friction::TimeWeakening::db_Tc = 
  pylith::friction::TimeWeakening::db_coefD + 1;
const int pylith::friction::TimeWeakening::db_cohesion =
  pylith::friction::TimeWeakening::db_Tc + 1;

// Indices of state variables.
const int pylith::friction::TimeWeakening::s_time = 0;

// Indices of database values (order must match dbProperties)
const int pylith::friction::TimeWeakening::db_time = 0;

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::TimeWeakening::TimeWeakening(void) :
  FrictionModel(materials::Metadata(_TimeWeakening::properties,
				    _TimeWeakening::numProperties,
				    _TimeWeakening::dbProperties,
				    _TimeWeakening::numDBProperties,
				    _TimeWeakening::stateVars,
				    _TimeWeakening::numStateVars,
				    _TimeWeakening::dbStateVars,
				    _TimeWeakening::numDBStateVars))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::TimeWeakening::~TimeWeakening(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::friction::TimeWeakening::_dbToProperties(
					   PylithScalar* const propValues,
					   const scalar_array& dbValues) const
{ // _dbToProperties
  assert(propValues);
  const int numDBValues = dbValues.size();
  assert(_TimeWeakening::numDBProperties == numDBValues);

  const PylithScalar db_static = dbValues[db_coefS];
  const PylithScalar db_dynamic = dbValues[db_coefD];
  const PylithScalar db_To = dbValues[db_Tc];
  const PylithScalar db_c = dbValues[db_cohesion];

  if (db_static < 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned negative value for static coefficient "
	<< "of friction.\n"
	<< "static coefficient of friction: " << db_static << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (db_dynamic < 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned negative value for dynamic coefficient "
	<< "of friction.\n"
	<< "dynamic coefficient of friction: " << db_dynamic << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (db_To <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for time weakening parameter "
	<< "of friction.\n"
	<< "time weakening parameter of friction: " << db_To << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[p_coefS] = db_static;
  propValues[p_coefD] = db_dynamic;
  propValues[p_Tc] = db_To;
  propValues[p_cohesion] = db_c;

} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::TimeWeakening::_nondimProperties(PylithScalar* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _TimeWeakening::numProperties);

  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_Tc] /= timeScale;
  values[p_cohesion] /= pressureScale;
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::TimeWeakening::_dimProperties(PylithScalar* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _TimeWeakening::numProperties);

  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_Tc] *= timeScale;
  values[p_cohesion] *= pressureScale;
} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
pylith::friction::TimeWeakening::_dbToStateVars(
					   PylithScalar* const stateValues,
					   const scalar_array& dbValues) const
{ // _dbToStateVars
  assert(stateValues);
  const int numDBValues = dbValues.size();
  assert(_TimeWeakening::numDBStateVars == numDBValues);

  const PylithScalar timeElasped = dbValues[db_time];
  
  stateValues[s_time] = timeElasped;
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::friction::TimeWeakening::_nondimStateVars(PylithScalar* const values,
						  const int nvalues) const
{ // _nondimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _TimeWeakening::numStateVars);

  const PylithScalar timeScale = _normalizer->timeScale();

  values[s_time] /= timeScale;
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::friction::TimeWeakening::_dimStateVars(PylithScalar* const values,
					       const int nvalues) const
{ // _dimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _TimeWeakening::numStateVars);

  const PylithScalar timeScale = _normalizer->timeScale();

  values[s_time] *= timeScale;
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
PylithScalar
pylith::friction::TimeWeakening::_calcFriction(const PylithScalar t,
					       const PylithScalar slip,
					       const PylithScalar slipRate,
					       const PylithScalar normalTraction,
					       const PylithScalar* properties,
					       const int numProperties,
					       const PylithScalar* stateVars,
					       const int numStateVars)
{ // _calcFriction
  assert(properties);
  assert(_TimeWeakening::numProperties == numProperties);
  assert(numStateVars);
  assert(_TimeWeakening::numStateVars == numStateVars);

  PylithScalar friction = 0.0;
  PylithScalar mu_f = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression
    if (stateVars[s_time] < properties[p_Tc]) {
	// if/else linear time-weakening form of mu_f 
	mu_f = properties[p_coefS] -
	  (properties[p_coefS] - properties[p_coefD]) * 
	  stateVars[s_time] / properties[p_Tc];
      } else {
	mu_f = properties[p_coefD];
      } // if/else
    friction = - mu_f * normalTraction + properties[p_cohesion];
  } else {
    friction = properties[p_cohesion];
  } // if/else

  PetscLogFlops(6);

  return friction;
} // _calcFriction

// ----------------------------------------------------------------------
// Compute derivative of friction with slip from properties and
// state variables.
PylithScalar
pylith::friction::TimeWeakening::_calcFrictionDeriv(const PylithScalar t,
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


// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::TimeWeakening::_updateStateVars(const PylithScalar t,
						  const PylithScalar slip,
						  const PylithScalar slipRate,
						  const PylithScalar normalTraction,
						  PylithScalar* const stateVars,
						  const int numStateVars,
						  const PylithScalar* properties,
						  const int numProperties)
{ // _updateStateVars
  assert(properties);
  assert(_TimeWeakening::numProperties == numProperties);
  assert(numStateVars);
  assert(_TimeWeakening::numStateVars == numStateVars);

  const PylithScalar tolerance = 1.0e-12;
  if (slipRate > tolerance) {
    const PylithScalar dt = _dt;

    stateVars[s_time] += dt;
  } else {
    stateVars[s_time] = 0.0;
  } // else

} // _updateStateVars


// End of file 
