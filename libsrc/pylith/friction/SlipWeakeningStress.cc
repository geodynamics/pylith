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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "SlipWeakeningStress.hh" // implementation of object methods

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
    namespace _SlipWeakeningStress {

      // Number of physical properties.
      const int numProperties = 4;

      // Physical properties.
      const pylith::materials::Metadata::ParamDescription properties[4] = {
	{ "static_stress", 1, pylith::topology::FieldBase::SCALAR },
	{ "dynamic_stress", 1, pylith::topology::FieldBase::SCALAR },
        { "slip_weakening_parameter", 1, pylith::topology::FieldBase::SCALAR },
	{ "weakening_time", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Number of State Variables.
      const int numStateVars = 2;

      // State Variables.
      const pylith::materials::Metadata::ParamDescription stateVars[] = {
	{ "cumulative_slip", 1, pylith::topology::FieldBase::SCALAR },
	{ "previous_slip", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 4;
      const char* dbProperties[4] = { "static-stress",
				      "dynamic-stress",
				      "slip-weakening-parameter",
				      "weakening-time",
      };

      const int numDBStateVars = 2;
      const char* dbStateVars[2] = { "cumulative-slip",
				     "previous-slip",
      };      
      
    } // _SlipWeakeningStress
  } // friction
} // pylith

// Indices of physical properties
const int pylith::friction::SlipWeakeningStress::p_stressStatic = 0;
const int pylith::friction::SlipWeakeningStress::p_stressDyn = 
  pylith::friction::SlipWeakeningStress::p_stressStatic + 1;
const int pylith::friction::SlipWeakeningStress::p_d0 = 
  pylith::friction::SlipWeakeningStress::p_stressDyn + 1;
const int pylith::friction::SlipWeakeningStress::p_weaktime =
  pylith::friction::SlipWeakeningStress::p_d0 + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::SlipWeakeningStress::db_stressStatic = 0;
const int pylith::friction::SlipWeakeningStress::db_stressDyn = 
  pylith::friction::SlipWeakeningStress::db_stressStatic + 1;
const int pylith::friction::SlipWeakeningStress::db_d0 = 
  pylith::friction::SlipWeakeningStress::db_stressDyn + 1;
const int pylith::friction::SlipWeakeningStress::db_weaktime =
  pylith::friction::SlipWeakeningStress::db_d0 + 1;

// Indices of state variables.
const int pylith::friction::SlipWeakeningStress::s_slipCum = 0;
const int pylith::friction::SlipWeakeningStress::s_slipPrev = 
  pylith::friction::SlipWeakeningStress::s_slipCum + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::SlipWeakeningStress::db_slipCum = 0;
const int pylith::friction::SlipWeakeningStress::db_slipPrev = 
  pylith::friction::SlipWeakeningStress::db_slipCum + 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::SlipWeakeningStress::SlipWeakeningStress(void) :
  FrictionModel(materials::Metadata(_SlipWeakeningStress::properties,
				    _SlipWeakeningStress::numProperties,
				    _SlipWeakeningStress::dbProperties,
				    _SlipWeakeningStress::numDBProperties,
				    _SlipWeakeningStress::stateVars,
				    _SlipWeakeningStress::numStateVars,
				    _SlipWeakeningStress::dbStateVars,
				    _SlipWeakeningStress::numDBStateVars))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::SlipWeakeningStress::~SlipWeakeningStress(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::friction::SlipWeakeningStress::_dbToProperties(PylithScalar* const propValues,
						       const scalar_array& dbValues) const
{ // _dbToProperties
  assert(propValues);
  const int numDBValues = dbValues.size();
  assert(_SlipWeakeningStress::numDBProperties == numDBValues);

  const PylithScalar db_static = dbValues[db_stressStatic];
  const PylithScalar db_dyn = dbValues[db_stressDyn];
  const PylithScalar db_do = dbValues[db_d0];
  const PylithScalar db_t = dbValues[db_weaktime];

  if (db_static < 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned negative value for static stress.\n"
	<< "static stress: " << db_static << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (db_dyn < 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned negative value for dynamic stress.\n"
	<< "dynamic stress: " << db_dyn << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (db_d0 <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for slip weakening parameter.\n"
	<< "slip weakening parameter: " << db_d0 << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (db_t < 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned negative value for weakening time.\n"
	<< "weakening time: " << db_t << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[p_stressStatic] = db_static;
  propValues[p_stressDyn] = db_dyn;
  propValues[p_d0] = db_do;
  propValues[p_weaktime] = db_t;

} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::SlipWeakeningStress::_nondimProperties(PylithScalar* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakeningStress::numProperties);

  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();

  values[p_d0] /= lengthScale;
  values[p_stressStatic] /= pressureScale;
  values[p_stressDyn] /= pressureScale;
  values[p_weaktime] /= timeScale;
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::SlipWeakeningStress::_dimProperties(PylithScalar* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakeningStress::numProperties);

  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();

  values[p_d0] *= lengthScale;
  values[p_stressStatic] *= pressureScale;
  values[p_stressDyn] *= pressureScale;
  values[p_weaktime] *= timeScale;
} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
pylith::friction::SlipWeakeningStress::_dbToStateVars(PylithScalar* const stateValues,
						      const scalar_array& dbValues) const
{ // _dbToStateVars
  assert(stateValues);
  const int numDBValues = dbValues.size();
  assert(_SlipWeakeningStress::numDBStateVars == numDBValues);

  const PylithScalar cumulativeSlip = dbValues[db_slipCum];
  const PylithScalar previousSlip = dbValues[db_slipPrev];
 
  stateValues[s_slipCum] = cumulativeSlip;
  stateValues[s_slipPrev] = previousSlip;
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::friction::SlipWeakeningStress::_nondimStateVars(PylithScalar* const values,
							const int nvalues) const
{ // _nondimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakeningStress::numStateVars);

  const PylithScalar lengthScale = _normalizer->lengthScale();

  values[s_slipCum] /= lengthScale;
  values[s_slipPrev] /= lengthScale;
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::friction::SlipWeakeningStress::_dimStateVars(PylithScalar* const values,
					       const int nvalues) const
{ // _dimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakeningStress::numStateVars);

  const PylithScalar lengthScale = _normalizer->lengthScale();

  values[s_slipCum] *= lengthScale;
  values[s_slipPrev] *= lengthScale;
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
PylithScalar
pylith::friction::SlipWeakeningStress::_calcFriction(const PylithScalar t,
						     const PylithScalar slip,
						     const PylithScalar slipRate,
						     const PylithScalar normalTraction,
						     const PylithScalar* properties,
						     const int numProperties,
						     const PylithScalar* stateVars,
						     const int numStateVars)
{ // _calcFriction
  assert(properties);
  assert(_SlipWeakeningStress::numProperties == numProperties);
  assert(stateVars);
  assert(_SlipWeakeningStress::numStateVars == numStateVars);

  PylithScalar friction = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression
    if (stateVars[s_slipCum] < properties[p_d0] &&
	t < properties[p_weaktime]) {
	// if/else linear slip-weakening form
	friction = properties[p_stressStatic] -
	  (properties[p_stressStatic] - properties[p_stressDyn]) * 
	  stateVars[s_slipCum] / properties[p_d0];
      } else {
	friction = properties[p_stressDyn];
      } // if/else
  } // if

  PetscLogFlops(6);

  return friction;
} // _calcFriction

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::SlipWeakeningStress::_updateStateVars(const PylithScalar t,
							const PylithScalar slip,
							const PylithScalar slipRate,
							const PylithScalar normalTraction,
							PylithScalar* const stateVars,
							const int numStateVars,
							const PylithScalar* properties,
							const int numProperties)
{ // _updateStateVars
  assert(properties);
  assert(_SlipWeakeningStress::numProperties == numProperties);
  assert(stateVars);
  assert(_SlipWeakeningStress::numStateVars == numStateVars);

  const PylithScalar tolerance = 1.0e-12;
  if (slipRate > tolerance) {
    const PylithScalar slipPrev = stateVars[s_slipPrev];

    stateVars[s_slipPrev] = slip;
    stateVars[s_slipCum] += fabs(slip - slipPrev);
  } else {
    // Sliding has stopped, so reset state variables.
    stateVars[s_slipPrev] = slip;
    stateVars[s_slipCum] = 0.0;
  } // else
} // _updateStateVars


// End of file 
