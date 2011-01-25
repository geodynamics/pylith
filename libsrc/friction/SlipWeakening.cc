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

#include "SlipWeakening.hh" // implementation of object methods

#include "pylith/materials/Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
#include <iostream>
// ----------------------------------------------------------------------
namespace pylith {
  namespace friction {
    namespace _SlipWeakening {

      // Number of physical properties.
      const int numProperties = 4;

      // Physical properties.
      const pylith::materials::Metadata::ParamDescription properties[] = {
	{ "static_coefficient", 1, pylith::topology::FieldBase::SCALAR },
	{ "dynamic_coefficient", 1, pylith::topology::FieldBase::SCALAR },
        { "slip_weakening_parameter", 1, pylith::topology::FieldBase::SCALAR },
	{ "cohesion", 1, pylith::topology::FieldBase::SCALAR },
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
      const char* dbProperties[4] = { "static-coefficient",
				      "dynamic-coefficient",
				      "slip-weakening-parameter",
				      "cohesion",
      };

      const int numDBStateVars = 2;
      const char* dbStateVars[2] = { "cumulative-slip",
				     "previous-slip",
      };      
      
    } // _SlipWeakening
  } // friction
} // pylith

// Indices of physical properties
const int pylith::friction::SlipWeakening::p_coefS = 0;
const int pylith::friction::SlipWeakening::p_coefD = 
  pylith::friction::SlipWeakening::p_coefS + 1;
const int pylith::friction::SlipWeakening::p_d0 = 
  pylith::friction::SlipWeakening::p_coefD + 1;
const int pylith::friction::SlipWeakening::p_cohesion =
  pylith::friction::SlipWeakening::p_d0 + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::SlipWeakening::db_coefS = 0;
const int pylith::friction::SlipWeakening::db_coefD = 
  pylith::friction::SlipWeakening::db_coefS + 1;
const int pylith::friction::SlipWeakening::db_d0 = 
  pylith::friction::SlipWeakening::db_coefD + 1;
const int pylith::friction::SlipWeakening::db_cohesion =
  pylith::friction::SlipWeakening::db_d0 + 1;

// Indices of state variables.
const int pylith::friction::SlipWeakening::s_slipCum = 0;
const int pylith::friction::SlipWeakening::s_slipPrev = 
  pylith::friction::SlipWeakening::s_slipCum + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::SlipWeakening::db_slipCum = 0;
const int pylith::friction::SlipWeakening::db_slipPrev = 
  pylith::friction::SlipWeakening::db_slipCum + 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::SlipWeakening::SlipWeakening(void) :
  FrictionModel(materials::Metadata(_SlipWeakening::properties,
				    _SlipWeakening::numProperties,
				    _SlipWeakening::dbProperties,
				    _SlipWeakening::numDBProperties,
				    _SlipWeakening::stateVars,
				    _SlipWeakening::numStateVars,
				    _SlipWeakening::dbStateVars,
				    _SlipWeakening::numDBStateVars))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::SlipWeakening::~SlipWeakening(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::friction::SlipWeakening::_dbToProperties(
					   double* const propValues,
					   const double_array& dbValues) const
{ // _dbToProperties
  assert(propValues);
  const int numDBValues = dbValues.size();
  assert(_SlipWeakening::numDBProperties == numDBValues);

  const double db_static = dbValues[db_coefS];
  const double db_dynamic = dbValues[db_coefD];
  const double db_do = dbValues[db_d0];
  const double db_c = dbValues[db_cohesion];

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

  if (db_d0 <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for slip weakening parameter "
	<< "of friction.\n"
	<< "slip weakening parameter of friction: " << db_d0 << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[p_coefS] = db_static;
  propValues[p_coefD] = db_dynamic;
  propValues[p_d0] = db_do;
  propValues[p_cohesion] = db_c;

} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::SlipWeakening::_nondimProperties(double* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakening::numProperties);

  const double lengthScale = _normalizer->lengthScale();
  const double pressureScale = _normalizer->pressureScale();

  values[p_d0] /= lengthScale;
  values[p_cohesion] /= pressureScale;
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::SlipWeakening::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakening::numProperties);

  const double lengthScale = _normalizer->lengthScale();
  const double pressureScale = _normalizer->pressureScale();

  values[p_d0] *= lengthScale;
  values[p_cohesion] *= pressureScale;
} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
pylith::friction::SlipWeakening::_dbToStateVars(
					   double* const stateValues,
					   const double_array& dbValues) const
{ // _dbToStateVars
  assert(stateValues);
  const int numDBValues = dbValues.size();
  assert(_SlipWeakening::numDBStateVars == numDBValues);

  const double cumulativeSlip = dbValues[db_slipCum];
  const double previousSlip = dbValues[db_slipPrev];
 
  stateValues[s_slipCum] = cumulativeSlip;
  stateValues[s_slipPrev] = previousSlip;
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::friction::SlipWeakening::_nondimStateVars(double* const values,
						  const int nvalues) const
{ // _nondimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakening::numStateVars);

  const double lengthScale = _normalizer->lengthScale();

  values[s_slipCum] /= lengthScale;
  values[s_slipPrev] /= lengthScale;
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::friction::SlipWeakening::_dimStateVars(double* const values,
					       const int nvalues) const
{ // _dimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _SlipWeakening::numStateVars);

  const double lengthScale = _normalizer->lengthScale();

  values[s_slipCum] *= lengthScale;
  values[s_slipPrev] *= lengthScale;
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
double
pylith::friction::SlipWeakening::_calcFriction(const double slip,
						const double slipRate,
						const double normalTraction,
						const double* properties,
						const int numProperties,
						const double* stateVars,
						const int numStateVars)
{ // _calcFriction
  assert(properties);
  assert(_SlipWeakening::numProperties == numProperties);
  assert(stateVars);
  assert(_SlipWeakening::numStateVars == numStateVars);

  double friction = 0.0;
  double mu_f = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression
    if (stateVars[s_slipCum] < properties[p_d0]) {
	// if/else linear slip-weakening form of mu_f 
	mu_f = properties[p_coefS] -
	  (properties[p_coefS] - properties[p_coefD]) * 
	  stateVars[s_slipCum] / properties[p_d0];
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
// Update state variables (for next time step).
void
pylith::friction::SlipWeakening::_updateStateVars(const double slip,
						  const double slipRate,
						  const double normalTraction,
						  double* const stateVars,
						  const int numStateVars,
						  const double* properties,
						  const int numProperties)
{ // _updateStateVars
  assert(properties);
  assert(_SlipWeakening::numProperties == numProperties);
  assert(stateVars);
  assert(_SlipWeakening::numStateVars == numStateVars);

  // Don't reset state variables when sliding stops
  // (SCEC dynamic rupture branch only).
  const double slipPrev = stateVars[s_slipPrev];
  stateVars[s_slipPrev] = stateVars[s_slipCum];
  stateVars[s_slipCum] += fabs(slip - slipPrev);
} // _updateStateVars


// End of file 
