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

#include "RateStateAgeing.hh" // implementation of object methods

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
    namespace _RateStateAgeing {

      // Number of physical properties.
      const int numProperties = 6;

      // Physical properties.
      const pylith::materials::Metadata::ParamDescription properties[] = {
        { "reference_friction_coefficient", 1, pylith::topology::FieldBase::SCALAR },
        { "reference_slip_rate", 1, pylith::topology::FieldBase::SCALAR },
        { "characteristic_slip_distance", 1, pylith::topology::FieldBase::SCALAR },
        { "constitutive_parameter_a", 1, pylith::topology::FieldBase::SCALAR },
        { "constitutive_parameter_b", 1, pylith::topology::FieldBase::SCALAR },
	{ "cohesion", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Number of State Variables.
      const int numStateVars = 1;

      // State Variables.
      const pylith::materials::Metadata::ParamDescription stateVars[] = {
        { "state_variable", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 6;
      const char* dbProperties[6] = {
	"reference-friction-coefficient",
	"reference-slip-rate",
	"characteristic-slip-distance",
	"constitutive-parameter-a",
	"constitutive-parameter-b",
	"cohesion",
      };

      const int numDBStateVars = 1;
      const char* dbStateVars[1] = {
	"state-variable",
      };
      
    } // _RateStateAgeing
  } // friction
} // pylith

// Indices of physical properties
const int pylith::friction::RateStateAgeing::p_coef = 0;
const int pylith::friction::RateStateAgeing::p_slipRate0 = 
  pylith::friction::RateStateAgeing::p_coef + 1;
const int pylith::friction::RateStateAgeing::p_L = 
  pylith::friction::RateStateAgeing::p_slipRate0 + 1;
const int pylith::friction::RateStateAgeing::p_a = 
  pylith::friction::RateStateAgeing::p_L + 1;
const int pylith::friction::RateStateAgeing::p_b = 
  pylith::friction::RateStateAgeing::p_a + 1;
const int pylith::friction::RateStateAgeing::p_cohesion =
  pylith::friction::RateStateAgeing::p_b + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::RateStateAgeing::db_coef = 0;
const int pylith::friction::RateStateAgeing::db_slipRate0 = 
  pylith::friction::RateStateAgeing::db_coef + 1;
const int pylith::friction::RateStateAgeing::db_L = 
  pylith::friction::RateStateAgeing::db_slipRate0 + 1;
const int pylith::friction::RateStateAgeing::db_a = 
  pylith::friction::RateStateAgeing::db_L + 1;
const int pylith::friction::RateStateAgeing::db_b = 
  pylith::friction::RateStateAgeing::db_a + 1;
const int pylith::friction::RateStateAgeing::db_cohesion =
  pylith::friction::RateStateAgeing::db_b + 1;

// Indices of state variables.
const int pylith::friction::RateStateAgeing::s_state = 0;

// Indices of database values (order must match dbProperties)
const int pylith::friction::RateStateAgeing::db_state = 0;

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::RateStateAgeing::RateStateAgeing(void) :
  FrictionModel(materials::Metadata(_RateStateAgeing::properties,
				    _RateStateAgeing::numProperties,
				    _RateStateAgeing::dbProperties,
				    _RateStateAgeing::numDBProperties,
				    _RateStateAgeing::stateVars,
				    _RateStateAgeing::numStateVars,
				    _RateStateAgeing::dbStateVars,
				    _RateStateAgeing::numDBStateVars)),
  _linearSlipRate(1.0e-12)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::RateStateAgeing::~RateStateAgeing(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set nondimensional slip rate below which friction varies
//  linearly with slip rate.
void
pylith::friction::RateStateAgeing::linearSlipRate(const PylithScalar value)
{ // linearSlipRate
  if (value < 0.0) {
    std::ostringstream msg;
    msg << "Nondimensional linear slip rate threshold (" << value << ") for rate state friction model "
	<< label() << " must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if
  
  _linearSlipRate = value;
} // linearSlipRate

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::friction::RateStateAgeing::_dbToProperties(PylithScalar* const propValues,
						   const scalar_array& dbValues) const
{ // _dbToProperties
  assert(propValues);
  const int numDBValues = dbValues.size();
  assert(_RateStateAgeing::numDBProperties == numDBValues);

  const PylithScalar frictionCoef = dbValues[db_coef];
  const PylithScalar slipRate0 = dbValues[db_slipRate0];
  const PylithScalar dc = dbValues[db_L];
  const PylithScalar a = dbValues[db_a];
  const PylithScalar b = dbValues[db_b];
  const PylithScalar cohesion = dbValues[db_cohesion];
 
  if (frictionCoef < 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned negative value for reference coefficient "
	<< "of Rate and State friction Ageing Law.\n"
	<< "reference coefficient of friction: " << frictionCoef << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (dc <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for characteristic"
	<< "slip distance of Rate and State friction Ageing Law.\n"
	<< "characteristic slip distance: " << dc << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (a <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
	<< "parameter 'a' of Rate and State friction Ageing Law.\n"
	<< "Rate and State parameter 'a' of Ageing Law of friction: " << a << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (b <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
	<< "parameter 'b' of Rate and State friction Ageing Law.\n"
	<< "Rate and State parameter 'b' of Ageing Law of friction: " << b << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[p_coef] = frictionCoef;
  propValues[p_slipRate0] = slipRate0;
  propValues[p_L] = dc;
  propValues[p_a] = a;
  propValues[p_b] = b;
  propValues[p_cohesion] = cohesion;

} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::RateStateAgeing::_nondimProperties(PylithScalar* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeing::numProperties);

  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_slipRate0] /= lengthScale / timeScale;
  values[p_L] /= lengthScale;
  values[p_cohesion] /= pressureScale;
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::RateStateAgeing::_dimProperties(PylithScalar* const values,
						  const int nvalues) const
{ // _dimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeing::numProperties);

  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_slipRate0] *= lengthScale / timeScale;
  values[p_L] *= lengthScale;
  values[p_cohesion] *= pressureScale;
} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
pylith::friction::RateStateAgeing::_dbToStateVars(PylithScalar* const stateValues,
						  const scalar_array& dbValues) const
{ // _dbToStateVars
  assert(stateValues);
  const int numDBValues = dbValues.size();
  assert(_RateStateAgeing::numDBStateVars == numDBValues);

  stateValues[s_state] = dbValues[db_state];
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::friction::RateStateAgeing::_nondimStateVars(PylithScalar* const values,
						    const int nvalues) const
{ // _nondimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeing::numStateVars);

  const PylithScalar timeScale = _normalizer->timeScale();

  values[s_state] /= timeScale;
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::friction::RateStateAgeing::_dimStateVars(PylithScalar* const values,
						 const int nvalues) const
{ // _dimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeing::numStateVars);

  const PylithScalar timeScale = _normalizer->timeScale();

  values[s_state] *= timeScale;
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
PylithScalar
pylith::friction::RateStateAgeing::_calcFriction(const PylithScalar t,
						 const PylithScalar slip,
						 const PylithScalar slipRate,
						 const PylithScalar normalTraction,
						 const PylithScalar* properties,
						 const int numProperties,
						 const PylithScalar* stateVars,
						 const int numStateVars)
{ // _calcFriction
  assert(properties);
  assert(_RateStateAgeing::numProperties == numProperties);
  assert(numStateVars);
  assert(_RateStateAgeing::numStateVars == numStateVars);

  PylithScalar friction = 0.0;
  PylithScalar mu_f = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression

    const PylithScalar slipRateLinear = _linearSlipRate;

    const PylithScalar f0 = properties[p_coef];
    const PylithScalar a = properties[p_a];
    const PylithScalar b = properties[p_b];
    const PylithScalar L = properties[p_L];
    const PylithScalar slipRate0 = properties[p_slipRate0];

    // Prevent zero value for theta, reasonable value is L / slipRate0
    const PylithScalar theta = (stateVars[s_state] > 0.0) ? stateVars[s_state] : L / slipRate0;

    if (slipRate >= slipRateLinear) {
      mu_f = f0 + a*log(slipRate / slipRate0) + b*log(slipRate0*theta/L);
    } else {
      mu_f = f0 + a*log(slipRateLinear / slipRate0) + b*log(slipRate0*theta/L) -
	a*(1.0 - slipRate/slipRateLinear);
    } // else
    friction = -mu_f * normalTraction + properties[p_cohesion];
    
  } else {
    friction = properties[p_cohesion];
  } // if/else

  PetscLogFlops(12);

  return friction;
} // _calcFriction


// ----------------------------------------------------------------------
// Compute derivative of friction with slip from properties and
// state variables.
PylithScalar
pylith::friction::RateStateAgeing::_calcFrictionDeriv(const PylithScalar t,
						      const PylithScalar slip,
						      const PylithScalar slipRate,
						      const PylithScalar normalTraction,
						      const PylithScalar* properties,
						      const int numProperties,
						      const PylithScalar* stateVars,
						      const int numStateVars)
{ // _calcFrictionDeriv
  assert(properties);
  assert(_RateStateAgeing::numProperties == numProperties);
  assert(numStateVars);
  assert(_RateStateAgeing::numStateVars == numStateVars);

  PylithScalar frictionDeriv = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression

    const PylithScalar slipRateLinear = _linearSlipRate;

    const PylithScalar a = properties[p_a];

    if (slipRate >= slipRateLinear) {
      frictionDeriv = -normalTraction * a / (slipRate * _dt);
    } else {
      frictionDeriv = -normalTraction * a / (slipRateLinear * _dt);
    } // else
  } // if    

  PetscLogFlops(12);

  return frictionDeriv;
} // _calcFrictionDeriv


// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::RateStateAgeing::_updateStateVars(const PylithScalar t,
						    const PylithScalar slip,
						    const PylithScalar slipRate,
						    const PylithScalar normalTraction,
						    PylithScalar* const stateVars,
						    const int numStateVars,
						    const PylithScalar* properties,
						    const int numProperties)
{ // _updateStateVars
  assert(properties);
  assert(_RateStateAgeing::numProperties == numProperties);
  assert(numStateVars);
  assert(_RateStateAgeing::numStateVars == numStateVars);

  // d(theta)/dt = (1 - slipRate * theta / L)
  //
  // Use separation of variables to integrate above ODE from t->t+dt,
  // keeping slip rate constant.
  //
  // thetaTpdt = thetaT * exp(-slipRate/L * dt)
  //             + L/slipRate * (1 -  exp(-slipRate/L * dt))
  //
  // As slipRate --> 0, L/sliprate --> infinity and
  // exp(-sliprate/L*dt) --> 1.  To determine, d(theta)/dt near
  // sliprate == 0, we expand the exponential term in a Taylor series:
  //
  // exp(-x) = 1 - x +1/2*x**2 + 1/6*x**3
  //
  // This leads to (in the vicinity of slipRate == 0):
  //
  // thetaTpdt = thetaT * exp(-slipRate/L * dt)
  //             + dt - 0.5*(sliprate/L)*dt**2 + 1.0/6.0*(slipRate/L)*dt**3;

  const PylithScalar dt = _dt;
  const PylithScalar thetaTVertex = stateVars[s_state];
  const PylithScalar L = properties[p_L];
  const PylithScalar vDtL = slipRate * dt / L;
  const PylithScalar expTerm = exp(-vDtL);

  PylithScalar thetaTpdtVertex = 0.0;
  if (vDtL > 1.0e-20) {
    thetaTpdtVertex = thetaTVertex * expTerm + L / slipRate * (1 - expTerm);
    PetscLogFlops(7);
  } else {
    thetaTpdtVertex = thetaTVertex * expTerm + dt - 0.5 * slipRate/L * dt*dt;
    PetscLogFlops(9);
  } // if/else
  
  stateVars[s_state] = thetaTpdtVertex;

} // _updateStateVars


// End of file 
