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

#include "RateStateAgeing.hh" // implementation of object methods

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
  _minSlipRate(1.0e-12)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::RateStateAgeing::~RateStateAgeing(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set floor for slip rate used in computing friction. Used to
void
pylith::friction::RateStateAgeing::minSlipRate(const double value)
{ // minSlipRate
  if (value < 0.0) {
    std::ostringstream msg;
    msg << "Minimum slip rate (" << value << ") for rate state friction model "
	<< label() << " must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  _minSlipRate = value;
} // minSlipRate

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::friction::RateStateAgeing::_dbToProperties(
					   double* const propValues,
					   const double_array& dbValues) const
{ // _dbToProperties
  assert(propValues);
  const int numDBValues = dbValues.size();
  assert(_RateStateAgeing::numDBProperties == numDBValues);

  const double frictionCoef = dbValues[db_coef];
  const double slipRate0 = dbValues[db_slipRate0];
  const double dc = dbValues[db_L];
  const double a = dbValues[db_a];
  const double b = dbValues[db_b];
  const double cohesion = dbValues[db_cohesion];
 
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
pylith::friction::RateStateAgeing::_nondimProperties(double* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeing::numProperties);

  const double lengthScale = _normalizer->lengthScale();
  const double timeScale = _normalizer->timeScale();
  const double pressureScale = _normalizer->pressureScale();

  values[p_slipRate0] /= lengthScale / timeScale;
  values[p_L] /= lengthScale;
  values[p_cohesion] /= pressureScale;
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::RateStateAgeing::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeing::numProperties);

  const double lengthScale = _normalizer->lengthScale();
  const double timeScale = _normalizer->timeScale();
  const double pressureScale = _normalizer->pressureScale();

  values[p_slipRate0] *= lengthScale / timeScale;
  values[p_L] *= lengthScale;
  values[p_cohesion] *= pressureScale;
} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
pylith::friction::RateStateAgeing::_dbToStateVars(
					   double* const stateValues,
					   const double_array& dbValues) const
{ // _dbToStateVars
  assert(stateValues);
  const int numDBValues = dbValues.size();
  assert(_RateStateAgeing::numDBStateVars == numDBValues);

  stateValues[s_state] = dbValues[db_state];
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::friction::RateStateAgeing::_nondimStateVars(double* const values,
						    const int nvalues) const
{ // _nondimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeing::numStateVars);

  const double timeScale = _normalizer->timeScale();

  values[s_state] /= timeScale;
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::friction::RateStateAgeing::_dimStateVars(double* const values,
						      const int nvalues) const
{ // _dimStateVars
  assert(_normalizer);
  assert(values);
  assert(nvalues == _RateStateAgeing::numStateVars);

  const double timeScale = _normalizer->timeScale();

  values[s_state] *= timeScale;
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
double
pylith::friction::RateStateAgeing::_calcFriction(const double slip,
						const double slipRate,
						const double normalTraction,
						const double* properties,
						const int numProperties,
						const double* stateVars,
						const int numStateVars)
{ // _calcFriction
  assert(properties);
  assert(_RateStateAgeing::numProperties == numProperties);
  assert(numStateVars);
  assert(_RateStateAgeing::numStateVars == numStateVars);

  double friction = 0.0;
  double mu_f = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression

#if 0
    // regularized rate and state equation
    const double f0 = properties[p_coef];

    // Since regulatized friction -> 0 as slipRate -> 0, limit slip
    // rate to some minimum value
    const double slipRateEff = std::max(_minSlipRate, slipRate);

    const double slipRate0 = properties[p_slipRate0];
    const double a = properties[p_a];

    const double theta = stateVars[s_state];

    const double L = properties[p_L];
    const double b = properties[p_b];
    const double bLnTerm = b * log(slipRate0 * theta / L);

    const double expTerm = exp((f0 + bLnTerm)/a);
    const double sinhArg = 0.5 * slipRateEff / slipRate0 * expTerm;

    mu_f = a * asinh(sinhArg);
    friction = -mu_f * normalTraction + properties[p_cohesion];

#else

    const double slipRateLinear = _minSlipRate;

    const double f0 = properties[p_coef];
    const double a = properties[p_a];
    const double b = properties[p_b];
    const double L = properties[p_L];
    const double slipRate0 = properties[p_slipRate0];
    const double theta = stateVars[s_state];

    if (slipRate >= slipRateLinear) {
      mu_f = f0 + a*log(slipRate / slipRate0) + b*log(slipRate0*theta/L);
    } else {
      mu_f = f0 + a*log(slipRateLinear / slipRate0) + b*log(slipRate0*theta/L) -
	a*(1.0 - slipRate/slipRateLinear);
    } // else

    friction = -mu_f * normalTraction + properties[p_cohesion];

#if 0
    std::cout << "slip: " << slip
	      << ", slipRate: " << slipRate
	      << ", stateVar: " << theta
	      << ", mu_f: " << mu_f
	      << std::endl;
#endif

#endif
  } // if

  PetscLogFlops(12);

  return friction;
} // _calcFriction

// ----------------------------------------------------------------------
// Compute change in friction for a change in slip (Jacobian).
double
pylith::friction::RateStateAgeing::_calcFrictionSlope(const double slip,
						      const double slipRate,
						      const double normalTraction,
						      const double* properties,
						      const int numProperties,
						      const double* stateVars,
						      const int numStateVars)
{ // _calcFrictionSlope
  assert(properties);
  assert(_RateStateAgeing::numProperties == numProperties);
  assert(numStateVars);
  assert(_RateStateAgeing::numStateVars == numStateVars);

  double slope = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression

    const double a = properties[p_a];
    const double slipRate0 = properties[p_slipRate0];

    const double slipRateLinear = _minSlipRate;
    const double slipRateEff = std::max(slipRate, slipRateLinear);

    //slope = -slipRateEff / (slipRate0 * normalTraction * a);
    slope = -normalTraction * a; // log space
  } // if

  PetscLogFlops(5);

  return slope;
} // _calcFrictionSlope


// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::RateStateAgeing::_updateStateVars(const double slip,
						  const double slipRate,
						  const double normalTraction,
						  double* const stateVars,
						  const int numStateVars,
						  const double* properties,
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

  // Since regulatized friction -> 0 as slipRate -> 0, limit slip
  // rate to some minimum value
    const double slipRateEff = std::max(_minSlipRate, slipRate);

  const double dt = _dt;
  const double thetaTVertex = stateVars[s_state];
  const double L = properties[p_L];
  const double vDtL = slipRateEff * dt / L;
  const double expTerm = exp(-vDtL);

  double thetaTpdtVertex = 0.0;
  if (vDtL > 1.0e-20) {
    thetaTpdtVertex = thetaTVertex * expTerm + L / slipRateEff * (1 - expTerm);
    PetscLogFlops(7);
  } else {
    thetaTpdtVertex = thetaTVertex * expTerm + dt - 0.5 * slipRateEff/L * dt*dt;
    PetscLogFlops(9);
  } // if/else
  
  stateVars[s_state] = thetaTpdtVertex;

} // _updateStateVars


// End of file 
