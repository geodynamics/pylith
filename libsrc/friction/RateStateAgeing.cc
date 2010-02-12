// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
	{ "reference_slip_velocity", 1, pylith::topology::FieldBase::SCALAR },
        { "characteristic_slip_distance", 1, pylith::topology::FieldBase::SCALAR },
	{ "initial_slip_velocity", 1, pylith::topology::FieldBase::SCALAR },
	{ "constitutive_parameter_a", 1, pylith::topology::FieldBase::SCALAR },
        { "constitutive_parameter_b", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Number of State Variables.
      const int numStateVars = 2;

      // State Variables.
      const pylith::materials::Metadata::ParamDescription stateVars[] = {
	{ "state-variable", 1, pylith::topology::FieldBase::SCALAR },
	{ "initial-state-variable", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 6;
      const char* dbProperties[] = { "reference-friction-coefficient"
				     "reference-slip-velocity"
				     "characteristic-slip-distance"
				     "initial slip velocity"
				     "constitutive parameter a"
				     "constitutive parameter b"
 };      

      const int numDBStateVars = 2;
      const char* dbStateVars[] = { "state-variable"
				     "initial state variable"
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
const int pylith::friction::RateStateAgeing::p_slipRateIn = 
  pylith::friction::RateStateAgeing::p_L + 1;
const int pylith::friction::RateStateAgeing::p_a = 
  pylith::friction::RateStateAgeing::p_slipRateIn + 1;
const int pylith::friction::RateStateAgeing::p_b = 
  pylith::friction::RateStateAgeing::p_a + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::RateStateAgeing::db_coef = 0;
const int pylith::friction::RateStateAgeing::db_slipRate0 = 
  pylith::friction::RateStateAgeing::db_coef + 1;
const int pylith::friction::RateStateAgeing::db_L = 
  pylith::friction::RateStateAgeing::db_slipRate0 + 1;
const int pylith::friction::RateStateAgeing::db_slipRateIn = 
  pylith::friction::RateStateAgeing::db_L + 1;
const int pylith::friction::RateStateAgeing::db_a = 
  pylith::friction::RateStateAgeing::db_slipRateIn + 1;
const int pylith::friction::RateStateAgeing::db_b = 
  pylith::friction::RateStateAgeing::db_a + 1;

// Indices of state variables.
const int pylith::friction::RateStateAgeing::s_state = 0;
const int pylith::friction::RateStateAgeing::s_stateIn = 
  pylith::friction::RateStateAgeing::s_state + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::RateStateAgeing::db_state = 0;
const int pylith::friction::RateStateAgeing::db_stateIn = 
  pylith::friction::RateStateAgeing::db_stateIn + 1;

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
				    _RateStateAgeing::numDBStateVars))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::RateStateAgeing::~RateStateAgeing(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::friction::RateStateAgeing::_dbToProperties(
					   double* const propValues,
					   const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_RateStateAgeing::numDBProperties == numDBValues);

  const double db_fricCoef = dbValues[db_coef];
  const double db_slipVel0 = dbValues[db_slipRate0];
  const double db_dC = dbValues[db_L];
  const double db_slipVelIn = dbValues[db_slipRateIn];
  const double db_parA = dbValues[db_a];
  const double db_parB = dbValues[db_b];
 
  if (db_fricCoef <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for reference coefficient "
	<< "of Rate and State friction Ageing Law.\n"
	<< "reference coefficient of friction: " << db_fricCoef << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (db_dC <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for characteristic"
	<< "slip distance of Rate and State friction Ageing Law.\n"
	<< "characteristic slip distance: " << db_dC << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (db_parA <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
	<< "parameter 'a' of Rate and State friction Ageing Law.\n"
	<< "Rate and State parameter 'a' of Ageing Law of friction: " << db_parA << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (db_parB <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for constitutive "
	<< "parameter 'b' of Rate and State friction Ageing Law.\n"
	<< "Rate and State parameter 'a' of Ageing Law of friction: " << db_parB << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[p_coef] = db_fricCoef;
  propValues[p_slipRate0] = db_slipVel0;
  propValues[p_L] = db_dC;
  propValues[p_slipRateIn] = db_slipVelIn;
  propValues[p_a] = db_parA;
  propValues[p_b] = db_parB;

} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::RateStateAgeing::_nondimProperties(double* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _RateStateAgeing::numProperties);

  const double lengthScale = _normalizer->lengthScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();

  values[p_slipRate0] /= lengthScale / timeScale;
  values[p_L] /= lengthScale;
  values[p_slipRateIn] /= lengthScale / timeScale;
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::RateStateAgeing::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _RateStateAgeing::numProperties);

  const double lengthScale = _normalizer->lengthScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();

  values[p_slipRate0] *= lengthScale / timeScale;
  values[p_L] *= lengthScale;
  values[p_slipRateIn] *= lengthScale / timeScale;
} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
pylith::friction::RateStateAgeing::_dbToStateVars(
					   double* const stateValues,
					   const double_array& dbValues) const
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_RateStateAgeing::numDBStateVars == numDBValues);

  const double stateVariable = dbValues[db_state];
  const double initialStateVariable = dbValues[db_stateIn];
 
  stateValues[s_state] = stateVariable;
  stateValues[s_stateIn] = initialStateVariable;
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::friction::RateStateAgeing::_nondimStateVars(double* const values,
						    const int nvalues) const
{ // _nondimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _RateStateAgeing::numStateVars);

  const double timeScale = _normalizer->timeScale();

  values[s_state] /= timeScale;
  values[s_stateIn] /= timeScale;
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::friction::RateStateAgeing::_dimStateVars(double* const values,
						      const int nvalues) const
{ // _dimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _RateStateAgeing::numStateVars);

  const double timeScale = _normalizer->timeScale();

  values[s_state] *= timeScale;
  values[s_stateIn] *= timeScale;
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
  assert(0 != properties);
  assert(_numPropsVertex == numProperties);
  assert(0 != numStateVars);
  assert(_numVarsVertex == numStateVars);

  // :TODO: slipRate has to be at (n+1)
  double friction = 0.0;
  double mu_f = 0.0;
  if (normalTraction < 0.0) {
    // if fault is in compression
    const double f0 = properties[p_coef];

    const double slipRate0 = properties[p_slipRate0];
    const double a = properties[p_a];
    const double aLnTerm = a * log(slipRate / slipRate0);

    const double theta = stateVars[s_state];
    const double L = properties[p_L];
    const double b = properties[p_b];
    const double bLnTerm = b * log(slipRate0 * theta / L);

    mu_f = f0 + aLnTerm + bLnTerm;
    friction = -mu_f * normalTraction;
  } // if

  PetscLogFlops(5);

  return friction;
} // _calcFriction

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::RateStateAgeing::_updateStateVars(const double slip,
						  const double slipRate,
						  double* const stateVars,
						  const int numStateVars,
						  const double* properties,
						  const int numProperties)
{ // _updateStateVars

  assert(0 != numStateVars);
  assert(0 != numProperties);

  stateVars[s_stateIn] = stateVars[s_state];

  const double deltaT = _dt;
  const double thetaN = stateVars[s_stateIn];
  const double L = properties[p_L];
  const double expTerm = exp(-slipRate * deltaT / L);
 
  stateVars[s_state] = thetaN * expTerm +
                       L / slipRate * (1 - expTerm);
    
} // _updateStateVars


// End of file 
