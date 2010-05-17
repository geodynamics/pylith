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

#include "ViscousFriction.hh" // implementation of object methods

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
    namespace _ViscousFriction {

      // Number of physical properties.
      const int numProperties = 3;

      // Friction model parameters.
      const pylith::materials::Metadata::ParamDescription properties[] = {
	{ "static_coefficient", 1, pylith::topology::FieldBase::SCALAR },
	{ "reference_slip_rate", 1, pylith::topology::FieldBase::SCALAR },
	{ "cohesion", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Number of state variables.
      const int numStateVars = 1;

      // State variables.
      const pylith::materials::Metadata::ParamDescription stateVars[] = {
	{ "slip_rate", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 3;
      const char* dbProperties[3] = { "viscous-coefficient",
				      "reference-slip-rate",
				      "cohesion",
      };

      const int numDBStateVars = 1;
      const char* dbStateVars[1] = { "slip-rate",
      };      
      
    } // _ViscousFriction
  } // friction
} // pylith

// Indices of physical properties
const int pylith::friction::ViscousFriction::p_coefS = 0;
const int pylith::friction::ViscousFriction::p_v0 = 
  pylith::friction::ViscousFriction::p_coefS + 1;
const int pylith::friction::ViscousFriction::p_cohesion =
  pylith::friction::ViscousFriction::p_v0 + 1;

// Indices of database values (order must match dbProperties)
const int pylith::friction::ViscousFriction::db_coefS = 0;
const int pylith::friction::ViscousFriction::db_v0 =
  pylith::friction::ViscousFriction::db_coefS + 1;
const int pylith::friction::ViscousFriction::db_cohesion =
  pylith::friction::ViscousFriction::db_v0 + 1;

// Indices of state variables.
const int pylith::friction::ViscousFriction::s_slipRate = 0;

// Indices of database values (order must match dbProperties)
const int pylith::friction::ViscousFriction::db_slipRate = 0;

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::ViscousFriction::ViscousFriction(void) :
  FrictionModel(materials::Metadata(_ViscousFriction::properties,
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
pylith::friction::ViscousFriction::~ViscousFriction(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::friction::ViscousFriction::_dbToProperties(
					   double* const propValues,
					   const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_ViscousFriction::numDBProperties == numDBValues);

  const double coefS = dbValues[db_coefS];
  const double v0 = dbValues[db_v0];
  const double cohesion = dbValues[db_cohesion];

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

  propValues[p_coefS] = coefS;
  propValues[p_v0] = db_v0;
  propValues[p_cohesion] = cohesion;

} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::ViscousFriction::_nondimProperties(double* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ViscousFriction::numProperties);

  const double lengthScale = _normalizer->lengthScale();
  const double timeScale = _normalizer->timeScale();
  const double pressureScale = _normalizer->pressureScale();
  const double velocityScale = lengthScale / timeScale;

  values[p_v0] /= velocityScale;
  values[p_cohesion] /= pressureScale;
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::ViscousFriction::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ViscousFriction::numProperties);

  const double lengthScale = _normalizer->lengthScale();
  const double timeScale = _normalizer->timeScale();
  const double pressureScale = _normalizer->pressureScale();
  const double velocityScale = lengthScale / timeScale;

  values[p_v0] *= velocityScale;
  values[p_cohesion] *= pressureScale;
} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
pylith::friction::ViscousFriction::_dbToStateVars(
					   double* const stateValues,
					   const double_array& dbValues) const
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_ViscousFriction::numDBStateVars == numDBValues);

  const double slipRate = dbValues[db_slipRate];
 
  stateValues[s_slipRate] = slipRate;
} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::friction::ViscousFriction::_nondimStateVars(double* const values,
						  const int nvalues) const
{ // _nondimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ViscousFriction::numStateVars);

  const double lengthScale = _normalizer->lengthScale();
  const double timeScale = _normalizer->timeScale();
  const double velocityScale = lengthScale / timeScale;

  values[s_slipRate] /= velocityScale;
} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::friction::ViscousFriction::_dimStateVars(double* const values,
					       const int nvalues) const
{ // _dimStateVars
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ViscousFriction::numStateVars);

  const double lengthScale = _normalizer->lengthScale();
  const double timeScale = _normalizer->timeScale();
  const double velocityScale = lengthScale / timeScale;

  values[s_slipRate] *= velocityScale;
} // _dimStateVars

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
double
pylith::friction::ViscousFriction::_calcFriction(const double slip,
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

  double friction = 0.0;
  double mu_f = 0.0;
  if (normalTraction <= 0.0) {
    // if fault is in compression
    mu_f = properties[p_coefS] * (1.0 + fabs(slipRate) / properties[p_v0]);
    friction = - mu_f * normalTraction + properties[p_cohesion];
  } // if

  PetscLogFlops(6);

  return friction;
} // _calcFriction

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::ViscousFriction::_updateStateVars(const double slip,
						  const double slipRate,
						  const double normalTraction,
						  double* const stateVars,
						  const int numStateVars,
						  const double* properties,
						  const int numProperties)
{ // _updateStateVars
  assert(0 != numStateVars);
  assert(0 != numProperties);

  stateVars[s_slipRate] = stateVars[s_slipRate]; 
} // _updateStateVars


// End of file 
