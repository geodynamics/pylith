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

#include "SlipWeakening.hh" // implementation of object methods

#include "pylith/materials/Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace friction {
    namespace _SlipWeakening {

      // Number of physical properties.
      const int numProperties = 3;

      // Physical properties.
      const pylith::materials::Metadata::ParamDescription properties[] = {
	{ "static-coefficient", 1, pylith::topology::FieldBase::SCALAR },
	{ "dynamic-coefficient", 2, pylith::topology::FieldBase::SCALAR },
        { "slip-weakening-parameter", 3, pylith::topology::FieldBase::SCALAR },
      };

      // Number of State Variables.
      const int numStateVars = 2;

      // State Variables.
      const pylith::materials::Metadata::ParamDescription stateVars[] = {
	{ "cumulative-slip", 1, pylith::topology::FieldBase::SCALAR },
	{ "previous-slip", 2, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 3;
      const char* dbProperties[] = { "static-coefficient"
				     "dynamic-coefficient"
				     "slip-weakening-parameter"
 };      

      const int numDBStateVars = 2;
      const char* dbStateVars[] = { "cumulative-slip"
				     "previous-slip"
};      
      
    } // _SlipWeakening
  } // friction
} // pylith

// Indices of physical properties
const int pylith::friction::SlipWeakening::p_coef = 0;

// Indices of database values (order must match dbProperties)
const int pylith::friction::SlipWeakening::db_coef = 0;

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::SlipWeakening::SlipWeakening(void) :
  FrictionModel(materials::Metadata(_SlipWeakening::properties,
				    _SlipWeakening::numProperties,
				    _SlipWeakening::dbProperties,
				    _SlipWeakening::numDBProperties,
				    0, 0,
				    0, 0))
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
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_SlipWeakening::numDBProperties == numDBValues);

  const double db_static = dbValues[db_coef];
  const double db_dynamic = dbValues[db_coef+1];
  const double db_d0 = dbValues[db_coef+2];
 
  if (db_static <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for static coefficient "
	<< "of friction.\n"
	<< "static coefficient of friction: " << db_static << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (db_dynamic <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for dynamic coefficient "
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

  propValues[p_coef] = db_static;
  propValues[p_coef+1] = db_dynamic;
  propValues[p_coef+2] = db_d0;
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::SlipWeakening::_nondimProperties(double* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);

  const double lengthScale = _normalizer->lengthScale();

  values[nvalues-1] = values[nvalues-1] / lengthScale;

  assert(nvalues == _SlipWeakening::numProperties);

} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::SlipWeakening::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);

  const double lengthScale = _normalizer->lengthScale();

  values[nvalues-1] = values[nvalues-1] * lengthScale;

  assert(nvalues == _SlipWeakening::numProperties);

} // _dimProperties

// ----------------------------------------------------------------------
// Compute state variables from values in spatial database.
void
pylith::friction::SlipWeakening::_dbToStateVars(
					   double* const stateValues,
					   const double_array& dbValues) const
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_SlipWeakening::numDBStateVars == numDBValues);

  const double cumulativeSlip = dbValues[db_coef+3];
  const double previousSlip = dbValues[db_coef+4];
 
  stateValues[0] = cumulativeSlip;
  stateValues[1] = previousSlip;

} // _dbToStateVars

// ----------------------------------------------------------------------
// Nondimensionalize state variables.
void
pylith::friction::SlipWeakening::_nondimStateVars(double* const values,
						    const int nvalues) const
{ // _nondimStateVars
  assert(0 != _normalizer);
  assert(0 != values);

  const double lengthScale = _normalizer->lengthScale();

  values[nvalues-1] = values[nvalues-1] / lengthScale;
  values[nvalues-2] = values[nvalues-2] / lengthScale;

  assert(nvalues == _SlipWeakening::numStateVars);

} // _nondimStateVars

// ----------------------------------------------------------------------
// Dimensionalize state variables.
void
pylith::friction::SlipWeakening::_dimStateVars(double* const values,
						      const int nvalues) const
{ // _dimStateVars
  assert(0 != _normalizer);
  assert(0 != values);

  const double lengthScale = _normalizer->lengthScale();

  values[nvalues-1] = values[nvalues-1] * lengthScale;
  values[nvalues-2] = values[nvalues-2] * lengthScale;

  assert(nvalues == _SlipWeakening::numStateVars);

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
  assert(0 != properties);
  assert(_numPropsVertex == numProperties);
  assert(0 != numStateVars);
  assert(_numVarsVertex == numStateVars);

  _updateStateVars(slip,slipRate,const_cast<double*>(&stateVars[0]),
		   numStateVars,&properties[0],numProperties);

  const double friction = (normalTraction < 0) ?
    ((stateVars[0] < properties[p_coef+2]) ?
     properties[p_coef]-(properties[p_coef]-properties[p_coef+1]) *
     stateVars[0] / properties[p_coef+2] * normalTraction : 
     properties[p_coef+1] * normalTraction) : 0.0;

  //  PetscLogFlops(1);

  return friction;
} // _calcFriction

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::SlipWeakening::_updateStateVars(const double slip,
						  const double slipRate,
						  double* const stateVars,
						  const int numStateVars,
						  const double* properties,
						  const int numProperties)
{ // _updateStateVars

  assert(0 != numStateVars);
  assert(0 != numProperties);

  const double tmpPreviousSlip = stateVars[1];
 
  stateVars[1] = stateVars[0];
  stateVars[0] += fabs(slip - tmpPreviousSlip);
    
} // _updateStateVars


// End of file 
