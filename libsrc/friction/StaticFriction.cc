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

#include "StaticFriction.hh" // implementation of object methods

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
    namespace _StaticFriction {

      // Number of physical properties.
      const int numProperties = 1;

      // Physical properties.
      const pylith::materials::Metadata::ParamDescription properties[] = {
	{ "friction-coefficient", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 1;
      const char* dbProperties[] = { "friction-coefficient" };      
      
    } // _StaticFriction
  } // friction
} // pylith

// Indices of physical properties
const int pylith::friction::StaticFriction::p_coef = 0;

// Indices of database values (order must match dbProperties)
const int pylith::friction::StaticFriction::db_coef = 0;

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
pylith::friction::StaticFriction::_dbToProperties(
					   double* const propValues,
					   const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_StaticFriction::numDBProperties == numDBValues);

  const double coef = dbValues[db_coef];
 
  if (coef <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for coefficient "
	<< "of friction.\n"
	<< "coefficient of friction: " << coef << "\n";
    throw std::runtime_error(msg.str());
  } // if

  propValues[p_coef] = coef;
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::friction::StaticFriction::_nondimProperties(double* const values,
						    const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _StaticFriction::numProperties);

  // No dimensions
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::friction::StaticFriction::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _StaticFriction::numProperties);

  // No dimensions
} // _dimProperties

// ----------------------------------------------------------------------
// Compute friction from properties and state variables.
double
pylith::friction::StaticFriction::_calcFriction(const double slip,
						const double slipRate,
						const double normalTraction,
						const double* properties,
						const int numProperties,
						const double* stateVars,
						const int numStateVars)
{ // _calcFriction
  assert(0 != properties);
  assert(_numProps == numProperties);
  assert(0 == numStateVars);

  const double friction = (normalTraction < 0) ?
    properties[p_coef] * normalTraction : 0.0;

  PetscLogFlops(1);
} // _calcFriction


// End of file 
