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

#include "SolutionFactory.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <typeinfo> // USES typeid()
#include <cassert>

// ----------------------------------------------------------------------
const char* pylith::problems::SolutionFactory::_genericComponent = "solutionfactory";

// ----------------------------------------------------------------------
// Default constructor.
pylith::problems::SolutionFactory::SolutionFactory(pylith::topology::Field& solution,
                                                   const spatialdata::units::Nondimensional& normalizer) :
    _solution(solution),
    _normalizer(normalizer),
    _spaceDim(solution.spaceDim())
{ // constructor
    GenericComponent::name(_genericComponent);
    assert(1 <= _spaceDim && _spaceDim <= 3);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::problems::SolutionFactory::~SolutionFactory(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Add time derivative of displacement field to solution field.
void
pylith::problems::SolutionFactory::displacement(const pylith::topology::Field::Discretization& discretization)
{ // displacement
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("displacement(discretization=typeid(discretization).name())");

    const char* fieldName = "displacement";
    const char* componentNames[3] = { "displacement_x", "displacement_y", "displacement_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = _normalizer.lengthScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // displacement


// ----------------------------------------------------------------------
// Add time derivative of velocity field to solution field.
void
pylith::problems::SolutionFactory::velocity(const pylith::topology::Field::Discretization& discretization)
{ // velocity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("velocity(discretization=typeid(discretization).name())");

    const char* fieldName = "velocity";
    const char* componentNames[3] = { "velocity_x", "velocity_y", "velocity_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = _normalizer.lengthScale() / _normalizer.timeScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // velocity


// ----------------------------------------------------------------------
// Add time derivative of pressure field to solution field.
void
pylith::problems::SolutionFactory::pressure(const pylith::topology::Field::Discretization& discretization)
{ // pressure
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("pressure(discretization=typeid(discretization).name())");

    const char* fieldName = "pressure";
    const char* componentNames[1] = { "pressure" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.pressureScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // pressure


// ----------------------------------------------------------------------
// Add fluid pressure field to solution field.
void
pylith::problems::SolutionFactory::fluidPressure(const pylith::topology::Field::Discretization& discretization)
{ // fluidPressure
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("fluidPressure(discretization=typeid(discretization).name())");

    const char* fieldName = "fluid_pressure";
    const char* componentNames[1] = { "fluid_pressure" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.pressureScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // fluidPressure


// ----------------------------------------------------------------------
// Add temperature field to solution field.
void
pylith::problems::SolutionFactory::temperature(const pylith::topology::Field::Discretization& discretization)
{ // temperature
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("temperature(discretization=typeid(discretization).name())");

    const char* fieldName = "temperature";
    const char* componentNames[1] = { "temperature" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.temperatureScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // temperature


// ----------------------------------------------------------------------
// Add time derivative of displacement field to solution field.
void
pylith::problems::SolutionFactory::displacementDot(const pylith::topology::Field::Discretization& discretization)
{ // displacementDot
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("displacementDot(discretization=typeid(discretization).name())");

    const char* fieldName = "displacement_dot";
    const char* componentNames[3] = { "displacement_dot_x", "displacement_dot_y", "displacement_dot_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = _normalizer.lengthScale() / _normalizer.timeScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // displacementDot


// ----------------------------------------------------------------------
// Add time derivative of velocity field to solution field.
void
pylith::problems::SolutionFactory::velocityDot(const pylith::topology::Field::Discretization& discretization)
{ // velocityDot
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("velocityDot(discretization=typeid(discretization).name())");

    const char* fieldName = "velocity_dot";
    const char* componentNames[3] = { "velocity_dot_x", "velocity_dot_y", "velocity_dot_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = _normalizer.lengthScale() / (_normalizer.timeScale() * _normalizer.timeScale());
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // velocityDot


// ----------------------------------------------------------------------
// Add time derivative of pressure field to solution field.
void
pylith::problems::SolutionFactory::pressureDot(const pylith::topology::Field::Discretization& discretization)
{ // pressureDot
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("pressureDot(discretization=typeid(discretization).name())");

    const char* fieldName = "pressure_dot";
    const char* componentNames[1] = { "pressure_dot" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.pressureScale() / _normalizer.timeScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // pressureDot


// ----------------------------------------------------------------------
// Add time derivative of fluid pressure field to solution field.
void
pylith::problems::SolutionFactory::fluidPressureDot(const pylith::topology::Field::Discretization& discretization)
{ // fluidPressureDot
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("fluidPressureDot(discretization=typeid(discretization).name())");

    const char* fieldName = "fluid_pressure_dot";
    const char* componentNames[1] = { "fluid_pressure_dot" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.pressureScale() / _normalizer.timeScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // fluidPressureDot


// ----------------------------------------------------------------------
// Add time derivative of temperature field to solution field.
void
pylith::problems::SolutionFactory::temperatureDot(const pylith::topology::Field::Discretization& discretization)
{ // temperatureDot
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("temperatureDot(discretization=typeid(discretization).name())");

    const char* fieldName = "temperature_dot";
    const char* componentNames[1] = { "temperature_dot" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.temperatureScale() / _normalizer.timeScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // temperatureDot


// ----------------------------------------------------------------------
void
pylith::problems::SolutionFactory::setValues(spatialdata::spatialdb::SpatialDB* db)
{ // setValues
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setValues(db="<<typeid(*db).name()<<")");

    // Allocate solution field
    _solution.zeroLocal();

    pylith::topology::FieldQuery query(_solution);
    const pylith::string_vector& subfieldNames = _solution.subfieldNames();
    const size_t numSubfields = subfieldNames.size();
    for (size_t i = 0; i < numSubfields; ++i) {
        query.queryFn(subfieldNames[i].c_str(), pylith::topology::FieldQuery::dbQueryGeneric);
    } // for
    query.openDB(db, _normalizer.lengthScale());
    query.queryDB();
    query.closeDB(db);

    PYLITH_METHOD_END;
} // setValues


// End of file
