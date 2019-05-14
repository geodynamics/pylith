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

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::problems::SolutionFactory::SolutionFactory(pylith::topology::Field& solution,
                                                   const spatialdata::units::Nondimensional& normalizer) :
    _solution(solution),
    _normalizer(normalizer),
    _spaceDim(solution.spaceDim()) {
    GenericComponent::setName("solutionfactory");
    assert(1 <= _spaceDim && _spaceDim <= 3);
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::problems::SolutionFactory::~SolutionFactory(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add displacement subfield to solution field.
void
pylith::problems::SolutionFactory::addDisplacement(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("displacement(discretization=typeid(discretization).name())");

    const char* fieldName = "displacement";
    const char* componentNames[3] = { "displacement_x", "displacement_y", "displacement_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
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
} // addDisplacement


// ---------------------------------------------------------------------------------------------------------------------
// Add time derivative of velocity subfield to solution field.
void
pylith::problems::SolutionFactory::addVelocity(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("velocity(discretization=typeid(discretization).name())");

    const char* fieldName = "velocity";
    const char* componentNames[3] = { "velocity_x", "velocity_y", "velocity_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
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
} // addVelocity


// ---------------------------------------------------------------------------------------------------------------------
// Add time derivative of pressure subfield to solution field.
void
pylith::problems::SolutionFactory::addPressure(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("pressure(discretization=typeid(discretization).name())");

    const char* fieldName = "pressure";
    const char* componentNames[1] = { "pressure" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.pressureScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addPressure


// ---------------------------------------------------------------------------------------------------------------------
// Add fluid pressure subfield to solution field.
void
pylith::problems::SolutionFactory::addFluidPressure(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("fluidPressure(discretization=typeid(discretization).name())");

    const char* fieldName = "fluid_pressure";
    const char* componentNames[1] = { "fluid_pressure" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.pressureScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addFluidPressure


// ---------------------------------------------------------------------------------------------------------------------
// Add fault Lagrange multiplier subfield to solution field.
void
pylith::problems::SolutionFactory::addLagrangeMultiplierFault(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addLagrangeMultiplierFault(discretization=typeid(discretization).name())");

    const char* fieldName = "lagrange_multiplier_fault";
    const char* componentNames[3] = {
        "lagrange_multiplier_fault_x",
        "lagrange_multiplier_fault_y",
        "lagrange_multiplier_fault_z",
    };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = 3;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = _normalizer.pressureScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addFluidPressure


// ---------------------------------------------------------------------------------------------------------------------
// Add temperature subfield to solution field.
void
pylith::problems::SolutionFactory::addTemperature(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("temperature(discretization=typeid(discretization).name())");

    const char* fieldName = "temperature";
    const char* componentNames[1] = { "temperature" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.temperatureScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addTemperature


// ---------------------------------------------------------------------------------------------------------------------
// Add time derivative of displacement subfield to solution field.
void
pylith::problems::SolutionFactory::addDisplacementDot(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("displacementDot(discretization=typeid(discretization).name())");

    const char* fieldName = "displacement_dot";
    const char* componentNames[3] = { "displacement_dot_x", "displacement_dot_y", "displacement_dot_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
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
} // addDisplacementDot


// ---------------------------------------------------------------------------------------------------------------------
// Add time derivative of velocity subfield to solution field.
void
pylith::problems::SolutionFactory::addVelocityDot(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("velocityDot(discretization=typeid(discretization).name())");

    const char* fieldName = "velocity_dot";
    const char* componentNames[3] = { "velocity_dot_x", "velocity_dot_y", "velocity_dot_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
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
} // addVelocityDot


// ---------------------------------------------------------------------------------------------------------------------
// Add time derivative of pressure subfield to solution field.
void
pylith::problems::SolutionFactory::addPressureDot(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("pressureDot(discretization=typeid(discretization).name())");

    const char* fieldName = "pressure_dot";
    const char* componentNames[1] = { "pressure_dot" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.pressureScale() / _normalizer.timeScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addPressureDot


// ---------------------------------------------------------------------------------------------------------------------
// Add time derivative of fluid pressure subfield to solution field.
void
pylith::problems::SolutionFactory::addFluidPressureDot(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("fluidPressureDot(discretization=typeid(discretization).name())");

    const char* fieldName = "fluid_pressure_dot";
    const char* componentNames[1] = { "fluid_pressure_dot" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.pressureScale() / _normalizer.timeScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addFluidPressureDot


// ---------------------------------------------------------------------------------------------------------------------
// Add time derivative of temperature subfield to solution field.
void
pylith::problems::SolutionFactory::addTemperatureDot(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("temperatureDot(discretization=typeid(discretization).name())");

    const char* fieldName = "temperature_dot";
    const char* componentNames[1] = { "temperature_dot" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _normalizer.temperatureScale() / _normalizer.timeScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addTemperatureDot


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::problems::SolutionFactory::setValues(spatialdata::spatialdb::SpatialDB* db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setValues(db="<<typeid(*db).name()<<")");

    // Set solution field to zero.
    _solution.zeroLocal();

    pylith::topology::FieldQuery query(_solution);
    query.initializeWithDefaultQueryFns();
    query.openDB(db, _normalizer.lengthScale());
    query.queryDB();
    query.closeDB(db);

    PYLITH_METHOD_END;
} // setValues


// End of file
