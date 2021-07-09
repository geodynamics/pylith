// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "SolutionFactory.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
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
    _spaceDim(solution.getSpaceDim()) {
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
    description.scale = _normalizer.getLengthScale();
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
    description.scale = _normalizer.getLengthScale() / _normalizer.getTimeScale();
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
    description.scale = _normalizer.getPressureScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addPressure


// ---------------------------------------------------------------------------------------------------------------------
// Add trace strain subfield to solution field.
void
pylith::problems::SolutionFactory::addTraceStrain(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("traceStrain(discretization=typeid(discretization).name())");

    const char* fieldName = "trace_strain";
    const char* componentNames[1] = { "trace_strain" };
    const PylithReal noScale = 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = noScale;
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addTraceStrain


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
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = _normalizer.getPressureScale();
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
    description.scale = _normalizer.getTemperatureScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addTemperature


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::problems::SolutionFactory::setValues(spatialdata::spatialdb::SpatialDB* db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setValues(db="<<typeid(*db).name()<<")");

    // Set solution field to zero.
    _solution.zeroLocal();

    pylith::topology::FieldQuery query(_solution);
    query.initializeWithDefaultQueries();
    query.openDB(db, _normalizer.getLengthScale());
    query.queryDB();
    query.closeDB(db);

    PYLITH_METHOD_END;
} // setValues


// End of file
