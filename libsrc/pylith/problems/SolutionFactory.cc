// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/problems/SolutionFactory.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Scales.hh" // USES Scales
#include "spatialdata/units/ElasticityScales.hh" // USES ElasticityScales
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <typeinfo> // USES "<<typeid()
#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::problems::SolutionFactory::SolutionFactory(pylith::topology::Field& solution,
                                                   const spatialdata::units::Scales& scales) :
    _solution(solution),
    _scales(scales),
    _spaceDim(solution.getSpaceDim()) {
    GenericComponent::setName("solutionfactory");
    assert(1 <= _spaceDim && _spaceDim <= 3);
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::problems::SolutionFactory::~SolutionFactory(void) {}


// ------------------------------------------------------------------------------------------------
// Add displacement subfield to solution field.
void
pylith::problems::SolutionFactory::addDisplacement(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addDisplacement(discretization="<<typeid(discretization).name()<<")");

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
    description.scale = _scales.getDisplacementScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addDisplacement


// ------------------------------------------------------------------------------------------------
// Add velocity subfield to solution field.
void
pylith::problems::SolutionFactory::addVelocity(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addVelocity(discretization="<<typeid(discretization).name()<<")");

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
    description.scale = _scales.getDisplacementScale() / _scales.getTimeScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addVelocity


// ------------------------------------------------------------------------------------------------
// Add pressure subfield to solution field.
void
pylith::problems::SolutionFactory::addPressure(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPressure(discretization="<<typeid(discretization).name()<<")");

    const char* fieldName = "pressure";
    const char* componentNames[1] = { "pressure" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = spatialdata::units::ElasticityScales::getStressScale(_scales);
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addPressure


// ------------------------------------------------------------------------------------------------
// Add time derivative of pressure subfield to solution field.
void
pylith::problems::SolutionFactory::addPressureDot(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPressureDot(discretization="<<typeid(discretization).name()<<")");

    const char* fieldName = "pressure_t";
    const char* componentNames[1] = { "pressure_t" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _scales.getPressureScale() / _scales.getTimeScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addPressureDot


// ------------------------------------------------------------------------------------------------
// Add trace strain subfield to solution field.
void
pylith::problems::SolutionFactory::addTraceStrain(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTraceStrain(discretization="<<typeid(discretization).name()<<")");

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


// ------------------------------------------------------------------------------------------------
// Add time derivative of trace strain subfield to solution field.
void
pylith::problems::SolutionFactory::addTraceStrainDot(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTraceStrain(discretization="<<typeid(discretization).name()<<")");

    const char* fieldName = "trace_strain_t";
    const char* componentNames[1] = { "trace_strain_t" };
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
} // addTraceStrainDot


// ------------------------------------------------------------------------------------------------
// Add fault Lagrange multiplier subfield to solution field.
void
pylith::problems::SolutionFactory::addLagrangeMultiplierFault(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addLagrangeMultiplierFault(discretization="<<typeid(discretization).name()<<")");

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
    description.scale = spatialdata::units::ElasticityScales::getStressScale(_scales);
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addLagrangeMultiplierFault


// ------------------------------------------------------------------------------------------------
// Add temperature subfield to solution field.
void
pylith::problems::SolutionFactory::addTemperature(const pylith::topology::Field::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTemperature(discretization="<<typeid(discretization).name()<<")");

    const char* fieldName = "temperature";
    const char* componentNames[1] = { "temperature" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = componentNames[0];
    description.scale = _scales.getTemperatureScale();
    description.validator = NULL;

    _solution.subfieldAdd(description, discretization);

    PYLITH_METHOD_END;
} // addTemperature


// ------------------------------------------------------------------------------------------------
void
pylith::problems::SolutionFactory::setValues(spatialdata::spatialdb::SpatialDB* db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setValues(db="<<typeid(*db).name()<<")");

    // Set solution field to zero.
    _solution.zeroLocal();

    pylith::topology::FieldQuery query(_solution);
    query.initializeWithDefaultQueries();
    query.openDB(db, _scales.getLengthScale());
    query.queryDB();
    query.closeDB(db);

    PYLITH_METHOD_END;
} // setValues


// End of file
