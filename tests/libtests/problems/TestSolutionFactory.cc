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

#include "TestSolutionFactory.hh" // Implementation of class methods

#include "pylith/problems/SolutionFactory.hh" // Test subject
#include "tests/src/FieldTester.hh" // USES FieldTester

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/error.hh" // USES PyLITH_METHOD*

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "pylith/scales/Scales.hh" // USES Scales

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
// Setup testing data.
pylith::problems::TestSolutionFactory::TestSolutionFactory(TestSolutionFactory_Data* data) :
    _data(data) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    assert(_data->scales);
    _data->scales->setLengthScale(1.0);
    _data->scales->setTimeScale(2.0);
    _data->scales->setRigidityScale(2.25e+6);

    pylith::topology::Field::SubfieldInfo info;
    pylith::string_vector componentNames;

    // displacement
    componentNames.resize(3);
    componentNames[0] = "displacement_x";
    componentNames[1] = "displacement_y";
    componentNames[2] = "displacement_z";
    info.description = pylith::topology::Field::Description(
        "displacement",
        "displacement",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::VECTOR,
        _data->scales->getLengthScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, -1, -1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["displacement"] = info;
    _data->subfields["displacement"].description.numComponents = _data->dimension;

    // velocity
    componentNames.resize(3);
    componentNames[0] = "velocity_x";
    componentNames[1] = "velocity_y";
    componentNames[2] = "velocity_z";
    info.description = pylith::topology::Field::Description(
        "velocity",
        "velocity",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::VECTOR,
        _data->scales->getLengthScale() / _data->scales->getTimeScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 3, -1, -1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, false
        );
    info.index = 1;
    _data->subfields["velocity"] = info;
    _data->subfields["velocity"].description.numComponents = _data->dimension;

    // pressure
    componentNames.resize(1);
    componentNames[0] = "pressure";
    info.description = pylith::topology::Field::Description(
        "pressure",
        "pressure",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->scales->getRigidityScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, -1, -1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["pressure"] = info;

    // trace_strain
    componentNames.resize(1);
    componentNames[0] = "trace_strain";
    info.description = pylith::topology::Field::Description(
        "trace_strain",
        "trace_strain",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->scales->getLengthScale() / _data->scales->getLengthScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, -1, -1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 1;
    _data->subfields["trace_strain"] = info;

    // pressure
    componentNames.resize(3);
    componentNames[0] = "lagrange_multiplier_fault_x";
    componentNames[1] = "lagrange_multiplier_fault_y";
    componentNames[2] = "lagrange_multiplier_fault_z";
    info.description = pylith::topology::Field::Description(
        "lagrange_multiplier_fault",
        "lagrange_multiplier_fault",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::VECTOR,
        _data->scales->getRigidityScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, -1, -1, true, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 1;
    _data->subfields["lagrange_multiplier_fault"] = info;
    _data->subfields["lagrange_multiplier_fault"].description.numComponents = _data->dimension;

    // temperature
    componentNames.resize(1);
    componentNames[0] = "temperature";
    info.description = pylith::topology::Field::Description(
        "temperature",
        "temperature",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->scales->getTemperatureScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 3, -1, -1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POINT_SPACE, true
        );
    info.index = 1;
    _data->subfields["temperature"] = info;

    _initialize();

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::problems::TestSolutionFactory::~TestSolutionFactory(void) {
    PYLITH_METHOD_BEGIN;

    delete _factory;_factory = NULL;
    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _solution;_solution = NULL;

    PYLITH_METHOD_END;
} // destructor


// ------------------------------------------------------------------------------------------------
// Test adding displacement and velocity subfields.
void
pylith::problems::TestSolutionFactory::testDispVel(void) {
    PYLITH_METHOD_BEGIN;
    assert(_factory);
    assert(_data);

    CHECK(!_solution->hasSubfield("displacement"));
    CHECK(!_solution->hasSubfield("velocity"));

    _factory->addDisplacement(_data->subfields["displacement"].fe);
    _factory->addVelocity(_data->subfields["velocity"].fe);

    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["displacement"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["velocity"]);

    PYLITH_METHOD_END;
} // testDispVel


// ------------------------------------------------------------------------------------------------
// Test adding displacement and fault Lagrange multiplier subfields.
void
pylith::problems::TestSolutionFactory::testDispLagrangeFault(void) {
    PYLITH_METHOD_BEGIN;
    assert(_factory);
    assert(_data);

    CHECK(!_solution->hasSubfield("displacement"));
    CHECK(!_solution->hasSubfield("lagrange_multiplier_fault"));

    _factory->addDisplacement(_data->subfields["displacement"].fe);
    _factory->addLagrangeMultiplierFault(_data->subfields["lagrange_multiplier_fault"].fe);

    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["displacement"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["lagrange_multiplier_fault"]);

    PYLITH_METHOD_END;
} // testDispVel


// ------------------------------------------------------------------------------------------------
// Test adding pressure and trace strain subfields.
void
pylith::problems::TestSolutionFactory::testPressure(void) {
    PYLITH_METHOD_BEGIN;
    assert(_factory);
    assert(_data);

    CHECK(!_solution->hasSubfield("pressure"));
    CHECK(!_solution->hasSubfield("trace_strain"));

    _factory->addPressure(_data->subfields["pressure"].fe);
    _factory->addTraceStrain(_data->subfields["trace_strain"].fe);

    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["pressure"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["trace_strain"]);

    PYLITH_METHOD_END;
} // testPressure


// ------------------------------------------------------------------------------------------------
// Test adding temperature subfields.
void
pylith::problems::TestSolutionFactory::testDispTemp(void) {
    PYLITH_METHOD_BEGIN;
    assert(_factory);
    assert(_data);

    CHECK(!_solution->hasSubfield("displacement"));
    CHECK(!_solution->hasSubfield("temperature"));

    _factory->addDisplacement(_data->subfields["displacement"].fe);
    _factory->addTemperature(_data->subfields["temperature"].fe);

    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["displacement"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["temperature"]);

    PYLITH_METHOD_END;
} // testTemperature


// ------------------------------------------------------------------------------------------------
// Test setValues().
void
pylith::problems::TestSolutionFactory::testSetValues(void) {
    PYLITH_METHOD_BEGIN;
    assert(_factory);
    assert(_data);

    _factory->addDisplacement(_data->subfields["displacement"].fe);
    _factory->addPressure(_data->subfields["pressure"].fe);
    _solution->subfieldsSetup();
    _solution->createDiscretization();
    _solution->allocate();

    assert(_data->solutionDB);
    _factory->setValues(_data->solutionDB);
    pylith::testing::FieldTester::checkFieldWithDB(*_solution, _data->solutionDB, _data->scales->getLengthScale());

    PYLITH_METHOD_END;
} // testSetValues


// ------------------------------------------------------------------------------------------------
// Initialize mesh, coordinate system, solution, and factory.
void
pylith::problems::TestSolutionFactory::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    assert(_data);

    pylith::meshio::MeshIOAscii iohandler;
    assert(_data->meshFilename);
    iohandler.setFilename(_data->meshFilename);
    _mesh = new pylith::topology::Mesh();assert(_mesh);
    iohandler.read(_mesh);

    assert(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    assert(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    // Setup coordinates.
    _mesh->setCoordSys(_data->cs);
    assert(_data->scales);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->scales);

    _solution = new pylith::topology::Field(*_mesh);assert(_solution);
    _solution->setLabel("solution");
    _factory = new SolutionFactory(*_solution, *_data->scales);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
pylith::problems::TestSolutionFactory_Data::TestSolutionFactory_Data(void) :
    meshFilename(NULL),
    cs(NULL),
    scales(new pylith::scales::Scales),
    solutionDB(new spatialdata::spatialdb::UserFunctionDB) {}


// ------------------------------------------------------------------------------------------------
pylith::problems::TestSolutionFactory_Data::~TestSolutionFactory_Data(void) {
    delete cs;cs = NULL;
    delete scales;scales = NULL;
    delete solutionDB;solutionDB = NULL;
}


// End of file
