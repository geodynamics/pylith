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

#include "TestSolutionFactory.hh" // Implementation of class methods

#include "pylith/problems/SolutionFactory.hh" // Test subject

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/testing/FieldTester.hh" // USES FieldTester
#include "pylith/utils/error.hh" // USES PyLITH_METHOD*

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::problems::TestSolutionFactory::setUp(void) {
    PYLITH_METHOD_BEGIN;
    _data = new TestSolutionFactory_Data();CPPUNIT_ASSERT(_data);

    CPPUNIT_ASSERT(_data->normalizer);
    _data->normalizer->setLengthScale(1.0e+03);
    _data->normalizer->setTimeScale(2.0);
    _data->normalizer->setDensityScale(3.0e+3);
    _data->normalizer->setPressureScale(2.25e+10);

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
        _data->normalizer->getLengthScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, -1, -1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["displacement"] = info;

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
        _data->normalizer->getLengthScale() / _data->normalizer->getTimeScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 3, -1, -1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, false
        );
    info.index = 1;
    _data->subfields["velocity"] = info;

    // pressure
    componentNames.resize(1);
    componentNames[0] = "pressure";
    info.description = pylith::topology::Field::Description(
        "pressure",
        "pressure",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getPressureScale()
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
        _data->normalizer->getLengthScale() / _data->normalizer->getLengthScale()
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
        _data->normalizer->getPressureScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, -1, -1, true, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 1;
    _data->subfields["lagrange_multiplier_fault"] = info;

    // temperature
    componentNames.resize(1);
    componentNames[0] = "temperature";
    info.description = pylith::topology::Field::Description(
        "temperature",
        "temperature",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getTemperatureScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 3, -1, -1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POINT_SPACE, true
        );
    info.index = 1;
    _data->subfields["temperature"] = info;

    PYLITH_METHOD_END;
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::problems::TestSolutionFactory::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _factory;_factory = NULL;
    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _solution;_solution = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Test adding displacement and velocity subfields.
void
pylith::problems::TestSolutionFactory::testDispVel(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_factory);
    CPPUNIT_ASSERT(_data);

    CPPUNIT_ASSERT(!_solution->hasSubfield("displacement"));
    CPPUNIT_ASSERT(!_solution->hasSubfield("velocity"));

    _factory->addDisplacement(_data->subfields["displacement"].fe);
    _factory->addVelocity(_data->subfields["velocity"].fe);

    CPPUNIT_ASSERT(_data->normalizer);

    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["displacement"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["velocity"]);

    PYLITH_METHOD_END;
} // testDispVel


// ---------------------------------------------------------------------------------------------------------------------
// Test adding displacement and fault Lagrange multiplier subfields.
void
pylith::problems::TestSolutionFactory::testDispLagrangeFault(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_factory);
    CPPUNIT_ASSERT(_data);

    CPPUNIT_ASSERT(!_solution->hasSubfield("displacement"));
    CPPUNIT_ASSERT(!_solution->hasSubfield("lagrange_multiplier_fault"));

    _factory->addDisplacement(_data->subfields["displacement"].fe);
    _factory->addLagrangeMultiplierFault(_data->subfields["lagrange_multiplier_fault"].fe);

    CPPUNIT_ASSERT(_data->normalizer);

    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["displacement"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["lagrange_multiplier_fault"]);

    PYLITH_METHOD_END;
} // testDispVel


// ---------------------------------------------------------------------------------------------------------------------
// Test adding pressure and trace strain subfields.
void
pylith::problems::TestSolutionFactory::testPressure(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_factory);

    CPPUNIT_ASSERT(!_solution->hasSubfield("pressure"));
    CPPUNIT_ASSERT(!_solution->hasSubfield("trace_strain"));

    _factory->addPressure(_data->subfields["pressure"].fe);
    _factory->addTraceStrain(_data->subfields["trace_strain"].fe);

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);

    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["pressure"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["trace_strain"]);

    PYLITH_METHOD_END;
} // testPressure


// ---------------------------------------------------------------------------------------------------------------------
// Test adding temperature subfields.
void
pylith::problems::TestSolutionFactory::testDispTemp(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_factory);

    CPPUNIT_ASSERT(!_solution->hasSubfield("displacement"));
    CPPUNIT_ASSERT(!_solution->hasSubfield("temperature"));

    _factory->addDisplacement(_data->subfields["displacement"].fe);
    _factory->addTemperature(_data->subfields["temperature"].fe);

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);

    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["displacement"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["temperature"]);

    PYLITH_METHOD_END;
} // testTemperature


// ---------------------------------------------------------------------------------------------------------------------
// Test setValues().
void
pylith::problems::TestSolutionFactory::testSetValues(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_factory);

    _factory->addDisplacement(_data->subfields["displacement"].fe);
    _factory->addPressure(_data->subfields["pressure"].fe);
    _solution->subfieldsSetup();
    _solution->createDiscretization();
    _solution->allocate();

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);

    CPPUNIT_ASSERT(_data->solutionDB);
    _factory->setValues(_data->solutionDB);
    pylith::testing::FieldTester::checkFieldWithDB(*_solution, _data->solutionDB, _data->normalizer->getLengthScale());

    PYLITH_METHOD_END;
} // testSetValues


// ---------------------------------------------------------------------------------------------------------------------
// Initialze mesh, coordinate system, solution, and factory.
void
pylith::problems::TestSolutionFactory::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_data);

    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(_data->meshFilename);
    iohandler.filename(_data->meshFilename);
    _mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    iohandler.read(_mesh);

    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any cells.",
                           pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any vertices.",
                           pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    // Setup coordinates.
    _mesh->setCoordSys(_data->cs);
    CPPUNIT_ASSERT(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    _solution = new pylith::topology::Field(*_mesh);CPPUNIT_ASSERT(_solution);
    _solution->setLabel("solution");
    _factory = new SolutionFactory(*_solution, *_data->normalizer);

    PYLITH_METHOD_END;
} // _initialize


// ---------------------------------------------------------------------------------------------------------------------
pylith::problems::TestSolutionFactory_Data::TestSolutionFactory_Data(void) :
    meshFilename(NULL),
    cs(NULL),
    normalizer(new spatialdata::units::Nondimensional),
    solutionDB(new spatialdata::spatialdb::UserFunctionDB)
{}


// ---------------------------------------------------------------------------------------------------------------------
pylith::problems::TestSolutionFactory_Data::~TestSolutionFactory_Data(void) {
    delete cs;cs = NULL;
    delete normalizer;normalizer = NULL;
    delete solutionDB;solutionDB = NULL;
}


// End of file
