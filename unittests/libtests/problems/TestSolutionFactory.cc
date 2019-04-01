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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
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
    _data->normalizer->lengthScale(1.0e+03);
    _data->normalizer->timeScale(2.0);
    _data->normalizer->densityScale(3.0e+3);
    _data->normalizer->pressureScale(2.25e+10);

    pylith::topology::Field::SubfieldInfo info;
    info.dm = NULL;
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
        _data->normalizer->lengthScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, -1, true, pylith::topology::Field::POLYNOMIAL_SPACE
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
        _data->normalizer->lengthScale() / _data->normalizer->timeScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 3, -1, false, pylith::topology::Field::POLYNOMIAL_SPACE
        );
    info.index = 1;
    _data->subfields["velocity"] = info;

    // displacement_dot
    componentNames.resize(3);
    componentNames[0] = "displacement_dot_x";
    componentNames[1] = "displacement_dot_y";
    componentNames[2] = "displacement_dot_z";
    info.description = pylith::topology::Field::Description(
        "displacement_dot",
        "displacement_dot",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::VECTOR,
        _data->normalizer->lengthScale() / _data->normalizer->timeScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, -1, true, pylith::topology::Field::POLYNOMIAL_SPACE
        );
    info.index = 2;
    _data->subfields["displacement_dot"] = info;

    // velocity_dot
    componentNames.resize(3);
    componentNames[0] = "velocity_dot_x";
    componentNames[1] = "velocity_dot_y";
    componentNames[2] = "velocity_dot_z";
    info.description = pylith::topology::Field::Description(
        "velocity_dot",
        "velocity_dot",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::VECTOR,
        _data->normalizer->lengthScale() / pow(_data->normalizer->timeScale(), 2)
        );
    info.fe = pylith::topology::Field::Discretization(
        3, 3, -1, true, pylith::topology::Field::POINT_SPACE
        );
    info.index = 3;
    _data->subfields["velocity_dot"] = info;

    // pressure
    componentNames.resize(1);
    componentNames[0] = "pressure";
    info.description = pylith::topology::Field::Description(
        "pressure",
        "pressure",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->pressureScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, -1, true, pylith::topology::Field::POLYNOMIAL_SPACE
        );
    info.index = 0;
    _data->subfields["pressure"] = info;

    // fluid_pressure
    componentNames.resize(1);
    componentNames[0] = "fluid_pressure";
    info.description = pylith::topology::Field::Description(
        "fluid_pressure",
        "fluid_pressure",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->pressureScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 3, -1, true, pylith::topology::Field::POINT_SPACE
        );
    info.index = 1;
    _data->subfields["fluid_pressure"] = info;

    // pressure_dot
    componentNames.resize(1);
    componentNames[0] = "pressure_dot";
    info.description = pylith::topology::Field::Description(
        "pressure_dot",
        "pressure_dot",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->pressureScale() / _data->normalizer->timeScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, -1, true, pylith::topology::Field::POINT_SPACE
        );
    info.index = 2;
    _data->subfields["pressure_dot"] = info;

    // fluid_pressure_dot
    componentNames.resize(1);
    componentNames[0] = "fluid_pressure_dot";
    info.description = pylith::topology::Field::Description(
        "fluid_pressure_dot",
        "fluid_pressure_dot",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->pressureScale() / _data->normalizer->timeScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 3, -1, true, pylith::topology::Field::POINT_SPACE
        );
    info.index = 3;
    _data->subfields["fluid_pressure_dot"] = info;

    // temperature
    componentNames.resize(1);
    componentNames[0] = "temperature";
    info.description = pylith::topology::Field::Description(
        "temperature",
        "temperature",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->temperatureScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 3, -1, true, pylith::topology::Field::POINT_SPACE
        );
    info.index = 1;
    _data->subfields["temperature"] = info;

    // temperature_dot
    componentNames.resize(1);
    componentNames[0] = "temperature_dot";
    info.description = pylith::topology::Field::Description(
        "temperature_dot",
        "temperature_dot",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->temperatureScale() / _data->normalizer->timeScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 3, -1, true, pylith::topology::Field::POINT_SPACE
        );
    info.index = 3;
    _data->subfields["temperature_dot"] = info;

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
    _factory->addDisplacementDot(_data->subfields["displacement_dot"].fe);
    _factory->addVelocityDot(_data->subfields["velocity_dot"].fe);

    CPPUNIT_ASSERT(_data->normalizer);

    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["displacement"]);
    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["velocity"]);
    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["displacement_dot"]);
    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["velocity_dot"]);

    PYLITH_METHOD_END;
} // testDispVel


// ---------------------------------------------------------------------------------------------------------------------
// Test adding pressure and fluid pressure subfields.
void
pylith::problems::TestSolutionFactory::testPressure(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_factory);

    CPPUNIT_ASSERT(!_solution->hasSubfield("pressure"));
    CPPUNIT_ASSERT(!_solution->hasSubfield("fluid_pressure"));

    _factory->addPressure(_data->subfields["pressure"].fe);
    _factory->addFluidPressure(_data->subfields["fluid_pressure"].fe);
    _factory->addPressureDot(_data->subfields["pressure_dot"].fe);
    _factory->addFluidPressureDot(_data->subfields["fluid_pressure_dot"].fe);

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);

    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["pressure"]);
    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["fluid_pressure"]);
    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["pressure_dot"]);
    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["fluid_pressure_dot"]);

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
    _factory->addDisplacementDot(_data->subfields["displacement_dot"].fe);
    _factory->addTemperatureDot(_data->subfields["temperature_dot"].fe);

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);

    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["displacement"]);
    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["temperature"]);
    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["displacement_dot"]);
    pylith::topology::FieldTester::checkSubfieldInfo(*_solution, _data->subfields["temperature_dot"]);

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
    _solution->allocate();

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);

    CPPUNIT_ASSERT(_data->solutionDB);
    _factory->setValues(_data->solutionDB);
    pylith::topology::FieldTester::checkFieldWithDB(*_solution, _data->solutionDB, _data->normalizer->lengthScale());

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

    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any cells.", _mesh->numCells() > 0);
    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any vertices.", _mesh->numVertices() > 0);

    // Setup coordinates.
    _mesh->coordsys(_data->cs);
    CPPUNIT_ASSERT(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    _solution = new pylith::topology::Field(*_mesh);CPPUNIT_ASSERT(_solution);
    _solution->label("solution");
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
