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

#include "TestMeshIOCubit.hh" // Implementation of class methods

#include "pylith/meshio/MeshIOCubit.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES JournalingComponent

#include <strings.h> // USES strcasecmp()

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestMeshIOCubit::setUp(void) {
    TestMeshIO::setUp();
    _io = new MeshIOCubit();CPPUNIT_ASSERT(_io);
    _data = NULL;

    _io->PyreComponent::setIdentifier("TestMeshIOCubit");
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::meshio::TestMeshIOCubit::tearDown(void) {
    const char* journalName = _io->PyreComponent::getName();
    pythia::journal::debug_t debug(journalName);
    debug.deactivate(); // DEBUGGING

    TestMeshIO::tearDown();

    delete _io;_io = NULL;
    delete _data;_data = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestMeshIOCubit::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    MeshIOCubit iohandler;

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test debug()
void
pylith::meshio::TestMeshIOCubit::testDebug(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);
    _testDebug(*_io);

    PYLITH_METHOD_END;
} // testDebug


// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestMeshIOCubit::testFilename(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);

    const char* filename = "hi.txt";
    _io->filename(filename);
    CPPUNIT_ASSERT(0 == strcasecmp(filename, _io->filename()));

    PYLITH_METHOD_END;
} // testFilename


// ----------------------------------------------------------------------
// Test read().
void
pylith::meshio::TestMeshIOCubit::testRead(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);
    CPPUNIT_ASSERT(_data);

    _io->filename(_data->filename);
    _io->useNodesetNames(true);

    // Read mesh
    delete _mesh;_mesh = new topology::Mesh;CPPUNIT_ASSERT(_mesh);
    _io->read(_mesh);

    // Make sure mesh matches data
    _checkVals();

    PYLITH_METHOD_END;
} // testRead


// ----------------------------------------------------------------------
// Get test data.
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit::_getData(void) {
    return _data;
} // _data


// ----------------------------------------------------------------------
// Constructor
pylith::meshio::TestMeshIOCubit_Data::TestMeshIOCubit_Data(void) :
    filename(NULL) {} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestMeshIOCubit_Data::~TestMeshIOCubit_Data(void) {}


// End of file
