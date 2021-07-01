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

#include "TestMeshIOAscii.hh" // Implementation of class methods

#include "pylith/meshio/MeshIOAscii.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES JournalingComponent

#include <strings.h> // USES strcasecmp()

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestMeshIOAscii::setUp(void) {
    TestMeshIO::setUp();
    _io = new MeshIOAscii();CPPUNIT_ASSERT(_io);
    _data = NULL;

    _io->PyreComponent::setIdentifier("TestMeshIOAscii");
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::meshio::TestMeshIOAscii::tearDown(void) {
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
pylith::meshio::TestMeshIOAscii::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    MeshIOAscii iohandler;

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test debug()
void
pylith::meshio::TestMeshIOAscii::testDebug(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);

    _testDebug(*_io);

    PYLITH_METHOD_END;
} // testDebug


// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestMeshIOAscii::testFilename(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);

    const char* filename = "hi.txt";
    _io->filename(filename);
    CPPUNIT_ASSERT(0 == strcasecmp(filename, _io->filename()));

    PYLITH_METHOD_END;
} // testFilename


// ----------------------------------------------------------------------
// Test write() and read().
void
pylith::meshio::TestMeshIOAscii::testWriteRead(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);
    CPPUNIT_ASSERT(_data);

    TestMeshIO::_createMesh();CPPUNIT_ASSERT(_mesh);

    // Write mesh
    CPPUNIT_ASSERT(_data->filename);
    _io->filename(_data->filename);
    _io->write(_mesh);

    // Read mesh
    delete _mesh;_mesh = new pylith::topology::Mesh;
    _io->read(_mesh);

    // Make sure meshIn matches data
    TestMeshIO::_checkVals();

    PYLITH_METHOD_END;
} // testWriteRead1D


// ----------------------------------------------------------------------
// Test read().
void
pylith::meshio::TestMeshIOAscii::testRead(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);
    CPPUNIT_ASSERT(_data);

    // Read mesh
    delete _mesh;_mesh = new pylith::topology::Mesh;
    CPPUNIT_ASSERT(_data->filename);
    _io->filename(_data->filename);
    _io->read(_mesh);

    // Make sure mesh matches data
    _checkVals();

    PYLITH_METHOD_END;
} // testRead3DIndexOne


// ----------------------------------------------------------------------
// Get test data.
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii::_getData(void) {
    return _data;
} // _data


// ----------------------------------------------------------------------
// Constructor
pylith::meshio::TestMeshIOAscii_Data::TestMeshIOAscii_Data(void) :
    filename(NULL) {} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestMeshIOAscii_Data::~TestMeshIOAscii_Data(void) {}


// End of file
