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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMeshIOPetsc.hh" // Implementation of class methods

#include "pylith/meshio/MeshIOPetsc.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES JournalingComponent

#include <strings.h> // USES strcasecmp()

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestMeshIOPetsc::setUp(void) {
    TestMeshIO::setUp();
    _io = new MeshIOPetsc();CPPUNIT_ASSERT(_io);
    _data = NULL;

    _io->PyreComponent::setIdentifier("TestMeshIOPetsc");
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::meshio::TestMeshIOPetsc::tearDown(void) {
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
pylith::meshio::TestMeshIOPetsc::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    MeshIOPetsc iohandler;

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test debug()
void
pylith::meshio::TestMeshIOPetsc::testDebug(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);
    _testDebug(*_io);

    PYLITH_METHOD_END;
} // testDebug


// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestMeshIOPetsc::testFilename(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);

    const std::string& filename = "hi.txt";
    _io->setFilename(filename.c_str());
    CPPUNIT_ASSERT_EQUAL(filename, std::string(_io->getFilename()));

    PYLITH_METHOD_END;
} // testFilename


// ----------------------------------------------------------------------
// Test read().
void
pylith::meshio::TestMeshIOPetsc::testRead(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);
    CPPUNIT_ASSERT(_data);

    _io->setFilename(_data->filename);

    // Read mesh
    delete _mesh;_mesh = new topology::Mesh;CPPUNIT_ASSERT(_mesh);
    _io->read(_mesh);

    pythia::journal::debug_t debug("TestMeshIOPetsc");
    if (debug.state()) {
        _mesh->view();
    } // if

    // Make sure mesh matches data
    _checkVals();

    PYLITH_METHOD_END;
} // testRead


// ----------------------------------------------------------------------
// Get test data.
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc::_getData(void) {
    return _data;
} // _data


// ----------------------------------------------------------------------
// Constructor
pylith::meshio::TestMeshIOPetsc_Data::TestMeshIOPetsc_Data(void) :
    filename(NULL) {} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestMeshIOPetsc_Data::~TestMeshIOPetsc_Data(void) {}


// End of file
