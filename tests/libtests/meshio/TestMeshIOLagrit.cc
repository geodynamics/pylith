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

#include "TestMeshIOLagrit.hh" // Implementation of class methods

#include "pylith/meshio/MeshIOLagrit.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES JournalingComponent

#include <strings.h> // USES strcasecmp()

#if defined(WORDS_BIGENDIAN)
#define NATIVE_BIG_ENDIAN
#else
#define NATIVE_LITTLE_ENDIAN
#endif

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestMeshIOLagrit::setUp(void) {
    TestMeshIO::setUp();
    _io = new MeshIOLagrit();CPPUNIT_ASSERT(_io);
    _data = NULL;

    _io->PyreComponent::setIdentifier("TestMeshIOLagrit");
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::meshio::TestMeshIOLagrit::tearDown(void) {
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
pylith::meshio::TestMeshIOLagrit::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    MeshIOLagrit iohandler;

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test debug()
void
pylith::meshio::TestMeshIOLagrit::testDebug(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);

    _testDebug(*_io);

    PYLITH_METHOD_END;
} // testDebug


// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestMeshIOLagrit::testFilename(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);

    const char* filenameGmv = "hi.txt";
    const char* filenamePset = "hi2.txt";
    _io->filenameGmv(filenameGmv);
    _io->filenamePset(filenamePset);

    CPPUNIT_ASSERT(0 == strcasecmp(filenameGmv, _io->filenameGmv()));
    CPPUNIT_ASSERT(0 == strcasecmp(filenamePset, _io->filenamePset()));

    PYLITH_METHOD_END;
} // testFilename


// ----------------------------------------------------------------------
// Test read().
void
pylith::meshio::TestMeshIOLagrit::testRead(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_io);
    CPPUNIT_ASSERT(_data);

    _io->filenameGmv(_data->filenameGmv);
    _io->filenamePset(_data->filenamePset);
    _io->ioInt32(_data->ioInt32);
    _io->isRecordHeader32Bit(_data->isRecordHeader32Bit);

    // LaGriT file was created on little endian machine, so flip endian if
    // running test on big endian machine.
#if defined(NATIVE_LITTLE_ENDIAN)
    _io->flipEndian(false);
#else
    _io->flipEndian(true);
#endif

    // Read mesh
    delete _mesh;_mesh = new topology::Mesh;
    _io->read(_mesh);

    // Make sure mesh matches data
    _checkVals();

    PYLITH_METHOD_END;
} // testReadTetAscii


// ----------------------------------------------------------------------
// Get test data.
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOLagrit::_getData(void) {
    return _data;
} // _data


// ----------------------------------------------------------------------
// Constructor
pylith::meshio::TestMeshIOLagrit_Data::TestMeshIOLagrit_Data(void) :
    filenameGmv(NULL),
    filenamePset(NULL) {} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestMeshIOLagrit_Data::~TestMeshIOLagrit_Data(void) {}


// End of file
