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

#include "TestMeshIOCubit.hh" // Implementation of class methods

#include "pylith/meshio/MeshIOCubit.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES JournalingComponent

#include "catch2/catch_test_macros.hpp"

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor.
pylith::meshio::TestMeshIOCubit::TestMeshIOCubit(TestMeshIO_Data* data) :
    TestMeshIO(data) {
    _io = new MeshIOCubit();assert(_io);
    _io->PyreComponent::setIdentifier("TestMeshIOCubit");
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::meshio::TestMeshIOCubit::~TestMeshIOCubit(void) {
    delete _io;_io = nullptr;
} // destructor


// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestMeshIOCubit::testFilename(void) {
    PYLITH_METHOD_BEGIN;
    assert(_io);

    const std::string& filename = "hi.txt";
    _io->setFilename(filename.c_str());
    CHECK(filename == std::string(_io->getFilename()));

    PYLITH_METHOD_END;
} // testFilename


// ----------------------------------------------------------------------
// Test read().
void
pylith::meshio::TestMeshIOCubit::testRead(void) {
    PYLITH_METHOD_BEGIN;
    assert(_io);
    assert(_data);

    _io->setFilename(_data->filename.c_str());
    _io->setUseNodesetNames(true);

    // Read mesh
    delete _mesh;_mesh = new topology::Mesh;assert(_mesh);
    _io->read(_mesh);

    pythia::journal::debug_t debug("TestMeshIOCubit");
    if (debug.state()) {
        _mesh->view();
    } // if

    // Make sure mesh matches data
    _checkVals();

    PYLITH_METHOD_END;
} // testRead


// End of file
