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

#include "TestBoundaryMesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorSubmesh.hh" // USES SubmeshIS
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/testing/FaultCohesiveStub.hh" // USES FaultCohesiveStub

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestBoundaryMesh::setUp(void) { // setUp
    PYLITH_METHOD_BEGIN;

    _data = NULL;

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::bc::TestBoundaryMesh::tearDown(void) { // tearDown
    PYLITH_METHOD_BEGIN;

    delete _data;_data = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
// Test submesh() without fault.
void
pylith::bc::TestBoundaryMesh::testSubmesh(void) { // testSubmesh
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_data);
    PetscErrorCode err;

    pylith::topology::Mesh mesh;
    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(_data->filename);
    iohandler.filename(_data->filename);
    iohandler.read(&mesh);

    // Set up coordinates
    spatialdata::geocoords::CSCart cs;
    spatialdata::units::Nondimensional normalizer;
    cs.setSpaceDim(mesh.getDimension());
    mesh.setCoordSys(&cs);
    pylith::topology::MeshOps::nondimensionalize(&mesh, normalizer);

    // Create submesh
    CPPUNIT_ASSERT(_data->bcLabel);
    pylith::topology::Mesh submesh(mesh, _data->bcLabel);
    PetscDM dmMesh = submesh.getDM();CPPUNIT_ASSERT(dmMesh);

    // Check vertices
    pylith::topology::Stratum verticesStratum(dmMesh, pylith::topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();
    CPPUNIT_ASSERT_EQUAL(_data->numVerticesNoFault, verticesStratum.size());

    // Check cells
    pylith::topology::Stratum cellsStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 1);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    CPPUNIT_ASSERT_EQUAL(_data->numCells, cellsStratum.size());
    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscInt *closure = NULL;
        PetscInt closureSize, numVertices = 0;

        err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        for (PetscInt p = 0; p < closureSize*2; p += 2) {
            if ((closure[p] >= vStart) && (closure[p] < vEnd)) { numVertices++; }
        } // for
        err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        CPPUNIT_ASSERT_EQUAL(_data->numCorners, numVertices);
    } // for

    PYLITH_METHOD_END;
} // testSubmesh


// ----------------------------------------------------------------------
// Test submesh() with fault.
void
pylith::bc::TestBoundaryMesh::testSubmeshFault(void) { // testSubmeshFault
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_data);
    PetscErrorCode err;

    pylith::topology::Mesh mesh;
    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(_data->filename);
    iohandler.filename(_data->filename);
    iohandler.read(&mesh);

    // Set up coordinates
    spatialdata::geocoords::CSCart cs;
    spatialdata::units::Nondimensional normalizer;
    cs.setSpaceDim(mesh.getDimension());
    mesh.setCoordSys(&cs);
    pylith::topology::MeshOps::nondimensionalize(&mesh, normalizer);

    // Adjust topology
    CPPUNIT_ASSERT(_data->faultLabel);
    pylith::faults::FaultCohesiveStub fault;
    fault.setLabel(_data->faultLabel);
    fault.id(_data->faultId);
    fault.adjustTopology(&mesh);

    // Create submesh
    CPPUNIT_ASSERT(_data->bcLabel);
    pylith::topology::Mesh submesh(mesh, _data->bcLabel);
    PetscDM dmMesh = submesh.getDM();CPPUNIT_ASSERT(dmMesh);
#if 0 // DEBUGGING
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL);
    DMView(mesh.getDM(), PETSC_VIEWER_STDOUT_WORLD);
    PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
#endif

    // Check vertices
    pylith::topology::Stratum verticesStratum(dmMesh, pylith::topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();
    CPPUNIT_ASSERT_EQUAL(_data->numVerticesWithFault, verticesStratum.size());

    // Check cells
    topology::Stratum cellsStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 1);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    CPPUNIT_ASSERT_EQUAL(_data->numCells, cellsStratum.size());

    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscInt *closure = NULL;
        PetscInt closureSize, numVertices = 0;

        err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        for (PetscInt p = 0; p < closureSize*2; p += 2) {
            if ((closure[p] >= vStart) && (closure[p] < vEnd)) { numVertices++; }
        } // for
        err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        CPPUNIT_ASSERT_EQUAL(_data->numCorners, numVertices);
    } // for

    PYLITH_METHOD_END;
} // testSubmeshFault


// ----------------------------------------------------------------------
// Constructor
pylith::bc::TestBoundaryMesh_Data::TestBoundaryMesh_Data(void) :
    filename(NULL),
    bcLabel(NULL),
    faultLabel(NULL),
    faultId(0),
    numCorners(0),
    numCells(0),
    numVerticesNoFault(0),
    numVerticesWithFault(0) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::bc::TestBoundaryMesh_Data::~TestBoundaryMesh_Data(void) { // destructor
} // destructor


// End of file
