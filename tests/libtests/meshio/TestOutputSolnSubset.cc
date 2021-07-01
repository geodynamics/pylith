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

#include "TestOutputSolnSubset.hh" // Implementation of class methods

#include "pylith/meshio/OutputSolnSubset.hh"

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include <string.h> // USES strcmp()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::meshio::TestOutputSolnSubset);

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestOutputSolnSubset::testConstructor(void) { // testConstructor
    PYLITH_METHOD_BEGIN;

    OutputSolnSubset output;

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test getLabel()
void
pylith::meshio::TestOutputSolnSubset::testLabel(void) { // testLabel
    PYLITH_METHOD_BEGIN;

    OutputSolnSubset output;

    const char* label = "boundary";

    output.setLabel(label);
    CPPUNIT_ASSERT(0 == strcmp(label, output._label.c_str()));

    PYLITH_METHOD_END;
} // testLabel


// ----------------------------------------------------------------------
// Test subdomainMesh()
void
pylith::meshio::TestOutputSolnSubset::testSubdomainMesh(void) { // testSubdomainMesh
    PYLITH_METHOD_BEGIN;

    const char* label = "bc3";
    const int nvertices = 3;
    const int verticesE[nvertices] = { 2, 3, 4 };
    const int ncells = 2;
    const int ncorners = 2;
    const int cellsE[ncells*ncorners] = {
        2, 3,
        3, 4,
    };

    topology::Mesh mesh;
    MeshIOAscii iohandler;
    iohandler.filename("data/quad4.mesh");
    iohandler.read(&mesh);

    OutputSolnSubset output;
    output.setLabel(label);

    PetscDM dmMesh = output.subdomainMesh(mesh).getDM();CPPUNIT_ASSERT(dmMesh);

    // Check vertices
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();
    CPPUNIT_ASSERT_EQUAL(nvertices, verticesStratum.size());
    for (PetscInt v = vStart, index = 0; v < vEnd; ++v, ++index) {
        CPPUNIT_ASSERT_EQUAL(verticesE[index], v);
    } // for

    // Check cells
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 1);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    PetscErrorCode err = 0;
    CPPUNIT_ASSERT_EQUAL(ncells, cellsStratum.size());
    for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
        PetscInt *closure = NULL;
        PetscInt closureSize = 0;
        err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        int count = 0;
        for (int i = 0; i < closureSize; ++i) {
            const PetscInt p = closure[2*i];
            if (( p >= vStart) && ( p < vEnd) ) {
                CPPUNIT_ASSERT_EQUAL(cellsE[index], p);
                ++count;
                ++index;
            } // if
        } // for
        CPPUNIT_ASSERT_EQUAL(ncorners, count);
        err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // testSubdomainMesh


// End of file
