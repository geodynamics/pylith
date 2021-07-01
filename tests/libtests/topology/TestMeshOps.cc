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

#include "TestMeshOps.hh" // Implementation of class methods

#include "pylith/topology/MeshOps.hh" // USES MeshOps

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestMeshOps);

// ---------------------------------------------------------------------------------------------------------------------
// Test nondimensionalize().
void
pylith::topology::TestMeshOps::testNondimensionalize(void) {
    PYLITH_METHOD_BEGIN;

    const PylithScalar lengthScale = 2.0;
    const int spaceDim = 2;
    const int numVertices = 4;
    const PylithScalar coordinates[numVertices*spaceDim] = {
        -1.0, 0.0,
        0.0, -1.0,
        0.0, 1.0,
        1.0, 0.0,
    };

    Mesh mesh;
    meshio::MeshIOAscii iohandler;
    iohandler.filename("data/tri3.mesh");
    iohandler.read(&mesh);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(2);
    mesh.setCoordSys(&cs);
    spatialdata::units::Nondimensional normalizer;
    normalizer.setLengthScale(lengthScale);
    MeshOps::nondimensionalize(&mesh, normalizer);

    // Get vertices
    PetscDM dmMesh = mesh.getDM();CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    // Check nondimensional coordinates
    CoordsVisitor coordsVisitor(dmMesh);
    const PetscScalar* coordsArray = coordsVisitor.localArray();CPPUNIT_ASSERT(coordsArray);

    const PylithScalar tolerance = 1.0e-06;
    for (PetscInt v = vStart, i = 0; v < vEnd; ++v) {
        CPPUNIT_ASSERT_EQUAL(spaceDim, coordsVisitor.sectionDof(v));
        const PetscInt off = coordsVisitor.sectionOffset(v);
        for (int iDim = 0; iDim < spaceDim; ++iDim, ++i) {
            const PylithScalar coordE = coordinates[i] / lengthScale;
            if (fabs(coordE) < 1.0) {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(coordE, coordsArray[off+iDim], tolerance);
            } else {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, coordsArray[off+iDim]/coordE, tolerance);
            } // if/else
        } // for
    } // for

    PYLITH_METHOD_END;
} // testNondimensionalize


// ---------------------------------------------------------------------------------------------------------------------
// Test checkTopology().
void
pylith::topology::TestMeshOps::testCheckTopology(void) {
    PYLITH_METHOD_BEGIN;

    const int numFiles = 4;
    const char* filenames[numFiles] = {
        "data/tri3.mesh",
        "data/fourquad4.mesh",
        "data/twotet4.mesh",
        "data/twohex8.mesh",
    };

    for (int i = 0; i < numFiles; ++i) {
        const char* filename = filenames[i];
        Mesh mesh;
        meshio::MeshIOAscii iohandler;
        iohandler.filename(filename);
        iohandler.read(&mesh);
        MeshOps::checkTopology(mesh);
    } // for

    PYLITH_METHOD_END;
} // testCheckTopology


// ---------------------------------------------------------------------------------------------------------------------
// Test isSimplexMesh().
void
pylith::topology::TestMeshOps::testIsSimplexMesh(void) {
    PYLITH_METHOD_BEGIN;

    const int numFiles = 4;
    const char* filenames[numFiles] = {
        "data/tri3.mesh",
        "data/fourquad4.mesh",
        "data/twotet4.mesh",
        "data/twohex8.mesh",
    };
    const bool results[numFiles] = {
        true,
        false,
        true,
        false,
    };

    for (int i = 0; i < numFiles; ++i) {
        const char* filename = filenames[i];
        const bool isSimplex = results[i];
        Mesh mesh;
        meshio::MeshIOAscii iohandler;
        iohandler.filename(filename);
        iohandler.read(&mesh);
        CPPUNIT_ASSERT_EQUAL(isSimplex, MeshOps::isSimplexMesh(mesh));
    } // for

    PYLITH_METHOD_END;
} // testIsSimplexMesh


// ---------------------------------------------------------------------------------------------------------------------
// Test checkMaterialIds().
void
pylith::topology::TestMeshOps::testCheckMaterialIds(void) {
    PYLITH_METHOD_BEGIN;

    Mesh mesh;
    meshio::MeshIOAscii iohandler;
    iohandler.filename("data/tri3.mesh");
    iohandler.read(&mesh);

    pylith::int_array materialIds(2);
    materialIds[0] = 4;
    materialIds[1] = 3;

    MeshOps::checkMaterialIds(mesh, materialIds);

    materialIds[0] = 99;
    CPPUNIT_ASSERT_THROW(MeshOps::checkMaterialIds(mesh, materialIds), std::runtime_error);

    PYLITH_METHOD_END;
} // testCheckMaterialIds


// End of file
