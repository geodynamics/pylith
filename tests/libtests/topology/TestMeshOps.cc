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

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/MeshOps.hh" // USES MeshOps

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <stdexcept> // USES std::runtime_error

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "catch2/matchers/catch_matchers_exception.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace topology {
        class TestMeshOps;
    } // topology
} // pylith

class pylith::topology::TestMeshOps : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test nondimensionalize().
    static
    void testNondimensionalize(void);

    /// Test checkTopology().
    static
    void testCheckTopology(void);

    /// Test isSimplexMesh().
    static
    void testIsSimplexMesh(void);

    /// Test checkMaterialIds().
    static
    void testCheckMaterialIds(void);

}; // class TestMeshOps

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshOps::testNondimensionalize", "[TestMeshOps]") {
    pylith::topology::TestMeshOps().testNondimensionalize();
}
TEST_CASE("TestMeshOps::testCheckTopology", "[TestMeshOps]") {
    pylith::topology::TestMeshOps().testCheckTopology();
}
TEST_CASE("TestMeshOps::testIsSimplexMesh", "[TestMeshOps]") {
    pylith::topology::TestMeshOps().testIsSimplexMesh();
}
TEST_CASE("TestMeshOps::testCheckMaterialIds", "[TestMeshOps]") {
    pylith::topology::TestMeshOps().testCheckMaterialIds();
}

// ------------------------------------------------------------------------------------------------
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
    iohandler.setFilename("data/tri3.mesh");
    iohandler.read(&mesh);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(2);
    mesh.setCoordSys(&cs);
    spatialdata::units::Nondimensional normalizer;
    normalizer.setLengthScale(lengthScale);
    MeshOps::nondimensionalize(&mesh, normalizer);

    // Get vertices
    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    // Check nondimensional coordinates
    CoordsVisitor coordsVisitor(dmMesh);
    const PetscScalar* coordsArray = coordsVisitor.localArray();assert(coordsArray);

    const PylithScalar tolerance = 1.0e-06;
    for (PetscInt v = vStart, i = 0; v < vEnd; ++v) {
        REQUIRE(spaceDim == coordsVisitor.sectionDof(v));
        const PetscInt off = coordsVisitor.sectionOffset(v);
        for (int iDim = 0; iDim < spaceDim; ++iDim, ++i) {
            const PylithScalar coordE = coordinates[i] / lengthScale;
            const double toleranceV = std::max(tolerance, tolerance*coordE);
            CHECK_THAT(coordsArray[off+iDim], Catch::Matchers::WithinAbs(coordE, toleranceV));
        } // for
    } // for

    PYLITH_METHOD_END;
} // testNondimensionalize


// ------------------------------------------------------------------------------------------------
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
        iohandler.setFilename(filename);
        iohandler.read(&mesh);
        MeshOps::checkTopology(mesh);
    } // for

    PYLITH_METHOD_END;
} // testCheckTopology


// ------------------------------------------------------------------------------------------------
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
        iohandler.setFilename(filename);
        iohandler.read(&mesh);
        CHECK(isSimplex == MeshOps::isSimplexMesh(mesh));
    } // for

    PYLITH_METHOD_END;
} // testIsSimplexMesh


// ------------------------------------------------------------------------------------------------
// Test checkMaterialIds().
void
pylith::topology::TestMeshOps::testCheckMaterialIds(void) {
    Mesh mesh;
    meshio::MeshIOAscii iohandler;
    iohandler.setFilename("data/tri3.mesh");
    iohandler.read(&mesh);

    pylith::int_array materialValues(2);
    materialValues[0] = 4;
    materialValues[1] = 3;

    MeshOps::checkMaterialLabels(mesh, materialValues);

    materialValues[0] = 99;
    REQUIRE_THROWS_AS(MeshOps::checkMaterialLabels(mesh, materialValues), std::runtime_error);
} // testCheckMaterialIds


// End of file
