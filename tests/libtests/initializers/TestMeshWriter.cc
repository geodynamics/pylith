// =================================================================================================
// This code is part of SpatialData, developed through the Computatreadernal Infrastructure
// for Geodynamics (https://github.com/geodynamics/spatialdata).
//
// Copyright (c) 2010-2026, University of California, Davis and the SpatialData Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license informatreadern.
// =================================================================================================

#include <portinfo>

#include "pylith/initializers/MeshWriter.hh" // USES MeshWriter

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/problems/Problem.hh" // USES Problem

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cmath> // USES fabs()
#include <valarray> // USES std::valarray

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace initializers {
        class TestMeshWriter;
    } // initializers
} // pylith

class pylith::initializers::TestMeshWriter {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test constructors.
    static
    void testConstructors(void);

    /// Test accessors.
    static
    void testAccessors(void);

    /// Test run().
    static
    void testRun(void);

}; // class TestMeshWriter

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshWriter::testConstructors", "[TestMeshWriter]") {
    pylith::initializers::TestMeshWriter::testConstructors();
}
TEST_CASE("TestMeshWriter::testAccessors", "[TestMeshWriter]") {
    pylith::initializers::TestMeshWriter::testAccessors();
}
TEST_CASE("TestMeshWriter::testRun", "[TestMeshWriter]") {
    pylith::initializers::TestMeshWriter::testRun();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::initializers::TestMeshWriter::testConstructors(void) {
    MeshWriter initializer;
    CHECK(!initializer._writer);
} // testConstructors


// ------------------------------------------------------------------------------------------------
// Test accessors.
void
pylith::initializers::TestMeshWriter::testAccessors(void) {
    const std::string filename = "mesh.txt";

    pylith::meshio::MeshIOAscii writer;
    writer.setFilename(filename.c_str());

    MeshWriter initializer;
    initializer.setWriter(&writer);
    REQUIRE(&writer == initializer._writer);
    REQUIRE(filename == std::string(writer.getFilename()));
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test run().
void
pylith::initializers::TestMeshWriter::testRun(void) {
    const std::string filenameIn = "data/quad_small.mesh";
    const std::string filenameOut = "quad_small.mesh";
    const size_t dim = 2;
    const size_t numCells = 12;
    const size_t numVertices = 17;

    pylith::topology::Mesh mesh;
    pylith::meshio::MeshIOAscii io;
    io.setFilename(filenameIn.c_str());
    io.read(&mesh);

    io.setFilename(filenameOut.c_str());
    MeshWriter initializer;
    initializer.setWriter(&io);

    pylith::problems::Problem problem;
    pylith::topology::Mesh* meshNew = initializer.run(&mesh, problem);
    CHECK(meshNew);
    CHECK(dim == meshNew->getDimension());
    CHECK(numCells == pylith::topology::MeshOps::getNumCells(*meshNew));
    CHECK(numVertices == pylith::topology::MeshOps::getNumVertices(*meshNew));
    delete meshNew;meshNew = nullptr;

    io.read(&mesh);
    CHECK(dim == mesh.getDimension());
    CHECK(numCells == pylith::topology::MeshOps::getNumCells(mesh));
    CHECK(numVertices == pylith::topology::MeshOps::getNumVertices(mesh));
} // testRun


// End of file
