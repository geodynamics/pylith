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

#include "pylith/initializers/MeshReader.hh" // USES MeshReader

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/problems/Problem.hh" // USES Problem
#include "pylith/scales/Scales.hh" // USES Scales

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cmath> // USES fabs()
#include <valarray> // USES std::valarray

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace initializers {
        class TestMeshReader;
    } // initializers
} // pylith

class pylith::initializers::TestMeshReader {
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

}; // class TestMeshReader

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshReader::testConstructors", "[TestMeshReader]") {
    pylith::initializers::TestMeshReader::testConstructors();
}
TEST_CASE("TestMeshReader::testAccessors", "[TestMeshReader]") {
    pylith::initializers::TestMeshReader::testAccessors();
}
TEST_CASE("TestMeshReader::testRun", "[TestMeshReader]") {
    pylith::initializers::TestMeshReader::testRun();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::initializers::TestMeshReader::testConstructors(void) {
    MeshReader initializer;
    CHECK(!initializer._reader);
} // testConstructors


// ------------------------------------------------------------------------------------------------
// Test accessors.
void
pylith::initializers::TestMeshReader::testAccessors(void) {
    const std::string filename = "mesh.txt";

    pylith::meshio::MeshIOAscii reader;
    reader.setFilename(filename.c_str());

    MeshReader initializer;
    initializer.setReader(&reader);
    REQUIRE(&reader == initializer._reader);
    REQUIRE(filename == std::string(reader.getFilename()));
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test run().
void
pylith::initializers::TestMeshReader::testRun(void) {
    const std::string filename = "data/quad_small.mesh";
    const size_t dim = 2;
    const size_t numCells = 12;
    const size_t numVertices = 17;

    pylith::meshio::MeshIOAscii reader;
    reader.setFilename(filename.c_str());

    MeshReader initializer;
    initializer.setReader(&reader);

    pylith::topology::Mesh* nullMesh = nullptr;
    pylith::problems::Problem problem;
    pylith::scales::Scales scales;
    problem.setScales(scales);
    pylith::topology::Mesh* mesh = initializer.run(nullMesh, problem);

    REQUIRE(mesh);
    CHECK(dim == mesh->getDimension());
    CHECK(numCells == pylith::topology::MeshOps::getNumCells(*mesh));
    CHECK(numVertices == pylith::topology::MeshOps::getNumVertices(*mesh));

    delete mesh;mesh = nullptr;
} // testRun


// End of file
