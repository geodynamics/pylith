// =================================================================================================
// This code is part of SpatialData, developed through the Computatreadernal Infrastructure
// for Geodynamics (https://github.com/geodynamics/spatialdata).
//
// Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license informatreadern.
// =================================================================================================

#include <portinfo>

#include "pylith/initializers/MeshRefiner.hh" // USES MeshRefiner

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/RefineUniform.hh" // USES RefineUniform
#include "pylith/problems/Problem.hh" // USES Problem

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cmath> // USES fabs()
#include <valarray> // USES std::valarray

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace initializers {
        class TestMeshRefiner;
    } // initializers
} // pylith

class pylith::initializers::TestMeshRefiner {
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

}; // class TestMeshRefiner

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshRefiner::testConstructors", "[TestMeshRefiner]") {
    pylith::initializers::TestMeshRefiner::testConstructors();
}
TEST_CASE("TestMeshRefiner::testAccessors", "[TestMeshRefiner]") {
    pylith::initializers::TestMeshRefiner::testAccessors();
}
TEST_CASE("TestMeshRefiner::testRun", "[TestMeshRefiner]") {
    pylith::initializers::TestMeshRefiner::testRun();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::initializers::TestMeshRefiner::testConstructors(void) {
    MeshRefiner initializer;
    CHECK(!initializer._refiner);
} // testConstructors


// ------------------------------------------------------------------------------------------------
// Test accessors.
void
pylith::initializers::TestMeshRefiner::testAccessors(void) {
    const size_t numLevels = 2;

    pylith::topology::RefineUniform refiner;
    refiner.setNumLevels(numLevels);

    MeshRefiner initializer;
    initializer.setRefiner(&refiner);
    CHECK(&refiner == initializer._refiner);
    CHECK(numLevels == refiner.getNumLevels());
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test run().
void
pylith::initializers::TestMeshRefiner::testRun(void) {
    const std::string filename = "data/quad_small.mesh";
    const size_t numLevels = 2;
    const size_t dim = 2;
    const size_t numCells = 192;
    const size_t numVertices = 209;

    pylith::topology::Mesh meshOrig;
    pylith::meshio::MeshIOAscii reader;
    reader.setFilename(filename.c_str());
    reader.read(&meshOrig);

    pylith::topology::RefineUniform refiner;
    refiner.setNumLevels(numLevels);

    MeshRefiner initializer;
    initializer.setRefiner(&refiner);

    pylith::problems::Problem problem;
    pylith::topology::Mesh* meshNew = initializer.run(&meshOrig, problem);

    CHECK(meshNew);
    CHECK(dim == meshNew->getDimension());
    CHECK(numCells == pylith::topology::MeshOps::getNumCells(*meshNew));
    CHECK(numVertices == pylith::topology::MeshOps::getNumVertices(*meshNew));

    delete meshNew;meshNew = nullptr;
} // testRun


// End of file
