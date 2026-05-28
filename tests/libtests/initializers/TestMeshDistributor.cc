// =================================================================================================
// This code is part of SpatialData, developed through the Computatdistributornal Infrastructure
// for Geodynamics (https://github.com/geodynamics/spatialdata).
//
// Copyright (c) 2010-2026, University of California, Davis and the SpatialData Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license informatdistributorn.
// =================================================================================================

#include <portinfo>

#include "pylith/initializers/MeshDistributor.hh" // USES MeshDistributor

#include "pylith/topology/Distributor.hh" // USES Distributor
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
        class TestMeshDistributor;
    } // initializers
} // pylith

class pylith::initializers::TestMeshDistributor {
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

}; // class TestMeshDistributor

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshDistributor::testConstructors", "[TestMeshDistributor]") {
    pylith::initializers::TestMeshDistributor::testConstructors();
}
TEST_CASE("TestMeshDistributor::testAccessors", "[TestMeshDistributor]") {
    pylith::initializers::TestMeshDistributor::testAccessors();
}
TEST_CASE("TestMeshDistributor::testRun", "[TestMeshDistributor]") {
    pylith::initializers::TestMeshDistributor::testRun();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::initializers::TestMeshDistributor::testConstructors(void) {
    MeshDistributor initializer;
    CHECK(!initializer._distributor);
} // testConstructors


// ------------------------------------------------------------------------------------------------
// Test accessors.
void
pylith::initializers::TestMeshDistributor::testAccessors(void) {
    pylith::topology::Distributor distributor;

    MeshDistributor initializer;
    initializer.setDistributor(&distributor);
    CHECK(&distributor == initializer._distributor);
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test run().
void
pylith::initializers::TestMeshDistributor::testRun(void) {
    const std::string partitioner = "parmetis";
    const std::string filename = "data/quad_small.mesh";
    const size_t dim = 2;
    const size_t numCells = 12;
    const size_t numVertices = 17;

    pylith::topology::Mesh meshOrig;
    pylith::meshio::MeshIOAscii reader;
    reader.setFilename(filename.c_str());
    reader.read(&meshOrig);

    pylith::topology::Distributor distributor;
    distributor.setPartitioner(partitioner.c_str());

    MeshDistributor initializer;
    initializer.setDistributor(&distributor);

    pylith::problems::Problem problem;
    pylith::topology::Mesh* meshNew = initializer.run(&meshOrig, problem);

    REQUIRE(meshNew);
    CHECK(dim == meshNew->getDimension());
    CHECK(numCells == pylith::topology::MeshOps::getNumCells(*meshNew));
    CHECK(numVertices == pylith::topology::MeshOps::getNumVertices(*meshNew));

    delete meshNew;meshNew = nullptr;
} // testRun


// End of file
