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

#include "pylith/initializers/Initializer.hh" // USES Initializer

#include "pylith/initializers/MeshReader.hh" // USES MeshReader
#include "pylith/initializers/MeshReordering.hh" // USES MeshReordering
#include "pylith/initializers/MeshRefiner.hh" // USES MeshRefiner

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/RefineUniform.hh" // USES RefineUniform
#include "pylith/problems/Problem.hh" // USES Problem
#include "pylith/scales/Scales.hh" // USES Scales

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cmath> // USES fabs()
#include <valarray> // USES std::valarray

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace initializers {
        class TestInitializer;
    } // initializers
} // pylith

class pylith::initializers::TestInitializer {
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

}; // class TestInitializer

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestInitializer::testConstructors", "[TestInitializer]") {
    pylith::initializers::TestInitializer::testConstructors();
}
TEST_CASE("TestInitializer::testAccessors", "[TestInitializer]") {
    pylith::initializers::TestInitializer::testAccessors();
}
TEST_CASE("TestInitializer::testRun", "[TestInitializer]") {
    pylith::initializers::TestInitializer::testRun();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::initializers::TestInitializer::testConstructors(void) {
    Initializer initializer;
    CHECK(!initializer._phases.size());
} // testConstructors


// ------------------------------------------------------------------------------------------------
// Test accessors.
void
pylith::initializers::TestInitializer::testAccessors(void) {
    pylith::initializers::MeshReader reader;
    pylith::initializers::MeshReordering reordering;
    pylith::initializers::MeshRefiner refiner;

    const size_t numPhases = 3;
    pylith::initializers::InitializePhase* phases[numPhases] = {
        &reader,
        &reordering,
        &refiner,
    };

    Initializer initializer;
    initializer.setPhases(phases, numPhases);
    REQUIRE(3 == initializer._phases.size());
    REQUIRE(&reader == initializer._phases[0]);
    REQUIRE(&reordering == initializer._phases[1]);
    REQUIRE(&refiner == initializer._phases[2]);
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test run().
void
pylith::initializers::TestInitializer::testRun(void) {
    const std::string filename = "data/quad_small.mesh";
    const size_t numLevels = 2;
    const size_t dim = 2;
    const size_t numCells = 192;
    const size_t numVertices = 209;

    pylith::initializers::MeshReader phaseRead;
    pylith::topology::Mesh meshOrig;
    pylith::meshio::MeshIOAscii reader;
    reader.setFilename(filename.c_str());
    phaseRead.setReader(&reader);

    pylith::initializers::MeshReordering phaseReorder;

    pylith::initializers::MeshRefiner phaseRefine;
    pylith::topology::RefineUniform refiner;
    refiner.setNumLevels(numLevels);
    phaseRefine.setRefiner(&refiner);

    const size_t numPhases = 3;
    pylith::initializers::InitializePhase* phases[numPhases] = {
        &phaseRead,
        &phaseReorder,
        &phaseRefine,
    };

    Initializer initializer;
    initializer.setPhases(phases, numPhases);

    pylith::problems::Problem problem;
    pylith::scales::Scales scales;
    problem.setScales(scales);
    pylith::topology::Mesh* meshNew = initializer.runPhases(problem);

    CHECK(meshNew);
    CHECK(dim == meshNew->getDimension());
    CHECK(numCells == pylith::topology::MeshOps::getNumCells(*meshNew));
    CHECK(numVertices == pylith::topology::MeshOps::getNumVertices(*meshNew));

    delete meshNew;meshNew = nullptr;
} // testRun


// End of file
