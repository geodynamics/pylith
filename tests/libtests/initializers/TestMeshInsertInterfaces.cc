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

#include "pylith/initializers/MeshInsertInterfaces.hh" // USES MeshInsertInterfaces

#include "tests/src/FaultCohesiveStub.hh" // USES FaultsCohesiveStub
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
        class TestMeshInsertInterfaces;
    } // initializers
} // pylith

class pylith::initializers::TestMeshInsertInterfaces {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test constructors.
    static
    void testConstructors(void);

    /// Test run().
    static
    void testRun(void);

}; // class TestMeshInsertInterfaces

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshInsertInterfaces::testConstructors", "[TestMeshInsertInterfaces]") {
    pylith::initializers::TestMeshInsertInterfaces::testConstructors();
}
TEST_CASE("TestMeshInsertInterfaces::testRun", "[TestMeshInsertInterfaces]") {
    pylith::initializers::TestMeshInsertInterfaces::testRun();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::initializers::TestMeshInsertInterfaces::testConstructors(void) {
    MeshInsertInterfaces initializer;
} // testConstructors


// ------------------------------------------------------------------------------------------------
// Test run().
void
pylith::initializers::TestMeshInsertInterfaces::testRun(void) {
    const std::string filename = "data/quad_small.mesh";
    const std::string faultGroup = "fault";

    const size_t dim = 2;
    const size_t numCells = 12+4;
    const size_t numVertices = 17+5;

    pylith::topology::Mesh meshOrig;
    pylith::meshio::MeshIOAscii reader;
    reader.setFilename(filename.c_str());
    reader.read(&meshOrig);

    pylith::faults::FaultCohesiveStub fault;
    pylith::faults::FaultCohesive* faults[1] = { &fault };
    fault.setCohesiveLabelName(pylith::topology::Mesh::cells_label_name);
    fault.setSurfaceLabelName(faultGroup.c_str());

    pylith::problems::Problem problem;
    problem.setInterfaces(faults, 1);

    pylith::initializers::MeshInsertInterfaces initializer;
    pylith::topology::Mesh* meshNew = initializer.run(&meshOrig, problem);

    CHECK(meshNew);
    CHECK(dim == meshNew->getDimension());
    CHECK(numCells == pylith::topology::MeshOps::getNumCells(*meshNew));
    CHECK(numVertices == pylith::topology::MeshOps::getNumVertices(*meshNew));

    delete meshNew;meshNew = nullptr;
} // testRun


// End of file
