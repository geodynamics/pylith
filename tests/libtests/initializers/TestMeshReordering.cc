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

#include "pylith/initializers/MeshReordering.hh" // USES MeshReordering

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/problems/Problem.hh" // USES Problem
#include "pylith/utils/error.hh" // USES PylithCallPetscRequire()

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cmath> // USES fabs()
#include <valarray> // USES std::valarray

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace initializers {
        class TestMeshReordering;
    } // initializers
} // pylith

class pylith::initializers::TestMeshReordering {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test constructors.
    static
    void testConstructors(void);

    /// Test run().
    static
    void testRun(void);

}; // class TestMeshReordering

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshReordering::testConstructors", "[TestMeshReordering]") {
    pylith::initializers::TestMeshReordering::testConstructors();
}
TEST_CASE("TestMeshReordering::testRun", "[TestMeshReordering]") {
    pylith::initializers::TestMeshReordering::testRun();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::initializers::TestMeshReordering::testConstructors(void) {
    MeshReordering initializer;
} // testConstructors


// ------------------------------------------------------------------------------------------------
// Test run().
void
pylith::initializers::TestMeshReordering::testRun(void) {
    const std::string filename = "data/quad_small.mesh";
    const size_t dim = 2;
    const size_t numCells = 12;
    const size_t numVertices = 17;

    pylith::meshio::MeshIOAscii reader;
    reader.setFilename(filename.c_str());

    pylith::topology::Mesh meshOrig;
    reader.read(&meshOrig);

    // Verify reduction in Jacobian bandwidth
    pylith::topology::Field fieldOrig(meshOrig);
    pylith::topology::Field::Description description;
    description.label = "solution";
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = "field";
    description.scale = 1.0;
    description.validator = NULL;

    pylith::topology::Field::Discretization discretization;
    discretization.basisOrder = 1;
    discretization.quadOrder = 1;
    fieldOrig.subfieldAdd(description, discretization);
    fieldOrig.subfieldsSetup();
    fieldOrig.createDiscretization();
    fieldOrig.allocate();
    PetscMat matrix = NULL;
    PetscInt bandwidthOrig = 0;
    PylithCallPetscRequire(DMCreateMatrix(fieldOrig.getDM(), &matrix));
    PylithCallPetscRequire(MatComputeBandwidth(matrix, 0.0, &bandwidthOrig));
    PylithCallPetscRequire(MatDestroy(&matrix));

    MeshReordering initializer;
    pylith::problems::Problem problem;
    pylith::topology::Mesh* meshNew = initializer.run(&meshOrig, problem);

    REQUIRE(meshNew);
    CHECK(dim == meshNew->getDimension());
    CHECK(numCells == pylith::topology::MeshOps::getNumCells(*meshNew));
    CHECK(numVertices == pylith::topology::MeshOps::getNumVertices(*meshNew));

    pylith::topology::Field fieldNew(*meshNew);
    fieldNew.subfieldAdd(description, discretization);
    fieldNew.subfieldsSetup();
    fieldNew.createDiscretization();
    fieldNew.allocate();
    PetscInt bandwidth = 0;
    PylithCallPetscRequire(DMCreateMatrix(fieldNew.getDM(), &matrix));
    PylithCallPetscRequire(MatComputeBandwidth(matrix, 0.0, &bandwidth));
    PylithCallPetscRequire(MatDestroy(&matrix));

    CHECK(bandwidthOrig > 0);
    CHECK(bandwidth > 0);
    // CHECK(bandwidth <= bandwidthOrig);

    delete meshNew;meshNew = nullptr;
} // testRun


// End of file
