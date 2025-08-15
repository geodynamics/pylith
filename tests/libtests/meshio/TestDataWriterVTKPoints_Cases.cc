// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestDataWriterVTKPoints.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestDataWriterVTKPoints_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTKPoints_Cases {
public:

    // Data factory methods
    static TestDataWriterVTKPoints_Data* Tri(void);

    static TestDataWriterVTKPoints_Data* Quad(void);

    static TestDataWriterVTKPoints_Data* Tet(void);

    static TestDataWriterVTKPoints_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
#if 0
TEST_CASE("TestDataWriterVTKMesh::Tri::testTimeStep", "[DataWriter][VTK][Points][Tri][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKPoints(pylith::meshio::TestDataWriterVTKPoints_Cases::Tri()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKPoints::Tri::testWriteVertexField", "[DataWriter][VTK][Points][Tri][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKPoints(pylith::meshio::TestDataWriterVTKPoints_Cases::Tri()).testWriteVertexField();
}

TEST_CASE("TestDataWriterVTKPoints::Quad::testTimeStep", "[DataWriter][VTK][Points][Quad][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKPoints(pylith::meshio::TestDataWriterVTKPoints_Cases::Quad()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKPoints::Quad::testWriteVertexField", "[DataWriter][VTK][Points][Quad][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKPoints(pylith::meshio::TestDataWriterVTKPoints_Cases::Quad()).testWriteVertexField();
}

TEST_CASE("TestDataWriterVTKPoints::Tet::testTimeStep", "[DataWriter][VTK][Points][Tet][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKPoints(pylith::meshio::TestDataWriterVTKPoints_Cases::Tet()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKPoints::Tet::testWriteVertexField", "[DataWriter][VTK][Points][Tet][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKPoints(pylith::meshio::TestDataWriterVTKPoints_Cases::Tet()).testWriteVertexField();
}

TEST_CASE("TestDataWriterVTKPoints::Hex::testTimeStep", "[DataWriter][VTK][Points][Hex][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKPoints(pylith::meshio::TestDataWriterVTKPoints_Cases::Hex()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKPoints::Hex::testWriteVertexField", "[DataWriter][VTK][Points][Hex][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKPoints(pylith::meshio::TestDataWriterVTKPoints_Cases::Hex()).testWriteVertexField();
}
#endif

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKPoints_Data*
pylith::meshio::TestDataWriterVTKPoints_Cases::Tri(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKPoints_Data* data = new TestDataWriterVTKPoints_Data();assert(data);

    data->timestepFilename = "tri3_points.vtu";
    data->vertexFilename = "tri3_points_vertex.vtu";

    TestDataWriterPoints::setDataTri(data);

    PYLITH_METHOD_RETURN(data);
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKPoints_Data*
pylith::meshio::TestDataWriterVTKPoints_Cases::Quad(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKPoints_Data* data = new TestDataWriterVTKPoints_Data();assert(data);

    data->timestepFilename = "quad4_points.vtu";
    data->vertexFilename = "quad4_points_vertex.vtu";

    TestDataWriterPoints::setDataQuad(data);

    PYLITH_METHOD_RETURN(data);
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKPoints_Data*
pylith::meshio::TestDataWriterVTKPoints_Cases::Tet(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKPoints_Data* data = new TestDataWriterVTKPoints_Data();assert(data);

    data->timestepFilename = "tet4_points.vtu";
    data->vertexFilename = "tet4_points_vertex.vtu";

    TestDataWriterPoints::setDataTet(data);

    PYLITH_METHOD_RETURN(data);
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKPoints_Data*
pylith::meshio::TestDataWriterVTKPoints_Cases::Hex(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKPoints_Data* data = new TestDataWriterVTKPoints_Data();assert(data);

    data->timestepFilename = "hex8_points.vtu";
    data->vertexFilename = "hex8_points_vertex.vtu";

    TestDataWriterPoints::setDataHex(data);

    PYLITH_METHOD_RETURN(data);
} // Hex


// End of file
