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

#include "TestMeshIOAscii.hh" // Implementation of class methods

#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestMeshIOAscii_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestMeshIOAscii_Cases {
public:

    // Data factory methods
    static TestMeshIO_Data* Quad2D(void);

    static TestMeshIO_Data* Quad2D_Comments(void);

    static TestMeshIO_Data* Hex3D(void);

    static TestMeshIO_Data* Hex3D_Index1(void);

    static TestMeshIO_Data* Tri_OrphanVertex(void);

    static TestMeshIO_Data* Hex_OrphanVertex(void);

}; // TestMeshIOAscii_Cases

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshIOAscii::Quad2D::testFilename", "[TestMeshIOAscii][testFilename]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Quad2D()).testFilename();
}

TEST_CASE("TestMeshIOAscii::Quad2D::testWriteRead", "[TestMeshIOAscii][testWriteRead]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Quad2D()).testWriteRead();
}
TEST_CASE("TestMeshIOAscii::Hex3D::testWriteRead", "[TestMeshIOAscii][testWriteRead]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Hex3D()).testWriteRead();
}

TEST_CASE("TestMeshIOAscii::Quad2D_Comments::testRead", "[TestMeshIOAscii][testRead]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Quad2D_Comments()).testRead();
}
TEST_CASE("TestMeshIOAscii::Hex3D_Index1::testRead", "[TestMeshIOAscii][testRead]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Hex3D_Index1()).testRead();
}

TEST_CASE("TestMeshIOAscii::Tri_OrphanVertex::testReadError", "[TestMeshIOAscii][testReadError]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Tri_OrphanVertex()).testReadError();
}
TEST_CASE("TestMeshIOAscii::Hex_OrphanVertex::testReadError", "[TestMeshIOAscii][testReadError]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Hex_OrphanVertex()).testReadError();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii_Cases::Quad2D(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "quad.mesh";
    const size_t numVertices = 9;
    const size_t spaceDim = 2;
    const size_t cellDim = 2;
    const size_t numCells = 3;
    const size_t numCorners = 4;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::QUADRILATERAL;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -1.0, +3.0,
        +1.0, +3.3,
        -1.2, +0.9,
        +0.9, +1.0,
        +3.0, +2.9,
        +6.0, +1.2,
        +3.4, -0.2,
        +0.1, -1.1,
        +2.9, -3.1,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0,  2,  3,  1,
        4,  3,  6,  5,
        3,  7,  8,  6,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);
    static const PylithInt materialIds[numCells] = {
        1, 0, 1,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numVertexGroups = 2;
    static const PylithInt vertexGroupSizes[2] = { 5, 3, };
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[5+3] = {
        0, 2, 4, 5, 6,
        0, 1, 2,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[2] = {
        "vertex-group A",
        "vertex-group B",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);

    data->numFaceGroups = 2;
    data->numFaceVertices = 2;
    static const PylithInt faceGroupSizes[2] = { 3, 2, };
    data->faceGroupSizes = const_cast<PylithInt*>(faceGroupSizes);
    static const PylithInt faceGroups[3*(1+2)+2*(1+2)] = {
        0,  0, 2,
        1,  6, 5,
        1,  5, 4,

        0,  0, 2,
        0,  1, 0,
    };
    data->faceGroups = const_cast<PylithInt*>(faceGroups);
    static const char* faceGroupNames[2] = {
        "face-group A",
        "face-group B",
    };
    data->faceGroupNames = const_cast<char**>(faceGroupNames);

    return data;
} // Quad2D


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii_Cases::Quad2D_Comments(void) {
    TestMeshIO_Data* data = Quad2D();assert(data);

    data->filename = "data/quad_comments.mesh";

    return data;
} // Quad2D_Comments


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii_Cases::Hex3D(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "hex.mesh";
    const size_t numVertices = 14;
    const size_t spaceDim = 3;
    const size_t cellDim = 3;
    const size_t numCells = 2;
    const size_t numCorners = 8;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::HEXAHEDRON;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -3.0, -1.0, +0.2,
        -3.0, -1.0, +1.3,
        -1.0, -1.2, +0.1,
        -1.0, -1.2, +1.2,
        -3.0, +5.0, +1.3,
        -3.0, +5.0, +0.1,
        -0.5, +4.8, +0.2,
        -0.5, +4.8, +1.4,
        +0.5, +7.0, +1.2,
        +1.0, +3.1, +1.3,
        +3.0, +4.1, +1.4,
        +0.5, +7.0, -0.1,
        +1.0, +3.0, -0.2,
        +3.0, +4.2, +0.1
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        6, 12, 13, 11,  7,  9, 10,  8,
        0,  2,  6,  5,  1,  3,  7,  4
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);
    static const PylithInt materialIds[numCells] = {
        1, 0,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numVertexGroups = 3;
    static const PylithInt vertexGroupSizes[3] = { 4, 4, 7,};
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[4+4+7] = {
        0,  1,  2,  3,
        4,  5,  6,  7,
        1,  3,  4,  7,  8,  9, 10,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[3] = {
        "vertex-group A",
        "vertex-group B",
        "vertex-group C",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);

    data->numFaceGroups = 3;
    data->numFaceVertices = 4;
    static const PylithInt faceGroupSizes[3] = { 1, 1, 2,};
    data->faceGroupSizes = const_cast<PylithInt*>(faceGroupSizes);
    static const PylithInt faceGroups[(1+1+2)*(1+4)] = {
        1,   0,  2,  3,  1,

        1,   6,  5,  4,  7,

        0,   7,  9, 10,  8,
        1,   1,  3,  7,  4,
    };
    data->faceGroups = const_cast<PylithInt*>(faceGroups);
    static const char* faceGroupNames[3] = {
        "face-group A",
        "face-group B",
        "face-group C",
    };
    data->faceGroupNames = const_cast<char**>(faceGroupNames);

    return data;
} // Hex3D


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii_Cases::Hex3D_Index1(void) {
    TestMeshIO_Data* data = Hex3D();assert(data);

    data->filename = "data/hex_index1.mesh";

    return data;
} // Hex3D_Index1


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii_Cases::Tri_OrphanVertex(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/tri_orphan_vertex.mesh";

    return data;
} // Hex3D_Index1


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii_Cases::Hex_OrphanVertex(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/hex_orphan_vertex.mesh";

    return data;
} // Hex3D_Index1


// End of file
