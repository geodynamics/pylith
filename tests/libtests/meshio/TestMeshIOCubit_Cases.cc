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

#include "TestMeshIOCubit.hh" // Implementation of class methods

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestMeshIOCubit_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestMeshIOCubit_Cases {
public:

    // Data factory methods
    static TestMeshIO_Data* Tri_v12(void);

    static TestMeshIO_Data* Tri_v13(void);

    static TestMeshIO_Data* Tri_v16(void);

    static TestMeshIO_Data* Quad_v12(void);

    static TestMeshIO_Data* Quad_v13(void);

    static TestMeshIO_Data* Quad_v16(void);

    static TestMeshIO_Data* Tet_v12(void);

    static TestMeshIO_Data* Tet_v13(void);

    static TestMeshIO_Data* Tet_v16(void);

    static TestMeshIO_Data* Hex_v12(void);

    static TestMeshIO_Data* Hex_v13(void);

    static TestMeshIO_Data* Hex_v16(void);

}; // TestMeshIOCubit_Cases

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshIOCubit::Tri::testFilename", "[TestMeshIOCubit][testFilename]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tri_v12()).testFilename();
}

TEST_CASE("TestMeshIOCubit::Tri_v12::testRead", "[TestMeshIOCubit][Tri][v12][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tri_v12()).testRead();
}
TEST_CASE("TestMeshIOCubit::Tri_v13::testRead", "[TestMeshIOCubit][Tri][v13][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tri_v13()).testRead();
}
TEST_CASE("TestMeshIOCubit::Tri_v16::testRead", "[TestMeshIOCubit][Tri][v16][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tri_v16()).testRead();
}

TEST_CASE("TestMeshIOCubit::Quad_v12::testRead", "[TestMeshIOCubit][Quad][v12][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Quad_v12()).testRead();
}
TEST_CASE("TestMeshIOCubit::Quad_v13::testRead", "[TestMeshIOCubit][Quad][v13][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Quad_v13()).testRead();
}
TEST_CASE("TestMeshIOCubit::Quad_v16::testRead", "[TestMeshIOCubit][Quad][v16][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Quad_v16()).testRead();
}

TEST_CASE("TestMeshIOCubit::Tet_v12::testRead", "[TestMeshIOCubit][Tet][v12][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tet_v12()).testRead();
}
TEST_CASE("TestMeshIOCubit::Tet_v13::testRead", "[TestMeshIOCubit][Tet][v13][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tet_v13()).testRead();
}
TEST_CASE("TestMeshIOCubit::Tet_v16::testRead", "[TestMeshIOCubit][Tet][v16][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tet_v16()).testRead();
}

TEST_CASE("TestMeshIOCubit::Hex_v12::testRead", "[TestMeshIOCubit][Hex][v12][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Hex_v12()).testRead();
}
TEST_CASE("TestMeshIOCubit::Hex_v13::testRead", "[TestMeshIOCubit][Hex][v13][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Hex_v13()).testRead();
}
TEST_CASE("TestMeshIOCubit::Hex_v16::testRead", "[TestMeshIOCubit][Hex][v16][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Hex_v16()).testRead();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Tri_v12(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/twotri3_12.2.exo";

    const size_t numVertices = 4;
    const size_t spaceDim = 2;
    const size_t cellDim = 2;
    const size_t numCells = 2;
    const size_t numCorners = 3;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::TRIANGLE;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -1.0,  +0.0,
        +0.0,  -1.0,
        +0.0,  +1.0,
        +1.0,  +0.0
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0,  1,  2,
        2,  1,  3,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);
    static const PylithInt materialIds[numCells] = {
        2, 3,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numVertexGroups = 2;
    static const PylithInt vertexGroupSizes[2] = { 1,  2, };
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[1+2] = {
        0,
        2, 3,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[2] = {
        "left_vertex",
        "right_vertex",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);

    return data;
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Tri_v13(void) {
    TestMeshIO_Data* data = Tri_v12();assert(data);

    data->filename = "data/twotri3_13.0.exo";

    return data;
} // Tri_v13


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Tri_v16(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/smalltri_v16.exo";

    const size_t numVertices = 6;
    const size_t spaceDim = 2;
    const size_t cellDim = 2;
    const size_t numCells = 4;
    const size_t numCorners = 3;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::TRIANGLE;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -1.0,  +0.0,
        +0.0,  +0.0,
        +0.0,  +1.0,
        -1.0,  +1.0,
        +1.0,   0.0,
        +1.0,   1.0,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0,  1,  2,
        0,  2,  3,
        1,  4,  5,
        1,  5,  2,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);
    static const PylithInt materialIds[numCells] = {
        1, 1,  2, 2,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numVertexGroups = 2;
    static const PylithInt vertexGroupSizes[2] = { 2,  3, };
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[2+3] = {
        0, 3,
        0, 1, 4,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[2] = {
        "vertices_xneg",
        "vertices_yneg",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);

    data->numFaceGroups = 2;
    data->numFaceVertices = 2;
    static const PylithInt faceGroupSizes[2] = { 1,  2, };
    data->faceGroupSizes = const_cast<PylithInt*>(faceGroupSizes);
    static const PylithInt faceGroups[3*(1+2)] = {
        1,   3,  0,
        0,   0,  1,
        2,   1,  4,
    };
    data->faceGroups = const_cast<PylithInt*>(faceGroups);
    static const char* faceGroupNames[2] = {
        "boundary_xneg",
        "boundary_yneg",
    };
    data->faceGroupNames = const_cast<char**>(faceGroupNames);

    return data;
} // Tri_v16


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Quad_v12(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/twoquad4_12.2.exo";

    const size_t numVertices = 6;
    const size_t spaceDim = 2;
    const size_t cellDim = 2;
    const size_t numCells = 2;
    const size_t numCorners = 4;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::QUADRILATERAL;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        0.0,  0.0,
        1.0,  0.0,
        1.0,  1.0,
        0.0,  1.0,
        2.0,  0.0,
        2.0,  1.0,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0,  1,  2,  3,
        1,  4,  5,  2,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);
    static const PylithInt materialIds[numCells] = {
        10, 11,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numVertexGroups = 2;
    static const PylithInt vertexGroupSizes[2] = { 2,  3, };
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[2+3] = {
        0, 3,
        2, 3, 5,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[2] = {
        "left_edge",
        "top_edge",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);

    return data;
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Quad_v13(void) {
    TestMeshIO_Data* data = Quad_v12();assert(data);

    data->filename = "data/twoquad4_13.0.exo";

    return data;
} // Quad_v13


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Quad_v16(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/smallquad_v16.exo";

    const size_t numVertices = 6;
    const size_t spaceDim = 2;
    const size_t cellDim = 2;
    const size_t numCells = 2;
    const size_t numCorners = 4;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::QUADRILATERAL;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -1.0,  +1.0,
        -1.0,  +0.0,
        +0.0,  +0.0,
        +0.0,  +1.0,
        +1.0,  +0.0,
        +1.0,  +1.0,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0,  1,  2,  3,
        3,  2,  4,  5,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);
    static const PylithInt materialIds[numCells] = {
        1, 2,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numVertexGroups = 2;
    static const PylithInt vertexGroupSizes[2] = { 2,  3, };
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[2+3] = {
        0, 1,
        1, 2, 4,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[2] = {
        "vertices_xneg",
        "vertices_yneg",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);

    data->numFaceGroups = 2;
    data->numFaceVertices = 2;
    static const PylithInt faceGroupSizes[2] = { 1,  2, };
    data->faceGroupSizes = const_cast<PylithInt*>(faceGroupSizes);
    static const PylithInt faceGroups[3*(1+2)] = {
        0,   0,  1,
        0,   1,  2,
        1,   2,  4,
    };
    data->faceGroups = const_cast<PylithInt*>(faceGroups);
    static const char* faceGroupNames[2] = {
        "boundary_xneg",
        "boundary_yneg",
    };
    data->faceGroupNames = const_cast<char**>(faceGroupNames);

    return data;
} // Quad_v16


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Tet_v12(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/twotet4_12.2.exo";

    const size_t numVertices = 5;
    const size_t spaceDim = 3;
    const size_t cellDim = 3;
    const size_t numCells = 2;
    const size_t numCorners = 4;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::TETRAHEDRON;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -2.0,  0.0,  0.0,
        +0.0, -1.0,  0.0,
        +0.0,  1.0,  0.0,
        +0.0,  0.0,  2.0,
        +2.0,  0.0,  0.0
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0,  1,  2,  3,
        1,  4,  2,  3,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);
    static const PylithInt materialIds[numCells] = {
        7, 8,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numVertexGroups = 2;
    static const PylithInt vertexGroupSizes[2] = { 3,  4, };
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[3+4] = {
        1, 2, 3,
        0, 1, 2, 3,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[2] = {
        "mid_face",
        "bottom_face",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);

    return data;
} // Tet_v12


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Tet_v13(void) {
    TestMeshIO_Data* data = Tet_v12();assert(data);

    data->filename = "data/twotet4_13.0.exo";

    return data;
} // Tet_v13


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Tet_v16(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/smalltet_v16.exo";

    const size_t numVertices = 13;
    const size_t spaceDim = 3;
    const size_t cellDim = 3;
    const size_t numCells = 18;
    const size_t numCorners = 4;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::TETRAHEDRON;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        +0.0,  +0.0,  +0.0,
        +0.0,  +0.5,  +0.5,
        +0.0,  -0.5,  +0.5,
        +1.0,  +0.5,  +0.5,
        -1.0,  +0.5,  -0.5,
        +0.0,  +0.5,  -0.5,
        -1.0,  -0.5,  -0.5,
        -1.0,  -0.5,  +0.5,
        +0.0,  -0.5,  -0.5,
        +1.0,  +0.5,  -0.5,
        +1.0,  -0.5,  -0.5,
        -1.0,  +0.5,  +0.5,
        +1.0,  -0.5,  +0.5,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0, 1, 2, 3,
        0, 4, 5, 6,
        0, 4, 6, 7,
        0, 8, 6, 5,
        0, 3, 9, 5,
        0, 7, 2, 1,
        0, 1, 3, 5,
        0, 8, 5, 9,
        0, 10, 9, 3,
        0, 10, 2, 8,
        0, 1, 5, 4,
        0, 8, 7, 6,
        0, 8, 2, 7,
        0, 9, 10, 8,
        2, 3, 10, 0,
        4, 7, 1, 11,
        10, 3, 2, 12,
        1, 7, 4, 0,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);
    static const PylithInt materialIds[numCells] = {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numVertexGroups = 2;
    static const PylithInt vertexGroupSizes[2] = { 4,  6, };
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[4+6] = {
        4, 6, 7, 11,
        4, 5, 6, 8, 9, 10,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[2] = {
        "vertices_xneg",
        "vertices_zneg",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);

    data->numFaceGroups = 2;
    data->numFaceVertices = 3;
    static const PylithInt faceGroupSizes[2] = { 2,  4, };
    data->faceGroupSizes = const_cast<PylithInt*>(faceGroupSizes);
    static const PylithInt faceGroups[(2+4)*(1+3)] = {
        2,     4,  6,  7,
        15,    7, 11,  4,

        1,     4,  5,  6,
        3,     8,  6,  5,
        7,     8,  5,  9,
        13,    9, 10,  8,
    };
    data->faceGroups = const_cast<PylithInt*>(faceGroups);
    static const char* faceGroupNames[2] = {
        "boundary_xneg",
        "boundary_zneg",
    };
    data->faceGroupNames = const_cast<char**>(faceGroupNames);

    return data;
} // Tet_v16


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Hex_v12(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/twohex8_12.2.exo";

    const size_t numVertices = 12;
    const size_t spaceDim = 3;
    const size_t cellDim = 3;
    const size_t numCells = 2;
    const size_t numCorners = 8;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::HEXAHEDRON;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -2.0, -1.0,  1.0,
        -2.0, -1.0, -1.0,
        -2.0,  1.0, -1.0,
        -2.0,  1.0,  1.0,
        +0.0, -1.0,  1.0,
        +0.0, -1.0, -1.0,
        +0.0,  1.0, -1.0,
        +0.0,  1.0,  1.0,
        +2.0, -1.0,  1.0,
        +2.0, -1.0, -1.0,
        +2.0,  1.0, -1.0,
        +2.0,  1.0,  1.0,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0,  1,  2,  3,  4,  5,  6,  7,
        4,  5,  6,  7,  8,  9, 10, 11
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);
    static const PylithInt materialIds[numCells] = {
        7, 8,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numVertexGroups = 2;
    static const PylithInt vertexGroupSizes[2] = { 4,  6, };
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[4+6] = {
        8,  9, 10, 11,
        0,  3,  4,  7,  8, 11
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[2] = {
        "right_face",
        "top_face",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);

    return data;
} // Hex_v12


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Hex_v13(void) {
    TestMeshIO_Data* data = Hex_v12();assert(data);

    data->filename = "data/twohex8_13.0.exo";

    return data;
} // Hex_v13


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Hex_v16(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/smallhex_v16.exo";

    const size_t numVertices = 12;
    const size_t spaceDim = 3;
    const size_t cellDim = 3;
    const size_t numCells = 2;
    const size_t numCorners = 8;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::HEXAHEDRON;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -1.0,  -0.5,  +0.5,
        -1.0,  -0.5,  -0.5,
        -1.0,  +0.5,  -0.5,
        -1.0,  +0.5,  +0.5,
        +0.0,  -0.5,  +0.5,
        +0.0,  -0.5,  -0.5,
        +0.0,  +0.5,  -0.5,
        +0.0,  +0.5,  +0.5,
        +1.0,  -0.5,  +0.5,
        +1.0,  -0.5,  -0.5,
        +1.0,  +0.5,  -0.5,
        +1.0,  +0.5,  +0.5,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0, 1, 2, 3, 4, 5, 6, 7,
        4, 5, 6, 7, 8, 9, 10, 11,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);
    static const PylithInt materialIds[numCells] = {
        1, 1,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numVertexGroups = 2;
    static const PylithInt vertexGroupSizes[2] = { 4,  6, };
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[4+6] = {
        0,  1,  2,  3,
        1,  2,  5,  6,  9, 10,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[2] = {
        "vertices_xneg",
        "vertices_zneg",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);

    data->numFaceGroups = 2;
    data->numFaceVertices = 4;
    static const PylithInt faceGroupSizes[2] = { 1,  2, };
    data->faceGroupSizes = const_cast<PylithInt*>(faceGroupSizes);
    static const PylithInt faceGroups[(1+2)*(1+4)] = {
        0,   0,  3,  2,  1,

        0,   1,  2,  6,  5,
        1,   5,  6, 10,  9,
    };
    data->faceGroups = const_cast<PylithInt*>(faceGroups);
    static const char* faceGroupNames[2] = {
        "boundary_xneg",
        "boundary_zneg",
    };
    data->faceGroupNames = const_cast<char**>(faceGroupNames);

    return data;
} // Hex_v16


// End of file
