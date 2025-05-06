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

#include "TestMeshIOPetsc.hh" // Implementation of class methods

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestMeshIOPetsc_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestMeshIOPetsc_Cases {
public:

    // Data factory methods

    // Gmsh Tri
    static TestMeshIO_Data* GmshBoxTri(void);

    static TestMeshIO_Data* GmshBoxTriVerticesAscii(void);

    static TestMeshIO_Data* GmshBoxTriVerticesBinary(void);

    static TestMeshIO_Data* GmshBoxTriBoundaryAscii(void);

    static TestMeshIO_Data* GmshBoxTriBoundaryBinary(void);

    // Gmsh Quad
    static TestMeshIO_Data* GmshBoxQuad(void);

    static TestMeshIO_Data* GmshBoxQuadVerticesAscii(void);

    static TestMeshIO_Data* GmshBoxQuadVerticesBinary(void);

    static TestMeshIO_Data* GmshBoxQuadBoundaryAscii(void);

    static TestMeshIO_Data* GmshBoxQuadBoundaryBinary(void);

    // Gmsh Tet
    static TestMeshIO_Data* GmshBoxTet(void);

    static TestMeshIO_Data* GmshBoxTetVerticesAscii(void);

    static TestMeshIO_Data* GmshBoxTetVerticesBinary(void);

    static TestMeshIO_Data* GmshBoxTetBoundaryAscii(void);

    static TestMeshIO_Data* GmshBoxTetBoundaryBinary(void);

    // Gmsh Hex
    static TestMeshIO_Data* GmshBoxHex(void);

    static TestMeshIO_Data* GmshBoxHexVerticesAscii(void);

    static TestMeshIO_Data* GmshBoxHexVerticesBinary(void);

    static TestMeshIO_Data* GmshBoxHexBoundaryAscii(void);

    static TestMeshIO_Data* GmshBoxHexBoundaryBinary(void);

    // HDF5
    static TestMeshIO_Data* HDF5BoxTri(void);

    static TestMeshIO_Data* HDF5BoxQuad(void);

    static TestMeshIO_Data* HDF5BoxTet(void);

    static TestMeshIO_Data* HDF5BoxHex(void);

}; // TestMeshIOPetsc_Cases

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshIOPetsc::testFilename", "[TestMeshIOPetsc][testFilename]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTriVerticesAscii()).testFilename();
}

TEST_CASE("TestMeshIOPetsc::GmshBoxTriVerticesAscii::testRead", "[TestMeshIOPetsc][Gmsh][Tri][Vertices][ASCII][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTriVerticesAscii()).testRead(true);
}
TEST_CASE("TestMeshIOPetsc::GmshBoxTriVerticesBinary::testRead", "[TestMeshIOPetsc][Gmsh][Tri][Vertices][binary][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTriVerticesBinary()).testRead(true);
}
TEST_CASE("TestMeshIOPetsc::GmshBoxTriBoundaryAscii::testRead", "[TestMeshIOPetsc][Gmsh][Tri][Boundary][ASCII][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTriBoundaryAscii()).testRead();
}
TEST_CASE("TestMeshIOPetsc::GmshBoxTriBoundaryBinary::testRead", "[TestMeshIOPetsc][Gmsh][Tri][Boundary][binary][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTriBoundaryBinary()).testRead();
}

TEST_CASE("TestMeshIOPetsc::GmshBoxQuadVerticesAscii::testRead", "[TestMeshIOPetsc][Gmsh][Quad][Vertices][ASCII][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxQuadVerticesAscii()).testRead(true);
}
TEST_CASE("TestMeshIOPetsc::GmshBoxQuadVerticesBinary::testRead", "[TestMeshIOPetsc][Gmsh][Quad][Vertices][binary][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxQuadVerticesBinary()).testRead(true);
}
TEST_CASE("TestMeshIOPetsc::GmshBoxQuadBoundaryAscii::testRead", "[TestMeshIOPetsc][Gmsh][Quad][Boundary][ASCII][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxQuadBoundaryAscii()).testRead();
}
TEST_CASE("TestMeshIOPetsc::GmshBoxQuadBoundaryBinary::testRead", "[TestMeshIOPetsc][Gmsh][Quad][Boundary][binary][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxQuadBoundaryBinary()).testRead();
}

TEST_CASE("TestMeshIOPetsc::GmshBoxTetVerticesAscii::testRead", "[TestMeshIOPetsc][Gmsh][Tet][Vertices][ASCII][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTetVerticesAscii()).testRead(true);
}
TEST_CASE("TestMeshIOPetsc::GmshBoxTetVerticesBinary::testRead", "[TestMeshIOPetsc][Gmsh][Tet][Vertices][binary][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTetVerticesBinary()).testRead(true);
}
TEST_CASE("TestMeshIOPetsc::GmshBoxTetBoundaryAscii::testRead", "[TestMeshIOPetsc][Gmsh][Tet][Boundary][ASCII][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTetBoundaryAscii()).testRead();
}
TEST_CASE("TestMeshIOPetsc::GmshBoxTetBoundaryBinary::testRead", "[TestMeshIOPetsc][Gmsh][Tet][Boundary][binary][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTetBoundaryBinary()).testRead();
}

TEST_CASE("TestMeshIOPetsc::GmshBoxHexVerticesAscii::testRead", "[TestMeshIOPetsc][Gmsh][Hex][Vertices][ASCII][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxHexVerticesAscii()).testRead(true);
}
TEST_CASE("TestMeshIOPetsc::GmshBoxHexVerticesBinary::testRead", "[TestMeshIOPetsc][Gmsh][Hex][Vertices][binary][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxHexVerticesBinary()).testRead(true);
}
TEST_CASE("TestMeshIOPetsc::GmshBoxHexBoundaryAscii::testRead", "[TestMeshIOPetsc][Gmsh][Hex][Boundary][ASCII][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxHexBoundaryAscii()).testRead();
}
TEST_CASE("TestMeshIOPetsc::GmshBoxHexBoundaryBinary::testRead", "[TestMeshIOPetsc][Gmsh][Hex][Boundary][binary][testRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxHexBoundaryBinary()).testRead();
}

TEST_CASE("TestMeshIOPetsc::HDF5BoxTri::testWriteRead", "[TestMeshIOPetsc][HDF5][Tri][testWriteRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::HDF5BoxTri()).testWriteRead();
}
TEST_CASE("TestMeshIOPetsc::HDF5BoxQuad::testWriteRead", "[TestMeshIOPetsc][HDF5][Quad][testWriteRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::HDF5BoxQuad()).testWriteRead();
}
TEST_CASE("TestMeshIOPetsc::HDF5BoxTet::testWriteRead", "[TestMeshIOPetsc][HDF5][Tet][testWriteRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::HDF5BoxTet()).testWriteRead();
}
TEST_CASE("TestMeshIOPetsc::HDF5BoxHex::testWriteRead", "[TestMeshIOPetsc][HDF5][Hex][testWriteRead]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_Cases::HDF5BoxHex()).testWriteRead();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTri(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "NONE";

    const size_t numVertices = 9;
    const size_t spaceDim = 2;
    const size_t numCells = 8;
    const size_t cellDim = 2;
    const size_t numCorners = 3;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::TRIANGLE;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -4.00000000e+03,  -4.00000000e+03,
        +0.00000000e+00,  -4.00000000e+03,
        +4.00000000e+03,  -4.00000000e+03,
        +4.00000000e+03,  +4.00000000e+03,
        +0.00000000e+00,  +4.00000000e+03,
        -4.00000000e+03,  +4.00000000e+03,
        +4.00000000e+03,  -2.66481948e-10,
        -4.00000000e+03,  +2.66481948e-10,
        +0.00000000e+00,  -2.66481948e-10,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        8, 7, 0,
        7, 8, 4,
        1, 8, 0,
        7, 4, 5,
        3, 8, 6,
        4, 8, 3,
        2, 8, 1,
        8, 2, 6,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);

    static const PylithInt materialIds[numCells] = {
        2, 2, 2, 2,
        1, 1, 1, 1,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    return data;
} // GmshBoxTri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTriVerticesAscii(void) {
    TestMeshIO_Data* data = GmshBoxTri();assert(data);

    data->filename = "data/box_tri_vertices_ascii.msh";

    data->numVertexGroups = 6;
    static const PylithInt vertexGroupSizes[6] = {3, 3, 3, 3, 3, 1};
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[16] = {
        0, 5, 7,
        2, 3, 6,
        0, 1, 2,
        3, 4, 5,
        1, 4, 8,
        1,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[6] = {
        "boundary_xneg",
        "boundary_xpos",
        "boundary_yneg",
        "boundary_ypos",
        "fault",
        "fault_end",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);
    static const PylithInt vertexGroupTags[6] = {
        10, 11, 12, 13, 20, 21,
    };
    data->vertexGroupTags = const_cast<PylithInt*>(vertexGroupTags);

    return data;
} // GmshBoxTriVerticesAscii


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTriVerticesBinary(void) {
    TestMeshIO_Data* data = GmshBoxTriVerticesAscii();assert(data);

    data->filename = "data/box_tri_vertices_binary.msh";

    return data;
} // GmshBoxTriVerticesBinary


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTriBoundaryAscii(void) {
    TestMeshIO_Data* data = GmshBoxTri();assert(data);

    data->filename = "data/box_tri_boundary_binary.msh";

    data->numFaceGroups = 5;
    data->numFaceVertices = 2;
    static const PylithInt faceGroupSizes[5] = {2, 2, 2, 2, 2};
    data->faceGroupSizes = const_cast<PylithInt*>(faceGroupSizes);
    static const PylithInt faceGroups[(5*2)*(1+2)] = {
        0,   7, 0, // xneg
        3,   5, 7,

        4,   6, 3, // xpos
        7,   2, 6,

        2,   0, 1, // yneg
        6,   1, 2,

        3,   4, 5, // ypos
        5,   3, 4,

        1,   8, 4, // fault
        2,   1, 8,
    };
    data->faceGroups = const_cast<PylithInt*>(faceGroups);
    static const char* faceGroupNames[5] = {
        "boundary_xneg",
        "boundary_xpos",
        "boundary_yneg",
        "boundary_ypos",
        "fault",
    };
    data->faceGroupNames = const_cast<char**>(faceGroupNames);
    static const PylithInt faceGroupTags[5] = {
        10, 11, 12, 13, 20,
    };
    data->faceGroupTags = const_cast<PylithInt*>(faceGroupTags);

    return data;
} // GmshBoxTriBoundaryAscii


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTriBoundaryBinary(void) {
    TestMeshIO_Data* data = GmshBoxTriBoundaryAscii();assert(data);

    data->filename = "data/box_tri_boundary_binary.msh";

    return data;
} // GmshBoxTriBoundaryBinary


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxQuad(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/box_quad_ascii.msh";

    const size_t numVertices = 12;
    const size_t spaceDim = 2;
    const size_t numCells = 6;
    const size_t cellDim = 2;
    const size_t numCorners = 4;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::QUADRILATERAL;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -4.00000000e+03,  -4.00000000e+03,
        +0.00000000e+00,  -4.00000000e+03,
        +4.00000000e+03,  -4.00000000e+03,
        +4.00000000e+03,  +4.00000000e+03,
        +0.00000000e+00,  +4.00000000e+03,
        -4.00000000e+03,  +4.00000000e+03,
        +4.00000000e+03,  -1.33333333e+03,
        +4.00000000e+03,  +1.33333333e+03,
        -4.00000000e+03,  +1.33333333e+03,
        -4.00000000e+03,  -1.33333333e+03,
        +0.00000000e+00,  -1.33333333e+03,
        +0.00000000e+00,  +1.33333333e+03,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0, 1, 10, 9,
        9, 10, 11, 8,
        8, 11, 4, 5,
        1, 2, 6, 10,
        10, 6, 7, 11,
        11, 7, 3, 4,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);

    static const PylithInt materialIds[numCells] = {
        2, 2, 2,
        1, 1, 1,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    return data;
} // GmshBoxQuad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxQuadVerticesAscii(void) {
    TestMeshIO_Data* data = GmshBoxQuad();assert(data);

    data->filename = "data/box_quad_vertices_ascii.msh";

    data->numVertexGroups = 6;
    static const PylithInt vertexGroupSizes[6] = {4, 4, 3, 3, 4, 1};
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[19] = {
        0, 5, 8, 9,
        2, 3, 6, 7,
        0, 1, 2,
        3, 4, 5,
        1, 4, 10, 11,
        1,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[6] = {
        "boundary_xneg",
        "boundary_xpos",
        "boundary_yneg",
        "boundary_ypos",
        "fault",
        "fault_end",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);
    static const PylithInt vertexGroupTags[6] = {
        10, 11, 12, 13, 20, 21
    };
    data->vertexGroupTags = const_cast<PylithInt*>(vertexGroupTags);

    return data;
} // GmshBoxQuadVerticesAscii


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxQuadVerticesBinary(void) {
    TestMeshIO_Data* data = GmshBoxQuadVerticesAscii();assert(data);

    data->filename = "data/box_quad_vertices_binary.msh";

    return data;
} // GmshBoxQuadVerticesBinary


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxQuadBoundaryAscii(void) {
    TestMeshIO_Data* data = GmshBoxQuad();assert(data);

    data->filename = "data/box_quad_boundary_ascii.msh";

    data->numFaceGroups = 5;
    data->numFaceVertices = 2;
    static const PylithInt faceGroupSizes[5] = {3, 3, 2, 2, 3};
    data->faceGroupSizes = const_cast<PylithInt*>(faceGroupSizes);
    static const PylithInt faceGroups[(3+3+2+2+3)*(1+2)] = {
        0,   9,  0, // xneg
        1,   8,  9,
        2,   5,  8,

        3,   2, 6, // xpos
        4,   6, 7,
        5,   7, 3,

        0,   0,  1, // yneg
        3,   1,  2,

        2,   4, 5, // ypos
        5,   3, 4,

        0,    1, 10,
        1,   10, 11,
        2,   11,  4,
    };
    data->faceGroups = const_cast<PylithInt*>(faceGroups);
    static const char* faceGroupNames[5] = {
        "boundary_xneg",
        "boundary_xpos",
        "boundary_yneg",
        "boundary_ypos",
        "fault",
    };
    data->faceGroupNames = const_cast<char**>(faceGroupNames);
    static const PylithInt faceGroupTags[5] = {
        10, 11, 12, 13, 20,
    };
    data->faceGroupTags = const_cast<PylithInt*>(faceGroupTags);

    return data;
} // GmshBoxQuadBoundaryAscii


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxQuadBoundaryBinary(void) {
    TestMeshIO_Data* data = GmshBoxQuadBoundaryAscii();assert(data);

    data->filename = "data/box_quad_boundary_binary.msh";

    return data;
} // GmshBoxQuadBoundaryBinary


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTet(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "NONE";

    const size_t numVertices = 31;
    const size_t spaceDim = 3;
    const size_t numCells = 70;
    const size_t cellDim = 3;
    const size_t numCorners = 4;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::TETRAHEDRON;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -4.00000000e+03,  -4.00000000e+03,  +0.00000000e+00,
        -4.00000000e+03,  -4.00000000e+03,  -8.00000000e+03,
        -4.00000000e+03,  +4.00000000e+03,  +0.00000000e+00,
        -4.00000000e+03,  +4.00000000e+03,  -8.00000000e+03,
        +0.00000000e+00,  -4.00000000e+03,  -8.00000000e+03,
        -9.09494702e-13,  -4.00000000e+03,  +0.00000000e+00,
        -9.09494702e-13,  +4.00000000e+03,  +0.00000000e+00,
        +0.00000000e+00,  +4.00000000e+03,  -8.00000000e+03,
        +4.00000000e+03,  -4.00000000e+03,  -8.00000000e+03,
        +4.00000000e+03,  -4.00000000e+03,  +0.00000000e+00,
        +4.00000000e+03,  +4.00000000e+03,  -8.00000000e+03,
        +4.00000000e+03,  +4.00000000e+03,  +0.00000000e+00,
        -4.00000000e+03,  -4.00000000e+03,  -4.00000000e+03,
        -4.00000000e+03,  +0.00000000e+00,  +0.00000000e+00,
        -4.00000000e+03,  +4.00000000e+03,  -4.00000000e+03,
        -4.00000000e+03,  +0.00000000e+00,  -8.00000000e+03,
        -4.49418280e-13,  -4.00000000e+03,  -4.00000000e+03,
        -9.04165631e-13,  +2.27373675e-13,  +0.00000000e+00,
        -4.49418280e-13,  +4.00000000e+03,  -4.00000000e+03,
        -5.32907052e-15,  +2.27373675e-13,  -8.00000000e+03,
        +4.00000000e+03,  -4.00000000e+03,  -4.00000000e+03,
        +4.00000000e+03,  +0.00000000e+00,  -8.00000000e+03,
        +4.00000000e+03,  +0.00000000e+00,  +0.00000000e+00,
        +4.00000000e+03,  +4.00000000e+03,  -4.00000000e+03,
        -4.00000000e+03,  +0.00000000e+00,  -4.00000000e+03,
        -4.54747351e-13,  -6.43873785e-14,  -4.00000000e+03,
        +4.00000000e+03,  +0.00000000e+00,  -4.00000000e+03,
        -1.33333333e+03,  -2.00000000e+03,  -2.66666667e+03,
        +2.67363374e+03,  -2.67949949e+03,  -2.14618370e+03,
        +1.84301691e+03,  -1.40920021e+03,  -1.40568443e+03,
        +2.00000000e+03,  +2.66666667e+03,  -1.33333333e+03,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        25, 5, 17, 27,
        5, 25, 16, 27,
        24, 16, 25, 27,
        24, 12, 16, 27,
        13, 0, 12, 5,
        25, 17, 24, 27,
        15, 24, 16, 25,
        7, 15, 18, 3,
        4, 16, 15, 1,
        18, 19, 15, 7,
        19, 16, 15, 4,
        18, 15, 24, 3,
        18, 15, 25, 24,
        19, 18, 15, 25,
        19, 16, 25, 15,
        2, 17, 13, 14,
        17, 13, 14, 24,
        17, 18, 25, 24,
        2, 17, 14, 6,
        1, 24, 16, 15,
        17, 14, 18, 24,
        17, 14, 6, 18,
        24, 18, 3, 14,
        1, 24, 12, 16,
        16, 12, 5, 27,
        27, 5, 13, 12,
        13, 5, 27, 17,
        13, 24, 27, 12,
        27, 24, 13, 17,
        9, 20, 16, 28,
        5, 9, 16, 28,
        6, 18, 11, 30,
        11, 18, 23, 30,
        9, 22, 20, 28,
        26, 20, 22, 28,
        22, 6, 11, 30,
        17, 6, 22, 30,
        10, 26, 19, 21,
        19, 26, 8, 21,
        22, 26, 28, 29,
        11, 23, 22, 30,
        17, 26, 22, 29,
        26, 16, 20, 28,
        26, 23, 19, 25,
        18, 19, 23, 25,
        19, 16, 8, 25,
        26, 10, 19, 23,
        10, 19, 23, 7,
        16, 19, 8, 4,
        16, 26, 8, 25,
        19, 26, 25, 8,
        26, 16, 8, 20,
        18, 19, 7, 23,
        9, 22, 28, 29,
        6, 17, 18, 30,
        5, 9, 28, 29,
        5, 17, 9, 29,
        9, 17, 22, 29,
        28, 16, 25, 26,
        28, 25, 16, 5,
        17, 29, 25, 26,
        17, 25, 29, 5,
        29, 28, 25, 26,
        29, 25, 28, 5,
        30, 18, 25, 17,
        30, 25, 18, 23,
        26, 30, 25, 17,
        26, 25, 30, 23,
        22, 30, 26, 17,
        22, 26, 30, 23,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);

    static const PylithInt materialIds[70] = {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    return data;
} // GmshBoxTet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTetVerticesAscii(void) {
    TestMeshIO_Data* data = GmshBoxTet();assert(data);

    data->filename = "data/box_tet_vertices_ascii.msh";

    data->numVertexGroups = 8;
    static const PylithInt vertexGroupSizes[8] = {9, 9, 9, 9, 9, 9, 9, 3};
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[66] = {
        0, 1, 2, 3, 12, 13, 14, 15, 24,
        8, 9, 10, 11, 20, 21, 22, 23, 26,
        0, 1, 4, 5, 8, 9, 12, 16, 20,
        2, 3, 6, 7, 10, 11, 14, 18, 23,
        1, 3, 4, 7, 8, 10, 15, 19, 21,
        0, 2, 5, 6, 9, 11, 13, 17, 22,
        4, 5, 6, 7, 16, 17, 18, 19, 25,
        4, 5, 16,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[8] = {
        "boundary_xneg",
        "boundary_xpos",
        "boundary_yneg",
        "boundary_ypos",
        "boundary_zneg",
        "boundary_zpos",
        "fault",
        "fault_end",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);
    static const PylithInt vertexGroupTags[8] = {
        10, 11, 12, 13, 14, 15, 20, 21
    };
    data->vertexGroupTags = const_cast<PylithInt*>(vertexGroupTags);

    return data;
} // GmshBoxTetVerticesAscii


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTetVerticesBinary(void) {
    TestMeshIO_Data* data = GmshBoxTetVerticesAscii();assert(data);

    data->filename = "data/box_tet_vertices_binary.msh";

    return data;
} // GmshBoxTetVerticesBinary


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTetBoundaryAscii(void) {
    TestMeshIO_Data* data = GmshBoxTet();assert(data);

    data->filename = "data/box_tet_boundary_ascii.msh";

    data->numFaceGroups = 7;
    data->numFaceVertices = 3;
    static const PylithInt faceGroupSizes[7] = {8, 8, 8, 8, 8, 8, 8 };
    data->faceGroupSizes = const_cast<PylithInt*>(faceGroupSizes);
    static const PylithInt faceGroups[(7*8)*(1+3)] = {
        4,     0, 13, 12, // xneg
        11,   15, 24,  3,
        15,   13,  2, 14,
        16,   13, 14, 24,
        19,   24, 15,  1,
        22,    3, 24, 14,
        23,   24,  1, 12,
        27,   24, 12, 13,

        33,   22,  9, 20, // xpos
        34,   20, 26, 22,
        37,   26, 21, 10,
        38,   26,  8, 21,
        40,   23, 11, 22,
        46,   10, 23, 26,
        51,    8, 26, 20,
        69,   26, 23, 22,

        4,     0, 12,  5, // yneg
        8,    16,  1,  4,
        23,   12,  1, 16,
        24,   12, 16,  5,
        29,   20,  9, 16,
        30,    9,  5, 16,
        48,    8, 16,  4,
        51,   16,  8, 20,

        7,   18,  7,  3, // ypos
        18,  14,  2,  6,
        21,  14,  6, 18,
        22,  18,  3, 14,
        31,  18,  6, 11,
        32,  18, 11, 23,
        47,  23, 10,  7,
        52,   7, 18, 23,

        7,   15,  3,  7, // zneg
        8,   15,  4,  1,
        9,   19, 15,  7,
        10,  15, 19,  4,
        37,  19, 10, 21,
        38,   8, 19, 21,
        47,  19,  7, 10,
        48,  19,  8,  4,

        4,    0,  5, 13, // zpos
        15,  17,  2, 13,
        18,  17,  6,  2,
        26,   5, 17, 13,
        35,   6, 22, 11,
        36,   6, 17, 22,
        56,  17,  5,  9,
        57,  17,  9, 22,

        0,    5, 25, 17,// fault
        1,   25,  5, 16,
        9,   19,  7, 18,
        10,  16,  4, 19,
        13,  18, 25, 19,
        14,  16, 19, 25,
        17,  18, 17, 25,
        21,   6, 17, 18,
    };
    data->faceGroups = const_cast<PylithInt*>(faceGroups);
    static const char* faceGroupNames[8] = {
        "boundary_xneg",
        "boundary_xpos",
        "boundary_yneg",
        "boundary_ypos",
        "boundary_zneg",
        "boundary_zpos",
        "fault",
    };
    data->faceGroupNames = const_cast<char**>(faceGroupNames);
    static const PylithInt faceGroupTags[8] = {
        10, 11, 12, 13, 14, 15, 20, 21
    };
    data->faceGroupTags = const_cast<PylithInt*>(faceGroupTags);

    return data;
} // GmshBoxTetBoundaryAscii


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxTetBoundaryBinary(void) {
    TestMeshIO_Data* data = GmshBoxTetBoundaryAscii();assert(data);

    data->filename = "data/box_tet_boundary_binary.msh";

    return data;
} // GmshBoxTetBoundaryBinary


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxHex(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "NONE";

    const size_t numVertices = 48;
    const size_t spaceDim = 3;
    const size_t numCells = 18;
    const size_t cellDim = 3;
    const size_t numCorners = 8;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::HEXAHEDRON;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        -4.00000000e+03,  -4.00000000e+03,  +0.00000000e+00,
        -4.00000000e+03,  -4.00000000e+03,  -8.00000000e+03,
        -4.00000000e+03,  +4.00000000e+03,  +0.00000000e+00,
        -4.00000000e+03,  +4.00000000e+03,  -8.00000000e+03,
        +0.00000000e+00,  -4.00000000e+03,  -8.00000000e+03,
        -9.09494702e-13,  -4.00000000e+03,  +0.00000000e+00,
        -9.09494702e-13,  +4.00000000e+03,  +0.00000000e+00,
        +0.00000000e+00,  +4.00000000e+03,  -8.00000000e+03,
        +4.00000000e+03,  -4.00000000e+03,  -8.00000000e+03,
        +4.00000000e+03,  -4.00000000e+03,  +0.00000000e+00,
        +4.00000000e+03,  +4.00000000e+03,  -8.00000000e+03,
        +4.00000000e+03,  +4.00000000e+03,  +0.00000000e+00,
        -4.00000000e+03,  -4.00000000e+03,  -5.33333333e+03,
        -4.00000000e+03,  -4.00000000e+03,  -2.66666667e+03,
        -4.00000000e+03,  -1.33333333e+03,  +0.00000000e+00,
        -4.00000000e+03,  +1.33333333e+03,  +0.00000000e+00,
        -4.00000000e+03,  +4.00000000e+03,  -5.33333333e+03,
        -4.00000000e+03,  +4.00000000e+03,  -2.66666667e+03,
        -4.00000000e+03,  -1.33333333e+03,  -8.00000000e+03,
        -4.00000000e+03,  +1.33333333e+03,  -8.00000000e+03,
        -5.97448017e-13,  -4.00000000e+03,  -2.66666667e+03,
        -3.01388544e-13,  -4.00000000e+03,  -5.33333333e+03,
        -9.04165631e-13,  -1.33333333e+03,  +0.00000000e+00,
        -9.04165631e-13,  +1.33333333e+03,  +0.00000000e+00,
        -5.97448017e-13,  +4.00000000e+03,  -2.66666667e+03,
        -3.01388544e-13,  +4.00000000e+03,  -5.33333333e+03,
        -5.32907052e-15,  -1.33333333e+03,  -8.00000000e+03,
        -5.32907052e-15,  +1.33333333e+03,  -8.00000000e+03,
        +4.00000000e+03,  -4.00000000e+03,  -5.33333333e+03,
        +4.00000000e+03,  -4.00000000e+03,  -2.66666667e+03,
        +4.00000000e+03,  -1.33333333e+03,  -8.00000000e+03,
        +4.00000000e+03,  +1.33333333e+03,  -8.00000000e+03,
        +4.00000000e+03,  -1.33333333e+03,  +0.00000000e+00,
        +4.00000000e+03,  +1.33333333e+03,  +0.00000000e+00,
        +4.00000000e+03,  +4.00000000e+03,  -5.33333333e+03,
        +4.00000000e+03,  +4.00000000e+03,  -2.66666667e+03,
        -4.00000000e+03,  -1.33333333e+03,  -5.33333333e+03,
        -4.00000000e+03,  +1.33333333e+03,  -5.33333333e+03,
        -4.00000000e+03,  -1.33333333e+03,  -2.66666667e+03,
        -4.00000000e+03,  +1.33333333e+03,  -2.66666667e+03,
        -6.02777088e-13,  -1.33333333e+03,  -2.66666667e+03,
        -3.06717614e-13,  -1.33333333e+03,  -5.33333333e+03,
        -6.02777088e-13,  +1.33333333e+03,  -2.66666667e+03,
        -3.06717614e-13,  +1.33333333e+03,  -5.33333333e+03,
        +4.00000000e+03,  -1.33333333e+03,  -5.33333333e+03,
        +4.00000000e+03,  +1.33333333e+03,  -5.33333333e+03,
        +4.00000000e+03,  -1.33333333e+03,  -2.66666667e+03,
        +4.00000000e+03,  +1.33333333e+03,  -2.66666667e+03,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        36, 12, 1, 18, 41, 21, 4, 26,
        37, 36, 18, 19, 43, 41, 26, 27,
        16, 37, 19, 3, 25, 43, 27, 7,
        38, 13, 12, 36, 40, 20, 21, 41,
        39, 38, 36, 37, 42, 40, 41, 43,
        17, 39, 37, 16, 24, 42, 43, 25,
        14, 0, 13, 38, 22, 5, 20, 40,
        15, 14, 38, 39, 23, 22, 40, 42,
        2, 15, 39, 17, 6, 23, 42, 24,
        28, 8, 4, 21, 44, 30, 26, 41,
        44, 30, 26, 41, 45, 31, 27, 43,
        45, 31, 27, 43, 34, 10, 7, 25,
        29, 28, 21, 20, 46, 44, 41, 40,
        46, 44, 41, 40, 47, 45, 43, 42,
        47, 45, 43, 42, 35, 34, 25, 24,
        9, 29, 20, 5, 32, 46, 40, 22,
        32, 46, 40, 22, 33, 47, 42, 23,
        33, 47, 42, 23, 11, 35, 24, 6,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);

    static const PylithInt materialIds[numCells] = {
        1, 1, 1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2, 2,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    return data;
} // GmshBoxHex


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxHexVerticesAscii(void) {
    TestMeshIO_Data* data = GmshBoxHex();assert(data);

    data->filename = "data/box_hex_vertices_ascii.msh";

    data->numVertexGroups = 8;
    static const PylithInt vertexGroupSizes[8] = {16, 16, 12, 12, 12, 12, 16, 4};
    data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);
    static const PylithInt vertexGroups[100] = {
        0, 1, 2, 3, 12, 13, 14, 15, 16, 17, 18, 19, 36, 37, 38, 39,
        8, 9, 10, 11, 28, 29, 30, 31, 32, 33, 34, 35, 44, 45, 46, 47,
        0, 1, 4, 5, 8, 9, 12, 13, 20, 21, 28, 29,
        2, 3, 6, 7, 10, 11, 16, 17, 24, 25, 34, 35,
        1, 3, 4, 7, 8, 10, 18, 19, 26, 27, 30, 31,
        0, 2, 5, 6, 9, 11, 14, 15, 22, 23, 32, 33,
        4, 5, 6, 7, 20, 21, 22, 23, 24, 25, 26, 27, 40, 41, 42, 43,
        4, 5, 20, 21,
    };
    data->vertexGroups = const_cast<PylithInt*>(vertexGroups);
    static const char* vertexGroupNames[9] = {
        "boundary_xneg",
        "boundary_xpos",
        "boundary_yneg",
        "boundary_ypos",
        "boundary_zneg",
        "boundary_zpos",
        "fault",
        "fault_end",
    };
    data->vertexGroupNames = const_cast<char**>(vertexGroupNames);
    static const PylithInt vertexGroupTags[8] = {
        10, 11, 12, 13, 14, 15, 20, 21,
    };
    data->vertexGroupTags = const_cast<PylithInt*>(vertexGroupTags);

    return data;
} // GmshBoxHexVerticesAscii


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxHexVerticesBinary(void) {
    TestMeshIO_Data* data = GmshBoxHexVerticesAscii();assert(data);

    data->filename = "data/box_hex_vertices_binary.msh";

    return data;
} // GmshBoxHexVerticesBinary


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxHexBoundaryAscii(void) {
    TestMeshIO_Data* data = GmshBoxHex();assert(data);

    data->filename = "data/box_hex_boundary_ascii.msh";

    data->numFaceGroups = 7;
    data->numFaceVertices = 4;
    static const PylithInt faceGroupSizes[7] = {9, 9, 6, 6, 6, 6, 9};
    data->faceGroupSizes = const_cast<PylithInt*>(faceGroupSizes);
    static const PylithInt faceGroups[(9+9+6+6+6+6+9)*(1+4)] = {
        0,   36, 18,  1, 12, // xneg
        1,   37, 19, 18, 36,
        2,   16,  3, 19, 37,
        3,   38, 36, 12, 13,
        4,   39, 37, 36, 38,
        5,   17, 16, 37, 39,
        6,   14, 38, 13,  0,
        7,   15, 39, 38, 14,
        8,    2, 17, 39, 15,

        9,   28,  8, 30, 44, // xpos
        10,  44, 30, 31, 45,
        11,  45, 31, 10, 34,
        12,  29, 28, 44, 46,
        13,  46, 44, 45, 47,
        14,  47, 45, 34, 35,
        15,   9, 29, 46, 32,
        16,  32, 46, 47, 33,
        17,  33, 47, 35, 11,

        0,   12,  1,  4, 21, // yneg
        3,   13, 12, 21, 20,
        6,    0, 13, 20,  5,
        9,   28, 21,  4,  8,
        12,  29, 20, 21, 28,
        15,   9,  5, 20, 29,

        2,   16, 25,  7,  3, // ypos
        5,   17, 24, 25, 16,
        8,    2,  6, 24, 17,
        11,  34, 10,  7, 25,
        14,  35, 34, 25, 24,
        17,  11, 35, 24,  6,

        0,    1, 18, 26,  4, // zneg
        1,   18, 19, 27, 26,
        2,   19,  3,  7, 27,
        9,    8,  4, 26, 30,
        10,  30, 26, 27, 31,
        11,  31, 27,  7, 10,

        6,   14,  0,  5, 22, // zpos
        7,   15, 14, 22, 23,
        8,    2, 15, 23,  6,
        15,   9, 32, 22,  5,
        16,  32, 33, 23, 22,
        17,  33, 11,  6, 23,

        0,   41, 21,  4, 26, // fault
        1,   43, 41, 26, 27,
        2,   25, 43, 27,  7,
        3,   40, 20, 21, 41,
        4,   42, 40, 41, 43,
        5,   24, 42, 43, 25,
        6,   22,  5, 20, 40,
        7,   23, 22, 40, 42,
        8,    6, 23, 42, 24,
    };
    data->faceGroups = const_cast<PylithInt*>(faceGroups);
    static const char* faceGroupNames[7] = {
        "boundary_xneg",
        "boundary_xpos",
        "boundary_yneg",
        "boundary_ypos",
        "boundary_zneg",
        "boundary_zpos",
        "fault",
    };
    data->faceGroupNames = const_cast<char**>(faceGroupNames);
    static const PylithInt faceGroupTags[7] = {
        10, 11, 12, 13, 14, 15, 20,
    };
    data->faceGroupTags = const_cast<PylithInt*>(faceGroupTags);

    return data;
} // GmshBoxHexBoundaryAscii


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::GmshBoxHexBoundaryBinary(void) {
    TestMeshIO_Data* data = GmshBoxHexBoundaryAscii();assert(data);

    data->filename = "data/box_hex_boundary_binary.msh";

    return data;
} // GmshBoxHexBoundaryBinary


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::HDF5BoxTri(void) {
    TestMeshIO_Data* data = GmshBoxTriBoundaryBinary();assert(data);

    data->filename = "box_tri.h5";
    static const PylithInt faceGroupTags[5] = {
        1, 1, 1, 1, 1,
    };
    data->faceGroupTags = const_cast<PylithInt*>(faceGroupTags);

    return data;
} // HDF5BoxTri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::HDF5BoxQuad(void) {
    TestMeshIO_Data* data = GmshBoxQuadBoundaryBinary();assert(data);

    data->filename = "box_quad.h5";
    static const PylithInt faceGroupTags[5] = {
        1, 1, 1, 1, 1,
    };
    data->faceGroupTags = const_cast<PylithInt*>(faceGroupTags);

    return data;
} // HDF5BoxQuad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::HDF5BoxTet(void) {
    TestMeshIO_Data* data = GmshBoxTetBoundaryBinary();assert(data);

    data->filename = "box_tet.h5";
    static const PylithInt faceGroupTags[7] = {
        1, 1, 1, 1, 1, 1, 1,
    };
    data->faceGroupTags = const_cast<PylithInt*>(faceGroupTags);

    return data;
} // HDF5BoxTet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_Cases::HDF5BoxHex(void) {
    TestMeshIO_Data* data = GmshBoxHexBoundaryBinary();assert(data);

    data->filename = "box_hex.h5";
    static const PylithInt faceGroupTags[7] = {
        1, 1, 1, 1, 1, 1, 1,
    };
    data->faceGroupTags = const_cast<PylithInt*>(faceGroupTags);

    return data;
} // HDF5BoxHex


// End of file
