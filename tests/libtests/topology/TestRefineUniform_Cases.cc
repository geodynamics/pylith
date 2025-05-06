// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
//

#include <portinfo>

#include "TestRefineUniform.hh" // Implementation of class methods

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace topology {
        class TestRefineUniform_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestRefineUniform_Cases {
public:

    // Data factory methods
    static TestRefineUniform_Data* Tri_2xNoFault(void);

    static TestRefineUniform_Data* Tri_2xFault(void);

    static TestRefineUniform_Data* Quad_2xNoFault(void);

    static TestRefineUniform_Data* Quad_2xFault(void);

    static TestRefineUniform_Data* Tet_2xNoFault(void);

    static TestRefineUniform_Data* Tet_2xFault(void);

    static TestRefineUniform_Data* Hex_2xNoFault(void);

    static TestRefineUniform_Data* Hex_2xFault(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestRefineUniform::Tri_2xNoFault::testRefine", "[TestRefineUniform][Tri][NoFault]") {
    pylith::topology::TestRefineUniform(pylith::topology::TestRefineUniform_Cases::Tri_2xNoFault()).testRefine();
}
TEST_CASE("TestRefineUniform::Tri_2xFault::testRefine", "[TestRefineUniform][Tri][Fault]") {
    pylith::topology::TestRefineUniform(pylith::topology::TestRefineUniform_Cases::Tri_2xFault()).testRefine();
}

TEST_CASE("TestRefineUniform::Quad_2xNoFault::testRefine", "[TestRefineUniform][Quad][NoFault]") {
    pylith::topology::TestRefineUniform(pylith::topology::TestRefineUniform_Cases::Quad_2xNoFault()).testRefine();
}
TEST_CASE("TestRefineUniform::Quad_2xFault::testRefine", "[TestRefineUniform][Quad][Fault]") {
    pylith::topology::TestRefineUniform(pylith::topology::TestRefineUniform_Cases::Quad_2xFault()).testRefine();
}

TEST_CASE("TestRefineUniform::Tet_2xNoFault::testRefine", "[TestRefineUniform][Tet][NoFault]") {
    pylith::topology::TestRefineUniform(pylith::topology::TestRefineUniform_Cases::Tet_2xNoFault()).testRefine();
}
TEST_CASE("TestRefineUniform::Tet_2xFault::testRefine", "[TestRefineUniform][Tet][Fault]") {
    pylith::topology::TestRefineUniform(pylith::topology::TestRefineUniform_Cases::Tet_2xFault()).testRefine();
}

TEST_CASE("TestRefineUniform::Hex_2xNoFault::testRefine", "[TestRefineUniform][Hex][NoFault]") {
    pylith::topology::TestRefineUniform(pylith::topology::TestRefineUniform_Cases::Hex_2xNoFault()).testRefine();
}
TEST_CASE("TestRefineUniform::Hex_2xFault::testRefine", "[TestRefineUniform][Hex][Fault]") {
    pylith::topology::TestRefineUniform(pylith::topology::TestRefineUniform_Cases::Hex_2xFault()).testRefine();
}

// ------------------------------------------------------------------------------------------------
pylith::topology::TestRefineUniform_Data*
pylith::topology::TestRefineUniform_Cases::Tri_2xNoFault(void) {
    TestRefineUniform_Data* data = new TestRefineUniform_Data();assert(data);

    data->filename = "data/fourtri3.mesh";
    data->refineLevel = 1;
    data->faultA = NULL;
    data->faultB = NULL;
    data->isSimplexMesh = true;

    data->numVertices = 13;
    data->spaceDim = 2;
    data->numCells = 16;
    data->numCellsCohesive = 0;
    data->cellDim = 2;
    data->numCorners = 3;
    data->numCornersCohesive = 4;

    data->matIdSum = 8*1 + 8*2;

    data->numVertexGroups = 4;
    static const size_t _vertexGroupSizes[4] = {
        5, 3, 2, 5,
    };
    data->vertexGroupSizes = const_cast<size_t*>(_vertexGroupSizes);

    static const char* _vertexGroupNames[4] = {
        "bc1_vertices",
        "bc2_vertices",
        "endpoints",
        "fault",
    };
    data->vertexGroupNames = const_cast<const char**>(_vertexGroupNames);

    data->numFaceGroups = 3;
    static const size_t _faceGroupSizes[3] = {
        4, 2, 4,
    };
    data->faceGroupSizes = const_cast<size_t*>(_faceGroupSizes);

    static const char* _faceGroupNames[3] = {
        "bc1",
        "bc2",
        "fault_faces",
    };
    data->faceGroupNames = const_cast<const char**>(_faceGroupNames);

    return data;
} // Tri_2xNoFault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestRefineUniform_Data*
pylith::topology::TestRefineUniform_Cases::Tri_2xFault(void) {
    TestRefineUniform_Data* data = new TestRefineUniform_Data();assert(data);

    data->filename = "data/fourtri3.mesh";
    data->refineLevel = 1;
    data->faultA = "fault";
    data->faultB = NULL;
    data->isSimplexMesh = true;

    data->numVertices = 18;
    data->spaceDim = 2;
    data->numCells = 16;
    data->numCellsCohesive = 4;
    data->cellDim = 2;
    data->numCorners = 3;
    data->numCornersCohesive = 4;

    data->matIdSum = 8*1 + 8*2 + 4*100;

    static const char* _vertexGroupNames[4] = {
        "bc1_vertices",
        "bc2_vertices",
        "endpoints",
        "fault",
    };
    data->vertexGroupNames = const_cast<const char**>(_vertexGroupNames);

    data->numFaceGroups = 3;
    static const size_t _faceGroupSizes[3] = {
        4, 2, 8,
    };
    data->faceGroupSizes = const_cast<size_t*>(_faceGroupSizes);

    static const char* _faceGroupNames[3] = {
        "bc1",
        "bc2",
        "fault_faces",
    };
    data->faceGroupNames = const_cast<const char**>(_faceGroupNames);

    return data;
} // Tri_2xFault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestRefineUniform_Data*
pylith::topology::TestRefineUniform_Cases::Quad_2xNoFault(void) {
    TestRefineUniform_Data* data = new TestRefineUniform_Data();assert(data);

    data->filename = "data/fourquad4.mesh";
    data->refineLevel = 1;
    data->faultA = NULL;
    data->faultB = NULL;
    data->isSimplexMesh = false;

    data->numVertices = 25;
    data->spaceDim = 2;
    data->numCells = 16;
    data->numCellsCohesive = 0;
    data->cellDim = 2;
    data->numCorners = 4;
    data->numCornersCohesive = 4;

    data->matIdSum = 8*1 + 8*2;

    data->numVertexGroups = 3;
    static const size_t _vertexGroupSizes[3] = {
        5, 5, 5,
    };
    data->vertexGroupSizes = const_cast<size_t*>(_vertexGroupSizes);
    static const char* _vertexGroupNames[3] = {
        "fault",
        "endpoints",
        "bc_vertices",
    };
    data->vertexGroupNames = const_cast<const char**>(_vertexGroupNames);

    data->numFaceGroups = 2;
    static const size_t _faceGroupSizes[2] = {
        4, 4,
    };
    data->faceGroupSizes = const_cast<size_t*>(_faceGroupSizes);
    static const char* _faceGroupNames[2] = {
        "fault_faces",
        "bc",
    };
    data->faceGroupNames = const_cast<const char**>(_faceGroupNames);

    return data;
} // Quad_2dNoFault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestRefineUniform_Data*
pylith::topology::TestRefineUniform_Cases::Quad_2xFault(void) {
    TestRefineUniform_Data* data = new TestRefineUniform_Data();assert(data);

    data->filename = "data/fourquad4.mesh";
    data->refineLevel = 1;
    data->faultA = "fault";
    data->faultB = NULL;
    data->isSimplexMesh = false;

    data->numVertices = 30;
    data->spaceDim = 2;
    data->numCells = 16;
    data->numCellsCohesive = 4;
    data->cellDim = 2;
    data->numCorners = 4;
    data->numCornersCohesive = 4;

    data->matIdSum = 8*1 + 8*2 + 4*100;

    data->numVertexGroups = 3;
    static const size_t _vertexGroupSizes[3] = {
        10, 5, 5+1,
    };
    data->vertexGroupSizes = const_cast<size_t*>(_vertexGroupSizes);
    static const char* _vertexGroupNames[3] = {
        "fault",
        "endpoints",
        "bc_vertices",
    };
    data->vertexGroupNames = const_cast<const char**>(_vertexGroupNames);

    data->numFaceGroups = 2;
    static const size_t _faceGroupSizes[2] = {
        8, 4,
    };
    data->faceGroupSizes = const_cast<size_t*>(_faceGroupSizes);
    static const char* _faceGroupNames[2] = {
        "fault_faces",
        "bc",
    };
    data->faceGroupNames = const_cast<const char**>(_faceGroupNames);

    return data;
} // Quad_2xFault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestRefineUniform_Data*
pylith::topology::TestRefineUniform_Cases::Tet_2xNoFault(void) {
    TestRefineUniform_Data* data = new TestRefineUniform_Data();assert(data);

    data->filename = "data/twotet4.mesh";
    data->refineLevel = 1;
    data->faultA = NULL;
    data->faultB = NULL;
    data->isSimplexMesh = true;

    data->numVertices = 14;
    data->spaceDim = 3;
    data->numCells = 16;
    data->numCellsCohesive = 0;
    data->cellDim = 3;
    data->numCorners = 4;
    data->numCornersCohesive = 6;

    data->matIdSum = 8*1 + 8*2;

    data->numVertexGroups = 4;
    static const size_t _vertexGroupSizes[4] = {
        3, 3, 2, 6,
    };
    data->vertexGroupSizes = const_cast<size_t*>(_vertexGroupSizes);

    static const char* _vertexGroupNames[4] = {
        "bc1_vertices",
        "bc2_vertices",
        "endpoints",
        "fault",
    };
    data->vertexGroupNames = const_cast<const char**>(_vertexGroupNames);

    data->numFaceGroups = 3;
    static const size_t _faceGroupSizes[3] = {
        4, 4, 4,
    };
    data->faceGroupSizes = const_cast<size_t*>(_faceGroupSizes);
    static const char* _faceGroupNames[3] = {
        "bc1",
        "bc2",
        "fault_faces",
    };
    data->faceGroupNames = const_cast<const char**>(_faceGroupNames);

    return data;
} // Tet_2xNoFault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestRefineUniform_Data*
pylith::topology::TestRefineUniform_Cases::Tet_2xFault(void) {
    TestRefineUniform_Data* data = new TestRefineUniform_Data();assert(data);

    data->filename = "data/twotet4.mesh";
    data->refineLevel = 1;
    data->faultA = "fault";
    data->faultB = NULL;
    data->isSimplexMesh = true;

    data->numVertices = 20;
    data->spaceDim = 3;
    data->numCells = 16;
    data->numCellsCohesive = 4;
    data->cellDim = 3;
    data->numCorners = 4;
    data->numCornersCohesive = 6;

    data->matIdSum = 8*1 + 8*2 + 4*100;

    data->numVertexGroups = 4;

    data->numVertexGroups = 4;
    static const size_t _vertexGroupSizes[4] = {
        3+1, 3+1, 2, 12,
    };
    data->vertexGroupSizes = const_cast<size_t*>(_vertexGroupSizes);

    static const char* _vertexGroupNames[4] = {
        "bc1_vertices",
        "bc2_vertices",
        "endpoints",
        "fault",
    };
    data->vertexGroupNames = const_cast<const char**>(_vertexGroupNames);

    data->numFaceGroups = 3;
    static const size_t _faceGroupSizes[3] = {
        4, 4, 8,
    };
    data->faceGroupSizes = const_cast<size_t*>(_faceGroupSizes);
    static const char* _faceGroupNames[3] = {
        "bc1",
        "bc2",
        "fault_faces",
    };
    data->faceGroupNames = const_cast<const char**>(_faceGroupNames);

    return data;
} // Tet_2xFault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestRefineUniform_Data*
pylith::topology::TestRefineUniform_Cases::Hex_2xNoFault(void) {
    TestRefineUniform_Data* data = new TestRefineUniform_Data();assert(data);

    data->filename = "data/twohex8.mesh";
    data->refineLevel = 1;
    data->faultA = NULL;
    data->faultB = NULL;
    data->isSimplexMesh = true;

    data->numVertices = 45;
    data->spaceDim = 3;
    data->numCells = 16;
    data->numCellsCohesive = 0;
    data->cellDim = 3;
    data->numCorners = 8;
    data->numCornersCohesive = 8;

    data->matIdSum = 8*1 + 8*2;

    data->numVertexGroups = 4;
    static const size_t _vertexGroupSizes[4] = {
        9, 2, 9, 9,
    };
    data->vertexGroupSizes = const_cast<size_t*>(_vertexGroupSizes);
    static const char* _vertexGroupNames[4] = {
        "fault",
        "endpoints",
        "face1_vertices",
        "face2_vertices",
    };
    data->vertexGroupNames = const_cast<const char**>(_vertexGroupNames);

    data->numFaceGroups = 3;
    static const size_t _faceGroupSizes[3] = {
        4, 4, 4,
    };
    data->faceGroupSizes = const_cast<size_t*>(_faceGroupSizes);
    static const char* _faceGroupNames[3] = {
        "fault_faces",
        "face1",
        "face2",
    };
    data->faceGroupNames = const_cast<const char**>(_faceGroupNames);

    return data;
} // Hex2x_NoFault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestRefineUniform_Data*
pylith::topology::TestRefineUniform_Cases::Hex_2xFault(void) {
    TestRefineUniform_Data* data = new TestRefineUniform_Data();assert(data);

    data->filename = "data/twohex8.mesh";
    data->refineLevel = 1;
    data->faultA = "fault";
    data->faultB = NULL;
    data->isSimplexMesh = true;

    data->numVertices = 54;
    data->spaceDim = 3;
    data->numCells = 16;
    data->numCellsCohesive = 4;
    data->cellDim = 3;
    data->numCorners = 8;
    data->numCornersCohesive = 8;

    data->matIdSum = 8*1 + 8*2 + 4*100;

    data->numVertexGroups = 4;
    static const size_t _vertexGroupSizes[4] = {
        18, 2, 9, 12,
    };
    data->vertexGroupSizes = const_cast<size_t*>(_vertexGroupSizes);
    static const char* _vertexGroupNames[4] = {
        "fault",
        "endpoints",
        "face1_vertices",
        "face2_vertices",
    };
    data->vertexGroupNames = const_cast<const char**>(_vertexGroupNames);

    data->numFaceGroups = 3;
    static const size_t _faceGroupSizes[3] = {
        8, 4, 4,
    };
    data->faceGroupSizes = const_cast<size_t*>(_faceGroupSizes);
    static const char* _faceGroupNames[3] = {
        "fault_faces",
        "face1",
        "face2",
    };
    data->faceGroupNames = const_cast<const char**>(_faceGroupNames);

    return data;
} // Hex_2xFault


// End of file
