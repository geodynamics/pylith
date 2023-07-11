// -*- C++ -*-
//
// -----------------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// -----------------------------------------------------------------------------
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

    data->numGroups = 4;

    static const int _groupSizes[4] = {
        9, 5, 2, 9,
    };
    data->groupSizes = const_cast<int*>(_groupSizes);

    static const char* _groupNames[4] = {
        "edge 1",
        "edge 2",
        "end points",
        "fault",
    };
    data->groupNames = const_cast<const char**>(_groupNames);

    static const char* _groupTypes[4] = {
        "vertex",
        "vertex",
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<const char**>(_groupTypes);

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

    data->numGroups = 4;

    static const int _groupSizes[4] = {
        11, 6, 2, 18,
    };
    data->groupSizes = const_cast<int*>(_groupSizes);

    static const char* _groupNames[4] = {
        "edge 1",
        "edge 2",
        "end points",
        "fault",
    };
    data->groupNames = const_cast<const char**>(_groupNames);

    static const char* _groupTypes[4] = {
        "vertex",
        "vertex",
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<const char**>(_groupTypes);

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

    data->numGroups = 3;

    static const int _groupSizes[3] = {
        9, 9, 9,
    };
    data->groupSizes = const_cast<int*>(_groupSizes);

    static const char* _groupNames[3] = {
        "edge 1",
        "end points",
        "fault",
    };
    data->groupNames = const_cast<const char**>(_groupNames);

    static const char* _groupTypes[3] = {
        "vertex",
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<const char**>(_groupTypes);

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

    data->numGroups = 3;

    static const int _groupSizes[3] = {
        6+4, // vertices+edges
        5+4,
        10+8,
    };
    data->groupSizes = const_cast<int*>(_groupSizes);

    static const char* _groupNames[3] = {
        "edge 1",
        "end points",
        "fault",
    };
    data->groupNames = const_cast<const char**>(_groupNames);

    static const char* _groupTypes[3] = {
        "vertex",
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<const char**>(_groupTypes);

    return data;
} // Quad_2dFault


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

    data->numGroups = 4;

    static const int _groupSizes[4] = {
        3+2+0, // vertices+edges+faces
        3+2+0,
        2+0+0,
        6+9+4
    };
    data->groupSizes = const_cast<int*>(_groupSizes);

    static const char* _groupNames[4] = {
        "edge 1",
        "edge 2",
        "end points",
        "fault",
    };
    data->groupNames = const_cast<const char**>(_groupNames);

    static const char* _groupTypes[4] = {
        "vertex",
        "vertex",
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<const char**>(_groupTypes);

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

    data->numGroups = 4;

    static const int _groupSizes[4] = {
        6,
        6,
        2,
        38,
    };
    data->groupSizes = const_cast<int*>(_groupSizes);

    static const char* _groupNames[4] = {
        "edge 1",
        "edge 2",
        "end points",
        "fault",
    };
    data->groupNames = const_cast<const char**>(_groupNames);

    static const char* _groupTypes[4] = {
        "vertex",
        "vertex",
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<const char**>(_groupTypes);

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

    data->numGroups = 4;

    static const int _groupSizes[4] = {
        2+0+0, // vertices+edges+faces
        9+12+4,
        9+12+4,
        9+12+4,
    };
    data->groupSizes = const_cast<int*>(_groupSizes);

    static const char* _groupNames[4] = {
        "end points",
        "face 1",
        "face 2",
        "fault",
    };
    data->groupNames = const_cast<const char**>(_groupNames);

    static const char* _groupTypes[4] = {
        "vertex",
        "vertex",
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<const char**>(_groupTypes);

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

    data->numGroups = 4;

    static const int _groupSizes[4] = {
        2+0+0, // vertices+edges+faces
        9+12+4,
        12+14+4,
        18+24+8,
    };
    data->groupSizes = const_cast<int*>(_groupSizes);

    static const char* _groupNames[4] = {
        "end points",
        "face 1",
        "face 2",
        "fault",
    };
    data->groupNames = const_cast<const char**>(_groupNames);

    static const char* _groupTypes[4] = {
        "vertex",
        "vertex",
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<const char**>(_groupTypes);

    return data;
} // Hex_2xFault


// End of file
