// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

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

    static TestMeshIO_Data* Quad_v12(void);

    static TestMeshIO_Data* Quad_v13(void);

    static TestMeshIO_Data* Tet_v12(void);

    static TestMeshIO_Data* Tet_v13(void);

    static TestMeshIO_Data* Hex_v12(void);

    static TestMeshIO_Data* Hex_v13(void);

}; // TestMeshIOCubit_Cases

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshIOCubit::Tri::testFilename", "[TestMeshIOCubit][testFilename]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tri_v12()).testFilename();
}

TEST_CASE("TestMeshIOCubit::Tri_v12::testRead", "[TestMeshIOCubit][Tri][v12][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tri_v12()).testRead();
}
TEST_CASE("TestMeshIOCubit::Tri::testRead", "[TestMeshIOCubit][Tri][v13][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tri_v13()).testRead();
}

TEST_CASE("TestMeshIOCubit::Quad_v12::testRead", "[TestMeshIOCubit][Quad][v12][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Quad_v12()).testRead();
}
TEST_CASE("TestMeshIOCubit::Quad::testRead", "[TestMeshIOCubit][Quad][v13][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Quad_v13()).testRead();
}

TEST_CASE("TestMeshIOCubit::Tet_v12::testRead", "[TestMeshIOCubit][Tet][v12][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tet_v12()).testRead();
}
TEST_CASE("TestMeshIOCubit::Tet::testRead", "[TestMeshIOCubit][Tet][v13][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Tet_v13()).testRead();
}

TEST_CASE("TestMeshIOCubit::Hex_v12::testRead", "[TestMeshIOCubit][Hex][v12][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Hex_v12()).testRead();
}
TEST_CASE("TestMeshIOCubit::Hex::testRead", "[TestMeshIOCubit][Hex][v13][testRead]") {
    pylith::meshio::TestMeshIOCubit(pylith::meshio::TestMeshIOCubit_Cases::Hex_v13()).testRead();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Tri_v12(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/twotri3_12.2.exo";

    data->numVertices = 4;
    data->spaceDim = 2;
    data->numCells = 2;
    data->cellDim = 2;
    data->numCorners = 3;

    static const PylithScalar vertices[4*2] = {
        -1.0,  +0.0,
        +0.0,  -1.0,
        +0.0,  +1.0,
        +1.0,  +0.0
    };
    data->vertices = const_cast<PylithScalar*>(vertices);

    static const PylithInt cells[2*3] = {
        0,  1,  2,
        2,  1,  3,
    };
    data->cells = const_cast<PylithInt*>(cells);
    static const PylithInt materialIds[2] = {
        2, 3,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numGroups = 2;
    static const PylithInt groupSizes[2] = { 1,  2, };
    data->groupSizes = const_cast<PylithInt*>(groupSizes);
    static const PylithInt groups[1+2] = {
        0,
        2, 3,
    };
    data->groups = const_cast<PylithInt*>(groups);
    static const char* groupNames[2] = {
        "left_vertex",
        "right_vertex",
    };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[2] = {
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<char**>(groupTypes);

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
pylith::meshio::TestMeshIOCubit_Cases::Quad_v12(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/twoquad4_12.2.exo";

    data->numVertices = 6;
    data->spaceDim = 2;
    data->numCells = 2;
    data->cellDim = 2;
    data->numCorners = 4;

    static const PylithScalar vertices[6*2] = {
        0.0,  0.0,
        1.0,  0.0,
        1.0,  1.0,
        0.0,  1.0,
        2.0,  0.0,
        2.0,  1.0,
    };
    data->vertices = const_cast<PylithScalar*>(vertices);

    static const PylithInt cells[2*4] = {
        0,  1,  2,  3,
        1,  4,  5,  2,
    };
    data->cells = const_cast<PylithInt*>(cells);
    static const PylithInt materialIds[2] = {
        10, 11,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numGroups = 2;
    static const PylithInt groupSizes[2] = { 2,  3, };
    data->groupSizes = const_cast<PylithInt*>(groupSizes);
    static const PylithInt groups[2+3] = {
        0, 3,
        2, 3, 5,
    };
    data->groups = const_cast<PylithInt*>(groups);
    static const char* groupNames[2] = {
        "left_edge",
        "top_edge",
    };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[2] = {
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<char**>(groupTypes);

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
pylith::meshio::TestMeshIOCubit_Cases::Tet_v12(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/twotet4_12.2.exo";

    data->numVertices = 5;
    data->spaceDim = 3;
    data->numCells = 2;
    data->cellDim = 3;
    data->numCorners = 4;

    static const PylithScalar vertices[5*3] = {
        -2.0,  0.0,  0.0,
        +0.0, -1.0,  0.0,
        +0.0,  1.0,  0.0,
        +0.0,  0.0,  2.0,
        +2.0,  0.0,  0.0
    };
    data->vertices = const_cast<PylithScalar*>(vertices);

    static const PylithInt cells[2*4] = {
        0,  1,  2,  3,
        1,  4,  2,  3,
    };
    data->cells = const_cast<PylithInt*>(cells);
    static const PylithInt materialIds[2] = {
        7, 8,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numGroups = 2;
    static const PylithInt groupSizes[2] = { 3,  4, };
    data->groupSizes = const_cast<PylithInt*>(groupSizes);
    static const PylithInt groups[3+4] = {
        1, 2, 3,
        0, 1, 2, 3,
    };
    data->groups = const_cast<PylithInt*>(groups);
    static const char* groupNames[2] = {
        "mid_face",
        "bottom_face",
    };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[2] = {
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<char**>(groupTypes);

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
pylith::meshio::TestMeshIOCubit_Cases::Hex_v12(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/twohex8_12.2.exo";

    data->numVertices = 12;
    data->spaceDim = 3;
    data->numCells = 2;
    data->cellDim = 3;
    data->numCorners = 8;

    static const PylithScalar vertices[12*3] = {
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
    data->vertices = const_cast<PylithScalar*>(vertices);

    static const PylithInt cells[2*8] = {
        0,  1,  2,  3,  4,  5,  6,  7,
        4,  5,  6,  7,  8,  9, 10, 11
    };
    data->cells = const_cast<PylithInt*>(cells);
    static const PylithInt materialIds[2] = {
        7, 8,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numGroups = 2;
    static const PylithInt groupSizes[2] = { 4,  6, };
    data->groupSizes = const_cast<PylithInt*>(groupSizes);
    static const PylithInt groups[4+6] = {
        8,  9, 10, 11,
        0,  3,  4,  7,  8, 11
    };
    data->groups = const_cast<PylithInt*>(groups);
    static const char* groupNames[2] = {
        "right_face",
        "top_face",
    };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[2] = {
        "vertex",
        "vertex",
    };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // Hex_v12


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOCubit_Cases::Hex_v13(void) {
    TestMeshIO_Data* data = Hex_v12();assert(data);

    data->filename = "data/twohex8_13.0.exo";

    return data;
} // Hex_v13


// End of file
