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

#include "TestMeshIOAscii.hh" // Implementation of class methods

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

    static TestMeshIO_Data* Quad3D(void);

    static TestMeshIO_Data* Hex3D(void);

    static TestMeshIO_Data* Hex3D_Index1(void);

}; // TestMeshIOAscii_Cases

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshIOAscii::Quad2D::testFilename", "[TestMeshIOAscii][testFilename]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Quad2D()).testFilename();
}

TEST_CASE("TestMeshIOAscii::Quad2D::testWriteRead", "[TestMeshIOAscii][testWriteRead]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Quad2D()).testWriteRead();
}
TEST_CASE("TestMeshIOAscii::Quad3D::testWriteRead", "[TestMeshIOAscii][testWriteRead]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Quad3D()).testWriteRead();
}
TEST_CASE("TestMeshIOAscii::Hex3D::testWriteRead", "[TestMeshIOAscii][testWriteRead]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Hex3D()).testWriteRead();
}

TEST_CASE("TestMeshIOAscii::Quad2D_Comments::testRead", "[TestMeshIOAscii][testRead]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Quad2D_Comments()).testRead();
}
TEST_CASE("TestMeshIOAscii::Hex3D_Index1::testReas", "[TestMeshIOAscii][testRead]") {
    pylith::meshio::TestMeshIOAscii(pylith::meshio::TestMeshIOAscii_Cases::Hex3D_Index1()).testRead();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii_Cases::Quad2D(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "mesh2D.txt";
    data->numVertices = 9;
    data->spaceDim = 2;
    data->numCells = 3;
    data->cellDim = 2;
    data->numCorners = 4;

    static const PylithScalar vertices[9*2] = {
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
    data->vertices = const_cast<PylithScalar*>(vertices);

    static const PylithInt cells[3*4] = {
        0,  2,  3,  1,
        4,  3,  6,  5,
        3,  7,  8,  6,
    };
    data->cells = const_cast<PylithInt*>(cells);
    static const PylithInt materialIds[3] = {
        1, 0, 1,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numGroups = 3;
    static const PylithInt groupSizes[3] = { 5, 3, 2, };
    data->groupSizes = const_cast<PylithInt*>(groupSizes);
    static const PylithInt groups[5+3+2] = {
        0, 2, 4, 6, 8,
        1, 4, 7,
        0, 2,
    };
    data->groups = const_cast<PylithInt*>(groups);
    static const char* groupNames[3] = {
        "group A",
        "group B",
        "group C",
    };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[3] = {
        "vertex",
        "vertex",
        "cell",
    };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // Quad2D


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii_Cases::Quad2D_Comments(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/mesh2D_comments.txt";
    data->numVertices = 9;
    data->spaceDim = 2;
    data->numCells = 3;
    data->cellDim = 2;
    data->numCorners = 4;

    static const PylithScalar vertices[9*2] = {
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
    data->vertices = const_cast<PylithScalar*>(vertices);

    static const PylithInt cells[3*4] = {
        0,  2,  3,  1,
        4,  3,  6,  5,
        3,  7,  8,  6,
    };
    data->cells = const_cast<PylithInt*>(cells);
    static const PylithInt materialIds[3] = {
        1, 0, 1,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numGroups = 3;
    static const PylithInt groupSizes[3] = { 5, 3, 2, };
    data->groupSizes = const_cast<PylithInt*>(groupSizes);
    static const PylithInt groups[5+3+2] = {
        0, 2, 4, 6, 8,
        1, 4, 7,
        0, 2,
    };
    data->groups = const_cast<PylithInt*>(groups);
    static const char* groupNames[3] = {
        "group A",
        "group B",
        "group C",
    };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[3] = {
        "vertex",
        "vertex",
        "cell",
    };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // Quad2D_Comments


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii_Cases::Quad3D(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "mesh2Din3D.txt";
    data->numVertices = 9;
    data->spaceDim = 3;
    data->numCells = 3;
    data->cellDim = 2;
    data->numCorners = 4;

    static const PylithScalar vertices[9*3] = {
        -1.0, +3.0, +0.2,
        +1.0, +3.3, +0.5,
        -1.2, +0.9, +0.3,
        +0.9, +1.0, +0.4,
        +3.0, +2.9, -0.1,
        +6.0, +1.2, -0.2,
        +3.4, -0.2, +0.1,
        +0.1, -1.1, +0.9,
        +2.9, -3.1, +0.8
    };
    data->vertices = const_cast<PylithScalar*>(vertices);

    static const PylithInt cells[3*4] = {
        0,  2,  3,  1,
        4,  3,  6,  5,
        3,  7,  8,  6,
    };
    data->cells = const_cast<PylithInt*>(cells);
    static const PylithInt materialIds[3] = {
        0, 1, 0,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numGroups = 1;
    static const PylithInt groupSizes[3] = { 3, };
    data->groupSizes = const_cast<PylithInt*>(groupSizes);
    static const PylithInt groups[3] = {
        0, 3, 6,
    };
    data->groups = const_cast<PylithInt*>(groups);
    static const char* groupNames[1] = {
        "group A",
    };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[1] = {
        "vertex",
    };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // Quad3D


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii_Cases::Hex3D(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "mesh3D.txt";
    data->numVertices = 14;
    data->spaceDim = 3;
    data->numCells = 2;
    data->cellDim = 3;
    data->numCorners = 8;

    static const PylithScalar vertices[14*3] = {
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
    data->vertices = const_cast<PylithScalar*>(vertices);

    static const PylithInt cells[2*8] = {
        6, 12, 13, 11,  7,  9, 10,  8,
        0,  2,  6,  5,  1,  3,  7,  4
    };
    data->cells = const_cast<PylithInt*>(cells);
    static const PylithInt materialIds[2] = {
        1, 0,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numGroups = 3;
    static const PylithInt groupSizes[3] = { 5, 2, 4,};
    data->groupSizes = const_cast<PylithInt*>(groupSizes);
    static const PylithInt groups[5+2+4] = {
        0, 4, 6, 7, 10,
        0, 1,
        0, 4, 12, 13
    };
    data->groups = const_cast<PylithInt*>(groups);
    static const char* groupNames[3] = {
        "group A",
        "group B",
        "group C",
    };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[3] = {
        "vertex",
        "cell",
        "vertex",
    };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // Hex3D


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOAscii_Cases::Hex3D_Index1(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/mesh3D_index1.txt";
    data->numVertices = 14;
    data->spaceDim = 3;
    data->numCells = 2;
    data->cellDim = 3;
    data->numCorners = 8;

    static const PylithScalar vertices[14*3] = {
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
    data->vertices = const_cast<PylithScalar*>(vertices);

    static const PylithInt cells[2*8] = {
        6, 12, 13, 11,  7,  9, 10,  8,
        0,  2,  6,  5,  1,  3,  7,  4
    };
    data->cells = const_cast<PylithInt*>(cells);
    static const PylithInt materialIds[2] = {
        2, 1,
    };
    data->materialIds = const_cast<PylithInt*>(materialIds);

    data->numGroups = 2;
    static const PylithInt groupSizes[2] = { 5, 2, };
    data->groupSizes = const_cast<PylithInt*>(groupSizes);
    static const PylithInt groups[5+2] = {
        0, 4, 6, 7, 10,
        0, 1,
    };
    data->groups = const_cast<PylithInt*>(groups);
    static const char* groupNames[3] = {
        "group A",
        "group B",
    };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[3] = {
        "vertex",
        "cell",
    };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // Hex3D_Index1


// End of file
