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

#include "TestSubmesh.hh" // Implementation of class methods

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace topology {
        class TestSubmesh_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestSubmesh_Cases {
public:

    // Data factory methods
    static TestSubmesh_Data* Tri(void);

    static TestSubmesh_Data* Quad(void);

    static TestSubmesh_Data* Tet(void);

    static TestSubmesh_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestSubmesh::Tri::testAccessors", "[TestSubmesh][Tri][testAccessors]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Tri()).testAccessors();
}
TEST_CASE("TestSubmesh::Tri::testSizes", "[TestSubmesh][Tri][testSizes]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Tri()).testSizes();
}
TEST_CASE("TestSubmesh::Tri::testCreateLowerDimMesh", "[TestSubmesh][Tri][testCreateLowerDimMesh]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Tri()).testCreateLowerDimMesh();
}
TEST_CASE("TestSubmesh::Tri::testCreateSubdomainMesh", "[TestSubmesh][Tri][testCreateSubdomainMesh]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Tri()).testCreateSubdomainMesh();
}

TEST_CASE("TestSubmesh::Quad::testAccessors", "[TestSubmesh][Quad][testAccessors]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Quad()).testAccessors();
}
TEST_CASE("TestSubmesh::Quad::testSizes", "[TestSubmesh][Quad][testSizes]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Quad()).testSizes();
}
TEST_CASE("TestSubmesh::Quad::testCreateLowerDimMesh", "[TestSubmesh][Quad][testCreateLowerDimMesh]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Quad()).testCreateLowerDimMesh();
}
TEST_CASE("TestSubmesh::Quad::testCreateSubdomainMesh", "[TestSubmesh][Quad][testCreateSubdomainMesh]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Quad()).testCreateSubdomainMesh();
}

TEST_CASE("TestSubmesh::Tet::testAccessors", "[TestSubmesh][Tet][testAccessors]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Tet()).testAccessors();
}
TEST_CASE("TestSubmesh::Tet::testSizes", "[TestSubmesh][Tet][testSizes]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Tet()).testSizes();
}
TEST_CASE("TestSubmesh::Tet::testCreateLowerDimMesh", "[TestSubmesh][Tet][testCreateLowerDimMesh]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Tet()).testCreateLowerDimMesh();
}
TEST_CASE("TestSubmesh::Tet::testCreateSubdomainMesh", "[TestSubmesh][Tet][testCreateSubdomainMesh]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Tet()).testCreateSubdomainMesh();
}

TEST_CASE("TestSubmesh::Hex::testAccessors", "[TestSubmesh][Hex][testAccessors]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Hex()).testAccessors();
}
TEST_CASE("TestSubmesh::Hex::testSizes", "[TestSubmesh][Hex][testSizes]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Hex()).testSizes();
}
TEST_CASE("TestSubmesh::Hex::testCreateLowerDimMesh", "[TestSubmesh][Hex][testCreateLowerDimMesh]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Hex()).testCreateLowerDimMesh();
}
TEST_CASE("TestSubmesh::Hex::testCreateSubdomainMesh", "[TestSubmesh][Hex][testCreateSubdomainMesh]") {
    pylith::topology::TestSubmesh(pylith::topology::TestSubmesh_Cases::Hex()).testCreateSubdomainMesh();
}

// ------------------------------------------------------------------------------------------------
pylith::topology::TestSubmesh_Data*
pylith::topology::TestSubmesh_Cases::Tri(void) {
    TestSubmesh_Data* data = new TestSubmesh_Data();assert(data);

    data->cellDim = 2;
    data->numVertices = 4;
    data->numCells = 2;
    data->numCorners = 3;
    static const int _cells[2*3] = {
        0, 1, 3,
        0, 3, 2,
    };
    data->cells = const_cast<int*>(_cells);
    static const PylithScalar _coordinates[4*2] = {
        0.0, 0.0,
        1.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
    };
    data->coordinates = const_cast<PylithScalar*>(_coordinates);

    // Submesh data
    data->groupLabel = "bc";
    data->groupSize = 3;
    static const int _groupVertices[3] = {
        1, 2, 3,
    };
    data->groupVertices = const_cast<int*>(_groupVertices);

    data->submeshNumCorners = 2;
    data->submeshNumVertices = 3;
    static const int _submeshVertices[3] = {
        2, 3, 4,
    };
    data->submeshVertices = const_cast<int*>(_submeshVertices);
    data->submeshNumCells = 2;
    static const int _submeshCells[2] = {
        0, 1,
    };
    data->submeshCells = const_cast<int*>(_submeshCells);

    // Subdomain data
    data->subdomainLabel = "material-id";
    static const int _subdomainLabelValues[2] = {
        10, 20,
    };
    data->subdomainLabelValues = const_cast<int*>(_subdomainLabelValues);
    data->subdomainLabelValue = 10;
    data->subdomainNumCorners = 3;
    data->subdomainNumVertices = 3;
    static const int _subdomainVertices[3] = {
        1, 2, 3,
    };
    data->subdomainVertices = const_cast<int*>(_subdomainVertices);
    data->subdomainNumCells = 1;
    static const int _subdomainCells[1] = {
        0,
    };
    data->subdomainCells = const_cast<int*>(_subdomainCells);

    return data;
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::topology::TestSubmesh_Data*
pylith::topology::TestSubmesh_Cases::Quad(void) {
    TestSubmesh_Data* data = new TestSubmesh_Data();assert(data);

    data->cellDim = 2;
    data->numVertices = 6;
    data->numCells = 2;
    data->numCorners = 4;
    static const int _cells[2*4] = {
        0, 2, 3, 1,
        2, 4, 5, 3,
    };
    data->cells = const_cast<int*>(_cells);
    static const PylithScalar _coordinates[6*2] = {
        -1.0, -1.0,
        -1.0, +1.0,
        +0.0, -1.0,
        +0.0, +1.0,
        +1.0, -1.0,
        +1.0, +1.0,
    };
    data->coordinates = const_cast<PylithScalar*>(_coordinates);

    // Submesh data
    data->groupLabel = "bc";
    data->groupSize = 3;
    static const int _groupVertices[3] = {
        0, 2, 4,
    };
    data->groupVertices = const_cast<int*>(_groupVertices);
    data->submeshNumCorners = 2;
    data->submeshNumVertices = 3;
    static const int _submeshVertices[3] = {
        2, 3, 4,
    };
    data->submeshVertices = const_cast<int*>(_submeshVertices);
    data->submeshNumCells = 2;
    static const int _submeshCells[2] = {
        0, 1,
    };
    data->submeshCells = const_cast<int*>(_submeshCells);

    // Subdomain data
    data->subdomainLabel = "material-id";
    static const int _subdomainLabelValues[2] = {
        10, 20,
    };
    data->subdomainLabelValues = const_cast<int*>(_subdomainLabelValues);
    data->subdomainLabelValue = 10;
    data->subdomainNumCorners = 4;
    data->subdomainNumVertices = 4;
    static const int _subdomainVertices[4] = {
        1, 2, 3, 4,
    };
    data->subdomainVertices = const_cast<int*>(_subdomainVertices);
    data->subdomainNumCells = 1;
    static const int _subdomainCells[1] = {
        0,
    };
    data->subdomainCells = const_cast<int*>(_subdomainCells);

    return data;
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::topology::TestSubmesh_Data*
pylith::topology::TestSubmesh_Cases::Tet(void) {
    TestSubmesh_Data* data = new TestSubmesh_Data();assert(data);

    data->cellDim = 3;
    data->numVertices = 5;
    data->numCells = 2;
    data->numCorners = 4;
    static const int _cells[2*4] = {
        1, 2, 3, 0,
        1, 3, 2, 4,
    };
    data->cells = const_cast<int*>(_cells);
    static const PylithScalar _coordinates[5*3] = {
        -1.0, +0.0, +0.0,
        +0.0, -1.0, +0.0,
        +0.0, +0.0, +1.0,
        +0.0, +1.0, +0.0,
        +1.0, +0.0, +0.0,
    };
    data->coordinates = const_cast<PylithScalar*>(_coordinates);

    // Submesh data
    data->groupLabel = "bc";
    data->groupSize = 4;
    static const int _groupVertices[4] = {
        0, 1, 3, 4,
    };
    data->groupVertices = const_cast<int*>(_groupVertices);
    data->submeshNumCorners = 3;
    data->submeshNumVertices = 4;
    static const int _submeshVertices[4] = {
        2, 3, 4, 5,
    };
    data->submeshVertices = const_cast<int*>(_submeshVertices);
    data->submeshNumCells = 2;
    static const int _submeshCells[2] = {
        0, 1,
    };
    data->submeshCells = const_cast<int*>(_submeshCells);

    // Subdomain data
    data->subdomainLabel = "subdomain-id";
    static const int _subdomainLabelValues[2] = {
        20, 10,
    };
    data->subdomainLabelValues = const_cast<int*>(_subdomainLabelValues);
    data->subdomainLabelValue = 10;
    data->subdomainNumCorners = 4;
    data->subdomainNumVertices = 4;
    static const int _subdomainVertices[4] = {
        1, 2, 3, 4,
    };
    data->subdomainVertices = const_cast<int*>(_subdomainVertices);
    data->subdomainNumCells = 1;
    static const int _subdomainCells[1] = {
        0,
    };
    data->subdomainCells = const_cast<int*>(_subdomainCells);

    return data;
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::topology::TestSubmesh_Data*
pylith::topology::TestSubmesh_Cases::Hex(void) {
    TestSubmesh_Data* data = new TestSubmesh_Data();assert(data);

    data->cellDim = 3;
    data->numVertices = 12;
    data->numCells = 2;
    data->numCorners = 8;
    static const int _cells[2*8] = {
        0,  2,  3,  1,  6,  8,  9,  7,
        2,  4,  5,  3,  8, 10, 11,  9,
    };
    data->cells = const_cast<int*>(_cells);
    static const PylithScalar _coordinates[12*3] = {
        -1.0, -1.0, -1.0,
        -1.0, +1.0, -1.0,
        +0.0, -1.0, -1.0,
        +0.0,  1.0, -1.0,
        +1.0, -1.0, -1.0,
        +1.0,  1.0, -1.0,
        -1.0, -1.0, +1.0,
        -1.0,  1.0, +1.0,
        +0.0, -1.0, +1.0,
        +0.0, +1.0, +1.0,
        +1.0, -1.0, +1.0,
        +1.0, +1.0, +1.0,
    };
    data->coordinates = const_cast<PylithScalar*>(_coordinates);

    // Submesh data
    data->groupLabel = "bc";
    data->groupSize = 6;
    static const int _groupVertices[6] = {
        6, 7, 8, 9, 10, 11,
    };
    data->groupVertices = const_cast<int*>(_groupVertices);
    data->submeshNumCorners = 4;
    data->submeshNumVertices = 6;
    static const int _submeshVertices[6] = {
        2, 3, 4, 5, 6, 7,
    };
    data->submeshVertices = const_cast<int*>(_submeshVertices);
    data->submeshNumCells = 2;
    static const int _submeshCells[2] = {
        0, 1,
    };
    data->submeshCells = const_cast<int*>(_submeshCells);

    // Subdomain data
    data->subdomainLabel = "sub-id";
    static const int _subdomainLabelValues[2] = {
        20, 10,
    };
    data->subdomainLabelValues = const_cast<int*>(_subdomainLabelValues);
    data->subdomainLabelValue = 10;
    data->subdomainNumCorners = 4;
    data->subdomainNumVertices = 8;
    static const int _subdomainVertices[8] = {
        1, 2, 3, 4, 5, 6, 7, 8,
    };
    data->subdomainVertices = const_cast<int*>(_subdomainVertices);
    data->subdomainNumCells = 1;
    static const int _subdomainCells[1] = {
        0,
    };
    data->subdomainCells = const_cast<int*>(_subdomainCells);

    return data;
} // Hex


// End of file
