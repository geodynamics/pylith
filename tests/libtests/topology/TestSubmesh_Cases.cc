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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// -----------------------------------------------------------------------------
//

#include <portinfo>

#include "TestSubmesh.hh" // Implementation of class methods

namespace pylith {
    namespace topology {
        class TestSubmesh_Tri;
        class TestSubmesh_Quad;
        class TestSubmesh_Tet;
        class TestSubmesh_Hex;
    } // topology
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::topology::TestSubmesh_Tri : public pylith::topology::TestSubmesh {
    CPPUNIT_TEST_SUB_SUITE(TestSubmesh_Tri, TestSubmesh);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestSubmesh::setUp();

        _data->cellDim = 2;
        _data->numVertices = 4;
        _data->numCells = 2;
        _data->numCorners = 3;
        static const int _cells[2*3] = {
            0, 1, 3,
            0, 3, 2,
        };
        _data->cells = const_cast<int*>(_cells);
        static const PylithScalar _coordinates[4*2] = {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            1.0, 1.0,
        };
        _data->coordinates = const_cast<PylithScalar*>(_coordinates);

        // Submesh data
        _data->groupLabel = "bc";
        _data->groupSize = 3;
        static const int _groupVertices[3] = {
            1, 2, 3,
        };
        _data->groupVertices = const_cast<int*>(_groupVertices);

        _data->submeshNumCorners = 2;
        _data->submeshNumVertices = 3;
        static const int _submeshVertices[3] = {
            2, 3, 4,
        };
        _data->submeshVertices = const_cast<int*>(_submeshVertices);
        _data->submeshNumCells = 2;
        static const int _submeshCells[2] = {
            0, 1,
        };
        _data->submeshCells = const_cast<int*>(_submeshCells);

        // Subdomain data
        _data->subdomainLabel = "material-id";
        static const int _subdomainLabelValues[2] = {
            10, 20,
        };
        _data->subdomainLabelValues = const_cast<int*>(_subdomainLabelValues);
        _data->subdomainLabelValue = 10;
        _data->subdomainNumCorners = 3;
        _data->subdomainNumVertices = 3;
        static const int _subdomainVertices[3] = {
            1, 2, 3,
        };
        _data->subdomainVertices = const_cast<int*>(_subdomainVertices);
        _data->subdomainNumCells = 1;
        static const int _subdomainCells[1] = {
            0,
        };
        _data->subdomainCells = const_cast<int*>(_subdomainCells);
    } // setUp

}; // _TestSubmesh_Tri
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestSubmesh_Tri);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::topology::TestSubmesh_Quad : public pylith::topology::TestSubmesh {
    CPPUNIT_TEST_SUB_SUITE(TestSubmesh_Quad, TestSubmesh);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestSubmesh::setUp();

        _data->cellDim = 2;
        _data->numVertices = 6;
        _data->numCells = 2;
        _data->numCorners = 4;
        static const int _cells[2*4] = {
            0, 2, 3, 1,
            2, 4, 5, 3,
        };
        _data->cells = const_cast<int*>(_cells);
        static const PylithScalar _coordinates[6*2] = {
            -1.0, -1.0,
            -1.0, +1.0,
            +0.0, -1.0,
            +0.0, +1.0,
            +1.0, -1.0,
            +1.0, +1.0,
        };
        _data->coordinates = const_cast<PylithScalar*>(_coordinates);

        // Submesh data
        _data->groupLabel = "bc";
        _data->groupSize = 3;
        static const int _groupVertices[3] = {
            0, 2, 4,
        };
        _data->groupVertices = const_cast<int*>(_groupVertices);
        _data->submeshNumCorners = 2;
        _data->submeshNumVertices = 3;
        static const int _submeshVertices[3] = {
            2, 3, 4,
        };
        _data->submeshVertices = const_cast<int*>(_submeshVertices);
        _data->submeshNumCells = 2;
        static const int _submeshCells[2] = {
            0, 1,
        };
        _data->submeshCells = const_cast<int*>(_submeshCells);

        // Subdomain data
        _data->subdomainLabel = "material-id";
        static const int _subdomainLabelValues[2] = {
            10, 20,
        };
        _data->subdomainLabelValues = const_cast<int*>(_subdomainLabelValues);
        _data->subdomainLabelValue = 10;
        _data->subdomainNumCorners = 4;
        _data->subdomainNumVertices = 4;
        static const int _subdomainVertices[4] = {
            1, 2, 3, 4,
        };
        _data->subdomainVertices = const_cast<int*>(_subdomainVertices);
        _data->subdomainNumCells = 1;
        static const int _subdomainCells[1] = {
            0,
        };
        _data->subdomainCells = const_cast<int*>(_subdomainCells);
    } // setUp

}; // _TestSubmesh_Quad
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestSubmesh_Quad);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::topology::TestSubmesh_Tet : public pylith::topology::TestSubmesh {
    CPPUNIT_TEST_SUB_SUITE(TestSubmesh_Tet, TestSubmesh);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestSubmesh::setUp();

        _data->cellDim = 3;
        _data->numVertices = 5;
        _data->numCells = 2;
        _data->numCorners = 4;
        static const int _cells[2*4] = {
            1, 2, 3, 0,
            1, 3, 2, 4,
        };
        _data->cells = const_cast<int*>(_cells);
        static const PylithScalar _coordinates[5*3] = {
            -1.0, +0.0, +0.0,
            +0.0, -1.0, +0.0,
            +0.0, +0.0, +1.0,
            +0.0, +1.0, +0.0,
            +1.0, +0.0, +0.0,
        };
        _data->coordinates = const_cast<PylithScalar*>(_coordinates);

        // Submesh data
        _data->groupLabel = "bc";
        _data->groupSize = 4;
        static const int _groupVertices[4] = {
            0, 1, 3, 4,
        };
        _data->groupVertices = const_cast<int*>(_groupVertices);
        _data->submeshNumCorners = 3;
        _data->submeshNumVertices = 4;
        static const int _submeshVertices[4] = {
            2, 3, 4, 5,
        };
        _data->submeshVertices = const_cast<int*>(_submeshVertices);
        _data->submeshNumCells = 2;
        static const int _submeshCells[2] = {
            0, 1,
        };
        _data->submeshCells = const_cast<int*>(_submeshCells);

        // Subdomain data
        _data->subdomainLabel = "subdomain-id";
        static const int _subdomainLabelValues[2] = {
            20, 10,
        };
        _data->subdomainLabelValues = const_cast<int*>(_subdomainLabelValues);
        _data->subdomainLabelValue = 10;
        _data->subdomainNumCorners = 4;
        _data->subdomainNumVertices = 4;
        static const int _subdomainVertices[4] = {
            1, 2, 3, 4,
        };
        _data->subdomainVertices = const_cast<int*>(_subdomainVertices);
        _data->subdomainNumCells = 1;
        static const int _subdomainCells[1] = {
            0,
        };
        _data->subdomainCells = const_cast<int*>(_subdomainCells);
    } // setUp

}; // _TestSubmesh_Tet
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestSubmesh_Tet);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::topology::TestSubmesh_Hex : public pylith::topology::TestSubmesh {
    CPPUNIT_TEST_SUB_SUITE(TestSubmesh_Hex, TestSubmesh);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestSubmesh::setUp();

        _data->cellDim = 3;
        _data->numVertices = 12;
        _data->numCells = 2;
        _data->numCorners = 8;
        static const int _cells[2*8] = {
            0,  2,  3,  1,  6,  8,  9,  7,
            2,  4,  5,  3,  8, 10, 11,  9,
        };
        _data->cells = const_cast<int*>(_cells);
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
        _data->coordinates = const_cast<PylithScalar*>(_coordinates);

        // Submesh data
        _data->groupLabel = "bc";
        _data->groupSize = 6;
        static const int _groupVertices[6] = {
            6, 7, 8, 9, 10, 11,
        };
        _data->groupVertices = const_cast<int*>(_groupVertices);
        _data->submeshNumCorners = 4;
        _data->submeshNumVertices = 6;
        static const int _submeshVertices[6] = {
            2, 3, 4, 5, 6, 7,
        };
        _data->submeshVertices = const_cast<int*>(_submeshVertices);
        _data->submeshNumCells = 2;
        static const int _submeshCells[2] = {
            0, 1,
        };
        _data->submeshCells = const_cast<int*>(_submeshCells);

        // Subdomain data
        _data->subdomainLabel = "sub-id";
        static const int _subdomainLabelValues[2] = {
            20, 10,
        };
        _data->subdomainLabelValues = const_cast<int*>(_subdomainLabelValues);
        _data->subdomainLabelValue = 10;
        _data->subdomainNumCorners = 4;
        _data->subdomainNumVertices = 8;
        static const int _subdomainVertices[8] = {
            1, 2, 3, 4, 5, 6, 7, 8,
        };
        _data->subdomainVertices = const_cast<int*>(_subdomainVertices);
        _data->subdomainNumCells = 1;
        static const int _subdomainCells[1] = {
            0,
        };
        _data->subdomainCells = const_cast<int*>(_subdomainCells);
    } // setUp

}; // _TestSubmesh_Hex
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestSubmesh_Hex);

// End of file
