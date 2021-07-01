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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestBoundaryMesh.hh" // Implementation of cases

// ----------------------------------------------------------------------
namespace pylith {
    namespace bc {

        // --------------------------------------------------------------
        class TestBoundaryMesh_Tri : public TestBoundaryMesh {
            CPPUNIT_TEST_SUB_SUITE(TestBoundaryMesh_Tri, TestBoundaryMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestBoundaryMesh::setUp();
                _data = new TestBoundaryMesh_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "data/tri3.mesh";
                _data->bcLabel = "bc";
                _data->faultLabel = "fault";
                _data->faultId = 100;

                _data->numCorners = 2;
                _data->numCells = 1;
                _data->numVerticesNoFault = 2;
                _data->numVerticesWithFault = 2;
            } // setUp
        }; // class TestBoundaryMesh_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION(TestBoundaryMesh_Tri);


        // --------------------------------------------------------------
        class TestBoundaryMesh_Quad : public TestBoundaryMesh {
            CPPUNIT_TEST_SUB_SUITE(TestBoundaryMesh_Quad, TestBoundaryMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestBoundaryMesh::setUp();
                _data = new TestBoundaryMesh_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "data/quad4.mesh";
                _data->bcLabel = "bc3";
                _data->faultLabel = "fault";
                _data->faultId = 100;

                _data->numCorners = 2;
                _data->numCells = 2;
                _data->numVerticesNoFault = 3;
                _data->numVerticesWithFault = 4;
            } // setUp
        }; // class TestBoundaryMesh_Quad
        CPPUNIT_TEST_SUITE_REGISTRATION(TestBoundaryMesh_Quad);


        // --------------------------------------------------------------
        class TestBoundaryMesh_Tet : public TestBoundaryMesh {
            CPPUNIT_TEST_SUB_SUITE(TestBoundaryMesh_Tet, TestBoundaryMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestBoundaryMesh::setUp();
                _data = new TestBoundaryMesh_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "data/tet4.mesh";
                _data->bcLabel = "bc4";
                _data->faultLabel = "fault";
                _data->faultId = 100;

                _data->numCorners = 3;
                _data->numCells = 2;
                _data->numVerticesNoFault = 4;
                _data->numVerticesWithFault = 6;
            } // setUp
        }; // class TestBoundaryMesh_Tet
        CPPUNIT_TEST_SUITE_REGISTRATION(TestBoundaryMesh_Tet);


        // --------------------------------------------------------------
        class TestBoundaryMesh_Hex : public TestBoundaryMesh {
            CPPUNIT_TEST_SUB_SUITE(TestBoundaryMesh_Hex, TestBoundaryMesh);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestBoundaryMesh::setUp();
                _data = new TestBoundaryMesh_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "data/hex8.mesh";
                _data->bcLabel = "bc2";
                _data->faultLabel = "fault";
                _data->faultId = 100;

                _data->numCorners = 4;
                _data->numCells = 2;
                _data->numVerticesNoFault = 6;
                _data->numVerticesWithFault = 8;
            } // setUp
        }; // class TestBoundaryMesh_Hex
        CPPUNIT_TEST_SUITE_REGISTRATION(TestBoundaryMesh_Hex);


    } // namespace bc
} // namespace pylith


// End of file
