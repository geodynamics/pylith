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

#include "TestMeshIOCubit.hh" // Implementation of class methods

namespace pylith {
    namespace meshio {

        // --------------------------------------------------------------
        class TestMeshIOCubit_Tri : public TestMeshIOCubit {
public:
            void setUp(void) {
                TestMeshIOCubit::setUp();
                _data = new TestMeshIOCubit_Data();CPPUNIT_ASSERT(_data);
                _data->numVertices = 4;
                _data->spaceDim = 2;
                _data->numCells = 2;
                _data->cellDim = 2;
                _data->numCorners = 3;

                static const PylithScalar vertices[4*2] = {
                    -1.0,  +0.0,
                    +0.0,  -1.0,
                    +0.0,  +1.0,
                    +1.0,  +0.0
                };
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[2*3] = {
                    0,  1,  2,
                    2,  1,  3,
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[2] = {
                    2, 3,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 2;
                static const PylithInt groupSizes[2] = { 1,  2, };
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[1+2] = {
                    0,
                    2, 3,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[2] = {
                    "left_vertex",
                    "right_vertex",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[2] = {
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp
        }; // class TestMeshIOCubit_Tri


        // --------------------------------------------------------------
        class TestMeshIOCubit_Tri_v12 : public TestMeshIOCubit_Tri {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOCubit_Tri_v12, TestMeshIOCubit_Tri);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOCubit_Tri::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/twotri3_12.2.exo";
            } // setUp
        }; // class TestMeshIOCubit_Tri_v12
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOCubit_Tri_v12);

        // --------------------------------------------------------------
        class TestMeshIOCubit_Tri_v13 : public TestMeshIOCubit_Tri {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOCubit_Tri_v13, TestMeshIOCubit_Tri);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOCubit_Tri::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/twotri3_13.0.exo";
            } // setUp
        }; // class TestMeshIOCubit_Tri_v13
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOCubit_Tri_v13);

        // --------------------------------------------------------------
        class TestMeshIOCubit_Quad : public TestMeshIOCubit {
public:
            void setUp(void) {
                TestMeshIOCubit::setUp();
                _data = new TestMeshIOCubit_Data();CPPUNIT_ASSERT(_data);
                _data->numVertices = 6;
                _data->spaceDim = 2;
                _data->numCells = 2;
                _data->cellDim = 2;
                _data->numCorners = 4;

                static const PylithScalar vertices[6*2] = {
                    0.0,  0.0,
                    1.0,  0.0,
                    1.0,  1.0,
                    0.0,  1.0,
                    2.0,  0.0,
                    2.0,  1.0,
                };
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[2*4] = {
                    0,  1,  2,  3,
                    1,  4,  5,  2,
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[2] = {
                    10, 11,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 2;
                static const PylithInt groupSizes[2] = { 2,  3, };
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[2+3] = {
                    0, 3,
                    2, 3, 5,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[2] = {
                    "left_edge",
                    "top_edge",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[2] = {
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp
        }; // class TestMeshIOCubit_Quad


        // --------------------------------------------------------------
        class TestMeshIOCubit_Quad_v12 : public TestMeshIOCubit_Quad {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOCubit_Quad_v12, TestMeshIOCubit_Quad);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOCubit_Quad::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/twoquad4_12.2.exo";
            } // setUp
        }; // class TestMeshIOCubit_Quad_v12
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOCubit_Quad_v12);

        // --------------------------------------------------------------
        class TestMeshIOCubit_Quad_v13 : public TestMeshIOCubit_Quad {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOCubit_Quad_v13, TestMeshIOCubit_Quad);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOCubit_Quad::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/twoquad4_13.0.exo";
            } // setUp
        }; // class TestMeshIOCubit_Quad_v13
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOCubit_Quad_v13);

        // --------------------------------------------------------------
        class TestMeshIOCubit_Tet : public TestMeshIOCubit {
public:
            void setUp(void) {
                TestMeshIOCubit::setUp();
                _data = new TestMeshIOCubit_Data();CPPUNIT_ASSERT(_data);
                _data->numVertices = 5;
                _data->spaceDim = 3;
                _data->numCells = 2;
                _data->cellDim = 3;
                _data->numCorners = 4;

                static const PylithScalar vertices[5*3] = {
                    -2.0,  0.0,  0.0,
                    +0.0, -1.0,  0.0,
                    +0.0,  1.0,  0.0,
                    +0.0,  0.0,  2.0,
                    +2.0,  0.0,  0.0
                };
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[2*4] = {
                    0,  1,  2,  3,
                    1,  4,  2,  3,
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[2] = {
                    7, 8,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 2;
                static const PylithInt groupSizes[2] = { 3,  4, };
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[3+4] = {
                    1, 2, 3,
                    0, 1, 2, 3,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[2] = {
                    "mid_face",
                    "bottom_face",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[2] = {
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp
        }; // class TestMeshIOCubit_Tet


        // --------------------------------------------------------------
        class TestMeshIOCubit_Tet_v12 : public TestMeshIOCubit_Tet {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOCubit_Tet_v12, TestMeshIOCubit_Tet);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOCubit_Tet::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/twotet4_12.2.exo";
            } // setUp
        }; // class TestMeshIOCubit_Tet_v12
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOCubit_Tet_v12);

        // --------------------------------------------------------------
        class TestMeshIOCubit_Tet_v13 : public TestMeshIOCubit_Tet {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOCubit_Tet_v13, TestMeshIOCubit_Tet);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOCubit_Tet::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/twotet4_13.0.exo";
            } // setUp
        }; // class TestMeshIOCubit_Tet_v13
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOCubit_Tet_v13);

        // --------------------------------------------------------------
        class TestMeshIOCubit_Hex : public TestMeshIOCubit {
public:
            void setUp(void) {
                TestMeshIOCubit::setUp();
                _data = new TestMeshIOCubit_Data();CPPUNIT_ASSERT(_data);
                _data->numVertices = 12;
                _data->spaceDim = 3;
                _data->numCells = 2;
                _data->cellDim = 3;
                _data->numCorners = 8;

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
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[2*8] = {
                    0,  1,  2,  3,  4,  5,  6,  7,
                    4,  5,  6,  7,  8,  9, 10, 11
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[2] = {
                    7, 8,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 2;
                static const PylithInt groupSizes[2] = { 4,  6, };
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[4+6] = {
                    8,  9, 10, 11,
                    0,  3,  4,  7,  8, 11
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[2] = {
                    "right_face",
                    "top_face",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[2] = {
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp
        }; // class TestMeshIOCubit_Hex


        // --------------------------------------------------------------
        class TestMeshIOCubit_Hex_v12 : public TestMeshIOCubit_Hex {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOCubit_Hex_v12, TestMeshIOCubit_Hex);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOCubit_Hex::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/twohex8_12.2.exo";
            } // setUp
        }; // class TestMeshIOCubit_Hex_v12
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOCubit_Hex_v12);

        // --------------------------------------------------------------
        class TestMeshIOCubit_Hex_v13 : public TestMeshIOCubit_Hex {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOCubit_Hex_v13, TestMeshIOCubit_Hex);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOCubit_Hex::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/twohex8_13.0.exo";
            } // setUp
        }; // class TestMeshIOCubit_Hex_v13
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOCubit_Hex_v13);

    } // meshio
} // pylith


// End of file
