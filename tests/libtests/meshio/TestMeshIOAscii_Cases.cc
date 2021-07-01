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

#include "TestMeshIOAscii.hh" // Implementation of class methods

// ----------------------------------------------------------------------
namespace pylith {
    namespace meshio {
        // --------------------------------------------------------------
        class TestMeshIOAscii_Line1D : public TestMeshIOAscii {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOAscii_Line1D, TestMeshIOAscii);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOAscii::setUp();
                _data = new TestMeshIOAscii_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "mesh1D.txt";
                _data->numVertices = 3;
                _data->spaceDim = 1;
                _data->numCells = 2;
                _data->cellDim = 1;
                _data->numCorners = 2;

                static const PylithScalar vertices[3] = {
                    -1.2,
                    +2.1,
                    +0.3,
                };
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[2*2] = {
                    0,  2,
                    2,  1,
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[2] = {
                    2, 1,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 3;
                static const PylithInt groupSizes[3] = { 1,  1,  2, };
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[1+1+2] = {
                    1,
                    0,
                    0, 1,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[3] = {
                    "group A",
                    "group B",
                    "group C",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[3] = {
                    "vertex",
                    "cell",
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOAscii_Line1D
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOAscii_Line1D);

        // --------------------------------------------------------------
        class TestMeshIOAscii_Line2D : public TestMeshIOAscii {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOAscii_Line2D, TestMeshIOAscii);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOAscii::setUp();
                _data = new TestMeshIOAscii_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "mesh1Din2D.txt";
                _data->numVertices = 4;
                _data->spaceDim = 2;
                _data->numCells = 3;
                _data->cellDim = 1;
                _data->numCorners = 2;

                static const PylithScalar vertices[4*2] = {
                    -3.0, -1.2,
                    +1.0, -1.0,
                    +2.6,  3.1,
                    +1.8, -4.0
                };
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[3*2] = {
                    3,  1,
                    0,  1,
                    1,  2,
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[3] = {
                    1, 0, 1,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 2;
                static const PylithInt groupSizes[3] = { 2, 3, };
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[2+3] = {
                    0, 2,
                    0, 1, 3,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[2] = {
                    "group A",
                    "group B",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[2] = {
                    "cell",
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOAscii_Line2D
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOAscii_Line2D);

        // --------------------------------------------------------------
        class TestMeshIOAscii_Line3D : public TestMeshIOAscii {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOAscii_Line3D, TestMeshIOAscii);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOAscii::setUp();
                _data = new TestMeshIOAscii_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "mesh1Din3D.txt";
                _data->numVertices = 4;
                _data->spaceDim = 3;
                _data->numCells = 3;
                _data->cellDim = 1;
                _data->numCorners = 2;

                static const PylithScalar vertices[4*3] = {
                    -3.0, -1.2, +0.3,
                    +1.0, -1.0, +0.0,
                    +2.6, +3.1, -0.5,
                    +1.8, -4.0, +1.0
                };
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[3*2] = {
                    3,  1,
                    0,  1,
                    1,  2,
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[3] = {
                    1, 1, 0,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 2;
                static const PylithInt groupSizes[2] = { 1, 1, };
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[1+1] = {
                    2,
                    1,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[2] = {
                    "group A",
                    "group B",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[2] = {
                    "vertex",
                    "cell",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOAscii_Line3D
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOAscii_Line3D);

        // --------------------------------------------------------------
        class TestMeshIOAscii_Quad2D : public TestMeshIOAscii {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOAscii_Quad2D, TestMeshIOAscii);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOAscii::setUp();
                _data = new TestMeshIOAscii_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "mesh2D.txt";
                _data->numVertices = 9;
                _data->spaceDim = 2;
                _data->numCells = 3;
                _data->cellDim = 2;
                _data->numCorners = 4;

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
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[3*4] = {
                    0,  2,  3,  1,
                    4,  3,  6,  5,
                    3,  7,  8,  6,
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[3] = {
                    1, 0, 1,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 3;
                static const PylithInt groupSizes[3] = { 5, 3, 2, };
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[5+3+2] = {
                    0, 2, 4, 6, 8,
                    1, 4, 7,
                    0, 2,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[3] = {
                    "group A",
                    "group B",
                    "group C",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[3] = {
                    "vertex",
                    "vertex",
                    "cell",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOAscii_Quad2D
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOAscii_Quad2D);

        // --------------------------------------------------------------
        class TestMeshIOAscii_Quad2D_Comments : public TestMeshIOAscii {
            CPPUNIT_TEST_SUITE(TestMeshIOAscii_Quad2D_Comments);
            CPPUNIT_TEST(testRead);
            CPPUNIT_TEST_SUITE_END();

public:

            void setUp(void) {
                TestMeshIOAscii::setUp();
                _data = new TestMeshIOAscii_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "data/mesh2D_comments.txt";
                _data->numVertices = 9;
                _data->spaceDim = 2;
                _data->numCells = 3;
                _data->cellDim = 2;
                _data->numCorners = 4;

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
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[3*4] = {
                    0,  2,  3,  1,
                    4,  3,  6,  5,
                    3,  7,  8,  6,
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[3] = {
                    1, 0, 1,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 3;
                static const PylithInt groupSizes[3] = { 5, 3, 2, };
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[5+3+2] = {
                    0, 2, 4, 6, 8,
                    1, 4, 7,
                    0, 2,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[3] = {
                    "group A",
                    "group B",
                    "group C",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[3] = {
                    "vertex",
                    "vertex",
                    "cell",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOAscii_Quad2D_Comments
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOAscii_Quad2D_Comments);

        // --------------------------------------------------------------
        class TestMeshIOAscii_Quad3D : public TestMeshIOAscii {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOAscii_Quad3D, TestMeshIOAscii);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOAscii::setUp();
                _data = new TestMeshIOAscii_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "mesh2Din3D.txt";
                _data->numVertices = 9;
                _data->spaceDim = 3;
                _data->numCells = 3;
                _data->cellDim = 2;
                _data->numCorners = 4;

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
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[3*4] = {
                    0,  2,  3,  1,
                    4,  3,  6,  5,
                    3,  7,  8,  6,
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[3] = {
                    0, 1, 0,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 1;
                static const PylithInt groupSizes[3] = { 3, };
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[3] = {
                    0, 3, 6,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[1] = {
                    "group A",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[1] = {
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOAscii_Quad3D
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOAscii_Quad3D);

        // --------------------------------------------------------------
        class TestMeshIOAscii_Hex3D : public TestMeshIOAscii {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOAscii_Hex3D, TestMeshIOAscii);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOAscii::setUp();
                _data = new TestMeshIOAscii_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "mesh3D.txt";
                _data->numVertices = 14;
                _data->spaceDim = 3;
                _data->numCells = 2;
                _data->cellDim = 3;
                _data->numCorners = 8;

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
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[2*8] = {
                    6, 12, 13, 11,  7,  9, 10,  8,
                    0,  2,  6,  5,  1,  3,  7,  4
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[2] = {
                    1, 0,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 3;
                static const PylithInt groupSizes[3] = { 5, 2, 4,};
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[5+2+4] = {
                    0, 4, 6, 7, 10,
                    0, 1,
                    0, 4, 12, 13
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[3] = {
                    "group A",
                    "group B",
                    "group C",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[3] = {
                    "vertex",
                    "cell",
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOAscii_Hex3D
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOAscii_Hex3D);

        // --------------------------------------------------------------
        class TestMeshIOAscii_Hex3D_Index1 : public TestMeshIOAscii {
            CPPUNIT_TEST_SUITE(TestMeshIOAscii_Hex3D_Index1);
            CPPUNIT_TEST(testRead);
            CPPUNIT_TEST_SUITE_END();

public:

            void setUp(void) {
                TestMeshIOAscii::setUp();
                _data = new TestMeshIOAscii_Data();CPPUNIT_ASSERT(_data);
                _data->filename = "data/mesh3D_index1.txt";
                _data->numVertices = 14;
                _data->spaceDim = 3;
                _data->numCells = 2;
                _data->cellDim = 3;
                _data->numCorners = 8;

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
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[2*8] = {
                    6, 12, 13, 11,  7,  9, 10,  8,
                    0,  2,  6,  5,  1,  3,  7,  4
                };
                _data->cells = const_cast<PylithInt*>(cells);
                static const PylithInt materialIds[2] = {
                    2, 1,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 2;
                static const PylithInt groupSizes[2] = { 5, 2, };
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[5+2] = {
                    0, 4, 6, 7, 10,
                    0, 1,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[3] = {
                    "group A",
                    "group B",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[3] = {
                    "vertex",
                    "cell",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOAscii_Hex3D_Index1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOAscii_Hex3D_Index1);

    } // meshio
} // pylith

// End of file
