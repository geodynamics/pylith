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

#include "TestMeshIOPetsc.hh" // Implementation of class methods

namespace pylith {
    namespace meshio {
        // ----------------------------------------------------------------------------------------
        class TestMeshIOPetsc_GmshTri : public TestMeshIOPetsc {
public:

            void setUp(void) {
                TestMeshIOPetsc::setUp();
                _data = new TestMeshIOPetsc_Data();CPPUNIT_ASSERT(_data);
                _data->numVertices = 9;
                _data->spaceDim = 2;
                _data->numCells = 8;
                _data->cellDim = 2;
                _data->numCorners = 3;

                static const PylithScalar vertices[9*2] = {
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
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[8*3] = {
                    8, 7, 0,
                    7, 8, 4,
                    1, 8, 0,
                    7, 4, 5,
                    3, 8, 6,
                    4, 8, 3,
                    2, 8, 1,
                    8, 2, 6,
                };
                _data->cells = const_cast<PylithInt*>(cells);

                static const PylithInt materialIds[8] = {
                    1, 1, 1, 1, 2, 2, 2, 2,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 6;
                static const PylithInt groupSizes[6] = {3, 3, 3, 3, 3, 1};
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[16] = {
                    0, 5, 7,
                    2, 3, 6,
                    0, 1, 2,
                    3, 4, 5,
                    1, 4, 8,
                    1,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[6] = {
                    "boundary_xneg",
                    "boundary_xpos",
                    "boundary_yneg",
                    "boundary_ypos",
                    "fault",
                    "fault_end",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const PylithInt groupTags[6] = {
                    10, 11, 12, 13, 20, 21,
                };
                _data->groupTags = const_cast<PylithInt*>(groupTags);
                static const char* groupTypes[6] = {
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOPetsc_GmshTri

        // ----------------------------------------------------------------------------------------
        class TestMeshIOPetsc_GmshTri_Ascii : public TestMeshIOPetsc_GmshTri {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOPetsc_GmshTri_Ascii, TestMeshIOPetsc_GmshTri);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOPetsc_GmshTri::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/tri_ascii.msh";
            } // setUp

        }; // class TestMeshIOPetsc_GmshTri_Ascii
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOPetsc_GmshTri_Ascii);

        // ----------------------------------------------------------------------------------------
        class TestMeshIOPetsc_GmshTri_Binary : public TestMeshIOPetsc_GmshTri {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOPetsc_GmshTri_Binary, TestMeshIOPetsc_GmshTri);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOPetsc_GmshTri::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/tri_binary.msh";
            } // setUp

        }; // class TestMeshIOPetsc_GmshTri_Binary
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOPetsc_GmshTri_Binary);

        // ----------------------------------------------------------------------------------------
        class TestMeshIOPetsc_GmshQuad : public TestMeshIOPetsc {
public:

            void setUp(void) {
                TestMeshIOPetsc::setUp();
                _data = new TestMeshIOPetsc_Data();CPPUNIT_ASSERT(_data);
                _data->numVertices = 12;
                _data->spaceDim = 2;
                _data->numCells = 6;
                _data->cellDim = 2;
                _data->numCorners = 4;

                static const PylithScalar vertices[12*2] = {
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
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[6*4] = {
                    0, 1, 10, 9,
                    9, 10, 11, 8,
                    8, 11, 4, 5,
                    1, 2, 6, 10,
                    10, 6, 7, 11,
                    11, 7, 3, 4,
                };
                _data->cells = const_cast<PylithInt*>(cells);

                static const PylithInt materialIds[6] = {
                    1,
                    1,
                    1,
                    2,
                    2,
                    2,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 6;
                static const PylithInt groupSizes[6] = {4, 4, 3, 3, 4, 1};
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[19] = {
                    0, 5, 8, 9,
                    2, 3, 6, 7,
                    0, 1, 2,
                    3, 4, 5,
                    1, 4, 10, 11,
                    1,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[6] = {
                    "boundary_xneg",
                    "boundary_xpos",
                    "boundary_yneg",
                    "boundary_ypos",
                    "fault",
                    "fault_end",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const PylithInt groupTags[6] = {
                    10, 11, 12, 13, 20, 21,
                };
                _data->groupTags = const_cast<PylithInt*>(groupTags);
                static const char* groupTypes[6] = {
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOPetsc_GmshQuad

        // ----------------------------------------------------------------------------------------
        class TestMeshIOPetsc_GmshQuad_Ascii : public TestMeshIOPetsc_GmshQuad {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOPetsc_GmshQuad_Ascii, TestMeshIOPetsc_GmshQuad);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOPetsc_GmshQuad::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/quad_ascii.msh";
            } // setUp

        }; // class TestMeshIOPetsc_GmshQuad_Ascii
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOPetsc_GmshQuad_Ascii);

        // ----------------------------------------------------------------------------------------
        class TestMeshIOPetsc_GmshQuad_Binary : public TestMeshIOPetsc_GmshQuad {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOPetsc_GmshQuad_Binary, TestMeshIOPetsc_GmshQuad);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOPetsc_GmshQuad::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/quad_binary.msh";
            } // setUp

        }; // class TestMeshIOPetsc_GmshQuad_Binary
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOPetsc_GmshQuad_Binary);

        // ----------------------------------------------------------------------------------------
        class TestMeshIOPetsc_GmshTet : public TestMeshIOPetsc {
public:

            void setUp(void) {
                TestMeshIOPetsc::setUp();
                _data = new TestMeshIOPetsc_Data();CPPUNIT_ASSERT(_data);
                _data->numVertices = 30;
                _data->spaceDim = 3;
                _data->numCells = 64;
                _data->cellDim = 3;
                _data->numCorners = 4;

                static const PylithScalar vertices[30*3] = {
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
                    -4.54747351e-13,  +6.73713060e-14,  -4.00000000e+03,
                    +4.00000000e+03,  +0.00000000e+00,  -4.00000000e+03,
                    +2.67937404e+03,  -2.70131505e+03,  -2.14680403e+03,
                    +1.83785671e+03,  -1.83994241e+03,  -1.45460401e+03,
                    +2.00000000e+03,  +2.66666667e+03,  -1.33333333e+03,
                };
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[64*4] = {
                    13, 0, 12, 5,
                    7, 15, 18, 3,
                    19, 18, 15, 25,
                    18, 19, 15, 7,
                    17, 12, 24, 25,
                    17, 18, 25, 24,
                    17, 13, 12, 5,
                    17, 16, 12, 25,
                    18, 15, 24, 3,
                    18, 15, 25, 24,
                    2, 17, 13, 14,
                    17, 13, 14, 24,
                    2, 17, 14, 6,
                    1, 24, 25, 15,
                    16, 1, 12, 25,
                    12, 16, 17, 5,
                    25, 16, 1, 4,
                    17, 14, 18, 24,
                    1, 24, 12, 25,
                    24, 13, 12, 17,
                    17, 14, 6, 18,
                    24, 18, 3, 14,
                    1, 15, 25, 4,
                    19, 25, 15, 4,
                    9, 20, 16, 27,
                    5, 9, 16, 27,
                    6, 18, 11, 29,
                    11, 18, 23, 29,
                    26, 20, 22, 27,
                    9, 22, 20, 27,
                    22, 6, 11, 29,
                    17, 6, 22, 29,
                    10, 26, 19, 21,
                    19, 26, 8, 21,
                    22, 26, 27, 28,
                    11, 23, 22, 29,
                    17, 26, 22, 28,
                    26, 16, 20, 27,
                    26, 23, 19, 25,
                    16, 17, 26, 25,
                    16, 8, 4, 25,
                    26, 10, 19, 23,
                    18, 19, 23, 25,
                    10, 19, 23, 7,
                    16, 26, 8, 25,
                    19, 8, 25, 4,
                    26, 16, 8, 20,
                    19, 26, 25, 8,
                    18, 19, 7, 23,
                    9, 22, 27, 28,
                    6, 17, 18, 29,
                    5, 9, 27, 28,
                    5, 17, 9, 28,
                    9, 17, 22, 28,
                    16, 27, 28, 5,
                    28, 27, 16, 26,
                    28, 17, 16, 5,
                    16, 17, 28, 26,
                    29, 18, 25, 17,
                    29, 25, 18, 23,
                    26, 29, 25, 17,
                    26, 25, 29, 23,
                    22, 29, 26, 17,
                    22, 26, 29, 23,
                };
                _data->cells = const_cast<PylithInt*>(cells);

                static const PylithInt materialIds[64] = {
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 8;
                static const PylithInt groupSizes[8] = {9, 9, 9, 9, 9, 9, 9, 3};
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[66] = {
                    0, 1, 2, 3, 12, 13, 14, 15, 24,
                    8, 9, 10, 11, 20, 21, 22, 23, 26,
                    0, 1, 4, 5, 8, 9, 12, 16, 20,
                    2, 3, 6, 7, 10, 11, 14, 18, 23,
                    1, 3, 4, 7, 8, 10, 15, 19, 21,
                    0, 2, 5, 6, 9, 11, 13, 17, 22,
                    4, 5, 6, 7, 16, 17, 18, 19, 25,
                    4, 5, 16,
                };

                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[8] = {
                    "boundary_xneg",
                    "boundary_xpos",
                    "boundary_yneg",
                    "boundary_ypos",
                    "boundary_zneg",
                    "boundary_zpos",
                    "fault",
                    "fault_end",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const PylithInt groupTags[8] = {
                    10, 11, 12, 13, 14, 15, 20, 21,
                };
                _data->groupTags = const_cast<PylithInt*>(groupTags);
                static const char* groupTypes[8] = {
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOPetsc_GmshTet

        // ----------------------------------------------------------------------------------------
        class TestMeshIOPetsc_GmshTet_Ascii : public TestMeshIOPetsc_GmshTet {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOPetsc_GmshTet_Ascii, TestMeshIOPetsc_GmshTet);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOPetsc_GmshTet::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/tet_ascii.msh";
            } // setUp

        }; // class TestMeshIOPetsc_GmshTet_Ascii
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOPetsc_GmshTet_Ascii);

        // ----------------------------------------------------------------------------------------
        class TestMeshIOPetsc_GmshTet_Binary : public TestMeshIOPetsc_GmshTet {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOPetsc_GmshTet_Binary, TestMeshIOPetsc_GmshTet);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOPetsc_GmshTet::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/tet_binary.msh";
            } // setUp

        }; // class TestMeshIOPetsc_GmshTet_Binary
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOPetsc_GmshTet_Binary);

        // --------------------------------------------------------------
        class TestMeshIOPetsc_GmshHex : public TestMeshIOPetsc {
public:

            void setUp(void) {
                TestMeshIOPetsc::setUp();
                _data = new TestMeshIOPetsc_Data();CPPUNIT_ASSERT(_data);
                _data->numVertices = 48;
                _data->spaceDim = 3;
                _data->numCells = 18;
                _data->cellDim = 3;
                _data->numCorners = 8;

                static const PylithScalar vertices[48*3] = {
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
                _data->vertices = const_cast<PylithScalar*>(vertices);

                static const PylithInt cells[18*8] = {
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
                _data->cells = const_cast<PylithInt*>(cells);

                static const PylithInt materialIds[18] = {
                    1, 1, 1, 1, 1, 1, 1, 1, 1,
                    2, 2, 2, 2, 2, 2, 2, 2, 2,
                };
                _data->materialIds = const_cast<PylithInt*>(materialIds);

                _data->numGroups = 8;
                static const PylithInt groupSizes[8] = {16, 16, 12, 12, 12, 12, 16, 4};
                _data->groupSizes = const_cast<PylithInt*>(groupSizes);
                static const PylithInt groups[100] = {
                    0, 1, 2, 3, 12, 13, 14, 15, 16, 17, 18, 19, 36, 37, 38, 39,
                    8, 9, 10, 11, 28, 29, 30, 31, 32, 33, 34, 35, 44, 45, 46, 47,
                    0, 1, 4, 5, 8, 9, 12, 13, 20, 21, 28, 29,
                    2, 3, 6, 7, 10, 11, 16, 17, 24, 25, 34, 35,
                    1, 3, 4, 7, 8, 10, 18, 19, 26, 27, 30, 31,
                    0, 2, 5, 6, 9, 11, 14, 15, 22, 23, 32, 33,
                    4, 5, 6, 7, 20, 21, 22, 23, 24, 25, 26, 27, 40, 41, 42, 43,
                    4, 5, 20, 21,
                };
                _data->groups = const_cast<PylithInt*>(groups);
                static const char* groupNames[8] = {
                    "boundary_xneg",
                    "boundary_xpos",
                    "boundary_yneg",
                    "boundary_ypos",
                    "boundary_zneg",
                    "boundary_zpos",
                    "fault",
                    "fault_end",
                };
                _data->groupNames = const_cast<char**>(groupNames);
                static const PylithInt groupTags[8] = {
                    10, 11, 12, 13, 14, 15, 20, 21,
                };
                _data->groupTags = const_cast<PylithInt*>(groupTags);
                static const char* groupTypes[8] = {
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // class TestMeshIOPetsc_GmshHex

        // ----------------------------------------------------------------------------------------
        class TestMeshIOPetsc_GmshHex_Ascii : public TestMeshIOPetsc_GmshHex {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOPetsc_GmshHex_Ascii, TestMeshIOPetsc_GmshHex);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOPetsc_GmshHex::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/hex_ascii.msh";
            } // setUp

        }; // class TestMeshIOPetsc_GmshHex_Ascii
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOPetsc_GmshHex_Ascii);

        // ----------------------------------------------------------------------------------------
        class TestMeshIOPetsc_GmshHex_Binary : public TestMeshIOPetsc_GmshHex {
            CPPUNIT_TEST_SUB_SUITE(TestMeshIOPetsc_GmshHex_Binary, TestMeshIOPetsc_GmshHex);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestMeshIOPetsc_GmshHex::setUp();
                CPPUNIT_ASSERT(_data);
                _data->filename = "data/hex_binary.msh";
            } // setUp

        }; // class TestMeshIOPetsc_GmshHex_Binary
        CPPUNIT_TEST_SUITE_REGISTRATION(TestMeshIOPetsc_GmshHex_Binary);

    } // meshio
} // pylith

// End of file
