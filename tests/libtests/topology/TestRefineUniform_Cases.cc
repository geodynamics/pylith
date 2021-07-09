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

#include "TestRefineUniform.hh" // Implementation of class methods

// -----------------------------------------------------------------------------
namespace pylith {
    namespace topology {

        // ---------------------------------------------------------------------
        class TestRefineUniform_Tri_2xNofault : public TestRefineUniform {

            CPPUNIT_TEST_SUB_SUITE( TestRefineUniform_Tri_2xNofault, TestRefineUniform );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestRefineUniform::setUp();

                _data->filename = "data/fourtri3.mesh";
                _data->refineLevel = 1;
                _data->faultA = NULL;
                _data->faultB = NULL;
                _data->isSimplexMesh = true;

                _data->numVertices = 13;
                _data->spaceDim = 2;
                _data->numCells = 16;
                _data->numCellsCohesive = 0;
                _data->cellDim = 2;
                _data->numCorners = 3;
                _data->numCornersCohesive = 4;

                _data->matIdSum = 8*1 + 8*2;

                _data->numGroups = 4;

                static const int _groupSizes[4] = {
                    9, 5, 2, 9,
                };
                _data->groupSizes = const_cast<int*>(_groupSizes);

                static const char* _groupNames[4] = {
                    "edge 1",
                    "edge 2",
                    "end points",
                    "fault",
                };
                _data->groupNames = const_cast<const char**>(_groupNames);

                static const char* _groupTypes[4] = {
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<const char**>(_groupTypes);

            }   // setUp

        };  // _TestRefineUniform_Tri_2xNofault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestRefineUniform_Tri_2xNofault );


        // ---------------------------------------------------------------------
        class TestRefineUniform_Tri_2xFault : public TestRefineUniform {

            CPPUNIT_TEST_SUB_SUITE( TestRefineUniform_Tri_2xFault, TestRefineUniform );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestRefineUniform::setUp();

                _data->filename = "data/fourtri3.mesh";
                _data->refineLevel = 1;
                _data->faultA = "fault";
                _data->faultB = NULL;
                _data->isSimplexMesh = true;

                _data->numVertices = 18;
                _data->spaceDim = 2;
                _data->numCells = 16;
                _data->numCellsCohesive = 4;
                _data->cellDim = 2;
                _data->numCorners = 3;
                _data->numCornersCohesive = 4;

                _data->matIdSum = 8*1 + 8*2 + 4*100;

                _data->numGroups = 4;

                static const int _groupSizes[4] = {
                    11, 6, 2, 18,
                };
                _data->groupSizes = const_cast<int*>(_groupSizes);

                static const char* _groupNames[4] = {
                    "edge 1",
                    "edge 2",
                    "end points",
                    "fault",
                };
                _data->groupNames = const_cast<const char**>(_groupNames);

                static const char* _groupTypes[4] = {
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<const char**>(_groupTypes);

            }   // setUp

        };  // _TestRefineUniform_Tri_2xFault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestRefineUniform_Tri_2xFault );


        // ---------------------------------------------------------------------
        class TestRefineUniform_Quad_2xNofault : public TestRefineUniform {

            CPPUNIT_TEST_SUB_SUITE( TestRefineUniform_Quad_2xNofault, TestRefineUniform );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestRefineUniform::setUp();

                _data->filename = "data/fourquad4.mesh";
                _data->refineLevel = 1;
                _data->faultA = NULL;
                _data->faultB = NULL;
                _data->isSimplexMesh = false;

                _data->numVertices = 25;
                _data->spaceDim = 2;
                _data->numCells = 16;
                _data->numCellsCohesive = 0;
                _data->cellDim = 2;
                _data->numCorners = 4;
                _data->numCornersCohesive = 4;

                _data->matIdSum = 8*1 + 8*2;

                _data->numGroups = 3;

                static const int _groupSizes[3] = {
                    9, 9, 9,
                };
                _data->groupSizes = const_cast<int*>(_groupSizes);

                static const char* _groupNames[3] = {
                    "edge 1",
                    "end points",
                    "fault",
                };
                _data->groupNames = const_cast<const char**>(_groupNames);

                static const char* _groupTypes[3] = {
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<const char**>(_groupTypes);

            }   // setUp

        };  // _TestRefineUniform_Quad_2xNofault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestRefineUniform_Quad_2xNofault );


        // ---------------------------------------------------------------------
        class TestRefineUniform_Quad_2xFault : public TestRefineUniform {

            CPPUNIT_TEST_SUB_SUITE( TestRefineUniform_Quad_2xFault, TestRefineUniform );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestRefineUniform::setUp();

                _data->filename = "data/fourquad4.mesh";
                _data->refineLevel = 1;
                _data->faultA = "fault";
                _data->faultB = NULL;
                _data->isSimplexMesh = false;

                _data->numVertices = 30;
                _data->spaceDim = 2;
                _data->numCells = 16;
                _data->numCellsCohesive = 4;
                _data->cellDim = 2;
                _data->numCorners = 4;
                _data->numCornersCohesive = 4;

                _data->matIdSum = 8*1 + 8*2 + 4*100;

                _data->numGroups = 3;

                static const int _groupSizes[3] = {
                    6+4, // vertices+edges
                    5+4,
                    10+8,
                };
                _data->groupSizes = const_cast<int*>(_groupSizes);

                static const char* _groupNames[3] = {
                    "edge 1",
                    "end points",
                    "fault",
                };
                _data->groupNames = const_cast<const char**>(_groupNames);

                static const char* _groupTypes[3] = {
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<const char**>(_groupTypes);

            }   // setUp

        };  // _TestRefineUniform_Quad_2xFault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestRefineUniform_Quad_2xFault );


        // ---------------------------------------------------------------------
        class TestRefineUniform_Tet_2xNofault : public TestRefineUniform {

            CPPUNIT_TEST_SUB_SUITE( TestRefineUniform_Tet_2xNofault, TestRefineUniform );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestRefineUniform::setUp();

                _data->filename = "data/twotet4.mesh";
                _data->refineLevel = 1;
                _data->faultA = NULL;
                _data->faultB = NULL;
                _data->isSimplexMesh = true;

                _data->numVertices = 14;
                _data->spaceDim = 3;
                _data->numCells = 16;
                _data->numCellsCohesive = 0;
                _data->cellDim = 3;
                _data->numCorners = 4;
                _data->numCornersCohesive = 6;

                _data->matIdSum = 8*1 + 8*2;

                _data->numGroups = 4;

                static const int _groupSizes[4] = {
                    3+2+0, // vertices+edges+faces
                    3+2+0,
                    2+0+0,
                    6+9+4
                };
                _data->groupSizes = const_cast<int*>(_groupSizes);

                static const char* _groupNames[4] = {
                    "edge 1",
                    "edge 2",
                    "end points",
                    "fault",
                };
                _data->groupNames = const_cast<const char**>(_groupNames);

                static const char* _groupTypes[4] = {
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<const char**>(_groupTypes);

            }   // setUp

        };  // _TestRefineUniform_Tet_2xNofault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestRefineUniform_Tet_2xNofault );


        // ---------------------------------------------------------------------
        class TestRefineUniform_Tet_2xFault : public TestRefineUniform {

            CPPUNIT_TEST_SUB_SUITE( TestRefineUniform_Tet_2xFault, TestRefineUniform );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestRefineUniform::setUp();

                _data->filename = "data/twotet4.mesh";
                _data->refineLevel = 1;
                _data->faultA = "fault";
                _data->faultB = NULL;
                _data->isSimplexMesh = true;

                _data->numVertices = 20;
                _data->spaceDim = 3;
                _data->numCells = 16;
                _data->numCellsCohesive = 4;
                _data->cellDim = 3;
                _data->numCorners = 4;
                _data->numCornersCohesive = 6;

                _data->matIdSum = 8*1 + 8*2 + 4*100;

                _data->numGroups = 4;

                static const int _groupSizes[4] = {
                    6,
                    6,
                    2,
                    38,
                };
                _data->groupSizes = const_cast<int*>(_groupSizes);

                static const char* _groupNames[4] = {
                    "edge 1",
                    "edge 2",
                    "end points",
                    "fault",
                };
                _data->groupNames = const_cast<const char**>(_groupNames);

                static const char* _groupTypes[4] = {
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<const char**>(_groupTypes);

            }   // setUp

        };  // _TestRefineUniform_Tet_2xFault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestRefineUniform_Tet_2xFault );


        // ---------------------------------------------------------------------
        class TestRefineUniform_Hex_2xNofault : public TestRefineUniform {

            CPPUNIT_TEST_SUB_SUITE( TestRefineUniform_Hex_2xNofault, TestRefineUniform );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestRefineUniform::setUp();

                _data->filename = "data/twohex8.mesh";
                _data->refineLevel = 1;
                _data->faultA = NULL;
                _data->faultB = NULL;
                _data->isSimplexMesh = true;

                _data->numVertices = 45;
                _data->spaceDim = 3;
                _data->numCells = 16;
                _data->numCellsCohesive = 0;
                _data->cellDim = 3;
                _data->numCorners = 8;
                _data->numCornersCohesive = 8;

                _data->matIdSum = 8*1 + 8*2;

                _data->numGroups = 4;

                static const int _groupSizes[4] = {
                    2+0+0, // vertices+edges+faces
                    9+12+4,
                    9+12+4,
                    9+12+4,
                };
                _data->groupSizes = const_cast<int*>(_groupSizes);

                static const char* _groupNames[4] = {
                    "end points",
                    "face 1",
                    "face 2",
                    "fault",
                };
                _data->groupNames = const_cast<const char**>(_groupNames);

                static const char* _groupTypes[4] = {
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<const char**>(_groupTypes);

            }   // setUp

        };  // _TestRefineUniform_Hex_2xNofault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestRefineUniform_Hex_2xNofault );


        // ---------------------------------------------------------------------
        class TestRefineUniform_Hex_2xFault : public TestRefineUniform {

            CPPUNIT_TEST_SUB_SUITE( TestRefineUniform_Hex_2xFault, TestRefineUniform );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestRefineUniform::setUp();

                _data->filename = "data/twohex8.mesh";
                _data->refineLevel = 1;
                _data->faultA = "fault";
                _data->faultB = NULL;
                _data->isSimplexMesh = true;

                _data->numVertices = 54;
                _data->spaceDim = 3;
                _data->numCells = 16;
                _data->numCellsCohesive = 4;
                _data->cellDim = 3;
                _data->numCorners = 8;
                _data->numCornersCohesive = 8;

                _data->matIdSum = 8*1 + 8*2 + 4*100;

                _data->numGroups = 4;

                static const int _groupSizes[4] = {
                    2+0+0, // vertices+edges+faces
                    9+12+4,
                    12+14+4,
                    18+24+8,
                };
                _data->groupSizes = const_cast<int*>(_groupSizes);

                static const char* _groupNames[4] = {
                    "end points",
                    "face 1",
                    "face 2",
                    "fault",
                };
                _data->groupNames = const_cast<const char**>(_groupNames);

                static const char* _groupTypes[4] = {
                    "vertex",
                    "vertex",
                    "vertex",
                    "vertex",
                };
                _data->groupTypes = const_cast<const char**>(_groupTypes);

            }   // setUp

        };  // _TestRefineUniform_Hex_2xFault
        CPPUNIT_TEST_SUITE_REGISTRATION( TestRefineUniform_Hex_2xFault );


    }   // topology
}   // pylith


// End of file
