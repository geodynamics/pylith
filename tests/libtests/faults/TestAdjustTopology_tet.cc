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

#include "TestAdjustTopology.hh" // Implementation of class methods

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace faults {
        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TetA : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TetA, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tet_a.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 8;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 5 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4+1, 6+6+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TetA
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TetA);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TetB : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TetB, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tet_b.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 8;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 5 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4+1, 6+6+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TetB
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TetB);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TetC : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TetC, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tet_c.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 8;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 5 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4+1, 6+6+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TetC
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TetC);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TetD : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TetD, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tet_d.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 8;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 5 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4+1, 6+6+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TetD
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TetD);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TetF : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TetF, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tet_f.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 8;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 5 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4+1, 6+6+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TetF
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TetF);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TetG : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TetG, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tet_g.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 8;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 5 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4+1, 6+6+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TetG
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TetG);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TetH : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TetH, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tet_h.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 8;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 5 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4+1, 6+6+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TetH
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TetH);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TetI : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TetI, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tet_i.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 10;

                static const size_t numCells = 6;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 4, 4, 5, 5 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 2+0+0, 8+10+4 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TetI
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TetI);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TetJ : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TetJ, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tet_j.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 12;

                static const size_t numCells = 7;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 4, 4, 4, 4, 5 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 2, 2, 0, 2, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 8+10+4, 6+6+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TetJ
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TetJ);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TetK : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TetK, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tet_k.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { "fault_edge" };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 354;

                static const size_t numCells = 1427;
                _data->numCells = numCells;

                static const int numCorners[numCells] = {
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 1400
                    4, 4, 4, 4, 4,
                    5, 5, 5, 5, 5,
                    5, 5, 5, 5, 5,
                    5, 5, 5, 5, 5,
                    5, 5, 5, 5, 5,
                    5, 5,
                };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = {
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   10, 10, 10, 10, 10, 10, 10, 10, 10, 10, // 1410
                    10, 10, 10, 10, 10,
                    100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100,
                    100, 100,
                };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 3;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 62+152+90, 10+9, 27+94+22 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "face_zpos", "fault_edge", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TetK
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TetK);

    } // faults
} // pylith

// End of file
