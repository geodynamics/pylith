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
        class TestAdjustTopology_TriA : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TriA, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tri_a.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 6;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 3, 3, 4 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4, 4+2 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TriA
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TriA);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TriB : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TriB, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tri_b.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 6;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 3, 3, 4 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4, 4+2 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TriB
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TriB);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TriC : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TriC, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tri_c.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 6;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 3, 3, 4 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4, 4+2 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TriC
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TriC);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TriD : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TriD, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tri_d.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 9;

                static const size_t numCells = 6;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 3, 3, 3, 3, 4, 4, };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4, 6+4 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TriD
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TriD);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TriE : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TriE, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tri_e.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 9;

                static const size_t numCells = 6;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 3, 3, 3, 3, 4, 4, };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+4, 6+4 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TriE
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TriE);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TriF : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TriF, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tri_f.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 8;

                static const size_t numCells = 5;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 3, 3, 3, 3, 4, };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 2, 2, 100, };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 4+2, 4+2 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TriF
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TriF);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TriG : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TriG, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tri_g.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { "edge" };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 8;

                static const size_t numCells = 7;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 3, 3, 3, 3, 3, 3, 4, };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 2, 0, 2, 0, 2, 100, };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 3;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 5+3, 1, 3+2 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "edge", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TriG
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TriG);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TriH : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TriH, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tri_h.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { "edge" };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 11;

                static const size_t numCells = 13;
                _data->numCells = numCells;

                static const int numCorners[numCells] = {
                    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                    4, 4, };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = {
                    0, 0, 0, 0,
                    2, 2, 3, 3, 3, 2, 3,
                    100, 100,
                };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 3;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 3+2, 2, 4+4, }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "edge", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_TriH
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TriH);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_TriI : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_TriI, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/tri_i.mesh";
                _data->failureExpected = true;

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { "edge" };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

            } // setUp

        }; // TestAdjustTopology_TriI
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_TriI);

    } // faults
} // pylith

// End of file
