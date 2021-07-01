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
        class TestAdjustTopology_QuadA : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_QuadA, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/quad_a.mesh";

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

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 4 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
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

        }; // TestAdjustTopology_QuadA
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_QuadA);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_QuadB : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_QuadB, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/quad_b.mesh";

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

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 4 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
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

        }; // TestAdjustTopology_QuadB
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_QuadB);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_QuadC : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_QuadC, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/quad_c.mesh";

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

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 4 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
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

        }; // TestAdjustTopology_QuadC
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_QuadC);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_QuadD : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_QuadD, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/quad_d.mesh";

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

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 4 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
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

        }; // TestAdjustTopology_QuadD
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_QuadD);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_QuadE : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_QuadE, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/quad_e.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 12;

                static const size_t numCells = 6;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 4, 4, 4, 4 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 4+2, 6+4 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_QuadE
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_QuadE);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_QuadF : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_QuadF, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/quad_f.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 12;

                static const size_t numCells = 6;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 4, 4, 4, 4 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 4+2, 6+4 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_QuadF
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_QuadF);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_QuadG : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_QuadG, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/quad_g.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 14;

                static const size_t numCells = 7;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 4, 4, 4, 4, 4, 4, 4 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 0, 2, 2, 100, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 3+2, 6+4 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_QuadG
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_QuadG);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_QuadH : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_QuadH, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/quad_h.mesh";

                _data->numFaults = 2;
                static const char* const faultSurfaceLabels[2] = { "faultA", "faultB" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[2] = { NULL, "faultB-edge" };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[2] = { 101, 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 22;

                static const size_t numCells = 14;
                _data->numCells = numCells;

                static const int numCorners[numCells] = {
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = {
                    10, 10, 11, 10, 10, 11, 12, 12, 11,
                    100, 100, 101, 101, 101 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 3;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 6+4, 8+6, 1 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "faultB", "faultA", "faultB-edge" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_QuadH
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_QuadH);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_QuadI : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_QuadI, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/quad_i.mesh";

                static const size_t numFaults = 1;
                _data->numFaults = numFaults;
                static const char* const faultSurfaceLabels[numFaults] = { "fault", };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[numFaults] = { "edge" };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[numFaults] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 2;
                _data->spaceDim = 2;
                _data->numVertices = 14;

                static const size_t numCells = 8;
                _data->numCells = numCells;

                static const int numCorners[numCells] = {
                    4, 4, 4, 4, 4, 4,
                    4, 4, };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = {
                    10, 10, 10, 11, 10, 11,
                    100, 100, };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 4;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 3+2, 4+2, 1, 5+4 }; // vertices + edges
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output2", "output1", "edge", "fault"};
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_QuadI
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_QuadI);

    } // faults
} // pylith

// End of file
