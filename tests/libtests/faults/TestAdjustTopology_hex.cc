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
        class TestAdjustTopology_HexA : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_HexA, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/hex_a.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 16;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 6, 6, 6 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_HexA
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_HexA);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_HexB : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_HexB, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/hex_b.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 16;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 6, 6, 6 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_HexB
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_HexB);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_HexC : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_HexC, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/hex_c.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 16;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 6, 6, 6 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_HexC
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_HexC);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_HexD : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_HexD, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/hex_d.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 16;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 6, 6, 6 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_HexD
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_HexD);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_HexE : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_HexE, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/hex_e.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 16;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 6, 6, 6 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_HexE
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_HexE);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_HexF : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_HexF, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/hex_e.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 16;

                static const size_t numCells = 3;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 6, 6, 6 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_HexF
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_HexF);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_HexG : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_HexG, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/hex_g.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 24;

                static const size_t numCells = 6;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 6, 6, 6, 6, 6, 6 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 8+8+2, 12+14+4 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_HexG
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_HexG);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_HexH : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_HexH, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/hex_h.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 24;

                static const size_t numCells = 6;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 6, 6, 6, 6, 6, 6 };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 8+8+2, 12+14+4 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_HexH
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_HexH);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_HexI : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_HexI, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/hex_i.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { NULL };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 36;

                static const size_t numCells = 9;
                _data->numCells = numCells;

                static const int numCorners[numCells] = { 6, 6, 6, 6, 6, 6, 6, 6, 6, };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = { 0, 0, 0, 0, 0, 0, 2, 100, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 2;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 10+13+4, 12+14+4 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_HexI
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_HexI);

        // ----------------------------------------------------------------------------------------
        class TestAdjustTopology_HexJ : public TestAdjustTopology {
            CPPUNIT_TEST_SUB_SUITE(TestAdjustTopology_HexJ, TestAdjustTopology);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestAdjustTopology::setUp();

                _data->filename = "data/hex_j.mesh";

                _data->numFaults = 1;
                static const char* const faultSurfaceLabels[1] = { "fault" };
                _data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
                static const char* const faultEdgeLabels[1] = { "fault_edge" };
                _data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
                static const int interfaceIds[1] = { 100 };
                _data->interfaceIds = const_cast<const int*>(interfaceIds);

                _data->cellDim = 3;
                _data->spaceDim = 3;
                _data->numVertices = 61;

                static const size_t numCells = 22;
                _data->numCells = numCells;

                static const int numCorners[numCells] = {
                    6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                    6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                    6, 6,
                };
                _data->numCorners = const_cast<int*>(numCorners);
                static const int materialIds[numCells] = {
                    10, 10, 11, 11, 10, 10, 11, 11, 11, 11,
                    10, 10, 10, 10, 11, 11, 11, 11, 10, 10,
                    100, 100 };
                _data->materialIds = const_cast<int*>(materialIds);

                static const size_t numGroups = 3;
                _data->numGroups = numGroups;
                static const int groupSizes[numGroups] = { 21+31+10, 5+4, 6+7+2 + 1+2+3 }; // vertices + edges + faces
                _data->groupSizes = const_cast<int*>(groupSizes);
                static const char* groupNames[numGroups] = { "output", "fault_edge", "fault" };
                _data->groupNames = const_cast<char**>(groupNames);
                static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex" };
                _data->groupTypes = const_cast<char**>(groupTypes);
            } // setUp

        }; // TestAdjustTopology_HexJ
        CPPUNIT_TEST_SUITE_REGISTRATION(TestAdjustTopology_HexJ);

    } // faults
} // pylith

#if 0

// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8.hh" // USES CohesiveDataHex8

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8(void) { // testAdjustTopologyHex8
    PYLITH_METHOD_BEGIN;

    CohesiveDataHex8 data;
    FaultCohesiveTract fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyHex8


// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8b.hh" // USES CohesiveDataHex8b

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8b(void) { // testAdjustTopologyHex8b
    PYLITH_METHOD_BEGIN;

    CohesiveDataHex8b data;
    FaultCohesiveTract fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyHex8b


// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8c.hh" // USES CohesiveDataHex8c

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8c(void) { // testAdjustTopologyHex8c
    PYLITH_METHOD_BEGIN;

    CohesiveDataHex8c data;
    FaultCohesiveTract fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyHex8c


// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8d.hh" // USES CohesiveDataHex8d

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8d(void) { // testAdjustTopologyHex8d
    PYLITH_METHOD_BEGIN;

    CohesiveDataHex8d data;
    FaultCohesiveTract fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyHex8d


// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8e.hh" // USES CohesiveDataHex8e

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8e(void) { // testAdjustTopologyHex8e
    PYLITH_METHOD_BEGIN;

    CohesiveDataHex8e data;
    FaultCohesiveTract fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyHex8e


// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8f.hh" // USES CohesiveDataHex8f

// Test adjustTopology() with 3-D hexahedral element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8f(void) { // testAdjustTopologyHex8f
    PYLITH_METHOD_BEGIN;

    CohesiveDataHex8f data;
    FaultCohesiveTract fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyHex8f


// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8g.hh" // USES CohesiveDataHex8g

// Test adjustTopology() with 3-D hexahedral element (2 cells easy).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8g(void) { // testAdjustTopologyHex8g
    PYLITH_METHOD_BEGIN;

    CohesiveDataHex8g data;
    FaultCohesiveTract fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyHex8g


// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8h.hh" // USES CohesiveDataHex8h

// Test adjustTopology() with 3-D hexahedral element (2 cells difficult).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8h(void) { // testAdjustTopologyHex8h
    PYLITH_METHOD_BEGIN;

    CohesiveDataHex8h data;
    FaultCohesiveTract fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyHex8h


// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8i.hh" // USES CohesiveDataHex8i

// Test adjustTopology() with 3-D hexahedral element (vertex/edge on fault).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8i(void) { // testAdjustTopologyHex8i
    PYLITH_METHOD_BEGIN;

    CohesiveDataHex8i data;
    FaultCohesiveTract fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyHex8i


// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8j.hh" // USES CohesiveDataHex8j

// Test adjustTopology() with 3-D hexahedral element (embedded fault).
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8j(void) { // testAdjustTopologyHex8j
    PYLITH_METHOD_BEGIN;

    CohesiveDataHex8j data;
    FaultCohesiveTract fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyHex8j


// ----------------------------------------------------------------------
#include "data/CohesiveDataTri3Lagrange.hh" // USES CohesiveDataTri3Lagrange

// Test adjustTopology() with 2-D triangular element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTri3Lagrange(void) { // testAdjustTopologyTri3Lagrange
    PYLITH_METHOD_BEGIN;

    CohesiveDataTri3Lagrange data;
    FaultCohesiveKin fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyTri3Lagrange


// ----------------------------------------------------------------------
#include "data/CohesiveDataQuad4Lagrange.hh" // USES CohesiveDataQuad4Lagrange

// Test adjustTopology() with 2-D quadrilateral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyQuad4Lagrange(void) { // testAdjustTopologyQuad4Lagrange
    PYLITH_METHOD_BEGIN;

    CohesiveDataQuad4Lagrange data;
    FaultCohesiveKin fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyQuad4Lagrange


// ----------------------------------------------------------------------
#include "data/CohesiveDataTet4Lagrange.hh" // USES CohesiveDataTet4Lagrange

// Test adjustTopology() with 3-D tetrahedral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyTet4Lagrange(void) { // testAdjustTopologyTet4Lagrange
    PYLITH_METHOD_BEGIN;

    CohesiveDataTet4Lagrange data;
    FaultCohesiveKin fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyTet4Lagrange


// ----------------------------------------------------------------------
#include "data/CohesiveDataHex8Lagrange.hh" // USES CohesiveDataHex8Lagrange

// Test adjustTopology() with 3-D hexahedral element for Lagrange
// multipliers.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyHex8Lagrange(void) { // testAdjustTopologyHex8Lagrange
    PYLITH_METHOD_BEGIN;

    CohesiveDataHex8Lagrange data;
    FaultCohesiveKin fault;
    _testAdjustTopology(&fault, data);

    PYLITH_METHOD_END;
} // testAdjustTopologyHex8Lagrange


#endif

// End of file
