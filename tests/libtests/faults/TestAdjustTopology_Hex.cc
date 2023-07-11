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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestAdjustTopology.hh" // USES TestAdjustTopology_Data

namespace pylith {
    namespace faults {
        class TestAdjustTopology_Hex;
    }
}

class pylith::faults::TestAdjustTopology_Hex {
public:

    // Data factory methods
    static TestAdjustTopology_Data* caseA(void);

    static TestAdjustTopology_Data* caseB(void);

    static TestAdjustTopology_Data* caseC(void);

    static TestAdjustTopology_Data* caseD(void);

    static TestAdjustTopology_Data* caseE(void);

    static TestAdjustTopology_Data* caseF(void);

    static TestAdjustTopology_Data* caseG(void);

    static TestAdjustTopology_Data* caseH(void);

    static TestAdjustTopology_Data* caseI(void);

    static TestAdjustTopology_Data* caseJ(void);

private:

    TestAdjustTopology_Hex(void); ///< Not implemented
}; // TestAdjustTopology_Hex

// ------------------------------------------------------------------------------------------------
#include "catch2/catch_test_macros.hpp"

TEST_CASE("TestAdjustTopology_HexA", "[TestAdjustTopology][Hex]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Hex::caseA()).run();
}
TEST_CASE("TestAdjustTopology_HexB", "[TestAdjustTopology][Hex]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Hex::caseB()).run();
}
TEST_CASE("TestAdjustTopology_HexC", "[TestAdjustTopology][Hex]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Hex::caseC()).run();
}
TEST_CASE("TestAdjustTopology_HexD", "[TestAdjustTopology][Hex]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Hex::caseD()).run();
}
TEST_CASE("TestAdjustTopology_HexE", "[TestAdjustTopology][Hex]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Hex::caseE()).run();
}
TEST_CASE("TestAdjustTopology_HexF", "[TestAdjustTopology][Hex]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Hex::caseF()).run();
}
TEST_CASE("TestAdjustTopology_HexG", "[TestAdjustTopology][Hex]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Hex::caseG()).run();
}
TEST_CASE("TestAdjustTopology_HexH", "[TestAdjustTopology][Hex]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Hex::caseH()).run();
}
TEST_CASE("TestAdjustTopology_HexI", "[TestAdjustTopology][Hex]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Hex::caseI()).run();
}
TEST_CASE("TestAdjustTopology_HexJ", "[TestAdjustTopology][Hex]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Hex::caseJ()).run();
}

// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Hex::caseA(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/hex_a.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 16;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 6, 6, 6 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Hex::caseB(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/hex_b.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 16;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 6, 6, 6 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseB


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Hex::caseC(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/hex_c.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 16;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 6, 6, 6 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseC


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Hex::caseD(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/hex_d.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 16;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 6, 6, 6 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseD


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Hex::caseE(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/hex_e.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 16;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 6, 6, 6 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseE


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Hex::caseF(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/hex_e.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 16;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 6, 6, 6 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseF


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Hex::caseG(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/hex_g.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 24;

    static const size_t numCells = 6;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 6, 6, 6, 6, 6, 6 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 12+14+4 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseG


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Hex::caseH(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/hex_h.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 24;

    static const size_t numCells = 6;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 6, 6, 6, 6, 6, 6 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 12+14+4 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseH


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Hex::caseI(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/hex_i.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 36;

    static const size_t numCells = 9;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 6, 6, 6, 6, 6, 6, 6, 6, 6, };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 0, 0, 0, 0, 2, 100, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 10+13+4, 12+14+4 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseI


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Hex::caseJ(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/hex_j.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { "fault_edge" };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 61;

    static const size_t numCells = 22;
    data->numCells = numCells;

    static const int numCorners[numCells] = {
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6,
    };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = {
        10, 10, 11, 11, 10, 10, 11, 11, 11, 11,
        10, 10, 10, 10, 11, 11, 11, 11, 10, 10,
        100, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 3;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 21+31+10, 5+4, 6+7+2 + 1+2+3 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault_edge", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseJ


// End of file
