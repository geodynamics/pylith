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
        class TestAdjustTopology_Quad;
    }
}

class pylith::faults::TestAdjustTopology_Quad {
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

private:

    TestAdjustTopology_Quad(void); ///< Not implemented
}; // TestAdjustTopology_Quad

// ------------------------------------------------------------------------------------------------
#include "catch2/catch_test_macros.hpp"

TEST_CASE("TestAdjustTopology_QuadA", "[TestAdjustTopology][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseA()).run();
}
TEST_CASE("TestAdjustTopology_QuadB", "[TestAdjustTopology][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseB()).run();
}
TEST_CASE("TestAdjustTopology_QuadC", "[TestAdjustTopology][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseC()).run();
}
TEST_CASE("TestAdjustTopology_QuadD", "[TestAdjustTopology][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseD()).run();
}
TEST_CASE("TestAdjustTopology_QuadE", "[TestAdjustTopology][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseE()).run();
}
TEST_CASE("TestAdjustTopology_QuadF", "[TestAdjustTopology][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseF()).run();
}
TEST_CASE("TestAdjustTopology_QuadG", "[TestAdjustTopology][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseG()).run();
}
TEST_CASE("TestAdjustTopology_QuadH", "[TestAdjustTopology][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseH()).run();
}
TEST_CASE("TestAdjustTopology_QuadI", "[TestAdjustTopology][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseI()).run();
}

// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseA(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/quad_a.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 8;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 4 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 4+2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseB(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/quad_b.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 8;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 4 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 4+2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseC(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/quad_c.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 8;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 4 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 4+2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseD(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/quad_d.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 8;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 4 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 4+2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseE(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/quad_e.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 12;

    static const size_t numCells = 6;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 4, 4, 4, 4 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 6+4 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseF(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/quad_f.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 12;

    static const size_t numCells = 6;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 4, 4, 4, 4 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 6+4 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseG(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/quad_g.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 14;

    static const size_t numCells = 7;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 4, 4, 4, 4, 4 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 0, 2, 2, 100, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 3+2, 6+4 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseH(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/quad_h.mesh";

    data->numFaults = 2;
    static const char* const faultSurfaceLabels[2] = { "faultA", "faultB" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[2] = { NULL, "faultB-edge" };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[2] = { 101, 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 22;

    static const size_t numCells = 14;
    data->numCells = numCells;

    static const int numCorners[numCells] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = {
        10, 10, 11, 10, 10, 11, 12, 12, 11,
        100, 100, 101, 101, 101 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 3;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 6+4, 8+6, 1 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "faultB", "faultA", "faultB-edge" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseI(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/quad_i.mesh";

    static const size_t numFaults = 1;
    data->numFaults = numFaults;
    static const char* const faultSurfaceLabels[numFaults] = { "fault", };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[numFaults] = { "edge" };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[numFaults] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 14;

    static const size_t numCells = 8;
    data->numCells = numCells;

    static const int numCorners[numCells] = {
        4, 4, 4, 4, 4, 4,
        4, 4, };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = {
        10, 10, 10, 11, 10, 11,
        100, 100, };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 3+2, 4+2, 1, 5+4 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output2", "output1", "edge", "fault"};
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseI


// End of file
