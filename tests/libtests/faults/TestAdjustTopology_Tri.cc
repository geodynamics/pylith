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
        class TestAdjustTopology_Tri;
    }
}

class pylith::faults::TestAdjustTopology_Tri {
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

    TestAdjustTopology_Tri(void); ///< Not implemented
}; // TestAdjustTopology_Tri

// ------------------------------------------------------------------------------------------------
#include "catch2/catch_test_macros.hpp"

TEST_CASE("TestAdjustTopology_TriA", "[TestAdjustTopology][Tri]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tri::caseA()).run();
}
TEST_CASE("TestAdjustTopology_TriB", "[TestAdjustTopology][Tri]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tri::caseB()).run();
}
TEST_CASE("TestAdjustTopology_TriC", "[TestAdjustTopology][Tri]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tri::caseC()).run();
}
TEST_CASE("TestAdjustTopology_TriD", "[TestAdjustTopology][Tri]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tri::caseD()).run();
}
TEST_CASE("TestAdjustTopology_TriE", "[TestAdjustTopology][Tri]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tri::caseE()).run();
}
TEST_CASE("TestAdjustTopology_TriF", "[TestAdjustTopology][Tri]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tri::caseF()).run();
}
TEST_CASE("TestAdjustTopology_TriG", "[TestAdjustTopology][Tri]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tri::caseG()).run();
}
TEST_CASE("TestAdjustTopology_TriH", "[TestAdjustTopology][Tri]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tri::caseH()).run();
}
TEST_CASE("TestAdjustTopology_TriI", "[TestAdjustTopology][Tri]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tri::caseI()).run();
}

// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tri::caseA(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/tri_a.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 6;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 3, 3, 4 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4, 4+2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tri::caseB(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/tri_b.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 6;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 3, 3, 4 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4, 4+2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseB


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tri::caseC(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/tri_c.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 6;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 3, 3, 4 };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4, 4+2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseC


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tri::caseD(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/tri_d.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 9;

    static const size_t numCells = 6;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 3, 3, 3, 3, 4, 4, };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4, 6+4 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseD


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tri::caseE(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/tri_e.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 9;

    static const size_t numCells = 6;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 3, 3, 3, 3, 4, 4, };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 0, 0, 100, 100 };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 2;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4, 6+4 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseE


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tri::caseF(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/tri_f.mesh";

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

    static const size_t numCells = 5;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 3, 3, 3, 3, 4, };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 0, 2, 2, 100, };
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
} // caseF


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tri::caseG(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/tri_g.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { "edge" };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 8;

    static const size_t numCells = 7;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 3, 3, 3, 3, 3, 3, 4, };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = { 0, 2, 0, 2, 0, 2, 100, };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 3;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+3, 1, 3+2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "edge", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseG


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tri::caseH(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/tri_h.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { "edge" };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 11;

    static const size_t numCells = 13;
    data->numCells = numCells;

    static const int numCorners[numCells] = {
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        4, 4, };
    data->numCorners = const_cast<int*>(numCorners);
    static const int materialIds[numCells] = {
        0, 0, 0, 0,
        2, 2, 3, 3, 3, 2, 3,
        100, 100,
    };
    data->materialIds = const_cast<int*>(materialIds);

    static const size_t numGroups = 3;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 3+2, 2, 4+4, }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output", "edge", "fault" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseH


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tri::caseI(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    assert(data);

    data->filename = "data/tri_i.mesh";
    data->failureExpected = true;

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { "edge" };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    return data;
} // caseI


// End of file
