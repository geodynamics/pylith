// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestTransformTopology.hh" // USES TestTransformTopology_Data

namespace pylith {
    namespace faults {
        class TestTransformTopology_Tet;
    }
}

class pylith::faults::TestTransformTopology_Tet {
public:

    // Data factory methods
    static TestTransformTopology_Data* caseA(void);

    static TestTransformTopology_Data* caseB(void);

    static TestTransformTopology_Data* caseC(void);

    static TestTransformTopology_Data* caseD(void);

    static TestTransformTopology_Data* caseF(void);

    static TestTransformTopology_Data* caseG(void);

    static TestTransformTopology_Data* caseH(void);

    static TestTransformTopology_Data* caseI(void);

    static TestTransformTopology_Data* caseJ(void);

    static TestTransformTopology_Data* caseK(void);

private:

    TestTransformTopology_Tet(void); ///< Not implemented
}; // TestTransformTopology_Tet

// ------------------------------------------------------------------------------------------------
#include "catch2/catch_test_macros.hpp"

TEST_CASE("TestTransform_TetA", "[TestTransformTopology][Tet]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tet::caseA()).run();
}
TEST_CASE("TestTransform_TetB", "[TestTransformTopology][Tet]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tet::caseB()).run();
}
TEST_CASE("TestTransform_TetC", "[TestTransformTopology][Tet]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tet::caseC()).run();
}
TEST_CASE("TestTransform_TetD", "[TestTransformTopology][Tet]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tet::caseD()).run();
}
TEST_CASE("TestTransform_TetF", "[TestTransformTopology][Tet]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tet::caseF()).run();
}
TEST_CASE("TestTransform_TetG", "[TestTransformTopology][Tet]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tet::caseG()).run();
}
TEST_CASE("TestTransform_TetH", "[TestTransformTopology][Tet]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tet::caseH()).run();
}
TEST_CASE("TestTransform_TetI", "[TestTransformTopology][Tet]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tet::caseI()).run();
}
TEST_CASE("TestTransform_TetJ", "[TestTransformTopology][Tet]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tet::caseJ()).run();
}
TEST_CASE("TestTransform_TetK", "[TestTransformTopology][Tet]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tet::caseK()).run();
}

// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tet::caseA(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tet_a.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 8;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 5 };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4+1, 6+6+2, 1, 6 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tet::caseB(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tet_b.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 8;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 5 };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4+1, 6+6+2, 1, 6 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseB


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tet::caseC(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tet_c.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 8;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 5 };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4+1, 6+6+2, 1, 6 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseC


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tet::caseD(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tet_d.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 8;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 5 };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4+1, 6+6+2, 1, 6 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseD


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tet::caseF(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tet_f.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 8;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 5 };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4+1, 6+6+2, 1, 6 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseF


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tet::caseG(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tet_g.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 8;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 5 };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4+1, 6+6+2, 1, 6 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseG


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tet::caseH(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tet_h.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 8;

    static const size_t numCells = 3;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 5 };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4+1, 6+6+2, 1, 6 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseH


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tet::caseI(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tet_i.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 10;

    static const size_t numCells = 6;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 4, 4, 5, 5 };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 2+0+0, 8+10+4, 1, 10 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces"};
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseI


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tet::caseJ(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tet_j.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 12;

    static const size_t numCells = 7;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 4, 4, 4, 4, 4, 4, 5 };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = [](const int cell,
                        const int numNoncohesiveCells,
                        const double* xyz) {
                         int value = 100;
                         if (cell < numNoncohesiveCells) {
                             value = (xyz[0] < 0.0) ? 0 : 2;
                         }
                         return value;
                     };

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+10+4, 6+6+2, 2, 6 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseJ


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tet::caseK(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tet_k.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 234+16+4;

    static const size_t numCells = 833+40;
    data->numCells = numCells;

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
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   4, 4, 4,
        5, 5, 5, 5, 5,
        5, 5, 5, 5, 5,
        5, 5, 5, 5, 5,
        5, 5, 5, 5, 5,
        5, 5, 5, 5, 5,
        5, 5, 5, 5, 5,
        5, 5, 5, 5, 5,
        5, 5, 5, 5, 5,
    };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = [](const int cell,
                        const int numNoncohesiveCells,
                        const double* xyz) {
                         return (cell < numNoncohesiveCells) ? 10 : 100;
                     };

    static const size_t numGroups = 6;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 223+5+4, 10+8, 2*139-8-10, 66, 101+21, 8 }; // vertices + edges +
                                                                                           // faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "vertices_zpos", "fault_edge_vertices", "fault_vertices", "boundary_zpos", "fault_faces", "fault_faces_edge_auto" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex", "face", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseK


// End of file
