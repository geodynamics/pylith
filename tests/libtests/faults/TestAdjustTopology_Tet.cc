// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestAdjustTopology.hh" // USES TestAdjustTopology_Data

namespace pylith {
    namespace faults {
        class TestAdjustTopology_Tet;
    }
}

class pylith::faults::TestAdjustTopology_Tet {
public:

    // Data factory methods
    static TestAdjustTopology_Data* caseA(void);

    static TestAdjustTopology_Data* caseB(void);

    static TestAdjustTopology_Data* caseC(void);

    static TestAdjustTopology_Data* caseD(void);

    static TestAdjustTopology_Data* caseF(void);

    static TestAdjustTopology_Data* caseG(void);

    static TestAdjustTopology_Data* caseH(void);

    static TestAdjustTopology_Data* caseI(void);

    static TestAdjustTopology_Data* caseJ(void);

    static TestAdjustTopology_Data* caseK(void);

private:

    TestAdjustTopology_Tet(void); ///< Not implemented
}; // TestAdjustTopology_Tet

// ------------------------------------------------------------------------------------------------
#include "catch2/catch_test_macros.hpp"

#if 0
TEST_CASE("TestAdjustTopology_TetA", "[TestAdjustTopology][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseA()).run();
}
TEST_CASE("TestAdjustTopology_TetB", "[TestAdjustTopology][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseB()).run();
}
TEST_CASE("TestAdjustTopology_TetC", "[TestAdjustTopology][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseC()).run();
}
TEST_CASE("TestAdjustTopology_TetD", "[TestAdjustTopology][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseD()).run();
}
TEST_CASE("TestAdjustTopology_TetF", "[TestAdjustTopology][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseF()).run();
}
TEST_CASE("TestAdjustTopology_TetG", "[TestAdjustTopology][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseG()).run();
}
TEST_CASE("TestAdjustTopology_TetH", "[TestAdjustTopology][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseH()).run();
}
TEST_CASE("TestAdjustTopology_TetI", "[TestAdjustTopology][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseI()).run();
}
TEST_CASE("TestAdjustTopology_TetJ", "[TestAdjustTopology][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseJ()).run();
}
TEST_CASE("TestAdjustTopology_TetK", "[TestAdjustTopology][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseK()).run();
}
#endif
TEST_CASE("TestTransform_TetA", "[TestTransform][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseA()).run_transform();
}
TEST_CASE("TestTransform_TetB", "[TestTransform][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseB()).run_transform();
}
TEST_CASE("TestTransform_TetC", "[TestTransform][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseC()).run_transform();
}
TEST_CASE("TestTransform_TetD", "[TestTransform][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseD()).run_transform();
}
TEST_CASE("TestTransform_TetF", "[TestTransform][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseF()).run_transform();
}
TEST_CASE("TestTransform_TetG", "[TestTransform][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseG()).run_transform();
}
TEST_CASE("TestTransform_TetH", "[TestTransform][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseH()).run_transform();
}
TEST_CASE("TestTransform_TetI", "[TestTransform][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseI()).run_transform();
}
TEST_CASE("TestTransform_TetJ", "[TestTransform][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseJ()).run_transform();
}
TEST_CASE("TestTransform_TetK", "[TestTransform][Tet]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Tet::caseK()).run_transform();
}

// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tet::caseA(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

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
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tet::caseB(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

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
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tet::caseC(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

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
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tet::caseD(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

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
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tet::caseF(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

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
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tet::caseG(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

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
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tet::caseH(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

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
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tet::caseI(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

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
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tet::caseJ(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
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
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Tet::caseK(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
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
