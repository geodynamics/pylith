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

#include "TestTransformTopology.hh" // USES TestTransformTopology_Data

namespace pylith {
    namespace faults {
        class TestTransformTopology_Hex;
    }
}

class pylith::faults::TestTransformTopology_Hex {
public:

    // Data factory methods
    static TestTransformTopology_Data* caseA(void);

    static TestTransformTopology_Data* caseB(void);

    static TestTransformTopology_Data* caseC(void);

    static TestTransformTopology_Data* caseD(void);

    static TestTransformTopology_Data* caseE(void);

    static TestTransformTopology_Data* caseF(void);

    static TestTransformTopology_Data* caseG(void);

    static TestTransformTopology_Data* caseH(void);

    static TestTransformTopology_Data* caseI(void);

    static TestTransformTopology_Data* caseJ(void);

private:

    TestTransformTopology_Hex(void); ///< Not implemented
}; // TestTransformTopology_Hex

// ------------------------------------------------------------------------------------------------
#include "catch2/catch_test_macros.hpp"

TEST_CASE("TestTransform_HexA", "[TestTransform][Hex]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Hex::caseA()).run();
}
TEST_CASE("TestTransform_HexB", "[TestTransform][Hex]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Hex::caseB()).run();
}
TEST_CASE("TestTransform_HexC", "[TestTransform][Hex]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Hex::caseC()).run();
}
TEST_CASE("TestTransform_HexD", "[TestTransform][Hex]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Hex::caseD()).run();
}
TEST_CASE("TestTransform_HexE", "[TestTransform][Hex]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Hex::caseE()).run();
}
TEST_CASE("TestTransform_HexF", "[TestTransform][Hex]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Hex::caseF()).run();
}
TEST_CASE("TestTransform_HexG", "[TestTransform][Hex]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Hex::caseG()).run();
}
TEST_CASE("TestTransform_HexH", "[TestTransform][Hex]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Hex::caseH()).run();
}
TEST_CASE("TestTransform_HexI", "[TestTransform][Hex]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Hex::caseI()).run();
}
TEST_CASE("TestTransform_HexJ", "[TestTransform][Hex]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Hex::caseJ()).run();
}

// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Hex::caseA(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/hex_a.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2, 2, 8 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Hex::caseB(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/hex_b.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2, 2, 8 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseB


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Hex::caseC(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/hex_c.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2, 2, 8 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseC


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Hex::caseD(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/hex_d.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2, 2, 8 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseD


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Hex::caseE(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/hex_e.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2, 2, 8 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseE


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Hex::caseF(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/hex_e.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 8+8+2, 2, 8 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseF


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Hex::caseG(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/hex_g.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 12+14+4, 2, 14 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseG


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Hex::caseH(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/hex_h.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+8+2, 12+14+4, 2, 14 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseH


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Hex::caseI(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/hex_i.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 39+3;

    static const size_t numCells = 14;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = [](const int cell,
                        const int numNoncohesiveCells,
                        const double* xyz) {
                         int value = 100;
                         if (cell < numNoncohesiveCells) {
                             value = (xyz[2] < 0.0) ? 0 : 2;
                         }
                         return value;
                     };

    static const size_t numGroups = 5;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 10+13+4, 9+12+4, 4, 10, 2 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces", "fault_faces_edge_auto" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseI


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Hex::caseJ(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/hex_j.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 3;
    data->spaceDim = 3;
    data->numVertices = 62;

    static const size_t numCells = 22;
    data->numCells = numCells;

    static const int numCorners[numCells] = {
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6,
    };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = [](const int cell,
                        const int numNoncohesiveCells,
                        const double* xyz) {
                         int value = 100;
                         if (cell < numNoncohesiveCells) {
                             value = (xyz[2] > -10.0) ? 10 : 11;
                         }
                         return value;
                     };

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 22+31+10, 10, 8, 3 }; // vertices + edges + faces
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "output_faces", "fault_faces", "fault_faces_edge_auto" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "face", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseJ


// End of file
