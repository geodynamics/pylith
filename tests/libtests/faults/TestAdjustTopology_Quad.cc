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

#if 0
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
#endif
TEST_CASE("TestTransform_QuadA", "[TestTransform][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseA()).run_transform();
}
TEST_CASE("TestTransform_QuadB", "[TestTransform][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseB()).run_transform();
}
TEST_CASE("TestTransform_QuadC", "[TestTransform][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseC()).run_transform();
}
TEST_CASE("TestTransform_QuadD", "[TestTransform][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseD()).run_transform();
}
TEST_CASE("TestTransform_QuadE", "[TestTransform][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseE()).run_transform();
}
TEST_CASE("TestTransform_QuadF", "[TestTransform][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseF()).run_transform();
}
TEST_CASE("TestTransform_QuadG", "[TestTransform][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseG()).run_transform();
}
TEST_CASE("TestTransform_QuadH", "[TestTransform][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseH()).run_transform();
}
TEST_CASE("TestTransform_QuadI", "[TestTransform][Quad]") {
    pylith::faults::TestAdjustTopology(pylith::faults::TestAdjustTopology_Quad::caseI()).run_transform();
}

// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseA(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    REQUIRE(data);

    data->filename = "data/quad_a.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 4+2, 2, 2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseB(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    REQUIRE(data);

    data->filename = "data/quad_b.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 4+2, 2, 2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseB


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseC(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    REQUIRE(data);

    data->filename = "data/quad_c.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 4+2, 2, 2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseC


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseD(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    REQUIRE(data);

    data->filename = "data/quad_d.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 4+2, 2, 2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseD


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseE(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    REQUIRE(data);

    data->filename = "data/quad_e.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 6+4, 2, 4 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseE


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseF(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    REQUIRE(data);

    data->filename = "data/quad_f.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestAdjustTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 4+2, 6+4, 2, 4 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseG(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    REQUIRE(data);

    data->filename = "data/quad_g.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = [](const int cell,
                        const int numNoncohesiveCells,
                        const double* xyz) {
                         int value = 100;
                         if (cell < numNoncohesiveCells) {
                             value = (xyz[1] > 0.0) ? 0 : 2;
                         }
                         return value;
                     };

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 3+2, 6+4, 2, 4 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseH(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    REQUIRE(data);

    data->filename = "data/quad_h.mesh";

    data->numFaults = 2;
    static const char* const faultSurfaceLabels[2] = { "faultA_faces", "faultB_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[2] = { NULL, NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[2] = { 100, 101 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    data->cellDim = 2;
    data->spaceDim = 2;
    data->numVertices = 22;

    static const size_t numCells = 14;
    data->numCells = numCells;

    static const int numCorners[numCells] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, };
    data->numCorners = const_cast<int*>(numCorners);
    data->getMatId = [](const int cell,
                        const int numNoncohesiveCells,
                        const double* xyz) {
                         int value = -999;
                         if (cell < numNoncohesiveCells) {
                             if (xyz[0] > 1.0) {
                                 value = 11;
                             } else if (xyz[1] > -1.0) {
                                 value = 10;
                             } else {
                                 value = 12;
                             }
                         } else {
                             const int numCohesiveCellsFaultA = 2;
                             value = (cell > numNoncohesiveCells + numCohesiveCellsFaultA) ? interfaceIds[1] : interfaceIds[0];
                         }
                         return value;
                     };

    static const size_t numGroups = 6;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 8+6, 6+4, 2, 6, 4, 1 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "faultA_vertices", "faultB_vertices", "faultB-edge_vertices", "faultA_faces", "faultB_faces", "faultB_faces_edge_auto" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex", "face", "face", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestAdjustTopology_Data*
pylith::faults::TestAdjustTopology_Quad::caseI(void) {
    pylith::faults::TestAdjustTopology_Data* data = new pylith::faults::TestAdjustTopology_Data();
    REQUIRE(data);

    data->filename = "data/quad_i.mesh";

    static const size_t numFaults = 1;
    data->numFaults = numFaults;
    static const char* const faultSurfaceLabels[numFaults] = { "fault_faces", };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[numFaults] = { NULL };
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
    data->getMatId = [](const int cell,
                        const int numNoncohesiveCells,
                        const double* xyz) {
                         int value = 100;
                         if (cell < numNoncohesiveCells) {
                             value = (xyz[0] < 0.0 || xyz[1] > 0.0) ? 10 : 11;
                         }
                         return value;
                     };

    static const size_t numGroups = 8;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 3+2, 4+2, 1, 5+4, 2, 2, 4, 1 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output2_vertices", "output1_vertices", "edge_vertices", "fault_vertices", "output2_faces", "output1_faces", "fault_faces", "fault_faces_edge_auto"};
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex", "vertex", "face", "face", "face", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseI


// End of file
