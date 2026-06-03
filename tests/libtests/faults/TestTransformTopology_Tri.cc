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
        class TestTransformTopology_Tri;
    }
}

class pylith::faults::TestTransformTopology_Tri {
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

private:

    TestTransformTopology_Tri(void); ///< Not implemented
}; // TestTransformTopology_Tri

// ------------------------------------------------------------------------------------------------
#include "catch2/catch_test_macros.hpp"

TEST_CASE("TestTransform_TriA", "[TestTransformTopology][Tri]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tri::caseA()).run();
}
TEST_CASE("TestTransform_TriB", "[TestTransformTopology][Tri]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tri::caseB()).run();
}
TEST_CASE("TestTransform_TriC", "[TestTransformTopology][Tri]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tri::caseC()).run();
}
TEST_CASE("TestTransform_TriD", "[TestTransformTopology][Tri]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tri::caseD()).run();
}
TEST_CASE("TestTransform_TriE", "[TestTransformTopology][Tri]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tri::caseE()).run();
}
TEST_CASE("TestTransform_TriF", "[TestTransformTopology][Tri]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tri::caseF()).run();
}
TEST_CASE("TestTransform_TriG", "[TestTransformTopology][Tri]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tri::caseG()).run();
}
TEST_CASE("TestTransform_TriH", "[TestTransformTopology][Tri]") {
    pylith::faults::TestTransformTopology(pylith::faults::TestTransformTopology_Tri::caseH()).run();
}

// ------------------------------------------------------------------------------------------------

pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tri::caseA(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tri_a.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4, 4+2, 2, 2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseA


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tri::caseB(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tri_b.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4, 4+2, 2, 2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseB


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tri::caseC(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tri_c.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4, 4+2, 2, 2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseC


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tri::caseD(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tri_d.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4, 6+4, 2, 4 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseD


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tri::caseE(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tri_e.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
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
    data->getMatId = pylith::faults::TestTransformTopology_Data::getMatIdDefault;

    static const size_t numGroups = 4;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+4, 6+4, 2, 4 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseE


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tri::caseF(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tri_f.mesh";

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

    static const size_t numCells = 5;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 3, 3, 3, 3, 4, };
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
    static const int groupSizes[numGroups] = { 4+2, 4+2, 2, 2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "fault_vertices", "output_faces", "fault_faces" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "face", "face" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseF


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tri::caseG(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tri_g.mesh";

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

    static const size_t numCells = 7;
    data->numCells = numCells;

    static const int numCorners[numCells] = { 3, 3, 3, 3, 3, 3, 4, };
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

    static const size_t numGroups = 6;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 5+3, 1, 3+2, 3, 2, 1 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "edge_vertices", "fault_vertices", "output_faces", "fault_faces", "fault_faces_edge_auto" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex", "face", "face", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseG


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tri::caseH(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tri_h.mesh";

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_faces" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
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
    data->getMatId = [](const int cell,
                        const int numNoncohesiveCells,
                        const double* xyz) {
                         int value = 100;
                         if (cell < numNoncohesiveCells) {
                             if (xyz[1] > 0.0) {
                                 value = 0;
                             } else if (xyz[1] > -1.0) {
                                 value = 2;
                             } else {
                                 value = 3;
                             } // if/else
                         } // if
                         return value;
                     };

    static const size_t numGroups = 6;
    data->numGroups = numGroups;
    static const int groupSizes[numGroups] = { 3+2, 2, 4+4, 2, 4, 2 }; // vertices + edges
    data->groupSizes = const_cast<int*>(groupSizes);
    static const char* groupNames[numGroups] = { "output_vertices", "edge_vertices", "fault_vertices", "output_faces", "fault_faces", "fault_faces_edge_auto" };
    data->groupNames = const_cast<char**>(groupNames);
    static const char* groupTypes[numGroups] = { "vertex", "vertex", "vertex", "face", "face", "vertex" };
    data->groupTypes = const_cast<char**>(groupTypes);

    return data;
} // caseH


// ------------------------------------------------------------------------------------------------
pylith::faults::TestTransformTopology_Data*
pylith::faults::TestTransformTopology_Tri::caseI(void) {
    pylith::faults::TestTransformTopology_Data* data = new pylith::faults::TestTransformTopology_Data();
    REQUIRE(data);

    data->filename = "data/tri_i.mesh";
    data->failureExpected = true;

    data->numFaults = 1;
    static const char* const faultSurfaceLabels[1] = { "fault_vertices" };
    data->faultSurfaceLabels = const_cast<const char**>(faultSurfaceLabels);
    static const char* const faultEdgeLabels[1] = { NULL };
    data->faultEdgeLabels = const_cast<const char**>(faultEdgeLabels);
    static const int interfaceIds[1] = { 100 };
    data->interfaceIds = const_cast<const int*>(interfaceIds);

    return data;
} // caseI


// End of file
