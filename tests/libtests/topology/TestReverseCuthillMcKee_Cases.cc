// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
//

#include <portinfo>

#include "TestReverseCuthillMcKee.hh" // Implementation of class methods

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace topology {
        class TestReverseCuthillMcKee_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestReverseCuthillMcKee_Cases {
public:

    // Data factory methods
    static TestReverseCuthillMcKee_Data* Tri_NoFault(void);

    static TestReverseCuthillMcKee_Data* Tri_Fault(void);

    static TestReverseCuthillMcKee_Data* Quad_NoFault(void);

    static TestReverseCuthillMcKee_Data* Quad_Fault(void);

    static TestReverseCuthillMcKee_Data* Tet_NoFault(void);

    static TestReverseCuthillMcKee_Data* Tet_Fault(void);

    static TestReverseCuthillMcKee_Data* Hex_NoFault(void);

    static TestReverseCuthillMcKee_Data* Hex_Fault(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestReverseCuthillMcKee::Tri_NoFault::testReorder", "[TestReverseCuthillMcKee][Tri][NoFault][testReorder]") {
    pylith::topology::TestReverseCuthillMcKee(pylith::topology::TestReverseCuthillMcKee_Cases::Tri_NoFault()).testReorder();
}
TEST_CASE("TestReverseCuthillMcKee::Tri_NoFault::testReorder", "[TestReverseCuthillMcKee][Tri][Fault][testReorder]") {
    pylith::topology::TestReverseCuthillMcKee(pylith::topology::TestReverseCuthillMcKee_Cases::Tri_Fault()).testReorder();
}

TEST_CASE("TestReverseCuthillMcKee::Quad_NoFault::testReorder", "[TestReverseCuthillMcKee][Quad][NoFault][testReorder]") {
    pylith::topology::TestReverseCuthillMcKee(pylith::topology::TestReverseCuthillMcKee_Cases::Quad_NoFault()).testReorder();
}
TEST_CASE("TestReverseCuthillMcKee::Quad_NoFault::testReorder", "[TestReverseCuthillMcKee][Quad][Fault][testReorder]") {
    pylith::topology::TestReverseCuthillMcKee(pylith::topology::TestReverseCuthillMcKee_Cases::Quad_Fault()).testReorder();
}

TEST_CASE("TestReverseCuthillMcKee::Tet_NoFault::testReorder", "[TestReverseCuthillMcKee][Tet][NoFault][testReorder]") {
    pylith::topology::TestReverseCuthillMcKee(pylith::topology::TestReverseCuthillMcKee_Cases::Tet_NoFault()).testReorder();
}
TEST_CASE("TestReverseCuthillMcKee::Tet_NoFault::testReorder", "[TestReverseCuthillMcKee][Tet][Fault][testReorder]") {
    pylith::topology::TestReverseCuthillMcKee(pylith::topology::TestReverseCuthillMcKee_Cases::Tet_Fault()).testReorder();
}

TEST_CASE("TestReverseCuthillMcKee::Hex_NoFault::testReorder", "[TestReverseCuthillMcKee][Hex][NoFault][testReorder]") {
    pylith::topology::TestReverseCuthillMcKee(pylith::topology::TestReverseCuthillMcKee_Cases::Hex_NoFault()).testReorder();
}
TEST_CASE("TestReverseCuthillMcKee::Hex_NoFault::testReorder", "[TestReverseCuthillMcKee][Hex][Fault][testReorder]") {
    pylith::topology::TestReverseCuthillMcKee(pylith::topology::TestReverseCuthillMcKee_Cases::Hex_Fault()).testReorder();
}

// ------------------------------------------------------------------------------------------------
pylith::topology::TestReverseCuthillMcKee_Data*
pylith::topology::TestReverseCuthillMcKee_Cases::Tri_NoFault(void) {
    TestReverseCuthillMcKee_Data* data = new TestReverseCuthillMcKee_Data();assert(data);

    data->filename = "data/reorder_tri3.mesh";
    data->faultLabel = NULL;

    return data;
}   // Tri_NoFault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestReverseCuthillMcKee_Data*
pylith::topology::TestReverseCuthillMcKee_Cases::Tri_Fault(void) {
    TestReverseCuthillMcKee_Data* data = new TestReverseCuthillMcKee_Data();assert(data);

    data->filename = "data/reorder_tri3.mesh";
    data->faultLabel = "fault";

    return data;
}   // Tri_Fault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestReverseCuthillMcKee_Data*
pylith::topology::TestReverseCuthillMcKee_Cases::Quad_NoFault(void) {
    TestReverseCuthillMcKee_Data* data = new TestReverseCuthillMcKee_Data();assert(data);

    data->filename = "data/reorder_quad4.mesh";
    data->faultLabel = NULL;

    return data;
}   // Quad_NoFault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestReverseCuthillMcKee_Data*
pylith::topology::TestReverseCuthillMcKee_Cases::Quad_Fault(void) {
    TestReverseCuthillMcKee_Data* data = new TestReverseCuthillMcKee_Data();assert(data);

    data->filename = "data/reorder_quad4.mesh";
    data->faultLabel = "fault";

    return data;
}   // Quad_Fault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestReverseCuthillMcKee_Data*
pylith::topology::TestReverseCuthillMcKee_Cases::Tet_NoFault(void) {
    TestReverseCuthillMcKee_Data* data = new TestReverseCuthillMcKee_Data();assert(data);

    data->filename = "data/reorder_tet4.mesh";
    data->faultLabel = NULL;

    return data;
}   // Tet_NoFault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestReverseCuthillMcKee_Data*
pylith::topology::TestReverseCuthillMcKee_Cases::Tet_Fault(void) {
    TestReverseCuthillMcKee_Data* data = new TestReverseCuthillMcKee_Data();assert(data);

    data->filename = "data/reorder_tet4.mesh";
    data->faultLabel = "fault";

    return data;
}   // Tet_Fault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestReverseCuthillMcKee_Data*
pylith::topology::TestReverseCuthillMcKee_Cases::Hex_NoFault(void) {
    TestReverseCuthillMcKee_Data* data = new TestReverseCuthillMcKee_Data();assert(data);

    data->filename = "data/reorder_hex8.mesh";
    data->faultLabel = NULL;

    return data;
}   // Hex_NoFault


// ------------------------------------------------------------------------------------------------
pylith::topology::TestReverseCuthillMcKee_Data*
pylith::topology::TestReverseCuthillMcKee_Cases::Hex_Fault(void) {
    TestReverseCuthillMcKee_Data* data = new TestReverseCuthillMcKee_Data();assert(data);

    data->filename = "data/reorder_hex8.mesh";
    data->faultLabel = "fault";

    return data;
}   // Hex_fault


// End of file
