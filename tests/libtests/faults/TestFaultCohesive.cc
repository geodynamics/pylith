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

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "tests/src/FaultCohesiveStub.hh" // USES FaultCohesiveStub
#include "pylith/topology/Mesh.hh" // USES Mesh::cells_label_name
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

namespace pylith {
    namespace faults {
        class TestFaultCohesive;
    } // faults
} // pylith

class pylith::faults::TestFaultCohesive : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Test set/getCohesiveLabelValue(), set/getSurfaceLabelName(), set/getBuriedEdgesLabelName(), setRefDir1/2().
    static
    void testAccessors(void);

}; // class TestFaultCohesive

// ------------------------------------------------------------------------------------------------
// Test set/getCohesiveLabelValue(), set/getSurfaceLabelName(), set/getBuriedEdgesLabelName(), setRefDir1/2().
void
pylith::faults::TestFaultCohesive::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    FaultCohesiveStub fault;

    const std::string& cohesiveLabelName = "cohesive-id";
    CHECK(std::string(pylith::topology::Mesh::cells_label_name) == std::string(fault.getCohesiveLabelName()));
    fault.setCohesiveLabelName(cohesiveLabelName.c_str());
    CHECK(cohesiveLabelName == std::string(fault.getCohesiveLabelName()));

    const int cohesiveLabelValue = 23;
    CHECK(100 == fault.getCohesiveLabelValue());
    fault.setCohesiveLabelValue(cohesiveLabelValue);
    CHECK(cohesiveLabelValue == fault.getCohesiveLabelValue());

    const std::string& surfaceLabel = "surface label";
    CHECK(std::string("") == std::string(fault.getSurfaceLabelName()));
    fault.setSurfaceLabelName(surfaceLabel.c_str());
    CHECK(surfaceLabel == std::string(fault.getSurfaceLabelName()));

    const int surfaceValue = 4;
    CHECK(1 == fault.getSurfaceLabelValue());
    fault.setSurfaceLabelValue(surfaceValue);
    CHECK(surfaceValue == fault.getSurfaceLabelValue());

    const std::string& edgeLabel = "edge label";
    CHECK(std::string("") == std::string(fault.getBuriedEdgesLabelName()));
    fault.setBuriedEdgesLabelName(edgeLabel.c_str());
    CHECK(edgeLabel == std::string(fault.getBuriedEdgesLabelName()));

    const int edgeValue = 4;
    CHECK(1 == fault.getBuriedEdgesLabelValue());
    fault.setBuriedEdgesLabelValue(edgeValue);
    CHECK(edgeValue == fault.getBuriedEdgesLabelValue());

    const int dim = 3;
    const PylithReal default1[3] = { 0.0, 0.0, 1.0 };
    for (int i = 0; i < dim; ++i) {
        CHECK(default1[i] == fault._refDir1[i]);
    } // for

    const double tolerance = 1.0e-6;

    const PylithReal dir1[3] = { 0.0, 2.0, -1.0 };
    fault.setRefDir1(dir1);
    for (int i = 0; i < dim; ++i) {
        CHECK_THAT(fault._refDir1[i], Catch::Matchers::WithinAbs(dir1[i]/sqrt(2*2+1*1), tolerance));
    } // for

    const PylithReal default2[3] = { 0.0, 1.0, 0.0 };
    for (int i = 0; i < dim; ++i) {
        CHECK_THAT(fault._refDir2[i], Catch::Matchers::WithinAbs(default2[i], tolerance));
    } // for
    const PylithReal dir2[3] = { 2.0, 1.0, 2.0 };
    fault.setRefDir2(dir2);
    for (int i = 0; i < dim; ++i) {
        CHECK_THAT(fault._refDir2[i], Catch::Matchers::WithinAbs(dir2[i]/sqrt(2*2+1*1+2*2), tolerance));
    } // for

    PYLITH_METHOD_END;
} // testAccessors


// ------------------------------------------------------------------------------------------------
TEST_CASE("FaultCohesive::testAccessors", "[FaultCohesive]") {
    pylith::faults::TestFaultCohesive::testAccessors();
}

// End of file
