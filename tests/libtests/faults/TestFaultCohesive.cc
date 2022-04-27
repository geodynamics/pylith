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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/Mesh.hh" // USES Mesh::cells_label_name
#include "pylith/testing/FaultCohesiveStub.hh" // USES FaultCohesiveStub
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

namespace pylith {
    namespace faults {
        class TestFaultCohesive;
    } // faults
} // pylith

class pylith::faults::TestFaultCohesive : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestFaultCohesive);

    CPPUNIT_TEST(testAccessors);
    // adjust topology tested by TestAdjustTopology

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Test set/getCohesiveLabelValue(), set/getSurfaceLabelName(), set/getBuriedEdgesLabelName(), setRefDir1/2().
    void testAccessors(void);

}; // class TestFaultCohesive

// ------------------------------------------------------------------------------------------------
// Test set/getCohesiveLabelValue(), set/getSurfaceLabelName(), set/getBuriedEdgesLabelName(), setRefDir1/2().
void
pylith::faults::TestFaultCohesive::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    FaultCohesiveStub fault;

    const std::string& cohesiveLabelName = "cohesive-id";
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default cohesive label name.", std::string(pylith::topology::Mesh::cells_label_name), std::string(fault.getCohesiveLabelName()));
    fault.setCohesiveLabelName(cohesiveLabelName.c_str());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in cohesive label name.", cohesiveLabelName, std::string(fault.getCohesiveLabelName()));

    const int cohesiveLabelValue = 23;
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default cohesive label value.", 100, fault.getCohesiveLabelValue());
    fault.setCohesiveLabelValue(cohesiveLabelValue);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in cohesive label value.", cohesiveLabelValue, fault.getCohesiveLabelValue());

    const std::string& surfaceLabel = "surface label";
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default surface label name.",
                                 std::string(""), std::string(fault.getSurfaceLabelName()));
    fault.setSurfaceLabelName(surfaceLabel.c_str());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in surface label name.",
                                 surfaceLabel, std::string(fault.getSurfaceLabelName()));

    const int surfaceValue = 4;
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default surface label value.", 1, fault.getSurfaceLabelValue());
    fault.setSurfaceLabelValue(surfaceValue);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in surface label value.", surfaceValue, fault.getSurfaceLabelValue());

    const std::string& edgeLabel = "edge label";
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default edge label name.",
                                 std::string(""), std::string(fault.getBuriedEdgesLabelName()));
    fault.setBuriedEdgesLabelName(edgeLabel.c_str());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in edge label name.",
                                 edgeLabel, std::string(fault.getBuriedEdgesLabelName()));

    const int edgeValue = 4;
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default edge label value.", 1, fault.getBuriedEdgesLabelValue());
    fault.setBuriedEdgesLabelValue(edgeValue);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in edge label value.", edgeValue, fault.getBuriedEdgesLabelValue());

    const int dim = 3;
    const PylithReal default1[3] = { 0.0, 0.0, 1.0 };
    for (int i = 0; i < dim; ++i) {
        CPPUNIT_ASSERT_EQUAL(default1[i], fault._refDir1[i]);
    } // for
    const PylithReal dir1[3] = { 0.0, 2.0, -1.0 };
    fault.setRefDir1(dir1);
    for (int i = 0; i < dim; ++i) {
        CPPUNIT_ASSERT_EQUAL(dir1[i]/sqrt(2*2+1*1), fault._refDir1[i]);
    } // for

    const PylithReal default2[3] = { 0.0, 1.0, 0.0 };
    for (int i = 0; i < dim; ++i) {
        CPPUNIT_ASSERT_EQUAL(default2[i], fault._refDir2[i]);
    } // for
    const PylithReal dir2[3] = { 2.0, 1.0, 2.0 };
    fault.setRefDir2(dir2);
    for (int i = 0; i < dim; ++i) {
        CPPUNIT_ASSERT_EQUAL(dir2[i]/sqrt(2*2+1*1+2*2), fault._refDir2[i]);
    } // for

    PYLITH_METHOD_END;
} // testAccessors


// End of file
