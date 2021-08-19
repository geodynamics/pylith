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

    /// Test set/getInterfaceId(), set/getSurfaceMarkerLabel(), set/getBuriedEdgesMarkerLabel(), setRefDir1/2().
    void testAccessors(void);

}; // class TestFaultCohesive

// ------------------------------------------------------------------------------------------------
// Test set/getInterfaceId(), set/getSurfaceMarkerLabel(), set/getBuriedEdgesMarkerLabel(), setRefDir1/2().
void
pylith::faults::TestFaultCohesive::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    FaultCohesiveStub fault;

    const int interfaceId = 23;
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default interface id.", 100, fault.getInterfaceId());
    fault.setInterfaceId(interfaceId);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in interface id.", interfaceId, fault.getInterfaceId());

    const std::string& surfaceLabel = "surface label";
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default surface marker label.",
                                 std::string(""), std::string(fault.getSurfaceMarkerLabel()));
    fault.setSurfaceMarkerLabel(surfaceLabel.c_str());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in surface marker label.",
                                 surfaceLabel, std::string(fault.getSurfaceMarkerLabel()));

    const std::string& edgeLabel = "edge label";
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default edge marker label.",
                                 std::string(""), std::string(fault.getBuriedEdgesMarkerLabel()));
    fault.setBuriedEdgesMarkerLabel(edgeLabel.c_str());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in edge marker label.",
                                 edgeLabel, std::string(fault.getBuriedEdgesMarkerLabel()));

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
