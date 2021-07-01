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

/**
 * @file tests/libtests/problems/TestProgressMonitor.hh
 *
 * @brief C++ TestProgressMonitor object.
 *
 * C++ unit testing for ProgressMonitor.
 */

#if !defined(pylith_problems_testprogressmonitor_hh)
#define pylith_problems_testprogressmonitor_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/testing/testingfwd.hh" // HOLDSA ProgressMonitorStub

/// Namespace for pylith package
namespace pylith {
    namespace problems {
        class TestProgressMonitor;
    } // problems
} // pylith

class pylith::problems::TestProgressMonitor : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestProgressMonitor);

    CPPUNIT_TEST(testAccessors);
    CPPUNIT_TEST(testOpenClose);
    CPPUNIT_TEST(testUpdate);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test get/setUpdatePercent() and get/setFilename().
    void testAccessors(void);

    /// Test open() and close().
    void testOpenClose(void);

    /// Test update().
    void testUpdate(void);

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::problems::ProgressMonitorStub* _monitor; ///< Test subject.

}; // class TestProgressMonitor

#endif // pylith_problems_testprogressmonitor_hh

// End of file
