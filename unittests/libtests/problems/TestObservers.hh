// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/feassemble/TestObservers.hh
 *
 * @brief C++ TestObservers object.
 *
 * C++ unit testing for Observers.
 */

#if !defined(pylith_feassemble_testobservers_hh)
#define pylith_feassemble_testobservers_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/feassemble/feassemblefwd.hh" // HOLDSA Observers

/// Namespace for pylith package
namespace pylith {
    namespace feassemble {
        class TestObservers;
    } // feassemble
} // pylith

class pylith::feassemble::TestObservers : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestObservers);

    CPPUNIT_TEST(testRegisterObserver);
    CPPUNIT_TEST(testRemoveObserver);
    CPPUNIT_TEST(testVerifyObservers);
    CPPUNIT_TEST(testNotifyObservers);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test registerObserver().
    void testRegisterObserver(void);

    /// Test removeObserver().
    void testRemoveObserver(void);

    /// Test verifyObservers().
    void testVerifyObservers(void);

    /// Test notifyObservers().
    void testNotifyObservers(void);

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::feassemble::Observers* _observers; ///< Test subject.

}; // class TestObservers

#endif // pylith_feassemble_testobservers_hh

// End of file
