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
 * @file tests/libtests/problems/TestObserversPhysics.hh
 *
 * @brief C++ TestObserversPhysics object.
 *
 * C++ unit testing for Observers of solution.
 */

#if !defined(pylith_problems_testobserversphysics_hh)
#define pylith_problems_testobserversphysics_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/problems/problemsfwd.hh" // HOLDSA ObserversPhysics

/// Namespace for pylith package
namespace pylith {
    namespace problems {
        class TestObserversPhysics;
    } // problems
} // pylith

class pylith::problems::TestObserversPhysics : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestObserversPhysics);

    CPPUNIT_TEST(testRegisterObserver);
    CPPUNIT_TEST(testRemoveObserver);
    CPPUNIT_TEST(testSetPhysicsImplementation);
    CPPUNIT_TEST(testTimeScale);
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

    /// Test setPhysicsImplemetation().
    void testSetPhysicsImplementation(void);

    /// Test setgetTimeScale().
    void testTimeScale(void);

    /// Test verifyObservers().
    void testVerifyObservers(void);

    /// Test notifyObservers().
    void testNotifyObservers(void);

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::problems::ObserversPhysics* _observers; ///< Test subject.

}; // class TestObserversPhysics

#endif // pylith_problems_testobserversphysics_hh

// End of file
