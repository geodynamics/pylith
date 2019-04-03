// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/problems/MMSTest.hh
 *
 * @brief Object for using the Method of Manufactured Solutions to test residual and Jacobian calculations.
 */

#if !defined(pylith_testing_mmstest_hh)
#define pylith_testing_mmstest_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/testing/testingfwd.hh" // forward declaration

#include "pylith/problems/problemsfwd.hh" // HOLDSA TimeDependent
#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh

class pylith::testing::MMSTest : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(MMSTest);

    CPPUNIT_TEST(testDiscretization);
    CPPUNIT_TEST(testResidual);
    CPPUNIT_TEST(testJacobianTaylorSeries);
    CPPUNIT_TEST(testJacobianFiniteDiff);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Verify discretization can represent solution field.
    void testDiscretization(void);

    /// Verify residual evaluated for solution is below specified tolerance.
    void testResidual(void);

    /** Verify Jacobian via Taylor series.
     *
     * || F(\vec{s} + \epsilon \vec{v}) - F(\vec{s} - \epsilon J \vec{v} || < \epsilon**2
     */
    void testJacobianTaylorSeries(void);

    /** Test Jacobian using finite differences.
     *
     * Compare computed Jacobian against one computed via finite-differences.
     */
    void testJacobianFiniteDiff(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::topology::TimeDependent* _problem; ///< Time-dependent problem.
    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    MMSTest(const MMSTest&); ///< Not implemented
    const MMSTest& operator=(const MMSTest&); ///< Not implemented

}; // MMSTest

#endif // pylith_testing_mmstest_hh

// End of file
