// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file libsrc/problems/MMSTest.hh
 *
 * @brief Object for using the Method of Manufactured Solutions to test residual and Jacobian calculations.
 */

#if !defined(pylith_testing_mmstest_hh)
#define pylith_testing_mmstest_hh

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent
#include <cppunit/extensions/HelperMacros.h> /// ISA Cppunit::TestFixture

#include "pylith/testing/testingfwd.hh" // forward declaration

#include "pylith/problems/problemsfwd.hh" // HOLDSA TimeDependent
#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh

#include "pylith/utils/petscfwd.h" // HASA PetscVec

class pylith::testing::MMSTest : public pylith::utils::GenericComponent,
    public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(MMSTest);

    CPPUNIT_TEST(testDiscretization);
    CPPUNIT_TEST(testResidual);
    CPPUNIT_TEST(testJacobianTaylorSeries);
    CPPUNIT_TEST(testJacobianFiniteDiff);

    CPPUNIT_TEST_SUITE_END_ABSTRACT();

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    virtual void setUp(void);

    /// Tear down testing data.
    virtual void tearDown(void);

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

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Initialize objects for test.
    virtual void _initialize(void);

    /// Set exact solution in domain.
    virtual void _setExactSolution(void) = 0;

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::problems::TimeDependent *_problem; ///< Time-dependent problem.
    pylith::topology::Mesh *_mesh; ///< Finite-element mesh.
    pylith::topology::Field *_solution; ///< Solution field.
    PetscVec _solutionExactVec; ///< Global vector to use for exact solution.
    PetscVec _solutionDotExactVec; ///< Global vector to use for time derivative of exact solution.
    PylithReal _jacobianConvergenceRate; ///< Expected convergence rate for Jacobiab (when not linear).
    bool _isJacobianLinear; ///< Jacobian is should be linear.
    bool _disableFiniteDifferenceCheck; ///< Flag to indicate not to perform finite-difference check of Jacobian.
    bool _allowZeroResidual; ///< Allow residual to be exactly zero.

}; // MMSTest

#endif // pylith_testing_mmstest_hh

// End of file
