// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/testing/testingfwd.hh" // forward declaration

#include "pylith/problems/problemsfwd.hh" // HOLDSA TimeDependent
#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh

#include "pylith/utils/petscfwd.h" // HASA PetscVec

class pylith::testing::MMSTest : public pylith::utils::GenericComponent {
    // PUBLIC TYPEDEFS ////////////////////////////////////////////////////////////////////////////
public:

    typedef PetscErrorCode (*solution_fn)(PetscInt /* dim */,
                                          PetscReal /* t */,
                                          const PetscReal /* x */[],
                                          PetscInt /* Nc */,
                                          PetscScalar /* u */[],
                                          void* /* ctx */);

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    MMSTest(void);

    /// Default destructor.
    ~MMSTest(void);

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

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Initialize objects for test.
    virtual void _initialize(void);

    /// Set exact solution and time derivative of solution in domain.
    virtual void _setExactSolution(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    pylith::problems::TimeDependent* _problem; ///< Time-dependent problem.
    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::topology::Field* _solution; ///< Solution field.
    PetscVec _solutionExactVec; ///< Global vector to use for exact solution.
    PetscVec _solutionDotExactVec; ///< Global vector to use for time derivative of exact solution.
    PylithReal _jacobianConvergenceRate; ///< Expected convergence rate for Jacobiab (when not linear).
    PylithReal _tolerance; ///< Tolerance for discretization and residual test.
    PylithReal _avgCellDim; ///< Average cell dimension.
    bool _isJacobianLinear; ///< Jacobian is should be linear.
    bool _allowZeroResidual; ///< Allow residual to be exactly zero.

}; // MMSTest

// End of file
