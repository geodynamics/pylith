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

#include <portinfo>

#include "MMSTest.hh" // implementation of class methods

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "petscts.h" // USES PetscTS

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include <cassert> // USES assert()

extern "C" {
    PetscErrorCode DMSNESCheckFromOptions_Internal(PetscSNES snes,
                                                   PetscDM dm,
                                                   PetscVec u,
                                                   PetscErrorCode(**exactFuncs)(PetscInt, PetscReal, const PetscReal x[], PetscInt, PetscScalar *u, void *ctx),
                                                   void **ctxs);

} // extern "C"

// ---------------------------------------------------------------------------------------------------------------------
/// Setup testing data.
void
pylith::testing::MMSTest::setUp(void) {
    _problem = new pylith::problems::TimeDependent;CPPUNIT_ASSERT(_problem);
    _mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    _solution = NULL;
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
/// Tear down testing data.
void
pylith::testing::MMSTest::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    journal::debug_t debug(_problem->PyreComponent::getName());
    debug.deactivate(); // DEBUGGING

    delete _problem;_problem = NULL;
    delete _mesh;_mesh = NULL;
    delete _solution;_solution = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
/// Verify discretization can represent solution field.
void
pylith::testing::MMSTest::testDiscretization(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();

#if 0
    // Call function for discretization test refactored from DMSNESCheckFromOptions_Internal().
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad @matt Implement test.", false);
#else
    CPPUNIT_ASSERT(_problem);
    CPPUNIT_ASSERT(_problem->_ts);
    PetscTS ts = _problem->_ts;
    PetscDM dm = NULL;
    PetscSNES snes = NULL;
    PetscErrorCode err = 0;
    err = TSGetDM(ts, &dm);CPPUNIT_ASSERT(!err);
    err = TSGetSNES(ts, &snes);CPPUNIT_ASSERT(!err);

    CPPUNIT_ASSERT(_solution);
    err = DMSNESCheckFromOptions_Internal(snes, dm, _solution->scatterVector("global"), NULL, NULL);CPPUNIT_ASSERT(!err);
#endif

    PYLITH_METHOD_END;
} // testDiscretization


// ---------------------------------------------------------------------------------------------------------------------
/// Verify residual evaluated for solution is below specified tolerance.
void
pylith::testing::MMSTest::testResidual(void) {
    PYLITH_METHOD_BEGIN;

    // Call function for residual test refactored from DMSNESCheckFromOptions_Internal().
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad @matt Implement test.", false);

    PYLITH_METHOD_END;
} // testResidual


// ---------------------------------------------------------------------------------------------------------------------
// Verify Jacobian via Taylor series.
//
// || F(\vec{s} + \epsilon \vec{v}) - F(\vec{s} - \epsilon J \vec{v} || < \epsilon**2
void
pylith::testing::MMSTest::testJacobianTaylorSeries(void) {
    PYLITH_METHOD_BEGIN;

    // Call function for Jacobian test refactored from DMSNESCheckFromOptions_Internal().
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad @matt Implement test.", false);

    PYLITH_METHOD_END;
} // testJacobianTaylorSeries


// ---------------------------------------------------------------------------------------------------------------------
// Test Jacobian using finite differences.
void
pylith::testing::MMSTest::testJacobianFiniteDiff(void) {
    PYLITH_METHOD_BEGIN;

    // Call SNESSolve with appropriate options.
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad @matt Implement test.", false);

#if 0
    err = SNESSetFromOptions("-snes_test_jacobian");
    err = SNESComputeJacobian(snes, solution, jacobian, preconditioner);
#endif

    PYLITH_METHOD_END;
} // testJacobianFiniteDiff


// ---------------------------------------------------------------------------------------------------------------------
// Initialize objects for test.
void
pylith::testing::MMSTest::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_problem);
    CPPUNIT_ASSERT(_solution);

    _problem->setSolverType(pylith::problems::Problem::NONLINEAR);
    _problem->setMaxTimeSteps(1);
    _problem->preinitialize(*_mesh);
    _problem->verifyConfiguration();

    _setExactSolutionBC();
    _problem->initialize();
    _setExactSolution();

#if 0
    err = SNESSetSolution(snes, u);CPPUNIT_ASSERT(!err);

    err = TSGetDM(ts, &dm);CPPUNIT_ASSERT(!err);
    err = TSGetSNES(ts, &snes);CPPUNIT_ASSERT(!err);
    err = DMSNESCheckFromOptions_Internal(snes, dm, sol, exactFuncs, ctxs);CPPUNIT_ASSERT(!err);
#endif

    PYLITH_METHOD_END;
} // _initialize


// End of file
