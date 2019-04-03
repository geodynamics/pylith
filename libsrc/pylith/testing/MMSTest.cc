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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "petscts.h" // USES PetscTS

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
/// Setup testing data.
void
pylith::testing::MMSTest::setUp(void) {
    _problem = NULL;
    _mesh = NULL;

    err = VecDuplicate(u, &sol);PYLITH_CHECK_ERROR(err);
    err = VecCopy(u, sol);PYLITH_CHECK_ERROR(err);
    err = TSSetSolution(ts, u);PYLITH_CHECK_ERROR(err);
    err = TSSetUp(ts);PYLITH_CHECK_ERROR(err);
    err = SNESSetSolution(snes, u);PYLITH_CHECK_ERROR(err);

#if 0
    err = TSGetDM(ts, &dm);PYLITH_CHECK_ERROR(err);
    err = TSGetSNES(ts, &snes);PYLITH_CHECK_ERROR(err);
    err = DMSNESCheckFromOptions_Internal(snes, dm, sol, exactFuncs, ctxs);PYLITH_CHECK_ERROR(err);

https: // bitbucket.org/petsc/petsc/src/4db401632e1443d930e54234abe1842a0bc30a81/src/ts/utils/dmplexts.c?fileviewer=file-view-default#lines-266
https: // bitbucket.org/petsc/petsc/src/4db401632e1443d930e54234abe1842a0bc30a81/src/snes/utils/dmplexsnes.c?fileviewer=file-view-default#lines-2470
https: // bitbucket.org/petsc/petsc/src/4db401632e1443d930e54234abe1842a0bc30a81/src/snes/examples/tutorials/ex17.c?fileviewer=file-view-default#lines-352
#endif

} // setUp


// ---------------------------------------------------------------------------------------------------------------------
/// Tear down testing data.
void
pylith::testing::MMSTest::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    err = VecDestroy(&sol);PYLITH_CHECK_ERROR(err);

    delete _problem;_problem = NULL;
    delete _mesh;_mesh = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
/// Verify discretization can represent solution field.
void
pylith::testing::MMSTest::testDiscretization(void) {
    PYLITH_METHOD_BEGIN;

    // Call function for discretization test refactored from DMSNESCheckFromOptions_Internal().
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad @matt Implement test.");

    PYLITH_METHOD_END;
} // testDiscretization


// ---------------------------------------------------------------------------------------------------------------------
/// Verify residual evaluated for solution is below specified tolerance.
void
pylith::testing::MMSTest::testResidual(void) {
    PYLITH_METHOD_BEGIN;

    // Call function for residual test refactored from DMSNESCheckFromOptions_Internal().
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad @matt Implement test.");

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
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad @matt Implement test.");

    PYLITH_METHOD_END;
} // testJacobianTaylorSeries


// ---------------------------------------------------------------------------------------------------------------------
// Test Jacobian using finite differences.
void
pylith::testing::MMSTest::testJacobianFiniteDiff(void) {
    PYLITH_METHOD_BEGIN;

    // Call SNESSolve with appropriate options.
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad @matt Implement test.");

#if 0
    err = SNESSetFromOptions("-snes_test_jacobian");
    err = SNESComputeJacobian(snes, solution, jacobian, preconditioner);
#endif

    PYLITH_METHOD_END;
} // testJacobianFiniteDiff


// ---------------------------------------------------------------------------------------------------------------------
// Initialize.
void
pylith::problems::MMSTest::initialize(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize()");

    Problem::initialize();

    // Set solution to initial conditions.
    assert(_solution);
    _solution->zeroLocal();
    const size_t numIC = _ic.size();
    for (size_t i = 0; i < numIC; ++i) {
        assert(_ic[i]);
        _ic[i]->setValues(_solution, *_normalizer);
    } // for

#if USE_SNES
    PetscErrorCode err = SNESDestroy(&_snes);PYLITH_CHECK_ERROR(err);assert(!_snes);
    const pylith::topology::Mesh& mesh = _solution->mesh();
    err = SNESCreate(mesh.comm(), &_snes);PYLITH_CHECK_ERROR(err);assert(_snes);
    err = SNESSetType(_snes, TSBEULER);PYLITH_CHECK_ERROR(err); // Backward Euler is default time stepping method.
    err = SNESSetFromOptions(_snes);PYLITH_CHECK_ERROR(err);
    err = SNESSetApplicationContext(_snes, (void*)this);PYLITH_CHECK_ERROR(err);

    PYLITH_COMPONENT_DEBUG("Setting SNES initial conditions using global vector for solution.");
    PetscVec solutionVec = NULL;
    err = DMCreateGlobalVector(_solution->dmMesh(), &solutionVec);PYLITH_CHECK_ERROR(err);
    _solution->scatterLocalToVector(solutionVec);
    err = SNESSetSolution(_snes, solutionVec);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&solutionVec);PYLITH_CHECK_ERROR(err);

    PYLITH_COMPONENT_DEBUG("Setting PetscTS callbacks prestep(), poststep(), computeRHSJacobian(), and computeRHSFunction().");
    err = SNESSetRHSJacobian(_snes, NULL, NULL, computeRHSJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
    err = SNESSetRHSFunction(_snes, NULL, computeRHSResidual, (void*)this);PYLITH_CHECK_ERROR(err);

    if (IMPLICIT == _formulationType) {
        PYLITH_COMPONENT_DEBUG("Setting PetscTS callbacks computeLHSJacobian(), and computeLHSFunction().");
        err = TSSetIFunction(_snes, NULL, computeLHSResidual, (void*)this);PYLITH_CHECK_ERROR(err);
        err = TSSetIJacobian(_snes, NULL, NULL, computeLHSJacobian, (void*)this);PYLITH_CHECK_ERROR(err);
    } // if

    // Setup time stepper.
    err = TSSetUp(_snes);PYLITH_CHECK_ERROR(err);

    // Setup field to hold inverse of lumped LHS Jacobian (if explicit).
    if (EXPLICIT == _formulationType) {
        PYLITH_COMPONENT_DEBUG("Setting up field for inverse of lumped LHS Jacobian.");

        delete _jacobianLHSLumpedInv;_jacobianLHSLumpedInv = new pylith::topology::Field(_solution->mesh());assert(_jacobianLHSLumpedInv);
        _jacobianLHSLumpedInv->cloneSection(*_solution);
    } // if
#endif

    PYLITH_METHOD_END;
} // initialize


// End of file
