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
#include "pylith/utils/array.hh" // USES real_array, string_vector
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include <cassert> // USES assert()

#if 0
extern "C" {
    PetscErrorCode DMSNESCheckDiscretization(SNES snes,
                                             DM dm,
                                             Vec u,
                                             PetscErrorCode(**exactFuncs)(PetscInt, PetscReal, const PetscReal x[], PetscInt, PetscScalar *u, void *ctx),
                                             void **ctxs,
                                             PetscReal tol);

    PetscErrorCode DMSNESCheckResidual(SNES snes,
                                       DM dm,
                                       Vec u,
                                       PetscReal tol);

    PetscErrorCode DMSNESCheckJacobian(SNES snes,
                                       DM dm,
                                       Vec u,
                                       PetscReal tol);

} // extern "C"
#endif

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

    CPPUNIT_ASSERT(_problem);
    CPPUNIT_ASSERT(_solution);
    PetscErrorCode err = 0;
    const PylithReal tolerance = -1.0;
    const pylith::string_vector subfieldNames = _solution->subfieldNames();
    const size_t numSubfields = subfieldNames.size();
    pylith::real_array error(numSubfields);
    err = DMSNESCheckDiscretization(_problem->getPetscSNES(), _problem->getPetscDM(), _solution->scatterVector("mmstest"),
                                    NULL, NULL, tolerance, &error[0]);CPPUNIT_ASSERT(!err);

    bool fail = false;
    std::ostringstream msg;
    for (size_t i_field = 0; i_field < numSubfields; ++i_field) {
        msg << "Discretization test failed for subfield(s): ";
        if (error[i_field] > 1.0e-6) {
            fail = true;
            msg << " " << subfieldNames[i_field] << " (" << error[i_field] << ")";
        } // if
    } // for
    if (fail) {
        CPPUNIT_ASSERT_MESSAGE(msg.str(), fail);
    } // if

    PYLITH_METHOD_END;
} // testDiscretization


// ---------------------------------------------------------------------------------------------------------------------
/// Verify residual evaluated for solution is below specified tolerance.
void
pylith::testing::MMSTest::testResidual(void) {
    PYLITH_METHOD_BEGIN;

    PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2"); // :DEBUG:
    PetscOptionsSetValue(NULL, "-dm_plex_print_l2", "2"); // :DEBUG:

    _initialize();

    CPPUNIT_ASSERT(_problem);
    CPPUNIT_ASSERT(_solution);
    PetscErrorCode err = 0;
    const PylithReal tolerance = -1.0;
    PylithReal norm = 0.0;
    err = DMSNESCheckResidual(_problem->getPetscSNES(), _problem->getPetscDM(), _solution->scatterVector("mmstest"),
                              tolerance, &norm);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_MESSAGE("L2 normal of residual is exactly zero, which suggests suspicious case with all residuals "
                           "entries exactly zero.", norm > 0.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test of F(s) - G(s) == 0 failed.", 0.0, norm, 1.0e-6);

    PYLITH_METHOD_END;
} // testResidual


// ---------------------------------------------------------------------------------------------------------------------
// Verify Jacobian via Taylor series.
//
// || F(\vec{s} + \epsilon \vec{v}) - F(\vec{s} - \epsilon J \vec{v} || < \epsilon**2
void
pylith::testing::MMSTest::testJacobianTaylorSeries(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();

    CPPUNIT_ASSERT(_problem);
    CPPUNIT_ASSERT(_solution);
    PetscErrorCode err = 0;
    const PylithReal tolerance = -1.0;
    PetscBool isLinear = PETSC_FALSE;
    PylithReal convergenceRate = 0.0;
    err = DMSNESCheckJacobian(_problem->getPetscSNES(), _problem->getPetscDM(), _solution->scatterVector("mmstest"),
                              tolerance, &isLinear, &convergenceRate);CPPUNIT_ASSERT(!err);

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

    _problem->initialize();
    _setExactSolution();

    // Global vector to use for solution in MMS tests.
    _solution->createScatter(_solution->dmMesh(), "mmstest");

    PYLITH_METHOD_END;
} // _initialize


// End of file
