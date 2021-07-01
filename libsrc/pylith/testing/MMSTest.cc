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
#include <iostream> // USES std::out

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::testing::MMSTest::setUp(void) {
    GenericComponent::setName("mmstest"); // Override in child class for finer control of journal output.
    _problem = new pylith::problems::TimeDependent;
    CPPUNIT_ASSERT(_problem);
    _mesh = new pylith::topology::Mesh();
    CPPUNIT_ASSERT(_mesh);
    _solution = NULL;
    _solutionExactVec = NULL;
    _solutionDotExactVec = NULL;
    _jacobianConvergenceRate = 0.0;
    _isJacobianLinear = false;
    _disableFiniteDifferenceCheck = false;
    _allowZeroResidual = false;
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::testing::MMSTest::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    pythia::journal::debug_t debug(GenericComponent::getName());
    debug.deactivate(); // DEBUGGING

    delete _problem;_problem = NULL;
    delete _mesh;_mesh = NULL;
    delete _solution;_solution = NULL;

    PetscErrorCode err;
    err = VecDestroy(&_solutionExactVec);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_solutionDotExactVec);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Verify discretization can represent solution field.
void
pylith::testing::MMSTest::testDiscretization(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();

    CPPUNIT_ASSERT(_problem);
    CPPUNIT_ASSERT(_solution);
    PetscErrorCode err = 0;
    const PylithReal tolerance = -1.0, t = 0.0;
    const pylith::string_vector subfieldNames = _solution->getSubfieldNames();
    const size_t numSubfields = subfieldNames.size();
    pylith::real_array error(numSubfields);
    err = DMSNESCheckDiscretization(_problem->getPetscSNES(), _problem->getPetscDM(), t, _solutionExactVec,
                                    tolerance, &error[0]);
    CPPUNIT_ASSERT(!err);

    bool fail = false;
    std::ostringstream msg;
    for (size_t i_field = 0; i_field < numSubfields; ++i_field) {
        msg << "Discretization test failed for subfield(s): ";
        if (error[i_field] > 1.0e-10) {
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
// Verify residual evaluated for solution is below specified tolerance.
void
pylith::testing::MMSTest::testResidual(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = 0;

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        err = PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2");
        CPPUNIT_ASSERT(!err);
        err = PetscOptionsSetValue(NULL, "-dm_plex_print_l2", "2");
        CPPUNIT_ASSERT(!err);
    } // if

    _initialize();

    CPPUNIT_ASSERT(_solution);
    if (debug.state()) {
        _solution->view("Solution field layout", pylith::topology::Field::VIEW_LAYOUT);
    } // if

    CPPUNIT_ASSERT(_problem);
    CPPUNIT_ASSERT(_solutionExactVec);
    CPPUNIT_ASSERT(_solutionDotExactVec);
    const PylithReal tolerance = -1.0;
    PylithReal norm = 0.0;
    err = DMTSCheckResidual(_problem->getPetscTS(), _problem->getPetscDM(), _problem->getStartTime(), _solutionExactVec,
                            _solutionDotExactVec, tolerance, &norm);
    if (debug.state()) {
        _solution->view("Solution field");
    } // if
    if (!_allowZeroResidual) {
        CPPUNIT_ASSERT_MESSAGE("L2 normal of residual is exactly zero, which suggests suspicious case with all residuals "
                               "entries exactly zero.",
                               norm > 0.0);
    } // if
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test of F(s) - G(s) == 0 failed.", 0.0, norm, 1.0e-10);

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
    CPPUNIT_ASSERT(_solutionExactVec);
    CPPUNIT_ASSERT(_solutionDotExactVec);
    PetscErrorCode err = 0;
    const PylithReal tolerance = -1.0;
    PetscBool isLinear = PETSC_FALSE;
    PylithReal convergenceRate = 0.0;
    err = DMTSCheckJacobian(_problem->getPetscTS(), _problem->getPetscDM(), _problem->getStartTime(), _solutionExactVec,
                            _solutionDotExactVec, tolerance, &isLinear, &convergenceRate);
    CPPUNIT_ASSERT(!err);

    if (_isJacobianLinear) {
        CPPUNIT_ASSERT_EQUAL_MESSAGE("Expected linear Jacobian.", PETSC_TRUE, isLinear);
    } else {
        const PylithReal tolerance = 1.0e-3;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Error in convergence rate for Jacobian.",
                                             _jacobianConvergenceRate, convergenceRate, tolerance);
    } // if/else

    PYLITH_METHOD_END;
} // testJacobianTaylorSeries


// ---------------------------------------------------------------------------------------------------------------------
// Test Jacobian using finite differences.
void
pylith::testing::MMSTest::testJacobianFiniteDiff(void) {
    PYLITH_METHOD_BEGIN;

    if (_disableFiniteDifferenceCheck) {
        PYLITH_JOURNAL_ERROR("Skipping Jacobian finite-difference check. Test disabled.");
        PYLITH_METHOD_END;
    } // if

    PetscErrorCode err = 0;
    err = PetscOptionsSetValue(NULL, "-ts_max_snes_failures", "1");CPPUNIT_ASSERT(!err);
    err = PetscOptionsSetValue(NULL, "-ts_error_if_step_fails", "false");CPPUNIT_ASSERT(!err);
    _initialize();

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        err = PetscOptionsSetValue(NULL, "-snes_test_jacobian_view", "");CPPUNIT_ASSERT(!err);
    } // if
    err = PetscOptionsSetValue(NULL, "-snes_test_jacobian", "1.0e-6");CPPUNIT_ASSERT(!err);
    err = PetscOptionsSetValue(NULL, "-snes_error_if_not_converged", "false");CPPUNIT_ASSERT(!err);
    err = SNESSetFromOptions(_problem->getPetscSNES());CPPUNIT_ASSERT(!err);

    CPPUNIT_ASSERT(_problem);
    CPPUNIT_ASSERT(_solution);

    _problem->solve();
    std::cout << "IMPORTANT: You must check the Jacobian values printed here manually!\n"
              << "           They should be O(1.0e-6) or smaller." << std::endl;
    err = PetscOptionsClearValue(NULL, "-snes_test_jacobian");CPPUNIT_ASSERT(!err);
    err = PetscOptionsClearValue(NULL, "-snes_test_jacobian_view");CPPUNIT_ASSERT(!err);

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
    TSSetUp(_problem->getPetscTS());
    _setExactSolution();

    // Global vectors to use for analytical solution in MMS tests.
    PetscErrorCode err = VecDuplicate(_solution->getGlobalVector(), &_solutionExactVec);CPPUNIT_ASSERT(!err);
    err = VecDuplicate(_solutionExactVec, &_solutionDotExactVec);CPPUNIT_ASSERT(!err);

    PYLITH_METHOD_END;
} // _initialize


// End of file
