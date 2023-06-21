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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "tests/src/MMSTest.hh" // implementation of class methods
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData
#include "pylith/utils/PetscOptions.hh" // USES PetscOptions

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "petscts.h" // USES PetscTS

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/array.hh" // USES real_array, string_vector
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::testing::MMSTest::MMSTest(void) :
    _problem(new pylith::problems::TimeDependent),
    _mesh(new pylith::topology::Mesh()),
    _solution(NULL),
    _solutionExactVec(NULL),
    _solutionDotExactVec(NULL),
    _jacobianConvergenceRate(0.0),
    _tolerance(1.0e-9),
    _isJacobianLinear(false),
    _allowZeroResidual(false) {
    GenericComponent::setName("mmstest"); // Override in child class for finer control of journal output.

    assert(_problem);
    assert(_mesh);

    _problem->setPetscDefaults(pylith::utils::PetscDefaults::TESTING | pylith::utils::PetscDefaults::SOLVER);
} // setUp


// ------------------------------------------------------------------------------------------------
// Tear down testing data.
pylith::testing::MMSTest::~MMSTest(void) {
    VecDestroy(&_solutionExactVec);
    VecDestroy(&_solutionDotExactVec);

    delete _problem;_problem = NULL;
    delete _mesh;_mesh = NULL;
    delete _solution;_solution = NULL;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Verify discretization can represent solution field.
void
pylith::testing::MMSTest::testDiscretization(void) {
    PYLITH_METHOD_BEGIN;
    assert(_problem);

    _initialize();

    pythia::journal::debug_t debug(GenericComponent::getName());
    const pylith::topology::Field* solution = _problem->getSolution();assert(solution);
    if (debug.state()) {
        solution->view("Solution field layout", pylith::topology::Field::VIEW_LAYOUT);
    } // if

    PetscErrorCode err = 0;
    const PylithReal tolerance = -1.0;
    const pylith::string_vector subfieldNames = solution->getSubfieldNames();
    const size_t numSubfields = subfieldNames.size();
    pylith::real_array error(numSubfields);
    err = DMSNESCheckDiscretization(_problem->getPetscSNES(), _problem->getPetscDM(), _problem->getStartTime(),
                                    _solutionExactVec, tolerance, &error[0]);PYLITH_CHECK_ERROR(err);

    if (debug.state()) {
        solution->view("Solution field");
    } // if

    bool fail = false;
    std::ostringstream msg;
    for (size_t i_field = 0; i_field < numSubfields; ++i_field) {
        msg << "Discretization test failed for subfield(s): ";
        if (error[i_field] > _tolerance) {
            fail = true;
            msg << " " << subfieldNames[i_field] << " (" << error[i_field] << ")";
        } // if
    } // for
    if (fail) {
        FAIL(msg.str());
    } // if

    PYLITH_METHOD_END;
} // testDiscretization


// ---------------------------------------------------------------------------------------------------------------------
// Verify residual evaluated for solution is below specified tolerance.
void
pylith::testing::MMSTest::testResidual(void) {
    PYLITH_METHOD_BEGIN;
    assert(_problem);

    PetscErrorCode err = 0;
    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        err = PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2");PYLITH_CHECK_ERROR(err);
        err = PetscOptionsSetValue(NULL, "-dm_plex_print_l2", "2");PYLITH_CHECK_ERROR(err);
    } // if

    _initialize();
    if (_problem->getFormulation() == pylith::problems::Physics::DYNAMIC) {
        _problem->_integrationData->removeField(pylith::feassemble::IntegrationData::lumped_jacobian_inverse);
    } // if

    assert(_solutionExactVec);
    assert(_solutionDotExactVec);
    PylithReal ignoreTolerance = -1.0;
    PylithReal norm = 0.0;
    err = DMTSCheckResidual(_problem->getPetscTS(), _problem->getPetscDM(), _problem->getStartTime(), _solutionExactVec,
                            _solutionDotExactVec, ignoreTolerance, &norm);PYLITH_CHECK_ERROR(err);
    if (!_allowZeroResidual && (0.0 == norm)) {
        FAIL("L2 normal of residual is exactly zero, which suggests suspicious case with all residual "
             "entries exactly zero.");
    } // if

    INFO("|G(s) - F(s)| == " << norm);
    REQUIRE_THAT(norm, Catch::Matchers::WithinAbs(0.0, _tolerance));

    PYLITH_METHOD_END;
} // testResidual


// ---------------------------------------------------------------------------------------------------------------------
// Verify Jacobian via Taylor series.
//
// || F(\vec{s} + \epsilon \vec{v}) - F(\vec{s} - \epsilon J \vec{v} || < \epsilon**2
void
pylith::testing::MMSTest::testJacobianTaylorSeries(void) {
    PYLITH_METHOD_BEGIN;
    assert(_problem);
    _initialize();

    assert(_solutionExactVec);
    assert(_solutionDotExactVec);
    PetscErrorCode err = 0;
    const PylithReal tolerance = -1.0;
    PetscBool isLinear = PETSC_FALSE;
    PylithReal convergenceRate = 0.0;
    err = DMTSCheckJacobian(_problem->getPetscTS(), _problem->getPetscDM(), _problem->getStartTime(), _solutionExactVec,
                            _solutionDotExactVec, tolerance, &isLinear, &convergenceRate);PYLITH_CHECK_ERROR(err);

    if (_isJacobianLinear) {
        REQUIRE(isLinear == PETSC_TRUE);
    } else {
        INFO("Convergence rate for Jacobian is " << convergenceRate);
        REQUIRE_THAT(convergenceRate, Catch::Matchers::WithinAbs(_jacobianConvergenceRate, 1.0e-3));
    } // if/else

    PYLITH_METHOD_END;
} // testJacobianTaylorSeries


// ---------------------------------------------------------------------------------------------------------------------
// Test Jacobian using finite differences.
void
pylith::testing::MMSTest::testJacobianFiniteDiff(void) {
    PYLITH_METHOD_BEGIN;
    assert(_problem);

    PetscErrorCode err = 0;
    err = PetscOptionsSetValue(NULL, "-ts_max_snes_failures", "1");PYLITH_CHECK_ERROR(err);
    err = PetscOptionsSetValue(NULL, "-ts_error_if_step_fails", "false");PYLITH_CHECK_ERROR(err);
    _initialize();

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        err = PetscOptionsSetValue(NULL, "-snes_test_jacobian_view", "");PYLITH_CHECK_ERROR(err);
    } // if
    err = PetscOptionsSetValue(NULL, "-snes_test_jacobian", "1.0e-6");PYLITH_CHECK_ERROR(err);
    err = PetscOptionsSetValue(NULL, "-snes_error_if_not_converged", "false");PYLITH_CHECK_ERROR(err);
    err = SNESSetFromOptions(_problem->getPetscSNES());PYLITH_CHECK_ERROR(err);

    _problem->solve();
    INFO("IMPORTANT: You must check the Jacobian values printed here manually!\n"
         << "           They should be O(1.0e-6) or smaller.\n");
    err = PetscOptionsClearValue(NULL, "-snes_test_jacobian");PYLITH_CHECK_ERROR(err);
    err = PetscOptionsClearValue(NULL, "-snes_test_jacobian_view");PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // testJacobianFiniteDiff


// ---------------------------------------------------------------------------------------------------------------------
// Initialize objects for test.
void
pylith::testing::MMSTest::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    assert(_problem);

    _problem->setSolverType(pylith::problems::Problem::NONLINEAR);
    _problem->setMaxTimeSteps(1);
    _problem->preinitialize(*_mesh);
    _problem->verifyConfiguration();

    _problem->initialize();
    TSSetUp(_problem->getPetscTS());
    _setExactSolution();

    PetscErrorCode err = PETSC_SUCCESS;
    if (_problem->getFormulation() == pylith::problems::Physics::DYNAMIC) {
        err = TSSetIFunction(_problem->getPetscTS(), NULL, pylith::problems::TimeDependent::computeLHSResidual,
                             (void*)_problem);PYLITH_CHECK_ERROR(err);
    } // if

    // Global vectors to use for analytical solution in MMS tests.
    const pylith::topology::Field* solution = _problem->getSolution();assert(solution);
    err = VecDuplicate(solution->getGlobalVector(), &_solutionExactVec);PYLITH_CHECK_ERROR(err);
    err = VecDuplicate(_solutionExactVec, &_solutionDotExactVec);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        const pylith::topology::Mesh& mesh = solution->getMesh();
        mesh.view();
        std::string options = std::string(":") + std::string(GenericComponent::getName()) + std::string("_mesh.tex:ascii_latex");
        mesh.view(options.c_str());
    } // if

    PYLITH_METHOD_END;
} // _initialize


// End of file
