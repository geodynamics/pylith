// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "tests/src/MMSTest.hh" // implementation of class methods
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData
#include "pylith/utils/PetscOptions.hh" // USES PetscOptions

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field

#include "petscts.h" // USES PetscTS

#include "pylith/utils/error.hh" // USES PylithCallPetsc()
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
    _tolerance(1.0e-8),
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

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PylithCallPetsc(PetscOptionsSetValue(NULL, "-dm_plex_print_l2", "2"));
    } // if

    _initialize();
    const pylith::topology::Field* solution = _problem->getSolution();assert(solution);
    if (debug.state()) {
        solution->view("Solution field layout", pylith::topology::Field::VIEW_LAYOUT);
    } // if

    const PylithReal ignoreTolerance = -1.0;
    const pylith::string_vector subfieldNames = solution->getSubfieldNames();
    const size_t numSubfields = subfieldNames.size();
    pylith::real_array error(numSubfields);
    PylithCallPetsc(DMSNESCheckDiscretization(_problem->getPetscSNES(), _problem->getPetscDM(), _problem->getStartTime(),
                                              _solutionExactVec, ignoreTolerance, &error[0]));

    if (debug.state()) {
        solution->view("Solution field");
    } // if

    bool fail = false;
    std::ostringstream msg;
    msg << "Discretization test failed for subfield(s):\n";
    for (size_t i_field = 0; i_field < numSubfields; ++i_field) {
        if (error[i_field] > _tolerance) {
            fail = true;
            msg << "    " << subfieldNames[i_field] << " (error=" << error[i_field] << ", tolerance="<< _tolerance << ")\n";
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

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PylithCallPetsc(PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2"));
        PylithCallPetsc(PetscOptionsSetValue(NULL, "-dm_plex_print_l2", "2"));
    } // if

    _initialize();
    if (_problem->getFormulation() == pylith::problems::Physics::DYNAMIC) {
        _problem->_integrationData->removeField(pylith::feassemble::IntegrationData::lumped_jacobian_inverse);
    } // if

    assert(_solutionExactVec);
    assert(_solutionDotExactVec);
    PylithReal ignoreTolerance = -1.0;
    PylithReal norm = 0.0;
    PylithCallPetsc(DMTSCheckResidual(_problem->getPetscTS(), _problem->getPetscDM(), _problem->getStartTime(), _solutionExactVec,
                                      _solutionDotExactVec, ignoreTolerance, &norm));
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

    PetscBool isLinear = PETSC_FALSE;
    PylithReal convergenceRate = 0.0;
    PylithCallPetsc(DMTSCheckJacobian(_problem->getPetscTS(), _problem->getPetscDM(), _problem->getStartTime(), _solutionExactVec,
                                      _solutionDotExactVec, _tolerance, &isLinear, &convergenceRate));

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

    PylithCallPetsc(PetscOptionsSetValue(NULL, "-ts_max_snes_failures", "1"));
    PylithCallPetsc(PetscOptionsSetValue(NULL, "-ts_error_if_step_fails", "false"));
    _initialize();

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PylithCallPetsc(PetscOptionsSetValue(NULL, "-snes_test_jacobian_view", ""));
    } // if
    PylithCallPetsc(PetscOptionsSetValue(NULL, "-snes_error_if_not_converged", "false"));
    PylithCallPetsc(SNESSetFromOptions(_problem->getPetscSNES()));

    PetscReal jacobianNorm = 1.0e+20, diffNorm = 1.0e+20;
    PylithCallPetsc(TSSetUp(_problem->getPetscTS()));
    PetscMat jacobian = PETSC_NULLPTR, jacobianPreconditioner = PETSC_NULLPTR;
    PylithCallPetsc(TSGetIJacobian(_problem->getPetscTS(), &jacobian, &jacobianPreconditioner, NULL, NULL));
    PylithCallPetsc(SNESComputeJacobian(_problem->getPetscSNES(), _solutionExactVec, jacobian, jacobianPreconditioner));
    PylithCallPetsc(SNESTestJacobian(_problem->getPetscSNES(), &jacobianNorm, &diffNorm));

    INFO("jacobianNorm=" << jacobianNorm << ", ||Code Jacobian - Finite Diff Jacobain||="<<diffNorm);
    const PetscReal jacobianTolerance = 20 * jacobianNorm * _tolerance;
    REQUIRE_THAT(jacobianNorm, !Catch::Matchers::WithinAbs(0.0, _tolerance));
    REQUIRE_THAT(diffNorm, Catch::Matchers::WithinAbs(0.0, jacobianTolerance));
    PylithCallPetsc(PetscOptionsClearValue(NULL, "-snes_test_jacobian_view"));

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

    if (_problem->getFormulation() == pylith::problems::Physics::DYNAMIC) {
        PylithCallPetsc(TSSetIFunction(_problem->getPetscTS(), NULL, pylith::problems::TimeDependent::computeLHSResidual,
                                       (void*)_problem));
    } // if

    // Global vectors to use for analytical solution in MMS tests.
    const pylith::topology::Field* solution = _problem->getSolution();assert(solution);
    PylithCallPetsc(VecDuplicate(solution->getGlobalVector(), &_solutionExactVec));
    PylithCallPetsc(VecDuplicate(_solutionExactVec, &_solutionDotExactVec));

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
