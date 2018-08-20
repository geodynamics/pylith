// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
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

#include "Integrator.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/Observers.hh" // USES Observers
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> \
    // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::feassemble::Integrator::Integrator(pylith::problems::Physics* const physics) :
    _physics(physics),
    _auxField(NULL),
    _derivedField(NULL),
    _observers(NULL),
    _logger(NULL),
    _needNewRHSJacobian(true),
    _needNewLHSJacobian(true)
{}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::feassemble::Integrator::~Integrator(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::Integrator::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete _auxField;_auxField = NULL;
    delete _derivedField;_derivedField = NULL;
    delete _observers;_observers = NULL;
    delete _logger;_logger = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary field.
const pylith::topology::Field*
pylith::feassemble::Integrator::getAuxiliaryField(void) const {
    return _auxField;
} // getAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Get derived field.
const pylith::topology::Field*
pylith::feassemble::Integrator::getDerivedField(void) const {
    return _derivedField;
} // getDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Check whether RHS Jacobian needs to be recomputed.
bool
pylith::feassemble::Integrator::needNewRHSJacobian(void) const {
    return _needNewRHSJacobian;
} // needNewRHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Check whether LHS Jacobian needs to be recomputed.
bool
pylith::feassemble::Integrator::needNewLHSJacobian(void) const {
    return _needNewLHSJacobian;
} // needNewLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::Integrator::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("intialize(solution="<<solution.label()<<")");

    const pylith::topology::Mesh& integrationDomainMesh = getIntegrationDomainMesh();

    delete _auxField;_auxField = _physics->createAuxiliaryField(solution, integrationDomainMesh);
    delete _derivedField;_derivedField = _physics->createDerivedField(solution, integrationDomainMesh);
    _observers = _physics->getObservers(); // Memory managed by Python

    //_auxField->view("MATERIAL AUXILIARY FIELD"); // :DEBUG: TEMPORARY
    const bool infoOnly = true;
    _observers->notifyObservers(0.0, 0, solution, infoOnly);

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::feassemble::Integrator::prestep(const PylithReal t,
                                        const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("prestep(t="<<t<<", dt="<<dt<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // prestep


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary fields at end of time step.
void
pylith::feassemble::Integrator::poststep(const PylithReal t,
                                         const PylithInt tindex,
                                         const PylithReal dt,
                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("poststep(t="<<t<<", dt="<<dt<<") empty method");

    _updateStateVars(t, dt, solution);

    const bool infoOnly = false;
    assert(_observers);
    _observers->notifyObservers(t, tindex, solution, infoOnly);

    PYLITH_METHOD_END;
} // poststep


// ---------------------------------------------------------------------------------------------------------------------
// Set constants used in finite-element kernels (point-wise functions).
void
pylith::feassemble::Integrator::_setKernelConstants(const pylith::topology::Field& solution,
                                                    const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setFEConstants(solution="<<solution.label()<<", dt="<<dt<<")");

    assert(_physics);
    const pylith::real_array& constants = _physics->getKernelConstants(dt);

    PetscDS prob = NULL;
    PetscDM dmSoln = solution.dmMesh();assert(dmSoln);

    PetscErrorCode err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);assert(prob);
    if (constants.size() > 0) {
        err = PetscDSSetConstants(prob, constants.size(), const_cast<double*>(&constants[0]));PYLITH_CHECK_ERROR(err);
    } else {
        err = PetscDSSetConstants(prob, 0, NULL);PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // _setKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::Integrator::_updateStateVars(const PylithReal t,
                                                 const PylithReal dt,
                                                 const pylith::topology::Field& solution) {
#if 0
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("updateStateVars(t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    if (0 == _updateStateVarsKernels.size()) {
        PYLITH_METHOD_END;
    } // if

    assert(_auxField);

    PetscErrorCode err;
#if 0
    PetscDM auxDM = _auxField->dmMesh(), stateVarDM, dms[2], superDM;
    PetscIS stateVarIS, *superIS;
    PetscVec stateVarVec, A, locA;

    // HAPPEN ONCE
    // Create subDM holding only the state vars
    err = DMCreateSubDM(auxDM, numStateSubfields, stateSubfields, &stateVarIS, &stateVarDM);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(stateVarDM, &stateVarVec);PYLITH_CHECK_ERROR(err);
    // Create superDM of {state vars, solution}
    dms[0] = stateVarDM;dms[1] = solution->dmMesh();
    err = DMCreateSuperDM(dms, 2, &superIS, &superDM);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(superDM, &A);PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(superDM, &locA);PYLITH_CHECK_ERROR(err);

    // Copy state vars from auxiliary vars
    err = VecISCopy(stateVarVec, stateVarIS, SCATTER_FORWARD, _auxField->globalVector());PYLITH_CHECK_ERROR(err);
    // Copy current state vars and solution into superDM space
    err = VecISCopy(A, subis[0], SCATTER_FORWARD, stateVarVec);PYLITH_CHECK_ERROR(err);
    err = VecISCopy(A, subis[1], SCATTER_FORWARD, solution->globalVector());PYLITH_CHECK_ERROR(err);
    // Move superDM data to a local vector
    err = DMGlobalToLocalBegin(superDM, A, INSERT_VALUES, locA);PYLITH_CHECK_ERROR(err);
    err = DMGlobalToLocalEnd(superDM, A, INSERT_VALUES, locA);PYLITH_CHECK_ERROR(err);
    // Attach superDM data as auxiliary for updateStateVar kernel
    err = PetscObjectCompose((PetscObject) auxDM, "dmAux", (PetscObject) superDM);CHKERRQ(ierr);
    err = PetscObjectCompose((PetscObject) auxDM, "A",     (PetscObject) locA);CHKERRQ(ierr);
#else
    PetscDM dmState = _auxField->dmMesh();
    err = PetscObjectCompose((PetscObject) dmState, "dmAux", (PetscObject) solution.dmMesh());PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmState, "A", (PetscObject) solution.localVector());PYLITH_CHECK_ERROR(err);
#endif

    _setFEConstants(*_auxField, dt);

    // Set update kernel for each auxiliary subfield.
    const pylith::string_vector& subfieldNames = _auxField->subfieldNames();
    const size_t numSubfields = subfieldNames.size();
    PetscPointFunc* stateVarsKernels = (numSubfields > 0) ? new PetscPointFunc[numSubfields] : NULL;
    // By default, set all auxiliary subfield update kernels to NULL.
    for (size_t i = 0; i < numSubfields; ++i) {
        stateVarsKernels[i] = NULL;
    } // for
    for (UpdateStateVarsMap::const_iterator iter = _updateStateVarsKernels.begin(); iter != _updateStateVarsKernels.end(); ++iter) {
        const pylith::topology::Field::SubfieldInfo& sinfo = _auxField->subfieldInfo(iter->first.c_str());
        stateVarsKernels[sinfo.index] = iter->second;
    } // for

    err = DMProjectFieldLocal(dmState, t, _auxField->localVector(), stateVarsKernels, INSERT_VALUES, _auxField->localVector());PYLITH_CHECK_ERROR(err);
    delete[] stateVarsKernels;stateVarsKernels = NULL;

#if 0
    // HAPPEN ONCE
    // Destroy superDM stuff
    for (i = 0; i < 2; ++i) {
        err = ISDestroy(&superIS[i]);PYLITH_CHECK_ERROR(err);
    }
    err = DMDestroy(&superDM);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&stateVarIS);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&stateVarDM);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&locA);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&A);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&stateVarVec);PYLITH_CHECK_ERROR(err);
#endif

    PYLITH_METHOD_END;
#endif
} // _updateStateVars


// ---------------------------------------------------------------------------------------------------------------------
// Compute field derived from solution and auxiliary field.
void
pylith::feassemble::Integrator::_computeDerivedFields(const PylithReal t,
                                                      const PylithReal dt,
                                                      const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_computeDerivedFields(t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    PYLITH_METHOD_END;
} // _computeDerivedFields


// End of file
