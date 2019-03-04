// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/feassemble/Constraint.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::Constraint::Constraint(pylith::problems::Physics* const physics) :
    PhysicsImplementation(physics),
    _constraintLabel(""),
    _subfieldName(""),
    _kernelConstraint(NULL)
{}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::Constraint::~Constraint(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Set indices of constrained degrees of freedom at each location.
void
pylith::feassemble::Constraint::setConstrainedDOF(const int* flags,
                                                  const int size) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setConstrainedDOF(flags="<<flags<<", size="<<size<<")");

    assert((size > 0 && flags) || (!size && !flags));

    _constrainedDOF.resize(size);
    for (int i = 0; i < size; ++i) {
        if (flags[i] < 0) {
            std::ostringstream msg;
            assert(_physics);
            msg << "Constrained DOF '" << flags[i] << "' must be nonnegative in constraint component '" << _physics->getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        _constrainedDOF[i] = flags[i];
    } // for

    PYLITH_METHOD_END;
} // setConstrainedDOF


// ---------------------------------------------------------------------------------------------------------------------
// Get indices of constrained degrees of freedom.
const pylith::int_array&
pylith::feassemble::Constraint::getConstrainedDOF(void) const {
    return _constrainedDOF;
} // getConstrainedDOF


// ---------------------------------------------------------------------------------------------------------------------
// Set label marking constrained degrees of freedom.
void
pylith::feassemble::Constraint::setMarkerLabel(const char* value) {
    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition integrator label.");
    } // if

    _constraintLabel = value;
} // setMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get label marking constrained degrees of freedom.
const char*
pylith::feassemble::Constraint::getMarkerLabel(void) const {
    return _constraintLabel.c_str();
} // getMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Set name of constrained solution subfield.
void
pylith::feassemble::Constraint::setSubfieldName(const char* value) {
    if (!value || (0 == strlen(value))) {
        std::ostringstream msg;
        assert(_physics);
        msg << "Empty string given for name of solution subfield for constraint '" << _physics->getIdentifier()
            <<"'.";
        throw std::runtime_error(msg.str());
    } // if
    _subfieldName = value;
} // setSubfieldName


// ---------------------------------------------------------------------------------------------------------------------
// Get name of constrained solution subfield.
const char*
pylith::feassemble::Constraint::getSubfieldName(void) const {
    return _subfieldName.c_str();
} // getSubfieldName


// ---------------------------------------------------------------------------------------------------------------------
// Set constraint kernel.
void
pylith::feassemble::Constraint::setKernelConstraint(const PetscPointFunc kernel) {
    _kernelConstraint = kernel;
} // setKernelConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Initialize constraint domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::Constraint::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("intialize(solution="<<solution.label()<<")");

    assert(_physics);

    const pylith::topology::Mesh& physicsDomainMesh = getPhysicsDomainMesh();

    delete _auxiliaryField;_auxiliaryField = _physics->createAuxiliaryField(solution, physicsDomainMesh);
    delete _derivedField;_derivedField = _physics->createDerivedField(solution, physicsDomainMesh);
    _observers = _physics->getObservers();assert(_observers); // Memory managed by Python
    _observers->setPhysicsImplementation(this);

    const bool infoOnly = true;
    _observers->notifyObservers(0.0, 0, solution, infoOnly);

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.
    PetscDS prob = NULL;
    PetscDM dmSoln = solution.dmMesh();assert(dmSoln);
    PetscErrorCode err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);assert(prob);

    void* context = NULL;
    const int labelId = 1;
    const PylithInt numConstrained = _constrainedDOF.size();
    const PetscInt i_field = solution.subfieldInfo(_subfieldName.c_str()).index;
    err = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL_FIELD, _constraintLabel.c_str(), _constraintLabel.c_str(), i_field,
                             numConstrained, &_constrainedDOF[0], (void (*)())_kernelConstraint, 1, &labelId, context);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Update at beginning of time step.
void
pylith::feassemble::Constraint::prestep(const double t,
                                        const double dt) { // prestep
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("prestep(t="<<t<<", dt="<<dt<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // prestep


// ---------------------------------------------------------------------------------------------------------------------
// Update at end of time step.
void
pylith::feassemble::Constraint::poststep(const PylithReal t,
                                         const PylithInt tindex,
                                         const PylithReal dt,
                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("poststep(t="<<t<<", dt="<<dt<<") empty method");

    const bool infoOnly = false;
    assert(_observers);
    _observers->notifyObservers(t, tindex, solution, infoOnly);

    PYLITH_METHOD_END;
} // poststep


// ---------------------------------------------------------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::feassemble::Constraint::setSolution(pylith::topology::Field* solution,
                                            const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setSolution(solution="<<solution->label()<<", t="<<t<<")");

    assert(solution);
    assert(_auxiliaryField);
    assert(_physics);

    PetscErrorCode err = 0;
    PetscDM dmSoln = solution->dmMesh();
    PetscDM dmAux = _auxiliaryField->dmMesh();

    // Get label for constraint.
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmSoln, _constraintLabel.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);

    // Set auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int labelId = 1;
    const int fieldIndex = solution->subfieldInfo(_subfieldName.c_str()).index;
    const PylithInt numConstrained = _constrainedDOF.size();
    assert(solution->localVector());
    err = DMPlexLabelAddCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssentialField(dmSoln, t, solution->localVector(), fieldIndex,
                                                   numConstrained, &_constrainedDOF[0], dmLabel, 1, &labelId,
                                                   _kernelConstraint, context, solution->localVector());PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);

    // solution->view("SOLUTION at end of setSolution()"); // :DEBUG: TEMPORARY

    PYLITH_METHOD_END;
} // setSolution


// ---------------------------------------------------------------------------------------------------------------------
// Set constants used in finite-element kernels (point-wise functions).
void
pylith::feassemble::Constraint::_setKernelConstants(const pylith::topology::Field& solution,
                                                    const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setKernelConstants(solution="<<solution.label()<<", dt="<<dt<<")");

    assert(_physics);
    const pylith::real_array& constants = _physics->getKernelConstants(dt);

    PetscDS prob = NULL;
    PetscDM dmSoln = solution.dmMesh();assert(dmSoln);

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.
    PetscErrorCode err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);assert(prob);
    if (constants.size() > 0) {
        err = PetscDSSetConstants(prob, constants.size(), const_cast<double*>(&constants[0]));PYLITH_CHECK_ERROR(err);
    } else {
        err = PetscDSSetConstants(prob, 0, NULL);PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // _setKernelConstants


// End of file
