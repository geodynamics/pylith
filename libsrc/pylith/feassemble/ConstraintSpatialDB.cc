// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/feassemble/ConstraintSpatialDB.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintSpatialDB::ConstraintSpatialDB(pylith::problems::Physics* const physics) :
    Constraint(physics),
    _kernelConstraint(NULL) {
    GenericComponent::setName("constraintspatialdb");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::ConstraintSpatialDB::~ConstraintSpatialDB(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Set constraint kernel.
void
pylith::feassemble::ConstraintSpatialDB::setKernelConstraint(const PetscBdPointFunc kernel) {
    _kernelConstraint = kernel;
} // setkernelConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Initialize constraint domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::ConstraintSpatialDB::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" initialize(solution="<<solution.getLabel()<<")");

    Constraint::initialize(solution);

    const pylith::topology::Mesh& physicsDomainMesh = getPhysicsDomainMesh();
    delete _auxiliaryField;_auxiliaryField = _physics->createAuxiliaryField(solution, physicsDomainMesh);
    delete _diagnosticField;_diagnosticField = _physics->createDiagnosticField(solution, physicsDomainMesh);
    _computeDiagnosticField();
    const pylith::problems::Observer::NotificationType notification = pylith::problems::ObserverPhysics::DIAGNOSTIC;
    _observers->notifyObservers(0.0, 0, solution, notification);
    if (_observers) {
        const pylith::problems::Observer::NotificationType notification = pylith::problems::ObserverPhysics::DIAGNOSTIC;
        _observers->notifyObservers(0.0, 0, solution, notification);
    } // if
    delete _diagnosticField;_diagnosticField = NULL;

    delete _derivedField;_derivedField = _physics->createDerivedField(solution, physicsDomainMesh);

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.
    PetscDS prob = NULL;
    PetscDMLabel label = NULL;
    PetscDM dmSoln = solution.getDM();assert(dmSoln);
    PetscErrorCode err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);assert(prob);

    void* context = NULL;
    const PylithInt numConstrained = _constrainedDOF.size();
    const PetscInt i_field = solution.getSubfieldInfo(_subfieldName.c_str()).index;
    err = DMGetLabel(dmSoln, _labelName.c_str(), &label);PYLITH_CHECK_ERROR(err);
    err = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL_BD_FIELD, _labelName.c_str(), label, 1, &_labelValue, i_field,
                             numConstrained, &_constrainedDOF[0], (void (*)()) _kernelConstraint, NULL, context, NULL);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Set auxiliary field values for current time.
void
pylith::feassemble::ConstraintSpatialDB::setState(const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setState(t="<<t<<")");

    assert(_physics);
    _physics->updateAuxiliaryField(_auxiliaryField, t);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        assert(_auxiliaryField);
        debug << pythia::journal::at(__HERE__)
              << "Constraint component '" << GenericComponent::getName() << "' for '"
              <<_physics->getIdentifier()<<"': viewing auxiliary field." << pythia::journal::endl;
        _auxiliaryField->view("Constraint auxiliary field", pylith::topology::Field::VIEW_ALL);
    } // if

    PYLITH_METHOD_END;
} // setState


// ---------------------------------------------------------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::feassemble::ConstraintSpatialDB::setSolution(pylith::feassemble::IntegrationData* integrationData) {
    assert(integrationData);
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setSolution(integrationData="<<integrationData->str()<<")");

    assert(_auxiliaryField);
    assert(_physics);

    const pylith::topology::Field* solution = integrationData->getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const PylithReal t = integrationData->getScalar(pylith::feassemble::IntegrationData::time);

    PetscErrorCode err = 0;
    PetscDM dmSoln = solution->getDM();

    // Get label for constraint.
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmSoln, _labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);

    // Set auxiliary data
    const PetscInt part = 0;
    err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, part, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int fieldIndex = solution->getSubfieldInfo(_subfieldName.c_str()).index;
    const PylithInt numConstrained = _constrainedDOF.size();
    assert(solution->getLocalVector());
    err = DMPlexLabelAddFaceCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssentialBdField(dmSoln, t, solution->getLocalVector(), fieldIndex,
                                                     numConstrained, &_constrainedDOF[0], dmLabel, 1, &_labelValue,
                                                     _kernelConstraint, context, solution->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        debug << pythia::journal::at(__HERE__)
              << "Constraint component '" << GenericComponent::getName() << "' for '"
              <<_physics->getIdentifier()<<"''." << pythia::journal::endl;
        _auxiliaryField->view("Auxiliary field", pylith::topology::Field::VIEW_ALL);
        solution->view("Solution field after setting constrained values", pylith::topology::Field::VIEW_ALL);
    } // if

    PYLITH_METHOD_END;
} // setSolution


// End of file
