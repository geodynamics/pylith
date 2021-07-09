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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/feassemble/Constraint.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::Constraint::Constraint(pylith::problems::Physics* const physics) :
    PhysicsImplementation(physics),
    _constraintLabel(""),
    _subfieldName(""),
    _boundaryMesh(NULL) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::Constraint::~Constraint(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::Constraint::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete _boundaryMesh;_boundaryMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate


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
    PYLITH_JOURNAL_DEBUG("setMarkerLabel(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition constraint label.");
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
    PYLITH_JOURNAL_DEBUG("setSubfieldName(value="<<value<<")");

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
// Get mesh associated with constrained boundary.
const pylith::topology::Mesh&
pylith::feassemble::Constraint::getPhysicsDomainMesh(void) const {
    assert(_boundaryMesh);
    return *_boundaryMesh;
} // getPhysicsDomainMesh


// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::feassemble::Constraint::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("initialize(solution="<<solution.getLabel()<<")");

    delete _boundaryMesh;_boundaryMesh = pylith::topology::MeshOps::createLowerDimMesh(solution.getMesh(), _constraintLabel.c_str());assert(_boundaryMesh);
    PetscDM dmBoundary = _boundaryMesh->getDM();assert(dmBoundary);
    pylith::topology::CoordsVisitor::optimizeClosure(dmBoundary);

    assert(_physics);
    _observers = _physics->getObservers(); // Memory managed by Physics
    if (_observers) {
        _observers->setPhysicsImplementation(this);
        _observers->setTimeScale(_physics->getNormalizer().getTimeScale());
    } // if

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Update at end of time step.
void
pylith::feassemble::Constraint::poststep(const PylithReal t,
                                         const PylithInt tindex,
                                         const PylithReal dt,
                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("poststep(t="<<t<<", dt="<<dt<<")");

    notifyObservers(t, tindex, solution);

    PYLITH_METHOD_END;
} // poststep


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary field values to current time.
void
pylith::feassemble::Constraint::updateState(const PylithReal t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("updateState(t="<<t<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // updateState


// ---------------------------------------------------------------------------------------------------------------------
// Set solution time derivative
void
pylith::feassemble::Constraint::setSolutionDot(pylith::topology::Field* solutionDot,
                                               const PylithReal t) {
    PYLITH_METHOD_BEGIN;
    // Default is to do nothing.

    PYLITH_METHOD_END;
} // updateState


// End of file
