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

#include "ConstraintBoundary.hh" // implementation of object methods

#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> \
    // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::ConstraintBoundary::ConstraintBoundary(void) :
    _boundaryMesh(NULL),
    _label(""),
    _field("")
{ // constructor
    _description.label = "unknown";
    _description.vectorFieldType = pylith::topology::FieldBase::OTHER;
    _description.numComponents = 0;
    _description.scale = 1.0;
    _description.validator = NULL;
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::bc::ConstraintBoundary::~ConstraintBoundary(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::ConstraintBoundary::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    ConstraintPointwise::deallocate();

    delete _boundaryMesh; _boundaryMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set mesh label associated with boundary condition surface.
void
pylith::bc::ConstraintBoundary::label(const char* value) {
    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition label.");
    } // if

    _label = value;
} // label


// ----------------------------------------------------------------------
// Get mesh label associated with boundary condition surface.
const char*
pylith::bc::ConstraintBoundary::label(void) const {
    return _label.c_str();
} // label


// ----------------------------------------------------------------------
// Set name of field in solution to constrain.
void
pylith::bc::ConstraintBoundary::field(const char* value) {
    PYLITH_METHOD_BEGIN;

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for name of solution field for boundary condition.");
    } // if
    _field = value;

    PYLITH_METHOD_END;
}  // field


// ----------------------------------------------------------------------
// Get name of field in solution to constrain.
const char*
pylith::bc::ConstraintBoundary::field(void) const {
    journal::debug_t debug("boundarycondition");
    debug << journal::at(__HERE__)
          << "BoundaryCondition::field()" << journal::endl;

    return _field.c_str();
} // field


// ----------------------------------------------------------------------
// Get mesh associated with integrator domain.
const pylith::topology::Mesh&
pylith::bc::ConstraintBoundary::domainMesh(void) const {
    assert(_boundaryMesh);
    return *_boundaryMesh;
} // domainMesh


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::ConstraintBoundary::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    if (!solution.hasSubfield(_field.c_str())) {
        std::ostringstream msg;
        msg << "Cannot constrain field '"<< _field
            << "' in component '" << PyreComponent::identifier() << "'"
            << "; field is not in solution.";
        throw std::runtime_error(msg.str());
    } // if

    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    const int numComponents = info.description.numComponents;
    const int numConstrained = _constrainedDOF.size();
    for (int iConstrained = 0; iConstrained < numConstrained; ++iConstrained) {
        if (_constrainedDOF[iConstrained] >= numComponents) {
            std::ostringstream msg;
            msg << "Cannot constrain degree of freedom '" << _constrainedDOF[iConstrained] << "'"
                << " in component '" << PyreComponent::identifier() << "'"
                << "; solution field '" << _field << "' contains only " << numComponents << " components.";
            throw std::runtime_error(msg.str());
        } // if
    } // for

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::ConstraintBoundary::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label()<<")");

    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    _description = info.description;

    delete _boundaryMesh; _boundaryMesh = new pylith::topology::Mesh(solution.mesh(), _label.c_str()); assert(_boundaryMesh);
    PetscDM dmBoundary = _boundaryMesh->dmMesh(); assert(dmBoundary);
    pylith::topology::CoordsVisitor::optimizeClosure(dmBoundary);

    delete _auxField; _auxField = new pylith::topology::Field(*_boundaryMesh); assert(_auxField);
    _auxField->label("Dirichlet auxiliary");
    _auxFieldSetup(solution);
    _auxField->subfieldsSetup();
    _auxField->allocate();
    _auxField->zeroLocal();

    assert(_normalizer);
    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    factory->initializeSubfields();

    //_auxField->view("AUXILIARY FIELD"); // :DEBUG: TEMPORARY
    const bool infoOnly = true;
    notifyObservers(0.0, 0, solution, infoOnly);

    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    _setFEKernelConstraint(solution);
    PetscPointFunc bcKernel = _getFEKernelConstraint(); assert(bcKernel);

    void* context = NULL;
    const int labelId = 1;
    const PylithInt numConstrained = _constrainedDOF.size();
    err = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL_FIELD, label(), label(), info.index, numConstrained, &_constrainedDOF[0],
                             (void (*)())bcKernel, 1, &labelId, context); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::bc::ConstraintBoundary::setSolution(pylith::topology::Field* solution,
                                            const double t)
{ // setSolution
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setSolution(solution="<<solution->label()<<", t="<<t<<")");

    assert(solution);
    assert(_auxField);

    PetscErrorCode err = 0;
    PetscDM dmSoln = solution->dmMesh();
    PetscDM dmAux = _auxField->dmMesh();

    // Get label for constraint.
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmSoln, _label.c_str(), &dmLabel); PYLITH_CHECK_ERROR(err);

    // Set auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxField->localVector()); PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int labelId = 1;
    const int fieldIndex = solution->subfieldInfo(_field.c_str()).index;
    PetscPointFunc bcKernel = _getFEKernelConstraint(); assert(bcKernel);
    const PylithInt numConstrained = _constrainedDOF.size();
    assert(solution->localVector());
    err = DMPlexLabelAddCells(dmSoln, dmLabel); PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssentialField(dmSoln, t,
                                                   solution->localVector(), fieldIndex, numConstrained, &_constrainedDOF[0], dmLabel, 1, &labelId, bcKernel,
                                                   context, solution->localVector()); PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmSoln, dmLabel); PYLITH_CHECK_ERROR(err);

    //solution->view("SOLUTION at end of setSolution()"); // :DEBUG: TEMPORARY

    PYLITH_METHOD_END;
} // setSolution


// End of file
