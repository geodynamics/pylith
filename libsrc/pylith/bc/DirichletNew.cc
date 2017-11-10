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

#include "pylith/bc/DirichletNew.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream>

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletNew::DirichletNew(void) :
    _boundaryMesh(NULL),
    _bcKernel(NULL)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletNew::~DirichletNew(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletNew::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    ConstraintPointwise::deallocate();

    delete _boundaryMesh; _boundaryMesh = NULL;
    _bcKernel = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::DirichletNew::verifyConfiguration(const pylith::topology::Field& solution) const {
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
pylith::bc::DirichletNew::initialize(const pylith::topology::Field& solution)
{ // initialize
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label()<<")");

    _setFEKernelsConstraint(solution);

    _boundaryMesh = new pylith::topology::Mesh(solution.mesh(), _label.c_str()); assert(_boundaryMesh);
    PetscDM dmBoundary = _boundaryMesh->dmMesh(); assert(dmBoundary);
    pylith::topology::CoordsVisitor::optimizeClosure(dmBoundary);

    delete _auxField; _auxField = new pylith::topology::Field(*_boundaryMesh); assert(_auxField);
    _auxField->label("auxiliary subfields");
    _auxFieldSetup(solution);
    _auxField->subfieldsSetup();
    _auxField->allocate();
    _auxField->zeroLocal();

    assert(_normalizer);
    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    factory->initializeSubfields();

    //_auxField->view("AUXILIARY FIELD"); // :DEBUG: TEMPORARY

    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int labelId = 1;
    const PylithInt numConstrained = _constrainedDOF.size();
    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    err = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL_FIELD, label(), label(), info.index, numConstrained, &_constrainedDOF[0],
                             (void (*)())_bcKernel, 1, &labelId, context); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::bc::DirichletNew::setSolution(pylith::topology::Field* solution,
                                      const double t)
{ // setSolution
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setSolution(solution="<<solution->label()<<", t="<<t<<")");

    assert(solution);
    assert(_auxField);

    PetscErrorCode err;
    PetscDM dmSoln = solution->dmMesh();
    PetscDM dmAux = _auxField->dmMesh();

    // Get label for constraint.
    PetscDMLabel dmLabel;
    err = DMGetLabel(dmSoln, _label.c_str(), &dmLabel); PYLITH_CHECK_ERROR(err);

    // Set auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxField->localVector()); PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int labelId = 1;
    const int fieldIndex = solution->subfieldInfo(_field.c_str()).index;
    const PylithInt numConstrained = _constrainedDOF.size();
    assert(solution->localVector());
    err = DMPlexLabelAddCells(dmSoln, dmLabel); PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssentialField(dmSoln, t,
                                                   solution->localVector(), fieldIndex, numConstrained, &_constrainedDOF[0], dmLabel, 1, &labelId, _bcKernel,
                                                   context, solution->localVector());
    PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmSoln, dmLabel); PYLITH_CHECK_ERROR(err);

    //solution->view("SOLUTION"); // :DEBUG: TEMPORARY

    PYLITH_METHOD_END;
} // setSolution


// End of file
