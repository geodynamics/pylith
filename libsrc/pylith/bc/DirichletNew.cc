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

#include "DirichletNew.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
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
    _description.label = "unknown";
    _description.vectorFieldType = pylith::topology::Field::OTHER;
    _description.numComponents = 0;
    _description.scale = 1.0;
    _description.validator = NULL;
    const pylith::topology::FieldBase::Discretization defaultInfo = {-1, -1, true, pylith::topology::FieldBase::POLYNOMIAL_SPACE};
    _auxFieldsFEInfo["default"] = defaultInfo;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletNew::~DirichletNew(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletNew::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    ConstraintPointwise::deallocate();

    delete _boundaryMesh; _boundaryMesh = NULL;
    _bcKernel = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletNew::initialize(const pylith::topology::Field& solution)
{ // initialize
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label()<<")");

    const topology::Field::SubfieldInfo& info = solution.subfieldInfo(_field.c_str());
    _description = info.description;

    _setFEKernelsConstraint(solution);

    _boundaryMesh = new pylith::topology::Mesh(solution.mesh(), _label.c_str()); assert(_boundaryMesh);
    PetscDM dmBoundary = _boundaryMesh->dmMesh(); assert(dmBoundary);
    pylith::topology::CoordsVisitor::optimizeClosure(dmBoundary);

    delete _auxFields; _auxFields = new pylith::topology::Field(*_boundaryMesh); assert(_auxFields);
    delete _auxFieldsQuery; _auxFieldsQuery = new pylith::topology::FieldQuery(*_auxFields); assert(_auxFieldsQuery);
    _auxFields->label("auxiliary fields");
    _auxFieldsSetup(solution);
    _auxFields->subfieldsSetup();
    _auxFields->allocate();
    _auxFields->zeroLocal();

    if (_auxFieldsDB) {
        assert(_normalizer);
        _auxFieldsQuery->openDB(_auxFieldsDB, _normalizer->lengthScale());
        _auxFieldsQuery->queryDB();
        _auxFieldsQuery->closeDB(_auxFieldsDB);
    } else { // else
        PYLITH_COMPONENT_ERROR("Unknown case for setting up auxiliary fields.");
        throw std::logic_error("Unknown case for setting up auxiliary fields.");
    } // if/else
    _auxFields->view("AUXILIARY FIELDS"); // :DEBUGGING: TEMPORARY

    // Get label for constraint.
    PetscDMLabel dmLabel;
    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscErrorCode err = DMGetLabel(dmSoln, _label.c_str(), &dmLabel); PYLITH_CHECK_ERROR(err);

    void* context = NULL;
    const int labelId = 1;
    const PylithInt numConstrained = _constrainedDOF.size();
    err = DMAddBoundary(dmSoln, DM_BC_ESSENTIAL_FIELD, label(), label(), info.index, numConstrained, &_constrainedDOF[0],
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
    assert(_auxFields);

    PetscErrorCode err;
    PetscDM dmSoln = solution->dmMesh();
    PetscDM dmAux = _auxFields->dmMesh();

    // Get label for constraint.
    PetscDMLabel dmLabel;
    err = DMGetLabel(dmSoln, _label.c_str(), &dmLabel); PYLITH_CHECK_ERROR(err);

    // Set auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxFields->localVector()); PYLITH_CHECK_ERROR(err);

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

    solution->view("SOLUTION"); // :DEBUGGING: TEMPORARY

    PYLITH_METHOD_END;
} // setSolution


// End of file
