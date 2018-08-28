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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> \
    // USES std::ostringstream

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintBoundary::ConstraintBoundary(pylith::problems::Physics* const physics) :
    Constraint(physics),
    _boundaryMesh(NULL) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::ConstraintBoundary::~ConstraintBoundary(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ConstraintBoundary::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    Constraint::deallocate();

    delete _boundaryMesh;_boundaryMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Get mesh associated with constrained domain.
const pylith::topology::Mesh&
pylith::feassemble::ConstraintBoundary::getConstraintDomainMesh(void) const {
    assert(_boundaryMesh);
    return *_boundaryMesh;
} // domainMesh


// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::feassemble::ConstraintBoundary::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("initialize(solution="<<solution.label()<<")");

    delete _boundaryMesh;_boundaryMesh = new pylith::topology::Mesh(solution.mesh(), _constraintLabel.c_str());assert(_boundaryMesh);
    PetscDM dmBoundary = _boundaryMesh->dmMesh();assert(dmBoundary);
    pylith::topology::CoordsVisitor::optimizeClosure(dmBoundary);

    Constraint::initialize(solution);

    PYLITH_METHOD_END;
} // initialize


// End of file
