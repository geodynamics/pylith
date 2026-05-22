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

#include "pylith/initializers/MeshDistributor.hh" // implementation of class methods

#include "pylith/topology/Distributor.hh" // HOLDSA Distributor
#include "pylith/problems/Problem.hh" // USES Problem
#include "pylith/utils/error.hh" // USES PylithCallPetsc()
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::initializers::MeshDistributor::MeshDistributor(void) : _distributor(nullptr) {}


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::initializers::MeshDistributor::~MeshDistributor(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::initializers::MeshDistributor::deallocate(void) {
    _distributor = nullptr; // :TODO: Use shared pointer

    InitializePhase::deallocate();
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set distributor.
void
pylith::initializers::MeshDistributor::setDistributor(pylith::topology::Distributor* distributor) {
    assert(distributor);
    _distributor = distributor;
}


// ------------------------------------------------------------------------------------------------
// Run initialization phase.
pylith::topology::Mesh*
pylith::initializers::MeshDistributor::run(pylith::topology::Mesh* mesh,
                                           const pylith::problems::Problem& problem) {
    PYLITH_METHOD_BEGIN;
    assert(_distributor);
    assert(mesh);

    pylith::topology::Mesh* newMesh = _distributor->distribute(*mesh, problem.getInterfaces());

    PYLITH_METHOD_RETURN(newMesh);
} // run


// End of file
