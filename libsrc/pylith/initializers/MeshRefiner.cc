// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/initializers/MeshRefiner.hh" // implementation of class methods

#include "pylith/topology/RefineMesh.hh" // HOLDSA RefineMesh

#include "pylith/utils/error.hh" // USES PylithCallPetsc()
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::initializers::MeshRefiner::MeshRefiner(void) : _refiner(nullptr) {}


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::initializers::MeshRefiner::~MeshRefiner(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::initializers::MeshRefiner::deallocate(void) {
    _refiner = nullptr; // :TODO: Use shared pointer

    InitializePhase::deallocate();
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set refiner.
void
pylith::initializers::MeshRefiner::setRefiner(pylith::topology::RefineMesh* refiner) {
    assert(refiner);
    _refiner = refiner;
} // setRefiner


// ------------------------------------------------------------------------------------------------
// Run initialization phase.
pylith::topology::Mesh*
pylith::initializers::MeshRefiner::run(pylith::topology::Mesh* mesh,
                                       const pylith::problems::Problem& problem) {
    PYLITH_METHOD_BEGIN;
    assert(_refiner);
    assert(mesh);

    pylith::topology::Mesh* newMesh = _refiner->refine(*mesh);

    PYLITH_METHOD_RETURN(newMesh);
} // run


// End of file
