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

#include "pylith/initializers/MeshReordering.hh" // implementation of class methods

#include "pylith/topology/ReverseCuthillMcKee.hh" // USES ReverseCuthillMcKee
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::initializers::MeshReordering::MeshReordering(void) {}


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::initializers::MeshReordering::~MeshReordering(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::initializers::MeshReordering::deallocate(void) {
    InitializePhase::deallocate();
} // deallocate


// ------------------------------------------------------------------------------------------------
// Run initialization phase.
pylith::topology::Mesh*
pylith::initializers::MeshReordering::run(pylith::topology::Mesh* mesh,
                                          const pylith::problems::Problem& problem) {
    PYLITH_METHOD_BEGIN;
    assert(mesh);

    pylith::topology::Mesh* meshNew = pylith::topology::ReverseCuthillMcKee::reorder(*mesh);

    PYLITH_METHOD_RETURN(meshNew);
} // run


// End of file
