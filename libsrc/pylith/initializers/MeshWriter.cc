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

#include "pylith/initializers/MeshWriter.hh" // implementation of class methods

#include "pylith/meshio/MeshIO.hh" // HOLDSA MeshIO
#include "pylith/utils/error.hh" // USES PylithCallPetsc()
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::initializers::MeshWriter::MeshWriter(void) : _writer(nullptr) {}


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::initializers::MeshWriter::~MeshWriter(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::initializers::MeshWriter::deallocate(void) {
    _writer = nullptr; // :TODO: Use shared pointer

    InitializePhase::deallocate();
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set writer.
void
pylith::initializers::MeshWriter::setWriter(pylith::meshio::MeshIO* writer) {
    assert(writer);
    _writer = writer;
}


// ------------------------------------------------------------------------------------------------
// Run initialization phase.
pylith::topology::Mesh*
pylith::initializers::MeshWriter::run(pylith::topology::Mesh* mesh,
                                      const pylith::problems::Problem& problem) {
    PYLITH_METHOD_BEGIN;
    assert(_writer);
    assert(mesh);

    _writer->write(mesh);

    PetscDM dmOrig = mesh->getDM();assert(dmOrig);
    PetscObjectReference(PetscObject(dmOrig));
    pylith::topology::Mesh* meshNew = new pylith::topology::Mesh(dmOrig, *mesh);assert(meshNew);

    PYLITH_METHOD_RETURN(meshNew);
} // run


// End of file
