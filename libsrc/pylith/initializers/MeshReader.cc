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

#include "pylith/initializers/MeshReader.hh" // implementation of class methods

#include "pylith/meshio/MeshIO.hh" // HOLDSA MeshIO
#include "pylith/problems/Problem.hh" // USES Problem
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/utils/error.hh" // USES PylithCallPetsc()
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::initializers::MeshReader::MeshReader(void) : _reader(nullptr) {}


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::initializers::MeshReader::~MeshReader(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::initializers::MeshReader::deallocate(void) {
    _reader = nullptr; // :TODO: Use shared pointer

    InitializePhase::deallocate();
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set reader.
void
pylith::initializers::MeshReader::setReader(pylith::meshio::MeshIO* reader) {
    assert(reader);
    _reader = reader;
}


// ------------------------------------------------------------------------------------------------
// Run initialization phase.
pylith::topology::Mesh*
pylith::initializers::MeshReader::run(pylith::topology::Mesh* mesh,
                                      const pylith::problems::Problem& problem) {
    PYLITH_METHOD_BEGIN;
    assert(_reader);

    pylith::topology::Mesh* meshNew = new pylith::topology::Mesh();
    _reader->read(meshNew);
    pylith::topology::MeshOps::nondimensionalize(meshNew, problem.getScales());

    PYLITH_METHOD_RETURN(meshNew);
} // run


// End of file
