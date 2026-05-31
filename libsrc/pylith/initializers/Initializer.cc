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

#include "pylith/initializers/Initializer.hh" // implementation of class methods

#include "pylith/initializers/InitializePhase.hh" // HASA InitializePhase
#include "pylith/topology/Mesh.hh" // HASA Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/utils/error.hh" // USES PylithCallPetsc()
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*


// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::initializers::Initializer::Initializer(void) {}


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::initializers::Initializer::~Initializer(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::initializers::Initializer::deallocate(void) {
    _phases.resize(0); // :TODO: Use shared pointer
}


// ------------------------------------------------------------------------------------------------
// Set phases.
void
pylith::initializers::Initializer::setPhases(pylith::initializers::InitializePhase* phases[],
                                             const size_t numPhases) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "setPhases(numPhases="<<numPhases<<")");

    assert( (!phases && 0 == numPhases) || (phases && 0 < numPhases) );

    _phases.resize(numPhases);
    for (size_t i = 0; i < numPhases; ++i) {
        _phases[i] = phases[i];
    } // for

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Run initialization phase.
pylith::topology::Mesh*
pylith::initializers::Initializer::runPhases(const pylith::problems::Problem& problem) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG(pylith::journal::application_flow, "Running mesh initialization.");

    pylith::topology::Mesh* meshNew = nullptr;
    pylith::topology::Mesh* meshPhase = nullptr;
    const size_t numPhases = _phases.size();
    for (size_t i = 0; i < numPhases; ++i) {
        assert(_phases[i]);
        meshNew = _phases[i]->run(meshPhase, problem);
        pylith::topology::MeshOps::checkTopology(*meshNew); // TEMPORARY
        delete meshPhase;meshPhase = meshNew;
    } // for

    pythia::journal::debug_t debug(pylith::journal::mesh);
    if (debug.state()) {
        meshNew->view(":mesh_domain_after_initialize.txt:ascii_info_detail");
        pylith::topology::Mesh* meshExploded = pylith::topology::MeshOps::explode(*meshNew);
        meshExploded->view(":mesh_domain_after_initialize.tex:ascii_latex");
        delete meshExploded;meshExploded = nullptr;
    } // if
    pylith::topology::MeshOps::checkTopology(*meshNew);

    PYLITH_METHOD_RETURN(meshNew);
} // runPhases


// End of file
