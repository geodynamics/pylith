// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/initializers/InitializePhase.hh" // ISA InitializePhase


class pylith::initializers::MeshReordering : public pylith::initializers::InitializePhase {
    friend class TestMeshReordering;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    MeshReordering(void);

    /// Default destructor
    ~MeshReordering(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Run initialization phase.
     *
     * @param[in] mesh Input mesh.
     * @param[in] problem Problem specification.
     * @returns Mesh after initialization phase.
     */
    pylith::topology::Mesh* run(pylith::topology::Mesh* mesh,
                                const pylith::problems::Problem& problem);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    MeshReordering(const MeshReordering&); ///< Not implemented
    const MeshReordering& operator=(const MeshReordering&); ///< Not implemented

}; // MeshReordering

// End of file
