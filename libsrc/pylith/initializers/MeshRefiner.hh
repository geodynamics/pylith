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

#include "pylith/topology/topologyfwd.hh" // HASA MeshRefiner

class pylith::initializers::MeshRefiner : public pylith::initializers::InitializePhase {
    friend class TestMeshRefiner;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    MeshRefiner(void);

    /// Default destructor
    ~MeshRefiner(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set refiner.
     *
     * @param[in] refiner Mesh refiner;
     */
    void setRefiner(pylith::topology::RefineMesh* refiner);

    /** Run initialization phase.
     *
     * @param[in] mesh Input mesh.
     * @param[in] problem Problem specification.
     * @returns Mesh after initialization phase.
     */
    pylith::topology::Mesh* run(pylith::topology::Mesh* mesh,
                                const pylith::problems::Problem& problem);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::topology::RefineMesh* _refiner; ///< Mesh refiner.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    MeshRefiner(const MeshRefiner&); ///< Not implemented
    const MeshRefiner& operator=(const MeshRefiner&); ///< Not implemented

}; // MeshRefiner

// End of file
