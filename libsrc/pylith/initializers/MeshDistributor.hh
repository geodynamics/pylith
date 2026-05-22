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

#include "pylith/topology/topologyfwd.hh" // HOLDSA Distributor

class pylith::initializers::MeshDistributor : public pylith::initializers::InitializePhase {
    friend class TestMeshDistributor;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    MeshDistributor(void);

    /// Default destructor
    ~MeshDistributor(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set distributor.
     *
     * @param[in] distributor Mesh distributor;
     */
    void setDistributor(pylith::topology::Distributor* distributor);

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

    pylith::topology::Distributor* _distributor; ///< Mesh distributor.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    MeshDistributor(const MeshDistributor&); ///< Not implemented
    const MeshDistributor& operator=(const MeshDistributor&); ///< Not implemented

}; // MeshDistributor


// End of file
