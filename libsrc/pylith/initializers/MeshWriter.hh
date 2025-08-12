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

#include "pylith/meshio/meshiofwd.hh" // HOLDSA MeshIO

class pylith::initializers::MeshWriter : public pylith::initializers::InitializePhase {
    friend class TestMeshWriter;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    MeshWriter(void);

    /// Default destructor
    ~MeshWriter(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set writer.
     *
     * @param[in] writer Mesh writer;
     */
    void setWriter(pylith::meshio::MeshIO* writer);

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

    pylith::meshio::MeshIO* _writer; ///< Mesh writer.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    MeshWriter(const MeshWriter&); ///< Not implemented
    const MeshWriter& operator=(const MeshWriter&); ///< Not implemented

}; // MeshWriter

// End of file
