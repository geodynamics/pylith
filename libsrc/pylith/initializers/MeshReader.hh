// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/initializers/InitializePhase.hh" // ISA InitializePhase

#include "pylith/meshio/meshiofwd.hh" // HOLDSA MeshIO

class pylith::initializers::MeshReader : public pylith::initializers::InitializePhase {
    friend class TestMeshReader;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    MeshReader(void);

    /// Default destructor
    ~MeshReader(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set reader.
     *
     * @param[in] reader Mesh reader;
     */
    void setReader(pylith::meshio::MeshIO* reader);

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

    pylith::meshio::MeshIO* _reader; ///< Mesh reader.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    MeshReader(const MeshReader&); ///< Not implemented
    const MeshReader& operator=(const MeshReader&); ///< Not implemented

}; // MeshReader

// End of file
