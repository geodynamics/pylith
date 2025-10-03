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

#include "pylith/meshio/meshiofwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder::GroupPtType

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/scales/scalesfwd.hh" // forward declarations
#include "pylith/utils/arrayfwd.hh" // USES scalar_array, int_array, string_vector

// ------------------------------------------------------------------------------------------------
class pylith::meshio::MeshIO : public pylith::utils::PyreComponent {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    MeshIO(void);

    /// Destructor
    virtual ~MeshIO(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Read mesh from file.
     *
     * @param[in] mesh PyLith finite-element mesh.
     * @param[in] check Check topology of mesh.
     */
    void read(pylith::topology::Mesh* mesh,
              const bool check=true);

    /** Write mesh to file.
     *
     * @param mesh PyLith finite-element mesh.
     */
    void write(pylith::topology::Mesh* const mesh);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    /// Write mesh
    virtual
    void _write(void) const = 0;

    /// Read mesh
    virtual
    void _read(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _mesh; ///< Pointer to finite-element mesh.

}; // MeshIO

// End of file
