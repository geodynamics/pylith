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

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

class pylith::topology::RefineMesh : public pylith::utils::PyreComponent {
    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Constructor
    RefineMesh(void);

    /// Destructor
    virtual ~RefineMesh(void);

    /// Deallocate data structures.
    void deallocate(void);

    /** Refine mesh.
     *
     * @param mesh Mesh to refine.
     * @returns Mesh after refinement.
     */
    virtual
    pylith::topology::Mesh* refine(const pylith::topology::Mesh& mesh) = 0;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    RefineMesh(const RefineMesh&); ///< Not implemented
    const RefineMesh& operator=(const RefineMesh&); ///< Not implemented

}; // RefineMesh

// End of file
