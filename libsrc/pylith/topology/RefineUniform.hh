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

#include "pylith/topology/topologyfwd.hh" // forward declarations

// RefineUniform --------------------------------------------------------
/// Object for managing uniform global mesh refinement.
class pylith::topology::RefineUniform { // RefineUniform
    friend class TestRefineUniform; // unit testing

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Constructor
    RefineUniform(void);

    /// Destructor
    ~RefineUniform(void);

    /// Deallocate data structures.
    void deallocate(void);

    /** Refine mesh.
     *
     * @param newMesh Refined mesh (result).
     * @param mesh Mesh to refine.
     * @param levels Number of levels to refine.
     */
    void refine(Mesh* const newMesh,
                const Mesh& mesh,
                const int levels=1);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    RefineUniform(const RefineUniform&); ///< Not implemented
    const RefineUniform& operator=(const RefineUniform&); ///< Not implemented

}; // RefineUniform

// End of file
