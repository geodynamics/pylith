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

#include "pylith/topology/RefineMesh.hh" // ISA RefineMesh

class pylith::topology::RefineUniform : public pylith::topology::RefineMesh {
    friend class TestRefineUniform; // unit testing

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Constructor
    RefineUniform(void);

    /// Destructor
    ~RefineUniform(void);

    /// Deallocate data structures.
    void deallocate(void);

    /** Set number of levels of refinement.
     *
     * @param[in] numLevels Number of levels.
     */
    void setNumLevels(const size_t numLevels);

    /** Get number of levels of refinement.
     *
     * @returns Number of levels.
     */
    size_t getNumLevels(void) const;

    /** Refine mesh.
     *
     * @param mesh Mesh to refine.
     * @returns Mesh after refinement.
     */
    pylith::topology::Mesh* refine(const pylith::topology::Mesh& mesh);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    int _numLevels; ///< Number of levels of refinement

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    RefineUniform(const RefineUniform&); ///< Not implemented
    const RefineUniform& operator=(const RefineUniform&); ///< Not implemented

}; // RefineUniform

// End of file
