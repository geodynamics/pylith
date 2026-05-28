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

#include "pylith/initializers/initializersfwd.hh" // forward declarations

#include "pylith/problems/problemsfwd.hh"
#include "pylith/topology/topologyfwd.hh"

class pylith::initializers::InitializePhase {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    InitializePhase(void);

    /// Default destructor
    virtual ~InitializePhase(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Run initialization phase.
     *
     * @param[in] mesh Input mesh.
     * @param[in] problem Problem specification.
     * @returns Mesh after initialization phase.
     */
    virtual pylith::topology::Mesh* run(pylith::topology::Mesh* mesh,
                                        const pylith::problems::Problem& problem) = 0;

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    InitializePhase(const InitializePhase&); ///< Not implemented
    const InitializePhase& operator=(const InitializePhase&); ///< Not implemented

}; // InitializePhase

// End of file
