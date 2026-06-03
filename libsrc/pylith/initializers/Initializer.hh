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

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/problems/problemsfwd.hh"
#include "pylith/topology/topologyfwd.hh"

#include <vector> // HASA std::vector

class pylith::initializers::Initializer : public pylith::utils::PyreComponent {
    friend class TestInitializer;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    Initializer(void);

    /// Default destructor
    ~Initializer(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set phases.
     *
     * @param[in] phases Initialization phases.
     */
    void setPhases(pylith::initializers::InitializePhase* phases[],
                   const size_t numPhases);

    /** Run initialization phase.
     *
     * @param[in] problem Problem specification.
     * @returns Mesh after initialization phase.
     */
    pylith::topology::Mesh* runPhases(const pylith::problems::Problem& problem);


    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    std::vector<pylith::initializers::InitializePhase*> _phases; ///< Array of phases.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    Initializer(const Initializer&); ///< Not implemented
    const Initializer& operator=(const Initializer&); ///< Not implemented

}; // Initializer

// End of file
