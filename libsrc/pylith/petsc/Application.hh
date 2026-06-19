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

#include "petscfwd.hh" // forward declarations

#include "pylith/petsc/petsc_types.h" // USES PetscErrorCode

#include <mpi.h>
#include <petscerror.h>

class pylith::petsc::Application {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Initialize
     *
     * @param[in] argc Number of command line arguments.
     * @param[in] argv Command line argument.
     */
    static
    void initialize(int argc,
                    char* argv[]);

    /// Finalize.
    static
    void finalize(void);

    /** Register citations.
     *
     * @param[in] citation BiBTeX citation string.
     */
    static
    void registerCitation(const char* citation);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    Application(void) = delete;
    Application(const Application&) = delete;
    const Application& operator=(const Application&) = delete;

}; // Application

// End of file
