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
#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

class pylith::petsc::TSMonitor : public pylith::utils::PyreComponent {
    friend class TestTSMonitor; // unit testing

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TSMonitor(void);

    /// Destructor
    ~TSMonitor(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Register monitor with PETSc.
     *
     * @param[inout] ts PETSc TS.
     */
    void registerMonitor(PetscTS ts);

    /** PyLith TS monitor function.
     *
     * @param[in] ts PETSc TS
     * @param[in] i_step Iteration number.
     * @param[in] time Current time.
     * @param[in] u Current iterate.
     * @param[inout] context User defined context.
     */
    static
    PetscErrorCode PylithTSMonitor(PetscTS ts,
                                   PetscInt i_step,
                                   PetscReal time,
                                   PetscVec u,
                                   PetscCtx context);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    TSMonitor(const TSMonitor&) = delete;

    const TSMonitor& operator=(const TSMonitor&) = delete;

}; // TSMonitor

// End of file
