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

class pylith::petsc::KSPMonitor : public pylith::utils::PyreComponent {
    friend class TestKSPMonitor; // unit testing

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    KSPMonitor(void);

    /// Destructor
    ~KSPMonitor(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Register monitor with PETSc.
     *
     * @param[inout] ts PETSc TS.
     */
    void registerMonitor(PetscKSP KSP);

    /** PyLith KSP monitor function.
     *
     * @param[in] KSP PETSc KSP
     * @param[in] iteration Iteration number.
     * @param[in] residualNorm Relative norma.
     * @param[inout] context User defined context.
     */
    static
    PetscErrorCode PylithKSPMonitor(PetscKSP KSP,
                                    PetscInt iteration,
                                    PetscReal residualNorm,
                                    PetscCtx context);

    /** PyLith KSP converged reason function.
     *
     * @param[in] KSP PETSc KSP
     * @param[inout] context User defined context.
     */
    static
    PetscErrorCode PylithKSPConvergedReason(PetscKSP KSP,
                                            PetscCtx context);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    KSPMonitor(const KSPMonitor&) = delete;

    const KSPMonitor& operator=(const KSPMonitor&) = delete;

}; // KSPMonitor

// End of file
