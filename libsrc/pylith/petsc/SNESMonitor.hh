// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


#pragma once

#include "petscfwd.hh" // forward declarations

#include "pylith/petsc/petsc_types.h" // USES PetscErrorCode
#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent


class pylith::petsc::SNESMonitor : public pylith::utils::PyreComponent {
    friend class TestSNESMonitor; // unit testing

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
  public:
    /// Constructor
    SNESMonitor(void);

    /// Destructor
    ~SNESMonitor(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Register monitor with PETSc.
     *
     * @param[inout] ts PETSc TS.
     */
    void registerMonitor(PetscSNES snes);

    /** PyLith SNES monitor function.
     *
     * @param[in] snes PETSc SNES
     * @param[in] iteration Iteration number.
     * @param[in] residualNorm Relative norma.
     * @param[inout] context User defined context.
     */
    static PetscErrorCode PylithSNESMonitor(PetscSNES snes,
                                            PetscInt iteration,
                                            PetscReal residualNorm,
                                            PetscCtx context);

    /** PyLith SNES converged reason function.
     *
     * @param[in] snes PETSc SNES
     * @param[inout] context User defined context.
     */
    static PetscErrorCode PylithSNESConvergedReason(PetscSNES snes,
                                                    PetscCtx context);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
  private:
    SNESMonitor(const SNESMonitor&) = delete;

    const SNESMonitor& operator=(const SNESMonitor&) = delete;

}; // SNESMonitor


// End of file
