// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


#include <portinfo>

#include "pylith/petsc/SNESMonitor.hh" // implementation of class methods

<<<<<<< HEAD
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
    == == ==
    =
>>>>>>> dbd6ccd44 (cxx: Switch from uncrustify to clang-format for C/C++ formatting.)
#include "pylith/exceptions/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iomanip> // USES setw()

        // ------------------------------------------------------------------------------------------------
        // Constructor
    pylith::petsc::SNESMonitor::SNESMonitor(void) {
}

// ------------------------------------------------------------------------------------------------
// Destructor
pylith::petsc::SNESMonitor::~SNESMonitor(void) {
}

// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::petsc::SNESMonitor::deallocate(void) {
}

// ------------------------------------------------------------------------------------------------
// Register monitor with PETSc.
void
pylith::petsc::SNESMonitor::registerMonitor(PetscSNES snes) {
    PYLITH_METHOD_BEGIN;

    PylithCallPetsc(SNESMonitorSet(snes, PylithSNESMonitor, nullptr, nullptr));
    PylithCallPetsc(SNESConvergedReasonViewSet(snes, PylithSNESConvergedReason, nullptr, nullptr));

    PYLITH_METHOD_END;
}

// ------------------------------------------------------------------------------------------------
PetscErrorCode
pylith::petsc::SNESMonitor::PylithSNESMonitor(PetscSNES snes,
                                              PetscInt iteration,
                                              PetscReal residualNorm,
                                              PetscCtx context) {
    PYLITH_METHOD_BEGIN;
    PetscErrorCode err = PETSC_SUCCESS;

    PYLITH_INFO_ROOT(pylith::journal::solver,
                     "PETSc SNES: Iteration " << iteration << ", nondimensional residual norm=" << std::scientific
                                              << std::setprecision(6) << std::setw(12) << residualNorm);

    PYLITH_METHOD_RETURN(err);
}

// ------------------------------------------------------------------------------------------------
// PyLith SNES converged reason function.
PetscErrorCode
pylith::petsc::SNESMonitor::PylithSNESConvergedReason(PetscSNES snes,
                                                      PetscCtx context) {
    PYLITH_METHOD_BEGIN;
    PetscErrorCode err = PETSC_SUCCESS;

    const char* reason = nullptr;
    PylithCallPetsc(SNESGetConvergedReasonString(snes, &reason));

    PetscInt numIterations = 0;
    PylithCallPetsc(SNESGetIterationNumber(snes, &numIterations));
    PYLITH_INFO_ROOT(pylith::journal::solver,
                     "PETSc SNES nonlinear solver converged in " << numIterations << " iteration(s) due to " << reason);

    PYLITH_METHOD_RETURN(err);
}


// End of file
