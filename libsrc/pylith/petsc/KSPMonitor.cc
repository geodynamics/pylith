// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/petsc/KSPMonitor.hh" // implementation of class methods

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include <iomanip> // USES setw()

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::petsc::KSPMonitor::KSPMonitor(void) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::petsc::KSPMonitor::~KSPMonitor(void) {}


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::petsc::KSPMonitor::deallocate(void) {
}


// ------------------------------------------------------------------------------------------------
// Register monitor with PETSc.
void
pylith::petsc::KSPMonitor::registerMonitor(PetscKSP KSP) {
    PYLITH_METHOD_BEGIN;

    PylithCallPetsc(KSPMonitorSet(KSP, PylithKSPMonitor, nullptr, nullptr));
    PylithCallPetsc(KSPConvergedReasonViewSet(KSP, PylithKSPConvergedReason, nullptr, nullptr));

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
PetscErrorCode
pylith::petsc::KSPMonitor::PylithKSPMonitor(PetscKSP KSP,
                                            PetscInt iteration,
                                            PetscReal residualNorm,
                                            PetscCtx context) {
    PYLITH_METHOD_BEGIN;
    PetscErrorCode err = PETSC_SUCCESS;

    PYLITH_INFO_ROOT(pylith::journal::solver_detail3, "PETSc KSP: Iteration "<< iteration << ", nondimensional preconditioned residual norm=" << std::scientific << std::setprecision(6) << std::setw(12) << residualNorm);

    PYLITH_METHOD_RETURN(err);
}


// ------------------------------------------------------------------------------------------------
// PyLith KSP converged reason function.
PetscErrorCode
pylith::petsc::KSPMonitor::PylithKSPConvergedReason(PetscKSP KSP,
                                                    PetscCtx context) {
    PYLITH_METHOD_BEGIN;
    PetscErrorCode err = PETSC_SUCCESS;

    const char* reason = nullptr;
    PylithCallPetsc(KSPGetConvergedReasonString(KSP, &reason));

    PetscInt numIterations = 0;
    PylithCallPetsc(KSPGetIterationNumber(KSP, &numIterations));
    PYLITH_INFO_ROOT(pylith::journal::solver, "PETSc KSP linear solver converged in " << numIterations << " interation due to " << reason);

    PYLITH_METHOD_RETURN(err);
}
