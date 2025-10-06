// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/utils/error.hh"

#include "pylith/utils/TSAdaptImpulse.hh" // implementation of class methods

#include "petscts.h"
#include "petsc/private/tsimpl.h"

// ------------------------------------------------------------------------------------------------
double pylith::utils::TSAdaptImpulse::_impulseStep = 1.0e-6;

// ------------------------------------------------------------------------------------------------
// Set time step for impulse portion.
void
pylith::utils::TSAdaptImpulse::setImpulseTimeStep(const double timestep) {
    _impulseStep = timestep;
}


// ------------------------------------------------------------------------------------------------
// Get time step for impulse portion.
double
pylith::utils::TSAdaptImpulse::getImpulseTimeStep(void) {
    return _impulseStep;
}


// ------------------------------------------------------------------------------------------------
// Set PETSc TS adapt.
void
pylith::utils::TSAdaptImpulse::set(PetscTS ts) {
    TSAdapt adapt;

    PetscErrorCode err = TSGetAdapt(ts, &adapt);PYLITH_CHECK_ERROR(err);
    adapt->ops->choose = TSAdaptChoose;
}


// ------------------------------------------------------------------------------------------------
// PETSc adaptive time stepper.
PetscErrorCode
pylith::utils::TSAdaptImpulse::TSAdaptChoose(TSAdapt adapt,
                                             TS ts,
                                             PetscReal h,
                                             PetscInt *next_sc,
                                             PetscReal *next_h,
                                             PetscBool *accept,
                                             PetscReal *wlte,
                                             PetscReal *wltea,
                                             PetscReal *wlter) {
    static PetscReal dtTarget = -1.0;
    PetscReal dtInitial = TSAdaptImpulse::getImpulseTimeStep();
    PetscInt step;

    PetscFunctionBeginUser;
    PetscCall(TSGetStepNumber(ts, &step));
    if (!step) {
        if (PetscAbsReal(dtInitial - h) > PETSC_SMALL) {
            *accept = PETSC_FALSE;
            *next_h = dtInitial;
            dtTarget = h;
        } else {
            *accept = PETSC_TRUE;
            *next_h = dtTarget < 0.0 ? dtInitial : dtTarget;
            dtTarget = -1.0;
        }
    } else {
        *accept = PETSC_TRUE;
        *next_h = h;
    }
    *next_sc = 0; /* Reuse the same order scheme */
    *wlte = -1; /* Weighted local truncation error was not evaluated */
    *wltea = -1; /* Weighted absolute local truncation error was not evaluated */
    *wlter = -1; /* Weighted relative local truncation error was not evaluated */
    PetscFunctionReturn(PETSC_SUCCESS);
}
