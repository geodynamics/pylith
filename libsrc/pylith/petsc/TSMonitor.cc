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

#include "pylith/petsc/TSMonitor.hh" // implementation of class methods

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/exceptions/error.hh" // USES PYLITH_METHOD_*


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::petsc::TSMonitor::TSMonitor(void) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::petsc::TSMonitor::~TSMonitor(void) {}


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::petsc::TSMonitor::deallocate(void) {
}


// ------------------------------------------------------------------------------------------------
// Register monitor with PETSc.
void
pylith::petsc::TSMonitor::registerMonitor(PetscTS ts) {
    PYLITH_METHOD_BEGIN;

    PylithCallPetsc(TSMonitorSet(ts, PylithTSMonitor, nullptr, nullptr));

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
PetscErrorCode
pylith::petsc::TSMonitor::PylithTSMonitor(PetscTS ts,
                                          PetscInt i_step,
                                          PetscReal time,
                                          PetscVec u,
                                          PetscCtx context) {
    PYLITH_METHOD_BEGIN;
    PetscErrorCode err = PETSC_SUCCESS;

    PetscReal dt;
    PylithCallPetsc(TSGetTimeStep(ts, &dt));
    PYLITH_INFO_ROOT(pylith::journal::application_flow, "PETSc TS: Time step "<< i_step << ", advancing from nondimensional time=" << time << " by dt=" << dt);

    PYLITH_METHOD_RETURN(err);
}
