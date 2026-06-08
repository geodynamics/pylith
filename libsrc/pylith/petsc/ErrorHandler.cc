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

#include "pylith/petsc/ErrorHandler.hh" // implementation of class methods

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include <petscsys.h>

#include <iomanip> // USES setw()

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::petsc::ErrorHandler::ErrorHandler(void) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::petsc::ErrorHandler::~ErrorHandler(void) {}


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::petsc::ErrorHandler::deallocate(void) {
}


// ------------------------------------------------------------------------------------------------
// Register monitor with PETSc.
void
pylith::petsc::ErrorHandler::registerErrorHandler(void) {
    PYLITH_METHOD_BEGIN;

    PylithCallPetsc(PetscPushErrorHandler(PylithErrorHandler, nullptr));

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
PetscErrorCode
pylith::petsc::ErrorHandler::PylithErrorHandler(MPI_Comm comm,
                                                int sourceLine,
                                                const char *sourceFunction,
                                                const char *sourceFilename,
                                                PetscErrorCode errorCode,
                                                PetscErrorType errorType,
                                                const char *message,
                                                PetscCtx context) {
    PYLITH_METHOD_BEGIN;

    throw pylith::PetscError(message, sourceFunction, sourceFilename, sourceLine);

    PYLITH_METHOD_RETURN(PETSC_SUCCESS);
}
