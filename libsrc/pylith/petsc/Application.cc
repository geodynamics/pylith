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

#include "pylith/petsc/Application.hh" // implementation of class methods

#include "pylith/petsc/ErrorHandler.hh" // USES ErrorHandler
#include "pylith/utils/error.hh" // USES PylithCallPetsc

#include <petscsys.h>

// ------------------------------------------------------------------------------------------------
// Initialize
void
pylith::petsc::Application::initialize(int argc,
                                       char* argv[]) {
    PylithCallPetsc(PetscInitialize(&argc, &argv, nullptr, nullptr));
    PylithCallPetsc(PetscPushErrorHandler(ErrorHandler::PylithErrorHandler, nullptr));
}


// ------------------------------------------------------------------------------------------------
// Finalize.
void
pylith::petsc::Application::finalize(void) {
    PetscPopErrorHandler();
    PylithCallPetsc(PetscFinalize());
}


// ------------------------------------------------------------------------------------------------
void
pylith::petsc::Application::registerCitation(const char* citation) {
    PylithCallPetsc(PetscOptionsSetValue(NULL, "-citations", ""));
    PetscBool set = PetscBool(PETSC_FALSE);
    PylithCallPetsc(PetscCitationsRegister(citation, &set));

}
