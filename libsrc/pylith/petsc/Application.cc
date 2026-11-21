// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


#include <portinfo>

#include "pylith/petsc/Application.hh" // implementation of class methods

#include "pylith/exceptions/error.hh" // USES PylithCallPetsc
#include "pylith/petsc/ErrorHandler.hh" // USES ErrorHandler

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


// End of file
