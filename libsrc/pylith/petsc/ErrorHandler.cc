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
#include "pylith/exceptions/error.hh" // USES PYLITH_METHOD_*
#include "pylith/exceptions/Exceptions.hh" // USES PetscError

#include <petscsys.h>

thread_local std::exception_ptr pylith::petsc::ErrorHandler::_exception;

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

    try {
        // Only capture the first (deepest) error; ExceptionSlot::set
        // already guards against overwriting, but constructing the
        // exception object isn't free, so skip the work if we can.
        if (!exists()) {
            auto exPtr = std::make_exception_ptr(pylith::exceptions::PetscError(message, sourceFilename, sourceLine, sourceFunction));
            set(exPtr);
        }
    } catch (...) {
        // Constructing PetscError or the exception_ptr failed
        // (almost certainly bad_alloc). Fall back to a pre-built
        // sentinel so the C++ side still sees *something*.
        try {
            set(std::make_exception_ptr(std::bad_alloc{}));
        } catch (...) {
            // Truly out of options. Swallow and let the error code
            // propagate.
        } // try/catcg
    } // try/catch

    PYLITH_METHOD_RETURN(errorCode);
}


// ------------------------------------------------------------------------------------------------
void
pylith::petsc::ErrorHandler::set(std::exception_ptr exceptionPtr) noexcept {
    if (!_exception) {
        _exception = std::move(exceptionPtr);
    } // if
}


// ------------------------------------------------------------------------------------------------
std::exception_ptr
pylith::petsc::ErrorHandler::release(void) noexcept {
    auto local = std::move(_exception);
    _exception = nullptr;

    return local;
}


// ------------------------------------------------------------------------------------------------
bool
pylith::petsc::ErrorHandler::exists(void) noexcept {
    return static_cast<bool>(_exception);
}


// ------------------------------------------------------------------------------------------------
void
pylith::petsc::ErrorHandler::clear(void) noexcept {
    _exception = nullptr;
}
