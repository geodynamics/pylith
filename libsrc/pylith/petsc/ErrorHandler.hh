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

#include <mpi.h>
#include <petscerror.h>
#include <stdexcept> // USES std::exception_ptr

class pylith::petsc::ErrorHandler : public pylith::utils::PyreComponent {
    friend class TestErrorHandler; // unit testing

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Capture the first exception.
     *
     * Ignore overwrites from other frames.
     *
     * @param[in] exceptionPtr Exception pointer
     */
    static
    void set(std::exception_ptr exceptionPtr) noexcept;

    static
    std::exception_ptr release(void) noexcept;

    static
    bool exists(void) noexcept;

    static
    void clear(void) noexcept;

    /** PyLith error handling function.
     *
     * @param[in] comm MPI communicator
     * @param[in] sourceLine Line number in source file.
     * @param[in] sourceFunction Function in source file.
     * @param[in] sourceFilename Name of source file.
     * @param[in] errorCode Error code.
     * @param[in] errorType Error type.
     * @param[in] errorMessage Error message.
     * @param[inout] context User defined context.
     */
    static
    PetscErrorCode PylithErrorHandler(MPI_Comm comm,
                                      int sourceLine,
                                      const char *sourceFunction,
                                      const char *sourceFilename,
                                      PetscErrorCode errorCode,
                                      PetscErrorType errorType,
                                      const char *message,
                                      PetscCtx context);


    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    static thread_local std::exception_ptr _exception;

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    ErrorHandler(void) = delete;
    ErrorHandler(const ErrorHandler&) = delete;
    const ErrorHandler& operator=(const ErrorHandler&) = delete;

}; // ErrorHandler

// End of file
