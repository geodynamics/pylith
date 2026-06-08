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

class pylith::petsc::ErrorHandler : public pylith::utils::PyreComponent {
    friend class TestErrorHandler; // unit testing

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    ErrorHandler(void);

    /// Destructor
    ~ErrorHandler(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Register monitor with PETSc.
     *
     * @param[inout] ts PETSc TS.
     */
    void registerErrorHandler(void);

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


    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    ErrorHandler(const ErrorHandler&) = delete;

    const ErrorHandler& operator=(const ErrorHandler&) = delete;

}; // ErrorHandler

// End of file
