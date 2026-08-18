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

#include "Exceptions.hh"
#include "pylith/petsc/ErrorHandler.hh"

#include "petsc.h"

#include <cassert>
#include <sstream>

#undef __FUNCT__
#if defined(__FUNCTION_NAME__)
#define __FUNCT__ __FUNCTION_NAME__
#else
#define __FUNCT__ __func__
#endif

#define PYLITH_METHOD_BEGIN PetscFunctionBeginUser
#define PYLITH_METHOD_END PetscFunctionReturnVoid()
#define PYLITH_METHOD_RETURN(v) PetscFunctionReturn(v)

#define PYLITH_ERROR_RETURN(comm,error,msg) SETERRQ(comm,error,"%s",msg)

#define PylithCallPetsc(...) \
        do { \
            PetscStackUpdateLine; \
            PetscErrorCode ierr_petsc_call = __VA_ARGS__; \
            if (PetscUnlikely(ierr_petsc_call != PETSC_SUCCESS)) { \
                auto exPtr = pylith::petsc::ErrorHandler::release(); \
                if (exPtr) { \
                    std::rethrow_exception(exPtr); \
                } \
                throw pylith::exceptions::PetscError("PETSc error with no handler-captured diagnostics.", __FILE__, __LINE__, __FUNCT__); \
            } \
        } while (0)

#define PylithCallPetscRequire(...) \
        do { \
            PetscStackUpdateLine; \
            PetscErrorCode ierr_petsc_call = __VA_ARGS__; \
            REQUIRE(!ierr_petsc_call); \
        } while (0)

// End of file
