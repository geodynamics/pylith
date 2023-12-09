// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file libsrc/utils/pylithtypes.h
 *
 * @brief Type definitions for PyLith.
 */

#if !defined(pylith_utils_pylithtypes_h)
#define pylith_utils_pylithtypes_h

#include "petsc.h"

typedef PetscScalar PylithScalar;
typedef PetscReal PylithReal;
typedef PetscInt PylithInt;

typedef PetscErrorCode (*PetscUserFieldFunc)(PetscInt,
                                             PetscReal,
                                             const PetscReal x[],
                                             PetscInt,
                                             PetscScalar *u,
                                             void *ctx);

#endif // pylith_utils_pylithtypes_h

// End of file
