// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

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
