// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

/*
 * @brief Declarations for LAPACK routines that we use. We rely on
 * PETSc's interface to these routines.
 */

#include <petscblaslapack.h>
#define lapack_dgesvd LAPACKgesvd_

// End of file
