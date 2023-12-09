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
 * @file libsrc/utils/lapack.h
 *
 * @brief Declarations for LAPACK routines that we use. We rely on
 * PETSc's interface to these routines.
 */

#if !defined(pylith_utils_lapack_h)
#define pylith_utils_lapack_h

#include <petscblaslapack.h>
#define lapack_dgesvd LAPACKgesvd_

#endif // pylith_utils_lapack_h

// End of file
