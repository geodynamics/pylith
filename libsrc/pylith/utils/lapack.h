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
