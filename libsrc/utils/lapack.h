// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
