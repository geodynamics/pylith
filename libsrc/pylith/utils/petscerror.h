// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/utils/petscerror.h
 *
 * @brief Wrappers around PETSc error handling routines.
 */

#if !defined(pylith_utils_petscerror_h)
#define pylith_utils_petscerror_h

#undef __FUNCT__
#if defined(__FUNCTION_NAME__)
#define __FUNCT__ __FUNCTION_NAME__
#else
#define __FUNCT__ "<unknown>"
#endif

#define CHECK_PETSC_ERROR(err) CHKERRXX(err)

#define CHECK_PETSC_ERROR_MSG(err, msg) \
  if (err) { \
    PetscError(PETSC_COMM_SELF, __LINE__, __FUNCT__, __FILE__, __SDIR__, err, PETSC_ERROR_IN_CXX, 0, " "); \
    throw std::runtime_error(msg); }

#endif // pylith_utils_petscerror_h

// End of file
