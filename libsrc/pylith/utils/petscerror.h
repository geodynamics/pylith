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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/utils/error.h
 *
 * @brief Wrappers around PETSc error handling routines.
 */

#if !defined(pylith_utils_error.h)
#define pylith_utils_error.h

#undef __FUNCT__
#if defined(__FUNCTION_NAME__)
#define __FUNCT__ __FUNCTION_NAME__
#undef PETSC_FUNCTION_NAME
#define PETSC_FUNCTION_NAME __FUNCT__
#else
#define __FUNCT__ __func__
#endif

#define PYLITH_METHOD_BEGIN PetscFunctionBegin
#define PYLITH_METHOD_END PetscFunctionReturnVoid()
#define PYLITH_METHOD_RETURN(v) PetscFunctionReturn(v)

#define PYLITH_CHECK_ERROR(err) CHKERRXX(err)

#define PYLITH_CHECK_ERROR_MSG(err, msg) \
  if (err) { \
    PetscError(PETSC_COMM_SELF, __LINE__, __FUNCT__, __FILE__, __SDIR__, err, PETSC_ERROR_IN_CXX, 0, " "); \
    throw std::runtime_error(msg); }

#define PYLITH_CHECK_ERROR(err) CHKERRXX(err)

#define PYLITH_CHECK_ERROR_MSG(err, msg) \
  if (err) { \
    PetscError(PETSC_COMM_SELF, __LINE__, __FUNCT__, __FILE__, __SDIR__, err, PETSC_ERROR_IN_CXX, 0, " "); \
    throw std::runtime_error(msg); }

#endif // pylith_utils_error.h

// End of file
