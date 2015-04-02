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
// Copyright (c) 2010-2015 University of California, Davis
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

#if !defined(pylith_utils_error_h)
#define pylith_utils_error_h

#include <cassert>
#include <stdexcept>
#include <sstream>

#undef __FUNCT__
#if defined(__FUNCTION_NAME__)
#define __FUNCT__ __FUNCTION_NAME__
#undef PETSC_FUNCTION_NAME
#define PETSC_FUNCTION_NAME __FUNCT__
#else
#define __FUNCT__ __func__
#endif

#define PYLITH_METHOD_BEGIN PetscFunctionBeginUser
#define PYLITH_METHOD_END PetscFunctionReturnVoid()
#define PYLITH_METHOD_RETURN(v) PetscFunctionReturn(v)

#define PYLITH_CHECK_ERROR(err) do {if (PetscUnlikely(err)) {PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,err,PETSC_ERROR_REPEAT,0);throw std::runtime_error("Error detected while in PETSc function.");}} while(0)

#define PYLITH_CHECK_ERROR_MSG(err, msg) \
  if (err) { \
    PetscError(PETSC_COMM_SELF,__LINE__,__FUNCT__,__FILE__,err,PETSC_ERROR_REPEAT, 0, " "); \
    throw std::runtime_error(msg); }

#endif // pylith_utils_error_h

// End of file
