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
 * @file pylith/utils/petscerror.h
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

#define CHECK_PETSC_ERROR(err) \
  if (err) { \
    PetscError(__LINE__, __FUNCT__, __FILE__, __SDIR__, err, 0, " "); \
    throw std::runtime_error("PETSc error."); }

#define CHECK_PETSC_ERROR_MSG(err, msg) \
  if (err) { \
    PetscError(__LINE__, __FUNCT__, __FILE__, __SDIR__, err, 0, " "); \
    throw std::runtime_error(msg); }

#endif // pylith_utils_petscerror_h

// End of file
