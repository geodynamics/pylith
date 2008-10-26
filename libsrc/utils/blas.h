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
 * @file pylith/utils/blas.h
 *
 * @brief Declarations for BLAS routines that we use. We rely on
 * PETSc's interface to these routines.
 */

#if !defined(pylith_utils_blas_h)
#define pylith_utils_blas_h

extern "C" {
  EXTERN void dgesvd(const char*,
		     const char*,
		     PetscBLASInt*,
		     PetscBLASInt*,
		     PetscScalar*,
		     PetscBLASInt*,
		     PetscReal*,
		     PetscScalar*,
		     PetscBLASInt*,
		     PetscScalar*,
		     PetscBLASInt*,
		     PetscScalar*,
		     PetscBLASInt*,
		     PetscBLASInt*);
}

#endif // pylith_utils_blas_h


// End of file
