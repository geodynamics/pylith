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
 * @file pylith/utils/petscfwd.hh
 *
 * @brief Forward declarations for Petsc objects.
 */

#if !defined(pylith_feassemble_petscfwd_hh)
#define pylith_feassemble_petscfwd_hh


/// forward declaration for PETSc Mat
typedef struct _p_Mat* PetscMat;

/// forward declaration for PETSc Vec
typedef struct _p_Vec* PetscVec;

/// forward declaration for PETSc ISLocalToGlobalMapping
typedef struct _p_ISLocalToGlobalMapping* PetscISLocalToGlobalMapping;

/// forward declaration for PETSc PetscErrorCode
typedef int PetscErrorCode;


#endif

// End of file
