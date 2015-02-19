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
 * @file libsrc/utils/petscfwd.h
 *
 * @brief Forward declarations for Petsc objects.
 */

#if !defined(pylith_utils_petscfwd_h)
#define pylith_utils_petscfwd_h

/// forward declaration for PETSc PetscErrorCode
typedef int PetscErrorCode;

/// forward declaration for PETSc Mat
typedef struct _p_Mat* PetscMat;

/// forward declaration for PETSc Vec
typedef struct _p_Vec* PetscVec;

/// forward declaration for PETSc VecScatter
typedef struct _p_VecScatter* PetscVecScatter;

/// forward declaration for PETSc KSP
typedef struct _p_KSP* PetscKSP;

/// forward declaration for PETSc SNES
typedef struct _p_SNES* PetscSNES;

/// forward declatation for PETSc line search
typedef struct _p_LineSearch* PetscSNESLineSearch;

/// forward declaration for PETSc PC
typedef struct _p_PC* PetscPC;

/// forward declaration for PETSc DM
typedef struct _p_DM* PetscDM;

/// forward declaration for PETSc DMLabel
typedef struct _n_DMLabel* PetscDMLabel;

/// forward declaration for PETSc IS
typedef struct _p_IS* PetscIS;

/// forward declaration for PETSc ISLocalToGlobalMapping
typedef struct _p_ISLocalToGlobalMapping* PetscISLocalToGlobalMapping;

/// forward declaration for PETSc DMMeshInterpolationInfo
typedef struct _DMMeshInterpolationInfo* PetscDMMeshInterpolationInfo;


#endif // pylith_utils_petscfwd_h

// End of file
