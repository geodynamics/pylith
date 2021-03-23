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

/// forward declaration for PETSc TS
typedef struct _p_TS* PetscTS;

/// forward declaration for PETSc PC
typedef struct _p_PC* PetscPC;

/// forward declaration for PETSc DM
typedef struct _p_DM* PetscDM;

/// forward declaration for PETSc DMLabel
typedef struct _p_DMLabel* PetscDMLabel;

/// forward declaration for PETSc IS
typedef struct _p_IS* PetscIS;

/// forward declaration for PETSc ISLocalToGlobalMapping
typedef struct _p_ISLocalToGlobalMapping* PetscISLocalToGlobalMapping;

/// forward declaration for PETSc DMMeshInterpolationInfo
typedef struct _DMMeshInterpolationInfo* PetscDMMeshInterpolationInfo;

/// forward declaration for PETSc DS
typedef struct _p_PetscDS* PetscDS;

/// forward declaration for PETSc FE
typedef struct _p_PetscFE* PetscFE;

/// forward declaration for PETSc weak form
typedef struct _p_PetscWeakForm* PetscWeakForm;

#endif // pylith_utils_petscfwd_h

// End of file
