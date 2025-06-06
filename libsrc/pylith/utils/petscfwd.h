// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

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

/// Forward declaration for point function
typedef PetscPointFn* PetscPointFnWrapper;
typedef PetscPointJacFn* PetscPointJacFnWrapper;
typedef PetscBdPointFn* PetscBdPointFnWrapper;

// End of file
