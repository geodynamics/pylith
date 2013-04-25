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

#include <portinfo>

#include "Solver.hh" // implementation of class methods

#include "Formulation.hh" // USES Formulation

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

#include <petscdmmesh_solvers.hh> // USES constructFieldSplit()

#include <cassert> // USES assert()


// ----------------------------------------------------------------------
EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MyMatGetSubMatrix"
PetscErrorCode  MyMatGetSubMatrix(PetscMat mat,
				  PetscIS isrow,
				  PetscIS iscol,
				  MatReuse reuse,
				  PetscMat *newmat) {
  FaultPreconCtx *ctx = NULL;
  PetscIS faultIS = NULL;
  PetscBool isFaultRow, isFaultCol;
  PetscErrorCode  ierr = 0;

  PetscFunctionBegin;
  ierr = MatShellGetContext(mat, (void **) &ctx);CHKERRQ(ierr);
  ierr = PCFieldSplitGetIS(ctx->pc, ctx->faultFieldName, &faultIS);CHKERRQ(ierr);assert(faultIS);
  ierr = ISEqual(isrow, faultIS, &isFaultRow);CHKERRQ(ierr);
  ierr = ISEqual(iscol, faultIS, &isFaultCol);CHKERRQ(ierr);
  if (isFaultRow && isFaultCol) {
    if (reuse == MAT_INITIAL_MATRIX) {
      ierr = PetscObjectReference((PetscObject) ctx->faultA);CHKERRQ(ierr);
      *newmat = ctx->faultA;
    } // if
  } else {
    ierr = MatGetSubMatrix(ctx->A, isrow, iscol, reuse, newmat);CHKERRQ(ierr);
  } // if/else

  PetscFunctionReturn(0);
} // MyMatGetSubMatrix
EXTERN_C_END

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Solver::Solver(void) :
  _formulation(0),
  _logger(0),
  _jacobianPC(0),
  _jacobianPCFault(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Solver::~Solver(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Solver::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  _formulation = 0; // Handle only, do not manage memory.
  delete _logger; _logger = 0;

  PetscErrorCode err = 0;
  err = MatDestroy(&_jacobianPC);PYLITH_CHECK_ERROR(err);
  err = MatDestroy(&_jacobianPCFault);PYLITH_CHECK_ERROR(err);

  _ctx.pc = 0; // KSP PC (managed separately)
  _ctx.A = 0; // Jacobian (managed separately)
  _ctx.faultA  = 0; // Handle to _jacobianPCFault

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::Solver::initialize(const topology::SolutionFields& fields,
				     const topology::Jacobian& jacobian,
				     Formulation* const formulation)
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(formulation);
  _formulation = formulation;

  // Make global preconditioner matrix
  PetscMat jacobianMat = jacobian.matrix();

  PetscSection solutionSection = fields.solution().petscSection();assert(solutionSection);
  PetscDM dmMesh = fields.solution().mesh().dmMesh();assert(dmMesh);
  PetscInt numFields;
  PetscErrorCode err;

  err = DMGetNumFields(dmMesh, &numFields);PYLITH_CHECK_ERROR(err);
  if (formulation->splitFields() && formulation->useCustomConstraintPC() &&
      numFields > 1) {
    // We have split fields with a custom constraint preconditioner
    // and constraints exist.

    PylithInt M, N, m, n;
    err = MatDestroy(&_jacobianPC);PYLITH_CHECK_ERROR(err);
    err = MatGetSize(jacobianMat, &M, &N);PYLITH_CHECK_ERROR(err);
    err = MatGetLocalSize(jacobianMat, &m, &n);PYLITH_CHECK_ERROR(err);
    err = MatCreateShell(fields.mesh().comm(), m, n, M, N, &_ctx, &_jacobianPC);PYLITH_CHECK_ERROR(err);
    err = MatShellSetOperation(_jacobianPC, MATOP_GET_SUBMATRIX, 
                               (void (*)(void)) MyMatGetSubMatrix);PYLITH_CHECK_ERROR(err);
  } else {
    _jacobianPC = jacobianMat;
    err = PetscObjectReference((PetscObject) jacobianMat);PYLITH_CHECK_ERROR(err);
  } // if/else

  PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Setup preconditioner for preconditioning with split fields.
void
pylith::problems::Solver::_setupFieldSplit(PetscPC* const pc,
					   Formulation* const formulation,
					   const topology::Jacobian& jacobian,
					   const topology::SolutionFields& fields)
{ // _setupFieldSplit
  PYLITH_METHOD_BEGIN;

  assert(pc);
  assert(formulation);

  PetscDM dmMesh = fields.mesh().dmMesh();assert(dmMesh);
  MPI_Comm comm;
  PetscSection solutionSection = fields.solution().petscSection();assert(solutionSection);
  PetscVec solutionVec = fields.solution().localVector();assert(solutionVec);
  PetscVec solutionGlobalVec = fields.solution().globalVector();assert(solutionGlobalVec);
  PetscInt spaceDim, numFields;
  PetscErrorCode err;

  const bool separateComponents = formulation->splitFieldComponents();

  err = PetscObjectGetComm((PetscObject) dmMesh, &comm);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetDimension(dmMesh, &spaceDim);PYLITH_CHECK_ERROR(err);
  err = PetscSectionGetNumFields(solutionSection, &numFields);PYLITH_CHECK_ERROR(err);

  err = PCSetDM(*pc, dmMesh);PYLITH_CHECK_ERROR(err);
  err = PCSetType(*pc, PCFIELDSPLIT);PYLITH_CHECK_ERROR(err);
  err = PCSetOptionsPrefix(*pc, "fs_");PYLITH_CHECK_ERROR(err);
  err = PCSetFromOptions(*pc);PYLITH_CHECK_ERROR(err);

  if (formulation->splitFields() && formulation->useCustomConstraintPC() &&
      numFields > 1) {
    // We have split fields with a custom constraint preconditioner
    // and constraints exist.

    // Get total number of DOF associated with constraints field split
    PetscDM lagrangeDM = NULL;
    PetscSection lagrangeSection = NULL;
    PetscInt lagrangeFields[1] = {1}, nrows, ncols;

    err = MatDestroy(&_jacobianPCFault);PYLITH_CHECK_ERROR(err);
    err = DMCreateSubDM(dmMesh, 1, lagrangeFields, NULL, &lagrangeDM);PYLITH_CHECK_ERROR(err);
    err = DMGetDefaultGlobalSection(lagrangeDM, &lagrangeSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetStorageSize(lagrangeSection, &nrows);PYLITH_CHECK_ERROR(err);
    ncols = nrows;

    err = MatCreate(comm, &_jacobianPCFault);PYLITH_CHECK_ERROR(err);
    err = MatSetSizes(_jacobianPCFault, nrows, ncols, PETSC_DECIDE, PETSC_DECIDE);PYLITH_CHECK_ERROR(err);
    err = MatSetType(_jacobianPCFault, MATAIJ);PYLITH_CHECK_ERROR(err);
    err = MatSetFromOptions(_jacobianPCFault);PYLITH_CHECK_ERROR(err);
    
    // Allocate just the diagonal.
    err = MatSeqAIJSetPreallocation(_jacobianPCFault, 1, NULL);PYLITH_CHECK_ERROR(err);
    err = MatMPIAIJSetPreallocation(_jacobianPCFault, 1, NULL, 0, NULL);PYLITH_CHECK_ERROR(err);
    // Set preconditioning matrix in formulation
    formulation->customPCMatrix(_jacobianPCFault);assert(_jacobianPCFault);

    _ctx.pc = *pc;
    _ctx.A = jacobian.matrix();
    switch (spaceDim) {
    case 1 :
      _ctx.faultFieldName = "1";
      break;
    case 2 :
      _ctx.faultFieldName = (separateComponents) ? "2" : "1";
      break;
    case 3 :
      _ctx.faultFieldName = (separateComponents) ? "3" : "1";
      break;
    default:
      assert(0);
      throw std::logic_error("Unknown space dimension in "
			     "Problems::_setupFieldSplit().");
    } // switch
    _ctx.faultA = _jacobianPCFault;
  } // if

  if (separateComponents) {
    throw std::logic_error("Separate components is no longer supported");
  } else {
    // Create rigid body null space.
    MatNullSpace nullsp = NULL;    
    PetscObject  field = NULL;
    PetscSection coordinateSection;
    PetscVec  coordinateVec;
    PetscInt dim = spaceDim, vStart, vEnd;
    PetscVec mode[6];

    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetCoordinateSection(dmMesh, &coordinateSection);PYLITH_CHECK_ERROR(err);assert(coordinateSection);
    err = DMGetCoordinatesLocal(dmMesh, &coordinateVec);PYLITH_CHECK_ERROR(err);assert(coordinateVec);
    if (dim > 1) {
      const int m = (dim * (dim + 1)) / 2;
      for(int i = 0; i < m; ++i) {
        err = VecDuplicate(solutionGlobalVec, &mode[i]);PYLITH_CHECK_ERROR(err);
      } // for
      // :KLUDGE: Assume P1
      for(int d = 0; d < dim; ++d) {
        PetscScalar values[3] = {0.0, 0.0, 0.0};
	
        values[d] = 1.0;
        err = VecSet(solutionVec, 0.0);PYLITH_CHECK_ERROR(err);
        for(PetscInt v = vStart; v < vEnd; ++v) {
          err = DMPlexVecSetClosure(dmMesh, solutionSection, solutionVec, v, values, INSERT_VALUES);PYLITH_CHECK_ERROR(err);
        } // for
        err = DMLocalToGlobalBegin(dmMesh, solutionVec, INSERT_VALUES, mode[d]);PYLITH_CHECK_ERROR(err);
        err = DMLocalToGlobalEnd(dmMesh, solutionVec, INSERT_VALUES, mode[d]);PYLITH_CHECK_ERROR(err);
      } // for
      for(int d = dim; d < m; ++d) {
        PetscInt k = (dim > 2) ? d - dim : d;
	
        err = VecSet(solutionVec, 0.0);PYLITH_CHECK_ERROR(err);
        for(PetscInt v = vStart; v < vEnd; ++v) {
          PetscScalar  values[3] = {0.0, 0.0, 0.0};
          PetscScalar *coords = NULL;

          err = DMPlexVecGetClosure(dmMesh, coordinateSection, coordinateVec, v, NULL, &coords);PYLITH_CHECK_ERROR(err);
          for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
              values[j] += _epsilon(i, j, k)*coords[i];
            } // for
          } // for
          err = DMPlexVecRestoreClosure(dmMesh, coordinateSection, coordinateVec, v, NULL, &coords);PYLITH_CHECK_ERROR(err);
          err = DMPlexVecSetClosure(dmMesh, solutionSection, solutionVec, v, values, INSERT_VALUES);PYLITH_CHECK_ERROR(err);
        }
        err = DMLocalToGlobalBegin(dmMesh, solutionVec, INSERT_VALUES, mode[d]);PYLITH_CHECK_ERROR(err);
        err = DMLocalToGlobalEnd(dmMesh, solutionVec, INSERT_VALUES, mode[d]);PYLITH_CHECK_ERROR(err);
      } // for
      for(int i = 0; i < dim; ++i) {
        err = VecNormalize(mode[i], NULL);PYLITH_CHECK_ERROR(err);
      } // for
      /* Orthonormalize system */
      for(int i = dim; i < m; ++i) {
        PetscScalar dots[6];

        err = VecMDot(mode[i], i, mode, dots);PYLITH_CHECK_ERROR(err);
        for(int j = 0; j < i; ++j) dots[j] *= -1.0;
        err = VecMAXPY(mode[i], i, dots, mode);PYLITH_CHECK_ERROR(err);
        err = VecNormalize(mode[i], NULL);PYLITH_CHECK_ERROR(err);
      } // for
      err = MatNullSpaceCreate(comm, PETSC_FALSE, m, mode, &nullsp);PYLITH_CHECK_ERROR(err);
      for(int i = 0; i< m; ++i) {err = VecDestroy(&mode[i]);PYLITH_CHECK_ERROR(err);}
    } // if

    err = DMSetNumFields(dmMesh, 2);PYLITH_CHECK_ERROR(err);
    err = DMGetField(dmMesh, 0, &field);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose(field, "nearnullspace", (PetscObject) nullsp);PYLITH_CHECK_ERROR(err);
    err = MatNullSpaceDestroy(&nullsp);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose(field, "pmat", (PetscObject) _jacobianPCFault);PYLITH_CHECK_ERROR(err);
  } // if/else

  PYLITH_METHOD_END;
} // _setupFieldSplit

// ----------------------------------------------------------------------
int
pylith::problems::Solver::_epsilon(int i, 
				   int j, 
				   int k)
{ // _epsilon
  switch(i) {
  case 0:
    switch(j) {
    case 0: return 0;
    case 1:
      switch(k) {
      case 0: return 0;
      case 1: return 0;
      case 2: return 1;
      }
    case 2:
      switch(k) {
      case 0: return 0;
      case 1: return -1;
      case 2: return 0;
      }
    }
  case 1:
    switch(j) {
    case 0:
      switch(k) {
      case 0: return 0;
      case 1: return 0;
      case 2: return -1;
      }
    case 1: return 0;
    case 2:
      switch(k) {
      case 0: return 1;
      case 1: return 0;
      case 2: return 0;
      }
    }
  case 2:
    switch(j) {
    case 0:
      switch(k) {
      case 0: return 0;
      case 1: return 1;
      case 2: return 0;
      }
    case 1:
      switch(k) {
      case 0: return -1;
      case 1: return 0;
      case 2: return 0;
      }
    case 2: return 0;
    }
  }
  return 0;
} // _epsilon


// End of file
