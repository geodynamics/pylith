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

#include <portinfo>

#include "Solver.hh" // implementation of class methods

#include "Formulation.hh" // USES Formulation

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()


// ----------------------------------------------------------------------
EXTERN_C_BEGIN
PetscErrorCode  MyMatGetSubMatrix(PetscMat mat,
				  PetscIS isrow,
				  PetscIS iscol,
				  MatReuse reuse,
				  PetscMat *newmat) {
  FaultPreconCtx *ctx = NULL;
  PetscIS faultIS = NULL;
  PetscBool isFaultRow, isFaultCol;
  PetscErrorCode  ierr = 0;

  PYLITH_METHOD_BEGIN;

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

  PYLITH_METHOD_RETURN(0);
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

  PetscSection solutionSection = fields.solution().localSection();assert(solutionSection);
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
// Create null space.
void
pylith::problems::Solver::_createNullSpace(const topology::SolutionFields& fields)
{ // _createNullSpace
  PYLITH_METHOD_BEGIN;

  PetscDM dmMesh = fields.solution().dmMesh();assert(dmMesh);
  PetscErrorCode err;

  const spatialdata::geocoords::CoordSys* cs = fields.mesh().coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  MPI_Comm comm;
  err = PetscObjectGetComm((PetscObject) dmMesh, &comm);PYLITH_CHECK_ERROR(err);

  PetscSection solutionSection = fields.solution().localSection();assert(solutionSection);
  PetscVec solutionVec = fields.solution().localVector();assert(solutionVec);
  PetscVec solutionGlobalVec = fields.solution().globalVector();assert(solutionGlobalVec);
  MatNullSpace nullsp = NULL;    
  PetscSection coordinateSection = NULL;
  PetscVec coordinateVec = NULL;
  PetscInt vStart, vEnd;
  PetscVec mode[6];
  
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
  err = DMGetCoordinateSection(dmMesh, &coordinateSection);PYLITH_CHECK_ERROR(err);assert(coordinateSection);
  err = DMGetCoordinatesLocal(dmMesh, &coordinateVec);PYLITH_CHECK_ERROR(err);assert(coordinateVec);
  if (spaceDim > 1) {
    const int m = (spaceDim * (spaceDim + 1)) / 2;
    for(int i = 0; i < m; ++i) {
      err = VecDuplicate(solutionGlobalVec, &mode[i]);PYLITH_CHECK_ERROR(err);
      // This is necessary to avoid circular references when we compose this MatNullSpace with a field in the DM
      err = VecSetDM(mode[i], NULL);PYLITH_CHECK_ERROR(err);
    } // for
    // :KLUDGE: Assume P1
    for(int d = 0; d < spaceDim; ++d) {
      PetscScalar values[3] = {0.0, 0.0, 0.0};      
      values[d] = 1.0;

      err = VecSet(solutionVec, 0.0);PYLITH_CHECK_ERROR(err);
      for(PetscInt v = vStart; v < vEnd; ++v) {
        err = DMPlexVecSetClosure(dmMesh, solutionSection, solutionVec, v, values, INSERT_VALUES);PYLITH_CHECK_ERROR(err);
      } // for
      err = DMLocalToGlobalBegin(dmMesh, solutionVec, INSERT_VALUES, mode[d]);PYLITH_CHECK_ERROR(err);
      err = DMLocalToGlobalEnd(dmMesh, solutionVec, INSERT_VALUES, mode[d]);PYLITH_CHECK_ERROR(err);
    } // for
    for(int d = spaceDim; d < m; ++d) {
      PetscInt k = (spaceDim > 2) ? d - spaceDim : d;
      
      err = VecSet(solutionVec, 0.0);PYLITH_CHECK_ERROR(err);
      for(PetscInt v = vStart; v < vEnd; ++v) {
        PetscScalar  values[3] = {0.0, 0.0, 0.0};
        PetscScalar *coords = NULL;
	
        err = DMPlexVecGetClosure(dmMesh, coordinateSection, coordinateVec, v, NULL, &coords);PYLITH_CHECK_ERROR(err);
        for(int i = 0; i < spaceDim; ++i) {
          for(int j = 0; j < spaceDim; ++j) {
            values[j] += _epsilon(i, j, k)*coords[i];
          } // for
        } // for
        err = DMPlexVecRestoreClosure(dmMesh, coordinateSection, coordinateVec, v, NULL, &coords);PYLITH_CHECK_ERROR(err);
        err = DMPlexVecSetClosure(dmMesh, solutionSection, solutionVec, v, values, INSERT_VALUES);PYLITH_CHECK_ERROR(err);
      } // for
      err = DMLocalToGlobalBegin(dmMesh, solutionVec, INSERT_VALUES, mode[d]);PYLITH_CHECK_ERROR(err);
      err = DMLocalToGlobalEnd(dmMesh, solutionVec, INSERT_VALUES, mode[d]);PYLITH_CHECK_ERROR(err);
    } // for
    for(int i = 0; i < spaceDim; ++i) {
      err = VecNormalize(mode[i], NULL);PYLITH_CHECK_ERROR(err);
    } // for
    // Orthonormalize system
    for(int i = spaceDim; i < m; ++i) {
      PetscScalar dots[6];
      
      err = VecMDot(mode[i], i, mode, dots);PYLITH_CHECK_ERROR(err);
      for(int j = 0; j < i; ++j) {
        dots[j] *= -1.0;
      } // for
      err = VecMAXPY(mode[i], i, dots, mode);PYLITH_CHECK_ERROR(err);
      err = VecNormalize(mode[i], NULL);PYLITH_CHECK_ERROR(err);
    } // for
    err = MatNullSpaceCreate(comm, PETSC_FALSE, m, mode, &nullsp);PYLITH_CHECK_ERROR(err);
    for(int i = 0; i< m; ++i) {
      err = VecDestroy(&mode[i]);PYLITH_CHECK_ERROR(err);
    } // for
  } // if
  
  PetscObject field = NULL;
  PetscInt    numFields;
  err = DMGetNumFields(dmMesh, &numFields);PYLITH_CHECK_ERROR(err);
  if (!numFields) {
    err = DMSetNumFields(dmMesh, 1);PYLITH_CHECK_ERROR(err);
  } // if
  err = DMGetField(dmMesh, 0, &field);PYLITH_CHECK_ERROR(err);
  err = PetscObjectCompose(field, "nearnullspace", (PetscObject) nullsp);PYLITH_CHECK_ERROR(err);
  err = MatNullSpaceDestroy(&nullsp);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _createNullSpace

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

  PetscDM dmMesh = fields.solution().dmMesh();assert(dmMesh);
  MPI_Comm comm;
  PetscSection solutionSection = fields.solution().localSection();assert(solutionSection);
  PetscVec solutionVec = fields.solution().localVector();assert(solutionVec);
  PetscVec solutionGlobalVec = fields.solution().globalVector();assert(solutionGlobalVec);
  PetscInt numFields;
  PetscErrorCode err;

  err = PetscObjectGetComm((PetscObject) dmMesh, &comm);PYLITH_CHECK_ERROR(err);

  err = PCSetDM(*pc, dmMesh);PYLITH_CHECK_ERROR(err);
  err = PCSetType(*pc, PCFIELDSPLIT);PYLITH_CHECK_ERROR(err);
  err = PCSetOptionsPrefix(*pc, "fs_");PYLITH_CHECK_ERROR(err);
  err = PCSetFromOptions(*pc);PYLITH_CHECK_ERROR(err);

  if (formulation->splitFields() && formulation->useCustomConstraintPC()) {
    err = PetscSectionGetNumFields(solutionSection, &numFields);PYLITH_CHECK_ERROR(err);
    if (numFields < 2) {
      throw std::logic_error("Cannot setup solution field for split fields. Solution field has only one field.");
    } // if

    // We have split fields with a custom constraint preconditioner
    // and constraints exist.

    // Get total number of DOF associated with constraints field split
    PetscInt nrows = 0;
    PetscSection lagrangeSection = NULL;

    err = MatDestroy(&_jacobianPCFault);PYLITH_CHECK_ERROR(err);

    PetscDM lagrangeDM = fields.solution().subfieldInfo("lagrange_multiplier").dm;assert(lagrangeDM);
    err = DMGetDefaultGlobalSection(lagrangeDM, &lagrangeSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetStorageSize(lagrangeSection, &nrows);PYLITH_CHECK_ERROR(err);
    PetscInt ncols = nrows;

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
    _ctx.faultFieldName = "1";
    _ctx.faultA = _jacobianPCFault;
  } // if
    
  PetscObject field = NULL;
  err = DMSetNumFields(dmMesh, 2);PYLITH_CHECK_ERROR(err);
  err = DMGetField(dmMesh, 1, &field);PYLITH_CHECK_ERROR(err);
  err = PetscObjectCompose(field, "pmat", (PetscObject) _jacobianPCFault);PYLITH_CHECK_ERROR(err);

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
