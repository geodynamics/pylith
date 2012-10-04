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
// Copyright (c) 2010-2012 University of California, Davis
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
#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

#include <petscdmmesh_solvers.hh> // USES constructFieldSplit()

#include <cassert> // USES assert()


// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MyMatGetSubMatrix"
PetscErrorCode  MyMatGetSubMatrix(Mat mat, IS isrow, IS iscol, MatReuse reuse, Mat *newmat) {
  FaultPreconCtx *ctx;
  IS              faultIS;
  PetscBool       isFaultRow, isFaultCol;
  PetscErrorCode  ierr = 0;

  PetscFunctionBegin;
  ierr = MatShellGetContext(mat, (void **) &ctx);CHKERRQ(ierr);
  ierr = PCFieldSplitGetIS(ctx->pc, ctx->faultFieldName, &faultIS);CHKERRQ(ierr);
  assert(faultIS);
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
  _formulation = 0; // Handle only, do not manage memory.
  delete _logger; _logger = 0;

  PetscErrorCode err = 0;
  err = MatDestroy(&_jacobianPC);CHECK_PETSC_ERROR(err);
  err = MatDestroy(&_jacobianPCFault);CHECK_PETSC_ERROR(err);

  _ctx.pc = 0; // KSP PC (managed separately)
  _ctx.A = 0; // Jacobian (managed separately)
  _ctx.faultA  = 0; // Handle to _jacobianPCFault
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::Solver::initialize(const topology::SolutionFields& fields,
				     const topology::Jacobian& jacobian,
				     Formulation* const formulation)
{ // initialize
  assert(formulation);
  _formulation = formulation;

  // Make global preconditioner matrix
  PetscMat jacobianMat = jacobian.matrix();

  PetscSection   solutionSection = fields.solution().petscSection();
  DM             dmMesh          = fields.solution().mesh().dmMesh();
  PetscInt       numFields;
  PetscErrorCode err;

  assert(solutionSection);assert(dmMesh);
  err = DMGetNumFields(dmMesh, &numFields);CHECK_PETSC_ERROR(err);
  if (formulation->splitFields() && formulation->useCustomConstraintPC() &&
      numFields > 1) {
    // We have split fields with a custom constraint preconditioner
    // and constraints exist.

    PylithInt M, N, m, n;
    PetscErrorCode err = 0;
    err = MatDestroy(&_jacobianPC);CHECK_PETSC_ERROR(err);
    err = MatGetSize(jacobianMat, &M, &N);CHECK_PETSC_ERROR(err);
    err = MatGetLocalSize(jacobianMat, &m, &n);CHECK_PETSC_ERROR(err);
    err = MatCreateShell(fields.mesh().comm(), m, n, M, N, &_ctx, &_jacobianPC);
    CHECK_PETSC_ERROR(err);
    err = MatShellSetOperation(_jacobianPC, MATOP_GET_SUBMATRIX, 
                               (void (*)(void)) MyMatGetSubMatrix);
    CHECK_PETSC_ERROR(err);
  } else {
    _jacobianPC = jacobianMat;
    PetscErrorCode err = PetscObjectReference((PetscObject) jacobianMat);
    CHECK_PETSC_ERROR(err);
  } // if/else
} // initialize


static inline PetscInt epsilon(PetscInt i, PetscInt j, PetscInt k)
{
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
}

// ----------------------------------------------------------------------
// Setup preconditioner for preconditioning with split fields.
void
pylith::problems::Solver::_setupFieldSplit(PetscPC* const pc,
					   Formulation* const formulation,
					   const topology::Jacobian& jacobian,
					   const topology::SolutionFields& fields)
{ // _setupFieldSplit
  assert(pc);
  assert(formulation);

  DM dmMesh = fields.mesh().dmMesh();
  assert(dmMesh);
  PetscSection   solutionSection   = fields.solution().petscSection();
  Vec            solutionVec       = fields.solution().localVector();
  Vec            solutionGlobalVec = fields.solution().globalVector();
  PetscInt       spaceDim, numFields;
  PetscErrorCode err;


  const bool separateComponents = formulation->splitFieldComponents();

  assert(solutionSection);
  err = DMComplexGetDimension(dmMesh, &spaceDim);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetNumFields(solutionSection, &numFields);CHECK_PETSC_ERROR(err);

  err = PCSetDM(*pc, dmMesh);CHECK_PETSC_ERROR(err);
  err = PCSetType(*pc, PCFIELDSPLIT);CHECK_PETSC_ERROR(err);
  err = PCSetOptionsPrefix(*pc, "fs_");CHECK_PETSC_ERROR(err);
  err = PCSetFromOptions(*pc);CHECK_PETSC_ERROR(err);

  if (formulation->splitFields() && formulation->useCustomConstraintPC() &&
      numFields > 1) {
    // We have split fields with a custom constraint preconditioner
    // and constraints exist.

    // Get total number of DOF associated with constraints field split
    DM           lagrangeDM;
    PetscSection lagrangeSection;
    PetscInt     lagrangeFields[1] = {1}, nrows, ncols;

    err = MatDestroy(&_jacobianPCFault);CHECK_PETSC_ERROR(err);
    err = DMCreateSubDM(dmMesh, 1, lagrangeFields, PETSC_NULL, &lagrangeDM);CHECK_PETSC_ERROR(err);
    err = DMGetDefaultGlobalSection(lagrangeDM, &lagrangeSection);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetStorageSize(lagrangeSection, &nrows);CHECK_PETSC_ERROR(err);
    ncols = nrows;

    err = MatCreate(((PetscObject) dmMesh)->comm, &_jacobianPCFault);CHECK_PETSC_ERROR(err);
    err = MatSetSizes(_jacobianPCFault, nrows, ncols, 
                      PETSC_DECIDE, PETSC_DECIDE);CHECK_PETSC_ERROR(err);
    err = MatSetType(_jacobianPCFault, MATAIJ);CHECK_PETSC_ERROR(err);
    err = MatSetFromOptions(_jacobianPCFault);CHECK_PETSC_ERROR(err);
    
    // Allocate just the diagonal.
    err = MatSeqAIJSetPreallocation(_jacobianPCFault, 1, 
                                    PETSC_NULL);CHECK_PETSC_ERROR(err);
    err = MatMPIAIJSetPreallocation(_jacobianPCFault, 1, PETSC_NULL, 
                                    0, PETSC_NULL);CHECK_PETSC_ERROR(err);
    // Set preconditioning matrix in formulation
    formulation->customPCMatrix(_jacobianPCFault);

    assert(_jacobianPCFault);

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
    throw std::runtime_error("Separate components is no longer supported");
    // constructFieldSplit(solutionSection, PETSC_DETERMINE, PETSC_NULL, PETSC_NULL, 
	//		sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", solutionSection), precon, PETSC_NULL, solution.vector(), *pc);
  } else {
    // Create rigid body null space.
    MatNullSpace nullsp;    
    PetscObject  field;
    PetscSection coordinateSection;
    Vec          coordinateVec;
    PetscInt     dim = spaceDim, vStart, vEnd;
    PetscVec     mode[6];

    err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
    err = DMComplexGetCoordinateSection(dmMesh, &coordinateSection);CHECK_PETSC_ERROR(err);
    err = DMGetCoordinatesLocal(dmMesh, &coordinateVec);CHECK_PETSC_ERROR(err);
    assert(coordinateSection);assert(coordinateVec);
    if (dim > 1) {
      const int m = (dim * (dim + 1)) / 2;
      for(int i = 0; i < m; ++i) {
        err = VecDuplicate(solutionGlobalVec, &mode[i]);CHECK_PETSC_ERROR(err);
      } // for
      // :KLUDGE: Assume P1
      for(int d = 0; d < dim; ++d) {
        PetscScalar values[3] = {0.0, 0.0, 0.0};
	
        values[d] = 1.0;
        err = VecSet(solutionVec, 0.0);CHECK_PETSC_ERROR(err);
        for(PetscInt v = vStart; v < vEnd; ++v) {
          err = DMComplexVecSetClosure(dmMesh, solutionSection, solutionVec, v, values, INSERT_VALUES);CHECK_PETSC_ERROR(err);
        } // for
        err = DMLocalToGlobalBegin(dmMesh, solutionVec, INSERT_VALUES, mode[d]);CHECK_PETSC_ERROR(err);
        err = DMLocalToGlobalEnd(dmMesh, solutionVec, INSERT_VALUES, mode[d]);CHECK_PETSC_ERROR(err);
      } // for
      for(int d = dim; d < m; ++d) {
        PetscInt k = (dim > 2) ? d - dim : d;
	
        err = VecSet(solutionVec, 0.0);CHECK_PETSC_ERROR(err);
        for(PetscInt v = vStart; v < vEnd; ++v) {
          PetscScalar values[3] = {0.0, 0.0, 0.0};
          const PetscScalar *coords;

          err = DMComplexVecGetClosure(dmMesh, coordinateSection, coordinateVec, v, PETSC_NULL, &coords);CHECK_PETSC_ERROR(err);
          for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
              values[j] += _epsilon(i, j, k)*coords[i];
            } // for
          } // for
          err = DMComplexVecRestoreClosure(dmMesh, coordinateSection, coordinateVec, v, PETSC_NULL, &coords);CHECK_PETSC_ERROR(err);
          err = DMComplexVecSetClosure(dmMesh, solutionSection, solutionVec, v, values, INSERT_VALUES);CHECK_PETSC_ERROR(err);
        }
        err = DMLocalToGlobalBegin(dmMesh, solutionVec, INSERT_VALUES, mode[d]);CHECK_PETSC_ERROR(err);
        err = DMLocalToGlobalEnd(dmMesh, solutionVec, INSERT_VALUES, mode[d]);CHECK_PETSC_ERROR(err);
      } // for
      for(int i = 0; i < dim; ++i) {
        err = VecNormalize(mode[i], PETSC_NULL);CHECK_PETSC_ERROR(err);
      } // for
      /* Orthonormalize system */
      for(int i = dim; i < m; ++i) {
        PetscScalar dots[6];

        err = VecMDot(mode[i], i, mode, dots);CHECK_PETSC_ERROR(err);
        for(int j = 0; j < i; ++j) dots[j] *= -1.0;
        err = VecMAXPY(mode[i], i, dots, mode);CHECK_PETSC_ERROR(err);
        err = VecNormalize(mode[i], PETSC_NULL);CHECK_PETSC_ERROR(err);
      } // for
      err = MatNullSpaceCreate(((PetscObject) dmMesh)->comm, PETSC_FALSE, m, mode, &nullsp);CHECK_PETSC_ERROR(err);
      for(int i = 0; i< m; ++i) {err = VecDestroy(&mode[i]);CHECK_PETSC_ERROR(err);}
    } // if

    err = DMSetNumFields(dmMesh, 2);CHECK_PETSC_ERROR(err);
    err = DMGetField(dmMesh, 0, &field);CHECK_PETSC_ERROR(err);
    err = PetscObjectCompose(field, "nearnullspace", (PetscObject) nullsp);CHECK_PETSC_ERROR(err);
    err = MatNullSpaceDestroy(&nullsp);CHECK_PETSC_ERROR(err);
    err = PetscObjectCompose(field, "pmat", (PetscObject) _jacobianPCFault);CHECK_PETSC_ERROR(err);
  } // if/else
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
