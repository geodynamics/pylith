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
// Copyright (c) 2010 University of California, Davis
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

#include <cassert> // USES assert()

#define FIELD_SPLIT

#if defined(FIELD_SPLIT)
#include <petscdmmesh_solvers.hh> // USES constructFieldSplit()
#endif

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

EXTERN_C_BEGIN
PetscErrorCode  MyMatGetSubMatrix(Mat mat, IS isrow, IS iscol, MatReuse reuse, Mat *newmat) {
  FaultPreconCtx *ctx;
  IS              faultIS;
  PetscBool       isFaultRow, isFaultCol;
  PetscErrorCode  ierr;

  ierr = MatShellGetContext(mat, (void **) &ctx);CHKERRQ(ierr);
  ierr = PCFieldSplitGetIS(ctx->pc, ctx->faultFieldName, &faultIS);CHKERRQ(ierr);
  ierr = ISEqual(isrow, faultIS, &isFaultRow);CHKERRQ(ierr);
  ierr = ISEqual(iscol, faultIS, &isFaultCol);CHKERRQ(ierr);
  if (isFaultRow && isFaultCol) {
    if (reuse == MAT_INITIAL_MATRIX) {
      ierr = PetscObjectReference((PetscObject) ctx->faultA);CHKERRQ(ierr);
      *newmat = ctx->faultA;
    }
  } else {
    ierr = MatGetSubMatrix(ctx->A, isrow, iscol, reuse, newmat);CHKERRQ(ierr);
  }
}
EXTERN_C_END

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Solver::Solver(void) :
  _formulation(0),
  _logger(0),
  _precondMatrix(0),
  _jacobianPre(0)
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
  if (0 != _jacobianPre) {
    PetscErrorCode err = MatDestroy(_jacobianPre); _jacobianPre = 0;
    CHECK_PETSC_ERROR(err);
  } // if
  if (0 != _precondMatrix) {
    PetscErrorCode err = MatDestroy(_precondMatrix); _precondMatrix = 0;
    CHECK_PETSC_ERROR(err);
  } // if
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::Solver::initialize(const topology::SolutionFields& fields,
				     const topology::Jacobian& jacobian,
				     Formulation* const formulation)
{ // initialize
  assert(0 != formulation);

  PetscMat jacobianMat = jacobian.matrix();
  _jacobianPre = jacobianMat;
  // Make global preconditioner matrix
  const ALE::Obj<RealSection>& solutionSection = fields.solution().section();
  assert(!solutionSection.isNull());
  const ALE::Obj<SieveMesh>& sieveMesh = fields.solution().mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  if (solutionSection->getNumSpaces() > sieveMesh->getDimension() && 0 != _precondMatrix) {
    PetscInt M, N, m, n;
    PetscErrorCode err;

    err = MatGetSize(jacobianMat, &M, &N);CHECK_PETSC_ERROR(err);
    err = MatGetLocalSize(jacobianMat, &m, &n);CHECK_PETSC_ERROR(err);
    err = MatCreateShell(fields.mesh().comm(), m, n, M, N, &_ctx, &_jacobianPre);CHECK_PETSC_ERROR(err);
    err = MatShellSetOperation(_jacobianPre, MATOP_GET_SUBMATRIX, (void (*)(void)) MyMatGetSubMatrix);
    _ctx.A              = jacobianMat;
    _ctx.faultFieldName = "3";
    _ctx.faultA         = _precondMatrix;
  }

  _formulation = formulation;
} // initialize

// ----------------------------------------------------------------------
// Setup preconditioner for preconditioning with split fields.
void
pylith::problems::Solver::_setupFieldSplit(PetscPC* const pc,
					   PetscMat* const precondMatrix,
					   Formulation* const formulation,
					   const topology::SolutionFields& fields)
{ // _setupFieldSplit
  assert(0 != pc);
  assert(0 != precondMatrix);
  assert(0 != formulation);

  PetscErrorCode err = 0;

  const ALE::Obj<SieveMesh>& sieveMesh = fields.mesh().sieveMesh();
  const topology::Field<topology::Mesh>& solution = fields.solution();
  const ALE::Obj<RealSection>& solutionSection = solution.section();
  assert(!solutionSection.isNull());

  err = PCSetType(*pc, PCFIELDSPLIT); CHECK_PETSC_ERROR(err);
  err = PCSetOptionsPrefix(*pc, "fs_"); CHECK_PETSC_ERROR(err);
  err = PCSetFromOptions(*pc); CHECK_PETSC_ERROR(err);

#if defined(FIELD_SPLIT)
  constructFieldSplit(solutionSection, 
	      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
						      solutionSection), 
		      solution.vector(), *pc);

  const int spaceDim = sieveMesh->getDimension();
  if (solutionSection->getNumSpaces() > spaceDim &&
      formulation->useCustomConstraintPC()) {
    // Get total number of DOF associated with constraints field split
    const ALE::Obj<SieveMesh::order_type>& lagrangeGlobalOrder =
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "faultDefault",
                                              solutionSection, spaceDim);
    assert(!lagrangeGlobalOrder.isNull());

    if (0 != *precondMatrix) {
      err = MatDestroy(*precondMatrix); *precondMatrix = 0;
      CHECK_PETSC_ERROR(err);
    } // if
    PetscInt nrows = lagrangeGlobalOrder->getLocalSize();
    PetscInt ncols = nrows;

    err = MatCreate(sieveMesh->comm(), precondMatrix); CHECK_PETSC_ERROR(err);
    err = MatSetSizes(*precondMatrix, nrows, ncols, 
		      PETSC_DECIDE, PETSC_DECIDE); CHECK_PETSC_ERROR(err);
    err = MatSetType(*precondMatrix, MATAIJ);
    err = MatSetFromOptions(*precondMatrix); CHECK_PETSC_ERROR(err);
    
#if 1
    // Allocate just the diagonal.
    err = MatSeqAIJSetPreallocation(*precondMatrix, 1, 
				    PETSC_NULL); CHECK_PETSC_ERROR(err);
    err = MatMPIAIJSetPreallocation(*precondMatrix, 1, PETSC_NULL, 
				    0, PETSC_NULL); CHECK_PETSC_ERROR(err);
#else
    // Allocate full matrix (overestimate).
    err = MatSeqAIJSetPreallocation(*precondMatrix, ncols, 
				    PETSC_NULL); CHECK_PETSC_ERROR(err);
    err = MatMPIAIJSetPreallocation(*precondMatrix, ncols, PETSC_NULL, 
				    0, PETSC_NULL); CHECK_PETSC_ERROR(err);
#endif
    
    // Set preconditioning matrix in formulation
    formulation->customPCMatrix(*precondMatrix);
  } // if

#endif
} // _setupFieldSplit


// End of file
