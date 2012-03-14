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

  const ALE::Obj<RealSection>& solutionSection = fields.solution().section();
  assert(!solutionSection.isNull());
  const ALE::Obj<SieveMesh>& sieveMesh = fields.solution().mesh().sieveMesh();
  assert(!sieveMesh.isNull());

  if (formulation->splitFields() && 
      formulation->useCustomConstraintPC() &&
      solutionSection->getNumSpaces() > sieveMesh->getDimension()) {
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

  PetscErrorCode err = 0;

  const ALE::Obj<SieveMesh>& sieveMesh = fields.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const topology::Field<topology::Mesh>& solution = fields.solution();
  const ALE::Obj<RealSection>& solutionSection = solution.section();
  assert(!solutionSection.isNull());

  err = PCSetType(*pc, PCFIELDSPLIT); CHECK_PETSC_ERROR(err);
  err = PCSetOptionsPrefix(*pc, "fs_"); CHECK_PETSC_ERROR(err);
  err = PCSetFromOptions(*pc); CHECK_PETSC_ERROR(err);

  constructFieldSplit(solutionSection, 
		      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
							      solutionSection), 
		      solution.vector(), *pc);

  const int spaceDim = sieveMesh->getDimension();
  if (formulation->splitFields() && 
      formulation->useCustomConstraintPC() &&
      solutionSection->getNumSpaces() > sieveMesh->getDimension()) {
    // We have split fields with a custom constraint preconditioner
    // and constraints exist.

    // Get total number of DOF associated with constraints field split
    const ALE::Obj<SieveMesh::order_type>& lagrangeGlobalOrder =
      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "faultDefault",
                                              solutionSection, spaceDim);
    assert(!lagrangeGlobalOrder.isNull());

    err = MatDestroy(&_jacobianPCFault); CHECK_PETSC_ERROR(err);
    PylithInt nrows = lagrangeGlobalOrder->getLocalSize();
    PylithInt ncols = nrows;

    err = MatCreate(sieveMesh->comm(), &_jacobianPCFault); CHECK_PETSC_ERROR(err);
    err = MatSetSizes(_jacobianPCFault, nrows, ncols, 
		      PETSC_DECIDE, PETSC_DECIDE); CHECK_PETSC_ERROR(err);
    err = MatSetType(_jacobianPCFault, MATAIJ);
    err = MatSetFromOptions(_jacobianPCFault); CHECK_PETSC_ERROR(err);
    
    // Allocate just the diagonal.
    err = MatSeqAIJSetPreallocation(_jacobianPCFault, 1, 
				    PETSC_NULL); CHECK_PETSC_ERROR(err);
    err = MatMPIAIJSetPreallocation(_jacobianPCFault, 1, PETSC_NULL, 
				    0, PETSC_NULL); CHECK_PETSC_ERROR(err);
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
      _ctx.faultFieldName = "2";
      break;
    case 3 :
      _ctx.faultFieldName = "3";
      break;
    default:
      assert(0);
      throw std::logic_error("Unknown space dimension in "
			     "Problems::_setupFieldSplit().");
    } // switch
    _ctx.faultA = _jacobianPCFault;
  } // if
} // _setupFieldSplit


// End of file
