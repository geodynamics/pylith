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
  const int spaceDim = sieveMesh->getDimension();
  const int numSpaces = solutionSection->getNumSpaces();
  const bool separateComponents = formulation->splitFieldComponents();

  err = PCSetType(*pc, PCFIELDSPLIT); CHECK_PETSC_ERROR(err);
  err = PCSetOptionsPrefix(*pc, "fs_"); CHECK_PETSC_ERROR(err);
  err = PCSetFromOptions(*pc); CHECK_PETSC_ERROR(err);

  if (formulation->splitFields() && 
      formulation->useCustomConstraintPC() &&
      numSpaces > spaceDim) {
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
    PetscMat* precon = new PetscMat[numSpaces];
    for (int i=0; i < numSpaces; ++i) {
      precon[i] = PETSC_NULL;
    } // for
    precon[numSpaces-1] = _jacobianPCFault;
    constructFieldSplit(solutionSection, PETSC_DETERMINE, PETSC_NULL, PETSC_NULL, 
			sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", solutionSection), precon, PETSC_NULL, solution.vector(), *pc);
    delete[] precon; precon = 0;
  } else {
    int numFields[2] = {spaceDim, (numSpaces > spaceDim) ? 1 : 0};
    MatNullSpace nullsp[2] = {PETSC_NULL, PETSC_NULL};
    PetscMat precon[2] = {PETSC_NULL, _jacobianPCFault};
    int* fields = new int[numSpaces];
    
    // Create rigid body null space.
    const ALE::Obj<RealSection>& coordinatesSection = sieveMesh->getRealSection("coordinates");
    assert(!coordinatesSection.isNull());
    const ALE::Obj<SieveMesh::label_sequence>& vertices = sieveMesh->depthStratum(0);
    assert(!vertices.isNull());
    const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
    const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
    PetscInt dim = spaceDim;
    PetscVec mode[6];

    if (dim > 1) {
      PetscInt n;
      err = VecGetLocalSize(solution.vector(), &n);CHECK_PETSC_ERROR(err);
      const int m = (dim * (dim + 1)) / 2;
      err = VecCreate(sieveMesh->comm(), &mode[0]);CHECK_PETSC_ERROR(err);
      err = VecSetSizes(mode[0], n, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
      err = VecSetUp(mode[0]);CHECK_PETSC_ERROR(err);
      for (int i = 1; i < m; ++i) {
	err = VecDuplicate(mode[0], &mode[i]);CHECK_PETSC_ERROR(err);
      } // for
      // :KLUDGE: Assume P1
      for (int d=0; d < dim; ++d) {
        PetscScalar values[3] = {0.0, 0.0, 0.0};
	
        values[d] = 1.0;
        solutionSection->zero();
        for(SieveMesh::label_sequence::iterator v_iter=verticesBegin; v_iter != verticesEnd; ++v_iter) {
          solutionSection->updatePoint(*v_iter, values);
        } // for
        solution.scatterSectionToVector();
        err = VecCopy(solution.vector(), mode[d]);CHECK_PETSC_ERROR(err);
      } // for
      for (int d = dim; d < m; ++d) {
        PetscInt k = (dim > 2) ? d - dim : d;
	
        solutionSection->zero();
        for (SieveMesh::label_sequence::iterator v_iter=verticesBegin; v_iter != verticesEnd; ++v_iter) {
          PetscScalar values[3] = {0.0, 0.0, 0.0};
          const PylithScalar* coords  = coordinatesSection->restrictPoint(*v_iter);

          for (int i=0; i < dim; ++i) {
            for (int j=0; j < dim; ++j) {
              values[j] += _epsilon(i, j, k)*coords[i];
            } // for
          } // for
          solutionSection->updatePoint(*v_iter, values);
        }
        solution.scatterSectionToVector();
        err = VecCopy(solution.vector(), mode[d]);CHECK_PETSC_ERROR(err);
      } // for
      for (int i=0; i < dim; ++i) {
	err = VecNormalize(mode[i], PETSC_NULL);CHECK_PETSC_ERROR(err);
      } // for
      /* Orthonormalize system */
      for (int i = dim; i < m; ++i) {
        PetscScalar dots[6];

        err = VecMDot(mode[i], i, mode, dots);CHECK_PETSC_ERROR(err);
        for(int j=0; j < i; ++j) dots[j] *= -1.0;
        err = VecMAXPY(mode[i], i, dots, mode);CHECK_PETSC_ERROR(err);
        err = VecNormalize(mode[i], PETSC_NULL);CHECK_PETSC_ERROR(err);
      } // for
      err = MatNullSpaceCreate(sieveMesh->comm(), PETSC_FALSE, m, mode, &nullsp[0]);CHECK_PETSC_ERROR(err);
      for(int i=0; i< m; ++i) {err = VecDestroy(&mode[i]);CHECK_PETSC_ERROR(err);}
    } // if

    for (int f=0; f < numSpaces; ++f) {
      fields[f] = f;
    } // for
    constructFieldSplit(solutionSection, 2, numFields, fields,
                        sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", solutionSection), precon, nullsp,
                        solution.vector(), *pc);
    err = MatNullSpaceDestroy(&nullsp[0]);CHECK_PETSC_ERROR(err);
    delete[] fields;
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
