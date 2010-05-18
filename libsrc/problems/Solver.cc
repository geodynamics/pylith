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

#include <portinfo>

#include "Solver.hh" // implementation of class methods

#include "Formulation.hh" // USES Formulation

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

#include <cassert> // USES assert()

#define FIELD_SPLIT

#if defined(FIELD_SPLIT)
#include <petscmesh_solvers.hh> // USES constructFieldSplit()
#endif

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Solver::Solver(void) :
  _formulation(0),
  _logger(0)
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
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::Solver::initialize(const topology::SolutionFields& fields,
				     const topology::Jacobian& jacobian,
				     Formulation* const formulation)
{ // initialize
  assert(0 != formulation);

  _formulation = formulation;
} // initialize

// ----------------------------------------------------------------------
// Setup preconditioner for preconditioning with split fields.
void
pylith::problems::Solver::_setupFieldSplit(PetscPC* const pc,
					   Formulation* const formulation,
					   const topology::SolutionFields& fields)
{ // _setupFieldSplit

  PetscErrorCode err = 0;

  const topology::Field<topology::Mesh>& solution = fields.solution();
  const ALE::Obj<topology::Mesh::SieveMesh>& sieveMesh = 
    fields.mesh().sieveMesh();
  err = PCSetType(*pc, PCFIELDSPLIT); CHECK_PETSC_ERROR(err);
  err = PCSetOptionsPrefix(*pc, "fs_"); CHECK_PETSC_ERROR(err);
  err = PCSetFromOptions(*pc); CHECK_PETSC_ERROR(err);

#if defined(FIELD_SPLIT)
  constructFieldSplit(solution.section(), 
	      sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
						      solution.section()), 
		      solution.vector(), *pc);

  if (formulation->useCustomFaultPC()) {
    PetscKSP *ksps = 0;
    PetscMat A, P;
    PetscInt num, m, n;

    err = PCFieldSplitGetSubKSP(*pc, &num, &ksps); CHECK_PETSC_ERROR(err);

    // Put in PC matrix for fault
    MatStructure flag;
    err = KSPGetOperators(ksps[num-1], &A, 
			  PETSC_NULL, &flag); CHECK_PETSC_ERROR(err);
    
    err = PetscObjectReference((PetscObject) A); CHECK_PETSC_ERROR(err);
    err = MatGetLocalSize(A, &m, &n); CHECK_PETSC_ERROR(err);
    err = MatCreate(sieveMesh->comm(), &P); CHECK_PETSC_ERROR(err);
    err = MatSetSizes(P, m, n, 
		      PETSC_DECIDE, PETSC_DECIDE); CHECK_PETSC_ERROR(err);

    // Allocate just the diagonal. This should really go in the Fault
    // object that calculates the preconditioner since it is
    // implementation dependent.
    err = MatSeqAIJSetPreallocation(P, 1, PETSC_NULL); CHECK_PETSC_ERROR(err);
    err = MatMPIAIJSetPreallocation(P, 1, PETSC_NULL, 
				    0, PETSC_NULL); CHECK_PETSC_ERROR(err);
    
    err = MatSetFromOptions(P); CHECK_PETSC_ERROR(err);
    err = KSPSetOperators(ksps[num-1], A, P, flag); CHECK_PETSC_ERROR(err);
    
    err = PetscFree(ksps); CHECK_PETSC_ERROR(err);
  } // if

#endif
} // _setupFieldSplit


// End of file
