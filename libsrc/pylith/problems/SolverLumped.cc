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

#include "SolverLumped.hh" // implementation of class methods

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/problems/Formulation.hh" // USES Formulation

#include "pylith/utils/EventLogger.hh" // USES EventLogger

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Constructor
pylith::problems::SolverLumped::SolverLumped(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::SolverLumped::~SolverLumped(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate data structures.
void
pylith::problems::SolverLumped::deallocate(void)
{ // deallocate
  Solver::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::SolverLumped::initialize(
				   const topology::SolutionFields& fields,
				   const topology::Field<topology::Mesh>& jacobian,
				   Formulation* formulation)
{ // initialize
  assert(0 != formulation);

  _initializeLogger();

  _formulation = formulation;
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverLumped::solve(
			      topology::Field<topology::Mesh>* solution,
			      const topology::Field<topology::Mesh>& jacobian,
			      const topology::Field<topology::Mesh>& residual)
{ // solve
  assert(0 != solution);
  assert(0 != _formulation);
  
  // solution = residual / jacobian
  
  const int setupEvent = _logger->eventId("SoLu setup");
  const int solveEvent = _logger->eventId("SoLu solve");
  const int adjustEvent = _logger->eventId("SoLu adjust");
  _logger->eventBegin(setupEvent);

  const spatialdata::geocoords::CoordSys* cs = solution->mesh().coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();
  
  // Get mesh vertices.
  DM             dmMesh = solution->mesh().dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  
  // Get sections.
  PetscSection solutionSection = solution->petscSection();
  Vec          solutionVec     = solution->localVector();
  PetscScalar *solutionArray;
  assert(solutionSection);assert(solutionVec);
	 
  PetscSection jacobianSection = jacobian.petscSection();
  Vec          jacobianVec     = jacobian.localVector();
  PetscScalar *jacobianArray;
  assert(jacobianSection);assert(jacobianVec);

  PetscSection residualSection = residual.petscSection();
  Vec          residualVec     = residual.localVector();
  PetscScalar *residualArray;
  assert(residualSection);assert(residualVec);
  
  _logger->eventEnd(setupEvent);
  _logger->eventBegin(solveEvent);

  err = VecGetArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(solutionVec, &solutionArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt jdof, joff;

    err = PetscSectionGetDof(jacobianSection, v, &jdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(jacobianSection, v, &joff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == jdof);
    PetscInt rdof, roff;

    err = PetscSectionGetDof(residualSection, v, &rdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(residualSection, v, &roff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == rdof);
    PetscInt sdof, soff;

    err = PetscSectionGetDof(solutionSection, v, &sdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(solutionSection, v, &soff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == sdof);

    for (int i=0; i < spaceDim; ++i) {
      assert(jacobianArray[joff+i] != 0.0);
      solutionArray[soff+i] = residualArray[roff+i] / jacobianArray[joff+i];
    } // for
  } // for
  PetscLogFlops((vEnd - vStart) * spaceDim);
  err = VecRestoreArray(jacobianVec, &jacobianArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(solutionVec, &solutionArray);CHECK_PETSC_ERROR(err);
  _logger->eventEnd(solveEvent);
  _logger->eventBegin(adjustEvent);

  // Update rate fields to be consistent with current solution.
  _formulation->calcRateFields();

  // Adjust solution to match constraints
  _formulation->adjustSolnLumped();

  // Update rate fields to be consistent with adjusted solution.
  _formulation->calcRateFields(); // :KLUDGE: Limit to only those changed?

  _logger->eventEnd(adjustEvent);
} // solve

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::problems::SolverLumped::_initializeLogger(void)
{ // initializeLogger
  delete _logger; _logger = new utils::EventLogger;
  assert(0 != _logger);
  _logger->className("SolverLumped");
  _logger->initialize();
  _logger->registerEvent("SoLu setup");
  _logger->registerEvent("SoLu solve");
  _logger->registerEvent("SoLu adjust");
} // initializeLogger


// End of file
