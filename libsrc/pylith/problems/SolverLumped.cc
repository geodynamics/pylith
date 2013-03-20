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
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/problems/Formulation.hh" // USES Formulation

#include "pylith/utils/EventLogger.hh" // USES EventLogger

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
pylith::problems::SolverLumped::initialize(const topology::SolutionFields& fields,
					   const topology::Field<topology::Mesh>& jacobian,
					   Formulation* formulation)
{ // initialize
  assert(formulation);

  _initializeLogger();

  _formulation = formulation;
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverLumped::solve(topology::Field<topology::Mesh>* solution,
				      const topology::Field<topology::Mesh>& jacobian,
				      const topology::Field<topology::Mesh>& residual)
{ // solve
  assert(solution);
  assert(_formulation);
  
  // solution = residual / jacobian
  
  const int setupEvent = _logger->eventId("SoLu setup");
  const int solveEvent = _logger->eventId("SoLu solve");
  const int adjustEvent = _logger->eventId("SoLu adjust");
  _logger->eventBegin(setupEvent);

  const spatialdata::geocoords::CoordSys* cs = solution->mesh().coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();
  
  // Get mesh vertices.
  PetscDM dmMesh = solution->mesh().dmMesh(); assert(dmMesh);
  topology::Stratum depthStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = depthStratum.begin();
  const PetscInt vEnd = depthStratum.end();
  
  // Get sections.
  topology::VecVisitorMesh solutionVisitor(*solution);
  PetscScalar* solutionArray = solutionVisitor.localArray();

  topology::VecVisitorMesh jacobianVisitor(jacobian);
  PetscScalar* jacobianArray = jacobianVisitor.localArray();

  topology::VecVisitorMesh residualVisitor(residual);
  PetscScalar* residualArray = residualVisitor.localArray();

  _logger->eventEnd(setupEvent);
  _logger->eventBegin(solveEvent);

  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt joff = jacobianVisitor.sectionOffset(v);
    assert(spaceDim == jacobianVisitor.sectionDof(v));

    const PetscInt roff = residualVisitor.sectionOffset(v);
    assert(spaceDim == residualVisitor.sectionDof(v));

    const PetscInt soff = solutionVisitor.sectionOffset(v);
    assert(spaceDim == solutionVisitor.sectionDof(v));

    for (int i=0; i < spaceDim; ++i) {
      assert(jacobianArray[joff+i] != 0.0);
      solutionArray[soff+i] = residualArray[roff+i] / jacobianArray[joff+i];
    } // for
  } // for
  PetscLogFlops((vEnd - vStart) * spaceDim);
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
  assert(_logger);
  _logger->className("SolverLumped");
  _logger->initialize();
  _logger->registerEvent("SoLu setup");
  _logger->registerEvent("SoLu solve");
  _logger->registerEvent("SoLu adjust");
} // initializeLogger


// End of file
