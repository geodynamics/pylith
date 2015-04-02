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

#include "SolverLumped.hh" // implementation of class methods

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/problems/Formulation.hh" // USES Formulation

#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

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
  PYLITH_METHOD_BEGIN;

  Solver::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::SolverLumped::initialize(const topology::SolutionFields& fields,
					   const topology::Field& jacobian,
					   Formulation* formulation)
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(formulation);

  _initializeLogger();

  _formulation = formulation;

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverLumped::solve(topology::Field* solution,
				      const topology::Field& jacobian,
				      const topology::Field& residual)
{ // solve
  PYLITH_METHOD_BEGIN;

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
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  
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

  PYLITH_METHOD_END;
} // solve

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::problems::SolverLumped::_initializeLogger(void)
{ // initializeLogger
  PYLITH_METHOD_BEGIN;

  delete _logger; _logger = new utils::EventLogger;assert(_logger);
  _logger->className("SolverLumped");
  _logger->initialize();
  _logger->registerEvent("SoLu setup");
  _logger->registerEvent("SoLu solve");
  _logger->registerEvent("SoLu adjust");

  PYLITH_METHOD_END;
} // initializeLogger


// End of file
