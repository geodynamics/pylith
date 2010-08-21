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
  const ALE::Obj<SieveMesh>& sieveMesh = solution->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = 
    vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  
  // Get sections.
  double_array solutionVertex(spaceDim);
  const ALE::Obj<RealSection>& solutionSection = solution->section();
  assert(!solutionSection.isNull());
	 
  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  assert(!jacobianSection.isNull());

  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());
  
  _logger->eventEnd(setupEvent);
  _logger->eventBegin(solveEvent);

  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin; 
       v_iter != verticesEnd;
       ++v_iter) {
    assert(spaceDim == jacobianSection->getFiberDimension(*v_iter));
    const double* jacobianVertex = jacobianSection->restrictPoint(*v_iter);

    assert(spaceDim == residualSection->getFiberDimension(*v_iter));
    const double* residualVertex = residualSection->restrictPoint(*v_iter);

    for (int i=0; i < spaceDim; ++i) {
      assert(jacobianVertex[i] != 0.0);
      solutionVertex[i] = residualVertex[i] / jacobianVertex[i];
    } // for
    
    assert(solutionSection->getFiberDimension(*v_iter) == spaceDim);
    solutionSection->updatePoint(*v_iter, &solutionVertex[0]);
  } // for
  PetscLogFlops(vertices->size() * spaceDim);
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
