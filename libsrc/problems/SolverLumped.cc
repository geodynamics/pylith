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
	 
  double_array jacobianVertex(spaceDim);
  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  assert(!jacobianSection.isNull());

  double_array residualVertex(spaceDim);
  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());
  
  _logger->eventEnd(setupEvent);
  _logger->eventBegin(solveEvent);

  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin; 
       v_iter != verticesEnd;
       ++v_iter) {
    jacobianSection->restrictPoint(*v_iter, &jacobianVertex[0],
				   jacobianVertex.size());
    residualSection->restrictPoint(*v_iter, &residualVertex[0],
				   residualVertex.size());

#if !defined(NDEBUG)
    for (int i=0; i < spaceDim; ++i)
      assert(jacobianVertex[i] != 0.0);
#endif

    solutionVertex = residualVertex / jacobianVertex;
    
    assert(solutionSection->getFiberDimension(*v_iter) == spaceDim);
    solutionSection->updatePoint(*v_iter, &solutionVertex[0]);
  } // for
  PetscLogFlops(vertices->size() * spaceDim);
  _logger->eventEnd(solveEvent);
  _logger->eventBegin(adjustEvent);

  // Adjust solution to match constraints
  _formulation->adjustSolnLumped();

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
