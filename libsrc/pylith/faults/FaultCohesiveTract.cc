// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesiveTract.hh" // implementation of object methods
#include "CohesiveTopology.hh" // USES CohesiveTopology

#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES StratumIS

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveTract::FaultCohesiveTract(void)
{ // constructor
  _useLagrangeConstraints = false;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveTract::~FaultCohesiveTract(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::FaultCohesiveTract::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  FaultCohesive::deallocate();

  PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveTract::initialize(const topology::Mesh& mesh,
					       const PylithScalar upDir[3])
{ // initialize
  PYLITH_METHOD_BEGIN;
  
  assert(upDir);
  assert(_quadrature);

  delete _faultMesh; _faultMesh = new topology::Mesh();
  CohesiveTopology::createFaultParallel(_faultMesh, mesh, id(), label(), useLagrangeConstraints());

  // Reset fields.
  delete _fields; _fields = new topology::Fields(*_faultMesh);assert(_fields);

  // Initialize quadrature geometry.
  _quadrature->initializeGeometry();

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveTract::integrateResidual(const topology::Field& residual,
						      const PylithScalar t,
						      topology::SolutionFields* const fields)
{ // integrateResidual
  throw std::logic_error("FaultCohesiveTract::integrateResidual() not implemented.");
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveTract::integrateJacobian(topology::Jacobian* jacobian,
						      const PylithScalar t,
						      topology::SolutionFields* const fields)
{ // integrateJacobian
  throw std::logic_error("FaultCohesiveTract::integrateJacobian() not implemented.");

  _needNewJacobian = false;
} // integrateJacobian
  
// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveTract::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;
  
  assert(_quadrature);

  const PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  PetscBool hasLabel = PETSC_FALSE;
  PetscErrorCode err = DMHasLabel(dmMesh, label(), &hasLabel);PYLITH_CHECK_ERROR(err);
  if (!hasLabel) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << label()
	<< " for boundary condition.";
    throw std::runtime_error(msg.str());
  } // if  

  // check compatibility of mesh and quadrature scheme
  const int dimension = mesh.dimension()-1;
  if (_quadrature->cellDim() != dimension) {
    std::ostringstream msg;
    msg << "Dimension of reference cell in quadrature scheme (" 
	<< _quadrature->cellDim() 
	<< ") does not match dimension of cells in mesh (" 
	<< dimension << ") for fault '" << label()
	<< "'.";
    throw std::runtime_error(msg.str());
  } // if

  const int numCorners = _quadrature->refGeometry().numCorners();
  const bool includeOnlyCells = true;
  topology::StratumIS cohesiveIS(dmMesh, "material-id", id(), includeOnlyCells);
  const PetscInt* cells = cohesiveIS.points();
  const PetscInt ncells = cohesiveIS.size();

  PetscInt coneSize = 0;
  for (PetscInt i=0; i < ncells; ++i) {
    err = DMPlexGetConeSize(dmMesh, cells[i], &coneSize);PYLITH_CHECK_ERROR(err);
    // TODO: Should be changed to Closure()
    if (2*numCorners != coneSize) {
      // No Lagrange vertices, just negative and positive sides of the
      // fault, so coneSize is 2*numCorners.
      std::ostringstream msg;
      msg << "Number of vertices in reference cell (" << numCorners 
	  << ") is not compatible with number of vertices (" << coneSize
	  << ") in cohesive cell " << cells[i] << " for fault '"
	  << label() << "'.";
      throw std::runtime_error(msg.str());
    } // if
  } // for

  PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field&
pylith::faults::FaultCohesiveTract::vertexField(const char* name,
						const topology::SolutionFields* fields)
{ // vertexField
  throw std::logic_error("FaultCohesiveTract::vertexField() not implemented.");
} // vertexField

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const pylith::topology::Field&
pylith::faults::FaultCohesiveTract::cellField(const char* name,
					      const topology::SolutionFields* fields)
{ // cellField
  throw std::logic_error("FaultCohesiveTract::cellField() not implemented.");
} // cellField


// End of file 
