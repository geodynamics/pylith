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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesive.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesive::FaultCohesive(void) :
  _fields(0),
  _useFaultMesh(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesive::~FaultCohesive(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::FaultCohesive::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  Fault::deallocate();
  feassemble::Integrator::deallocate();

  delete _fields; _fields = 0;

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set flag for using fault mesh or group of vertices to define
// fault surface.
void
pylith::faults::FaultCohesive::useFaultMesh(const bool flag)
{ // useFaultMesh
  _useFaultMesh = flag;
} // useFaultMesh

// ----------------------------------------------------------------------
// Get the number of vertices associated with the fault (before
// fault mesh exists).
int
pylith::faults::FaultCohesive::numVerticesNoMesh(const topology::Mesh& mesh) const
{ // numVerticesNoMesh
  PYLITH_METHOD_BEGIN;

  PetscInt nvertices = 0;

  if (!_useFaultMesh) {
    // Get group of vertices associated with fault
    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
    PetscBool hasLabel = PETSC_FALSE;
    PetscErrorCode err;

    assert(std::string("") != label());
    err = DMPlexHasLabel(dmMesh, label(), &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
      std::ostringstream msg;
      msg << "Mesh missing group of vertices '" << label()
          << "' for fault interface condition.";
      throw std::runtime_error(msg.str());
    } // if  
    err = DMPlexGetStratumSize(dmMesh, label(), 1, &nvertices);PYLITH_CHECK_ERROR(err);
  } else {
    assert(3 == mesh.dimension());
    nvertices = -1;
  } // else

  PYLITH_METHOD_RETURN(nvertices);
} // numVerticesNoMesh

// ----------------------------------------------------------------------
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesive::adjustTopology(topology::Mesh* const mesh,
                                              int *firstFaultVertex,
                                              int *firstLagrangeVertex,
                                              int *firstFaultCell,
                                              const bool flipFault)
{ // adjustTopology
  PYLITH_METHOD_BEGIN;

  assert(mesh);
  assert(std::string("") != label());
  
  try {
    topology::Mesh faultMesh;
    PetscDM faultBoundary = NULL;
  
    // Get group of vertices associated with fault
    PetscDM dmMesh = mesh->dmMesh();assert(dmMesh);
    
    if (!_useFaultMesh) {
      const char* charlabel = label();

      PetscDMLabel groupField;
      PetscBool hasLabel;
      PetscErrorCode err;
      err = DMPlexHasLabel(dmMesh, charlabel, &hasLabel);PYLITH_CHECK_ERROR(err);
      if (!hasLabel) {
        std::ostringstream msg;
        msg << "Mesh missing group of vertices '" << label()
            << "' for fault interface condition.";
        throw std::runtime_error(msg.str());
      } // if
      err = DMPlexGetLabel(dmMesh, charlabel, &groupField);PYLITH_CHECK_ERROR(err);
      CohesiveTopology::createFault(&faultMesh, faultBoundary, *mesh, groupField, flipFault);
      
      CohesiveTopology::create(mesh, faultMesh, faultBoundary, groupField, id(), *firstFaultVertex, *firstLagrangeVertex, *firstFaultCell, useLagrangeConstraints());
      err = DMDestroy(&faultBoundary);PYLITH_CHECK_ERROR(err);
    } else {
      const int faultDim = 2;
      assert(3 == mesh->dimension());
      throw std::logic_error("Support for UCD fault files no longer implemented."); 
    } // if/else
  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error occurred while adjusting topology to create cohesive cells for fault '" << label() << "'.\n"
	<< err.what();
    throw std::runtime_error(msg.str());
  }

  PYLITH_METHOD_END;
} // adjustTopology


// End of file 
