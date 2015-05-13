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

#include "FaultCohesive.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/MeshOps.hh" // USES MeshOps

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
    PetscMPIInt rank;
    PetscErrorCode err;

    assert(std::string("") != label());
    // We do not have labels on all ranks until after distribution
    err = MPI_Comm_rank(PetscObjectComm((PetscObject) dmMesh), &rank);PYLITH_CHECK_ERROR(err);
    err = DMPlexHasLabel(dmMesh, label(), &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel && !rank) {
      std::ostringstream msg;
      msg << "Mesh missing group of vertices '" << label() << "' for fault interface condition.";
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
                                              int *firstFaultCell)
{ // adjustTopology
  PYLITH_METHOD_BEGIN;

  assert(mesh);
  assert(std::string("") != label());
  
  try {
    topology::Mesh faultMesh;
  
    // Get group of vertices associated with fault
    PetscDM dmMesh = mesh->dmMesh();assert(dmMesh);
    
    if (!_useFaultMesh) {
      const char* charlabel = label();

      PetscDMLabel   groupField;
      PetscBool      hasLabel;
      PetscInt       depth, gdepth, dim;
      PetscMPIInt    rank;
      PetscErrorCode err;
      // We do not have labels on all ranks until after distribution
      err = MPI_Comm_rank(PetscObjectComm((PetscObject) dmMesh), &rank);PYLITH_CHECK_ERROR(err);
      err = DMPlexHasLabel(dmMesh, charlabel, &hasLabel);PYLITH_CHECK_ERROR(err);
      if (!hasLabel && !rank) {
        std::ostringstream msg;
        msg << "Mesh missing group of vertices '" << label()
            << "' for fault interface condition.";
        throw std::runtime_error(msg.str());
      } // if
      err = DMGetDimension(dmMesh, &dim);PYLITH_CHECK_ERROR(err);
      err = DMPlexGetDepth(dmMesh, &depth);PYLITH_CHECK_ERROR(err);
      err = MPI_Allreduce(&depth, &gdepth, 1, MPIU_INT, MPI_MAX, mesh->comm());PYLITH_CHECK_ERROR(err);
      err = DMPlexGetLabel(dmMesh, charlabel, &groupField);PYLITH_CHECK_ERROR(err);
      CohesiveTopology::createFault(&faultMesh, *mesh, groupField);
      PetscDMLabel faultBdLabel = NULL;

      // We do not have labels on all ranks until after distribution
      if (strlen(edge()) > 0 && !rank) {
	err = DMPlexGetLabel(dmMesh, edge(), &faultBdLabel);PYLITH_CHECK_ERROR(err);
	if (!faultBdLabel) {
	  std::ostringstream msg;
	  msg << "Could not find nodeset/pset '" << edge() << "' marking buried edges for fault '" << label() << "'.";
	  throw std::runtime_error(msg.str());
	} // if
      } // if
      CohesiveTopology::create(mesh, faultMesh, faultBdLabel, id(), *firstFaultVertex, *firstLagrangeVertex, *firstFaultCell, useLagrangeConstraints());
    } else {
      assert(3 == mesh->dimension());
      throw std::logic_error("Support for UCD fault files no longer implemented."); 
    } // if/else

    // Check consistency of mesh.
    topology::MeshOps::checkTopology(*mesh);
    topology::MeshOps::checkTopology(faultMesh);

  } catch (const std::exception& err) {
    std::ostringstream msg;
    msg << "Error occurred while adjusting topology to create cohesive cells for fault '" << label() << "'.\n"
	<< err.what();
    throw std::runtime_error(msg.str());
  }


  PYLITH_METHOD_END;
} // adjustTopology


// End of file 
