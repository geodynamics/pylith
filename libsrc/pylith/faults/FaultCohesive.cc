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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesive.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/meshio/UCDFaultFile.hh" // USES UCDFaultFile
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesive::FaultCohesive(void) :
  _fields(0),
  _useFaultMesh(false),
  _faultMeshFilename("fault.inp")
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
  Fault::deallocate();

  delete _fields; _fields = 0;
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
// Set filename of UCD file for fault mesh.
void
pylith::faults::FaultCohesive::faultMeshFilename(const char* filename)
{ // faultMeshFilename
  _faultMeshFilename = filename;
} // faultMeshFilename

// ----------------------------------------------------------------------
// Get number of vertices in fault.
int
pylith::faults::FaultCohesive::numVertices(const topology::Mesh& mesh) const
{ // numVertices
  int nvertices = 0;

  if (!_useFaultMesh) {
    // Get group of vertices associated with fault
    const ALE::Obj<topology::Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());
    if (!sieveMesh->hasIntSection(label())) {
      std::ostringstream msg;
      msg << "Mesh missing group of vertices '" << label()
          << "' for fault interface condition.";
      throw std::runtime_error(msg.str());
    } // if  

    assert(std::string("") != label());
    const ALE::Obj<topology::Mesh::IntSection>& groupField = 
      mesh.sieveMesh()->getIntSection(label());
    nvertices = groupField->size();
  } else {
    assert(3 == mesh.dimension());
    nvertices = meshio::UCDFaultFile::numVertices(_faultMeshFilename.c_str());
  } // else

  return nvertices;
} // numVertices

// ----------------------------------------------------------------------
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesive::adjustTopology(topology::Mesh* const mesh,
                                              int *firstFaultVertex,
                                              int *firstLagrangeVertex,
                                              int *firstFaultCell,
                                              const bool flipFault)
{ // adjustTopology
  assert(0 != mesh);
  assert(std::string("") != label());
  
  topology::SubMesh faultMesh;
  ALE::Obj<SieveFlexMesh> faultBoundary;
  
  // Get group of vertices associated with fault
  const ALE::Obj<topology::Mesh::SieveMesh>& sieveMesh = mesh->sieveMesh();
  assert(!sieveMesh.isNull());

  if (!_useFaultMesh) {
    if (!sieveMesh->hasIntSection(label())) {
      std::ostringstream msg;
      msg << "Mesh missing group of vertices '" << label()
          << "' for fault interface condition.";
      throw std::runtime_error(msg.str());
    } // if  
    const ALE::Obj<topology::Mesh::IntSection>& groupField = 
      sieveMesh->getIntSection(label());
    assert(!groupField.isNull());
    CohesiveTopology::createFault(&faultMesh, faultBoundary, *mesh, groupField, 
				  flipFault);

    CohesiveTopology::create(mesh, faultMesh, faultBoundary, groupField, id(), 
                             *firstFaultVertex, *firstLagrangeVertex, *firstFaultCell, useLagrangeConstraints());

  } else {
    //std::cout << "BEFORE ADJUSTING TOPOLOGY FOR FAULT '" << label() << "' firstFaultVertex: " << *firstFaultVertex << ", firstFaultCell: " << *firstFaultCell << std::endl;

    const int faultDim = 2;
    assert(3 == mesh->dimension());

    meshio::UCDFaultFile::read(_faultMeshFilename.c_str(),
			       &faultMesh, faultBoundary, *mesh);

    // Set coordinates in fault mesh
    const ALE::Obj<topology::SubMesh::SieveMesh>& faultSieveMesh = 
      faultMesh.sieveMesh();
    assert(!faultSieveMesh.isNull());
    faultSieveMesh->setRealSection("coordinates", 
				   sieveMesh->getRealSection("coordinates"));

    const ALE::Obj<topology::Mesh::IntSection>& groupField = 
      sieveMesh->getIntSection(label());
    assert(!groupField.isNull());
    CohesiveTopology::create(mesh, faultMesh, faultBoundary, groupField, id(),
                             *firstFaultVertex, *firstLagrangeVertex, *firstFaultCell, useLagrangeConstraints());

    //std::cout << "AFTER ADJUSTING TOPOLOGY FOR FAULT '" << label() << "' firstFaultVertex: " << *firstFaultVertex << ", firstFaultCell: " << *firstFaultCell << std::endl;
  } // if/else
} // adjustTopology


// End of file 
