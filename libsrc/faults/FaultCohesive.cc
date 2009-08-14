// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesive.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/meshio/UCDFaultFile.hh" // USES UCDFaultFile

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesive::FaultCohesive(void) :
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
                                              int *firstFaultCell,
                                              const bool flipFault)
{ // adjustTopology
  assert(0 != mesh);
  assert(std::string("") != label());
  
  topology::SubMesh faultMesh;
  ALE::Obj<ALE::Mesh> faultBoundary;
  
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
			     *firstFaultVertex, *firstFaultCell, useLagrangeConstraints());

  } else {
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

    if (!sieveMesh->hasIntSection(label())) {
      std::ostringstream msg;
      msg << "Mesh missing group of vertices '" << label()
          << "' for fault interface condition.";
      throw std::runtime_error(msg.str());
    } // if  
    const ALE::Obj<topology::Mesh::IntSection>& groupField = 
      sieveMesh->getIntSection(label());
    assert(!groupField.isNull());
    CohesiveTopology::create(mesh, faultMesh, faultBoundary, groupField, id(),
                  *firstFaultVertex, *firstFaultCell, useLagrangeConstraints());
  } // if/else
} // adjustTopology


// End of file 
