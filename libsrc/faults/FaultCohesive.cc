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

#include "CohesiveTopology.hh" // USES CohesiveTopology::create()

#include "pylith/meshio/UCDFaultFile.hh" // USES UCDFaultFile::read()
#include "pylith/utils/array.hh" // USES double_array

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
} // destructor

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
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesive::adjustTopology(const topology::Mesh& mesh,
					      const bool flipFault)
{ // adjustTopology
  assert(std::string("") != label());
  SubMesh faultMesh;
  Mesh fault
  Obj<SubMesh>   faultMesh = NULL;
  Obj<ALE::Mesh> faultBd   = NULL;

  if (_useFaultMesh) {
    const int faultDim = 2;
    assert(3 == mesh.dimension());

    //MPI_Bcast(&faultDim, 1, MPI_INT, 0, comm);
    faultMesh = new SubMesh(mesh->comm(), faultDim, mesh->debug());
    meshio::UCDFaultFile::readFault(_faultMeshFilename, mesh, 
				    faultMesh, faultBd);

    // Get group of vertices associated with fault
    const ALE::Obj<IntSection>& groupField = sieveMesh->getIntSection(label());
    faultMesh->setRealSection("coordinates", 
			      mesh->getRealSection("coordinates"));

    CohesiveTopology::create(faultMesh, faultBd, mesh, groupField, id(),
			     _useLagrangeConstraints());
  } else {
    if (!sieveMesh->hasIntSection(label())) {
      std::ostringstream msg;
      msg << "Mesh missing group of vertices '" << label()
          << " for fault interface condition.";
      throw std::runtime_error(msg.str());
    } // if  

    // Get group of vertices associated with fault
    const ALE::Obj<IntSection>& groupField = sieveMesh->getIntSection(label());
    assert(!groupField.isNull());

    faultMesh = 
      new SubMesh(mesh->comm(), mesh->getDimension()-1, mesh->debug());

    CohesiveTopology::createFault(faultMesh, faultBd, mesh, groupField, 
				  flipFault);

    CohesiveTopology::create(faultMesh, faultBd, mesh, groupField, id(), 
			     _useLagrangeConstraints());
  } // if/else
} // adjustTopology


// End of file 
