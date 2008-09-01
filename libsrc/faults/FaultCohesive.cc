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

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()
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
pylith::faults::FaultCohesive::adjustTopology(const ALE::Obj<Mesh>& mesh)
{ // adjustTopology
  assert(std::string("") != label());

  if (!_useFaultMesh) {
    // Use group of vertices to define fault.
    if (!mesh->hasIntSection(label())) {
      std::ostringstream msg;
      msg << "Mesh missing group of vertices '" << label()
	  << " for fault interface condition.";
      throw std::runtime_error(msg.str());
    } // if  
    
    // Get group of vertices associated with fault
    const ALE::Obj<int_section_type>& groupField = 
      mesh->getIntSection(label());
    assert(!groupField.isNull());
    
    ALE::Obj<Mesh> faultMesh;
    CohesiveTopology::create(&faultMesh, mesh, groupField, id(),
			     _useLagrangeConstraints());
  } else {
    // Use fault mesh to define fault.
    std::cout << "ADD FAULT MESH ADJUSTING TOPOLOGY STUFF HERE." << std::endl;
  } // else
} // adjustTopology


// End of file 
