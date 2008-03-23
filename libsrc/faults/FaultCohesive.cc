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
pylith::faults::FaultCohesive::FaultCohesive(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesive::~FaultCohesive(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesive::adjustTopology(const ALE::Obj<ALE::Mesh>& mesh)
{ // adjustTopology
  assert(std::string("") != label());

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
} // adjustTopology


// End of file 
