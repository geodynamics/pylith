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

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesive::FaultCohesive(void) :
  _faultMesh(new ALE::Obj<ALE::Mesh>)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesive::~FaultCohesive(void)
{ // destructor
  delete _faultMesh; _faultMesh = 0;
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::faults::FaultCohesive::FaultCohesive(const FaultCohesive& f) :
  Fault(f),
  _faultMesh(new ALE::Obj<ALE::Mesh>)
{ // copy constructor
  *_faultMesh = *f._faultMesh;
} // copy constructor

// ----------------------------------------------------------------------
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesive::adjustTopology(ALE::Obj<ALE::Mesh>* mesh) const
{ // adjustTopology
  assert(0 != mesh);
  assert("" != label());

  // Get group of vertices associated with fault
  const ALE::Obj<int_section_type>& groupField = 
    (*mesh)->getIntSection(label());
  assert(!groupField.isNull());

  CohesiveTopology::create(_faultMesh, *mesh, groupField, id());
} // adjustTopology


// End of file 
