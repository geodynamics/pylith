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

#include "BoundaryCondition.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh

#include <stdexcept> // USES std::runtime_error()

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::BoundaryCondition::BoundaryCondition(void) :
  _label(""),
  _db(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::BoundaryCondition::~BoundaryCondition(void)
{ // destructor
  _db = 0;
} // destructor

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::BoundaryCondition::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  const ALE::Obj<topology::Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!sieveMesh->hasIntSection(_label)) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << _label
	<< " for boundary condition.";
    throw std::runtime_error(msg.str());
  } // if
} // verifyConfiguration


// End of file 
