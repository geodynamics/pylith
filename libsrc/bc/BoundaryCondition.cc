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
pylith::bc::BoundaryCondition::verifyConfiguration(const ALE::Obj<Mesh>& mesh) const
{ // verifyConfiguration
  assert(!mesh.isNull());

  if (!mesh->hasIntSection(_label)) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << _label
	<< " for boundary condition.";
    throw std::runtime_error(msg.str());
  } // if
} // verifyConfiguration


// End of file 
