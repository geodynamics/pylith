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

#include "Fault.hh" // implementation of object methods

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::Fault::Fault(void) :
  _id(0),
  _label(""),
  _faultMesh(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::Fault::~Fault(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get mesh associated with fault fields.
const ALE::Obj<pylith::SubMesh>&
pylith::faults::Fault:: faultMesh(void) const
{ // faultMesh
  return _faultMesh;
} // faultMesh


// End of file 
