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

#include "pylith/topology/SubMesh.hh" // USES SubMesh

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
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::Fault::deallocate(void)
{ // deallocate
  delete _faultMesh; _faultMesh = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Get mesh associated with fault fields.
const pylith::topology::SubMesh&
pylith::faults::Fault:: faultMesh(void) const
{ // faultMesh
  return *_faultMesh;
} // faultMesh


// End of file 
