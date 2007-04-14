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
  _label("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::Fault::~Fault(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::faults::Fault::Fault(const Fault& f) :
  _id(f._id),
  _label(f._label)
{ // copy constructor
} // copy constructor


// End of file 
