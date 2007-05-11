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

#include "SlipTimeFn.hh" // implementation of object methods

#include "pylith/feassemble/ParameterManager.hh" // USES ParameterManager

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::SlipTimeFn::SlipTimeFn(void) :
  _parameters(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::SlipTimeFn::~SlipTimeFn(void)
{ // destructor
  delete _parameters; _parameters = 0;
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::faults::SlipTimeFn::SlipTimeFn(const SlipTimeFn& f) :
  _parameters(0)
{ // copy constructor
} // copy constructor


// End of file 
