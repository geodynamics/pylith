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

#include "EqKinSrc.hh" // implementation of object methods

#include "SlipTimeFn.hh" // USES SlipTimeFn

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::EqKinSrc::EqKinSrc(void) :
  _slipfn(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::EqKinSrc::~EqKinSrc(void)
{ // destructor
  delete _slipfn; _slipfn = 0;
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::faults::EqKinSrc::EqKinSrc(const EqKinSrc& s) :
  _slipfn(0)
{ // copy constructor
  if (0 != s._slipfn)
    _slipfn = s._slipfn->clone();
} // copy constructor


// End of file 
