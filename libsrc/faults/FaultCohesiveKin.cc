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

#include "FaultCohesiveKin.hh" // implementation of object methods

#include "EqKinSrc.hh" // USES EqKinSrc

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(void) :
  _eqsrc(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveKin::~FaultCohesiveKin(void)
{ // destructor
  delete _eqsrc; _eqsrc = 0;
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(const FaultCohesiveKin& f) :
  FaultCohesive(f),
  _eqsrc(0)
{ // copy constructor
  if (0 != f._eqsrc)
    _eqsrc = f._eqsrc->clone();
} // copy constructor

// ----------------------------------------------------------------------
// Set kinematic earthquake source.
void
pylith::faults::FaultCohesiveKin::eqsrc(EqKinSrc* src)
{ // eqsrc
  delete _eqsrc; _eqsrc = (0 != src) ? src->clone() : 0;
} // eqsrc


// End of file 
