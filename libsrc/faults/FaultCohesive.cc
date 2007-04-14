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
// Copy constructor.
pylith::faults::FaultCohesive::FaultCohesive(const FaultCohesive& f) :
  Fault(f)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesive::adjustTopology(ALE::Obj<ALE::Mesh>* mesh) const
{ // adjustTopology
} // adjustTopology


// End of file 
