// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "IntegratorExplicit.hh" // implementation of class methods

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorExplicit::IntegratorExplicit(void) :
  Integrator(),
  _dt(0.0),
  _dtm1(0.0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorExplicit::~IntegratorExplicit(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor
pylith::feassemble::IntegratorExplicit::IntegratorExplicit(const IntegratorExplicit& i) :
  Integrator(i),
  _dt(i._dt),
  _dtm1(i._dtm1)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::IntegratorExplicit::timeStep(const double dt)
{ // timeStep
  _dtm1 = _dt;
  _dt = dt;
  assert(_dt == _dtm1); // For now, don't allow variable time step
} // timeStep

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
double
pylith::feassemble::IntegratorExplicit::stableTimeStep(void) const
{ // stableTimeStep
  // Default is current time step
  return _dt;
} // stableTimeStep


// End of file 
