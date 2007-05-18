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

#include "IntegratorImplicit.hh" // implementation of class methods

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorImplicit::IntegratorImplicit(void) :
  Integrator(),
  _dt(-1.0),
  _dtm1(-1.0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorImplicit::~IntegratorImplicit(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor
pylith::feassemble::IntegratorImplicit::IntegratorImplicit(const IntegratorImplicit& i) :
  Integrator(i),
  _dt(i._dt),
  _dtm1(i._dtm1)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::IntegratorImplicit::timeStep(const double dt)
{ // timeStep
  if (_dt != -1.0)
    _dtm1 = _dt;
  else
    _dtm1 = dt;
  _dt = dt;
  assert(_dt == _dtm1); // For now, don't allow variable time step
} // timeStep

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
double
pylith::feassemble::IntegratorImplicit::stableTimeStep(void) const
{ // stableTimeStep
  // Default is current time step
  return _dt;
} // stableTimeStep


// End of file 
