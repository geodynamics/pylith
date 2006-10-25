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

#include "Integrator.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Integrator::Integrator(void) :
  _quadrature(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Integrator::~Integrator(void)
{ // destructor
  delete _quadrature; _quadrature = 0;
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor
pylith::feassemble::Integrator::Integrator(const Integrator& i) :
  _quadrature(0)
{ // copy constructor
  if (0 != i._quadrature)
    _quadrature = i._quadrature->clone();
} // copy constructor

// ----------------------------------------------------------------------
// Set quadrature for integration.
void
pylith::feassemble::Integrator::quadrature(const Quadrature* q)
{ // quadrature
  delete _quadrature;
  _quadrature = (0 != q) ? q->clone() : 0;
} // quadrature

// End of file 
