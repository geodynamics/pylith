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

#include "CellFilter.hh" // implementation of class methods

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::CellFilter::CellFilter(void) :
  _quadrature(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::CellFilter::~CellFilter(void)
{ // destructor
  delete _quadrature; _quadrature = 0;
} // destructor  

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::CellFilter::CellFilter(const CellFilter& f) :
  _quadrature(0)
{ // copy constructor
  if (0 != f._quadrature)
    _quadrature = f._quadrature->clone();
} // copy constructor

// ----------------------------------------------------------------------
// Set quadrature associated with cells.
void
pylith::meshio::CellFilter::quadrature(const feassemble::Quadrature* q)
{ // quadrature
    delete _quadrature; _quadrature = (0 != q) ? q->clone() : 0;
} // quadrature


// End of file
