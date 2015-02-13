// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "CellFilter.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/Stratum.hh" // USES StratumIS

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::CellFilter::CellFilter(void) :
  _quadrature(0),
  _cellsIS(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::CellFilter::~CellFilter(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::CellFilter::CellFilter(const CellFilter& f) :
  _quadrature(0),
  _cellsIS(0)
{ // copy constructor
  PYLITH_METHOD_BEGIN;

  if (f._quadrature)
    _quadrature = new feassemble::Quadrature(*f._quadrature);

  PYLITH_METHOD_END;
} // copy constructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::CellFilter::deallocate(void)
{ // deallocate
  delete _quadrature; _quadrature = 0;
  delete _cellsIS; _cellsIS = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Set quadrature associated with cells.
void
pylith::meshio::CellFilter::quadrature(const feassemble::Quadrature* q)
{ // quadrature
  PYLITH_METHOD_BEGIN;

  delete _quadrature; 
  _quadrature = (q) ? new feassemble::Quadrature(*q) : 0;

  PYLITH_METHOD_END;
} // quadrature


// End of file
