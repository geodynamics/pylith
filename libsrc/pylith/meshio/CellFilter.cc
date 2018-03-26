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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "CellFilter.hh" // Implementation of class methods

#include "pylith/topology/Stratum.hh" // USES StratumIS

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::CellFilter::CellFilter(void) :
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
  _cellsIS(0)
{ // copy constructor
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_END;
} // copy constructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::CellFilter::deallocate(void)
{ // deallocate
  delete _cellsIS; _cellsIS = 0;
} // deallocate
  

// End of file
