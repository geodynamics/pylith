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

#include "VertexFilter.hh" // Implementation of class methods

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::VertexFilter::VertexFilter(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::VertexFilter::~VertexFilter(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::VertexFilter::deallocate(void)
{ // deallocate
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::VertexFilter::VertexFilter(const VertexFilter& f)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// operator=.
const pylith::meshio::VertexFilter&
pylith::meshio::VertexFilter::operator=(const VertexFilter& f)
{ // operator=
  return f;
} // operator=


// End of file
