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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

// ----------------------------------------------------------------------
// Constructor
template<typename field_type>
pylith::meshio::VertexFilter<field_type>::VertexFilter(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename field_type>
pylith::meshio::VertexFilter<field_type>::~VertexFilter(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename field_type>
void
pylith::meshio::VertexFilter<field_type>::deallocate(void)
{ // deallocate
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
template<typename field_type>
pylith::meshio::VertexFilter<field_type>::VertexFilter(const VertexFilter& f)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// operator=.
template<typename field_type>
const pylith::meshio::VertexFilter<field_type>&
pylith::meshio::VertexFilter<field_type>::operator=(const VertexFilter& f)
{ // operator=
} // operator=


// End of file
