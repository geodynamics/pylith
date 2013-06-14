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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type, typename field_type>
pylith::meshio::CellFilter<mesh_type, field_type>::CellFilter(void) :
  _quadrature(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type, typename field_type>
pylith::meshio::CellFilter<mesh_type, field_type>::~CellFilter(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type, typename field_type>
pylith::meshio::CellFilter<mesh_type, field_type>::CellFilter(const CellFilter& f) :
  _quadrature(0)
{ // copy constructor
  if (f._quadrature)
    _quadrature = new feassemble::Quadrature<mesh_type>(*f._quadrature);
} // copy constructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type, typename field_type>
void
pylith::meshio::CellFilter<mesh_type, field_type>::deallocate(void)
{ // deallocate
  delete _quadrature; _quadrature = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Set quadrature associated with cells.
template<typename mesh_type, typename field_type>
void
pylith::meshio::CellFilter<mesh_type, field_type>::quadrature(const feassemble::Quadrature<mesh_type>* q)
{ // quadrature
  delete _quadrature; 
  _quadrature = (q) ? new feassemble::Quadrature<mesh_type>(*q) : 0;
} // quadrature


// End of file
