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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include "pylith/utils/petscerror.h" // USES PYLITH_METHOD_BEGIN/END

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
  PYLITH_METHOD_BEGIN;

  if (f._quadrature)
    _quadrature = new feassemble::Quadrature<mesh_type>(*f._quadrature);

  PYLITH_METHOD_END;
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
  PYLITH_METHOD_BEGIN;

  delete _quadrature; 
  _quadrature = (q) ? new feassemble::Quadrature<mesh_type>(*q) : 0;

  PYLITH_METHOD_END;
} // quadrature


// End of file
