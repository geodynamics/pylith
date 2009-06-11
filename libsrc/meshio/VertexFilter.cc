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
