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

#include "OutputFilter.hh" // implementation of class methods

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputFilter::OutputFilter(void) :
  _filterType(VERTEX_FILTER)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputFilter::~OutputFilter(void)
{ // destructor
} // destructor  

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::OutputFilter::OutputFilter(const OutputFilter& f) :
  _filterType(f._filterType)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// operator=.
const pylith::meshio::OutputFilter&
pylith::meshio::OutputFilter::operator=(const OutputFilter& f)
{ // operator=
  if (this != &f) {
    _filterType = f._filterType;
  } // if
} // operator=


// End of file
