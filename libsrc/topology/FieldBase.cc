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

#include "FieldBase.hh" // implementation of class methods

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::FieldBase::FieldBase(void) :
  _scale(1.0),
  _name("unknown"),
  _vecFieldType(OTHER),
  _dimensionsOkay(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::FieldBase::~FieldBase(void)
{ // destructor
} // destructor


// End of file 
