// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "BoundaryCondition.hh" // implementation of object methods

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::BoundaryCondition::BoundaryCondition(void) :
  _id(0),
  _label(""),
  _db(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::BoundaryCondition::~BoundaryCondition(void)
{ // destructor
  _db = 0;
} // destructor


// End of file 
