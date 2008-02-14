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

#include "DirichletData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::DirichletData::DirichletData(void) :
  tRef(0),
  valueRate(0),
  numDOF(0),
  numFixedDOF(0),
  numConstrainedPts(0),
  id(0),
  label(0),
  fixedDOF(0),
  constrainedPoints(0),
  valuesInitial(0),
  meshFilename(0),
  dbFilename(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::DirichletData::~DirichletData(void)
{ // destructor
} // destructor


// End of file
