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

#include "DirichletPointsData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::DirichletPointsData::DirichletPointsData(void) :
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
pylith::bc::DirichletPointsData::~DirichletPointsData(void)
{ // destructor
} // destructor


// End of file
