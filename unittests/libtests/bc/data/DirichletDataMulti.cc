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

#include "DirichletDataMulti.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::DirichletDataMulti::DirichletDataMulti(void) :
  numDOF(0),
  numFixedDOFA(0),
  numConstrainedPtsA(0),
  idA(0),
  labelA(0),
  fixedDOFA(0),
  constrainedPointsA(0),
  valuesA(0),
  dbFilenameA(0),
  numFixedDOFB(0),
  numConstrainedPtsB(0),
  idB(0),
  labelB(0),
  fixedDOFB(0),
  constrainedPointsB(0),
  valuesB(0),
  dbFilenameB(0),
  field(0),
  constraintSizes(0),
  constrainedDOF(0),
  meshFilename(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::DirichletDataMulti::~DirichletDataMulti(void)
{ // destructor
} // destructor


// End of file
