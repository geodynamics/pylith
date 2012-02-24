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

#include "DirichletDataMulti.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::DirichletDataMulti::DirichletDataMulti(void) :
  numDOF(0), // General
  numFixedDOFA(0), // BC A
  numConstrainedPtsA(0),
  idA(0),
  labelA(0),
  fixedDOFA(0),
  constrainedPointsA(0),
  dbFilenameA(0),
  dbFilenameARate(0),
  tRefA(0),
  numFixedDOFB(0), // BC B
  numConstrainedPtsB(0),
  idB(0),
  labelB(0),
  fixedDOFB(0),
  constrainedPointsB(0),
  dbFilenameB(0),
  dbFilenameBRate(0),
  tRefB(0),
  numFixedDOFC(0), // BC C
  numConstrainedPtsC(0),
  idC(0),
  labelC(0),
  fixedDOFC(0),
  constrainedPointsC(0),
  dbFilenameC(0),
  dbFilenameCRate(0),
  tRefC(0),
  field(0), // General
  fieldIncr(0), // General
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
