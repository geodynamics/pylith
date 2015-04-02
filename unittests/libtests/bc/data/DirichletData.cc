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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
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
  dbFilename(0),
  lengthScale(1.0e+3),
  pressureScale(2.25e+10),
  densityScale(1.0),
  timeScale(2.0)
{ // constructor
  const PylithScalar velScale = lengthScale / timeScale;
  densityScale = pressureScale / (velScale*velScale);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::DirichletData::~DirichletData(void)
{ // destructor
} // destructor


// End of file
