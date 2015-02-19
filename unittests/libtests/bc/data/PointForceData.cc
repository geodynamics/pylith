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

#include "PointForceData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::PointForceData::PointForceData(void) :
  tRef(0),
  forceRate(0),
  tResidual(0),
  numDOF(0),
  numForceDOF(0),
  numForcePts(0),
  id(0),
  label(0),
  forceDOF(0),
  forcePoints(0),
  forceInitial(0),
  residual(0),
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
pylith::bc::PointForceData::~PointForceData(void)
{ // destructor
} // destructor


// End of file
