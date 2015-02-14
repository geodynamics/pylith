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

#include "CohesiveImpulsesData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::CohesiveImpulsesData::CohesiveImpulsesData(void) :
  meshFilename(0),
  lengthScale(1.0e+3),
  pressureScale(2.25e+10),
  densityScale(1.0),
  timeScale(2.0),
  spaceDim(0),
  cellDim(0),
  numBasis(0),
  numQuadPts(0),
  quadPts(0),
  quadWts(0),
  basis(0),
  basisDeriv(0),
  verticesRef(0),
  id(0),
  label(0),
  impulseAmpFilename(0),
  impulseDOF(0),
  numComponents(0),
  fieldT(0),
  fieldIncr(0),
  orientation(0),
  area(0),
  amplitude(0),
  numImpulses(0),
  residual(0),
  constraintEdges(0),
  numConstraintEdges(0)
{ // constructor
  const PylithScalar velScale = lengthScale / timeScale;
  densityScale = pressureScale / (velScale*velScale);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::faults::CohesiveImpulsesData::~CohesiveImpulsesData(void)
{ // destructor
} // destructor


// End of file
