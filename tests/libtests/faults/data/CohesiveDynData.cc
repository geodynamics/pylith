// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include "CohesiveDynData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::CohesiveDynData::CohesiveDynData(void) :
  meshFilename(0),
  setLengthScale(1.0e+3),
  setPressureScale(2.25e+10),
  setDensityScale(1.0),
  setTimeScale(2.0),
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
  initialTractFilename(0),
  fieldT(0),
  fieldIncrStick(0),
  fieldIncrSlip(0),
  fieldIncrOpen(0),
  jacobian(0),
  orientation(0),
  initialTractions(0),
  area(0),
  slipStickE(0),
  fieldIncrSlipE(0),
  slipSlipE(0),
  fieldIncrOpenE(0),
  slipOpenE(0),
  constraintEdges(0),
  numConstraintEdges(0)
{ // constructor
  const PylithScalar velScale = lengthScale / timeScale;
  densityScale = pressureScale / (velScale*velScale);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::faults::CohesiveDynData::~CohesiveDynData(void)
{ // destructor
} // destructor


// End of file
