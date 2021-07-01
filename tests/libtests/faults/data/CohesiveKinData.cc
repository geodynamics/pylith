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

#include "CohesiveKinData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::CohesiveKinData::CohesiveKinData(void) :
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
  edge(0),
  finalSlipFilename(0),
  slipTimeFilename(0),
  riseTimeFilename(0),
  fieldT(0),
  fieldIncr(0),
  jacobianLumped(0),
  orientation(0),
  area(0),
  residual(0),
  jacobian(0),
  fieldIncrAdjusted(0),
  verticesFault(0),
  edgesLagrange(0),
  verticesPositive(0),
  verticesNegative(0),
  numFaultVertices(0),
  numCohesiveCells(0),
  cellMappingFault(0),
  cellMappingCohesive(0)
{ // constructor
  const PylithScalar velScale = lengthScale / timeScale;
  densityScale = pressureScale / (velScale*velScale);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::faults::CohesiveKinData::~CohesiveKinData(void)
{ // destructor
} // destructor


// End of file
