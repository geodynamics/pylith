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

#include "IntegratorData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorData::IntegratorData(void) :
  spaceDim(0),
  gravityVec(0),
  cellDim(0),
  numVertices(0),
  numCells(0),
  vertices(0),
  cells(0),
  verticesRef(0),
  numBasis(0),
  numQuadPts(0),
  quadPts(0),
  quadWts(0),
  basis(0),
  basisDerivRef(0),
  lengthScale(1.0e+3),
  pressureScale(2.25e+10),
  densityScale(1.0),
  timeScale(2.0),
  matType(0),
  matDBFilename(0),
  matId(0),
  matLabel(0),
  dt(0),
  fieldTIncr(0),
  fieldT(0),
  fieldTmdt(0),
  valsResidual(0),
  valsJacobian(0)
{ // constructor
  const PylithScalar velScale = lengthScale / timeScale;
  densityScale = pressureScale / (velScale*velScale);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorData::~IntegratorData(void)
{ // destructor
} // destructor


// End of file
