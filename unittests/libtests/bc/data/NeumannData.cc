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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "NeumannData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::NeumannData::NeumannData(void) :
  meshFilename(0),
  lengthScale(1.0e+3),
  pressureScale(2.25e+10),
  densityScale(1.0),
  timeScale(2.0),
  numBasis(0),
  numQuadPts(0),
  quadPts(0),
  quadWts(0),
  basis(0),
  basisDerivRef(0),
  spatialDBFilename(0),
  id(0),
  label(""),
  spaceDim(0),
  cellDim(0),
  numVertices(0),
  numCells(0),
  numCorners(0),
  cells(0),
  tractionsCell(0),
  valsResidual(0)
{ // constructor
  const PylithScalar velScale = lengthScale / timeScale;
  densityScale = pressureScale / (velScale*velScale);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::NeumannData::~NeumannData(void)
{ // destructor
} // destructor


// End of file
