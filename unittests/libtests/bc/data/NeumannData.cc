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
// Copyright (c) 2010 University of California, Davis
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
  numBoundaryVertices(0),
  numBoundaryCells(0),
  numCorners(0),
  cellVertices(0),
  tractionsCell(0),
  valsResidual(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::NeumannData::~NeumannData(void)
{ // destructor
} // destructor


// End of file
