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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "AbsorbingDampersData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::AbsorbingDampersData::AbsorbingDampersData(void) :
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
  dt(0),
  fieldTIncr(0),
  fieldT(0),
  fieldTmdt(0),
  spaceDim(0),
  cellDim(0),
  numVertices(0),
  numCells(0),
  numCorners(0),
  cells(0),
  dampingConsts(0),
  valsResidual(0),
  valsJacobian(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::AbsorbingDampersData::~AbsorbingDampersData(void)
{ // destructor
} // destructor


// End of file
