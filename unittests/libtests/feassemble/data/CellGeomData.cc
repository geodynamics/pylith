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

#include "CellGeomData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::CellGeomData::CellGeomData(void) :
  cellDim(0),
  spaceDim(0),
  numCorners(0),
  numLocs(0),
  gravityVec(0),
  vertices(0),
  locations(0),
  jacobian(0),
  jacobianDet(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::CellGeomData::~CellGeomData(void)
{ // destructor
} // destructor


// End of file
