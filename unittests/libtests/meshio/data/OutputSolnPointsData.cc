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

#include "OutputSolnPointsData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnPointsData::OutputSolnPointsData(void) :
  meshFilename(0),
  spaceDim(0),
  numPoints(0),
  points(0),
  numVertices(0),
  fiberDim(0),
  field(0),
  fieldInterp(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnPointsData::~OutputSolnPointsData(void)
{ // destructor
} // destructor


// End of file
