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

#include "DataWriterDataPoints.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriterDataPoints::DataWriterDataPoints(void) :
  numPoints(0),
  spaceDim(0),
  points(0),
  names(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriterDataPoints::~DataWriterDataPoints(void)
{ // destructor
} // destructor


// End of file
