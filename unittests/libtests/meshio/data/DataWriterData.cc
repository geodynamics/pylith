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

#include "DataWriterData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriterData::DataWriterData(void) :
  meshFilename(0),
  faultLabel(0),
  faultId(0),
  bcLabel(0),
  timestepFilename(0),
  vertexFilename(0),
  cellFilename(0),
  time(0),
  timeFormat(0),
  cellsLabel(0),
  labelId(0),
  numVertexFields(0),
  numVertices(0),
  vertexFieldsInfo(0),
  numCellFields(0),
  numCells(0),
  cellFieldsInfo(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriterData::~DataWriterData(void)
{ // destructor
} // destructor


// End of file
