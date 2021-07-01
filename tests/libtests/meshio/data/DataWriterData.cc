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

#include "DataWriterData.hh"

const int pylith::meshio::DataWriterData::DataWriterData::numVertexFields = 4;
const int pylith::meshio::DataWriterData::DataWriterData::numCellFields = 4;

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
  numVertices(0),
  vertexFieldsInfo(0),
  numCells(0),
  cellFieldsInfo(0)
{ // constructor
  for (int i=0; i < numVertexFields; ++i) {
    vertexFields[i] = 0;
  } // for

  for (int i=0; i < numCellFields; ++i) {
    cellFields[i] = 0;
  } // for
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriterData::~DataWriterData(void)
{ // destructor
} // destructor


// End of file
