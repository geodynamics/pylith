// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include "DataWriterVTKData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriterVTKData::DataWriterVTKData(void) :
  timestepFilename(0),
  vertexFilename(0),
  cellFilename(0),
  time(0),
  timeFormat(0),
  numVertexFields(0),
  numVertices(0),
  vertexFields(0),
  numCellFields(0),
  numCells(0),
  cellFields(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriterVTKData::~DataWriterVTKData(void)
{ // destructor
} // destructor


// End of file
