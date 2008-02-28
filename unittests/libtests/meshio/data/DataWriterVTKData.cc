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
pylith::meshio::DataWriterVTKData::~DataWriterVTKData(void)
{ // destructor
} // destructor


// End of file
