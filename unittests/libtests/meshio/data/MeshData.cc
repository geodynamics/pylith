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

#include "MeshData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshData::MeshData(void) :
  numVertices(0),
  spaceDim(0),
  numCells(0),
  cellDim(0),
  numCorners(0),
  vertices(0),
  cells(0),
  useIndexZero(true)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshData::~MeshData(void)
{ // destructor
} // destructor

// End of file
