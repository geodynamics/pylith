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

#include "BoundaryData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::BoundaryData::BoundaryData(void) :
  numVertices(0),
  spaceDim(0),
  numCells(0),
  cellDim(0),
  vertices(0),
  numCorners(0),
  cells(0),
  filename(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::faults::BoundaryData::~BoundaryData(void)
{ // destructor
} // destructor


// End of file
