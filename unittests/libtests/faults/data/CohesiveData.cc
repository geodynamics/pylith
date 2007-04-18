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

#include "CohesiveData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::CohesiveData::CohesiveData(void) :
  numVertices(0),
  spaceDim(0),
  numCells(0),
  cellDim(0),
  vertices(0),
  numCorners(0),
  cells(0),
  materialIds(0),
  groups(0),
  groupSizes(0),
  groupNames(0),
  groupTypes(0),
  numGroups(0),
  filename(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::faults::CohesiveData::~CohesiveData(void)
{ // destructor
} // destructor


// End of file
