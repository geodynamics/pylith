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

#include "BoundaryMeshData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::BoundaryMeshData::BoundaryMeshData(void) :
  filename(0),
  bcLabel(0),
  faultLabel(0),
  faultId(0),
  verticesNoFault(0),
  cellsNoFault(0),
  verticesFault(0),
  cellsFault(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::BoundaryMeshData::~BoundaryMeshData(void)
{ // destructor
} // destructor


// End of file
