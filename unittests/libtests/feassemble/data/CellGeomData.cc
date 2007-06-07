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

#include "CellGeomData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::CellGeomData::CellGeomData(void) :
  cellDim(0),
  spaceDim(0),
  numCorners(0),
  numLocs(0),
  vertices(0),
  locations(0),
  jacobian(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::CellGeomData::~CellGeomData(void)
{ // destructor
} // destructor


// End of file
