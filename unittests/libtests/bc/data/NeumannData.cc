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

#include "NeumannData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::NeumannData::NeumannData(void) :
  meshFilename(0),
  spaceDim(0),
  cellDim(0),
  numBasis(0),
  numQuadPts(0),
  quadPts(0),
  quadWts(0),
  basis(0),
  basisDeriv(0),
  verticesRef(0),
  id(0),
  label(0),
  dbFilename(0),
  numBoundaryCells(0),
  numVertices(0),
  numCorners(0),
  cells(0),
  tractionsCell(0),
  valsResidual(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::NeumannData::~NeumannData(void)
{ // destructor
} // destructor


// End of file
