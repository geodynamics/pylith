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

#include "AbsorbingDampersData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::AbsorbingDampersData::AbsorbingDampersData(void) :
  meshFilename(0),
  spaceDim(0),
  cellDim(0),
  numVertices(0),
  numCells(0),
  vertices(0),
  cells(0),
  verticesRef(0),
  numBasis(0),
  numQuadPts(0),
  quadPts(0),
  quadWts(0),
  basis(0),
  basisDerivRef(0),
  spatialDBFilename(0),
  dt(0),
  fieldTpdt(0),
  fieldT(0),
  fieldTmdt(0),
  valsResidual(0),
  valsJacobian(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::bc::AbsorbingDampersData::~AbsorbingDampersData(void)
{ // destructor
} // destructor


// End of file
