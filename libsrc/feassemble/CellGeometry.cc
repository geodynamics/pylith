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

#include <portinfo>

#include "CellGeometry.hh" // implementation of class methods

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::CellGeometry::CellGeometry(const int cellDim,
					     const int spaceDim,
					     const int numCorners) :
  _cellDim(cellDim),
  _spaceDim(spaceDim),
  _numCorners(numCorners)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::CellGeometry::~CellGeometry(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::CellGeometry::CellGeometry(const CellGeometry& g) :
  _cellDim(g._cellDim),
  _spaceDim(g._spaceDim),
  _numCorners(g._numCorners)
{ // copy constructor
} // copy constructor


// End of file
