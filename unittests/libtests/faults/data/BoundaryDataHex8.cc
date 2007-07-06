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

/* Original mesh
 *
 * Cells are 0-1 and vertices are 2-13.
 *
 *       2,3,4,5 -------- 6,7,8,9 -------- 10,11,12,13
 *
 *                        ^^^^^^^ Vertices forming fault
 *
 * Taking the boundary
 *
 */

#include "BoundaryDataHex8.hh"

const int pylith::faults::BoundaryDataHex8::_numVertices = 12;

const int pylith::faults::BoundaryDataHex8::_spaceDim = 3;

const int pylith::faults::BoundaryDataHex8::_numCells = 10;

const int pylith::faults::BoundaryDataHex8::_cellDim = 2;

const double pylith::faults::BoundaryDataHex8::_vertices[] = {
  -2.0, -1.0, -1.0,
  -2.0,  1.0, -1.0,
  -2.0, -1.0,  1.0,
  -2.0,  1.0,  1.0,
   0.0, -1.0, -1.0,
   0.0,  1.0, -1.0,
   0.0, -1.0,  1.0,
   0.0,  1.0,  1.0,
   2.0, -1.0, -1.0,
   2.0,  1.0, -1.0,
   2.0, -1.0,  1.0,
   2.0,  1.0,  1.0,
};

const int pylith::faults::BoundaryDataHex8::_numCorners[] = {
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4
};

const int pylith::faults::BoundaryDataHex8::_cells[] = {
  4,  5,  3,  2,
  3,  7,  6,  2,
  6,  8,  4,  2,
  5,  9,  7,  3,
  8,  9,  5,  4,
  7, 11, 10,  6,
 10, 12,  8,  6,
  9, 13, 11,  7,
 12, 13,  9,  8,
 11, 13, 12, 10,
};

const char* pylith::faults::BoundaryDataHex8::_filename = 
  "data/hex8traction.mesh";

pylith::faults::BoundaryDataHex8::BoundaryDataHex8(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  vertices = const_cast<double*>(_vertices);
  numCorners = const_cast<int*>(_numCorners);
  cells = const_cast<int*>(_cells);
  filename = const_cast<char*>(_filename);
} // constructor

pylith::faults::BoundaryDataHex8::~BoundaryDataHex8(void)
{}


// End of file
