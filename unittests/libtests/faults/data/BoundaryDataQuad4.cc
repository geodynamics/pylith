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
 * Cells are 0-1, vertices are 2-7.
 *
 *       3 -------- 5 -------- 7
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       2 -------- 4 -------- 6
 *
 * Taking the boundary
 *
 *       3 -------- 5 -------- 7
 *       |                     |
 *       |                     |
 *       |                     |
 *       |                     |
 *       |                     |
 *       |                     |
 *       |                     |
 *       |                     |
 *       2 -------- 4 -------- 6
 */

#include "BoundaryDataQuad4.hh"

const int pylith::faults::BoundaryDataQuad4::_numVertices = 6;

const int pylith::faults::BoundaryDataQuad4::_spaceDim = 2;

const int pylith::faults::BoundaryDataQuad4::_numCells = 6;

const int pylith::faults::BoundaryDataQuad4::_cellDim = 1;

const double pylith::faults::BoundaryDataQuad4::_vertices[] = {
  -2.0, -1.0,
  -2.0,  1.0,
   0.0, -1.0,
   0.0,  1.0,
   2.0, -1.0,
   2.0,  1.0,
};

const int pylith::faults::BoundaryDataQuad4::_numCorners[] = {
  2,
  2,
  2,
  2,
  2,
  2
};

const int pylith::faults::BoundaryDataQuad4::_cells[] = {
  3,  2,
  2,  4,
  5,  3,
  4,  6,
  7,  5,
  6,  7,
};

const char* pylith::faults::BoundaryDataQuad4::_filename = 
  "data/quad4traction.mesh";

pylith::faults::BoundaryDataQuad4::BoundaryDataQuad4(void)
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

pylith::faults::BoundaryDataQuad4::~BoundaryDataQuad4(void)
{}


// End of file
