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
 * Cells are 0-1, vertices are 2-6.
 *
 * 2   3,4,5  6
 *
 *     ^^^^^ Face in x-y plane
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,10, vertices are 2-9.
 *
 * 2   3,4,5  7,8,9   6
 *
 *     ^^^^^^^^^^^^ Cohesive element in x-y plane.
 */

#include "BoundaryDataTet4.hh"

const int pylith::faults::BoundaryDataTet4::_numVertices = 5;

const int pylith::faults::BoundaryDataTet4::_spaceDim = 3;

const int pylith::faults::BoundaryDataTet4::_numCells = 6;

const int pylith::faults::BoundaryDataTet4::_cellDim = 2;

const double pylith::faults::BoundaryDataTet4::_vertices[] = {
  -1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0,
   1.0,  0.0,  0.0,
};

const int pylith::faults::BoundaryDataTet4::_numCorners[] = {
  3,
  3,
  3,
  3,
  3,
  3
};

const int pylith::faults::BoundaryDataTet4::_cells[] = {
  3,  4,  2,
  2,  5,  3,
  4,  5,  2,
  6,  4,  3,
  3,  5,  6,
  5,  4,  6,
};

const char* pylith::faults::BoundaryDataTet4::_filename = "data/tet4traction.mesh";

pylith::faults::BoundaryDataTet4::BoundaryDataTet4(void)
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

pylith::faults::BoundaryDataTet4::~BoundaryDataTet4(void)
{}


// End of file
