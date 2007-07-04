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
 * Cells are 0-1, vertices are 2-5.
 *
 *              3
 *             /|\
 *            / | \
 *           /  |  \
 *          /   |   \
 *         2    |    5
 *          \   |   /
 *           \  |  /
 *            \ | /
 *             \|/
 *              4
 *
 *  The boundary is
 *
 *              3
 *             / \
 *            /   \
 *           /     \
 *          /       \
 *         2         5
 *          \       /
 *           \     /
 *            \   /
 *             \ /
 *              4
 */

#include "BoundaryDataTri3.hh"

const int pylith::faults::BoundaryDataTri3::_numVertices = 4;

const int pylith::faults::BoundaryDataTri3::_spaceDim = 2;

const int pylith::faults::BoundaryDataTri3::_numCells = 4;

const int pylith::faults::BoundaryDataTri3::_cellDim = 1;

const double pylith::faults::BoundaryDataTri3::_vertices[] = {
 -1.0,  0.0,
  0.0,  1.0,
  0.0, -1.0,
  1.0,  0.0,
};

const int pylith::faults::BoundaryDataTri3::_numCorners[] = {
  2,
  2,
  2,
  2
};

const int pylith::faults::BoundaryDataTri3::_cells[] = {
  3,  2,
  2,  4,
  5,  3,
  4,  5,
};

const char* pylith::faults::BoundaryDataTri3::_filename = "data/tri3traction.mesh";

pylith::faults::BoundaryDataTri3::BoundaryDataTri3(void)
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

pylith::faults::BoundaryDataTri3::~BoundaryDataTri3(void)
{}


// End of file
