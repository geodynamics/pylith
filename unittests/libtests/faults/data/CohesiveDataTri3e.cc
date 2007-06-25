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
 * Cells are 0-3, vertices are 4-9.
 *
 *
 *         9
 *        / \
 *       /   \
 *      /  2  \
 *     /       \
 *    8---------5
 *     \       /|\
 *      \  3  / | \
 *       \   /  |  \
 *        \ /   |   \
 *         4  0 | 1  7
 *          \   |   /
 *           \  |  /
 *            \ | /
 *             \|/
 *              6
 *
 *
 * After adding cohesive elements
 *
 * Cells are 0-3, 13-14, vertices are 4-12.
 *
 *         9
 *        / \
 *       /   \
 *      /     \
 *     /       \
 *    8---------  5
 *    |          /|
 *   12--------10 |\
 *     \       /| | \
 *      \     / | |  \
 *       \   /  | |   \
 *        \ /   | |    \
 *         4    | |    7
 *          \   | |   /
 *           \  | |  /
 *            \ | | /
 *             \| |/
 *             11-6
 */

#include "CohesiveDataTri3e.hh"

const int pylith::faults::CohesiveDataTri3e::_numVertices = 9;

const int pylith::faults::CohesiveDataTri3e::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3e::_numCells = 6;

const int pylith::faults::CohesiveDataTri3e::_cellDim = 2;

const double pylith::faults::CohesiveDataTri3e::_vertices[] = {
 -1.0,  0.0,
  0.0,  1.0,
  0.0, -1.0,
  1.0,  0.0,
 -2.0,  1.0,
 -1.0,  2.0,
  0.0,  1.0,
  0.0, -1.0,
 -2.0,  1.0,
};

const int pylith::faults::CohesiveDataTri3e::_numCorners[] = {
  3,
  3,
  3,
  3,
  4,
  4,
};

const int pylith::faults::CohesiveDataTri3e::_cells[] = {
  4, 11, 10,
  5,  6,  7,
  8,  5,  9,
 12,  4, 10,
  5,  6, 10, 11,
  8,  5, 12, 10,
};

const int pylith::faults::CohesiveDataTri3e::_materialIds[] = {
  0,  0,  0,  0,
  1,  1
};

const int pylith::faults::CohesiveDataTri3e::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3e::_groupSizes[] = 
  { 6, 5 };

const int pylith::faults::CohesiveDataTri3e::_groups[] = {
  5, 6, 8, 10, 11, 12,
  5, 6, 7, 10, 11
};

const char* pylith::faults::CohesiveDataTri3e::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataTri3e::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3e::_filename = "data/tri3e.mesh";

pylith::faults::CohesiveDataTri3e::CohesiveDataTri3e(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  vertices = const_cast<double*>(_vertices);
  numCorners = const_cast<int*>(_numCorners);
  cells = const_cast<int*>(_cells);
  materialIds = const_cast<int*>(_materialIds);
  groups = const_cast<int*>(_groups);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  filename = const_cast<char*>(_filename);
} // constructor

pylith::faults::CohesiveDataTri3e::~CohesiveDataTri3e(void)
{}


// End of file
