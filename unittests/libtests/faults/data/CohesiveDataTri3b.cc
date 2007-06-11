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
 *
 * After adding cohesive elements
 *
 * Cells are 0-1, 8, vertices are 2-7.
 *
 *              3 -- 6
 *             /|    |\
 *            / |    | \
 *           /  |    |  \
 *          /   |    |   \
 *         2    |    |    5
 *          \   |    |   /
 *           \  |    |  /
 *            \ |    | /
 *             \|    |/
 *              4 -- 7
 */

#include "CohesiveDataTri3b.hh"

const int pylith::faults::CohesiveDataTri3b::_numVertices = 6;

const int pylith::faults::CohesiveDataTri3b::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3b::_numCells = 3;

const int pylith::faults::CohesiveDataTri3b::_cellDim = 2;

const double pylith::faults::CohesiveDataTri3b::_vertices[] = {
 -1.0,  0.0,
  0.0,  1.0,
  0.0, -1.0,
  1.0,  0.0,
  0.0,  1.0,
  0.0, -1.0
};

const int pylith::faults::CohesiveDataTri3b::_numCorners[] = {
  3,
  3,
  4
};

const int pylith::faults::CohesiveDataTri3b::_cells[] = {
  4,  3,  2,
  7,  5,  6,
  3,  4,  6, 7
};

const int pylith::faults::CohesiveDataTri3b::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTri3b::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3b::_groupSizes[] = 
  { 4, 5 };

const int pylith::faults::CohesiveDataTri3b::_groups[] = {
  3, 4, 6, 7,
  3, 4, 5, 6, 7
};

const char* pylith::faults::CohesiveDataTri3b::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataTri3b::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3b::_filename = "data/tri3b.mesh";

pylith::faults::CohesiveDataTri3b::CohesiveDataTri3b(void)
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

pylith::faults::CohesiveDataTri3b::~CohesiveDataTri3b(void)
{}


// End of file
