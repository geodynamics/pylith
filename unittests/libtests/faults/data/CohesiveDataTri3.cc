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

#include "CohesiveDataTri3.hh"

const int pylith::faults::CohesiveDataTri3::_numVertices = 6;

const int pylith::faults::CohesiveDataTri3::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3::_numCells = 3;

const int pylith::faults::CohesiveDataTri3::_cellDim = 2;

const double pylith::faults::CohesiveDataTri3::_vertices[] = {
 -1.0,  0.0,
  0.0,  1.0,
  0.0, -1.0,
  1.0,  0.0,
  0.0,  1.0,
  0.0, -1.0
};

const int pylith::faults::CohesiveDataTri3::_numCorners[] = {
  3,
  3,
  4
};

const int pylith::faults::CohesiveDataTri3::_cells[] = {
  2,  4,  3,
  6,  7,  5,
  3,  4,  6, 7
};

const int pylith::faults::CohesiveDataTri3::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTri3::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3::_groupSizes[] = 
  { 4, 5 };

const int pylith::faults::CohesiveDataTri3::_groups[] = {
  3, 4, 6, 7,
  3, 4, 5, 6, 7
};

const char* pylith::faults::CohesiveDataTri3::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataTri3::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3::_filename = "data/meshTri3A.txt";

pylith::faults::CohesiveDataTri3::CohesiveDataTri3(void)
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

pylith::faults::CohesiveDataTri3::~CohesiveDataTri3(void)
{}


// End of file
