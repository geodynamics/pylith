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
 *              1
 *             /|\
 *            / | \
 *           /  |  \
 *          /   |   \
 *         0    |    3
 *          \   |   /
 *           \  |  /
 *            \ | /
 *             \|/
 *              2
 *
 *
 * After adding cohesive elements
 *
 *              1 -- 4
 *             /|    |\
 *            / |    | \
 *           /  |    |  \
 *          /   |    |   \
 *         0    |    |    3
 *          \   |    |   /
 *           \  |    |  /
 *            \ |    | /
 *             \|    |/
 *              2 -- 5
 */

#include "CohesiveDataTet4.hh"

const int pylith::faults::CohesiveDataTet4::_numVertices = 8;

const int pylith::faults::CohesiveDataTet4::_spaceDim = 3;

const int pylith::faults::CohesiveDataTet4::_numCells = 3;

const int pylith::faults::CohesiveDataTet4::_cellDim = 3;

const double pylith::faults::CohesiveDataTet4::_vertices[] = {
  -1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0,
   1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0
};

const int pylith::faults::CohesiveDataTet4::_numCorners[] = {
  4,
  4,
  6
};

const int pylith::faults::CohesiveDataTet4::_cells[] = {
  1,  2,  3,  0,
  5,  7,  6,  4,
  1,  3,  2,  5,  7,  6
};

const int pylith::faults::CohesiveDataTet4::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTet4::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4::_groupSizes[] = 
  { 3, 3 };

const int pylith::faults::CohesiveDataTet4::_groups[] = {
  1, 2, 3,
  0, 2, 3
};

const char* pylith::faults::CohesiveDataTet4::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataTet4::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4::_filename = "data/meshTet4A.txt";

pylith::faults::CohesiveDataTet4::CohesiveDataTet4(void)
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

pylith::faults::CohesiveDataTet4::~CohesiveDataTet4(void)
{}


// End of file
