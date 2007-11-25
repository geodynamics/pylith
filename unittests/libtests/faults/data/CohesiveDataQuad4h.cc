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
 * Cells are 0-8, vertices are 9-24.
 *
 * 9 ----10 ----11 ----12
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 *13 ----14 ----15 ----16
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 *17 ----18 ----19 ----20
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 * |      |      |      |
 *21 ----22 ----23 ----24
 *
 * After adding cohesive elements
 *
 * Cells are 0-8,16-17 vertices are 9-24.
 *
 * 9 ----10 ----25-11 ----12
 * |      |      |  |      |
 * |  0   |  1   |28|  2   |
 * |      |      |  |      |
 * |      |      |  |      |
 *13 ----14 ----26-15 ----16
 * |      |      |  |      |
 * |  3   |  4   |29|  5   |
 * |      |      |  |      |
 * |      |      |  |      |
 *30 ----31 ----27-19 ----20
 * |  32  |/       |       |
 *17 ----18        |       |
 * |      |   7    |   8   |
 * |  6   |        |       |
 * |      |        |       |
 * |      |        |       |
 *21 ----22 ------23------24
 *
 */

#include "CohesiveDataQuad4h.hh"

const int pylith::faults::CohesiveDataQuad4h::_numVertices = 21;

const int pylith::faults::CohesiveDataQuad4h::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4h::_numCells = 12;

const int pylith::faults::CohesiveDataQuad4h::_cellDim = 2;

const double pylith::faults::CohesiveDataQuad4h::_vertices[] = {
  -3.0,  3.0,
  -1.0,  3.0,
   1.0,  3.0,
   3.0,  3.0,
  -3.0,  1.0,
  -1.0,  1.0,
   1.0,  1.0,
   3.0,  1.0,
  -3.0, -1.0,
  -1.0, -1.0,
   1.0, -1.0,
   3.0, -1.0,
  -3.0, -3.0,
  -1.0, -3.0,
   1.0, -3.0,
   3.0, -3.0,
   1.0,  3.0,
   1.0,  1.0,
   1.0, -1.0,
  -3.0, -1.0,
  -1.0, -1.0,
};

const int pylith::faults::CohesiveDataQuad4h::_numCorners[] = {
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
  4,
};

const int pylith::faults::CohesiveDataQuad4h::_cells[] = {
   9, 13, 14, 10,
  10, 14, 26, 25,
  11, 15, 16, 12,
  13, 30, 31, 14,
  14, 18, 27, 26,
  15, 19, 20, 16,
  17, 21, 22, 18,
  18, 22, 23, 19,
  19, 23, 24, 20,
  11, 15, 25, 26,
  15, 19, 26, 27,
  18, 17, 31, 30,
};

const int pylith::faults::CohesiveDataQuad4h::_materialIds[] = {
  10, 10, 11, 10, 10, 11, 12, 12, 11,
  1, 1, 1, 2, 2,
};

const int pylith::faults::CohesiveDataQuad4h::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4h::_groupSizes[] = 
  { 6, 4 };

const int pylith::faults::CohesiveDataQuad4h::_groups[] = {
  11, 15, 19, 25, 26, 27,
  17, 18, 30, 31,
};

const char* pylith::faults::CohesiveDataQuad4h::_groupNames[] = {
  "faultA", "faultB"
};

const char* pylith::faults::CohesiveDataQuad4h::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataQuad4h::_filename = 
  "data/quad4h.mesh";

pylith::faults::CohesiveDataQuad4h::CohesiveDataQuad4h(void)
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

pylith::faults::CohesiveDataQuad4h::~CohesiveDataQuad4h(void)
{}


// End of file
