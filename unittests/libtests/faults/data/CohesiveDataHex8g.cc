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
 * After adding cohesive elements
 *
 * Cells are 0-1,16 and vertices are 4-15.
 *
 *       2,3,4,5 -------- 6,7,8,9 -- 14,15,16,17 -------- 10,11,12,13
 *
 *                        ^^^^^^^^^^^^^^^^^^^^^^ Cohesive element
 *
 */

#include "CohesiveDataHex8g.hh"

const int pylith::faults::CohesiveDataHex8g::_numVertices = 24;

const int pylith::faults::CohesiveDataHex8g::_spaceDim = 3;

const int pylith::faults::CohesiveDataHex8g::_numCells = 6;

const int pylith::faults::CohesiveDataHex8g::_cellDim = 3;

const double pylith::faults::CohesiveDataHex8g::_vertices[] = {
  -2.0, -1.0, -2.0,
  -2.0,  1.0, -2.0,
  -2.0, -1.0,  0.0,
  -2.0,  1.0,  0.0,
  -2.0, -1.0,  2.0,
  -2.0,  1.0,  2.0,
   0.0, -1.0, -2.0,
   0.0,  1.0, -2.0,
   0.0, -1.0,  0.0,
   0.0,  1.0,  0.0,
   0.0, -1.0,  2.0,
   0.0,  1.0,  2.0,
   2.0, -1.0, -2.0,
   2.0,  1.0, -2.0,
   2.0, -1.0,  0.0,
   2.0,  1.0,  0.0,
   2.0, -1.0,  2.0,
   2.0,  1.0,  2.0,
   0.0, -1.0, -2.0,
   0.0,  1.0, -2.0,
   0.0, -1.0,  0.0,
   0.0,  1.0,  0.0,
   0.0, -1.0,  2.0,
   0.0,  1.0,  2.0,
};

const int pylith::faults::CohesiveDataHex8g::_numCorners[] = {
  8,
  8,
  8,
  8,
  8,
  8
};

const int pylith::faults::CohesiveDataHex8g::_cells[] = {
   4, 10, 11,  5,  6, 12, 13,  7,
   6, 12, 13,  7,  8, 14, 15,  9,
  22, 16, 17, 23, 24, 18, 19, 25,
  24, 18, 19, 25, 26, 20, 21, 27,
  10, 11, 13, 12, 22, 23, 25, 24,
  12, 13, 15, 14, 24, 25, 27, 26,
};

const int pylith::faults::CohesiveDataHex8g::_materialIds[] = {
  0,  0,  0,  0,
  1,  1
};

const int pylith::faults::CohesiveDataHex8g::_numGroups = 2;

const int pylith::faults::CohesiveDataHex8g::_groupSizes[] = 
  { 12, 8 };

const int pylith::faults::CohesiveDataHex8g::_groups[] = {
  10, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 27,
  8, 9, 14, 15, 20, 21, 26, 27
};

const char* pylith::faults::CohesiveDataHex8g::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataHex8g::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataHex8g::_filename = 
  "data/hex8g.mesh";

pylith::faults::CohesiveDataHex8g::CohesiveDataHex8g(void)
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

pylith::faults::CohesiveDataHex8g::~CohesiveDataHex8g(void)
{}


// End of file
