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

#include "CohesiveDataHex8h.hh"

const int pylith::faults::CohesiveDataHex8h::_numVertices = 24;

const int pylith::faults::CohesiveDataHex8h::_spaceDim = 3;

const int pylith::faults::CohesiveDataHex8h::_numCells = 6;

const int pylith::faults::CohesiveDataHex8h::_cellDim = 3;

const double pylith::faults::CohesiveDataHex8h::_vertices[] = {
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

const int pylith::faults::CohesiveDataHex8h::_numCorners[] = {
  8,
  8,
  8,
  8,
  8,
  8
};

const int pylith::faults::CohesiveDataHex8h::_cells[] = {
   4, 10, 11,  5,  6, 12, 13,  7,
  19, 18, 20, 21, 25, 24, 26, 27,
  23, 17, 19, 25, 22, 16, 18, 24,
  15, 14,  8,  9, 13, 12,  6,  7,
  11, 13, 12, 10, 23, 25, 24, 22,
  15, 14, 12, 13, 27, 26, 24, 25, 
};

const int pylith::faults::CohesiveDataHex8h::_materialIds[] = {
  0,  0,  0,  0,
  1,  1
};

const int pylith::faults::CohesiveDataHex8h::_numGroups = 2;

const int pylith::faults::CohesiveDataHex8h::_groupSizes[] = 
  { 12, 8 };

const int pylith::faults::CohesiveDataHex8h::_groups[] = {
  10, 11, 12, 13, 14, 15, 22, 23, 24, 25, 26, 27,
  8, 9, 14, 15, 20, 21, 26, 27
};

const char* pylith::faults::CohesiveDataHex8h::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataHex8h::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataHex8h::_filename = 
  "data/hex8h.mesh";

pylith::faults::CohesiveDataHex8h::CohesiveDataHex8h(void)
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

pylith::faults::CohesiveDataHex8h::~CohesiveDataHex8h(void)
{}


// End of file
