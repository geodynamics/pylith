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

#include "CohesiveDataHex8e.hh"

const int pylith::faults::CohesiveDataHex8e::_numVertices = 16;

const int pylith::faults::CohesiveDataHex8e::_spaceDim = 3;

const int pylith::faults::CohesiveDataHex8e::_numCells = 3;

const int pylith::faults::CohesiveDataHex8e::_cellDim = 3;

const double pylith::faults::CohesiveDataHex8e::_vertices[] = {
  -2.0, -1.0, -1.0,
  -2.0,  1.0, -1.0,
  -2.0, -1.0,  1.0,
  -2.0,  1.0,  1.0,
   0.0, -1.0, -1.0,
   0.0,  1.0, -1.0,
   0.0, -1.0,  1.0,
   0.0,  1.0,  1.0,
   2.0, -1.0, -1.0,
   2.0,  1.0, -1.0,
   2.0, -1.0,  1.0,
   2.0,  1.0,  1.0,
   0.0, -1.0, -1.0,
   0.0,  1.0, -1.0,
   0.0, -1.0,  1.0,
   0.0,  1.0,  1.0,
};

const int pylith::faults::CohesiveDataHex8e::_numCorners[] = {
  8,
  8,
  8
};

const int pylith::faults::CohesiveDataHex8e::_cells[] = {
  5,  9,  8,  4,  3,  7,  6,  2,
 17, 13, 12, 16, 15, 11, 10, 14,
  6,  7,  9,  8, 14, 15, 17, 16,
};

const int pylith::faults::CohesiveDataHex8e::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataHex8e::_numGroups = 2;

const int pylith::faults::CohesiveDataHex8e::_groupSizes[] = 
  { 8, 8 };

const int pylith::faults::CohesiveDataHex8e::_groups[] = {
  6, 7, 8, 9, 14, 15, 16, 17,
  4, 5, 8, 9, 12, 13, 16, 17
};

const char* pylith::faults::CohesiveDataHex8e::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataHex8e::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataHex8e::_filename = 
  "data/hex8e.mesh";

pylith::faults::CohesiveDataHex8e::CohesiveDataHex8e(void)
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

pylith::faults::CohesiveDataHex8e::~CohesiveDataHex8e(void)
{}


// End of file
