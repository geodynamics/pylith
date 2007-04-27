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
 *       2,3,4,5 -------- 6,7,8,9 -- 14,16,18,20 -------- 10,11,12,13
 *                                    15,17,19,21
 *                        ^^^^^^^^^^^^^^^^^^^^^^ Cohesive element
 *
 */

#include "CohesiveLagrangeDataHex8.hh"

const int pylith::faults::CohesiveLagrangeDataHex8::_numVertices = 20;

const int pylith::faults::CohesiveLagrangeDataHex8::_spaceDim = 3;

const int pylith::faults::CohesiveLagrangeDataHex8::_numCells = 3;

const int pylith::faults::CohesiveLagrangeDataHex8::_cellDim = 3;

const double pylith::faults::CohesiveLagrangeDataHex8::_vertices[] = {
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
   0.0, -1.0, -1.0,
   0.0,  1.0, -1.0,
   0.0,  1.0, -1.0,
   0.0, -1.0,  1.0,
   0.0, -1.0,  1.0,
   0.0,  1.0,  1.0,
   0.0,  1.0,  1.0
};

const int pylith::faults::CohesiveLagrangeDataHex8::_numCorners[] = {
  8,
  8,
  12
};

const int pylith::faults::CohesiveLagrangeDataHex8::_cells[] = {
  2,  4,  5,  3,  6,  8,  9,  7,
 14, 18, 20, 16, 10, 12, 13, 11,
  6,  8,  9,  7, 14, 18, 20, 16, 15, 19, 21, 17
};

const int pylith::faults::CohesiveLagrangeDataHex8::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveLagrangeDataHex8::_numGroups = 2;

const int pylith::faults::CohesiveLagrangeDataHex8::_groupSizes[] = 
  { 12, 10 };

const int pylith::faults::CohesiveLagrangeDataHex8::_groups[] = {
  6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21,
  4, 5, 8, 9, 12, 13, 18, 19, 20, 21
};

const char* pylith::faults::CohesiveLagrangeDataHex8::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveLagrangeDataHex8::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveLagrangeDataHex8::_filename = "data/meshHex8A.txt";

pylith::faults::CohesiveLagrangeDataHex8::CohesiveLagrangeDataHex8(void)
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

pylith::faults::CohesiveLagrangeDataHex8::~CohesiveLagrangeDataHex8(void)
{}


// End of file
