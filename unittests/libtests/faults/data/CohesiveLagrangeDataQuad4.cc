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
 * Cells are 0-1, vertices are 2-7.
 *
 *       3 -------- 5 -------- 7
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       2 -------- 4 -------- 6
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,10 vertices are 2-9.
 *
 *       3 -------- 5 -11- 10 -------- 7
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       2 -------- 4 --9-- 8 -------- 6
 */

#include "CohesiveLagrangeDataQuad4.hh"

const int pylith::faults::CohesiveLagrangeDataQuad4::_numVertices = 10;

const int pylith::faults::CohesiveLagrangeDataQuad4::_spaceDim = 2;

const int pylith::faults::CohesiveLagrangeDataQuad4::_numCells = 3;

const int pylith::faults::CohesiveLagrangeDataQuad4::_cellDim = 2;

const double pylith::faults::CohesiveLagrangeDataQuad4::_vertices[] = {
  -2.0, -1.0,
  -2.0,  1.0,
   0.0, -1.0,
   0.0,  1.0,
   2.0, -1.0,
   2.0,  1.0,
   0.0, -1.0,
   0.0, -1.0,
   0.0,  1.0,
   0.0,  1.0,
};

const int pylith::faults::CohesiveLagrangeDataQuad4::_numCorners[] = {
  4,
  4,
  6
};

const int pylith::faults::CohesiveLagrangeDataQuad4::_cells[] = {
  2,  4,  5,  3,
  6,  7, 10,  8,
  4,  5,  8, 10,  9, 11
};

const int pylith::faults::CohesiveLagrangeDataQuad4::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveLagrangeDataQuad4::_numGroups = 2;

const int pylith::faults::CohesiveLagrangeDataQuad4::_groupSizes[] = 
  { 6, 5 };

const int pylith::faults::CohesiveLagrangeDataQuad4::_groups[] = {
  4, 5, 8, 9, 10, 11,
  3, 5, 7, 10, 11
};

const char* pylith::faults::CohesiveLagrangeDataQuad4::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveLagrangeDataQuad4::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveLagrangeDataQuad4::_filename = "data/meshQuad4A.txt";

pylith::faults::CohesiveLagrangeDataQuad4::CohesiveLagrangeDataQuad4(void)
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

pylith::faults::CohesiveLagrangeDataQuad4::~CohesiveLagrangeDataQuad4(void)
{}


// End of file
