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
 *   Cells are 0-1, vertices are 2-4.
 *
 *   2 -------- 3 -------- 4
 *
 * After adding cohesive elements
 *   2 -------- 3 -6- 5 -------- 4
 */

#include "CohesiveLagrangeDataLine2.hh"

const int pylith::faults::CohesiveLagrangeDataLine2::_numVertices = 5;

const int pylith::faults::CohesiveLagrangeDataLine2::_spaceDim = 1;

const int pylith::faults::CohesiveLagrangeDataLine2::_numCells = 3;

const int pylith::faults::CohesiveLagrangeDataLine2::_cellDim = 1;

const double pylith::faults::CohesiveLagrangeDataLine2::_vertices[] = {
  -1.0,
   0.0,
   1.0,
   0.0,
   0.0
};

const int pylith::faults::CohesiveLagrangeDataLine2::_numCorners[] = {
  2,
  2,
  3
};

const int pylith::faults::CohesiveLagrangeDataLine2::_cells[] = {
       2,  3,
       5,  4,
       3,  5,  6
};

const int pylith::faults::CohesiveLagrangeDataLine2::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveLagrangeDataLine2::_numGroups = 2;

const int pylith::faults::CohesiveLagrangeDataLine2::_groupSizes[] = 
  { 3, 4 };

const int pylith::faults::CohesiveLagrangeDataLine2::_groups[] = {
  3, 5, 6,
  2, 3, 5, 6
};

const char* pylith::faults::CohesiveLagrangeDataLine2::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveLagrangeDataLine2::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveLagrangeDataLine2::_filename = "data/meshLine.txt";

pylith::faults::CohesiveLagrangeDataLine2::CohesiveLagrangeDataLine2(void)
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

pylith::faults::CohesiveLagrangeDataLine2::~CohesiveLagrangeDataLine2(void)
{}


// End of file
