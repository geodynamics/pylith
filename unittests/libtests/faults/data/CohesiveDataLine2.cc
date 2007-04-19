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
 *   0 -------- 1 -------- 2
 *
 * After adding cohesive elements
 *   0 -------- 1 -- 3 -------- 2
 */

#include "CohesiveDataLine2.hh"

const int pylith::faults::CohesiveDataLine2::_numVertices = 4;

const int pylith::faults::CohesiveDataLine2::_spaceDim = 1;

const int pylith::faults::CohesiveDataLine2::_numCells = 3;

const int pylith::faults::CohesiveDataLine2::_cellDim = 1;

const double pylith::faults::CohesiveDataLine2::_vertices[] = {
  -1.0,
   0.0,
   1.0,
   0.0
};

const int pylith::faults::CohesiveDataLine2::_numCorners[] = {
  2,
  2,
  2
};

const int pylith::faults::CohesiveDataLine2::_cells[] = {
       0,  1,
       3,  2,
       1,  3,
};

const int pylith::faults::CohesiveDataLine2::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataLine2::_numGroups = 2;

const int pylith::faults::CohesiveDataLine2::_groupSizes[] = 
  { 1, 2 };

const int pylith::faults::CohesiveDataLine2::_groups[] = {
  1,
  0, 1
};

const char* pylith::faults::CohesiveDataLine2::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataLine2::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataLine2::_filename = "data/meshLine.txt";

pylith::faults::CohesiveDataLine2::CohesiveDataLine2(void)
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

pylith::faults::CohesiveDataLine2::~CohesiveDataLine2(void)
{}


// End of file
