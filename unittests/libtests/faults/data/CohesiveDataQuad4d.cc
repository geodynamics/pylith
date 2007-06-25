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
 *       3 -------- 9 -- 5 -------- 7
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       2 -------- 8 -- 4 -------- 6
 */

#include "CohesiveDataQuad4d.hh"

const int pylith::faults::CohesiveDataQuad4d::_numVertices = 8;

const int pylith::faults::CohesiveDataQuad4d::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4d::_numCells = 3;

const int pylith::faults::CohesiveDataQuad4d::_cellDim = 2;

const double pylith::faults::CohesiveDataQuad4d::_vertices[] = {
  -2.0, -1.0,
  -2.0,  1.0,
   0.0, -1.0,
   0.0,  1.0,
   2.0, -1.0,
   2.0,  1.0,
   0.0, -1.0,
   0.0,  1.0,
};

const int pylith::faults::CohesiveDataQuad4d::_numCorners[] = {
  4,
  4,
  4
};

const int pylith::faults::CohesiveDataQuad4d::_cells[] = {
  3,  2,  8,  9,
  4,  6,  7,  5,
  5,  4,  9,  8,
};

const int pylith::faults::CohesiveDataQuad4d::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataQuad4d::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4d::_groupSizes[] = 
  { 4, 4 };

const int pylith::faults::CohesiveDataQuad4d::_groups[] = {
  4, 5, 8, 9,
  3, 5, 7, 9
};

const char* pylith::faults::CohesiveDataQuad4d::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataQuad4d::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataQuad4d::_filename = 
  "data/quad4d.mesh";

pylith::faults::CohesiveDataQuad4d::CohesiveDataQuad4d(void)
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

pylith::faults::CohesiveDataQuad4d::~CohesiveDataQuad4d(void)
{}


// End of file
