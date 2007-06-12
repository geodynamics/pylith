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
 *       3 -------- 5 -- 9 -------- 7
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       |          |    |          |
 *       2 -------- 4 -- 8 -------- 6
 */

#include "CohesiveDataQuad4b.hh"

const int pylith::faults::CohesiveDataQuad4b::_numVertices = 8;

const int pylith::faults::CohesiveDataQuad4b::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4b::_numCells = 3;

const int pylith::faults::CohesiveDataQuad4b::_cellDim = 2;

const double pylith::faults::CohesiveDataQuad4b::_vertices[] = {
  -2.0, -1.0,
  -2.0,  1.0,
   0.0, -1.0,
   0.0,  1.0,
   2.0, -1.0,
   2.0,  1.0,
   0.0, -1.0,
   0.0,  1.0,
};

const int pylith::faults::CohesiveDataQuad4b::_numCorners[] = {
  4,
  4,
  4
};

const int pylith::faults::CohesiveDataQuad4b::_cells[] = {
  4,  5,  3,  2,
  7,  9,  8,  6,
  4,  5,  8,  9,
};

const int pylith::faults::CohesiveDataQuad4b::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataQuad4b::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4b::_groupSizes[] = 
  { 4, 4 };

const int pylith::faults::CohesiveDataQuad4b::_groups[] = {
  4, 5, 8, 9,
  3, 5, 7, 9
};

const char* pylith::faults::CohesiveDataQuad4b::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataQuad4b::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataQuad4b::_filename = 
  "data/quad4b.mesh";

pylith::faults::CohesiveDataQuad4b::CohesiveDataQuad4b(void)
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

pylith::faults::CohesiveDataQuad4b::~CohesiveDataQuad4b(void)
{}


// End of file
