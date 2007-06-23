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
 * Cells are 0-3, vertices are 4-12.
 *
 *      10 --------11 --------12
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       5 -------- 7 -------- 9
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       4 -------- 6 -------- 8
 *
 * After adding cohesive elements
 *
 * Cells are 0-3,16-17 vertices are 4-15.
 *
 *      10 --------11--15 --------12
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       5 -------- 7--14 -------- 9
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       4 -------- 6--13 -------- 8
 */

#include "CohesiveDataQuad4f.hh"

const int pylith::faults::CohesiveDataQuad4f::_numVertices = 12;

const int pylith::faults::CohesiveDataQuad4f::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4f::_numCells = 6;

const int pylith::faults::CohesiveDataQuad4f::_cellDim = 2;

const double pylith::faults::CohesiveDataQuad4f::_vertices[] = {
  -2.0, -1.0,
  -2.0,  1.0,
   0.0, -1.0,
   0.0,  1.0,
   2.0, -1.0,
   2.0,  1.0,
  -2.0,  3.0,
   0.0,  3.0,
   2.0,  3.0,
   0.0, -1.0,
   0.0,  1.0,
   0.0,  3.0,
};

const int pylith::faults::CohesiveDataQuad4f::_numCorners[] = {
  4,
  4,
  4,
  4,
  4,
  4,
};

const int pylith::faults::CohesiveDataQuad4f::_cells[] = {
  5,  7, 11, 10,
 14,  9, 12, 15,
 13,  8,  9, 14,
  4,  6,  7,  5,
  6,  7, 13, 14,
  7, 11, 14, 15,
};

const int pylith::faults::CohesiveDataQuad4f::_materialIds[] = {
  0,  0,  0, 0,
  1,  1,
};

const int pylith::faults::CohesiveDataQuad4f::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4f::_groupSizes[] = 
  { 6, 4 };

const int pylith::faults::CohesiveDataQuad4f::_groups[] = {
  6, 7, 11, 13, 14, 15,
  5, 7, 9, 14
};

const char* pylith::faults::CohesiveDataQuad4f::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataQuad4f::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataQuad4f::_filename = 
  "data/quad4f.mesh";

pylith::faults::CohesiveDataQuad4f::CohesiveDataQuad4f(void)
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

pylith::faults::CohesiveDataQuad4f::~CohesiveDataQuad4f(void)
{}


// End of file
