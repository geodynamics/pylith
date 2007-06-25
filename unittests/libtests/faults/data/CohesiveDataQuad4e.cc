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
 *      10 --------15--11 --------12
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       5 --------14-- 7 -------- 9
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       |          |   |          |
 *       4 --------13-- 6 -------- 8
 */

#include "CohesiveDataQuad4e.hh"

const int pylith::faults::CohesiveDataQuad4e::_numVertices = 12;

const int pylith::faults::CohesiveDataQuad4e::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4e::_numCells = 6;

const int pylith::faults::CohesiveDataQuad4e::_cellDim = 2;

const double pylith::faults::CohesiveDataQuad4e::_vertices[] = {
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

const int pylith::faults::CohesiveDataQuad4e::_numCorners[] = {
  4,
  4,
  4,
  4,
  4,
  4,
};

const int pylith::faults::CohesiveDataQuad4e::_cells[] = {
  4, 13, 14,  5,
  6,  8,  9,  7,
  5, 14, 15, 10,
  7,  9, 12, 11,
  7,  6, 14, 13,
 11,  7, 15, 14,
};

const int pylith::faults::CohesiveDataQuad4e::_materialIds[] = {
  0,  0,  0, 0,
  1,  1,
};

const int pylith::faults::CohesiveDataQuad4e::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4e::_groupSizes[] = 
  { 6, 4 };

const int pylith::faults::CohesiveDataQuad4e::_groups[] = {
  6, 7, 11, 13, 14, 15,
  5, 7, 9, 14
};

const char* pylith::faults::CohesiveDataQuad4e::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataQuad4e::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataQuad4e::_filename = 
  "data/quad4e.mesh";

pylith::faults::CohesiveDataQuad4e::CohesiveDataQuad4e(void)
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

pylith::faults::CohesiveDataQuad4e::~CohesiveDataQuad4e(void)
{}


// End of file
