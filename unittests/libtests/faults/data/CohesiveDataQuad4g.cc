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
 *       7 --------12 --------15
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       6 --------11 --------14
 *       |        / |          |
 *       |       /  |          |
 *       |      /   |          |
 *       5-----9    |          |
 *             |    |          |
 *             |    |          |
 *             |    |          |
 *             |    |          |
 *             8---10 --------13
 *
 * After adding cohesive elements
 *
 * Cells are 0-3,16-17 vertices are 4-15.
 *
 *       7 --------12 -- 18 --------15
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       |          |     |          |
 *       6 --------11 -- 17 --------14
 *       |        / |     |          |
 *       |       /  |     |          |
 *       |      /   |     |          |
 *       5-----9    |     |          |
 *             |    |     |          |
 *             |    |     |          |
 *             |    |     |          |
 *             |    |     |          |
 *             8---10 -- 16 --------13
 */

#include "CohesiveDataQuad4g.hh"

const int pylith::faults::CohesiveDataQuad4g::_numVertices = 14;

const int pylith::faults::CohesiveDataQuad4g::_spaceDim = 2;

const int pylith::faults::CohesiveDataQuad4g::_numCells = 7;

const int pylith::faults::CohesiveDataQuad4g::_cellDim = 2;

const double pylith::faults::CohesiveDataQuad4g::_vertices[] = {
  -2.0, -2.0,
  -2.0,  0.0,
  -2.0,  2.0,
  -1.0, -2.0,
  -1.0, -1.0,
   0.0, -2.0,
   0.0,  0.0,
   0.0,  2.0,
   2.0, -2.0,
   2.0,  0.0,
   2.0,  2.0,
   0.0, -2.0,
   0.0,  0.0,
   0.0,  2.0,
};

const int pylith::faults::CohesiveDataQuad4g::_numCorners[] = {
  4,
  4,
  4,
  4,
  4,
  4,
  4,
};

const int pylith::faults::CohesiveDataQuad4g::_cells[] = {
  6, 17, 18,  7,
 11, 10, 13, 14,
 15, 12, 11, 14,
  6,  5,  9, 17,
  9,  8, 16, 17,
 11, 10, 17, 16,
 12, 11, 18, 17,
};

const int pylith::faults::CohesiveDataQuad4g::_materialIds[] = {
  0,  0,  0,  2,  2,
  1,  1,
};

const int pylith::faults::CohesiveDataQuad4g::_numGroups = 2;

const int pylith::faults::CohesiveDataQuad4g::_groupSizes[] = 
  { 6, 3 };

const int pylith::faults::CohesiveDataQuad4g::_groups[] = {
  10, 11, 12, 16, 17, 18,
  13, 14, 15
};

const char* pylith::faults::CohesiveDataQuad4g::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataQuad4g::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataQuad4g::_filename = 
  "data/quad4g.mesh";

pylith::faults::CohesiveDataQuad4g::CohesiveDataQuad4g(void)
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

pylith::faults::CohesiveDataQuad4g::~CohesiveDataQuad4g(void)
{}


// End of file
