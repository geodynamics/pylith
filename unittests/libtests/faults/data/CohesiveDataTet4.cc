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
 * Cells are 0-1, vertices are 2-6.
 *
 * 2   3,4,5  6
 *
 *     ^^^^^ Face in x-y plane
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,10, vertices are 2-9.
 *
 * 2   3,4,5  7,8,9   6
 *
 *     ^^^^^^^^^^^^ Cohesive element in x-y plane.
 */

#include "CohesiveDataTet4.hh"

const int pylith::faults::CohesiveDataTet4::_numVertices = 8;

const int pylith::faults::CohesiveDataTet4::_spaceDim = 3;

const int pylith::faults::CohesiveDataTet4::_numCells = 3;

const int pylith::faults::CohesiveDataTet4::_cellDim = 3;

const double pylith::faults::CohesiveDataTet4::_vertices[] = {
  -1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0,
   1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0
};

const int pylith::faults::CohesiveDataTet4::_numCorners[] = {
  4,
  4,
  6
};

const int pylith::faults::CohesiveDataTet4::_cells[] = {
  3,  4,  5,  2,
  7,  9,  8,  6,
  3,  5,  4,  7,  9,  8
};

const int pylith::faults::CohesiveDataTet4::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTet4::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4::_groupSizes[] = 
  { 6, 5 };

const int pylith::faults::CohesiveDataTet4::_groups[] = {
  3, 4, 5, 7, 8, 9,
  2, 4, 5, 8, 9
};

const char* pylith::faults::CohesiveDataTet4::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataTet4::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4::_filename = "data/tet4.mesh";

pylith::faults::CohesiveDataTet4::CohesiveDataTet4(void)
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

pylith::faults::CohesiveDataTet4::~CohesiveDataTet4(void)
{}


// End of file
