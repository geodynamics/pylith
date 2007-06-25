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

#include "CohesiveDataTet4b.hh"

const int pylith::faults::CohesiveDataTet4b::_numVertices = 8;

const int pylith::faults::CohesiveDataTet4b::_spaceDim = 3;

const int pylith::faults::CohesiveDataTet4b::_numCells = 3;

const int pylith::faults::CohesiveDataTet4b::_cellDim = 3;

const double pylith::faults::CohesiveDataTet4b::_vertices[] = {
  -1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0,
   1.0,  0.0,  0.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  1.0,
   0.0,  1.0,  0.0
};

const int pylith::faults::CohesiveDataTet4b::_numCorners[] = {
  4,
  4,
  6
};

const int pylith::faults::CohesiveDataTet4b::_cells[] = {
  2,  7,  9,  8,
  3,  6,  5,  4,
  4,  5,  3,  8,  9,  7
};

const int pylith::faults::CohesiveDataTet4b::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTet4b::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4b::_groupSizes[] = 
  { 6, 5 };

const int pylith::faults::CohesiveDataTet4b::_groups[] = {
  3, 4, 5, 7, 8, 9,
  2, 4, 5, 8, 9
};

const char* pylith::faults::CohesiveDataTet4b::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataTet4b::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4b::_filename = "data/tet4b.mesh";

pylith::faults::CohesiveDataTet4b::CohesiveDataTet4b(void)
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

pylith::faults::CohesiveDataTet4b::~CohesiveDataTet4b(void)
{}


// End of file
