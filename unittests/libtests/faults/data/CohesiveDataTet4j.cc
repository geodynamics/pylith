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
 * Cells are 0-5, vertices are 6-17.
 */

#include "CohesiveDataTet4j.hh"

const int pylith::faults::CohesiveDataTet4j::_numVertices = 12;

const int pylith::faults::CohesiveDataTet4j::_spaceDim = 3;

const int pylith::faults::CohesiveDataTet4j::_numCells = 7;

const int pylith::faults::CohesiveDataTet4j::_cellDim = 3;

const double pylith::faults::CohesiveDataTet4j::_vertices[] = {
  -2.0, -1.0,  0.0,
  -2.0,  0.0,  0.0,
  -2.0,  0.0,  1.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  0.0,
   0.0,  0.0,  1.0,
   2.0, -1.0,  0.0,
   2.0,  0.0,  0.0,
   2.0,  0.0,  1.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  0.0,
   0.0,  0.0,  1.0,
};

const int pylith::faults::CohesiveDataTet4j::_numCorners[] = {
  4,
  4,
  4,
  4,
  4,
  4,
  6
};

const int pylith::faults::CohesiveDataTet4j::_cells[] = {
  7,  9,  8,  6,
  17, 15, 13, 14,
  16, 13, 17, 15,
  10,  9, 11,  7,
  15, 13, 14, 12,
  7, 11,  8,   9,
  11, 9, 10, 17, 15, 16
};

const int pylith::faults::CohesiveDataTet4j::_materialIds[] = {
  0, 2, 2, 0, 2, 0,
  1, 1
};

const int pylith::faults::CohesiveDataTet4j::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4j::_groupSizes[] = 
  { 6, 8 };

const int pylith::faults::CohesiveDataTet4j::_groups[] = {
  9, 10, 11, 15, 16, 17,
  6, 8, 9, 11, 12, 14, 15, 17
};

const char* pylith::faults::CohesiveDataTet4j::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataTet4j::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4j::_filename = "data/tet4j.mesh";

pylith::faults::CohesiveDataTet4j::CohesiveDataTet4j(void)
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

pylith::faults::CohesiveDataTet4j::~CohesiveDataTet4j(void)
{}


// End of file
