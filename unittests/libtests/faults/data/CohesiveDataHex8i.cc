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
 * Cells are 0-6 and vertices are 8-43
 */

#include "CohesiveDataHex8i.hh"

const int pylith::faults::CohesiveDataHex8i::_numVertices = 36;

const int pylith::faults::CohesiveDataHex8i::_spaceDim = 3;

const int pylith::faults::CohesiveDataHex8i::_numCells = 9;

const int pylith::faults::CohesiveDataHex8i::_cellDim = 3;

const double pylith::faults::CohesiveDataHex8i::_vertices[] = {
  -2.0, -2.0, -2.0,
  -2.0, -1.0, -2.0,
  -3.0,  0.0, -2.0,
  -2.0,  1.0, -2.0,
  -2.0,  2.0, -2.0,
  -2.0, -2.0,  0.0,
  -2.0, -1.0,  0.0,
  -3.0,  0.0,  0.0,
  -2.0,  1.0,  0.0,
  -2.0,  2.0,  0.0,
  -2.0, -1.0,  2.0,
  -3.0,  0.0,  2.0,
  -2.0,  1.0,  2.0,
   0.0, -2.0, -2.0,
   0.0,  0.0, -2.0,
   0.0,  2.0, -2.0,
   0.0, -2.0,  0.0,
   0.0,  0.0,  0.0,
   0.0,  2.0,  0.0,
   0.0,  0.0,  2.0,
   2.0, -2.0, -2.0,
   2.0, -1.0, -2.0,
   3.0,  0.0, -2.0,
   2.0,  1.0, -2.0,
   2.0,  2.0, -2.0,
   2.0, -2.0,  0.0,
   2.0, -1.0,  0.0,
   3.0,  0.0,  0.0,
   2.0,  1.0,  0.0,
   2.0,  2.0,  0.0,
   0.0, -2.0, -2.0,
   0.0,  0.0, -2.0,
   0.0,  2.0, -2.0,
   0.0, -2.0,  0.0,
   0.0,  0.0,  0.0,
   0.0,  2.0,  0.0,
};

const int pylith::faults::CohesiveDataHex8i::_numCorners[] = {
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
};

const int pylith::faults::CohesiveDataHex8i::_cells[] = {
    8,  7, 20, 21, 13, 12, 23, 24,
   14,  9, 10, 15, 13,  8, 21, 24,
   15, 10, 11, 16, 24, 21, 22, 25,
   30, 38, 28, 29, 35, 41, 33, 34,
   41, 35, 30, 38, 42, 36, 31, 39,
   27, 28, 38, 37, 32, 33, 41, 40,
   15, 14, 13, 24, 19, 18, 17, 26,
   24, 23, 20, 21, 41, 40, 37, 38,
   24, 21, 22, 25, 41, 38, 39, 42
};

const int pylith::faults::CohesiveDataHex8i::_materialIds[] = {
  0,  0,  0,  0,  0,  0,  2,
  1,  1
};

const int pylith::faults::CohesiveDataHex8i::_numGroups = 2;

const int pylith::faults::CohesiveDataHex8i::_groupSizes[] = 
  { 12, 10 };

const int pylith::faults::CohesiveDataHex8i::_groups[] = {
  20, 21, 22, 23, 24, 25, 37, 38, 39, 40, 41, 42,
  27, 28, 29, 30, 31, 32, 33, 34, 35, 36
};

const char* pylith::faults::CohesiveDataHex8i::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataHex8i::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataHex8i::_filename = 
  "data/hex8i.mesh";

pylith::faults::CohesiveDataHex8i::CohesiveDataHex8i(void)
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

pylith::faults::CohesiveDataHex8i::~CohesiveDataHex8i(void)
{}


// End of file
