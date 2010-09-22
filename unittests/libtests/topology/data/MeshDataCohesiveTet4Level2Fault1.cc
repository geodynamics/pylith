// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "MeshDataCohesiveTet4Level2Fault1.hh"

const char* pylith::topology::MeshDataCohesiveTet4Level2Fault1::_filename = 
  "data/twotet4.mesh";

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_refineLevel = 2;
const char* pylith::topology::MeshDataCohesiveTet4Level2Fault1::_faultA = 
  "fault";
const char* pylith::topology::MeshDataCohesiveTet4Level2Fault1::_faultB = 0;

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_numVertices = 26;

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_spaceDim = 3;

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_numCells = 16;

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_numCellsCohesive = 4;

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_cellDim = 3;

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_numCorners = 4;

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_numCornersCohesive = 9;

const double pylith::topology::MeshDataCohesiveTet4Level2Fault1::_vertices[] = {
  -1.0,  0.0, 0.0,
   0.0, -1.0, 0.0,
   0.0,  0.0, 1.0,
   0.0,  1.0, 0.0,
   1.0,  0.0, 0.0,
   0.0, -1.0, 0.0,
   0.0,  0.0, 1.0,
   0.0,  1.0, 0.0,
   0.0, -1.0, 0.0,
   0.0,  0.0, 1.0,
   0.0,  1.0, 0.0,
   0.0, -0.5, 0.5,
   0.0,  0.5, 0.5,
   0.0,  0.0, 0.0,
  -0.5, -0.5, 0.0,
  -0.5,  0.0, 0.5,
  -0.5,  0.5, 0.0,
   0.0,  0.0, 0.0,
   0.0,  0.5, 0.5,
   0.0, -0.5, 0.5,
   0.5, -0.5, 0.0,
   0.5,  0.5, 0.0,
   0.5,  0.0, 0.5,
   0.0, -0.5, 0.5,
   0.0,  0.0, 0.0,
   0.0,  0.5, 0.5,
};

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_cells[] = {
   1,  14,  11,  13, // original cell 1
  11,  12,  13,  14,
  11,  14,  15,  12,
   2,  15,  12,  11,
  13,  16,  14,  12,
   3,  16,  13,  12,
  12,  15,  16,  14,
   0,  14,  16,  15,
   5,  20,  17,  19, // original cell 2
  17,  18,  19,  20,
  17,  20,  21,  18,
   7,  21,  18,  17,
  19,  22,  20,  18,
   6,  22,  19,  18,
  18,  21,  22,  20,
   4,  20,  22,  21,
};
const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_cellsCohesive[] = {
   2,  11,  12,   6,  19,  18,   9,  23,  25,
  11,  13,  12,  19,  17,  18,  23,  24,  25,
   1,  13,  11,   5,  17,  19,   8,  24,  23,
   3,  12,  13,   7,  18,  17,  10,  25,  24,
};
const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_materialIds[] = {
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
  100, 100, 100, 100,
};

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_numGroups = 4;

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_groupSizes[] = {
  4, 4, 2, 18,
};

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_groups[] = {
  0, 1, 5, 14,
  2, 4, 6, 22,
  0, 4,
  1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18, 19, 23, 24, 25,
};

const char* pylith::topology::MeshDataCohesiveTet4Level2Fault1::_groupNames[] = {
  "edge 1",
  "edge 2",
  "end points",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveTet4Level2Fault1::_groupTypes[] = {
  "vertex", 
  "vertex",
  "vertex",
  "vertex",
};

pylith::topology::MeshDataCohesiveTet4Level2Fault1::MeshDataCohesiveTet4Level2Fault1(void)
{ // constructor
  filename = const_cast<char*>(_filename);
  refineLevel = _refineLevel;
  faultA = const_cast<char*>(_faultA);
  faultB = const_cast<char*>(_faultB);

  numVertices = _numVertices;
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numCells = _numCells;
  numCorners = _numCorners;
  numCellsCohesive = _numCellsCohesive;
  numCornersCohesive = _numCornersCohesive;
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  cellsCohesive = const_cast<int*>(_cellsCohesive);
  materialIds = const_cast<int*>(_materialIds);
  groups = const_cast<int*>(_groups);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
} // constructor

pylith::topology::MeshDataCohesiveTet4Level2Fault1::~MeshDataCohesiveTet4Level2Fault1(void)
{}


// End of file
