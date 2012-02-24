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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "MeshDataCohesiveHex8Level2Fault1.hh"

const char* pylith::topology::MeshDataCohesiveHex8Level2Fault1::_filename = 
  "data/twohex8.mesh";

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_refineLevel = 2;
const char* pylith::topology::MeshDataCohesiveHex8Level2Fault1::_faultA = 
  "fault";
const char* pylith::topology::MeshDataCohesiveHex8Level2Fault1::_faultB = 0;

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_numVertices = 63;

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_spaceDim = 3;

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_numCells = 16;

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_numCellsCohesive = 4;

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_cellDim = 3;

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_numCorners = 8;

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_numCornersCohesive = 12;

const double pylith::topology::MeshDataCohesiveHex8Level2Fault1::_vertices[] = {
  -2.0, -1.0, -1.0, // 16
  -2.0, +1.0, -1.0,
  -2.0, -1.0, +1.0,
  -2.0, +1.0, +1.0,
   0.0, -1.0, -1.0,
   0.0, +1.0, -1.0,
   0.0, -1.0, +1.0,
   0.0, +1.0, +1.0,
  +2.0, -1.0, -1.0,
  +2.0, +1.0, -1.0,
  +2.0, -1.0, +1.0,
  +2.0, +1.0, +1.0,
   0.0, -1.0, -1.0,
   0.0, +1.0, -1.0,
   0.0, -1.0, +1.0,
   0.0, +1.0, +1.0,
  -1.0, -1.0, -1.0, // 32 (edges)
   0.0,  0.0, -1.0,
  -1.0, +1.0, -1.0,
  -2.0,  0.0, -1.0,
  -1.0, -1.0, +1.0, // 36
   0.0,  0.0, +1.0,
  -1.0, +1.0, +1.0,
  -2.0,  0.0, +1.0,
  -2.0, -1.0, +0.0, // 40
  +0.0, -1.0, +0.0,
  +0.0, +1.0, +0.0,
  -2.0, +1.0, +0.0,
  -1.0, -1.0, +0.0, // 44 (faces)
  +0.0, +0.0, +0.0,
  -1.0, +1.0, +0.0,
  -2.0, +0.0, +0.0,
  -1.0, +0.0, -1.0,
  -1.0, +0.0, +1.0,
  -1.0, +0.0, +0.0, // 50 (volume)
  +2.0, +0.0, -1.0, // 51 (edges)
  +1.0, +1.0, -1.0,
  +0.0,  0.0, -1.0,
  +1.0, -1.0, -1.0,
  +2.0, +0.0, +1.0, // 55
  +1.0, +1.0, +1.0,
  +0.0,  0.0, +1.0,
  +1.0, -1.0, +1.0,
  +2.0, -1.0, +0.0, // 59
  +2.0, +1.0, +0.0,
  +0.0, +1.0, +0.0,
  +0.0, -1.0, +0.0,
  +2.0, +0.0, +0.0, // 63 (faces)
  +1.0, +1.0, +0.0,
   0.0, +0.0, +0.0,
  +1.0, -1.0, +0.0,
  +1.0, +0.0, -1.0,
  +1.0, +0.0, +1.0,
  +1.0, +0.0, +0.0, // 69 (volume)
   0.0, -1.0, -1.0, // 70 (Lagrange vertices)
   0.0, +1.0, -1.0,
   0.0, -1.0, +1.0,
   0.0, +1.0, +1.0,
   0.0,  0.0, -1.0, // 74 (edges)
   0.0, +1.0,  0.0,
   0.0,  0.0, +1.0,
   0.0, -1.0,  0.0,
   0.0,  0.0,  0.0, // 78 (face)
};

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_cells[] = {
  16, 32, 48, 35, 40, 44, 50, 47, 
  20, 33, 48, 32, 41, 45, 50, 44, 
  21, 34, 48, 33, 42, 46, 50, 45, 
  17, 35, 48, 34, 43, 47, 50, 46, 
  40, 44, 50, 47, 18, 36, 49, 39, 
  41, 45, 50, 44, 22, 37, 49, 36, 
  42, 46, 50, 45, 23, 38, 49, 37, 
  43, 47, 50, 46, 19, 39, 49, 38, 
  24, 51, 67, 54, 59, 63, 69, 66, 
  25, 52, 67, 51, 60, 64, 69, 63, 
  29, 53, 67, 52, 61, 65, 69, 64, 
  28, 54, 67, 53, 62, 66, 69, 65, 
  59, 63, 69, 66, 26, 55, 68, 58, 
  60, 64, 69, 63, 27, 56, 68, 55, 
  61, 65, 69, 64, 31, 57, 68, 56, 
  62, 66, 69, 65, 30, 58, 68, 57, 
};
const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_cellsCohesive[] = {
  20, 33, 45, 41, 28, 53, 65, 62, 70, 74, 78, 77, 
  21, 42, 45, 33, 29, 61, 65, 53, 71, 75, 78, 74, 
  23, 37, 45, 42, 31, 57, 65, 61, 73, 76, 78, 75, 
  22, 41, 45, 37, 30, 62, 65, 57, 72, 77, 78, 76, 
};
const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_materialIds[] = {
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
  100, 100, 100, 100,
};

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_numGroups = 4;

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_groupSizes[] = {
  2, 9, 12, 27,
};

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_groups[] = {
  16, 24,
  16, 17, 18, 19, 35, 39, 40, 43, 47,
  20, 22, 24, 26, 28, 30, 41, 54, 58, 59, 62, 66,
  20, 21, 22, 23, 28, 29, 30, 31, 33, 37, 41, 42, 45, 53, 57, 61, 62, 65,
  70, 71, 72, 73, 74, 75, 76, 77, 78,
};

const char* pylith::topology::MeshDataCohesiveHex8Level2Fault1::_groupNames[] = {
  "end points",
  "face 1",
  "face 2",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveHex8Level2Fault1::_groupTypes[] = {
  "vertex",
  "vertex",
  "vertex",
  "vertex",
};

pylith::topology::MeshDataCohesiveHex8Level2Fault1::MeshDataCohesiveHex8Level2Fault1(void)
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

pylith::topology::MeshDataCohesiveHex8Level2Fault1::~MeshDataCohesiveHex8Level2Fault1(void)
{}


// End of file
