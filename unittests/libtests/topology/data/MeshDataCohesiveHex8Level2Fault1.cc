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
  -2.0, -1.0, -1.0, // 21
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
   0.0, -1.0, -1.0,
   0.0, +1.0, -1.0,
   0.0, -1.0, +1.0,
   0.0, +1.0, +1.0,

  -1.0, -1.0, -1.0, // 41 (edges)
   0.0,  0.0, -1.0,
  -1.0, +1.0, -1.0,
  -2.0,  0.0, -1.0,
  -1.0, -1.0, +1.0, // 45
   0.0,  0.0, +1.0,
  -1.0, +1.0, +1.0,
  -2.0,  0.0, +1.0,
  -2.0, -1.0, +0.0, // 49
  +0.0, -1.0, +0.0,
  +0.0, +1.0, +0.0,
  -2.0, +1.0, +0.0,
  -1.0, -1.0, +0.0, // 53 (faces)
  +0.0, +0.0, +0.0,
  -1.0, +1.0, +0.0,
  -2.0, +0.0, +0.0,
  -1.0, +0.0, -1.0,
  -1.0, +0.0, +1.0,
  -1.0, +0.0, +0.0, // 59 (volume)
  +2.0, +0.0, -1.0, // 60 (edges)
  +1.0, +1.0, -1.0,
  +1.0, -1.0, -1.0,
  +2.0, +0.0, +1.0, // 63
  +1.0, +1.0, +1.0,
  +1.0, -1.0, +1.0,
  +2.0, -1.0, +0.0, // 66
  +2.0, +1.0, +0.0,
  +2.0, +0.0, +0.0, // 68 (faces)
  +1.0, +1.0, +0.0,
  +1.0, -1.0, +0.0,
  +1.0, +0.0, -1.0,
  +1.0, +0.0, +1.0,
  +1.0, +0.0, +0.0, // 73 (volume)
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
};

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_cells[] = {
  16,   28,   44,   31,   36,   40,   46,   43,
  20,   29,   44,   28,   37,   41,   46,   40, 
  21,   30,   44,   29,   38,   42,   46,   41, 
  17,   31,   44,   30,   39,   43,   46,   42, 
  36,   40,   46,   43,   18,   32,   45,   35, 
  37,   41,   46,   40,   22,   33,   45,   32, 
  38,   42,   46,   41,   23,   34,   45,   33, 
  39,   43,   46,   42,   19,   35,   45,   34, 
  24,   47,   58,   49,   53,   55,   60,   57, 
  25,   48,   58,   47,   54,   56,   60,   55, 
  21,   29,   58,   48,   38,   41,   60,   56, 
  20,   49,   58,   29,   37,   57,   60,   41, 
  53,   55,   60,   57,   26,   50,   59,   52, 
  54,   56,   60,   55,   27,   51,   59,   50, 
  38,   41,   60,   56,   23,   33,   59,   51, 
  37,   57,   60,   41,   22,   52,   59,   33, 
};
const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_cellsCohesive[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_materialIds[] = {
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
  10, 10, 10, 10,
};

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_numGroups = 4;

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_groupSizes[] = {
  2, 9, 9, 9,
};

const int pylith::topology::MeshDataCohesiveHex8Level2Fault1::_groups[] = {
  16, 24,
  16, 17, 18, 19, 31, 35, 36, 39, 43,
  20, 22, 24, 26, 37, 49, 52, 53, 57,
  20, 21, 22, 23, 29, 33, 37, 38, 41,
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
