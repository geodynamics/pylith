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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "MeshDataCohesiveHex8Level2.hh"

const char* pylith::topology::MeshDataCohesiveHex8Level2::_filename = 
  "data/twohex8.mesh";

const int pylith::topology::MeshDataCohesiveHex8Level2::_refineLevel = 2;
const char* pylith::topology::MeshDataCohesiveHex8Level2::_faultA = 0;
const char* pylith::topology::MeshDataCohesiveHex8Level2::_faultB = 0;

const int pylith::topology::MeshDataCohesiveHex8Level2::_numVertices = 45;

const int pylith::topology::MeshDataCohesiveHex8Level2::_spaceDim = 3;

const int pylith::topology::MeshDataCohesiveHex8Level2::_numCells = 16;

const int pylith::topology::MeshDataCohesiveHex8Level2::_numCellsCohesive = 0;

const int pylith::topology::MeshDataCohesiveHex8Level2::_cellDim = 3;

const int pylith::topology::MeshDataCohesiveHex8Level2::_numCorners = 8;

const int pylith::topology::MeshDataCohesiveHex8Level2::_numCornersCohesive = 6;

const PylithScalar pylith::topology::MeshDataCohesiveHex8Level2::_vertices[] = {
  -2.0, -1.0, -1.0, // 26
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
  -1.0, -1.0, -1.0, // 28 (edges)
   0.0,  0.0, -1.0,
  -1.0, +1.0, -1.0,
  -2.0,  0.0, -1.0,
  -1.0, -1.0, +1.0, // 32
   0.0,  0.0, +1.0,
  -1.0, +1.0, +1.0,
  -2.0,  0.0, +1.0,
  -2.0, -1.0, +0.0, // 36
  +0.0, -1.0, +0.0,
  +0.0, +1.0, +0.0,
  -2.0, +1.0, +0.0,
  -1.0, -1.0, +0.0, // 40 (faces)
  +0.0, +0.0, +0.0,
  -1.0, +1.0, +0.0,
  -2.0, +0.0, +0.0,
  -1.0, +0.0, -1.0,
  -1.0, +0.0, +1.0,
  -1.0, +0.0, +0.0, // 46 (volume)
  +2.0, +0.0, -1.0, // 47 (edges)
  +1.0, +1.0, -1.0,
  +1.0, -1.0, -1.0,
  +2.0, +0.0, +1.0, // 50
  +1.0, +1.0, +1.0,
  +1.0, -1.0, +1.0,
  +2.0, -1.0, +0.0, // 53
  +2.0, +1.0, +0.0,
  +2.0, +0.0, +0.0, // 55 (faces)
  +1.0, +1.0, +0.0,
  +1.0, -1.0, +0.0,
  +1.0, +0.0, -1.0,
  +1.0, +0.0, +1.0,
  +1.0, +0.0, +0.0, // 60 (volume)
};

const int pylith::topology::MeshDataCohesiveHex8Level2::_cells[] = {
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
const int pylith::topology::MeshDataCohesiveHex8Level2::_cellsCohesive[] = {
};
const int pylith::topology::MeshDataCohesiveHex8Level2::_materialIds[] = {
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
};

const int pylith::topology::MeshDataCohesiveHex8Level2::_numGroups = 4;

const int pylith::topology::MeshDataCohesiveHex8Level2::_groupSizes[] = {
  2, 9, 9, 9,
};

const int pylith::topology::MeshDataCohesiveHex8Level2::_groups[] = {
  16, 24,
  16, 17, 18, 19, 31, 35, 36, 39, 43,
  20, 22, 24, 26, 37, 49, 52, 53, 57,
  20, 21, 22, 23, 29, 33, 37, 38, 41,
};

const char* pylith::topology::MeshDataCohesiveHex8Level2::_groupNames[] = {
  "end points",
  "face 1",
  "face 2",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveHex8Level2::_groupTypes[] = {
  "vertex",
  "vertex",
  "vertex",
  "vertex",
};

pylith::topology::MeshDataCohesiveHex8Level2::MeshDataCohesiveHex8Level2(void)
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
  vertices = const_cast<PylithScalar*>(_vertices);
  cells = const_cast<int*>(_cells);
  cellsCohesive = const_cast<int*>(_cellsCohesive);
  materialIds = const_cast<int*>(_materialIds);
  groups = const_cast<int*>(_groups);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
} // constructor

pylith::topology::MeshDataCohesiveHex8Level2::~MeshDataCohesiveHex8Level2(void)
{}


// End of file
