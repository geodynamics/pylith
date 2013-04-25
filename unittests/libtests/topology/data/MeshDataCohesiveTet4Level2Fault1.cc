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
// Copyright (c) 2010-2013 University of California, Davis
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

const PylithScalar pylith::topology::MeshDataCohesiveTet4Level2Fault1::_vertices[] = {
  -1.0,  0.0, 0.0,
   0.0, -1.0, 0.0,
   0.0,  0.0, 1.0,
   0.0,  1.0, 0.0,
   1.0,  0.0, 0.0,
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
   0.0, -1.0, 0.0,
   0.0,  0.0, 1.0,
   0.0,  1.0, 0.0,
   0.0, -0.5, 0.5,
   0.0,  0.0, 0.0,
   0.0,  0.5, 0.5,
};

const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_cells[] = {
  17,  27,  24,  26,
  24,  25,  26,  27,
  24,  27,  28,  25,
  18,  28,  25,  24,
  26,  29,  27,  25,
  19,  29,  26,  25,
  25,  28,  29,  27,
  16,  27,  29,  28,
  21,  33,  30,  32,
  30,  31,  32,  33,
  30,  33,  34,  31,
  23,  34,  31,  30,
  32,  35,  33,  31,
  22,  35,  32,  31,
  31,  34,  35,  33,
  20,  33,  35,  34,
};
const int pylith::topology::MeshDataCohesiveTet4Level2Fault1::_cellsCohesive[] = {
  18,  24,  25,  22,  32,  31,  37,  39,  41,
  24,  26,  25,  32,  30,  31,  39,  40,  41,
  17,  26,  24,  21,  30,  32,  36,  40,  39,
  19,  25,  26,  23,  31,  30,  38,  41,  40,
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
 16, 17, 21, 27,
 18, 20, 22, 35,
 16, 20,
 17, 18, 19, 21, 22, 23, 24, 25, 26, 30, 31, 32, 36, 37, 38, 39, 40, 41,
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

pylith::topology::MeshDataCohesiveTet4Level2Fault1::~MeshDataCohesiveTet4Level2Fault1(void)
{}


// End of file
