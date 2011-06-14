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

#include "MeshDataCohesiveTri3Level2Fault1.hh"

const char* pylith::topology::MeshDataCohesiveTri3Level2Fault1::_filename = 
  "data/fourtri3.mesh";

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_refineLevel = 2;
const char* pylith::topology::MeshDataCohesiveTri3Level2Fault1::_faultA = 
  "fault";
const char* pylith::topology::MeshDataCohesiveTri3Level2Fault1::_faultB = 0;

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_numVertices = 23;

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_spaceDim = 2;

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_numCells = 16;

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_numCellsCohesive = 4;

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_cellDim = 2;

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_numCorners = 3;

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_numCornersCohesive = 6;

const double pylith::topology::MeshDataCohesiveTri3Level2Fault1::_vertices[] = {
  -1.0,  0.0,
   0.0, -1.0,
   0.0,  0.0,
   0.0,  1.0,
   1.0,  0.0,
   0.0, -1.0,
   0.0,  0.0,
   0.0,  1.0,
  -0.5, -0.5,
   0.0, -0.5,
  -0.5,  0.0,
   0.0,  0.5,
  -0.5,  0.5,
   0.0, -0.5,
   0.5, -0.5,
   0.5,  0.0,
   0.5,  0.5,
   0.0,  0.5,
   0.0, -1.0,
   0.0,  0.0,
   0.0,  1.0,
   0.0, -0.5,
   0.0,  0.5,
};

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_cells[] = {
  16,  24,  26,
  24,  25,  26,
  17,  25,  24,
  18,  26,  25,
  18,  27,  26,
  27,  28,  26,
  19,  28,  27,
  16,  26,  28,
  22,  29,  31,
  29,  30,  31,
  21,  30,  29,
  20,  31,  30,
  22,  31,  33,
  31,  32,  33,
  20,  32,  31,
  23,  33,  32,
};
const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_cellsCohesive[] = {
  17,  25,  21,  29,  34,  37,
  25,  18,  29,  22,  37,  35,
  18,  27,  22,  33,  35,  38,
  27,  19,  33,  23,  38,  36,
};
const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_materialIds[] = {
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
  100, 100, 100, 100,
};

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_numGroups = 4;

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_groupSizes[] = {
  7, 4, 2, 15,
};

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_groups[] = {
  16, 17, 19, 21, 23, 24, 28,
  17, 20, 21, 30,
  16, 20,
  17, 18, 19, 21, 22, 23, 25, 27, 29, 33, 34, 35, 36, 37, 38,
};

const char* pylith::topology::MeshDataCohesiveTri3Level2Fault1::_groupNames[] = {
  "edge 1",
  "edge 2",
  "end points",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveTri3Level2Fault1::_groupTypes[] = {
  "vertex", 
  "vertex",
  "vertex",
  "vertex",
};

pylith::topology::MeshDataCohesiveTri3Level2Fault1::MeshDataCohesiveTri3Level2Fault1(void)
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

pylith::topology::MeshDataCohesiveTri3Level2Fault1::~MeshDataCohesiveTri3Level2Fault1(void)
{}


// End of file
