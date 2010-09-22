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

#include "MeshDataCohesiveTet4Level2.hh"

const char* pylith::topology::MeshDataCohesiveTet4Level2::_filename = 
  "data/twotet4.mesh";

const int pylith::topology::MeshDataCohesiveTet4Level2::_refineLevel = 2;
const char* pylith::topology::MeshDataCohesiveTet4Level2::_faultA = 0;
const char* pylith::topology::MeshDataCohesiveTet4Level2::_faultB = 0;

const int pylith::topology::MeshDataCohesiveTet4Level2::_numVertices = 14;

const int pylith::topology::MeshDataCohesiveTet4Level2::_spaceDim = 3;

const int pylith::topology::MeshDataCohesiveTet4Level2::_numCells = 16;

const int pylith::topology::MeshDataCohesiveTet4Level2::_numCellsCohesive = 0;

const int pylith::topology::MeshDataCohesiveTet4Level2::_cellDim = 3;

const int pylith::topology::MeshDataCohesiveTet4Level2::_numCorners = 4;

const int pylith::topology::MeshDataCohesiveTet4Level2::_numCornersCohesive = 9;

const double pylith::topology::MeshDataCohesiveTet4Level2::_vertices[] = {
  -1.000000e+00,      0.000000e+00,      0.000000e+00,
   0.000000e+00,     -1.000000e+00,      0.000000e+00,
   0.000000e+00,      0.000000e+00,      1.000000e+00,
   0.000000e+00,      1.000000e+00,      0.000000e+00,
   1.000000e+00,      0.000000e+00,      0.000000e+00,
   0.000000e+00,     -5.000000e-01,      5.000000e-01,
   0.000000e+00,      5.000000e-01,      5.000000e-01,
   0.000000e+00,      0.000000e+00,      0.000000e+00,
  -5.000000e-01,     -5.000000e-01,      0.000000e+00,
  -5.000000e-01,      0.000000e+00,      5.000000e-01,
  -5.000000e-01,      5.000000e-01,      0.000000e+00,
   5.000000e-01,     -5.000000e-01,      0.000000e+00,
   5.000000e-01,      5.000000e-01,      0.000000e+00,
   5.000000e-01,      0.000000e+00,      5.000000e-01,
};

const int pylith::topology::MeshDataCohesiveTet4Level2::_cells[] = {
  1,   8,  5,  7,
  5,   6,  7,  8,
  5,   8,  9,  6,
  2,   9,  6,  5,
  7,  10,  8,  6,
  3,  10,  7,  6,
  6,   9, 10,  8,
  0,   8, 10,  9,
  1,  11,  7,  5,
  7,   6,  5, 11,
  7,  11, 12,  6,
  3,  12,  6,  7,
  5,  13, 11,  6,
  2,  13,  5,  6,
  6,  12, 13, 11,
  4,  11, 13, 12,
};
const int pylith::topology::MeshDataCohesiveTet4Level2::_cellsCohesive[] = {
};
const int pylith::topology::MeshDataCohesiveTet4Level2::_materialIds[] = {
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
};

const int pylith::topology::MeshDataCohesiveTet4Level2::_numGroups = 4;

const int pylith::topology::MeshDataCohesiveTet4Level2::_groupSizes[] = {
  3, 3, 2, 6,
};

const int pylith::topology::MeshDataCohesiveTet4Level2::_groups[] = {
  0, 1, 8,
  2, 4, 13,
  0, 4,
  1, 2, 3, 5, 6, 7,
};

const char* pylith::topology::MeshDataCohesiveTet4Level2::_groupNames[] = {
  "edge 1",
  "edge 2",
  "end points",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveTet4Level2::_groupTypes[] = {
  "vertex", 
  "vertex",
  "vertex",
  "vertex",
};

pylith::topology::MeshDataCohesiveTet4Level2::MeshDataCohesiveTet4Level2(void)
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

pylith::topology::MeshDataCohesiveTet4Level2::~MeshDataCohesiveTet4Level2(void)
{}


// End of file
