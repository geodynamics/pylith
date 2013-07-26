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

#include "MeshDataCohesiveTet4Level1.hh"

const char* pylith::topology::MeshDataCohesiveTet4Level1::_filename = 
  "data/twotet4.mesh";

const int pylith::topology::MeshDataCohesiveTet4Level1::_refineLevel = 2;
const char* pylith::topology::MeshDataCohesiveTet4Level1::_faultA = 0;
const char* pylith::topology::MeshDataCohesiveTet4Level1::_faultB = 0;

const int pylith::topology::MeshDataCohesiveTet4Level1::_numVertices = 14;

const int pylith::topology::MeshDataCohesiveTet4Level1::_spaceDim = 3;

const int pylith::topology::MeshDataCohesiveTet4Level1::_numCells = 16;

const int pylith::topology::MeshDataCohesiveTet4Level1::_numCellsCohesive = 0;

const int pylith::topology::MeshDataCohesiveTet4Level1::_cellDim = 3;

const int pylith::topology::MeshDataCohesiveTet4Level1::_numCorners = 4;

const int pylith::topology::MeshDataCohesiveTet4Level1::_numCornersCohesive = 9;

const PylithScalar pylith::topology::MeshDataCohesiveTet4Level1::_vertices[] = {
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

const int pylith::topology::MeshDataCohesiveTet4Level1::_cells[] = {
  17,  24,  21,  23,
  21,  22,  23,  24,
  21,  24,  25,  22,
  18,  25,  22,  21,
  23,  26,  24,  22,
  19,  26,  23,  22,
  22,  25,  26,  24,
  16,  24,  26,  25,
  17,  27,  23,  21,
  23,  22,  21,  27,
  23,  27,  28,  22,
  19,  28,  22,  23,
  21,  29,  27,  22,
  18,  29,  21,  22,
  22,  28,  29,  27,
  20,  27,  29,  28,
};
const int pylith::topology::MeshDataCohesiveTet4Level1::_cellsCohesive[] = {
};
const int pylith::topology::MeshDataCohesiveTet4Level1::_materialIds[] = {
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
};

const int pylith::topology::MeshDataCohesiveTet4Level1::_numGroups = 4;

const int pylith::topology::MeshDataCohesiveTet4Level1::_groupSizes[] = {
  3, 3, 2, 6,
};

const int pylith::topology::MeshDataCohesiveTet4Level1::_groups[] = {
  16, 17, 24,
  18, 20, 29,
  16, 20,
  17, 18, 19, 21, 22, 23,
};

const char* pylith::topology::MeshDataCohesiveTet4Level1::_groupNames[] = {
  "edge 1",
  "edge 2",
  "end points",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveTet4Level1::_groupTypes[] = {
  "vertex", 
  "vertex",
  "vertex",
  "vertex",
};

pylith::topology::MeshDataCohesiveTet4Level1::MeshDataCohesiveTet4Level1(void)
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

pylith::topology::MeshDataCohesiveTet4Level1::~MeshDataCohesiveTet4Level1(void)
{}


// End of file
