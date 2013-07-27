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

#include "MeshDataCohesiveTri3Level1.hh"

const char* pylith::topology::MeshDataCohesiveTri3Level1::_filename = 
  "data/fourtri3.mesh";

const int pylith::topology::MeshDataCohesiveTri3Level1::_refineLevel = 1;
const char* pylith::topology::MeshDataCohesiveTri3Level1::_faultA = 0;
const char* pylith::topology::MeshDataCohesiveTri3Level1::_faultB = 0;

const int pylith::topology::MeshDataCohesiveTri3Level1::_numVertices = 13;

const int pylith::topology::MeshDataCohesiveTri3Level1::_spaceDim = 2;

const int pylith::topology::MeshDataCohesiveTri3Level1::_numCells = 16;

const int pylith::topology::MeshDataCohesiveTri3Level1::_numCellsCohesive = 0;

const int pylith::topology::MeshDataCohesiveTri3Level1::_cellDim = 2;

const int pylith::topology::MeshDataCohesiveTri3Level1::_numCorners = 3;

const int pylith::topology::MeshDataCohesiveTri3Level1::_numCornersCohesive = 6;

const PylithScalar pylith::topology::MeshDataCohesiveTri3Level1::_vertices[] = {
  -1.0,  0.0,
   0.0, -1.0,
   0.0,  0.0,
   0.0,  1.0,
   1.0,  0.0,
  -0.5, -0.5,
   0.0, -0.5,
  -0.5,  0.0,
   0.0,  0.5,
  -0.5,  0.5,
   0.5, -0.5,
   0.5,  0.0,
   0.5,  0.5,
};

const int pylith::topology::MeshDataCohesiveTri3Level1::_cells[] = {
  16,  21,  23,
  21,  17,  22,
  23,  22,  18,
  21,  22,  23,
  18,  24,  23,
  24,  19,  25,
  23,  25, 16,
  24,  25,  23,
  18,  22,  27,
  22,  17,  26,
  27,  26,  20,
  22,  26,  27,
  18,  27,  24,
  27,  20,  28,
  24,  28,  19,
  27,  28,  24,
};
const int pylith::topology::MeshDataCohesiveTri3Level1::_cellsCohesive[] = {
};
const int pylith::topology::MeshDataCohesiveTri3Level1::_materialIds[] = {
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
};

const int pylith::topology::MeshDataCohesiveTri3Level1::_numGroups = 4;

const int pylith::topology::MeshDataCohesiveTri3Level1::_groupSizes[] = {
  5, 3, 2, 5,
};

const int pylith::topology::MeshDataCohesiveTri3Level1::_groups[] = {
  16, 17, 19, 21, 25,
  17, 20, 26,
  16, 20,
  17, 18, 19, 22, 24,
};

const char* pylith::topology::MeshDataCohesiveTri3Level1::_groupNames[] = {
  "edge 1",
  "edge 2",
  "end points",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveTri3Level1::_groupTypes[] = {
  "vertex", 
  "vertex",
  "vertex",
  "vertex",
};

pylith::topology::MeshDataCohesiveTri3Level1::MeshDataCohesiveTri3Level1(void)
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

pylith::topology::MeshDataCohesiveTri3Level1::~MeshDataCohesiveTri3Level1(void)
{}


// End of file
