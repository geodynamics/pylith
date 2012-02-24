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

#include "MeshDataCohesiveQuad4Level2.hh"

const char* pylith::topology::MeshDataCohesiveQuad4Level2::_filename = 
  "data/fourquad4.mesh";

const int pylith::topology::MeshDataCohesiveQuad4Level2::_refineLevel = 2;
const char* pylith::topology::MeshDataCohesiveQuad4Level2::_faultA = 0;
const char* pylith::topology::MeshDataCohesiveQuad4Level2::_faultB = 0;

const int pylith::topology::MeshDataCohesiveQuad4Level2::_numVertices = 25;

const int pylith::topology::MeshDataCohesiveQuad4Level2::_spaceDim = 2;

const int pylith::topology::MeshDataCohesiveQuad4Level2::_numCells = 16;

const int pylith::topology::MeshDataCohesiveQuad4Level2::_numCellsCohesive = 0;

const int pylith::topology::MeshDataCohesiveQuad4Level2::_cellDim = 2;

const int pylith::topology::MeshDataCohesiveQuad4Level2::_numCorners = 4;

const int pylith::topology::MeshDataCohesiveQuad4Level2::_numCornersCohesive = 6;

const double pylith::topology::MeshDataCohesiveQuad4Level2::_vertices[] = {
  -1.0, -1.0,
  -1.0,  0.0,
  -1.0,  1.0,
   0.0, -1.0,
   0.0,  0.0,
   0.0,  1.0,
   1.0, -1.0,
   1.0,  0.0,
   1.0,  1.0,
  -0.5, -1.0,
   0.0, -0.5,
  -0.5,  0.0,
  -1.0, -0.5,
  -0.5, -0.5,
   1.0, -0.5,
   0.5,  0.0,
   0.5, -1.0,
   0.5, -0.5,
  -1.0,  0.5,
   0.0,  0.5,
  -0.5,  1.0,
  -0.5,  0.5,
   0.5,  1.0,
   1.0,  0.5,
   0.5,  0.5,
};

const int pylith::topology::MeshDataCohesiveQuad4Level2::_cells[] = {
  16,  25,  29,  28,
  19,  26,  29,  25,
  17,  28,  29,  27,
  20,  27,  29,  26,
  22,  30,  33,  32,
  23,  31,  33,  30,
  19,  32,  33,  26,
  20,  26,  33,  31,
  18,  34,  37,  36,
  17,  27,  37,  34,
  21,  36,  37,  35,
  20,  35,  37,  27,
  24,  38,  40,  39,
  21,  35,  40,  38,
  23,  39,  40,  31,
  20,  31,  40,  35,
};
const int pylith::topology::MeshDataCohesiveQuad4Level2::_cellsCohesive[] = {
};
const int pylith::topology::MeshDataCohesiveQuad4Level2::_materialIds[] = {
  1, 1, 1, 1, 2, 2, 2, 2,
  1, 1, 1, 1, 2, 2, 2, 2,
};

const int pylith::topology::MeshDataCohesiveQuad4Level2::_numGroups = 3;

const int pylith::topology::MeshDataCohesiveQuad4Level2::_groupSizes[] = {
  5, 5, 5,
};

const int pylith::topology::MeshDataCohesiveQuad4Level2::_groups[] = {
  16, 19, 22, 25, 32,
  16, 17, 18, 28, 34,
  19, 20, 21, 26, 35,
};

const char* pylith::topology::MeshDataCohesiveQuad4Level2::_groupNames[] = {
  "edge 1",
  "end points",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveQuad4Level2::_groupTypes[] = {
  "vertex",
  "vertex",
  "vertex",
};

pylith::topology::MeshDataCohesiveQuad4Level2::MeshDataCohesiveQuad4Level2(void)
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

pylith::topology::MeshDataCohesiveQuad4Level2::~MeshDataCohesiveQuad4Level2(void)
{}


// End of file
