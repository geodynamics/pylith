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

#include "MeshDataCohesiveQuad4Level1.hh"

const char* pylith::topology::MeshDataCohesiveQuad4Level1::_filename = 
  "data/fourquad4.mesh";

const int pylith::topology::MeshDataCohesiveQuad4Level1::_refineLevel = 1;
const char* pylith::topology::MeshDataCohesiveQuad4Level1::_faultA = 0;
const char* pylith::topology::MeshDataCohesiveQuad4Level1::_faultB = 0;

const int pylith::topology::MeshDataCohesiveQuad4Level1::_numVertices = 25;

const int pylith::topology::MeshDataCohesiveQuad4Level1::_spaceDim = 2;

const int pylith::topology::MeshDataCohesiveQuad4Level1::_numCells = 16;

const int pylith::topology::MeshDataCohesiveQuad4Level1::_numCellsCohesive = 0;

const int pylith::topology::MeshDataCohesiveQuad4Level1::_cellDim = 2;

const int pylith::topology::MeshDataCohesiveQuad4Level1::_numCorners = 4;

const int pylith::topology::MeshDataCohesiveQuad4Level1::_numCornersCohesive = 6;

const PylithScalar pylith::topology::MeshDataCohesiveQuad4Level1::_vertices[25*2] = {
  -1.0, -1.0,
  -1.0,  0.0,
  -1.0,  1.0,
   0.0, -1.0,
   0.0,  0.0,
   0.0,  1.0,
   1.0, -1.0,
   1.0,  0.0,
   1.0,  1.0,
  -0.5, -1.0, // 9
   0.0, -0.5,
  -0.5,  0.0,
  -1.0, -0.5,
   1.0, -0.5, // 13
   0.5,  0.0,
   0.5, -1.0,
  -1.0,  0.5, // 16
   0.0,  0.5,
  -0.5,  1.0,
   0.5,  1.0, // 19
   1.0,  0.5,
  -0.5, -0.5, // 21
   0.5, -0.5,
  -0.5,  0.5,
   0.5,  0.5,
};

const int pylith::topology::MeshDataCohesiveQuad4Level1::_cells[16*4] = {
  16,  25,  37,  28,
  25,  19,  26,  37,
  37,  26,  20,  27,
  28,  37,  27,  17,
  22,  29,  38,  31,
  29,  23,  30,  38,
  38,  30,  20,  26,
  31,  38,  26,  19,
  18,  32,  39,  34,
  32,  17,  27,  39,
  39,  27,  20,  33,
  34,  39,  33,  21,
  24,  35,  40,  36,
  35,  21,  33,  40,
  40,  33,  20,  30,
  36,  40,  30,  23,
};
const int pylith::topology::MeshDataCohesiveQuad4Level1::_cellsCohesive[] = {
};
const int pylith::topology::MeshDataCohesiveQuad4Level1::_materialIds[16] = {
  1, 1, 1, 1, 2, 2, 2, 2,
  1, 1, 1, 1, 2, 2, 2, 2,
};

const int pylith::topology::MeshDataCohesiveQuad4Level1::_numGroups = 3;

const int pylith::topology::MeshDataCohesiveQuad4Level1::_groupSizes[3] = {
  9, 9, 9,
};

const int pylith::topology::MeshDataCohesiveQuad4Level1::_groups[9+9+9] = { // vertices, edges
  16, 19, 22, 25, 31,    41, 42, 53, 54,
  16, 17, 18, 28, 32,    47, 48, 55, 56,
  19, 20, 21, 26, 33,    43, 44, 57, 58,
};

const char* pylith::topology::MeshDataCohesiveQuad4Level1::_groupNames[3] = {
  "edge 1",
  "end points",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveQuad4Level1::_groupTypes[3] = {
  "vertex",
  "vertex",
  "vertex",
};

pylith::topology::MeshDataCohesiveQuad4Level1::MeshDataCohesiveQuad4Level1(void)
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

pylith::topology::MeshDataCohesiveQuad4Level1::~MeshDataCohesiveQuad4Level1(void)
{}


// End of file
