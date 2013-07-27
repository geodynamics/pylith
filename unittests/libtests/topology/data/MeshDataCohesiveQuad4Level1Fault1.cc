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

#include "MeshDataCohesiveQuad4Level1Fault1.hh"

const char* pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_filename = 
  "data/fourquad4.mesh";

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_refineLevel = 1;
const char* pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_faultA = "fault";
const char* pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_faultB = 0;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numVertices = 35;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_spaceDim = 2;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numCells = 16;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numCellsCohesive = 4;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_cellDim = 2;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numCorners = 4;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numCornersCohesive = 6;

const PylithScalar pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_vertices[] = {
  -1.0, -1.0,
  -1.0,  0.0,
  -1.0,  1.0,
   0.0, -1.0,
   0.0,  0.0,
   0.0,  1.0,
   1.0, -1.0,
   1.0,  0.0,
   1.0,  1.0,
   0.0, -1.0,
   0.0,  0.0,
   0.0,  1.0,
  -0.5, -1.0,
   0.0, -0.5,
  -0.5,  0.0,
  -1.0, -0.5,
  -0.5, -0.5,
   1.0, -0.5,
   0.5,  0.0,
   0.0, -0.5,
   0.5, -1.0,
   0.5, -0.5,
  -1.0,  0.5,
   0.0,  0.5,
  -0.5,  1.0,
  -0.5,  0.5,
   0.5,  1.0,
   0.0,  0.5,
   1.0,  0.5,
   0.5,  0.5,
   0.0, -1.0,
   0.0,  0.0,
   0.0,  1.0,
   0.0, -0.5,
   0.0,  0.5,
};

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_cells[] = {
  16,  28,  32,  31,
  19,  29,  32,  28,
  17,  31,  32,  30,
  20,  30,  32,  29,
  22,  33,  37,  36,
  23,  34,  37,  33,
  25,  36,  37,  35,
  26,  35,  37,  34,
  18,  38,  41,  40,
  17,  30,  41,  38,
  21,  40,  41,  39,
  20,  39,  41,  30,
  24,  42,  45,  44,
  27,  43,  45,  42,
  23,  44,  45,  34,
  26,  34,  45,  43,
};
const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_cellsCohesive[] = {
  19,  29,  25,  35,  46,  49,
  29,  20,  35,  26,  49,  47,
  20,  39,  26,  43,  47,  50,
  39,  21,  43,  27,  50,  48,
};
const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_materialIds[] = {
  1, 1, 1, 1, 2, 2, 2, 2,
  1, 1, 1, 1, 2, 2, 2, 2,
  100, 100, 100, 100,
};

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numGroups = 3;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_groupSizes[] = {
  6, 5, 15,
};

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_groups[] = {
  16, 19, 22, 25, 28, 36,
  16, 17, 18, 31, 38,
  19, 20, 21, 25, 26, 27, 29, 35, 39, 43, 46, 47, 48, 49, 50,
};

const char* pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_groupNames[] = {
  "edge 1",
  "end points",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_groupTypes[] = {
  "vertex",
  "vertex",
  "vertex",
};

pylith::topology::MeshDataCohesiveQuad4Level1Fault1::MeshDataCohesiveQuad4Level1Fault1(void)
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

pylith::topology::MeshDataCohesiveQuad4Level1Fault1::~MeshDataCohesiveQuad4Level1Fault1(void)
{}


// End of file
