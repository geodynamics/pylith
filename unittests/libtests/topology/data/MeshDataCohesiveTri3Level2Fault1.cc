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
   0.0, -0.5,
   0.0,  0.5,
};

const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_cells[] = {
     0, 11, 13,
    11, 12, 13,
     1, 12, 11, 
     2, 13, 12,
     2, 14, 13, 
    14, 15, 13, 
     3, 15, 14, 
     0, 13, 15, 
     6, 16, 18, 
    16, 17, 18, 
     5, 17, 16, 
     4, 18, 17, 
     6, 18, 20, 
    18, 19, 20, 
     4, 19, 18, 
     7, 20, 19, 

};
const int pylith::topology::MeshDataCohesiveTri3Level2Fault1::_cellsCohesive[] = {
   1,  12,   5,  16,   8,  21, 
  12,   2,  16,   6,  21,   9, 
   2,  14,   6,  20,   9,  22,
  14,   3,  20,   7,  22,  10, 

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
  0, 1, 3, 5, 7, 11, 15,
  1, 4, 5, 17,
  0, 4,
  1, 2, 3, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 21, 22,
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
