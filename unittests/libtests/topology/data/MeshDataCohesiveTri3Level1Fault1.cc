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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "MeshDataCohesiveTri3Level1Fault1.hh"

const char* pylith::topology::MeshDataCohesiveTri3Level1Fault1::_filename = 
  "data/fourtri3.mesh";

const int pylith::topology::MeshDataCohesiveTri3Level1Fault1::_refineLevel = 1;
const char* pylith::topology::MeshDataCohesiveTri3Level1Fault1::_faultA = 
  "fault";
const char* pylith::topology::MeshDataCohesiveTri3Level1Fault1::_faultB = 0;

const int pylith::topology::MeshDataCohesiveTri3Level1Fault1::_numVertices = 18;

const int pylith::topology::MeshDataCohesiveTri3Level1Fault1::_spaceDim = 2;

const int pylith::topology::MeshDataCohesiveTri3Level1Fault1::_numCells = 16;

const int pylith::topology::MeshDataCohesiveTri3Level1Fault1::_numCellsCohesive = 4;

const int pylith::topology::MeshDataCohesiveTri3Level1Fault1::_cellDim = 2;

const int pylith::topology::MeshDataCohesiveTri3Level1Fault1::_numCorners = 3;

const int pylith::topology::MeshDataCohesiveTri3Level1Fault1::_numCornersCohesive = 4;

const int pylith::topology::MeshDataCohesiveTri3Level1Fault1::_matIdSum =
  8*1 + 8*2 + 4*100;

const int pylith::topology::MeshDataCohesiveTri3Level1Fault1::_numGroups = 4;

const int pylith::topology::MeshDataCohesiveTri3Level1Fault1::_groupSizes[4] = {
  11, 6, 2, 18,
};

const char* pylith::topology::MeshDataCohesiveTri3Level1Fault1::_groupNames[4] = {
  "edge 1",
  "edge 2",
  "end points",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveTri3Level1Fault1::_groupTypes[4] = {
  "vertex", 
  "vertex",
  "vertex",
  "vertex",
};

pylith::topology::MeshDataCohesiveTri3Level1Fault1::MeshDataCohesiveTri3Level1Fault1(void)
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
  matIdSum = _matIdSum;
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
} // constructor

pylith::topology::MeshDataCohesiveTri3Level1Fault1::~MeshDataCohesiveTri3Level1Fault1(void)
{}


// End of file
