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

#include "MeshDataCohesiveQuad4Level1Fault1.hh"

const char* pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_filename = 
  "data/fourquad4.mesh";

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_refineLevel = 1;
const char* pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_faultA = "fault";
const char* pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_faultB = 0;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numVertices = 30;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_spaceDim = 2;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numCells = 16;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numCellsCohesive = 4;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_cellDim = 2;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numCorners = 4;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numCornersCohesive = 4;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_matIdSum = 
  8*1 + 8*2 + 4*100;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_numGroups = 3;

const int pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_groupSizes[3] = {
  6+4, // vertices+edges
  5+4, 
  10+8,
};

const char* pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_groupNames[3] = {
  "edge 1",
  "end points",
  "fault",
};

const char* pylith::topology::MeshDataCohesiveQuad4Level1Fault1::_groupTypes[3] = {
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
  matIdSum = _matIdSum;
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
} // constructor

pylith::topology::MeshDataCohesiveQuad4Level1Fault1::~MeshDataCohesiveQuad4Level1Fault1(void)
{}


// End of file
