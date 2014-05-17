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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/* Original mesh
 *
 * Cells are 0-6 and vertices are 8-43
 */

#include "CohesiveDataHex8j.hh"

const int pylith::faults::CohesiveDataHex8j::_numVertices = 61;

const int pylith::faults::CohesiveDataHex8j::_spaceDim = 3;

const int pylith::faults::CohesiveDataHex8j::_numCells = 22;

const int pylith::faults::CohesiveDataHex8j::_cellDim = 3;

const int pylith::faults::CohesiveDataHex8j::_numCorners[22] = {
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  5,
  5,
};

const int pylith::faults::CohesiveDataHex8j::_materialIds[22] = {
  10, 10, 11, 11, 10, 10, 11, 11, 11, 11, 10, 10, 10, 10, 11, 11, 11, 11, 10, 10,
  1,  1,
};

const int pylith::faults::CohesiveDataHex8j::_numGroups = 3;

const int pylith::faults::CohesiveDataHex8j::_groupSizes[3] = {
  21+31+10, 5+4, 6+7+2 + 1+2+3 // vertices+edges+faces + split versions
};

const char* pylith::faults::CohesiveDataHex8j::_groupNames[3] = {
  "output", "fault_edge", "fault"
};

const char* pylith::faults::CohesiveDataHex8j::_groupTypes[3] = {
  "vertex", "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataHex8j::_filename = "data/hex8j.mesh";
const char* pylith::faults::CohesiveDataHex8j::_fault = "fault";
const char* pylith::faults::CohesiveDataHex8j::_edge = "fault_edge";

pylith::faults::CohesiveDataHex8j::CohesiveDataHex8j(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numCorners = const_cast<int*>(_numCorners);
  materialIds = const_cast<int*>(_materialIds);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  filename = const_cast<char*>(_filename);
  fault = const_cast<char*>(_fault);
  edge = const_cast<char*>(_edge);
} // constructor

pylith::faults::CohesiveDataHex8j::~CohesiveDataHex8j(void)
{}


// End of file
