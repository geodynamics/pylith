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

/* Original mesh
 *
 * Cells are 0-6 and vertices are 8-43
 */

#include "CohesiveDataHex8i.hh"

const int pylith::faults::CohesiveDataHex8i::_numVertices = 36;

const int pylith::faults::CohesiveDataHex8i::_spaceDim = 3;

const int pylith::faults::CohesiveDataHex8i::_numCells = 9;

const int pylith::faults::CohesiveDataHex8i::_cellDim = 3;

const int pylith::faults::CohesiveDataHex8i::_numCorners[9] = {
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
  8,
};

const int pylith::faults::CohesiveDataHex8i::_materialIds[9] = {
  0,  0,  0,  0,  0,  0,  2,
  1,  1
};

const int pylith::faults::CohesiveDataHex8i::_numGroups = 2;

const int pylith::faults::CohesiveDataHex8i::_groupSizes[2] = {
  10+13+4, 12+14+4 // vertices+edges+faces 
};

const char* pylith::faults::CohesiveDataHex8i::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataHex8i::_groupTypes[2] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataHex8i::_filename = 
  "data/hex8i.mesh";

pylith::faults::CohesiveDataHex8i::CohesiveDataHex8i(void)
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
} // constructor

pylith::faults::CohesiveDataHex8i::~CohesiveDataHex8i(void)
{}


// End of file
